# ============================================================================
# R/baselines/refit_cglasso_core.R — CGLasso's INDEPENDENT copy of the LNM refit core.
#
# Engineering-discipline separation (2026-07-22, user-mandated): TV and CGLasso do NOT
# share one refit file. This is a cg_-PREFIXED copy of R/refit.R (refit_* -> cg_refit_*)
# so both cores coexist in one R session (methods=c('TV','CGLasso')) with NO name collision
# (the P0-1 lesson). CGLasso wiring: refit_cglasso.R sources THIS (not R/refit.R) and calls
# cg_refit_fixed_support. A numerical fix in R/refit.R (TV) must be MIRRORED here.
# ============================================================================

# ============================================================================
# GLASSO-NATIVE relaxed (de-biased) refit for CGLasso (Tian et al. 2023).
#
# CGLasso is a glasso method (Yuan's Compo_glasso base fit = NR-latent-Z + glasso-Omega);
# its relaxed refit stays glasso-native. Pass 1 (penalized Compo_glasso) SELECTS the support;
# pass 2 (this file) re-estimates the FULL precision matrix on that FIXED zero pattern with NO
# penalty, together with the latent Z. Given Z, the Omega-update per slice is the Gaussian-
# graphical MLE with the pass-1 zero pattern held fixed = glasso(S_k, rho=0, zero=non-edges):
#   Dempster (1972) "Covariance Selection" (Biometrics 28:157) / Speed & Kiiveri (1986,
#   Ann.Stat. 14:138) / Hastie-Tibshirani-Friedman ESL 2e Alg 17.1 / glasso `zero`+`rho=0`
#   (Friedman-Hastie-Tibshirani 2008, Biostatistics 9:432).
# This CONSTRAINED MLE RE-ESTIMATES the diagonal (FREE -- the fitted covariance matches S on the
# diagonal AND the edges, the moment-matching characterization) + the active off-diagonals, and
# forces the non-edges to exactly 0. The support (edges) is UNCHANGED from pass 1 => the
# lambda-sweep ROC is refit-invariant; only the deployed magnitudes/BIC change. The refit
# alternates Yuan's OWN latent-Z update (cg_NR = Compo_glasso's Newton-Raphson, via
# cg_refit_renew_z_nr) with the glasso-zero Omega-update to convergence -- i.e. it is LITERALLY
# Compo_glasso's iteration on the fixed support. There is NO ADAM / NO reduced-coordinate hand-
# optimizer, and the Z-update is the base method's own solver (NOT the TV-inherited L-BFGS): CGLasso does NOT inherit
# TV's optimizer -- the Omega solve is delegated to glasso, exactly as the base method does.
#
# Model-selection likelihood = the JOINT LNM density evaluated at the INFERRED latent Zhat
# (multinomial + Gaussian-graphical layers), via cg_refit_joint_nll_average. NOT the Z_0
# Gaussian-only loglik. cg_refit_ic_unbalanced sits next to cg_refit_information_criteria for
# real data with UNEQUAL per-slice n (the balanced version asserts equal n).
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(MASS)
  library(Matrix)
  library(splines)
  library(glasso)          # constrained-MLE Omega-update (glasso `zero`+`rho=0`)
})

if (!exists("main_function_final", mode = "function") ||
    !exists("G_beta_Rcpp", mode = "function")) {
  source(here::here("R", "tvcglasso.R"))
}

.cg_refit_atomic_save <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- sprintf("%s.tmp.%d", path, Sys.getpid())
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("atomic save failed: ", path)
  }
  invisible(path)
}

.cg_refit_log <- function(path, fmt, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  line <- sprintf(fmt, ...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), line),
      file = path, append = TRUE)
  message(line)
  invisible(line)
}

# Canonical pooled option-2 pseudocount used by main_function_final_0624.R.
# ★ P0-5 (2026-07-22): the refit's multinomial likelihood uses the RAW observed counts.
# The refit gets its latent Z from the base fit (Z_start), so it never needs a pseudocount
# for an ALR init -> return the counts unchanged. (Previously this added the option-2
# pseudocount, leaking it into the refit likelihood; that is now confined to the base
# engine's Z_0 init only.)
cg_refit_prepare_counts <- function(X) {
  stopifnot(is.list(X), length(X) > 0L)
  X
}

cg_refit_support_from_beta <- function(beta, zero_tol = 0) {
  lapply(beta, function(B) {
    idx <- which(upper.tri(B) & abs(B) > zero_tol, arr.ind = TRUE)
    if (length(idx) == 0L) matrix(integer(0), nrow = 0L, ncol = 2L,
                                  dimnames = list(NULL, c("row", "col"))) else idx
  })
}

cg_refit_support_df <- function(support) sum(vapply(support, nrow, integer(1)))

cg_refit_support_signature <- function(support) {
  paste(vapply(seq_along(support), function(h) {
    idx <- support[[h]]
    if (nrow(idx) == 0L) sprintf("h%d:", h) else
      sprintf("h%d:%s", h, paste(sprintf("%d-%d", idx[, 1], idx[, 2]), collapse = ","))
  }, character(1)), collapse = "|")
}

# ★ codex P0-3 (2026-07-22): a content signature of the raw counts, bound into the checkpoint
# config so a resume with DIFFERENT data (same file path / same support signature) is DETECTED
# and refused (never silently continued on stale state). Per-slice dims + sum + sum-of-squares
# + a position-weighted sum catch any change in the data.
cg_refit_data_fp <- function(X_work) {
  vapply(X_work, function(X) {
    X <- as.matrix(X)
    c(nrow(X), ncol(X), sum(X), sum(X^2), sum(as.numeric(X) * seq_along(X)))
  }, numeric(5))
}

# ★ glasso-native constrained-MLE Omega-update (2026-07-22, user-mandated). CGLasso is a GLASSO
# method (Yuan's Compo_glasso base fit = NR-Z + glasso-Omega); its refit must stay glasso-native,
# NOT the ADAM reduced-parameterization it had inherited from the shared TV core. Given the latent
# Z, the Omega-update on the FIXED pass-1 support is the Gaussian-graphical MLE with that zero
# pattern held fixed = glasso(S_k, rho=0, zero=non-edges):
#   Dempster (1972) "Covariance Selection" (Biometrics 28:157) / Speed & Kiiveri (1986, Ann.Stat.
#   14:138) / Hastie-Tibshirani-Friedman ESL 2e Alg 17.1 / the glasso `zero`+`rho=0` route
#   (Friedman-Hastie-Tibshirani 2008, Biostatistics 9:432).
# The constrained MLE RE-ESTIMATES the diagonal (FREE -- the fitted covariance matches S on the
# diagonal AND on the edges) + the active off-diagonals, forces the non-edges to exactly 0, and is
# PD by construction. `zero_k` = the fixed non-edge upper-triangular index pairs (glasso enforces
# symmetry, so each pair is given once). Verified: moment-match on edges+diagonal ~1e-11, non-edges
# exactly 0. The rho=0 not-full-rank warning (n<P) is muffled; a support whose constrained MLE does
# not exist returns a non-PD Omega -> caught by the .all_slices_pd gate (that rho is not deployable).
# thr = 1e-4 = glasso's OWN DEFAULT (2026-07-25). A tighter hand-picked thr (we had 1e-8) is NOT
# reachable on the near-singular problems that arise at n<P dense supports, so glasso ground to
# the iteration cap there (measured: minutes per call, vs 0.05 s at the package default, with a
# moment-match difference of only ~1e-4 on Omega entries of order 0.1-1 -- immaterial). Use the
# package's own convergence control; do not hand-tighten a reproduced method's tolerances.
# maxit=500 (glasso's default is 10000): with the default thr the solve converges naturally well
# below 500 sweeps (measured: maxit=100 and 500 give IDENTICAL results), so this is only a backstop.
# Historical note: when the constrained MLE EXISTS glasso converges far
# below 500 sweeps; when it does NOT exist (dense support at n<P) glasso never converges no matter
# the cap and each sweep is expensive on the near-singular problem (measured: one such solve at
# P=15/n=12 exceeds 2 MINUTES at maxit=10000, x m slices x 70 rho = days). Capping only bounds the
# wasted work; the outcome is unchanged (PD -> usable, BIC computable; non-PD -> no BIC, hence not
# a candidate).
# ★★ RIDGED SLICE COVARIANCE in pass 2 (`refit_ridge_eps`, default 0.01; 2026-07-25, user-approved).
# OURS, must be disclosed (like the cg_NR iteration cap) -- it is NOT part of Yuan's CGLasso.
#
# WHY. The pass-2 estimator above is the given-graph Gaussian MLE, and that estimator has a DOMAIN:
# it exists iff the graph-restricted S admits a POSITIVE-DEFINITE COMPLETION [Grone-Johnson-de Sa-
# Wolkowicz 1984 LAA 58:109-124; Dempster 1972 Biometrics 28:157-175; Uhler 2012 Ann.Statist.
# 40(1):238-261 Thm 2.1]. Our per-slice S_k violates that criterion at BOTH ends of the rho path,
# in two different ways, and BOTH make the alternation non-terminating:
#   (A) DENSE end -- S_k is centered so rank(S_k) <= n_k-1 < P; once the pass-1 support is too dense
#       for that rank the MLE does not exist, the likelihood is literally UNBOUNDED (Uhler 2017
#       arXiv:1707.04345 Sec.3), and glasso grinds forever. MEASURED on the real stall points
#       (valid16 seed 1): cell14 rho-index 1 (n=12, rank(S)=11, 245 edges) and cell08 rho-index 1
#       (n=20, rank(S)=19, 281 edges) BOTH fail to return in 120 s -- and cannot even be interrupted,
#       since glasso's Fortran does not poll for R interrupts.
#   (B) SPARSE end -- one latent coordinate collapses within a slice (S_jj -> 0; for PSD S a zero
#       diagonal forces the whole row/column to zero), so no PD completion exists either. The empty
#       graph's MLE is diag(1/S_jj): MEASURED cell04 rho-index 45, min diag(S_k) = 1.30e-23 ->
#       max diag(Omega) = 7.71e22, cond = 7.2e23 -- yet min_eig = +0.107, i.e. POSITIVE DEFINITE.
#       Every PD-based guard (including .all_slices_pd below) is STRUCTURALLY BLIND to this, which
#       is why the divergence slid past the Omega-update and hung inside cg_NR instead.
#       NB the EMPTY graph's maximum likelihood threshold is 1, so n<P is mathematically incapable of
#       being the cause here; the two ends share ONE criterion but bind through different mechanisms.
#
# THE FIX. Solve pass 2 on S_k + eps*mean(diag(S_k))*I. A strictly positive-definite matrix IS its
# own positive-definite completion, so the constrained MLE exists and is unique for EVERY graph and
# EVERY n -- one repair at the root of the criterion, not two patches. Equivalently (exactly, since
# Omega_jj>0 at any PD solution so eps*sum|Omega_jj| = tr(eps*I*Omega)) it is a penalty eps on the
# DIAGONAL OF Omega ONLY: the selected off-diagonals stay COMPLETELY unpenalized, so the refit's
# de-biasing purpose is untouched. This is the ridge precision estimator of van Wieringen & Peeters
# (2016, CSDA 103:284-303; R pkg rags2ridges, ridgePchordal(S, lambda, zeros)) restricted to the
# selected support, and it is the SAME protection the penalized pass-1 already enjoys -- huge/glasso
# penalize the diagonal, which is exactly what yields Banerjee, El Ghaoui & d'Aspremont (2008, JMLR
# 9:485-516) Thm 1's bound p/lambda >= ||Omega||_2 with NO condition on n vs p. As a MAP under a
# conjugate inverse-Wishart/G-Wishart prior it is also how mclust cures the identical degeneracy in
# Gaussian mixtures (Fraley & Raftery 2007, J.Classification 24:155-181). Note glasso's own source
# documents the failure this repairs: "With rho=0, there may be convergence problems if the input
# matrix is not of full rank".
#
# CHOICE OF eps = 0.01 (relative to mean(diag(S_k)), so it follows the data scale = Ledoit-Wolf's
# scaled-identity target mu*I, mu = tr(S)/p). MEASURED on our own geometry:
#   * fixes both ends -- dense: 120s-no-return -> 0.02 s / 0.01 s, PD, max|Omega| 13.7 / 8.0;
#     sparse: max|Omega| 7.71e22 -> 37.5 (cond 7.2e23 -> 352) and 1.74e10 -> 18.1.
#   * eps=1e-4 is NOT enough (still 18.5 s on the dense case); the 1/eps scaling is exact
#     (max diag(Omega) = 3.7489/eps to 5 s.f. over eps = 1e-4 .. 1e-1).
#   * COST: eps=0.01 retains 97.9% of the exact refit's off-diagonal magnitude (1.7% relative
#     Frobenius) against the ~30% L1 shrinkage the refit exists to remove; eps=0.1 retains only 82%.
#   * do NOT use an analytic shrinkage intensity: Ledoit-Wolf / Schaefer-Strimmer give 0.81-0.84 at
#     our n_k, two orders too large -- the refit would then recover only ~19% of the true |Omega_ij|.
# HONEST SCOPE. Every ingredient above is published; the COMPOSITION ("ridge the slice covariance
# before the given-graph refit") we did not find stated in any paper -> disclose it as ours. The
# bound Omega_jj <= 1/S~_jj holds only for the EMPTY graph; for a general support the valid
# statement is ||Omega|| <= (#cliques + #separators)/eps on decomposable G (measured ~1/(2eps)).
# Setting refit_ridge_eps = 0 restores the exact unridged given-graph MLE (and the stalls).
cg_refit_glasso_slice <- function(S_k, zero_k, ridge_eps = 0.01) {
  if (ridge_eps > 0) S_k <- S_k + ridge_eps * mean(diag(S_k)) * diag(nrow(S_k))
  fit <- withCallingHandlers(
    if (nrow(zero_k) > 0L)
      glasso::glasso(S_k, rho = 0, zero = zero_k, penalize.diagonal = FALSE,
                     thr = 1e-4, maxit = 500L)
    else
      glasso::glasso(S_k, rho = 0, penalize.diagonal = FALSE, thr = 1e-4, maxit = 500L),
    warning = function(w) if (grepl("rho=0|full rank", conditionMessage(w)))
      invokeRestart("muffleWarning"))
  (fit$wi + t(fit$wi)) / 2   # ONE canonical symmetric Omega (glasso BCD can be mildly asymmetric)
}

cg_refit_validate_fixed_coordinates <- function(beta, beta_template, support,
                                              tol = 0, free_diag = FALSE) {
  P <- nrow(beta_template[[1]])
  # ★ free_diag (2026-07-22): the constrained-MLE refit re-estimates the diagonal (Dempster 1972
  # covariance selection / ESL Alg 17.1 / glasso zero+rho=0), so the diagonal is NO LONGER bit-
  # identical to the pass-1 template; do not require it. (When free_diag, the diagonal is IN the
  # optimization support -> `inactive_ok`/`diag(allowed)` still hold; only the diag-bit-identical
  # requirement is dropped.)
  diag_ok <- all(vapply(seq_along(beta), function(h) {
    identical(as.numeric(diag(beta[[h]])), as.numeric(diag(beta_template[[h]])))
  }, logical(1)))
  symmetric_ok <- all(vapply(beta, function(B) max(abs(B - t(B))) <= tol, logical(1)))
  inactive_ok <- all(vapply(seq_along(beta), function(h) {
    allowed <- matrix(FALSE, P, P)
    idx <- support[[h]]
    if (nrow(idx) > 0L) {
      allowed[idx] <- TRUE
      allowed[cbind(idx[, 2], idx[, 1])] <- TRUE
    }
    diag(allowed) <- TRUE
    all(beta[[h]][!allowed] == 0)
  }, logical(1)))
  list(ok = (free_diag || diag_ok) && symmetric_ok && inactive_ok,
       diagonal_bit_identical = diag_ok,
       symmetric = symmetric_ok,
       inactive_exact_zero = inactive_ok)
}

cg_refit_active_gradient_blocks <- function(beta, Z, basis, support) {
  m <- length(Z)
  P <- nrow(beta[[1]])
  nk <- vapply(Z, nrow, integer(1)); w_slice <- nk / sum(nk)   # P0-6: n_k/N (=1/m for equal n)
  Om <- G_beta_Rcpp(beta, basis, m)
  S <- S_Z_t(Z)
  lapply(seq_along(beta), function(h) {
    idx <- support[[h]]
    if (nrow(idx) == 0L) return(numeric(0))
    G <- matrix(0, P, P)
    for (k in seq_len(m)) {
      G <- G + w_slice[k] * basis[k, h] * (S[[k]] - .inv_pd(Om[[k]]))
    }
    as.numeric(G[idx])
  })
}

cg_refit_diagonal_score <- function(beta, Z, basis) {
  m <- length(Z)
  P <- nrow(beta[[1]])
  nk <- vapply(Z, nrow, integer(1)); w_slice <- nk / sum(nk)   # P0-6: n_k/N (=1/m for equal n)
  Om <- G_beta_Rcpp(beta, basis, m)
  S <- S_Z_t(Z)
  scores <- lapply(seq_along(beta), function(h) {
    g <- numeric(P)
    for (k in seq_len(m)) {
      # A diagonal beta coefficient enters Omega only once, hence the 1/2.
      g <- g + w_slice[k] * 0.5 * basis[k, h] * diag(S[[k]] - .inv_pd(Om[[k]]))
    }
    g
  })
  list(by_basis = scores,
       max_abs = max(abs(unlist(scores))),
       l2 = sqrt(sum(unlist(scores)^2)))
}

# The unpenalized complete/joint LNM criterion, with constants common to all
# candidate supports omitted.  This is not the integrated observed-data
# likelihood; it is used only for same-data, same-nuisance path comparisons.
# Per-OBSERVATION-averaged joint LNM neg-loglik (multinomial + Gaussian-graphical).
# ★ P0-6 (2026-07-22): each slice k is weighted by n_k (its sample size), so
#   2*N*total == the TRUE joint -2loglik  Sigma_k Sigma_i q_ki  for ANY n_k (unequal-n
#   real data works by default). For EQUAL n it reduces EXACTLY to the old slice-mean
#   (n_k/N = 1/m), so all equal-n simulations are bit-unchanged.
# ★ P0-5: X_work must be the RAW observed counts (the caller passes raw X; the pseudocount
#   is used only for the base engine's Z_0 ALR init, never here).
cg_refit_joint_nll_average <- function(X_work, Z, beta, basis) {
  m <- length(X_work)
  P <- ncol(X_work[[1]]) - 1L
  n_k <- vapply(Z, nrow, integer(1)); N <- sum(n_k)
  Om <- G_beta_Rcpp(beta, basis, m)
  logdets <- vapply(Om, .logdet_chol, numeric(1))
  if (anyNA(logdets)) {
    return(list(total = Inf, multinomial = Inf, neg_logdet = Inf,
                trace = Inf, by_slice = rep(Inf, m), Omega = Om))
  }
  # ★ P0-6/per-slice (2026-07-22): per-slice joint NLL contributions (multinomial +
  # Gaussian-graphical), SUMMED over that slice's OWN observations (NOT divided by N).
  # sum(by_slice) == N * total, and 2 * by_slice[k] is slice k's EXACT contribution to the
  # joint -2loglik -> enables PER-SLICE INDEPENDENT model selection (CGLasso, basis = I_m,
  # slices conditionally separable). ADDITIVE: the pooled total/multinomial/neg_logdet/trace
  # below are computed bit-identically to before, so TV/JGL callers are unaffected.
  S <- S_Z_t(Z)
  mult_k <- vapply(seq_len(m), function(k) {
    Xi <- X_work[[k]]
    Zi <- Z[[k]]
    M <- rowSums(Xi)
    -sum(rowSums(Xi[, seq_len(P), drop = FALSE] * Zi) -
           M * log1p(rowSums(exp(Zi))))
  }, numeric(1))
  gauss_k <- -0.5 * n_k * logdets +
    0.5 * n_k * vapply(seq_len(m), function(k) sum(S[[k]] * Om[[k]]), numeric(1))
  by_slice <- mult_k + gauss_k
  multinomial <- sum(mult_k) / N
  neg_logdet <- -0.5 * sum(n_k * logdets) / N
  trace <- 0.5 * sum(n_k * vapply(seq_len(m), function(k) {
    sum(S[[k]] * Om[[k]])
  }, numeric(1))) / N
  list(total = multinomial + neg_logdet + trace,
       multinomial = multinomial,
       neg_logdet = neg_logdet,
       trace = trace,
       by_slice = by_slice,
       Omega = Om)
}

# Exact smooth Z-block for the same joint logistic-normal-multinomial
# criterion used above. Slices are conditionally separable when Omega is
# fixed; within a slice, centering by colMeans(Z) is handled analytically in
# both the objective and its gradient.
cg_refit_z_slice_value_gradient <- function(par, X, Omega) {
  n <- nrow(X)
  P <- ncol(X) - 1L
  Z <- matrix(par, nrow = n, ncol = P)
  M <- rowSums(X)
  row_max <- pmax(0, apply(Z, 1L, max))
  ez <- exp(sweep(Z, 1L, row_max, "-"))
  den <- exp(-row_max) + rowSums(ez)
  log_den <- row_max + log(den)
  prob <- ez / den
  centered <- sweep(Z, 2L, colMeans(Z), "-")
  value <- mean(-rowSums(X[, seq_len(P), drop = FALSE] * Z) +
                  M * log_den) +
    0.5 * mean(rowSums((centered %*% Omega) * centered))
  gradient <- (prob * M - X[, seq_len(P), drop = FALSE] +
                 centered %*% Omega) / n
  list(value = value, gradient = as.numeric(gradient))
}

cg_refit_z_gradient_stats <- function(X_work, Z, beta, basis) {
  Omega <- G_beta_Rcpp(beta, basis, length(Z))
  by_slice <- lapply(seq_along(Z), function(k) {
    g <- cg_refit_z_slice_value_gradient(as.numeric(Z[[k]]), X_work[[k]],
                                      Omega[[k]])$gradient
    matrix(g, nrow(Z[[k]]), ncol(Z[[k]])) / length(Z)
  })
  flat <- unlist(by_slice)
  list(by_slice = by_slice, max_abs = max(abs(flat)),
       l2 = sqrt(sum(flat^2)), rms = sqrt(mean(flat^2)))
}

# ★ Z-update = Yuan's OWN Compo_glasso Newton-Raphson (cg_NR), applied per slice (2026-07-23,
# user-mandated: the baseline comparison must respect the original method). This makes the refit
# iteration LITERALLY Compo_glasso on the fixed support (alternate cg_NR-latent-Z <-> glasso-Omega),
# not the TV-inherited L-BFGS. cg_NR (CompoGlasso.R) is a per-slice damped Newton with an Armijo
# line search that solves the LNM latent MAP given Omega; the subproblem is convex so the solution
# is the unique MAP (independent of the solver). Returns the outer machinery's Z-update CONTRACT
# (Z / usable / all_converged / max_abs_gradient / diagnostics); `...` is accepted + ignored.
cg_refit_renew_z_nr <- function(X_work, Z_start, Omega_list, ...) {
  m <- length(X_work)
  Z <- lapply(seq_len(m), function(k)
    cg_NR(as.matrix(X_work[[k]]), as.matrix(Z_start[[k]]), Omega_list[[k]]))
  usable <- all(vapply(Z, function(z) all(is.finite(z)), logical(1)))
  gs <- if (usable) cg_refit_z_gradient_stats(X_work, Z, Omega_list, diag(m))
        else list(by_slice = vector("list", m), max_abs = NA_real_)
  per_slice_grad <- if (usable) vapply(gs$by_slice, function(g) max(abs(g)), numeric(1))
                    else rep(NA_real_, m)
  list(Z = Z,
       diagnostics = data.frame(slice = seq_len(m), convergence = 0L,   # cg_NR iterates each sample to its own step threshold (Yuan's convention)
                                max_abs_gradient = per_slice_grad),
       usable = usable && is.finite(gs$max_abs),
       all_converged = TRUE,
       max_abs_gradient = gs$max_abs)
}

# ★ P0-6 (2026-07-22): `nll_average` from cg_refit_joint_nll_average is now the n_k-weighted
# per-observation joint neg-loglik, so 2*N*nll_average IS the true joint -2loglik for ANY
# n_k (equal OR unequal). This ONE function is correct for both sims (equal n) and real
# Zebrafish (unequal n) -- no assertion, no separate unbalanced variant. df = # nonzero
# off-diagonal beta coefficients; N = sum(n_per_slice).
cg_refit_information_criteria <- function(nll_average, n_per_slice, df, P) {
  N <- sum(n_per_slice)
  neg2loglik <- 2 * N * nll_average
  c(neg2loglik = neg2loglik,
    AIC = neg2loglik + 2 * df,
    BIC = neg2loglik + log(N) * df,
    eBIC = neg2loglik + (log(N) + log(P)) * df)
}
# Back-compat alias (the unequal-n case is now handled by cg_refit_information_criteria itself).
cg_refit_ic_unbalanced <- cg_refit_information_criteria

cg_refit_fixed_support <- function(X_work, Z_start, beta_start, basis,
                                checkpoint_file, progress_log,
                                max_outer = 50L, inner_max = 10L,   # 50 = Compo_glasso's own max_iter
                                z_align_max = 160L,
                                initial_learning_rate = 0.01,
                                conv_tol_Omega = 1e-4,
                                conv_tol_Z = 5e-5,
                                conv_tol_z_grad = 1e-4,
                                conv_tol_grad = 1e-4,
                                conv_sustain = 1L,   # Yuan's while-loop exits as soon as the condition holds (no 'sustained K' requirement)
                                # ★ 2026-07-24 (user-mandated: reproduce the ORIGINAL method's design).
                                # Yuan's Compo_glasso declares convergence with RELATIVE thresholds
                                # derived from the FIRST iteration's movement, on the PARAMETER CHANGE
                                # only (never on a gradient):
                                #     Omega.1 <- matrix(0,K,K); z.end <- matrix(0,n,K)   # <- ZERO at this point
                                #     O_thr <- mean((Omega.0-Omega.1)^2)/O_ratio   == mean(Omega_init^2)/1000
                                #     z_thr <- mean((z.start-z.end)^2)/z_ratio     == mean(z_init^2)/1000
                                #     while ((dOmega > O_thr || dZ > z_thr) && iter <= max_iter)  # max_iter = 50
                                # Because the "end" iterates are ZERO when the thresholds are formed, the rule is
                                # "mean-squared CHANGE <= mean-squared MAGNITUDE of the initial value / ratio",
                                # i.e. RELATIVE-TO-MAGNITUDE. In our normalized units (dZ = ||dZ||_F/||Z||_F)
                                # that is exactly dZ <= 1/sqrt(z_ratio) ~ 3.16% for ratio=1000 -- which is why
                                # Yuan's max_iter=50 suffices.
                                # Our previous ABSOLUTE tolerances (dZ<5e-5) PLUS a self-invented GRADIENT
                                # condition were far stricter than the method being reproduced, and were
                                # unreachable at low depth / large P (the latent-Z jitter floor) -> every rho
                                # rode the iteration cap. Adopting Yuan's own rule restores the original
                                # method's convergence semantics AND removes that pathology.
                                yuan_conv = TRUE,   # use Compo_glasso's relative parameter-change rule
                                z_ratio = 1000, O_ratio = 1000,
                                pd_max_backtrack = 25L,
                                pd_c1 = 1e-4,
                                pd_nonmono_K = 5L,
                                # ★ OURS, disclosed: ridge on the slice covariance in pass 2 ONLY
                                # (see the block above cg_refit_glasso_slice). 0 = exact unridged
                                # given-graph MLE (and the dense-/sparse-end non-termination).
                                refit_ridge_eps = 0.01,
                                free_diag = TRUE) {
  # NOTE (2026-07-22): the Omega-update is now the glasso-native constrained MLE (see
  # cg_refit_glasso_slice); `inner_max`/`initial_learning_rate`/`inner_grad_tol`/`pd_*` (the old
  # ADAM/line-search knobs) are retained for call-signature compatibility but UNUSED by this update.
  # free_diag = TRUE: the constrained MLE re-estimates the diagonal (Dempster/ESL) -> the diagonal
  # is no longer bit-identical to the pass-1 template, so the fixed-coordinate check does not
  # require it (only the non-edge off-diagonals stay exactly 0).
  m <- length(X_work)
  P <- nrow(beta_start[[1]])
  support <- cg_refit_support_from_beta(beta_start)        # ACTIVE off-diagonals (edges) = df; diagonal is NOT here (glasso frees it)
  # fixed non-edge pattern per slice (upper-tri off-diagonals forced to 0 by the pass-1 threshold);
  # given ONCE to glasso as its `zero` constraint at every refit iteration.
  zero_pattern <- lapply(beta_start, function(B) {
    ze <- which(upper.tri(B) & B == 0, arr.ind = TRUE)
    if (length(ze) == 0L) matrix(integer(0), nrow = 0L, ncol = 2L) else ze
  })
  signature <- cg_refit_support_signature(support)
  beta_template <- beta_start
  initial_check <- cg_refit_validate_fixed_coordinates(beta_start, beta_template, support, free_diag = free_diag)
  stopifnot(initial_check$ok,
            .all_slices_pd(G_beta_Rcpp(beta_start, basis, m)))

  config <- list(P = P, m = m, J_n = length(beta_start),
                 df = cg_refit_support_df(support), signature = signature,
                 data_fp = cg_refit_data_fp(X_work),   # codex P0-3: bind the RAW-COUNTS content -> a different-data resume is refused
                 max_outer = max_outer, inner_max = inner_max,
                 z_align_max = z_align_max,
                 z_method = "joint_slice_lbfgsb",
                 initial_learning_rate = initial_learning_rate,
                 conv_tol_Omega = conv_tol_Omega,
                 conv_tol_Z = conv_tol_Z,
                 conv_tol_z_grad = conv_tol_z_grad,
                 conv_tol_grad = conv_tol_grad,
                 conv_sustain = conv_sustain,
                 yuan_conv = yuan_conv, z_ratio = z_ratio, O_ratio = O_ratio,
                 pd_max_backtrack = pd_max_backtrack,
                 pd_c1 = pd_c1, pd_nonmono_K = pd_nonmono_K,
                 # in `config` ON PURPOSE: the ridge CHANGES the pass-2 estimator, so a checkpoint
                 # written under a different eps must be REFUSED, not silently resumed.
                 refit_ridge_eps = refit_ridge_eps)

  if (file.exists(checkpoint_file)) {
    state <- readRDS(checkpoint_file)
    if (!identical(state$config, config)) {
      extension_fields <- setdiff(names(config), "max_outer")
      max_outer_extension <- setequal(names(state$config), names(config)) &&
        all(vapply(extension_fields, function(nm) {
          identical(state$config[[nm]], config[[nm]])
        }, logical(1))) &&
        isTRUE(config$max_outer > state$config$max_outer) &&
        isTRUE(config$max_outer > state$cg_refit_iter)
      legacy_fields <- c(
        "P", "m", "J_n", "df", "signature", "data_fp", "max_outer", "inner_max",
        "initial_learning_rate", "conv_tol_Omega", "conv_tol_Z",
        "conv_tol_grad", "conv_sustain",
        "pd_max_backtrack", "pd_c1", "pd_nonmono_K")
      legacy_compatible <- identical(state$phase, "z_align") &&
        identical(state$cg_refit_iter, 0L) &&
        !"z_method" %in% names(state$config) &&
        all(vapply(legacy_fields, function(nm) {
          identical(state$config[[nm]], config[[nm]])
        }, logical(1))) && z_align_max > state$z_iter
      if (max_outer_extension) {
        old_max_outer <- state$config$max_outer
        state$config <- config
        .cg_refit_atomic_save(state, checkpoint_file)
        .cg_refit_log(progress_log,
                   "extend checkpoint max_outer from %d to %d at cg_refit_iter=%d",
                   old_max_outer, config$max_outer, state$cg_refit_iter)
      } else if (!legacy_compatible) {
        stop("checkpoint config mismatch: ", checkpoint_file)
      } else {
        old_iter <- state$z_iter
        state$config <- config
        state$version <- "2026-07-16-refit-v2"
        state$z_streak <- 0L
        state$z_align_converged <- FALSE
        if (nrow(state$trace) && !"max_z_grad" %in% names(state$trace)) {
          state$trace$max_z_grad <- NA_real_
        }
        if (is.null(state$z_optimizer_stats)) state$z_optimizer_stats <- list()
        .cg_refit_atomic_save(state, checkpoint_file)
        .cg_refit_log(progress_log,
                   paste0("migrate legacy Z checkpoint at iter=%d; retain Z and ",
                          "continue with exact joint-slice Z block"), old_iter)
      }
    }
    .cg_refit_log(progress_log, "resume phase=%s z_iter=%d cg_refit_iter=%d df=%d",
               state$phase, state$z_iter, state$cg_refit_iter, config$df)
  } else {
    state <- list(version = "2026-07-16-refit-v2", config = config,
                  phase = "z_align", Z = Z_start, beta = beta_start,
                  z_iter = 0L, cg_refit_iter = 0L,
                  z_streak = 0L, z_align_converged = FALSE,
                  cg_refit_streak = 0L, trace = data.frame(),
                  optimizer_stats = list(), z_optimizer_stats = list())
    .cg_refit_atomic_save(state, checkpoint_file)
    .cg_refit_log(progress_log, "start z-alignment df=%d", config$df)
  }

  Omega_fixed <- G_beta_Rcpp(beta_start, basis, m)
  if (identical(state$phase, "z_align")) {
    start <- state$z_iter + 1L
    if (start <= z_align_max) {
      for (iter in seq.int(start, z_align_max)) {
        z_update <- cg_refit_renew_z_nr(X_work, state$Z, Omega_fixed)   # Yuan's cg_NR (base method's own Z-update)
        if (!z_update$usable) {
          state$last_z_optimizer_failure <- z_update$diagnostics
          .cg_refit_atomic_save(state, checkpoint_file)
          stop("exact Z optimizer returned a non-finite alignment state: ",
               checkpoint_file)
        }
        if (!z_update$all_converged) {
          .cg_refit_log(progress_log,
                     "z-align iter=%d retains finite optimizer intermediate; nonzero slices=%s",
                     iter,
                     paste(z_update$diagnostics$slice[
                       z_update$diagnostics$convergence != 0L], collapse = ","))
        }
        Z_new <- z_update$Z
        dZ <- sqrt(sum(mapply(function(A, B) sum((A - B)^2), Z_new, state$Z))) /
          (sqrt(sum(vapply(state$Z, function(M) sum(M^2), numeric(1)))) + 1e-12)
        max_z_grad <- z_update$max_abs_gradient
        # ★ Yuan's LITERAL rule (Compo_glasso). In his code the thresholds are formed while the
        # "end" iterates are still ZERO matrices (z.end <- matrix(0,n,K); Omega.1 <- matrix(0,K,K)):
        #     z_thr <- mean((z.start - z.end)^2)/z_ratio   ==  mean(z_init^2)/z_ratio
        # so the rule is "MEAN-SQUARED CHANGE <= MEAN-SQUARED MAGNITUDE of the initial value / ratio",
        # i.e. a RELATIVE-TO-MAGNITUDE criterion, NOT relative to the first iteration's change.
        # In our normalized units (dZ = ||dZ||_F/||Z||_F) that is exactly dZ <= 1/sqrt(z_ratio)
        # (mean(dZ^2) <= mean(Z^2)/r  <=>  ||dZ||^2/||Z||^2 <= 1/r  <=>  dZ <= r^-1/2),
        # i.e. ~3.16% relative change for r=1000 -- which is why Yuan's max_iter=50 suffices.
        # Parameter change only: Compo_glasso never tests a gradient. (The absolute+gradient rule
        # is kept as a fallback for yuan_conv = FALSE.)
        z_thr_rel <- 1 / sqrt(z_ratio)
        state$z_streak <- if (isTRUE(yuan_conv)) {
          if (dZ <= z_thr_rel) state$z_streak + 1L else 0L
        } else if (dZ < conv_tol_Z && max_z_grad < conv_tol_z_grad) {
          state$z_streak + 1L
        } else 0L
        nll <- cg_refit_joint_nll_average(X_work, Z_new, beta_start, basis)$total
        state$trace <- rbind(state$trace,
                             data.frame(phase = "z_align", iter = iter,
                                        objective = nll, dOmega = 0, dZ = dZ,
                                        max_active_grad = NA_real_,
                                        max_z_grad = max_z_grad,
                                        min_eig = min(vapply(Omega_fixed, .min_eig, numeric(1))),
                                        streak = state$z_streak))
        state$Z <- Z_new
        state$z_iter <- iter
        state$z_optimizer_stats[[length(state$z_optimizer_stats) + 1L]] <-
          list(phase = "z_align", iter = iter,
               diagnostics = z_update$diagnostics)
        .cg_refit_atomic_save(state, checkpoint_file)
        .cg_refit_log(progress_log,
                   paste0("z-align iter=%d obj=%.8f dZ=%.3e ",
                          "max|g_Z|=%.3e streak=%d"),
                   iter, nll, dZ, max_z_grad, state$z_streak)
        if (state$z_streak >= conv_sustain) break
      }
    }
    state$z_align_converged <- state$z_streak >= conv_sustain
    .cg_refit_atomic_save(state, checkpoint_file)
    if (!isTRUE(state$z_align_converged)) {
      # ★ 2026-07-23 (user-mandated): do NOT hard-error on a strict-gate miss. At low sequencing
      # depth the latent-Z jitter floor sits above conv_tol_Z/conv_tol_z_grad, so the z-alignment
      # cannot reach the strict sustained gate within z_align_max -- but the Z it reached IS stable
      # (dZ tiny; only the gradient oscillates above the strict floor). FLAG z_align_converged=FALSE
      # and CONTINUE with that Z into the refit phase, letting the downstream deployment gate decide
      # what to do with the flag -- instead of throwing the whole rho's refit away (which also lost
      # the still-valid penalized pass-1 `pre` candidate). Previously this stop() made EVERY rho of
      # every low-depth cell error out -> empty/incomplete deployed points.
      .cg_refit_log(progress_log,
                 paste0("z-alignment did NOT reach the strict gate by iter=%d (low-depth ",
                        "Z-jitter floor); CONTINUING with the current Z, flagged z_align_converged=FALSE"),
                 state$z_iter)
    }
    state$pre <- list(Z = state$Z, beta = beta_start,
                      criterion = cg_refit_joint_nll_average(X_work, state$Z,
                                                           beta_start, basis))
    state$phase <- "refit"
    state$beta <- beta_start
    state$cg_refit_iter <- 0L
    state$cg_refit_streak <- 0L
    .cg_refit_atomic_save(state, checkpoint_file)
    .cg_refit_log(progress_log, "z-alignment complete iter=%d pre_obj=%.8f",
               state$z_iter, state$pre$criterion$total)
  }

  exit_reason <- if (state$cg_refit_streak >= conv_sustain) {
    "converged"
  } else "max_outer"
  post_not_pd <- FALSE   # set TRUE if the constrained MLE turns non-PD (support not estimable at this n)
  start <- state$cg_refit_iter + 1L
  if (!identical(exit_reason, "converged") && start <= max_outer) {
    for (iter in seq.int(start, max_outer)) {
      beta_old <- state$beta
      Z_old <- state$Z
      Omega_old <- G_beta_Rcpp(beta_old, basis, m)
      z_update <- cg_refit_renew_z_nr(X_work, Z_old, Omega_old)   # Yuan's cg_NR (base method's own Z-update)
      if (!z_update$usable) {
        state$last_z_optimizer_failure <- z_update$diagnostics
        .cg_refit_atomic_save(state, checkpoint_file)
        stop("exact Z optimizer returned a non-finite refit state: ",
             checkpoint_file)
      }
      if (!z_update$all_converged) {
        .cg_refit_log(progress_log,
                   "refit iter=%d retains finite optimizer intermediate; nonzero slices=%s",
                   iter,
                   paste(z_update$diagnostics$slice[
                     z_update$diagnostics$convergence != 0L], collapse = ","))
      }
      Z_new <- z_update$Z
      # glasso-native constrained-MLE Omega-update on the FIXED support (see cg_refit_glasso_slice):
      # per slice, glasso(S_k, rho=0, zero=non-edges) re-estimates the diagonal + active off-diagonals
      # and forces the non-edges to exactly 0. basis = I_m => beta_k = Omega_k.
      S_new <- S_Z_t(Z_new)
      Omega_new <- lapply(seq_len(m), function(k) cg_refit_glasso_slice(S_new[[k]], zero_pattern[[k]], refit_ridge_eps))
      beta_new <- Omega_new
      if (!.all_slices_pd(Omega_new)) {
        # ★ 2026-07-23 (user-mandated, same flag-not-throw policy as the z-align gate): the
        # unpenalized constrained MLE glasso(rho=0, zero=...) does NOT exist when the pass-1 support
        # is not estimable at this n (e.g. a DENSE support with n<P) -> glasso returns a non-PD Omega.
        # Do NOT throw: mark the post invalid, KEEP the valid penalized `pre`, and stop the loop.
        # Sparser rho whose MLE exists still succeed; the deployment gate excludes this rho's post
        # via converged=FALSE. Previously this stopifnot errored the whole rho (losing the pre too).
        post_not_pd <- TRUE
        exit_reason <- "post_not_pd"
        .cg_refit_log(progress_log,
                   "refit iter=%d: constrained MLE non-PD (support not estimable at this n); post invalid, pre retained",
                   iter)
        break
      }
      fixed_check <- cg_refit_validate_fixed_coordinates(beta_new, beta_template, support, free_diag = free_diag)
      stopifnot(fixed_check$ok)

      dOm <- sqrt(sum(mapply(function(A, B) sum((A - B)^2), Omega_new, Omega_old))) /
        (sqrt(sum(vapply(Omega_old, function(M) sum(M^2), numeric(1)))) + 1e-12)
      dZ <- sqrt(sum(mapply(function(A, B) sum((A - B)^2), Z_new, Z_old))) /
        (sqrt(sum(vapply(Z_old, function(M) sum(M^2), numeric(1)))) + 1e-12)
      active_grad <- unlist(cg_refit_active_gradient_blocks(beta_new, Z_new, basis, support))
      max_grad <- if (length(active_grad) == 0L) 0 else max(abs(active_grad))
      max_z_grad <- cg_refit_z_gradient_stats(X_work, Z_new, beta_new,
                                           basis)$max_abs
      # ★ Yuan's LITERAL rule (Compo_glasso), see the z-align block above for the derivation:
      # his thresholds are mean(Omega_init^2)/O_ratio and mean(z_init^2)/z_ratio (the "end"
      # iterates are zero matrices when the thresholds are formed), i.e. mean-squared CHANGE vs
      # mean-squared MAGNITUDE -> in our normalized units dOm <= 1/sqrt(O_ratio) AND
      # dZ <= 1/sqrt(z_ratio). Both blocks' parameter change, no gradient test, exactly as
      # `while ((dOmega > O_thr || dZ > z_thr) && iter <= max_iter)`.
      state$cg_refit_streak <- if (isTRUE(yuan_conv)) {
        if (dOm <= 1 / sqrt(O_ratio) && dZ <= 1 / sqrt(z_ratio))
          state$cg_refit_streak + 1L else 0L
      } else if (dOm < conv_tol_Omega && dZ < conv_tol_Z &&
                 max_grad < conv_tol_grad && max_z_grad < conv_tol_z_grad) {
        state$cg_refit_streak + 1L
      } else 0L
      criterion <- cg_refit_joint_nll_average(X_work, Z_new, beta_new, basis)
      state$trace <- rbind(state$trace,
                           data.frame(phase = "refit", iter = iter,
                                      objective = criterion$total,
                                      dOmega = dOm, dZ = dZ,
                                      max_active_grad = max_grad,
                                      max_z_grad = max_z_grad,
                                      min_eig = min(vapply(Omega_new, .min_eig, numeric(1))),
                                      streak = state$cg_refit_streak))
      state$Z <- Z_new
      state$beta <- beta_new
      state$cg_refit_iter <- iter
      state$optimizer_stats[[length(state$optimizer_stats) + 1L]] <-
        list(update = "glasso_zero_mle", max_active_grad = max_grad)
      state$z_optimizer_stats[[length(state$z_optimizer_stats) + 1L]] <-
        list(phase = "refit", iter = iter,
             diagnostics = z_update$diagnostics)
      .cg_refit_atomic_save(state, checkpoint_file)
      .cg_refit_log(progress_log,
                 paste0("refit iter=%d obj=%.8f dOm=%.3e dZ=%.3e ",
                        "max|g_active|=%.3e max|g_Z|=%.3e ",
                        "minEig=%.4g streak=%d"),
                 iter, criterion$total, dOm, dZ, max_grad, max_z_grad,
                 min(vapply(Omega_new, .min_eig, numeric(1))),
                 state$cg_refit_streak)
      if (state$cg_refit_streak >= conv_sustain) {
        exit_reason <- "converged"
        break
      }
    }
  }

  final_Z <- state$Z
  final_beta <- state$beta
  post_valid <- !isTRUE(post_not_pd)
  if (post_valid) {
    final_Omega <- G_beta_Rcpp(final_beta, basis, m)
    fixed_check <- cg_refit_validate_fixed_coordinates(final_beta, beta_template, support, free_diag = free_diag)
    stopifnot(fixed_check$ok, .all_slices_pd(final_Omega))
    active_grad <- unlist(cg_refit_active_gradient_blocks(final_beta, final_Z, basis, support))
    max_grad <- if (length(active_grad) == 0L) 0 else max(abs(active_grad))
    max_z_grad <- cg_refit_z_gradient_stats(X_work, final_Z, final_beta, basis)$max_abs
    diag_score <- cg_refit_diagonal_score(final_beta, final_Z, basis)
    post_out <- list(Z = final_Z, beta = final_beta,
                     criterion = cg_refit_joint_nll_average(X_work, final_Z, final_beta, basis))
  } else {
    # constrained MLE does not exist on this support at this n -> NO valid de-biased post (post=NULL);
    # the penalized `pre` is retained and returned. The deployment gate excludes this rho's post.
    fixed_check <- list(ok = TRUE)
    max_grad <- NA_real_; max_z_grad <- NA_real_; diag_score <- NA_real_
    post_out <- NULL
  }

  list(version = "2026-07-16-refit-v2", config = config,
       support = support, support_signature = signature,
       df = config$df, beta_template = beta_template,
       pre = state$pre,
       post = post_out,
       post_valid = post_valid,
       trace = state$trace,
       optimizer_stats = state$optimizer_stats,
       z_optimizer_stats = state$z_optimizer_stats,
       z_align_converged = isTRUE(state$z_align_converged),
       z_align_iterations = state$z_iter,
       exit_reason = exit_reason,
       converged = post_valid && identical(exit_reason, "converged"),
       max_active_grad = max_grad,
       max_z_grad = max_z_grad,
       diagonal_score = diag_score,
       fixed_coordinate_check = fixed_check,
       checkpoint_file = checkpoint_file,
       progress_log = progress_log)
}
