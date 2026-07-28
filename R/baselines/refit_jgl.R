# ============================================================================
# baselines/refit_jgl.R — JGL (Danaher et al. 2014) relaxed refit + Gaussian IC.
#
# JGL is a GAUSSIAN graphical method on the plug-in centered ALR Y (S_k = crossprod(Y_k)/n_k);
# it has NO LNM layer and NO latent Zhat. Its relaxed refit therefore stays in JGL's OWN
# Gaussian family (faithful — JGL never had a Zhat) and does NOT route through the LNM refit
# core used by TV/CGLasso. This file + method_jgl.R are JGL's OWN, INDEPENDENT set — editing
# them never touches any other method's files.
#
# ── THE REFIT — a STANDARD object, not an ad-hoc construction ─────────────────
# Two-pass relaxed graphical lasso / "refitted MLE in the selected model":
#   pass 1  fused JGL selects the per-slice support E_k (its nonzero off-diagonal pattern).
#   pass 2  drop BOTH penalties (lambda1 AND lambda2). Since lambda2=0 removes the fusion, the
#           JGL objective sum_k n_k[tr(S_k Theta_k) - logdet Theta_k] + pen DECOUPLES across
#           slices, and each slice becomes the UNPENALIZED CONSTRAINED GAUSSIAN MLE (Dempster
#           1972 covariance selection) on E_k:
#               Theta_k^refit = argmax_{Theta>0, Theta[off-E_k]=0} [ logdet Theta - tr(S_k Theta) ]
#           computed by  glasso(S_k, rho = 0, zero = off-support)  at GLASSO'S OWN DEFAULT
#           TOLERANCE (thr=1e-4, maxit=500 backstop — the SAME settings as the CGLasso refit;
#           see .jgl_slice_refit). rho = 0 is the TRUE unpenalized MLE (NOT rho=1e-6, which is
#           a tiny-L1 fit masquerading as unpenalized).
#           The diagonal is FREE: the constraint is only on the OFF-support inverse-covariance
#           entries — the standard refitted-MLE / relaxed-lasso phi=0 form (Meinshausen 2007).
#   Refit does NOT move the support (edges stay E_k) -> the lambda-sweep ROC is UNCHANGED; its
#   ONLY role is the de-biased DEPLOYED operating point (it moves the IC argmin off the too-sparse
#   penalized point) + magnitude accuracy.
#
# ── EXISTENCE of the given-graph MLE, and the RIDGE that guarantees it ───────
# The unpenalized given-graph MLE has a DOMAIN: it exists (finite, unique, PD) iff the graph-
# restricted S_k admits a POSITIVE-DEFINITE COMPLETION [Grone-Johnson-de Sa-Wolkowicz 1984,
# LAA 58:109-124; Dempster 1972; Uhler 2012 Ann.Statist. 40(1):238-261 Thm 2.1]; equivalently
# the EFFECTIVE SAMPLE SIZE must reach the graph's Maximum Likelihood Threshold [Buhl 1993
# Scand.J.Statist.; Uhler 2012; Gross-Sullivant 2018], with clique#(E_k) <= MLT(E_k) <=
# treewidth(E_k)+1. So it is the effective n vs the SELECTED GRAPH, not n_k vs P: a sparse support
# has a small MLT and the MLE exists even at n<P.
#   ★ THE EFFECTIVE SAMPLE SIZE HERE IS n_k - 1, NOT n_k. method_jgl.R centers each slice before
#     forming S_k = crossprod(Y_k)/n_k, so rank(S_k) <= n_k - 1 and the correct criterion is
#         n_k - 1  >=  MLT(E_k).
#     Stating it as "n_k >= MLT" (as an earlier version of this file did) is off by one IN THE
#     UNSAFE DIRECTION -- it declares a class of supports estimable that are not. Counterexample
#     (measured, P=15, n=12, support = complete graph on 12 of the 15 nodes, so MLT = 12): the
#     "n_k >= MLT" rule says 12 >= 12 => exists, but glasso returns max|Theta| = 3.4e5 and a
#     NON-PD matrix, so .jgl_slice_refit returns NULL. The true boundary is q = n_k - 1 = 11
#     (PD, max|Theta| = 120, moment-match residual 3.4e-04). The sibling CGLasso refit states the
#     rank bound correctly (refit_cglasso_core.R:146).
# Our per-slice S_k VIOLATES the criterion at the DENSE END of the lambda1 path: rank(S_k) <=
# n_k - 1 < P, and once the fused-selected support is too dense for that rank the MLE does not
# exist and the likelihood is literally UNBOUNDED [Uhler 2017 arXiv:1707.04345 Sec.3].
#
# TWO EARLIER CLAIMS IN THIS FILE WERE WRONG AND ARE CORRECTED HERE (measured 2026-07-25):
#   (i)  "maxit only bounds the wasted work". FALSE in cost. When the MLE does not exist glasso
#        never converges at ANY cap and each sweep is expensive on the near-singular problem.
#        MEASURED at P=15/n=12 (cell01, lambda1-index 30, support density 0.74): ONE slice's
#        glasso(rho=0, zero=) call did not return in 2 h 34 min, while the base JGL fit at the
#        same grid point took 4.9 s and the other 6 slices refitted in ~0 s. On the cluster the
#        n<P tasks advanced ZERO grid points in 3 h; the n>P cells finished all 300 in 12-16 min.
#        Worse, glasso's Fortran does not poll for R interrupts, so such a call cannot even be
#        interrupted -- only killed.
#   (ii) "a non-existent MLE diverges / goes indefinite -> Cholesky fails", i.e. PD detection is
#        sufficient. NOT SAFE. A PD-and-finite return can still be numerically degenerate, and
#        every PD-based guard is structurally blind to it (measured on the sibling CGLasso refit:
#        min_eig = +0.107 with cond = 7.2e23 and max diag(Omega) = 7.7e22). Such a Theta passes
#        is.finite + chol, contributes a huge +logdet, hence a hugely NEGATIVE -2loglik and AIC,
#        and would WIN the grid argmin -- a fabricated deployed point, not merely a slow one.
#
# ★★ THE FIX -- RIDGED SLICE COVARIANCE IN PASS 2 (`ridge_eps`, cfg$jgl_refit_ridge_eps,
# default 1e-3). OURS, MUST BE DISCLOSED -- it is NOT part of Danaher's JGL. (Same repair, same
# rationale and the same user decision as the sibling CGLasso refit, R/baselines/refit_cglasso_core.R;
# JGL keeps its OWN copy and its OWN cfg key -- the two refits never share a file or a parameter.)
# Solve pass 2 on  S_k + eps * mean(diag(S_k)) * I. A strictly positive-definite matrix IS its own
# positive-definite completion, so the given-graph MLE exists and is unique for EVERY graph and
# EVERY n -- one repair at the ROOT of the existence criterion rather than a patch per failure
# mode. Equivalently (exactly, since Omega_jj > 0 at any PD solution, so eps*sum|Omega_jj| =
# tr(eps*I*Omega)) it is a penalty eps on the DIAGONAL OF Omega ONLY: the selected OFF-diagonals
# stay COMPLETELY unpenalized, so the refit's de-biasing purpose -- the whole reason pass 2 exists
# -- is untouched. This is the ridge precision estimator of van Wieringen & Peeters (2016, CSDA
# 103:284-303; R pkg rags2ridges, ridgePchordal(S, lambda, zeros)) restricted to the selected
# support, and it is the SAME protection the penalized pass-1 already enjoys. glasso's own source
# documents the failure it repairs: "With rho=0, there may be convergence problems if the input
# matrix is not of full rank". eps is RELATIVE to mean(diag(S_k)), so it follows the data scale and
# is invariant to rescaling Y. Setting jgl_refit_ridge_eps = 0 restores the exact unridged MLE
# (and the stalls) -- kept as an explicit switch so the estimator is reproducible either way.
# HONEST SCOPE: every ingredient is published; the COMPOSITION ("ridge the slice covariance before
# the given-graph refit") is ours -> disclose it.
#
# CHOICE OF eps = 1e-3 -- CALIBRATED ON JGL'S OWN GEOMETRY, NOT COPIED. The sibling CGLasso refit
# uses 1e-2, but that value was measured on the LNM latent-Z covariance; JGL's plug-in ALR S_k is
# markedly MORE ridge-sensitive, so the number does not transfer and we re-measured it here
# (cell01 = P15/n12/low and cell03 = P15/n20/low, seed 1, lambda2 = 0; off-diagonal magnitude
# RETAINED relative to the exact eps = 0 MLE, at path points where that MLE exists):
#     density      0.08     0.24*    0.40     0.56     0.70          (* ~ the deployed operating point)
#     eps 1e-4   99.95%   99.89%   99.85%   99.19%   98.84%
#     eps 1e-3   99.49%   98.91%   98.48%   92.50%   89.65%
#     eps 1e-2   95.16%   90.29%   86.71%   57.56%   48.91%
# TERMINATION -- the actual defect -- is fixed at ANY eps > 0: every solve returned in < 2 s at
# every eps tested, including the densest grid point (density 0.96). eps only controls whether the
# solve at the EXTREME dense end is numerically usable: at density 0.96 the per-slice MLE was
# recovered 0/7 at eps=1e-4, 6/7 at 1e-3, 7/7 at 1e-2 (cond ~5e2). We take 1e-3 because
#   (a) the refit exists ONLY to undo the ~30% L1 shrinkage, so paying ~10% of the off-diagonal
#       magnitude (eps=1e-2, at the deployed density) would forfeit a third of its purpose --
#       ~1% (eps=1e-3) is the same order as the 1.7% the CGLasso refit accepts;
#   (b) the grid points 1e-3 cannot refit sit at the extreme dense end, which the df-penalised IC
#       never deploys; they are FLAGGED refit_exists = FALSE, never silently substituted.
# Using a different eps from CGLasso is deliberate: each method's ridge is calibrated on its own
# covariance geometry and disclosed, rather than inheriting a constant measured on another one.
#
# THE RIDGE IS CONFINED TO THE SOLVER. It is applied to a LOCAL copy of S_k inside
# .jgl_slice_refit and NOWHERE else: jgl_gaussian_neg2ll below scores the resulting Theta against
# the ORIGINAL, UNRIDGED S_list. The ridge regularises the ESTIMATOR (so that it exists); the
# reported likelihood remains the honest Gaussian likelihood of the data given that Theta. (Same
# convention as the CGLasso refit: refit_cglasso_core.R:192 ridges locally, :294 scores on S_Z_t(Z).)
#
# ── LIKELIHOOD / IC ──────────────────────────────────────────────────────────
# Gaussian deviance  -2ell = sum_k n_k [ tr(S_k Theta_k) - logdet Theta_k ], logdet via CHOLESKY
# (returns NA if Theta_k is not PD -> that candidate's IC is NA -> it is NEVER selected). This
# PD check applies to BASE and REFIT alike: a base JGL fit that hit maxiter can return a non-PD
# auxiliary Z (the ADMM's sparse variable is only PD at convergence), so its likelihood must be
# PD-checked too. SLICE-IC (per-slice independent, matching the fully-relaxed per-slice refit):
#     AIC  = -2ell + 2 D
#     BIC  = -2ell + sum_k d_k * log(n_k)
#     eBIC = BIC + 4 * gamma * D * log(P)              (Foygel-Drton 2010; gamma default 0.25)
# D = sum_k d_k off-diagonal edges (the constant P diagonal params drop from the grid argmin;
# same convention as the CGLasso slice-IC). df d_per_slice is ALWAYS the pass-1 support (never
# re-thresholded on refit magnitudes).
#
# References: Dempster 1972 (Biometrics 28:157) covariance selection; Speed-Kiiveri 1986 (Ann.
# Statist. 14:138) IPS; Buhl 1993 MLE existence; Meinshausen 2007 (CSDA 52:374) relaxed lasso;
# Danaher-Wang-Witten 2014 (JRSS-B 76:373) JGL.
# ============================================================================
suppressPackageStartupMessages({ library(glasso) })

# Confusion (FPR/TPR/F1) pooled over eval_slices; off-diagonal edge = Omega != 0. (JGL's own copy.)
.jgl_confusion <- function(Om_list, true_list, eval_slices) {
  P <- nrow(Om_list[[1]]); ut <- upper.tri(matrix(0, P, P))
  pred <- unlist(lapply(eval_slices, function(k) { d <- Om_list[[k]]; diag(d) <- 0; d[ut] != 0 }))
  tru  <- unlist(lapply(eval_slices, function(k) { d <- true_list[[k]]; diag(d) <- 0; d[ut] != 0 }))
  TP <- sum(pred & tru); FP <- sum(pred & !tru); FN <- sum(!pred & tru); TN <- sum(!pred & !tru)
  TPR <- TP / max(TP + FN, 1); FPR <- FP / max(FP + TN, 1); PR <- TP / max(TP + FP, 1)
  F1  <- if (TPR + PR > 0) 2 * TPR * PR / (TPR + PR) else 0
  c(FPR = FPR, TPR = TPR, precision = PR, F1 = F1)      # precision reported too (shared `deployed` contract)
}

# PD-safe log-determinant: log det(Theta) on Theta>0 via Cholesky, else NA.
# (Fix: the old determinant()$modulus = log|det| ignored the sign, so an indefinite Theta with
# an even number of negative eigenvalues [det>0] passed with a finite, even spuriously optimal,
# "likelihood". chol() succeeds iff Theta is symmetric PD -> the correct precision-likelihood check.)
.jgl_logdet_pd <- function(Th) {
  Th <- (Th + t(Th)) / 2
  R <- tryCatch(chol(Th), error = function(e) NULL)
  if (is.null(R)) return(NA_real_)
  2 * sum(log(diag(R)))
}

# Unpenalized constrained Gaussian MLE (covariance selection) on ONE slice's support: the nonzero
# off-diagonal pattern of theta_pen. glasso(S, rho=0, zero=off-support) = the TRUE unpenalized MLE
# (free diagonal, i.e. the fitted covariance matches S on the diagonal AND on the edges).
#
# ★ SOLVER SETTINGS = the SAME as the CGLasso refit (R/baselines/refit_cglasso_core.R:130), which
# was calibrated first: thr = 1e-4 = GLASSO'S OWN DEFAULT, maxit = 500 as a mere backstop. Do NOT
# hand-tighten a reproduced method's tolerances: a tighter thr (this file previously had 1e-6) is
# NOT reachable on the near-singular problems that arise at n<P dense supports, while the Omega
# difference at the default thr is only ~1e-4 on entries of order 0.1-1 (immaterial). Use the
# package's own convergence control. With the ridge below the solve converges far under 500 sweeps,
# so maxit is a pure backstop; WITHOUT the ridge no cap suffices (see the header: the likelihood is
# unbounded where the MLE does not exist).
#
# EXISTENCE: guaranteed by the ridge (a strictly PD S has a PD completion for EVERY graph), so on
# the default eps > 0 the MLE exists at every grid point. The finite + symmetric-PD checks below are
# RETAINED as a belt-and-braces contract (they also cover eps = 0, and any pathological glasso
# return): on failure we return NULL and the caller FLAGS that grid point's refit invalid, keeps the
# penalized base, and excludes it from the deployment candidates — never a silent substitute.
# NB the PD check alone is NOT a sufficient degeneracy test (header item (ii)); the ridge, not the
# check, is what makes the estimator well posed.
.jgl_slice_refit <- function(S, theta_pen, maxit = 500L, thr = 1e-4, ridge_eps = 1e-3) {
  ut_zero <- which(upper.tri(theta_pen) & theta_pen == 0, arr.ind = TRUE)   # off-support upper-tri -> constrain to 0
  Ssym <- (S + t(S)) / 2
  # ★ OURS, disclosed: ridge the slice covariance for the SOLVE ONLY (local copy; the likelihood in
  #   jgl_gaussian_neg2ll is scored against the caller's ORIGINAL S). eps is relative to the data
  #   scale via mean(diag(S)). eps = 0 -> exact unridged MLE. See the header block.
  if (ridge_eps > 0) Ssym <- Ssym + ridge_eps * mean(diag(Ssym)) * diag(nrow(Ssym))
  fit <- tryCatch(withCallingHandlers({
    if (nrow(ut_zero) == 0L)
      glasso::glasso(Ssym, rho = 0, penalize.diagonal = FALSE, thr = thr, maxit = maxit)   # full support
    else
      glasso::glasso(Ssym, rho = 0, zero = ut_zero, penalize.diagonal = FALSE, thr = thr, maxit = maxit)
  }, warning = function(w) if (grepl("rho=0|full rank", conditionMessage(w))) invokeRestart("muffleWarning")),
  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  Th <- fit$wi
  if (!all(is.finite(Th))) return(NULL)                                   # diverged (MLE does not exist)
  Th <- (Th + t(Th)) / 2                                                  # ONE canonical symmetric Omega (glasso BCD can be mildly asymmetric)
  if (is.null(tryCatch(chol(Th), error = function(e) NULL))) return(NULL) # not PD -> MLE does not exist on this support
  Th
}

# Per-slice constrained MLE on each slice's pass-1 support. Returns a list of per-slice refit Theta
# (NULL where the MLE does not exist). The grid point's refit "exists" iff ALL slices exist (needed
# for the joint per-grid-point IC — JGL selects ONE (l1,l2) for the whole slice set).
jgl_refit_supports <- function(S_list, theta_pen_list, maxit = 500L, thr = 1e-4, ridge_eps = 1e-3) {
  lapply(seq_along(S_list), function(k)
    .jgl_slice_refit(S_list[[k]], theta_pen_list[[k]], maxit, thr, ridge_eps))
}

# Gaussian -2loglik on the plug-in ALR: -2ell = sum_k n_k (tr(S_k Theta_k) - logdet Theta_k),
# logdet via the PD-safe Cholesky. Returns NA if ANY slice is NULL (refit non-existent) or not PD.
jgl_gaussian_neg2ll <- function(S_list, theta_list, n_per_slice) {
  vals <- vapply(seq_along(theta_list), function(k) {
    Th <- theta_list[[k]]; if (is.null(Th)) return(NA_real_)
    ld <- .jgl_logdet_pd(Th); if (is.na(ld)) return(NA_real_)
    n_per_slice[k] * (sum(S_list[[k]] * Th) - ld)          # sum(S*Th) = tr(S Theta) for symmetric
  }, numeric(1))
  if (anyNA(vals)) NA_real_ else sum(vals)
}

# Gaussian SLICE-IC. eBIC = BIC + 4*gamma*D*log(P) (Foygel-Drton 2010; gamma explicit, default 0.25
# -> identical to the CGLasso slice-eBIC's +D*log(P)). df d_per_slice = pass-1 support (fixed).
jgl_slice_ic <- function(neg2ll, n_per_slice, d_per_slice, P, ebic_gamma = 0.25) {
  D <- sum(d_per_slice); slice_pen <- sum(d_per_slice * log(n_per_slice))
  c(neg2loglik = neg2ll, AIC = neg2ll + 2 * D,
    BIC_slice = neg2ll + slice_pen, eBIC_slice = neg2ll + slice_pen + 4 * ebic_gamma * D * log(P))
}
