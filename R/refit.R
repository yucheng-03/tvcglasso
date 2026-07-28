# ============================================================================
# refit.R — TV's fixed-support, FREE-diagonal relaxed (de-biased) refit.
#
# ★ TV-OWNED (2026-07-22 engineering-discipline separation): this file is OUR method's
#   (TV) refit core ONLY. CGLasso keeps an INDEPENDENT copy in
#   R/baselines/refit_cglasso_core.R (functions cg_-prefixed) — do NOT source THIS from
#   the CGLasso path, and do NOT couple the two via a shared cfg key (TV -> tv_refit_max_outer,
#   CGLasso -> cglasso_refit_max_outer). A numerical fix here must be MIRRORED into the CG copy.
#
# The publication estimator is TV(base fit) + refit. Pass 1 (the penalized
# tv_warm_path fit) SELECTS the off-diagonal support; pass 2 (this file) re-estimates,
# UNPENALIZED, the active off-diagonal magnitudes AND the DIAGONAL AND the latent Z on
# that fixed support. The optimization variables are GENUINELY reduced coordinates: the
# ADAM state holds the selected upper-triangular beta entries PLUS the P diagonal entries
# per basis block (the ACTIVE set = refit_active_support); inactive OFF-diagonal entries
# are reconstructed as exact zeros (this is the "reduced parameterization", NOT gradient-
# masking). The diagonal is FREE and unpenalized — the relaxed-lasso / de-biasing / Dempster
# convention (freed 2026-07-24; it was previously held fixed at the base-fit value). df/IC
# still count the OFF-diagonal support only (the always-free diagonal is a constant offset).
#
# Model-selection likelihood = the JOINT LNM density evaluated at the INFERRED
# latent Zhat (multinomial + Gaussian-graphical layers), via refit_joint_nll_average
# -> refit_information_criteria. It is NOT the Z_0 Gaussian-only loglik.
#
# Copied verbatim from code/refit_fixed_support_2026-07-16.R; the ONLY changes here
# are (a) the self-source guard points at the packaged engine R/tvcglasso.R, and
# (b) refit_ic_unbalanced is added next to refit_information_criteria for real data
# with UNEQUAL per-slice n (the balanced version asserts equal n).
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(MASS)
  library(Matrix)
  library(splines)
})

if (!exists("main_function_final", mode = "function") ||
    !exists("G_beta_Rcpp", mode = "function")) {
  source(here::here("R", "tvcglasso.R"))
}

.refit_atomic_save <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- sprintf("%s.tmp.%d", path, Sys.getpid())
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("atomic save failed: ", path)
  }
  invisible(path)
}

.refit_log <- function(path, fmt, ...) {
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
refit_prepare_counts <- function(X) {
  stopifnot(is.list(X), length(X) > 0L)
  X
}

refit_support_from_beta <- function(beta, zero_tol = 0) {
  lapply(beta, function(B) {
    idx <- which(upper.tri(B) & abs(B) > zero_tol, arr.ind = TRUE)
    if (length(idx) == 0L) matrix(integer(0), nrow = 0L, ncol = 2L,
                                  dimnames = list(NULL, c("row", "col"))) else idx
  })
}

refit_support_df <- function(support) sum(vapply(support, nrow, integer(1)))

# ★ FREE-DIAGONAL refit (2026-07-24): the diagonal is ALWAYS a free, unpenalized parameter
# (relaxed-lasso / de-biasing convention: refit re-estimates ALL retained entries — the
# selected off-diagonals AND the always-present diagonal — as the unpenalized MLE on the
# selected zero pattern; Meinshausen 2007 / Dempster 1972 / glasso). ACTIVE coordinates per
# basis block = the P diagonal entries (i,i) PLUS the selected off-diagonal support. df/IC
# still count the OFF-DIAGONAL support only (the always-free diagonal is a constant offset
# that does not change model selection).
refit_active_support <- function(support, P) {
  diag_idx <- cbind(seq_len(P), seq_len(P)); colnames(diag_idx) <- c("row", "col")
  lapply(support, function(idx) if (nrow(idx) == 0L) diag_idx else rbind(diag_idx, idx))
}

# Reduced-coordinate gradient extracted from the per-block Omega-gradient matrix Gfull =
# sum_k w_k B_h(t_k) (S_k - Omega_k^{-1}). An OFF-diagonal coordinate beta_ij uses Gfull[i,j]
# (its two symmetric appearances in Omega cancel the objective's 1/2); a DIAGONAL coordinate
# beta_ii uses 0.5*Gfull[i,i] (it enters Omega only once). Matches the base engine's dgrad_lt_h.
.refit_grad_at <- function(Gfull, idx) {
  if (nrow(idx) == 0L) return(numeric(0))
  g <- Gfull[idx]; on_diag <- idx[, 1] == idx[, 2]; g[on_diag] <- 0.5 * g[on_diag]; g
}

refit_support_signature <- function(support) {
  paste(vapply(seq_along(support), function(h) {
    idx <- support[[h]]
    if (nrow(idx) == 0L) sprintf("h%d:", h) else
      sprintf("h%d:%s", h, paste(sprintf("%d-%d", idx[, 1], idx[, 2]), collapse = ","))
  }, character(1)), collapse = "|")
}

refit_pack_blocks <- function(beta, support) {
  Map(function(B, idx) {
    if (nrow(idx) == 0L) numeric(0) else as.numeric(B[idx])
  }, beta, support)
}

refit_rebuild_beta <- function(theta_blocks, beta_template, support) {
  P <- nrow(beta_template[[1]])
  stopifnot(length(theta_blocks) == length(beta_template),
            length(support) == length(beta_template))
  lapply(seq_along(beta_template), function(h) {
    B <- matrix(0, P, P)
    diag(B) <- diag(beta_template[[h]])
    idx <- support[[h]]
    if (nrow(idx) > 0L) {
      stopifnot(length(theta_blocks[[h]]) == nrow(idx))
      B[idx] <- theta_blocks[[h]]
      B[cbind(idx[, 2], idx[, 1])] <- theta_blocks[[h]]
    }
    B
  })
}

# The refit must keep every INACTIVE entry (off the diagonal AND off the selected off-diagonal
# support) EXACTLY zero, and Omega symmetric. `active` = the diagonal + off-diagonal support (from
# refit_active_support). The diagonal is now a FREE, re-estimated parameter, so there is no longer
# a "diagonal held at the template" check (that was the fixed-diagonal refit; freed 2026-07-24).
refit_validate_fixed_coordinates <- function(beta, active, tol = 0) {
  P <- nrow(beta[[1]])
  symmetric_ok <- all(vapply(beta, function(B) max(abs(B - t(B))) <= tol, logical(1)))
  inactive_ok <- all(vapply(seq_along(beta), function(h) {
    allowed <- matrix(FALSE, P, P)
    idx <- active[[h]]
    if (nrow(idx) > 0L) {
      allowed[idx] <- TRUE
      allowed[cbind(idx[, 2], idx[, 1])] <- TRUE
    }
    all(beta[[h]][!allowed] == 0)
  }, logical(1)))
  list(ok = symmetric_ok && inactive_ok,
       symmetric = symmetric_ok,
       inactive_exact_zero = inactive_ok)
}

# ★ P0-6/codex: slice k weighted by w_slice[k]=n_k/N (=1/m for equal n -> bit-identical).
refit_beta_objective <- function(beta, basis, m, Z = NULL, S_list = NULL,
                                 Omega_list = NULL, w_slice = NULL) {
  if (is.null(w_slice)) w_slice <- if (!is.null(Z)) { nk <- vapply(Z, nrow, integer(1)); nk/sum(nk) } else rep(1/m, m)
  if (is.null(S_list)) S_list <- S_Z_t(Z)
  if (is.null(Omega_list)) Omega_list <- G_beta_Rcpp(beta, basis, m)
  logdets <- vapply(Omega_list, .logdet_chol, numeric(1))
  if (anyNA(logdets)) return(Inf)
  -0.5 * sum(w_slice * logdets) +
    0.5 * sum(vapply(seq_len(m), function(k) {
      w_slice[k] * sum(S_list[[k]] * Omega_list[[k]])
    }, numeric(1)))
}

refit_active_gradient_blocks <- function(beta, Z, basis, support) {
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
    .refit_grad_at(G, idx)   # off-diag: G[i,j]; diagonal: 0.5*G[i,i] (no-op when idx has no diagonal)
  })
}

refit_diagonal_score <- function(beta, Z, basis) {
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
refit_joint_nll_average <- function(X_work, Z, beta, basis) {
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
refit_z_slice_value_gradient <- function(par, X, Omega) {
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

refit_renew_z_exact <- function(X_work, Z_start, Omega_list,
                                maxit = 1000L, pgtol = 1e-8,
                                factr = 1e3) {
  stopifnot(length(X_work) == length(Z_start),
            length(Z_start) == length(Omega_list))
  fits <- lapply(seq_along(X_work), function(k) {
    cache <- new.env(parent = emptyenv())
    evaluate <- function(par) {
      if (!is.null(cache$par) && identical(par, cache$par)) return(cache$out)
      cache$par <- par
      cache$out <- refit_z_slice_value_gradient(par, X_work[[k]],
                                                 Omega_list[[k]])
      cache$out
    }
    optim(par = as.numeric(Z_start[[k]]),
          fn = function(par) evaluate(par)$value,
          gr = function(par) evaluate(par)$gradient,
          method = "L-BFGS-B",
          control = list(maxit = maxit, pgtol = pgtol, factr = factr))
  })
  Z <- Map(function(f, old) matrix(f$par, nrow(old), ncol(old)),
           fits, Z_start)
  diagnostics <- do.call(rbind, lapply(seq_along(fits), function(k) {
    f <- fits[[k]]
    g <- refit_z_slice_value_gradient(f$par, X_work[[k]], Omega_list[[k]])
    data.frame(slice = k, convergence = f$convergence,
               value = f$value, max_abs_gradient = max(abs(g$gradient)),
               iterations_function = unname(f$counts["function"]),
               iterations_gradient = unname(f$counts["gradient"]),
               message = if (is.null(f$message) || !length(f$message)) ""
                 else as.character(f$message)[1],
               stringsAsFactors = FALSE)
  }))
  ## ★ D-1 (2026-07-27): slice k enters the JOINT gradient with weight n_k/N,
  ## not 1/m. refit_z_slice_value_gradient already divides by n_k (see its
  ## `/ n` above), so the joint gradient w.r.t. Z_k is (n_k/N) * g_k. Dividing
  ## by length(X_work) = m instead is EXACT when every slice has the same n
  ## (n_k/N = 1/m) -- which is why every equal-n simulation is bit-unchanged --
  ## but wrong for unequal n. On the zebrafish days (n = 7..23) it makes the
  ## convergence test up to 40% too strict on the smallest day and 21% too
  ## loose on the largest. This is a STOPPING-RULE correction: it changes WHEN
  ## the refit declares convergence, not the estimate it converges to.
  nk_ <- vapply(Z, nrow, integer(1))
  list(Z = Z, diagnostics = diagnostics,
       usable = all(vapply(fits, function(f) {
         all(is.finite(f$par)) && is.finite(f$value)
       }, logical(1))) && all(is.finite(diagnostics$max_abs_gradient)),
       all_converged = all(diagnostics$convergence == 0L),
       max_abs_gradient = max(diagnostics$max_abs_gradient * nk_ / sum(nk_)))
}

refit_z_gradient_stats <- function(X_work, Z, beta, basis) {
  Omega <- G_beta_Rcpp(beta, basis, length(Z))
  N_ <- sum(vapply(Z, nrow, integer(1)))          # D-1: n_k/N, not 1/m (see above)
  by_slice <- lapply(seq_along(Z), function(k) {
    g <- refit_z_slice_value_gradient(as.numeric(Z[[k]]), X_work[[k]],
                                      Omega[[k]])$gradient
    matrix(g, nrow(Z[[k]]), ncol(Z[[k]])) * (nrow(Z[[k]]) / N_)
  })
  flat <- unlist(by_slice)
  list(by_slice = by_slice, max_abs = max(abs(flat)),
       l2 = sqrt(sum(flat^2)), rms = sqrt(mean(flat^2)))
}

# ★ P0-6 (2026-07-22): `nll_average` from refit_joint_nll_average is now the n_k-weighted
# per-observation joint neg-loglik, so 2*N*nll_average IS the true joint -2loglik for ANY
# n_k (equal OR unequal). This ONE function is correct for both sims (equal n) and real
# Zebrafish (unequal n) -- no assertion, no separate unbalanced variant. df = # nonzero
# off-diagonal beta coefficients; N = sum(n_per_slice).
refit_information_criteria <- function(nll_average, n_per_slice, df, P) {
  N <- sum(n_per_slice)
  neg2loglik <- 2 * N * nll_average
  c(neg2loglik = neg2loglik,
    AIC = neg2loglik + 2 * df,
    BIC = neg2loglik + log(N) * df,
    eBIC = neg2loglik + (log(N) + log(P)) * df)
}
# Back-compat alias (the unequal-n case is now handled by refit_information_criteria itself).
refit_ic_unbalanced <- refit_information_criteria

refit_add_direction <- function(beta, support, direction_blocks, scale) {
  out <- beta
  for (h in seq_along(out)) {
    idx <- support[[h]]
    if (nrow(idx) == 0L) next
    vals <- out[[h]][idx] + scale * direction_blocks[[h]]
    out[[h]][idx] <- vals
    out[[h]][cbind(idx[, 2], idx[, 1])] <- vals
  }
  out
}

refit_fd_check <- function(beta, Z, basis, support, n_direction = 3L,
                           eps = 1e-6, seed = 716L) {
  if (refit_support_df(support) == 0L) {
    return(data.frame(direction = integer(0), analytic = numeric(0),
                      numeric = numeric(0), rel_error = numeric(0)))
  }
  set.seed(seed)
  grad <- refit_active_gradient_blocks(beta, Z, basis, support)
  rows <- vector("list", n_direction)
  for (r in seq_len(n_direction)) {
    direction <- lapply(support, function(idx) rnorm(nrow(idx)))
    dn <- sqrt(sum(unlist(direction)^2))
    direction <- lapply(direction, function(x) x / dn)
    analytic <- sum(unlist(Map(`*`, grad, direction)))
    step <- eps
    repeat {
      bp <- refit_add_direction(beta, support, direction, step)
      bm <- refit_add_direction(beta, support, direction, -step)
      if (.all_slices_pd(G_beta_Rcpp(bp, basis, length(Z))) &&
          .all_slices_pd(G_beta_Rcpp(bm, basis, length(Z)))) break
      step <- step / 2
      if (step < 1e-12) stop("FD step could not remain in the PD cone")
    }
    fp <- refit_beta_objective(bp, basis, length(Z), Z = Z)
    fm <- refit_beta_objective(bm, basis, length(Z), Z = Z)
    numerical <- (fp - fm) / (2 * step)
    rel <- abs(analytic - numerical) /
      max(1e-10, abs(analytic) + abs(numerical))
    rows[[r]] <- data.frame(direction = r, analytic = analytic,
                            numeric = numerical, rel_error = rel,
                            step = step)
  }
  do.call(rbind, rows)
}

# One beta update with true reduced-coordinate ADAM.  The Armijo slope is the
# exact reduced-coordinate directional derivative g^T d.  If ADAM momentum is
# not a descent direction, the method falls back to -g.
refit_renew_beta_reduced <- function(Z, beta, beta_template, support, basis,
                                     initial_learning_rate = 0.01,
                                     max_iterations = 10L,
                                     grad_tol = 1e-5,
                                     pd_max_backtrack = 25L,
                                     pd_c1 = 1e-4,
                                     pd_nonmono_K = 5L) {
  m <- length(Z)
  S <- S_Z_t(Z)
  nk <- vapply(Z, nrow, integer(1)); w_slice <- nk / sum(nk)   # P0-6: n_k/N (=1/m for equal n)
  theta <- refit_pack_blocks(beta, support)
  objective <- function(bt, Om = NULL) {
    refit_beta_objective(bt, basis, m, S_list = S, Omega_list = Om, w_slice = w_slice)
  }
  stats <- list(adam_accept = 0L, gradient_accept = 0L,
                block_break = 0L, inner_iterations = 0L)

  for (h in seq_along(beta)) {
    if (length(theta[[h]]) == 0L) next
    m_adam <- numeric(length(theta[[h]]))
    v_adam <- numeric(length(theta[[h]]))
    F_hist <- rep(objective(beta), pd_nonmono_K)

    for (i in seq_len(max_iterations)) {
      stats$inner_iterations <- stats$inner_iterations + 1L
      Om <- G_beta_Rcpp(beta, basis, m)
      Gfull <- matrix(0, nrow(beta[[1]]), ncol(beta[[1]]))
      for (k in seq_len(m)) {
        Gfull <- Gfull + w_slice[k] * basis[k, h] * (S[[k]] - .inv_pd(Om[[k]]))   # P0-6: n_k/N
      }
      g <- .refit_grad_at(Gfull, support[[h]])   # `support` here = ACTIVE (diag+off-diag); diag gets 0.5
      if (max(abs(g)) < grad_tol) break

      m_adam <- 0.9 * m_adam + 0.1 * g
      v_adam <- 0.999 * v_adam + 0.001 * (g^2)
      m_hat <- m_adam / (1 - 0.9^i)
      v_hat <- v_adam / (1 - 0.999^i)
      d_adam <- -m_hat / (sqrt(v_hat) + 1e-8)
      slope_adam <- sum(g * d_adam)
      accepted <- FALSE
      F_ref <- max(F_hist)

      if (is.finite(slope_adam) && slope_adam < 0) {
        for (bt in 0:pd_max_backtrack) {
          alpha <- initial_learning_rate * 0.5^bt
          theta_cand <- theta
          theta_cand[[h]] <- theta[[h]] + alpha * d_adam
          cand <- refit_rebuild_beta(theta_cand, beta_template, support)
          Om_cand <- G_beta_Rcpp(cand, basis, m)
          if (.all_slices_pd(Om_cand)) {
            Fnew <- objective(cand, Om_cand)
            if (is.finite(Fnew) && Fnew <= F_ref + pd_c1 * alpha * slope_adam) {
              theta <- theta_cand
              beta <- cand
              accepted <- TRUE
              stats$adam_accept <- stats$adam_accept + 1L
              break
            }
          }
        }
      }

      if (!accepted) {
        d_grad <- -g
        slope_grad <- -sum(g^2)
        F_old <- objective(beta)
        for (bt in 0:pd_max_backtrack) {
          alpha <- initial_learning_rate * 0.5^bt
          theta_cand <- theta
          theta_cand[[h]] <- theta[[h]] + alpha * d_grad
          cand <- refit_rebuild_beta(theta_cand, beta_template, support)
          Om_cand <- G_beta_Rcpp(cand, basis, m)
          if (.all_slices_pd(Om_cand)) {
            Fnew <- objective(cand, Om_cand)
            if (is.finite(Fnew) && Fnew <= F_old + pd_c1 * alpha * slope_grad) {
              theta <- theta_cand
              beta <- cand
              accepted <- TRUE
              stats$gradient_accept <- stats$gradient_accept + 1L
              break
            }
          }
        }
      }

      if (!accepted) {
        stats$block_break <- stats$block_break + 1L
        break
      }
      F_hist <- c(F_hist[-1], objective(beta))
    }
  }
  list(beta = beta, stats = stats)
}

refit_fixed_support <- function(X_work, Z_start, beta_start, basis,
                                checkpoint_file, progress_log,
                                max_outer = 150L, inner_max = 10L,
                                z_align_max = 160L,
                                z_optim_maxit = 1000L,
                                z_optim_pgtol = 1e-8,
                                z_optim_factr = 1e3,
                                initial_learning_rate = 0.01,
                                conv_tol_Omega = 1e-4,
                                conv_tol_Z = 5e-5,
                                conv_tol_z_grad = 1e-4,
                                conv_tol_grad = 1e-4,
                                conv_sustain = 3L,
                                inner_grad_tol = 1e-5,
                                pd_max_backtrack = 25L,
                                pd_c1 = 1e-4,
                                pd_nonmono_K = 5L) {
  m <- length(X_work)
  P <- nrow(beta_start[[1]])
  support <- refit_support_from_beta(beta_start)          # selected OFF-diagonal edges (for df/IC)
  active  <- refit_active_support(support, P)              # + the always-free diagonal (re-estimated)
  signature <- refit_support_signature(support)
  beta_template <- beta_start
  theta0 <- refit_pack_blocks(beta_start, active)
  beta_start <- refit_rebuild_beta(theta0, beta_template, active)
  initial_check <- refit_validate_fixed_coordinates(beta_start, active)
  stopifnot(initial_check$ok,
            .all_slices_pd(G_beta_Rcpp(beta_start, basis, m)))

  config <- list(P = P, m = m, J_n = length(beta_start),
                 df = refit_support_df(support), signature = signature,
                 max_outer = max_outer, inner_max = inner_max,
                 z_align_max = z_align_max,
                 z_method = "joint_slice_lbfgsb",
                 z_optim_maxit = z_optim_maxit,
                 z_optim_pgtol = z_optim_pgtol,
                 z_optim_factr = z_optim_factr,
                 initial_learning_rate = initial_learning_rate,
                 conv_tol_Omega = conv_tol_Omega,
                 conv_tol_Z = conv_tol_Z,
                 conv_tol_z_grad = conv_tol_z_grad,
                 conv_tol_grad = conv_tol_grad,
                 conv_sustain = conv_sustain,
                 inner_grad_tol = inner_grad_tol,
                 pd_max_backtrack = pd_max_backtrack,
                 pd_c1 = pd_c1, pd_nonmono_K = pd_nonmono_K)

  if (file.exists(checkpoint_file)) {
    state <- readRDS(checkpoint_file)
    if (!identical(state$config, config)) {
      extension_fields <- setdiff(names(config), "max_outer")
      max_outer_extension <- setequal(names(state$config), names(config)) &&
        all(vapply(extension_fields, function(nm) {
          identical(state$config[[nm]], config[[nm]])
        }, logical(1))) &&
        isTRUE(config$max_outer > state$config$max_outer) &&
        isTRUE(config$max_outer > state$refit_iter)
      legacy_fields <- c(
        "P", "m", "J_n", "df", "signature", "max_outer", "inner_max",
        "initial_learning_rate", "conv_tol_Omega", "conv_tol_Z",
        "conv_tol_grad", "conv_sustain", "inner_grad_tol",
        "pd_max_backtrack", "pd_c1", "pd_nonmono_K")
      legacy_compatible <- identical(state$phase, "z_align") &&
        identical(state$refit_iter, 0L) &&
        !"z_method" %in% names(state$config) &&
        all(vapply(legacy_fields, function(nm) {
          identical(state$config[[nm]], config[[nm]])
        }, logical(1))) && z_align_max > state$z_iter
      if (max_outer_extension) {
        old_max_outer <- state$config$max_outer
        state$config <- config
        .refit_atomic_save(state, checkpoint_file)
        .refit_log(progress_log,
                   "extend checkpoint max_outer from %d to %d at refit_iter=%d",
                   old_max_outer, config$max_outer, state$refit_iter)
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
        .refit_atomic_save(state, checkpoint_file)
        .refit_log(progress_log,
                   paste0("migrate legacy Z checkpoint at iter=%d; retain Z and ",
                          "continue with exact joint-slice Z block"), old_iter)
      }
    }
    .refit_log(progress_log, "resume phase=%s z_iter=%d refit_iter=%d df=%d",
               state$phase, state$z_iter, state$refit_iter, config$df)
  } else {
    state <- list(version = "2026-07-16-refit-v2", config = config,
                  phase = "z_align", Z = Z_start, beta = beta_start,
                  z_iter = 0L, refit_iter = 0L,
                  z_streak = 0L, z_align_converged = FALSE,
                  refit_streak = 0L, trace = data.frame(),
                  optimizer_stats = list(), z_optimizer_stats = list())
    .refit_atomic_save(state, checkpoint_file)
    .refit_log(progress_log, "start z-alignment df=%d", config$df)
  }

  Omega_fixed <- G_beta_Rcpp(beta_start, basis, m)
  if (identical(state$phase, "z_align")) {
    start <- state$z_iter + 1L
    if (start <= z_align_max) {
      for (iter in seq.int(start, z_align_max)) {
        z_update <- refit_renew_z_exact(
          X_work, state$Z, Omega_fixed,
          maxit = z_optim_maxit, pgtol = z_optim_pgtol,
          factr = z_optim_factr)
        if (!z_update$usable) {
          state$last_z_optimizer_failure <- z_update$diagnostics
          .refit_atomic_save(state, checkpoint_file)
          stop("exact Z optimizer returned a non-finite alignment state: ",
               checkpoint_file)
        }
        if (!z_update$all_converged) {
          .refit_log(progress_log,
                     "z-align iter=%d retains finite optimizer intermediate; nonzero slices=%s",
                     iter,
                     paste(z_update$diagnostics$slice[
                       z_update$diagnostics$convergence != 0L], collapse = ","))
        }
        Z_new <- z_update$Z
        dZ <- sqrt(sum(mapply(function(A, B) sum((A - B)^2), Z_new, state$Z))) /
          (sqrt(sum(vapply(state$Z, function(M) sum(M^2), numeric(1)))) + 1e-12)
        max_z_grad <- z_update$max_abs_gradient
        state$z_streak <- if (dZ < conv_tol_Z &&
                                max_z_grad < conv_tol_z_grad) {
          state$z_streak + 1L
        } else 0L
        nll <- refit_joint_nll_average(X_work, Z_new, beta_start, basis)$total
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
        .refit_atomic_save(state, checkpoint_file)
        .refit_log(progress_log,
                   paste0("z-align iter=%d obj=%.8f dZ=%.3e ",
                          "max|g_Z|=%.3e streak=%d"),
                   iter, nll, dZ, max_z_grad, state$z_streak)
        if (state$z_streak >= conv_sustain) break
      }
    }
    state$z_align_converged <- state$z_streak >= conv_sustain
    .refit_atomic_save(state, checkpoint_file)
    if (!isTRUE(state$z_align_converged)) {
      .refit_log(progress_log,
                 paste0("STOP z-alignment did not converge by iter=%d; ",
                        "checkpoint retained and refit phase not entered"),
                 state$z_iter)
      stop("Z-only alignment failed its convergence gate after ",
           state$z_iter, " iterations: ", checkpoint_file)
    }
    state$pre <- list(Z = state$Z, beta = beta_start,
                      criterion = refit_joint_nll_average(X_work, state$Z,
                                                           beta_start, basis))
    state$phase <- "refit"
    state$beta <- beta_start
    state$refit_iter <- 0L
    state$refit_streak <- 0L
    .refit_atomic_save(state, checkpoint_file)
    .refit_log(progress_log, "z-alignment complete iter=%d pre_obj=%.8f",
               state$z_iter, state$pre$criterion$total)
  }

  exit_reason <- if (state$refit_streak >= conv_sustain) {
    "converged"
  } else "max_outer"
  start <- state$refit_iter + 1L
  if (!identical(exit_reason, "converged") && start <= max_outer) {
    for (iter in seq.int(start, max_outer)) {
      beta_old <- state$beta
      Z_old <- state$Z
      Omega_old <- G_beta_Rcpp(beta_old, basis, m)
      z_update <- refit_renew_z_exact(
        X_work, Z_old, Omega_old,
        maxit = z_optim_maxit, pgtol = z_optim_pgtol,
        factr = z_optim_factr)
      if (!z_update$usable) {
        state$last_z_optimizer_failure <- z_update$diagnostics
        .refit_atomic_save(state, checkpoint_file)
        stop("exact Z optimizer returned a non-finite refit state: ",
             checkpoint_file)
      }
      if (!z_update$all_converged) {
        .refit_log(progress_log,
                   "refit iter=%d retains finite optimizer intermediate; nonzero slices=%s",
                   iter,
                   paste(z_update$diagnostics$slice[
                     z_update$diagnostics$convergence != 0L], collapse = ","))
      }
      Z_new <- z_update$Z
      upd <- refit_renew_beta_reduced(
        Z = Z_new, beta = beta_old, beta_template = beta_template,
        support = active, basis = basis,   # ACTIVE = diag + off-diag support (diagonal re-estimated)
        initial_learning_rate = initial_learning_rate,
        max_iterations = inner_max, grad_tol = inner_grad_tol,
        pd_max_backtrack = pd_max_backtrack, pd_c1 = pd_c1,
        pd_nonmono_K = pd_nonmono_K)
      beta_new <- upd$beta
      Omega_new <- G_beta_Rcpp(beta_new, basis, m)
      stopifnot(.all_slices_pd(Omega_new))
      fixed_check <- refit_validate_fixed_coordinates(beta_new, active)
      stopifnot(fixed_check$ok)

      dOm <- sqrt(sum(mapply(function(A, B) sum((A - B)^2), Omega_new, Omega_old))) /
        (sqrt(sum(vapply(Omega_old, function(M) sum(M^2), numeric(1)))) + 1e-12)
      dZ <- sqrt(sum(mapply(function(A, B) sum((A - B)^2), Z_new, Z_old))) /
        (sqrt(sum(vapply(Z_old, function(M) sum(M^2), numeric(1)))) + 1e-12)
      active_grad <- unlist(refit_active_gradient_blocks(beta_new, Z_new, basis, active))
      max_grad <- if (length(active_grad) == 0L) 0 else max(abs(active_grad))
      max_z_grad <- refit_z_gradient_stats(X_work, Z_new, beta_new,
                                           basis)$max_abs
      state$refit_streak <- if (dOm < conv_tol_Omega && dZ < conv_tol_Z &&
                                 max_grad < conv_tol_grad &&
                                 max_z_grad < conv_tol_z_grad) {
        state$refit_streak + 1L
      } else 0L
      criterion <- refit_joint_nll_average(X_work, Z_new, beta_new, basis)
      state$trace <- rbind(state$trace,
                           data.frame(phase = "refit", iter = iter,
                                      objective = criterion$total,
                                      dOmega = dOm, dZ = dZ,
                                      max_active_grad = max_grad,
                                      max_z_grad = max_z_grad,
                                      min_eig = min(vapply(Omega_new, .min_eig, numeric(1))),
                                      streak = state$refit_streak))
      state$Z <- Z_new
      state$beta <- beta_new
      state$refit_iter <- iter
      state$optimizer_stats[[length(state$optimizer_stats) + 1L]] <- upd$stats
      state$z_optimizer_stats[[length(state$z_optimizer_stats) + 1L]] <-
        list(phase = "refit", iter = iter,
             diagnostics = z_update$diagnostics)
      .refit_atomic_save(state, checkpoint_file)
      .refit_log(progress_log,
                 paste0("refit iter=%d obj=%.8f dOm=%.3e dZ=%.3e ",
                        "max|g_active|=%.3e max|g_Z|=%.3e ",
                        "minEig=%.4g streak=%d"),
                 iter, criterion$total, dOm, dZ, max_grad, max_z_grad,
                 min(vapply(Omega_new, .min_eig, numeric(1))),
                 state$refit_streak)
      if (state$refit_streak >= conv_sustain) {
        exit_reason <- "converged"
        break
      }
    }
  }

  final_beta <- state$beta
  final_Z <- state$Z
  final_Omega <- G_beta_Rcpp(final_beta, basis, m)
  fixed_check <- refit_validate_fixed_coordinates(final_beta, active)
  stopifnot(fixed_check$ok, .all_slices_pd(final_Omega))
  active_grad <- unlist(refit_active_gradient_blocks(final_beta, final_Z, basis, active))
  max_grad <- if (length(active_grad) == 0L) 0 else max(abs(active_grad))
  max_z_grad <- refit_z_gradient_stats(X_work, final_Z, final_beta,
                                       basis)$max_abs
  diag_score <- refit_diagonal_score(final_beta, final_Z, basis)

  list(version = "2026-07-16-refit-v2", config = config,
       support = support, support_signature = signature,
       df = config$df, beta_template = beta_template,
       pre = state$pre,
       post = list(Z = final_Z, beta = final_beta,
                   criterion = refit_joint_nll_average(X_work, final_Z,
                                                        final_beta, basis)),
       trace = state$trace,
       optimizer_stats = state$optimizer_stats,
       z_optimizer_stats = state$z_optimizer_stats,
       z_align_converged = isTRUE(state$z_align_converged),
       z_align_iterations = state$z_iter,
       exit_reason = exit_reason,
       converged = identical(exit_reason, "converged"),
       max_active_grad = max_grad,
       max_z_grad = max_z_grad,
       diagonal_score = diag_score,
       fixed_coordinate_check = fixed_check,
       checkpoint_file = checkpoint_file,
       progress_log = progress_log)
}
