# ============================================================================
# baselines/refit_cglasso.R — CGLasso (Tian 2023) relaxed refit.
#
# Uses CGLasso's OWN refit core (R/baselines/refit_cglasso_core.R : cg_refit_fixed_support),
# a cg_-prefixed INDEPENDENT copy of the refit machinery — TV and CGLasso NEVER share a refit
# file (editing one never touches the other). basis = I_m (identity) => beta_k = Omega_k, each
# time slice independent (NO cross-time borrowing — the static control).
#
# Method-faithful choices:
#   * RAW counts in the joint-LNM working likelihood (CompoGlasso NR is on raw counts;
#     the option-2 pseudocount only INITIALISES Z). X_work = raw X (NOT TV's pooled-adjusted).
#   * S_k = (1/n_k) sum (z-zbar)(z-zbar)' — the 1/n MLE convention (matches the likelihood
#     and TV-refit's S_Z_t), NOT cov()'s 1/(n-1).
#   * Support E_{k,rho} = per-slice nonzero off-diagonal of the pass-1 Omega_k (NO cross-time
#     union). The refit is the glasso-native CONSTRAINED MLE on that FIXED zero pattern
#     (glasso(S_k, rho=0, zero=non-edges); Dempster 1972 / ESL Alg 17.1): it RE-ESTIMATES the
#     diagonal (FREE) + the active off-diagonals and forces the non-edges to 0; NO rho penalty;
#     Z re-estimated. df/support ALWAYS from pass-1 EDGES (never re-thresholded on refit values;
#     the freed diagonal is not counted as a selected df -- it is always present).
#   * BIC granularity = SLICE-BIC  -2*ell + sum_k d_k*log(n_k)  (each Omega_k fit independently
#     on n_k samples). The old pooled  D*log(sum n_k)  over-penalises by log(m) -> kept as AUDIT only.
#     (The slice penalty already handles unequal n_k; only the neg2ll term assumes equal n —
#     stopifnot in cglasso_slice_ic — so real-data unequal-n needs the same relaxation as
#     refit_ic_unbalanced, i.e. drop the stopifnot. Sims here are equal-n.)
#
# The refit CORE (R/baselines/refit_cglasso_core.R : cg_refit_fixed_support, glasso-native) provides
# the guarantees: constrained-MLE Omega per slice (glasso `zero`+`rho=0`, PD by construction, free
# diagonal), exact-zero non-edges, exact-Z (NR/L-BFGS), atomic checkpoint, hard validity gates.
# ============================================================================

suppressPackageStartupMessages({ library(here); library(MASS); library(Matrix); library(splines); library(glasso) })
if (!exists("cg_refit_fixed_support", mode = "function")) source(here::here("R", "baselines", "refit_cglasso_core.R"))   # CG's OWN refit core (cg_-prefixed) — INDEPENDENT of TV's R/refit.R (editing one never touches the other)
if (!exists("fit_cglasso_static", mode = "function")) source(here::here("R", "baselines", "cglasso_static.R"))  # -> fit_cglasso_static + CompoGlasso (z_hat_offset)

`%||%` <- function(a, b) if (is.null(a)) b else a

.cg_is_pd  <- function(M) !is.null(tryCatch(chol((M + t(M)) / 2), error = function(e) NULL))
.cg_mineig <- function(M) min(eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values)

# PD repair (P1-3, 2026-07-22; margin refinement per codex 2026-07-22): if zeroing sub-tol
# off-diagonals lost PD, scale the ALLOWED off-diagonals toward the FIXED diagonal until the
# matrix is PD WITH A RELATIVE EIGENVALUE MARGIN (min_eig >= margin_rel * min(diag)) -- NOT the
# bare Cholesky-just-succeeds boundary (which leaves an ill-conditioned min_eig ~ 1e-15 warm
# start). No ridge, no diagonal change, no filling back disabled edges; the support is preserved
# for ANY s>0 (scaling never zeros an edge), and the refit re-estimates magnitudes from here.
# min_eig(D+s*OFF) is concave in s with min_eig(0)=min(diag) >= floor, so {s: min_eig>=floor} is
# an interval [0,s*] -> bisection is exact. Hard-errors only if NO positive s meets the margin.
.cg_pd_repair <- function(Om, margin_rel = 1e-4) {
  Om <- (Om + t(Om)) / 2
  D <- diag(diag(Om)); OFF <- Om - D
  floor_eig <- margin_rel * min(diag(D))
  ok <- function(s) .cg_mineig(D + s * OFF) >= floor_eig
  if (ok(1)) return(Om)                              # already PD with margin -> keep pass-1 Omega as-is
  lo <- 0; hi <- 1; s_ok <- NA_real_                 # s=0 (diagonal) meets the margin; s=1 does not (we're here)
  for (it in seq_len(50L)) {
    s <- (lo + hi) / 2
    if (ok(s)) { s_ok <- s; lo <- s } else hi <- s
  }
  if (is.na(s_ok) || s_ok <= .Machine$double.eps)
    stop(".cg_pd_repair: no positive off-diagonal scale keeps Omega PD with margin (pass-1 support cannot be preserved)")
  out <- D + s_ok * OFF
  attr(out, "pd_repair_scale") <- s_ok
  out
}

# SLICE-BIC / IC. nll_average = refit_joint_nll_average(...)$total (per-observation avg NLL);
# -2loglik = 2*N*nll_average (N = sum n_k, equal-n asserted as in the TV refit).
cglasso_slice_ic <- function(nll_average, n_per_slice, d_per_slice, P) {
  stopifnot(length(unique(n_per_slice)) == 1L)   # equal-n sims; real n_k-varying needs the weighted form
  N <- sum(n_per_slice); D <- sum(d_per_slice)
  neg2ll <- 2 * N * nll_average
  slice_pen <- sum(d_per_slice * log(n_per_slice))
  c(neg2loglik = neg2ll,
    AIC        = neg2ll + 2 * D,
    BIC_slice  = neg2ll + slice_pen,                 # <- pooled-support SLICE-BIC (AUDIT)
    eBIC_slice = neg2ll + slice_pen + D * log(P),
    BIC_pooled_AUDIT = neg2ll + D * log(N))          # <- old pooled, AUDIT ONLY (over-penalises by log m)
}

# ★ PER-SLICE IC (2026-07-22): the FAITHFUL per-slice CGLasso selection input. Each slice k is
# Yuan's STATIC CGLasso fit on n_k samples, so its BIC = -2*ell_k + d_k*log(n_k) with
# -2*ell_k = 2*by_slice[k] (the core's exact per-slice joint -2loglik contribution). This is
# what "each time point tuned/refit by its OWN independent BIC" needs; it is naturally
# unequal-n-safe (uses n_per_slice[k] directly). Deployment then picks, PER SLICE, the rho
# minimizing that slice's own BIC (done at analysis / in method_cglasso across the rho path).
.cg_per_slice_ic <- function(by_slice, n_per_slice, d_per_slice, P) {
  neg2ll_k <- 2 * by_slice
  data.frame(slice = seq_along(by_slice), n = n_per_slice, df = d_per_slice,
             neg2loglik = neg2ll_k,
             AIC  = neg2ll_k + 2 * d_per_slice,
             BIC  = neg2ll_k + d_per_slice * log(n_per_slice),
             eBIC = neg2ll_k + d_per_slice * log(n_per_slice) + d_per_slice * log(P))
}

# Refit ONE rho: Om_slices_rho = list (length m) of the pass-1 per-slice Omega at this rho.
cglasso_refit_one_rho <- function(X_raw, Om_slices_rho, checkpoint_file, progress_log,
                                  option = 2L, max_outer = 150L, inner_max = 10L,
                                  z_align_max = 160L, ...) {
  m <- length(X_raw); P <- ncol(X_raw[[1]]) - 1L
  basis <- diag(m)                                   # I_m => beta_k = Omega_k, slices independent
  repaired <- lapply(Om_slices_rho, .cg_pd_repair)   # fixed support+diagonal; PD-repair (continuous) if thresholding broke PD
  pd_repair_scale <- vapply(repaired, function(M) { s <- attr(M, "pd_repair_scale"); if (is.null(s)) 1.0 else s }, numeric(1))
  beta_start <- lapply(repaired, function(M) { attr(M, "pd_repair_scale") <- NULL; M })  # pristine matrices into the core
  Z_start <- lapply(X_raw, function(Xk)
    as.matrix(z_hat_offset(as.matrix(Xk), offset = P + 1L, option = option)))  # pseudocount INIT only (raw counts used in the likelihood)

  fit <- cg_refit_fixed_support(X_work = lapply(X_raw, as.matrix), Z_start = Z_start,   # CG's OWN core (not TV's refit_fixed_support)
                             beta_start = beta_start, basis = basis,
                             checkpoint_file = checkpoint_file, progress_log = progress_log,
                             max_outer = max_outer, inner_max = inner_max,
                             z_align_max = z_align_max, ...)

  n_per_slice <- vapply(X_raw, nrow, integer(1))
  d_per_slice <- vapply(fit$support, nrow, integer(1))
  pre_ic  <- cglasso_slice_ic(fit$pre$criterion$total,  n_per_slice, d_per_slice, P)   # pooled AUDIT
  pre_ic_slice  <- .cg_per_slice_ic(fit$pre$criterion$by_slice,  n_per_slice, d_per_slice, P)   # the SELECTION input
  # fit$post is NULL when the constrained MLE does not exist on this support (support not estimable
  # at this n -- e.g. a DENSE support with n<P); the penalized `pre` is still valid. Guard the post IC.
  post_ic       <- if (!is.null(fit$post)) cglasso_slice_ic(fit$post$criterion$total,  n_per_slice, d_per_slice, P) else NULL
  post_ic_slice <- if (!is.null(fit$post)) .cg_per_slice_ic(fit$post$criterion$by_slice, n_per_slice, d_per_slice, P) else NULL
  c(fit, list(pre_ic = pre_ic, post_ic = post_ic,
              pre_ic_slice = pre_ic_slice, post_ic_slice = post_ic_slice,
              pd_repair_scale = pd_repair_scale,
              d_per_slice = d_per_slice, n_per_slice = n_per_slice,
              df_total = sum(d_per_slice)))
}
