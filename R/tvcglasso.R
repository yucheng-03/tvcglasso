# ============================================================================
# tvcglasso.R — TVCGLasso engine: sparse TIME-VARYING precision-matrix estimation
# for longitudinal compositional counts.
#
# Model. Counts X_i(t_k) are logistic-normal-multinomial: a latent Gaussian
# Z_i(t_k) ~ N(mu, Omega(t_k)) is mapped to multinomial probabilities by an
# additive-log-ratio (ALR) link. Each entry of the precision matrix is a smooth
# function of time, Omega_ij(t) = sum_h beta_h[i,j] * B_h(t), expanded in a
# B-spline basis B_h; an adaptive sliding-window group-lasso on the beta's makes
# each edge sparse and time-local. The fit alternates estimating Z and beta.
#
# PUBLIC ENTRY POINTS
#   main_function_final(X, lambda, x_sequence, ...)  — ONE fit at a single lambda.
#       Block-coordinate descent: alternate an NR update of the latent Z (given
#       Omega) with an ADAM update of beta (given Z) under a Cholesky positive-
#       definiteness (PD) guard, until Omega and Z stop changing; then a final
#       entrywise threshold selects the sparse edge set. Returns beta (+ the
#       pre-threshold beta_presel), Z, the adaptive weights, and PD / convergence
#       diagnostics. The reported Omega_hat(t_k) is PD on every slice or it errors.
#   tv_warm_path(X, lambda_grid, x_sequence, ...)    — a WHOLE lambda path, fit
#       from HIGH (sparse) to LOW (dense) lambda, each lambda warm-started from
#       the previous (sparser) solution. This continuation is the fix for the
#       n<P "upward stall". Returns an ascending-lambda list of the fits above.
#
# INTERNALS (in the order they appear below)
#   Z-update    : obj, NR, renew_z            — per-sample Newton-Raphson for the
#                                               LNM latent Gaussian Z.
#   beta-update : Ffinal_beta, dgrad_lt_h,
#                 renew_beta_final            — PD-guarded ADAM on the smooth
#                                               penalized objective (no in-loop
#                                               thresholding).
#   penalty     : tv_pen_value, tv_pen_grad,
#                 tv_compute_W                — adaptive OVERLAPPING-group lasso on
#                                               the B-spline coefficients, made
#                                               differentiable by Nesterov smoothing.
#   init/helpers: generate_beta, S_Z_t,
#                 calculate_percentile        — project a per-slice glasso fit onto
#                                               the spline basis, per-slice second
#                                               moment S(Z), default threshold.
#   PD utils    : .sym, .pd_chol_ok, .min_eig,
#                 .all_slices_pd, .logdet_chol
#
# CONVENTIONS. ALR zero-handling = proportional pseudocount pooled across slices,
# used for the Z INIT ONLY (the multinomial likelihood uses the RAW counts). The
# diagonal of Omega is freely estimated (free_diag = TRUE, unpenalized). Slices
# are weighted by n_k / N in the likelihood/gradient (= 1/m for equal n). Requires
# Matrix, MASS, splines, Rcpp, glasso, here; compiles tvcglasso.cpp (G_beta_Rcpp,
# which assembles Omega(t_k) = sum_h B_h(t_k) beta_h) on load. Run from the repo
# root so here::here() resolves.
# ============================================================================

suppressPackageStartupMessages({ library(Matrix); library(MASS); library(splines); library(Rcpp); library(glasso) })
Rcpp::sourceCpp(here::here("R", "tvcglasso.cpp"))   # -> G_beta_Rcpp (beta -> Omega(t) assembler)

# Diagnostic counter env (reset per main_function_final call): Path-B line-search
# outcomes (ADAM accepts / negative-gradient fallbacks / block breaks).
.pd_diag_env <- new.env(parent = emptyenv())
.pd_diag_reset <- function() { .pd_diag_env$adam_accept <- 0L; .pd_diag_env$fallback_accept <- 0L; .pd_diag_env$block_break <- 0L }
.pd_diag_reset()


# ============================================================================
# main_function_final — the single-lambda fit (block-coordinate descent).
#
# Two-stage adaptive: a preliminary estimate sets the frozen adaptive weights W,
# then the main solve runs BCD to convergence. Path-B keeps every iterate PD; a
# final entrywise threshold (per B-spline block, PD-guarded) selects the edges.
# ============================================================================
main_function_final <- function(X, lambda, x_sequence_realdata, N_n = 1, q = 2,
                                r = 1.0, mu = 1e-3, free_diag = TRUE,   # (C) ADOPTED 2026-07-22: diagonal freely estimated (unpenalized), the glasso/CGLasso/Xue-Shu-Qu standard; data-validated (free_diag experiment). FALSE = legacy frozen diagonal.
                                Max_iterations = 150, initial_learning_rate = 0.01,
                                conv_tol_Omega = 1e-4, conv_tol_Z = 5e-5, conv_sustain = 3L, inner_max = 10L,
                                pilot_max = 6L, pilot_inner = 6L,
                                pd_init_eps = 1e-3, sel_threshold = NULL, sel_type = "hard",
                                weight_mode = "glasso",
                                warm_beta = NULL,   # NULL => glasso init (bit-identical to 0608); else BCD warm-starts from this beta
                                two_stage = TRUE, record_pd_trace = TRUE, conv_trace = FALSE, verbose = FALSE) {
  # The engine is fully deterministic given its inputs (glasso / lm / ADAM / NR use no RNG), so
  # no seeding is needed. A legacy set.seed(7) was removed here (2026-07-24): it was a no-op that
  # seeded a stream the engine never draws from, and only risked clobbering the caller's RNG.
  m <- length(X); P <- dim(X[[1]])[2] - 1; J_n <- N_n + q + 1; offset_x <- P + 1
  n_per_slice <- vapply(X, nrow, integer(1)); w_slice <- n_per_slice / sum(n_per_slice)   # P0-6: n_k/N (=1/m for equal n)
  midpoints <- seq(x_sequence_realdata[1], x_sequence_realdata[m], length.out = N_n+2)[-c(1, N_n+2)]
  basis <- bs(x_sequence_realdata, degree=q, knots=midpoints, Boundary.knots=c(x_sequence_realdata[1], x_sequence_realdata[m]), intercept=TRUE)
  .pd_diag_reset()

  # ---- init: per-slice glasso -> project onto the B-spline basis (generate_beta) ----
  Z_0 <- vector("list", m); Omega_0 <- vector("list", m)
  # ALR zero-handling: Yuan z_hat_offset OPTION 2 (proportional pseudocount), POOLED across slices.
  # p.hat = mean relative abundance per taxon, pooled over ALL slices, so a globally-retained taxon
  # is never all-zero (p.hat>0 => no log(0)); offset mass = P+1, so total added mass per sample
  # equals the old flat-+1 (sum(p.hat)=1), just redistributed proportionally to abundance.
  # ★ P0-5 (2026-07-22): the pseudocount is used ONLY to build the initial ALR Z_0 (to avoid log(0));
  # X is NOT overwritten, so the multinomial LIKELIHOOD (renew_z/NR below) uses the RAW observed counts
  # -- the model as written. (Previously X was overwritten with augmented counts, leaking the pseudocount
  # into the likelihood; matters at low depth. Matches Yuan's CGLasso: pseudocount for init only.)
  Xpool <- do.call(rbind, X); off2 <- colMeans(Xpool / rowSums(Xpool)) * offset_x
  X_aug <- lapply(X, function(Xi) t(t(Xi) + off2))          # augmented — for the Z_0 ALR init ONLY
  for (i in 1:m) {
    Z_0[[i]] <- log(X_aug[[i]][, -(P+1)] / X_aug[[i]][, P+1])
    Omega_0[[i]] <- glasso(s = cov(Z_0[[i]]), rho = lambda)$wi
  }
  beta_0 <- generate_beta(Omega_0, P, m, q, J_n, basis)     # X stays RAW from here on (likelihood uses raw counts)

  # ---- Phase-I PD init repair (partition-of-unity ridge shift) ----
  init_repair_c <- 0
  Om0 <- G_beta_Rcpp(beta_0, basis, m); min_e <- min(vapply(Om0, .min_eig, numeric(1)))
  if (min_e < pd_init_eps) { init_repair_c <- pd_init_eps - min_e
    for (h in 1:J_n) diag(beta_0[[h]]) <- diag(beta_0[[h]]) + init_repair_c }

  # ---- warm-start BCD from a supplied beta (PD-repaired); NULL => glasso init (bit-identical) ----
  bcd_start <- NULL
  if (!is.null(warm_beta)) {
    bcd_start <- warm_beta
    Omw <- G_beta_Rcpp(bcd_start, basis, m); mew <- min(vapply(Omw, .min_eig, numeric(1)))
    if (mew < pd_init_eps) { sh <- pd_init_eps - mew; for (h in 1:J_n) diag(bcd_start[[h]]) <- diag(bcd_start[[h]]) + sh }
  }

  pd_trace <- list(); ctrace <- list(); iteration <- 0
  rec <- function(Om, tag, stage) { if (!record_pd_trace) return(invisible(NULL))
    pd_trace[[length(pd_trace)+1L]] <<- data.frame(stage=stage, iter=iteration, tag=tag,
      k=seq_along(Om), min_eig=vapply(Om,.min_eig,numeric(1)), pd=vapply(Om,.pd_chol_ok,logical(1))) }

  # ---- BCD solve for a given (frozen) W: alternate renew_z and renew_beta_final ----
  solve_bcd <- function(W, stage, max_outer, inner, reweight = FALSE, beta_init = NULL) {
    bi <- if (is.null(beta_init)) beta_0 else beta_init          # warm start (NULL => glasso init, bit-identical)
    Om_init <- if (is.null(beta_init)) Omega_0 else G_beta_Rcpp(bi, basis, m)
    Z_old <- Z_0; Omega_old <- Om_init; beta_old <- bi
    rec(G_beta_Rcpp(bi, basis, m), "init", stage)
    conv_streak <- 0L; iteration <<- 0
    repeat {
      if (reweight) W <- tv_compute_W(beta_old, q, N_n, r = r)   # iterative reweighting (LLA): W from current beta
      Z_new <- renew_z(X, Z_old, Omega_old, m)
      beta_new <- renew_beta_final(Z_new, beta_old, basis, J_n, m, P, q, lambda, W, N_n, mu,
                                   initial_learning_rate = initial_learning_rate,
                                   max_iterations = inner, free_diag = free_diag, w_slice = w_slice)
      Omega_new <- G_beta_Rcpp(beta_new, basis, m)
      rec(Omega_new, "post_beta", stage)
      # PARAMETER-CHANGE convergence (Yuan/CompoGlasso & CGLasso family): stop when the relative
      # change in Omega(t) (all slices) AND in Z are both small, sustained `conv_sustain` iters.
      # Replaces the old beta-only objective relative-decrease (which omitted the Z-likelihood term
      # and was unreachable at high lambda). Tols sit BELOW the per-config relative change observed
      # at edge-set stabilization, so the criterion never stops before the reported edge-set is
      # stable; Max_iterations is the backstop (set above the empirical max stabilization iter ~118).
      dOm <- sqrt(sum(mapply(function(A,B) sum((A-B)^2), Omega_new, Omega_old))) /
             (sqrt(sum(vapply(Omega_old, function(M) sum(M^2), numeric(1)))) + 1e-12)
      dZ  <- sqrt(sum(mapply(function(A,B) sum((A-B)^2), Z_new, Z_old))) /
             (sqrt(sum(vapply(Z_old, function(M) sum(M^2), numeric(1)))) + 1e-12)
      conv_streak <- if (iteration > 0 && dOm < conv_tol_Omega && dZ < conv_tol_Z) conv_streak + 1L else 0L
      if (conv_trace) {  # diagnostic only (default off); does not affect the result
        ctrace[[length(ctrace)+1L]] <<- list(stage=stage, iter=iteration, dOm_rel=dOm, dZ_rel=dZ,
          conv_streak=conv_streak, min_eig=min(vapply(Omega_new,.min_eig,numeric(1))), beta=beta_new)
      }
      if (verbose) cat(sprintf("  [%s] iter %d dOm %.2e dZ %.2e streak %d min_eig %.4f\n",
        stage, iteration, dOm, dZ, conv_streak, min(vapply(Omega_new,.min_eig,numeric(1)))))
      hit_conv <- conv_streak >= conv_sustain; hit_cap <- iteration >= max_outer
      if (hit_conv || hit_cap) break
      Z_old <- Z_new; Omega_old <- Omega_new; beta_old <- beta_new; iteration <<- iteration + 1
    }
    # convergence CONTRACT: converged = broke on the param-change criterion (NOT the max_outer cap),
    # so the analysis/validated-DONE layer can flag a path that hit the iteration cap. n_outer/final_d*
    # give the stationarity at exit. (Base-fit status gap flagged by the 2026-07-22 re-audit.)
    list(Z = Z_new, beta = beta_new,
         converged = isTRUE(hit_conv), exit_reason = if (isTRUE(hit_conv)) "param_change" else "max_outer",
         n_outer = iteration, final_dOm = dOm, final_dZ = dZ)
  }

  # ---- adaptive-weight scheme (weight_mode): "glasso" (default) | "pilot" | "iterative" ----
  # The frozen plug-in weights W = (||gamma_hat||+eps)^(-r) need a preliminary estimate gamma_hat:
  #   glasso    : from beta_0 (the glasso-init already used as the BCD start) -> consistent, fast (no pilot),
  #               clean monotone lambda-swept ROC. DEFAULT (chosen 2026-06-08 over pilot/iterative by ROC shape).
  #   pilot     : from a cheap 6-iter pilot fit (the old 0605 default). Non-converged -> unstable W -> zigzag ROC (OC-2).
  #   iterative : reweighted-L1 / LLA -- W recomputed from the current beta each outer iter (non-convex log penalty).
  W_one <- lapply(1:(N_n+1), function(.) matrix(1, P, P))
  if (weight_mode == "glasso") {
    W <- tv_compute_W(beta_0, q, N_n, r = r)
    fit <- solve_bcd(W, "adaptive", max_outer = Max_iterations, inner = inner_max, beta_init = bcd_start)
  } else if (weight_mode == "iterative") {
    fit <- solve_bcd(W_one, "adaptive", max_outer = Max_iterations, inner = inner_max, reweight = TRUE, beta_init = bcd_start)
    W <- tv_compute_W(fit$beta, q, N_n, r = r)   # final W, for reporting only
  } else {  # "pilot"
    pilot <- solve_bcd(W_one, "pilot", max_outer = pilot_max, inner = pilot_inner, beta_init = bcd_start)
    W <- tv_compute_W(pilot$beta, q, N_n, r = r)
    fit <- solve_bcd(W, "adaptive", max_outer = Max_iterations, inner = inner_max, beta_init = bcd_start)
  }

  # ---- final PD-guarded ENTRYWISE-THRESHOLD selection (A4-i); sel_type = "hard" (default) | "soft" ----
  # "hard" = keep-or-kill (Lan/XSQ convention; survivors keep FULL magnitude -> more faithful Omega_hat);
  # "soft" = sign(x)*max(|x|-thr,0) (shrinks survivors by thr). Both zero the SAME entries (|x|<=thr), so
  # the selected edge-set and the lambda-swept ROC are IDENTICAL; they differ only in surviving magnitude.
  thr_fun <- if (identical(sel_type, "hard"))
    function(x, thr) { y <- ifelse(abs(x) > thr, x, 0); diag(y) <- diag(x); y }
  else
    function(x, thr) { y <- sign(x) * pmax(abs(x) - thr, 0); diag(y) <- diag(x); y }
  if (is.null(sel_threshold)) { sel_threshold <- calculate_percentile(fit$beta, 0.2); sel_threshold <- max(0.01, ifelse(is.na(sel_threshold),0.01,sel_threshold)) }
  beta_sel <- fit$beta
  # ★ P1-5 (2026-07-22): the final threshold runs PER B-spline block; a block that would break PD
  # gets its threshold HALVED (or reverted to un-thresholded). So the EFFECTIVE threshold can differ
  # per block and depends on block order. RECORD it: threshold_per_block[h] = the tau actually applied
  # to block h (NA if the block was reverted un-thresholded); threshold_reverted[h] = TRUE if reverted.
  # (Empirically hard==soft and tau ~= its 0.01 floor -> revert is rare; kept for the pre-publication
  # health check, per the user's "record now, revisit after results".)
  threshold_per_block <- rep(sel_threshold, J_n); threshold_reverted <- logical(J_n)
  for (h in 1:J_n) {
    th <- sel_threshold; cand <- beta_sel; cand[[h]] <- thr_fun(beta_sel[[h]], th); tries <- 0L; rev_h <- FALSE
    # PD-1 fix: halve until PD, OR give up by reverting block h to its un-thresholded value (keeps
    # the assembly PD, since beta_sel was PD on entry to this block-iteration).
    while (!.all_slices_pd(G_beta_Rcpp(cand, basis, m))) {
      th <- th/2; tries <- tries + 1L
      if (th < 1e-12 || tries >= 80L) { cand[[h]] <- beta_sel[[h]]; rev_h <- TRUE; break }
      cand[[h]] <- thr_fun(beta_sel[[h]], th)
    }
    threshold_per_block[h] <- if (rev_h) NA_real_ else th; threshold_reverted[h] <- rev_h
    beta_sel <- cand
  }
  # PD-1 safety net: the REPORTED estimator must be PD on every slice.
  stopifnot("reported Omega_hat(t_k) not PD on all slices" = .all_slices_pd(G_beta_Rcpp(beta_sel, basis, m)))

  pd_trace_df <- if (record_pd_trace && length(pd_trace)>0) do.call(rbind, pd_trace) else NULL
  list(Z = fit$Z, beta = beta_sel, beta_presel = fit$beta, W = W,
       pd_trace = pd_trace_df, init_repair_c = init_repair_c, sel_threshold = sel_threshold,
       threshold_per_block = threshold_per_block, threshold_reverted = threshold_reverted,   # P1-5
       r = r, mu = mu, free_diag = free_diag, two_stage = two_stage, weight_mode = weight_mode,
       ls_adam = .pd_diag_env$adam_accept, ls_fallback = .pd_diag_env$fallback_accept, ls_break = .pd_diag_env$block_break,
       converged = isTRUE(fit$converged), exit_reason = fit$exit_reason,   # base convergence contract (re-audit gap)
       n_outer = fit$n_outer, final_dOm = fit$final_dOm, final_dZ = fit$final_dZ,
       conv_trace_data = if (conv_trace) ctrace else NULL)
}


# ============================================================================
# tv_warm_path — high->low-lambda continuation (homotopy / warm-start) driver.
#
# Fits a lambda grid from HIGH (sparse) to LOW (dense), each lambda warm-started
# from the previous (sparser) solution. This is the validated fix for the n<P
# "upward stall" (d18-d24). Returns an ASCENDING-lambda list of main_function_final
# fits. Optional per-lambda checkpointing (ckpt_file) makes it resumable.
#   init_mode "diag"   : first (highest) lambda starts from a sparse DIAGONAL beta
#                        (glasso(rho=lam_hi) per slice, off-diagonal zeroed).  [DEFAULT]
#   init_mode "pooled" : first lambda starts from a TIME-POOLED glasso (stack all
#                        slices, n_eff=n*m) held CONSTANT in time (static-CGLasso start).
# `...` is forwarded to main_function_final (sel_type, weight_mode, q, N_n, ...).
# ============================================================================
tv_warm_path <- function(X, lambda_grid, x_sequence_realdata, q = 2, N_n = 1,
                         init_mode = "diag", init_rho = NULL, ckpt_file = NULL, ...) {
  m <- length(X); P <- dim(X[[1]])[2] - 1; J_n <- N_n + q + 1; offset_x <- P + 1
  midpoints <- seq(x_sequence_realdata[1], x_sequence_realdata[m], length.out = N_n+2)[-c(1, N_n+2)]
  basis <- bs(x_sequence_realdata, degree=q, knots=midpoints,
              Boundary.knots=c(x_sequence_realdata[1], x_sequence_realdata[m]), intercept=TRUE)
  # ALR zero-handling: option2 (proportional), pooled p.hat across slices (matches main_function_final)
  off2 <- colMeans(do.call(rbind, X) / rowSums(do.call(rbind, X))) * offset_x
  Z0 <- lapply(X, function(Xi){ Xo <- t(t(Xi) + off2); log(Xo[, -(P+1)] / Xo[, P+1]) })
  if (init_mode == "pooled") {                                   # time-pooled (centered) glasso -> constant-in-time warm start
    Zc <- lapply(Z0, function(Z) sweep(Z, 2, colMeans(Z)))
    Om_pool <- glasso(cov(do.call(rbind, Zc)), rho = if (is.null(init_rho)) max(lambda_grid) else init_rho)$wi
    warm <- lapply(1:J_n, function(.) Om_pool)                   # Omega(t)=Om_pool for all t (partition-of-unity)
  } else {                                                       # "diag": sparse diagonal start
    Om_hi <- lapply(Z0, function(Z) glasso(cov(Z), rho = max(lambda_grid))$wi)
    beta0_hi <- generate_beta(Om_hi, P, m, q, J_n, basis)
    warm <- lapply(beta0_hi, function(M) diag(diag(M)))
  }
  lams <- sort(lambda_grid, decreasing = TRUE)
  out <- list(); start_i <- 1L
  if (!is.null(ckpt_file) && file.exists(ckpt_file)) {           # RESUME from per-lambda checkpoint
    cp <- readRDS(ckpt_file); out <- cp$out; warm <- cp$warm; start_i <- cp$next_i
  }
  idx <- if (start_i <= length(lams)) start_i:length(lams) else integer(0)
  for (i in idx) {                                              # ckpt_file=NULL => start_i=1 => bit-identical loop
    lam <- lams[i]
    f <- main_function_final(X = X, lambda = lam, x_sequence_realdata = x_sequence_realdata,
                             q = q, N_n = N_n, warm_beta = warm, ...)
    warm <- f$beta_presel                                        # warm-start next (lower) lambda
    out[[sprintf("%.6f", lam)]] <- f
    if (!is.null(ckpt_file)) saveRDS(list(out = out, warm = warm, next_i = i + 1L), ckpt_file)  # CHECKPOINT per lambda
  }
  out[order(as.numeric(names(out)))]   # ascending-lambda
}


# ============================================================================
# Z-UPDATE — Newton-Raphson for the LNM latent Gaussian (renew_z -> NR -> obj).
#
# Given Omega, infer the latent Gaussian Z from the multinomial counts. Each
# sample's Z-row is found by Newton-Raphson on the (negative) log-joint density
# obj(); NR() loops over samples in one slice; renew_z() maps NR over all slices.
# ============================================================================

# Per-sample objective: Gaussian-graphical quadratic in (z - mu) minus the
# multinomial log-likelihood. Used only as the line-search accept/shrink test.
obj <- function(x, z, Omega, K) {
  M <- sum(x)
  mu = mean(z)
  f = 1 / 2 * t(z - mu) %*% Omega %*% (z - mu) - (t(x) %*% z - M * log(as.numeric(t(rep(1, K)) %*% exp(z) + 1)))
  return(as.numeric(f))
}

# Newton-Raphson update of Z for one slice (all n samples), with a Levenberg-
# Marquardt damped Hessian and an Armijo backtracking line search.
NR <- function(x, z.0, Omega.0, alpha_0 = 1, delta = 5, epsilon = 0.001, threshold = 0.0001, lambda = 1e-4) {
  # Initialization
  n = dim(z.0)[1]
  K = dim(z.0)[2]
  M <- as.numeric(apply(x, 1, sum))
  mu.0 = apply(z.0, 2, mean)
  z.1 <- matrix(0, n, K)

  for (j in 1:n) {
    z.iter <- 0
    alpha <- alpha_0
    time_rep <- 0
    # Loop to update z for the j-th sample
    while (mean((z.0[j,] - z.1[j,])^2) > threshold) {
      if(time_rep > 100) break

      if (z.iter != 0 &&
          (obj(x[j, 1:K], z.1[j,], Omega.0, K) <= (obj(x[j, 1:K], z.0[j,], Omega.0, K) + epsilon * alpha * h_0))) {
        z.0[j,] <- z.1[j,]
      }

      z.iter <- z.iter + 1

      # Gradient of the objective function (equation 7) with respect to z
      dipi <- M[j] * exp(z.0[j,]) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1) -
        x[j, 1:K] + as.vector(Omega.0 %*% (z.0[j,] - mu.0))

      # Hessian of the objective function with respect to z
      tripi <- M[j] * diag(exp(z.0[j,])) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1) -
        M[j] * (exp(z.0[j,])) %*% t(exp(z.0[j,])) / (as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1))^2 +
        Omega.0

      # Levenberg-Marquardt modification: add damping term to Hessian
      tripi.reg <- tripi + lambda * diag(K)

      # Newton-Raphson's updating rule with regularized Hessian
      temp_matrix <- ginv(tripi.reg) %*% dipi

      # Optional step-size cap to further limit the update magnitude
       max_update <- 1.0
       update_norm <- sqrt(sum(temp_matrix^2))
       if (update_norm > max_update) {
         temp_matrix <- temp_matrix * (max_update / update_norm)
       }

      z.1[j,] <- z.0[j,] - alpha * temp_matrix

      # Shrink the step size using Armijo's Rule
      dk = (-1) * temp_matrix
      h_0 = t(z.0[j,] - mu.0) %*% Omega.0 %*% dk - x[j, 1:K] %*% dk +
        M[j] * as.numeric(t(dk) %*% exp(z.0[j,])) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1)

      if (obj(x[j, 1:K], z.1[j,], Omega.0, K) >
          (obj(x[j, 1:K], z.0[j,], Omega.0, K) + epsilon * alpha * h_0)) {
        alpha <- alpha / delta
      }

      time_rep <- time_rep + 1
    }
  }
  return(z.1)
}

# Update Z on all slices (one NR solve per slice).
renew_z <- function(x, z.0, Omega.0, m, alpha_0 = 1, delta = 5, epsilon_NR = 0.01, threshold_NR = 0.0001) {
  z_end <- lapply(1:m, function(i) {
    NR(x[[i]], z.0[[i]], Omega.0[[i]], alpha_0, delta, epsilon_NR, threshold_NR)
  })
  return(z_end)
}


# ============================================================================
# BETA-UPDATE — PD-guarded ADAM on the smooth penalized objective.
#
# Given Z, update the B-spline coefficients beta. Ffinal_beta is the (beta-part)
# objective value; dgrad_lt_h is its gradient w.r.t. one basis block; renew_beta_final
# runs one ADAM step per block with a nonmonotone-Armijo PD line search (the smoothed
# penalty does the shrinkage, so there is NO in-loop hard thresholding).
# ============================================================================

# Objective VALUE of the beta-controlled part, with a TRUE Cholesky log-det barrier:
#   result_two (-log det) + result_three (trace S*Omega) + result_four (R penalty).
# Returns +Inf if any assembled slice is not PD (the barrier). Slices are n_k/N-weighted.
Ffinal_beta <- function(beta, lambda, W, q, N_n, mu, m, x_b_spline_base,
                        Omega_list = NULL, S_list = NULL, w_slice = NULL) {
  if (is.null(w_slice)) w_slice <- rep(1/m, m)                 # P0-6: n_k/N weighting (=1/m for equal n)
  if (is.null(Omega_list)) Omega_list <- G_beta_Rcpp(beta, x_b_spline_base, m)
  logdets <- vapply(Omega_list, .logdet_chol, numeric(1))
  if (any(is.na(logdets))) return(Inf)
  result_two   <- -0.5 * sum(w_slice * logdets)
  result_three <- sum(vapply(seq_len(m), function(k) w_slice[k] * 0.5 * sum(S_list[[k]] * Omega_list[[k]]), numeric(1)))
  result_four  <- tv_pen_value(beta, lambda, W, q, N_n, mu)
  result_two + result_three + result_four
}

# PD-safe inverse: Path B keeps iterates PD, so chol2inv (true inverse) is exact
# and ~5-10x faster than ginv's SVD; fall back to ginv if a slice is not PD.
.inv_pd <- function(M) tryCatch(chol2inv(chol(.sym(M))), error = function(e) ginv(M))

# Gradient of (result_two + result_three) w.r.t. basis block beta_h. For slice k the
# per-slice term is B_h(t_k) * (-Omega(t_k)^{-1} + S(t_k)); slices are summed with
# weight w_slice[k] = n_k/N. free_diag=TRUE keeps the diagonal (halved, it enters once);
# free_diag=FALSE zeros the diagonal gradient (frozen diagonal).
# ★ P0-6/codex (2026-07-22): slice k weighted by w_slice[k] = n_k/N (the joint LNM likelihood),
# NOT the slice-uniform 1/m. For EQUAL n, w_slice = 1/m -> bit-identical (all sims unchanged);
# for unequal-n real data the beta ESTIMATOR now matches the reported n_k-weighted likelihood/IC.
dgrad_lt_h <- function(G_list, S_list, x_b_spline_base, m, P, h, free_diag, w_slice = NULL) {
  if (is.null(w_slice)) w_slice <- rep(1/m, m)
  temp <- matrix(0, P, P)
  for (k in 1:m) {
    Ginv <- .inv_pd(G_list[[k]])
    temp <- temp + w_slice[k] * x_b_spline_base[k,h] * (-Ginv + S_list[[k]])
  }
  if (free_diag) diag(temp) <- 0.5 * diag(temp) else diag(temp) <- 0
  temp
}

# One ADAM step per basis block on the SMOOTH corrected objective, block by block.
# For each block h: form the gradient (likelihood + smoothed penalty), take an ADAM
# step, and accept it only if (a) all assembled slices stay PD AND (b) a nonmonotone
# Armijo sufficient-decrease holds; else fall back to the negative-gradient direction
# (a true descent of Fcur); if neither is accepted after backtracking, stop this block.
renew_beta_final <- function(Z, beta, x_b_spline_base, J_n, m, P, q, lambda, W, N_n, mu,
                             initial_learning_rate, max_iterations = 20, tol = 1e-4,
                             free_diag = FALSE, pd_max_backtrack = 25L, pd_c1 = 1e-4,
                             pd_nonmono_K = 5L, w_slice = NULL) {
  if (is.null(w_slice)) w_slice <- rep(1/m, m)                 # P0-6: n_k/N slice weights (=1/m for equal n)
  temp_beta <- beta
  S_list <- S_Z_t(Z)
  Fcur <- function(bt, Om = NULL) Ffinal_beta(bt, lambda, W, q, N_n, mu, m, x_b_spline_base,
                                              Omega_list = Om, S_list = S_list, w_slice = w_slice)
  for (h in 1:J_n) {
    learning_rate <- initial_learning_rate
    m_adam <- matrix(0, P, P); v_adam <- matrix(0, P, P)
    F_hist <- rep(Fcur(temp_beta), pd_nonmono_K)
    for (i in 1:max_iterations) {
      G_list <- G_beta_Rcpp(temp_beta, x_b_spline_base, m)
      d1 <- dgrad_lt_h(G_list, S_list, x_b_spline_base, m, P, h, free_diag, w_slice = w_slice)
      d2 <- tv_pen_grad(temp_beta, lambda, W, q, N_n, mu)[[h]]
      grad <- d1 + d2
      if (!free_diag) diag(grad) <- 0
      grad <- (grad + t(grad)) / 2
      grad_norm <- sqrt(sum(grad^2)); if (grad_norm < tol) break
      m_adam <- 0.9*m_adam + 0.1*grad; v_adam <- 0.999*v_adam + 0.001*(grad^2)
      m_hat <- m_adam/(1-0.9^i); v_hat <- v_adam/(1-0.999^i)
      adam_step <- (learning_rate/(sqrt(v_hat)+1e-8))*m_hat
      adam_step <- (adam_step + t(adam_step))/2
      if (!free_diag) diag(adam_step) <- 0
      F_ref <- max(F_hist); accepted <- FALSE
      # (a) ADAM direction, PD + nonmonotone Armijo
      for (bt in 0:pd_max_backtrack) {
        alpha <- 0.5^bt
        cand_h <- (function(B){B<-B-alpha*adam_step;(B+t(B))/2})(temp_beta[[h]])
        cand <- temp_beta; cand[[h]] <- cand_h
        Om <- G_beta_Rcpp(cand, x_b_spline_base, m)
        if (.all_slices_pd(Om)) { Fnew <- Fcur(cand, Om=Om)
          if (is.finite(Fnew) && Fnew <= F_ref - pd_c1*alpha*grad_norm^2) { temp_beta <- cand; accepted <- TRUE
            .pd_diag_env$adam_accept <- .pd_diag_env$adam_accept + 1L; break } }
      }
      # (b) fallback: negative gradient (true descent of Fcur), strict Armijo
      if (!accepted) {
        F_old <- Fcur(temp_beta)
        for (bt in 0:pd_max_backtrack) {
          alpha <- initial_learning_rate*0.5^bt
          cand_h <- (function(B){B<-B-alpha*grad;(B+t(B))/2})(temp_beta[[h]])
          cand <- temp_beta; cand[[h]] <- cand_h
          Om <- G_beta_Rcpp(cand, x_b_spline_base, m)
          if (.all_slices_pd(Om)) { Fnew <- Fcur(cand, Om=Om)
            if (is.finite(Fnew) && Fnew <= F_old - pd_c1*alpha*grad_norm^2) { temp_beta <- cand; accepted <- TRUE
              .pd_diag_env$fallback_accept <- .pd_diag_env$fallback_accept + 1L; break } }
        }
      }
      if (!accepted) { .pd_diag_env$block_break <- .pd_diag_env$block_break + 1L; break }
      F_hist <- c(F_hist[-1], Fcur(temp_beta))
    }
  }
  temp_beta
}


# ============================================================================
# PENALTY — adaptive OVERLAPPING sliding-window group-lasso on the beta's.
#
# Each edge (i,j) has J_n B-spline coefficients (beta_1[i,j], ..., beta_{J_n}[i,j]),
# the coefficients of its Omega_ij(t) function. The penalty groups OVERLAPPING
# windows of q+1 adjacent coefficients:
#   gamma_g^{ij} = (beta_g[i,j], ..., beta_{g+q}[i,j]),   g = 1 .. N_n+1
# and penalizes each window's L2 norm ||gamma_g^{ij}|| (group lasso -> shrinks a
# whole time-local window to zero). Adaptive weights W_g[i,j] downweight already-
# large windows (oracle property). The non-smooth L2 norm is made differentiable
# by Nesterov smoothing psi_mu so the ADAM/gradient solver can use it:
#   psi_mu(t) = t - mu/2        if t >  mu
#             = t^2 / (2 mu)     if t <= mu       (psi_mu(t) -> t as mu -> 0)
#   d/dg psi_mu = g / max(||g||, mu)
# So the penalty value = sum_{i<j} sum_g lambda * W_g[i,j] * psi_mu(||gamma_g^{ij}||).
# Everything is vectorized over the upper-triangular edges (Bmat = n_edges x J_n).
# ============================================================================

.ut_idx <- function(P) which(upper.tri(matrix(0, P, P)))                        # linear indices of the upper triangle
.beta_to_Bmat <- function(beta, ut) vapply(beta, function(M) M[ut], numeric(length(ut)))  # n_edges x J_n coefficient matrix

# Penalty VALUE: sum over the N_n+1 overlapping windows of lambda * W_g * psi_mu(||gamma_g||).
tv_pen_value <- function(beta, lambda, W, q, N_n, mu) {
  P <- nrow(beta[[1]]); G <- N_n + 1L; ut <- .ut_idx(P)                         # G = number of overlapping windows
  Bmat <- .beta_to_Bmat(beta, ut); Wmat <- vapply(W, function(M) M[ut], numeric(length(ut)))
  tot <- 0
  for (g in 1:G) {
    gam <- Bmat[, g:(g+q), drop = FALSE]; ng <- sqrt(rowSums(gam^2))            # per-edge window norm ||gamma_g||
    psi <- ifelse(ng > mu, ng - mu/2, ng^2/(2*mu))                             # Nesterov-smoothed norm
    tot <- tot + lambda * sum(Wmat[, g] * psi)
  }
  tot
}

# Penalty GRADIENT: returns a list of J_n P x P matrices, d(penalty)/d(beta_h). A window g
# contributes to the q+1 blocks it spans; its gradient is coef * gamma_g / max(||gamma_g||, mu).
tv_pen_grad <- function(beta, lambda, W, q, N_n, mu) {
  P <- nrow(beta[[1]]); J_n <- length(beta); G <- N_n + 1L; ut <- .ut_idx(P)
  Bmat <- .beta_to_Bmat(beta, ut); Wmat <- vapply(W, function(M) M[ut], numeric(length(ut)))
  Gmat <- matrix(0, length(ut), J_n)                                          # accumulate per-edge x per-block gradient
  for (g in 1:G) {
    gam <- Bmat[, g:(g+q), drop = FALSE]; ng <- sqrt(rowSums(gam^2)); denom <- pmax(ng, mu)
    dpsi <- gam / denom; coef <- lambda * Wmat[, g]                           # smoothed-norm gradient x weight
    for (e in 0:q) Gmat[, g+e] <- Gmat[, g+e] + coef * dpsi[, e+1]            # scatter into the q+1 spanned blocks
  }
  lapply(1:J_n, function(h) { M <- matrix(0, P, P); M[ut] <- Gmat[, h]; M + t(M) })  # symmetric, diag 0
}

# FROZEN adaptive plug-in weights from a pilot beta: W_g[i,j] = (||gamma_g^{ij}(pilot)|| + eps)^(-r).
# Larger pilot windows get smaller weight (adaptive-lasso, exponent r); computed once and frozen.
tv_compute_W <- function(beta_pilot, q, N_n, r = 1.0, w_floor = 1e-6) {
  P <- nrow(beta_pilot[[1]]); G <- N_n + 1L; ut <- .ut_idx(P)
  Bmat <- .beta_to_Bmat(beta_pilot, ut)
  lapply(1:G, function(g) {
    ng <- sqrt(rowSums(Bmat[, g:(g+q), drop = FALSE]^2))
    Wg <- matrix(0, P, P); Wg[ut] <- (ng + w_floor)^(-r); Wg + t(Wg)
  })
}


# ============================================================================
# INIT / HELPERS — glasso->beta projection, per-slice S(Z), default threshold.
# ============================================================================

# Per-slice second moment S(Z) = centered sample covariance (used in the beta objective/gradient).
S_Z_t <- function(Z) {
  S_temp <- list()                     # output list, one S per slice
  m <- length(Z)
  P <- dim(Z[[1]])[2]
  for (i in 1:m) {
    Z_matrix <- Z[[i]]
    n_i <- dim(Z_matrix)[1]
    mu_t <- colMeans(Z_matrix)         # per-taxon mean
    centered_Z <- sweep(Z_matrix, 2, mu_t)   # center, then form the covariance
    S_t <- (t(centered_Z) %*% centered_Z) / n_i
    S_temp[[i]] <- S_t
  }
  return(S_temp)
}

# Project a per-slice Omega (from glasso) onto the B-spline basis to get an initial beta:
# for each (i,j), least-squares fit Omega_ij(t_v) ~ B(t_v) -> the J_n coefficients. Init only.
generate_beta <- function(Omega_t, P, m, q, J_n, x_b_spline_base) {
  # beta_temp: one coefficient vector per upper-triangular (i,j) entry
  beta_temp <- vector("list", P * (P + 1) / 2)   # store only the upper triangle (P*(P+1)/2 entries)
  idx <- 1
  for (i in 1:P) {                               # loop over each (i,j) entry
    for (j in i:P) {
      x_temp <- sapply(1:m, function(v){         # collect Omega_ij(t_v) across the m slices
        Omega_t[[v]][i, j]
      }
      )
      fit <- lm(x_temp ~ x_b_spline_base - 1)
      coefficients <- coef(fit)
      coefficients[is.na(coefficients)] <- 0
      beta_temp[[idx]] <- coefficients           # store the fitted coefficients
      idx <- idx + 1
    }
  }
  # scatter the coefficients into J_n P x P matrices (one matrix per basis function)
  beta_matrices <- lapply(1:J_n, function(k) {
    mat <- matrix(0, nrow = P, ncol = P)
    idx <- 1
    for (i in 1:P) {
      for (j in i:P) {
        mat[i, j] <- beta_temp[[idx]][k]
        mat[j, i] <- mat[i, j]                   # symmetry
        idx <- idx + 1
      }
    }
    return(mat)
  })
  return(beta_matrices)
}

# Percentile of |off-diagonal beta| (nonzero entries) — the default sel_threshold.
calculate_percentile <- function(matrices, percentile = 0.05) {
  # flatten all matrices, keep the nonzero off-diagonal magnitudes
  all_elements <- unlist(lapply(matrices, function(mat) {
    off_diag_elements <- abs(mat[upper.tri(mat) | lower.tri(mat)])   # absolute off-diagonal entries
    off_diag_elements <- off_diag_elements[off_diag_elements!= 0]    # drop exact zeros
    return(off_diag_elements)
  }))
  threshold <- quantile(all_elements, percentile)
  return(threshold)
}


# ============================================================================
# PD UTILITIES (Cholesky-based; a det>0 test misses ~1/3 of indefinite cases).
# ============================================================================
.sym <- function(M) (M + t(M)) / 2                                  # symmetrize

.pd_chol_ok <- function(M) {                                        # TRUE iff M is PD (Cholesky succeeds)
  tryCatch({ chol(.sym(M)); TRUE }, error = function(e) FALSE)
}

.min_eig <- function(M) {                                           # smallest eigenvalue (PD margin)
  min(eigen(.sym(M), symmetric = TRUE, only.values = TRUE)$values)
}

.all_slices_pd <- function(Omega_list) {                            # TRUE iff every slice is PD
  all(vapply(Omega_list, .pd_chol_ok, logical(1)))
}

# Cholesky log-det of an SPD matrix; NA if not PD (the barrier caller returns +Inf).
.logdet_chol <- function(M) {
  ch <- tryCatch(chol(.sym(M)), error = function(e) NULL)
  if (is.null(ch)) return(NA_real_)   # not PD
  2 * sum(log(diag(ch)))
}
