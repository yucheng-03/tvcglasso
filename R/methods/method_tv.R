# ============================================================================
# methods/method_tv.R — TVCGLasso, OUR method (base fit + optional refit).
#
# This is the publication estimator. `cfg$tv_refit` toggles the refit:
#   FALSE -> pass 1 only (penalized warm-start continuation): ROC/scores + the
#            per-lambda base intermediate.
#   TRUE  -> pass 1 + pass 2 (relaxed refit at EVERY lambda on the path): the
#            output ALSO carries, per lambda, the de-biased refit results with
#            the full BIC-tuning components (pre = penalized magnitudes on the
#            selected support with Z re-aligned; post = refitted magnitudes),
#            so any selector (AIC/BIC/eBIC, pre or post) is computable at
#            analysis WITHOUT re-fitting.
#
# refit does NOT change WHICH entries are nonzero (the support is fixed from
# pass 1), only their MAGNITUDES + the latent Z -> so the lambda-sweep ROC
# (edge = Omega_hat != 0) is identical with or without refit; refit's role is
# the de-biased DEPLOYED operating point (refit-BIC) + magnitude accuracy.
#
# Model-selection likelihood = the JOINT LNM density at the inferred Zhat
# (multinomial + Gaussian-graphical), via refit_joint_nll_average -> the IC
# helpers. NOT the Z_0 Gaussian-only loglik (a recorded statistical error).
#
# Requires (sourced ONCE by the master driver): R/tvcglasso.R, R/refit.R,
#   R/roc_utils.R; libraries splines, glasso, Matrix, MASS, Rcpp.
# ============================================================================

# Data-adaptive lambda grid (CGLasso recipe): lambda_max = 2 * max over slices of the max
# |off-diagonal| of cov(ALR Z_0), guaranteeing the empty graph at the sparse end; log-spaced
# down to lambda_max/1000 (dense end). Anchored on the max across slices so ONE grid spans
# empty->dense for every slice. Reaches both corners with REAL lambda, no wasted all-empty
# points. (Z_0 built with the option-2 pseudocount, matching the engine's init.)
tv_lambda_grid <- function(X, n_lambda) {
  P <- ncol(X[[1]]) - 1L
  Xpool <- do.call(rbind, X); off2 <- colMeans(Xpool / rowSums(Xpool)) * (P + 1L)
  maxc <- max(vapply(X, function(Xi) {
    Xa <- t(t(Xi) + off2); Z <- log(Xa[, -(P + 1)] / Xa[, P + 1]); S <- cov(Z)
    max(abs(S[upper.tri(S)]))
  }, numeric(1)))
  exp(seq(log(2 * maxc), log(maxc / 1000), length.out = n_lambda))
}

run_method_tv <- function(dat, cfg, eval_slices) {
  X <- dat$X; xseq <- dat$x_sequence; true <- dat$true_Omega_list
  m <- length(X); P <- ncol(X[[1]]) - 1
  q <- dat$tv_q %||% cfg$tv_q %||% 2L; N_n <- dat$tv_N_n %||% cfg$tv_N_n %||% 1L   # per-cell basis (scales with m)
  mid <- seq(xseq[1], xseq[m], length.out = N_n + 2)[-c(1, N_n + 2)]
  basis <- bs(xseq, degree = q, knots = mid, Boundary.knots = c(xseq[1], xseq[m]), intercept = TRUE)
  # lambda grid: use the config's fixed grid if given, else data-adaptive (reaches the corners)
  lam_grid <- cfg$lambda_grid_tv %||% tv_lambda_grid(X, cfg$tv_n_lambda %||% 24L)

  # ---- pass 1: penalized warm-start continuation (SELECTS the support per lambda) ----
  fits <- tv_warm_path(X, lam_grid, xseq, q = q, N_n = N_n,
                       init_mode   = cfg$tv_init_mode   %||% "diag",
                       sel_type    = cfg$tv_sel_type    %||% "hard",
                       weight_mode = cfg$tv_weight_mode %||% "glasso",
                       free_diag   = cfg$tv_free_diag   %||% TRUE,   # (C) ADOPTED: free (estimated) diagonal
                       Max_iterations = cfg$tv_max_iter %||% 150L,
                       sel_threshold  = cfg$tv_sel_threshold,   # NULL => engine 20th-pct/floor-0.01 default
                       ckpt_file      = dat$tv_ckpt_file)       # NULL => no checkpoint
  lam <- as.numeric(names(fits))                               # ASCENDING lambda
  Om_by_lambda <- lapply(fits, function(f) G_beta_Rcpp(f$beta, basis, m))

  # refit changes the COEFFICIENT support's magnitudes, not the lambda-sweep support -> this IS
  # the TV method's lambda-sweep ROC (base path; the ROC is unaffected by refit by construction).
  rp     <- roc_from_lambda_path(Om_by_lambda, lam, true, eval_slices)
  scores <- entry_lambda_scores(Om_by_lambda, lam, m)          # symmetric score-sweep

  # RAW counts (P0-5: pseudocount is init-only; the refit likelihood uses raw counts).
  # ONE IC (P0-6: refit_information_criteria is n_k-weighted -> correct for equal AND unequal n).
  X_work <- refit_prepare_counts(X)                            # = raw X
  n_per_slice <- vapply(X, nrow, integer(1)); N <- sum(n_per_slice)
  ic_fun <- refit_information_criteria

  # ---- per-lambda BASE geometry (df/edges/min_eig) + the LAMBDA-ALIGNED operating point (FPR/TPR).
  #      `detail` is in ASCENDING-lambda order, aligned INDEX-FOR-INDEX with fits / Om_by_lambda /
  #      refit[[i]] (lambda_index = i). `rp$roc` is FPR-SORTED and drops lambda, so it cannot map a
  #      SELECTED lambda back to its operating point; `detail$FPR`/`detail$TPR` here CAN -> the
  #      as-deployed point (e.g. refit-BIC argmin index) + any tuning table is a DIRECT LOOKUP
  #      `detail[<sel_index>, c("FPR","TPR")]`, with NO graph reconstruction. (Same edge rule as the
  #      ROC: edge = Omega_hat_ij(t_k) != 0, pooled over eval_slices.)
  #      The base "without-refit" model-selection IC is NOT computed here (P1-4): the correct base IC
  #      = the refit's Z-ALIGNED `pre` criterion (below), not post-threshold-beta + pre-threshold-Z. ----
  detail <- do.call(rbind, lapply(seq_along(fits), function(i) {
    f <- fits[[i]]; Oml <- Om_by_lambda[[i]]
    df    <- sum(vapply(f$beta, function(b) sum(b[upper.tri(b)] != 0), integer(1)))
    edges <- sum(vapply(Oml, function(M) { d <- M; diag(d) <- 0; sum(d[upper.tri(d)] != 0) }, integer(1)))
    ft  <- edge_fpr_tpr(Oml, true, eval_slices)                  # lambda-aligned operating point (edge = Omega!=0)
    tpb <- f$threshold_per_block                                 # P1-5: per-block effective threshold audit
    data.frame(lambda = lam[i], df = df, n_edges = edges,
               FPR = unname(ft["FPR"]), TPR = unname(ft["TPR"]),  # <- lambda-aligned: deployed point = direct lookup
               min_eig     = min(vapply(Oml, .min_eig, numeric(1))),
               ls_fallback = f$ls_fallback,
               converged = isTRUE(f$converged), exit_reason = f$exit_reason %||% NA_character_,   # base convergence contract
               n_outer = f$n_outer %||% NA_integer_,
               final_dOm = f$final_dOm %||% NA_real_, final_dZ = f$final_dZ %||% NA_real_,
               thr_min = if (all(is.na(tpb))) NA_real_ else min(tpb, na.rm = TRUE),
               thr_n_reverted = sum(f$threshold_reverted %||% FALSE),
               row.names = NULL)
  }))

  # ---- pass 2: refit at EVERY lambda (cfg$tv_refit) ----
  # For each lambda: support = f$beta nonzeros; re-estimate active magnitudes + Z on the FIXED
  # support (unpenalized). Store BOTH pre (penalized magnitudes, Z re-aligned) and post (refitted)
  # criteria with every numeric BIC component so the deployed point (refit-BIC argmin) is
  # recomputable at analysis. Each lambda's refit is wrapped so one failure does not lose the fit.
  refit <- NULL
  if (isTRUE(cfg$tv_refit %||% TRUE)) {
    ckdir <- dat$refit_ckpt_dir                                  # per-(cell,seed) dir; NULL => tempfiles
    .refit_one <- function(i) {
      f <- fits[[i]]
      ckf  <- if (!is.null(ckdir)) file.path(ckdir, sprintf("refit_lam%02d.rds", i)) else tempfile(fileext = ".rds")
      plog <- if (!is.null(ckdir)) file.path(ckdir, sprintf("refit_lam%02d.log", i)) else tempfile(fileext = ".log")
      fit <- tryCatch(
        refit_fixed_support(X_work = X_work, Z_start = f$Z, beta_start = f$beta, basis = basis,
                            checkpoint_file = ckf, progress_log = plog,
                            max_outer       = cfg$tv_refit_max_outer %||% cfg$refit_max_outer %||% 500L,  # TV-OWN key (decoupled from CGLasso's cap); 500 = 2026-07-24 default
                            inner_max       = cfg$refit_inner_max      %||% 10L,
                            z_align_max     = cfg$refit_z_align_max    %||% 160L,
                            conv_tol_Z      = cfg$refit_conv_tol_Z     %||% 5e-5,
                            conv_tol_z_grad = cfg$refit_conv_tol_z_grad %||% 1e-4),
        error = function(e) structure(list(msg = conditionMessage(e)), class = "refit_error"))
      if (inherits(fit, "refit_error"))
        return(list(lambda = lam[i], lambda_index = i, df = NA_integer_, error = fit$msg))
      df <- fit$df
      pre_ic  <- ic_fun(fit$pre$criterion$total,  n_per_slice, df, P)
      post_ic <- ic_fun(fit$post$criterion$total, n_per_slice, df, P)
      comp <- function(cr, ic) list(
        mult = cr$multinomial, gauss = cr$neg_logdet + cr$trace, nll = cr$total,
        neg2loglik = unname(ic["neg2loglik"]), AIC = unname(ic["AIC"]),
        BIC = unname(ic["BIC"]), eBIC = unname(ic["eBIC"]))
      list(lambda = lam[i], lambda_index = i, df = df,
           support = fit$support,
           pre  = comp(fit$pre$criterion,  pre_ic),
           post = comp(fit$post$criterion, post_ic),
           beta_post  = fit$post$beta,                           # refitted (de-biased) beta
           Omega_post = fit$post$criterion$Omega,                # refitted Omega(t_k) per slice (completeness)
           Z_post     = fit$post$Z,                              # refitted latent Z per slice (completeness)
           min_eig   = min(vapply(fit$post$criterion$Omega, .min_eig, numeric(1))),
           z_align_converged = isTRUE(fit$z_align_converged),
           converged = isTRUE(fit$converged), exit_reason = fit$exit_reason,
           max_active_grad = fit$max_active_grad, max_z_grad = fit$max_z_grad)
    }
    # ---- OPTIONAL PARALLEL REFIT (env-gated; NOT a cfg key, so provenance$config is byte-unchanged
    #      and the stored rds format is identical either way). The refits at different lambdas are
    #      INDEPENDENT: each starts from its OWN base fit fits[[i]] (no cross-lambda warm start), uses
    #      no RNG, and writes its own refit_lam%02d checkpoint/log -- so evaluating them concurrently
    #      yields the SAME list, in the SAME order, with the SAME values. TVCG_REFIT_CORES unset or 1
    #      => the plain lapply, i.e. the exact previous code path. ----
    .rc <- suppressWarnings(as.integer(Sys.getenv("TVCG_REFIT_CORES", "1")))
    if (is.na(.rc) || .rc < 1L) .rc <- 1L
    refit <- if (.rc > 1L)
      parallel::mclapply(seq_along(fits), .refit_one, mc.cores = .rc, mc.preschedule = FALSE)
    else lapply(seq_along(fits), .refit_one)
    # A forked worker that dies outside R's condition system (OOM/segfault) yields a try-error rather
    # than a list; normalise it to the SAME error element the serial tryCatch produces, so downstream
    # (ic_col / deployed / the stored payload) sees one and only one shape. No-op in the serial path.
    refit <- lapply(seq_along(refit), function(i) {
      x <- refit[[i]]
      if (is.list(x) && !is.null(x$lambda_index)) return(x)
      list(lambda = lam[i], lambda_index = i, df = NA_integer_,
           error = if (inherits(x, "try-error")) as.character(x) else "parallel refit worker failed")
    })
  }

  # ---- flatten the refit pre/post IC into `detail` (index-aligned with refit[[i]] / fits) so the
  #      tuning table + the "refit shifts BIC toward oracle" figure read a FLAT column, not the deep
  #      payload (SHELL contract: detail carries AIC/BIC/eBIC). AIC/BIC/eBIC = refit-POST (the native
  #      refit-BIC selector's IC); *_pre = penalized pre-refit. df is already a detail column (same support). ----
  if (!is.null(refit)) {
    ic_col <- function(pass, key) vapply(refit, function(x) {
      p <- x[[pass]]; if (is.null(p) || is.null(p[[key]])) NA_real_ else p[[key]] }, numeric(1))
    detail$AIC     <- ic_col("post", "AIC"); detail$BIC     <- ic_col("post", "BIC"); detail$eBIC     <- ic_col("post", "eBIC")
    detail$AIC_pre <- ic_col("pre",  "AIC"); detail$BIC_pre <- ic_col("pre",  "BIC"); detail$eBIC_pre <- ic_col("pre",  "eBIC")
  }

  # ---- as-deployed operating point = the NATIVE selector (refit-BIC argmin over CONVERGED lambdas).
  #      The refit does not change the support, so the deployed edge-set = the base graph at that lambda;
  #      report its FPR/TPR/F1/precision/n_edges from the lambda-aligned (edge x slice) confusion. ----
  deployed <- NULL
  if (!is.null(refit)) {
    bic_conv <- vapply(refit, function(x) if (!is.null(x$post) && isTRUE(x$converged)) x$post$BIC else NA_real_, numeric(1))
    bic_all  <- vapply(refit, function(x) if (!is.null(x$post)) x$post$BIC else NA_real_, numeric(1))
    di <- if (any(is.finite(bic_conv))) which.min(bic_conv)       # prefer a CONVERGED refit lambda
          else if (any(is.finite(bic_all)))  which.min(bic_all)   # else best BIC overall (flagged not-converged)
          else NA_integer_
    if (!is.na(di)) {
      dep_conv <- isTRUE(refit[[di]]$converged)
      P_ <- nrow(Om_by_lambda[[di]][[1]]); ut <- upper.tri(matrix(0, P_, P_)); tp <- fp <- fn <- tn <- 0
      for (k in eval_slices) {
        pe <- Om_by_lambda[[di]][[k]][ut] != 0; te <- true[[k]][ut] != 0
        tp <- tp + sum(pe & te); fp <- fp + sum(pe & !te); fn <- fn + sum(!pe & te); tn <- tn + sum(!pe & !te)
      }
      deployed <- list(selector = "refit-BIC", lambda_index = di, lambda = lam[di],
                       FPR = fp / max(1, fp + tn), TPR = tp / max(1, tp + fn),
                       precision = if (tp + fp == 0) NA_real_ else tp / (tp + fp),
                       F1 = if (2*tp + fp + fn == 0) NA_real_ else 2*tp / (2*tp + fp + fn),
                       n_edges = tp + fp, df = refit[[di]]$df,
                       refit_converged = dep_conv)               # FALSE => deployed lambda's refit hit the iter cap
    }
  }

  list(method = "tvcglasso", convention = "lambda-sweep",
       # ---- SHELL (shared across all 4 methods) ----
       roc = rp$roc, auc = rp$auc, scores = scores, detail = detail,
       deployed = deployed,                                      # as-deployed operating point (native selector)
       n_per_slice = n_per_slice, N = N, q = q, N_n = N_n,
       threshold_audit = lapply(fits, function(f) list(per_block = f$threshold_per_block, reverted = f$threshold_reverted)),  # P1-5
       # ---- PAYLOAD SLOT: this method's native fitted objects, BOTH passes, EVERY lambda (completeness) ----
       estimate = list(
         base  = list(beta  = lapply(fits, function(f) f$beta),  # base fit: beta (reconstruct Omega(t) at ANY t),
                      Omega = Om_by_lambda,                      #   Omega(t_k) per slice, and the latent Z per slice
                      Z     = lapply(fits, function(f) f$Z)),
         refit = refit))                                        # refit: beta_post/Omega_post/Z_post + pre/post IC
}
