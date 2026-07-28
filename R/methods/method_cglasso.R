# ============================================================================
# methods/method_cglasso.R — CGLasso static control (lambda-sweep), REAL solver.
# Requires (sourced ONCE by the master driver):
#   R/baselines/cglasso_static.R  (fit_cglasso_static -> vendored Compo_glasso)
#   R/baselines/stars_cglasso.R   (stars_select_cglasso, .cg_confusion)
#   R/baselines/refit_cglasso.R   (cglasso_refit_one_rho, .cg_per_slice_ic)
#   R/roc_utils.R                 (roc_from_lambda_path, entry_lambda_scores, edge_fpr_tpr)
#   libraries: glasso, huge, MASS, propagate
#
# This is the REAL Comp-gLASSO (Tian et al. 2023, JASA) — the LNM joint solver
# (NR z-update + glasso) — fit INDEPENDENTLY per time slice (the no-time-awareness
# static control). We only REPRODUCE their static method + ADD a relaxed refit.
#
# ★ PER-SLICE, DATA-ADAPTIVE, INDEPENDENT (2026-07-22, faithful reproduction):
#   * rho grid: EACH slice gets its OWN data-adaptive grid from its OWN max(cov(ALR-z))
#     (Yuan's Functions_ROC.R recipe applied per slice) — cglasso_static.R.
#   * ROC: swept by the common rho-INDEX (relative regularisation level), pooled over
#     eval_slices; support-based (edge = Omega != 0), so REFIT-INVARIANT.
#   * deployed point: PER-SLICE INDEPENDENT — each slice k picks the rho minimising ITS
#     OWN BIC = -2*ell_k + d_k*log(n_k) (Yuan's static BIC on n_k samples), pre AND post
#     refit; the deployed net is assembled from the per-slice selected rho. This is the
#     most faithful "each time point is Yuan's static CGLasso tuned by its own BIC".
#
# refit (pass 2, cfg$cglasso_refit): relaxed refit at EVERY rho index via the CGLasso
# adapter (refit_cglasso.R) — the TV refit core with basis = I_m (=> beta_k = Omega_k,
# slices independent), RAW counts in the working likelihood, joint-LNM-at-Zhat per-slice
# BIC. Refit re-estimates the active per-slice magnitudes + Z on the FIXED pass-1 support
# -> it does NOT move the support, so the ROC is unchanged; its role is the de-biased
# deployed operating point. The output carries, per rho index, the pre/post per-slice BIC
# components + diagnostics (免re-fit selection at analysis).
# ============================================================================

run_method_cglasso <- function(dat, cfg, eval_slices) {
  m <- length(dat$X); P <- ncol(dat$X[[1]]) - 1L
  # per-slice base checkpoint dir + tail-able progress log (P0-4: resume-safe base path)
  base_ckpt <- if (!is.null(dat$refit_ckpt_dir)) file.path(dat$refit_ckpt_dir, "cglasso_base") else NULL
  base_plog <- if (!is.null(dat$refit_ckpt_dir)) file.path(dat$refit_ckpt_dir, "cglasso_base_progress.log") else NULL

  # rho grid: per-slice data-adaptive by default (paper's Functions_ROC.R per slice).
  # rho-grid resolution: 40 for TESTING (covers low-FPR, ~40% cheaper than the paper's 70;
  # bump back to 70 for the final publication run). cfg$cglasso_rho_list forces one common grid.
  f <- fit_cglasso_static(dat$X,
                          length_rholist = cfg$cglasso_length_rholist %||% 40L,
                          rho_hi_mult    = cfg$cglasso_rho_hi_mult %||% 2,
                          rho_lo_div     = cfg$cglasso_rho_lo_div  %||% 1000,  # larger -> denser dense-end -> FPR->1
                          option         = cfg$cglasso_option   %||% 2,
                          max_iter       = cfg$cglasso_max_iter %||% 50L,
                          rho_list       = cfg$cglasso_rho_list,   # NULL -> per-slice data-adaptive
                          verbose        = FALSE,
                          ckpt_dir       = base_ckpt,
                          progress_log   = base_plog)

  # ROC by rho-INDEX; rp$detail carries index-aligned FPR/TPR (P0-8: kept, joinable).
  rp     <- roc_from_lambda_path(f$Om_by_lambda, f$lambda, dat$true_Omega_list, eval_slices)
  scores <- entry_lambda_scores(f$Om_by_lambda, f$lambda, m)

  # per-index path table (P0-8): index, per-slice rho, index-aligned FPR/TPR, n_edges, +
  # the per-index per-slice support mask so ANY rho reverse-maps to its graph and (FPR,TPR).
  n_rho <- length(f$Om_by_lambda)
  ft_by_index <- lapply(f$Om_by_lambda, function(Oml) edge_fpr_tpr(Oml, dat$true_Omega_list, eval_slices))
  support_by_index <- lapply(f$Om_by_lambda, function(Oml)
    lapply(Oml, function(M) { d <- M; diag(d) <- 0; which(d != 0 & upper.tri(d), arr.ind = TRUE) }))
  detail <- data.frame(
    index   = seq_len(n_rho),
    FPR     = vapply(ft_by_index, function(x) unname(x["FPR"]), numeric(1)),
    TPR     = vapply(ft_by_index, function(x) unname(x["TPR"]), numeric(1)),
    n_edges = sapply(f$Om_by_lambda, function(Oml)
                sum(sapply(Oml, function(M) { d <- M; diag(d) <- 0; sum(d[upper.tri(d)] != 0) }))))
  # LEGACY diagnostic ONLY (Z_0 Gaussian-only; NOT the CGLasso likelihood; NOT used for any
  # selection — the deployed point uses the joint-LNM-at-Zhat per-slice BIC below). Kept + renamed
  # so nothing downstream mistakes it for the CGLasso likelihood/AIC (codex P1-5).
  detail$loglik_Z0_LEGACY <- sapply(f$Om_by_lambda, function(Oml) sum(sapply(seq_along(Oml), function(k) {
    Om <- (Oml[[k]] + t(Oml[[k]])) / 2; Z <- dat$Z_0[[k]]; nk <- nrow(Z); S <- crossprod(Z) / nk
    ev <- eigen(Om, symmetric = TRUE, only.values = TRUE)$values
    nk * (sum(log(pmax(ev, 1e-8))) - sum(S * Om)) / 2 })))

  # StARS (CGLasso's NATIVE selector) is DEFERRED to the real-data leg. It assumes ONE common
  # rho grid across slices; the per-slice data-adaptive rho path (2026-07-22) needs a per-slice
  # StARS rework, so hard-error if enabled rather than silently mis-selecting on the index.
  sel_point <- NULL
  if (isTRUE(cfg$cglasso_stars))
    stop("cglasso_stars is DEFERRED (real-data leg): per-slice data-adaptive rho needs a per-slice StARS rework; do not enable it for the sim.")

  # ---- refit (pass 2): relaxed refit at EVERY rho index (cfg$cglasso_refit) ----
  # basis = I_m, raw counts, per-slice fixed support; store per-index pre/post per-slice BIC
  # + diagnostics so the deployed point is recomputable at analysis. Each index wrapped so one
  # failure does not lose the fit. Refit does NOT move the support -> ROC unchanged.
  refit <- NULL; deployed <- NULL
  n_per_slice <- vapply(dat$X, nrow, integer(1)); N <- sum(n_per_slice)
  if (isTRUE(cfg$cglasso_refit %||% TRUE)) {
    ckdir <- if (!is.null(dat$refit_ckpt_dir)) file.path(dat$refit_ckpt_dir, "cglasso") else NULL
    if (!is.null(ckdir)) dir.create(ckdir, showWarnings = FALSE, recursive = TRUE)
    # ★ PARALLEL (2026-07-22): the refit at each rho is INDEPENDENT (own pass-1 Omega, own checkpoint +
    # log file) -> fit them across SLURM_CPUS_PER_TASK cores. This is the run's bottleneck (70 rho x up
    # to max_outer iters); mclapply gives ~ncores speedup. Serial (1 core) if ncores==1.
    .ncr <- { s <- Sys.getenv("SLURM_CPUS_PER_TASK"); if (nzchar(s)) max(1L, as.integer(s)) else 1L }
    .refit_one <- function(i) {
      ckf  <- if (!is.null(ckdir)) file.path(ckdir, sprintf("refit_lam%02d.rds", i)) else tempfile(fileext = ".rds")
      plog <- if (!is.null(ckdir)) file.path(ckdir, sprintf("refit_lam%02d.log", i)) else tempfile(fileext = ".log")
      r <- tryCatch(
        cglasso_refit_one_rho(dat$X, f$Om_by_lambda[[i]], ckf, plog,
                              option      = cfg$cglasso_option    %||% 2L,
                              max_outer   = cfg$cglasso_refit_max_outer %||% cfg$refit_max_outer %||% 50L,   # CGLasso-OWN key; 50 = Compo_glasso's own max_iter
                              inner_max   = cfg$refit_inner_max   %||% 10L,
                              z_align_max = cfg$refit_z_align_max %||% 160L,
                              conv_tol_Z      = cfg$refit_conv_tol_Z      %||% 5e-5,
                              conv_tol_z_grad = cfg$refit_conv_tol_z_grad %||% 1e-4,
                              # OURS, disclosed: ridge on the pass-2 slice covariance (CGLasso-OWN
                              # key). See the block above cg_refit_glasso_slice; 0 = unridged.
                              refit_ridge_eps = cfg$cglasso_refit_ridge_eps %||% 0.01),
        error = function(e) structure(list(msg = conditionMessage(e)), class = "refit_error"))
      if (inherits(r, "refit_error"))
        return(list(index = i, rho_per_slice = f$rho_matrix[i, ], df = NA_integer_, error = r$msg))
      comp <- function(cr, ic) list(
        mult = cr$multinomial, gauss = cr$neg_logdet + cr$trace, nll = cr$total,
        neg2loglik = unname(ic["neg2loglik"]), AIC = unname(ic["AIC"]),
        BIC_slice = unname(ic["BIC_slice"]), eBIC_slice = unname(ic["eBIC_slice"]),
        BIC_pooled = unname(ic["BIC_pooled_AUDIT"]))
      post_ok <- !is.null(r$post)   # NULL when the constrained MLE does not exist (support not estimable at this n); pre still valid
      list(index = i, rho_per_slice = f$rho_matrix[i, ], df = r$df_total, d_per_slice = r$d_per_slice,
           pre  = comp(r$pre$criterion,  r$pre_ic),
           post = if (post_ok) comp(r$post$criterion, r$post_ic) else NULL,
           pre_ic_slice = r$pre_ic_slice, post_ic_slice = r$post_ic_slice,   # the PER-SLICE selection input
           Omega_post = if (post_ok) r$post$criterion$Omega else NULL,       # refitted per-slice Omega (basis=I_m => = beta_post[[k]])
           pd_repair_scale = r$pd_repair_scale,
           min_eig   = if (post_ok) min(vapply(r$post$criterion$Omega, .min_eig, numeric(1))) else NA_real_,
           z_align_converged = isTRUE(r$z_align_converged),
           converged = isTRUE(r$converged), exit_reason = r$exit_reason,
           max_active_grad = r$max_active_grad, max_z_grad = r$max_z_grad)
    }
    refit <- if (.ncr > 1L) parallel::mclapply(seq_along(f$Om_by_lambda), .refit_one, mc.cores = .ncr, mc.preschedule = FALSE)
             else lapply(seq_along(f$Om_by_lambda), .refit_one)
    .re <- which(vapply(refit, function(x) inherits(x, "try-error"), logical(1)))   # mclapply returns try-error
    if (length(.re)) for (i in .re) refit[[i]] <- list(index = i, df = NA_integer_, error = as.character(attr(refit[[i]], "condition")$message))

    # ---- PER-SLICE INDEPENDENT deployed point (pre & post refit) ----
    # each slice k picks the rho index minimising ITS OWN BIC (-2*ell_k + d_k*log n_k); the
    # deployed net = assembly of the per-slice selected Omega_k. The point sits OFF the pooled
    # index-ROC (each slice at its own regularisation level) — the honest as-deployed graph.
    #
    # ★ DEPLOYMENT SELECTION (2026-07-23, user-mandated re-design): the deployed point is the
    # per-slice BIC-minimiser over the refit path (Yuan's native BIC selector applied to the
    # de-biased fit). A rho is a legitimate candidate for a slice iff its refit produced a WELL-
    # DEFINED criterion there -- i.e. the per-slice Omega is POSITIVE-DEFINITE and its BIC is
    # FINITE. Nothing else is required: whether the refit strictly CONVERGED (z-align / outer) is
    # RECORDED for transparency (deployed$*$sel_z_align_converged / sel_converged) but is NOT used
    # to FILTER candidates. Rationale (user's fairness principle: honest + consistent, respect the
    # original method): at low depth NOTHING strictly converges (latent-Z jitter floor), so a
    # convergence FILTER would leave NO candidate; but a non-PD refit (the constrained MLE does not
    # exist, e.g. a dense support at n<P) genuinely has no defined BIC and is naturally excluded.
    # The empty support is always PD + finite -> there is ALWAYS >= 1 candidate, so the deployed
    # point is always DEFINED (honestly empty/sparse at low depth = CGLasso's real static-baseline
    # weakness there), with its convergence status disclosed rather than the point suppressed.
    pd_floor <- cfg$cglasso_pd_floor %||% 1e-8
    fin_ic   <- function(ic) !is.null(ic) && all(is.finite(ic$BIC))
    ok_post <- vapply(refit, function(r) is.null(r$error) && !is.null(r$post_ic_slice) &&
                        fin_ic(r$post_ic_slice) && is.finite(r$min_eig) && r$min_eig >= pd_floor,
                      logical(1))   # candidate iff Omega PD + BIC finite; convergence NOT required (recorded only)
    ok_pre  <- vapply(refit, function(r) is.null(r$error) && !is.null(r$pre_ic_slice) &&
                        fin_ic(r$pre_ic_slice), logical(1))   # pre = penalized pass-1 (PD-repaired); needs only a finite BIC
    Pdim <- nrow(dat$true_Omega_list[[1]])
    pick_perslice <- function(ic_field, source, ok_idx) {
      sel_idx <- rep(NA_integer_, m); sel_rho <- rep(NA_real_, m); Om_sel <- vector("list", m); no_valid <- 0L
      sel_zac <- rep(NA, m); sel_conv <- rep(NA, m)   # RECORD (not filter) the selected rho's convergence status
      for (k in seq_len(m)) {
        bic_k <- if (length(ok_idx)) vapply(ok_idx, function(j) refit[[j]][[ic_field]]$BIC[k], numeric(1)) else numeric(0)
        good  <- is.finite(bic_k)
        if (!any(good)) { Om_sel[[k]] <- matrix(0, Pdim, Pdim); no_valid <- no_valid + 1L; next }
        jbest <- ok_idx[good][which.min(bic_k[good])]   # index j ascending = DENSEST->sparsest; ties (measure-zero) -> denser
        sel_idx[k] <- jbest; sel_rho[k] <- f$rho_matrix[jbest, k]
        sel_zac[k] <- isTRUE(refit[[jbest]]$z_align_converged); sel_conv[k] <- isTRUE(refit[[jbest]]$converged)
        Om_sel[[k]] <- if (source == "pass1") f$Om_by_lambda[[jbest]][[k]] else refit[[jbest]]$Omega_post[[k]]
      }
      cm <- .cg_confusion(Om_sel[eval_slices], dat$true_Omega_list[eval_slices])
      list(sel_index = sel_idx, sel_rho = sel_rho, Omega = Om_sel, n_no_valid = no_valid,
           sel_z_align_converged = sel_zac, sel_converged = sel_conv,   # DISCLOSED (transparency), not used to filter
           FPR = unname(cm["FPR"]), TPR = unname(cm["TPR"]), F1 = unname(cm["F1"]), edges = unname(cm["edges"]))
    }
    idx_pre <- which(ok_pre); idx_post <- which(ok_post)
    deployed <- list(
      selector = "per-slice independent BIC (-2*ell_k + d_k*log n_k) over the refit path; candidate iff Omega PD + BIC finite; convergence DISCLOSED not filtered",
      n_valid_pre = length(idx_pre), n_valid_post = length(idx_post), n_rho_total = length(refit),
      pre  = pick_perslice("pre_ic_slice",  "pass1", idx_pre),
      post = pick_perslice("post_ic_slice", "refit", idx_post))
    deployed$incomplete <- (deployed$post$n_no_valid > 0L) || (deployed$pre$n_no_valid > 0L)
  }

  # ---- SHARED `detail` CONTRACT (2026-07-25): flatten the per-rho quantities that already live
  # in `refit` into the standard per-path table, so a client reads df / min_eig / convergence /
  # AIC,BIC,eBIC (post = refit, *_pre = penalized pass-1) in the SAME place for every method
  # (cf. TV's detail). Pure OUTPUT ASSEMBLY -- values are copied, nothing is recomputed. ----
  if (!is.null(refit)) {
    g <- function(fn) vapply(refit, function(r) {
      v <- tryCatch(fn(r), error = function(e) NA_real_)
      if (is.null(v) || length(v) != 1L) NA_real_ else as.numeric(v)
    }, numeric(1))
    gl <- function(fn) vapply(refit, function(r) {
      v <- tryCatch(fn(r), error = function(e) NA); if (is.null(v)) NA else isTRUE(v)
    }, logical(1))
    gc_ <- function(fn) vapply(refit, function(r) {
      v <- tryCatch(fn(r), error = function(e) NA_character_)
      if (is.null(v)) NA_character_ else as.character(v)[1]
    }, character(1))
    detail$df              <- g(function(r) r$df)
    detail$min_eig         <- g(function(r) r$min_eig)
    detail$converged       <- gl(function(r) r$converged)
    detail$z_align_converged <- gl(function(r) r$z_align_converged)
    detail$exit_reason     <- gc_(function(r) r$exit_reason)
    detail$max_active_grad <- g(function(r) r$max_active_grad)
    detail$max_z_grad      <- g(function(r) r$max_z_grad)
    # POST (refit, de-biased) and PRE (penalized pass-1) pooled ICs, both from the joint-LNM-at-Zhat
    # likelihood. The DEPLOYED selector uses the PER-SLICE BIC in refit[[i]]$post_ic_slice; these
    # pooled columns are the path-level summary that mirrors the other methods' detail table.
    detail$AIC       <- g(function(r) r$post$AIC);        detail$AIC_pre  <- g(function(r) r$pre$AIC)
    detail$BIC       <- g(function(r) r$post$BIC_slice);  detail$BIC_pre  <- g(function(r) r$pre$BIC_slice)
    detail$eBIC      <- g(function(r) r$post$eBIC_slice); detail$eBIC_pre <- g(function(r) r$pre$eBIC_slice)
    detail$neg2loglik     <- g(function(r) r$post$neg2loglik)
    detail$neg2loglik_pre <- g(function(r) r$pre$neg2loglik)
    detail$refit_error    <- gc_(function(r) r$error)
  }

  # SAVE ALL PER-SLICE INTERMEDIATE (免re-fit for any future ROC/tuning change): the FULL
  # pass-1 per-index per-slice Omega (magnitudes, not just the support) + the refit Omega_post
  # (in `refit`) + per-slice rho + per-slice pre/post NLL/AIC/BIC/eBIC (in `refit`). So a later
  # switch of the ROC rule (support- vs magnitude-based) OR the tuning method needs NO re-fit.
  # `estimate` = the shared PAYLOAD-slot name (same content as Om_by_lambda + refit, which are
  # kept under their original names so existing analysis scripts keep working).
  list(method = "CGLasso", convention = "lambda-sweep",
       roc = rp$roc, auc = rp$auc, scores = scores, detail = detail,
       Om_by_lambda = f$Om_by_lambda, support_by_index = support_by_index, rho_matrix = f$rho_matrix,
       sel_point = sel_point, deployed = deployed,
       estimate = list(base = f$Om_by_lambda, refit = refit, rho_matrix = f$rho_matrix,
                       support_by_index = support_by_index),
       is_proxy = FALSE, refit = refit, n_per_slice = n_per_slice, N = N)
}
