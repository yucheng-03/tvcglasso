# ============================================================================
# code/baselines/cglasso_static.R — convenient MULTI-TIMEPOINT wrapper around the
# vendored canonical Compo_glasso (CGLasso, Tian et al. 2023).
#
# CGLasso is a STATIC single-network solver. As the "no-time-awareness static control"
# (CLAUDE.md §2), we fit it INDEPENDENTLY on each time slice and assemble the per-rho
# precision matrices across slices — exactly the comparison TV must beat by being
# time-aware. Returns the structure the lambda-sweep ROC machinery expects.
#
# RHO GRID = the paper's DATA-ADAPTIVE simulation-ROC grid (Comp-gLASSO-JASA-main/
# Simulation/Functions_ROC.R, random network):
#     Sigma = cov(ALR-z);  rho.list = exp(seq(log(max(Sigma)*2), log(max(Sigma)/1000), len=70))
# It auto-scales to the data so it spans dense -> EMPTY (covers the low-FPR corner) for
# any cell. We anchor on the LARGEST max(Sigma_k) across slices so ONE common grid spans
# empty->dense for every slice (needed for the common lambda-sweep ROC). length=70, the
# paper's resolution. Pass an explicit `rho_list` to override.
#
# fit_cglasso_static(X_list, length_rholist=70, rho_hi_mult=2, rho_lo_div=1000,
#                    option=2, max_iter=50, rho_list=NULL, verbose=FALSE)
#   X_list -> list (length m) of count matrices, each n_k x (P+1), reference taxon = last col.
#   -> list(Om_by_lambda, lambda)  with lambda ASCENDING and
#        Om_by_lambda : list (over ascending rho) of lists (over slices) of P x P precision.
#   (glasso exact-zero off-diagonals kept as 0 -> the lambda-sweep edge rule |Omega|!=0.)
# ============================================================================
source(here::here("R", "baselines", "CompoGlasso.R"))   # -> Compo_glasso + z_hat_offset/obj/NR

fit_cglasso_static <- function(X_list, length_rholist = 70L, rho_hi_mult = 2, rho_lo_div = 1000,
                               option = 2, max_iter = 50L, rho_list = NULL, verbose = FALSE,
                               ckpt_dir = NULL, progress_log = NULL, tol = 1e-8) {
  m <- length(X_list); K <- ncol(X_list[[1]]) - 1
  n_rho <- if (!is.null(rho_list)) length(rho_list) else length_rholist   # explicit grid may differ from length_rholist (codex minor-1)
  if (!is.null(ckpt_dir)) dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)
  .plog <- function(...) if (!is.null(progress_log)) cat(sprintf(...), file = progress_log, append = TRUE)

  # ★ PER-SLICE data-adaptive rho (2026-07-22): each time slice gets its OWN rho grid built
  # from ITS OWN max(cov(ALR-z)) -- the faithful reproduction of Yuan's data-adaptive recipe
  # (Functions_ROC.R) applied to each slice INDEPENDENTLY (= running the static CGLasso control
  # on that slice alone). No cross-slice coupling via a shared absolute grid; each slice spans
  # its own empty->dense range so the pooled INDEX ROC reaches both corners for every slice.
  # An explicit `rho_list` (one common grid) still overrides, for diagnostics only.
  per_slice_rho_desc <- if (!is.null(rho_list)) {
    replicate(m, sort(rho_list, decreasing = TRUE), simplify = FALSE)
  } else {
    lapply(seq_len(m), function(k) {
      zk <- z_hat_offset(as.matrix(X_list[[k]]), offset = K + 1, option = option)
      maxSig <- max(max(cov(zk)), .Machine$double.eps)
      sort(exp(seq(log(maxSig * rho_hi_mult), log(maxSig / rho_lo_div), length = n_rho)),
           decreasing = TRUE)   # sparse -> dense warm-start order (Compo_glasso warm-starts each rho from the previous)
    })
  }

  # per slice: one Compo_glasso call over that slice's OWN rho grid (warm-started sparse->dense).
  # SLICE-LEVEL checkpoint (resume-safe, preserves the warm-start chain bit-exact) + per-rho
  # progress redirected to progress_log (tail-able) instead of being swallowed.
  # ★ FINGERPRINT (codex P0-3, 2026-07-22): the checkpoint is bound to a content+config signature
  # (slice counts, rho grid endpoints/length, option, max_iter); a filename hit with a MISMATCHED
  # fingerprint (different data/config re-using the same tag) is treated as STALE -> re-fit, never
  # silently returned. Saved ATOMICally (tmp + rename) so a mid-write kill can't leave a torn file.
  .fp <- function(Xk, rho_desc_k) c(dim(Xk), sum(Xk), sum(Xk^2), sum(as.numeric(Xk) * seq_along(Xk)),
                                    option, max_iter, length(rho_desc_k),
                                    signif(rho_desc_k[1], 8), signif(rho_desc_k[length(rho_desc_k)], 8))
  # ★ PARALLEL (2026-07-22): slices are INDEPENDENT (own data, own checkpoint) -> fit them across
  # SLURM_CPUS_PER_TASK cores (like tvmgm's sweep). Each slice writes its OWN per-slice progress log
  # (progress_log.sliceKK) so parallel forks never share a file. Serial (1 core) if ncores==1.
  ncores <- { s <- Sys.getenv("SLURM_CPUS_PER_TASK"); if (nzchar(s)) max(1L, as.integer(s)) else 1L }
  .fit_slice <- function(k) {
    Xk <- as.matrix(X_list[[k]]); rho_desc_k <- per_slice_rho_desc[[k]]; fp_k <- .fp(Xk, rho_desc_k)
    plog_k <- if (!is.null(progress_log)) sprintf("%s.slice%02d", progress_log, k) else NULL
    plk <- function(...) if (!is.null(plog_k)) cat(sprintf(...), file = plog_k, append = TRUE)
    ckf <- if (!is.null(ckpt_dir)) file.path(ckpt_dir, sprintf("cgbase_slice%02d.rds", k)) else NULL
    if (!is.null(ckf) && file.exists(ckf)) {
      ck <- tryCatch(readRDS(ckf), error = function(e) NULL)
      if (is.list(ck) && !is.null(ck$fp) && identical(ck$fp, fp_k)) { plk("[cgbase] slice %d resumed (fingerprint ok)\n", k); return(ck$Oms) }
      plk("[cgbase] slice %d checkpoint STALE -> re-fitting\n", k)
    }
    call_cg <- function() Compo_glasso(Xk, rho.list = rho_desc_k, option = option, para_NR = FALSE, max_iter = max_iter)
    plk("[cgbase] slice %d start: %d rho in [%.4g, %.4g]\n", k, length(rho_desc_k), min(rho_desc_k), max(rho_desc_k))
    res <- withCallingHandlers(
      if (verbose) call_cg()
      else if (!is.null(plog_k)) { tmp <- NULL; capture.output(tmp <- call_cg(), file = plog_k, append = TRUE); tmp }
      else { tmp <- NULL; invisible(capture.output(tmp <- call_cg())); tmp },
      warning = function(w) if (grepl("not a multiple of split", conditionMessage(w))) invokeRestart("muffleWarning"))
    Oms <- res$Omegas.1
    if (dim(Oms)[3] != n_rho) stop(sprintf("slice %d: Compo_glasso returned %d rho, expected %d", k, dim(Oms)[3], n_rho))
    if (!is.null(ckf)) { tmp_ckf <- paste0(ckf, ".tmp"); saveRDS(list(Oms = Oms, fp = fp_k), tmp_ckf); file.rename(tmp_ckf, ckf) }
    plk("[cgbase] slice %d done\n", k); Oms
  }
  per_slice <- if (ncores > 1L) parallel::mclapply(seq_len(m), .fit_slice, mc.cores = min(ncores, m), mc.preschedule = FALSE)
               else lapply(seq_len(m), .fit_slice)
  .e <- which(vapply(per_slice, function(x) inherits(x, "try-error"), logical(1)))   # mclapply returns try-error, not stop
  if (length(.e)) stop("cgbase slice ", .e[1], " failed: ", attr(per_slice[[.e[1]]], "condition")$message)

  # reorganize to ASCENDING-INDEX list: index j = each slice's j-th ASCENDING rho (a common
  # RELATIVE regularization level; the ABSOLUTE rho differs per slice). j-th ascending == the
  # (n_rho-j+1)-th of the descending path (exact index, no float match).
  Om_by_lambda <- lapply(seq_len(n_rho), function(j) {
    idx <- n_rho - j + 1L
    lapply(seq_len(m), function(k) {
      Mk <- per_slice[[k]][, , idx]
      Mk <- (Mk + t(Mk)) / 2                 # ONE canonical SYMMETRIC Omega (codex P0-2): huge's BCD can
      Mk[abs(Mk) <= tol] <- 0                # return a mildly asymmetric matrix; symmetrize BEFORE thresholding
      Mk                                     # so base ROC/support AND the refit (which also symmetrizes) agree.
    })
  })
  rho_asc_mat <- vapply(per_slice_rho_desc, rev, numeric(n_rho))   # n_rho x m, ascending per slice
  list(Om_by_lambda = Om_by_lambda, lambda = seq_len(n_rho),
       rho_matrix = rho_asc_mat, per_slice_rho_desc = per_slice_rho_desc)
}
