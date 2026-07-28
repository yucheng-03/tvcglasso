# ============================================================================
# methods/stars_cglasso.R — StARS selection for the static CGLasso deployed point.
# Faithful to Yuan/Tian's Comp-gLASSO-JASA-main/补充材料/TARA_StARS.R rule:
#   - B = iter_sub subsamples (Yuan: 100), instability threshold beta (Yuan: 0.05);
#   - D_b(rho) = mean over edge-slots of 2*p*(1-p), p = selection frequency across subsamples;
#   - monotonized sup-instability; pick the DENSEST rho whose sup_D_b <= beta (sparse->dense scan, take last stable).
# ADAPTATIONS (documented for the supervisor): (1) subsample size n_sub = floor(frac*n) (Yuan's
#   10*sqrt(n) EXCEEDS n at our small n=12/20, so cannot be used); (2) per-slice time-varying:
#   CGLasso is fit independently per slice, so edge-slots are POOLED over the m slices' off-diagonals
#   (same pooling the lambda-sweep ROC uses), giving one instability per rho over all (i,j,k) slots;
#   (3) all subsamples use the SAME rho grid as the full fit so column r indexes the same rho.
# Requires fit_cglasso_static (code/baselines/cglasso_static.R) sourced.
# ============================================================================

stars_select_cglasso <- function(X_list, rho_grid, frac = 0.8, B = 100L, beta = 0.05,
                                  cglasso_args = list(), seed = NULL) {
  m <- length(X_list); n <- nrow(X_list[[1]]); P <- ncol(X_list[[1]]) - 1L
  n_sub <- max(3L, floor(frac * n))
  ut <- upper.tri(matrix(0, P, P))
  # one subsample -> (n_edge_slots x n_rho) 0/1 matrix: pooled off-diag edges over slices, per rho
  sub_edge_mat <- function(Om_by_rho) sapply(Om_by_rho, function(Oml)
    unlist(lapply(Oml, function(Mt) { d <- Mt; diag(d) <- 0; as.integer(d[ut] != 0) })))
  fit_args <- c(list(rho_list = rho_grid, verbose = FALSE), cglasso_args)

  subs <- lapply(seq_len(B), function(b) {
    if (!is.null(seed)) set.seed(as.integer(seed) * 100000L + b)   # reproducible, distinct per (seed,b)
    idx <- sample.int(n, n_sub)
    Xs  <- lapply(X_list, function(Xk) Xk[idx, , drop = FALSE])
    fb  <- tryCatch(do.call(fit_cglasso_static, c(list(X_list = Xs), fit_args)), error = function(e) NULL)
    if (is.null(fb)) NULL else sub_edge_mat(fb$Om_by_lambda)
  })
  subs <- Filter(Negate(is.null), subs); Bok <- length(subs)
  if (Bok < 2L) return(list(idx = NA_integer_, rho = NA_real_, Db = rep(NA_real_, length(rho_grid)),
                            Bok = Bok, n_sub = n_sub, note = "StARS failed: <2 usable subsamples"))
  nrho <- length(rho_grid)
  Db <- sapply(seq_len(nrho), function(r) {
    sel <- sapply(subs, function(s) s[, r]); p <- rowMeans(sel); mean(2 * p * (1 - p)) })   # instability per rho
  # Yuan rule: scan rho SPARSE -> DENSE (rho high -> low), monotone sup, keep last rho with sup <= beta
  ord <- order(rho_grid, decreasing = TRUE)
  supD <- 0; pick <- ord[1]
  for (r in ord) { supD <- max(supD, Db[r]); if (supD > beta) break; pick <- r }
  list(idx = pick, rho = rho_grid[pick], Db = Db, Bok = Bok, n_sub = n_sub, B_requested = B)
}

# deployed-point confusion vs truth, pooled off-diag over slices (matches the ROC edge rule)
.cg_confusion <- function(Om_list, true_list) {
  P <- nrow(Om_list[[1]]); ut <- upper.tri(matrix(0, P, P))
  pred <- unlist(lapply(Om_list,  function(M) { d <- M; diag(d) <- 0; d[ut] != 0 }))
  tru  <- unlist(lapply(true_list, function(M) { d <- M; diag(d) <- 0; d[ut] != 0 }))
  TP <- sum(pred & tru); FP <- sum(pred & !tru); FN <- sum(!pred & tru); TN <- sum(!pred & !tru)
  TPR <- TP / max(TP + FN, 1); FPR <- FP / max(FP + TN, 1); PR <- TP / max(TP + FP, 1)
  c(FPR = FPR, TPR = TPR, F1 = if (TPR + PR > 0) 2 * TPR * PR / (TPR + PR) else 0, edges = TP + FP)
}
