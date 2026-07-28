# ============================================================================
# methods/_roc_utils.R — shared ROC machinery + the per-method OUTPUT CONTRACT.
#
# Every method file (method_tv.R, method_cglasso.R, method_tvmgm.R, method_lupine.R)
# exposes ONE function  run_method_X(dat, cfg, eval_slices)  that returns:
#
#   list(
#     method     = "TV" | "CGLasso" | "tvmgm" | "LUPINE_single" | "LUPINE_long",
#     convention = "lambda-sweep" | "score-sweep",
#     roc        = data.frame(FPR, TPR)  ordered by FPR  (the method's NATIVE ROC),
#     auc        = <numeric, trapezoid AUC with (0,0)/(1,1) anchors>,
#     scores     = list length m of P x P score matrices (higher = more edge-like),
#                  OR NULL. Used by the MASTER for a future METHOD-SYMMETRIC score-sweep.
#                  For lambda-sweep methods this is the per-edge ENTRY-LAMBDA (largest
#                  lambda at which the edge survives = regularization-path entry, the
#                  glmnet-style continuous score) -> all 4 methods become score-comparable.
#     detail     = <method-specific extras: lambda grid, n_edges/lambda, bw_sel, ...>
#   )
#
#   dat         = output of generate_simulation_data_mirrorexp() (X, true_Omega_list, x_sequence)
#   cfg         = list of run parameters (lambda grids, q, N_n, ncomp, bwSeq, ...)
#   eval_slices = integer vector of time slices to evaluate on (e.g. 2:m when LUPINE_long
#                 is in the comparison, since it cannot estimate slice 1) — applied to EVERY
#                 method so the comparison is slice-consistent.
# ============================================================================

# null-coalescing helper so cfg can omit optional keys (shared by all method files).
`%||%` <- function(a, b) if (is.null(a)) b else a

# ★ P1-7 (2026-07-22): Trapezoid AUC over the ROC's REAL operating points. Anchors (0,0)
# (the empty graph IS a real, reachable operating point at high lambda -> data-adaptive
# lambda_max guarantees it) but does NOT fabricate the (1,1) corner: a thresholded sparse
# estimator (TV) caps its density and may never reach FPR=1; extrapolating to (1,1) invents
# area the data never demonstrated (criticized -- Muschelli 2019; RJafroc). The REPORTED
# comparison is the as-deployed tuning table + the ROC curves (Xue-Shu-Qu/SPACE/CGLasso
# practice; the sparse-GGM field does NOT report a full-range scalar AUC). This scalar is
# for INTERNAL ranking (e.g. JGL lambda2 / mgm bw) + diagnostics, over the covered FPR span.
auc_trap <- function(df) {
  d <- df[order(df$FPR), c("FPR", "TPR")]
  if (!any(d$FPR == 0 & d$TPR == 0)) d <- rbind(data.frame(FPR = 0, TPR = 0), d)
  sum(diff(d$FPR) * (head(d$TPR, -1) + tail(d$TPR, -1)) / 2)
}

# ★ P1-7: partial AUC over [0, fmax] on REAL points, standardized to [0,1] by /fmax. The pAUC
# CAPABILITY is kept in the code but is NOT reported by default (no committed story yet).
pauc_trap <- function(df, fmax = 0.2) {
  d <- df[order(df$FPR), c("FPR", "TPR")]
  if (!any(d$FPR == 0 & d$TPR == 0)) d <- rbind(data.frame(FPR = 0, TPR = 0), d)
  d <- d[d$FPR <= fmax + 1e-12, , drop = FALSE]
  if (nrow(d) < 2) return(NA_real_)
  if (max(d$FPR) < fmax) d <- rbind(d, data.frame(FPR = fmax, TPR = tail(d$TPR, 1)))
  sum(diff(d$FPR) * (head(d$TPR, -1) + tail(d$TPR, -1)) / 2) / fmax
}

# --- ROC SELF-CHECK (always run after producing an ROC) -----------------------
# ★ P1-7: roc_anchor ensures the curve STARTS at (0,0) (the empty graph = a real high-lambda
# operating point) but does NOT add (1,1) -- a sparse estimator may never reach the dense
# corner, so we plot/integrate over REAL coverage only (no fabricated top-right corner).
roc_anchor <- function(roc) {
  d <- unique(roc[order(roc$FPR, roc$TPR), c("FPR", "TPR")])
  if (!any(d$FPR == 0 & d$TPR == 0)) d <- rbind(data.frame(FPR = 0, TPR = 0), d)
  d[order(d$FPR, d$TPR), ]
}

# roc_sanity_report: flag when a method's REAL data does not reach the corners, so the
# AUC over the un-reached span is just straight-line interpolation (e.g. a too-narrow
# lambda grid). dense_gap = 1 - maxFPR (interpolated region at the dense end);
# sparse_gap = minFPR. A large gap (>0.1) means "extend the grid before trusting the AUC".
roc_sanity_report <- function(roc, method = "") {
  maxF <- max(roc$FPR); minF <- min(roc$FPR); maxT <- max(roc$TPR); minT <- min(roc$TPR)
  dgap <- 1 - maxF; sgap <- minF
  flag <- character()
  if (dgap > 0.10) flag <- c(flag, sprintf("DENSE-END GAP: data stops at FPR=%.2f -> AUC interpolates top %.0f%% (extend dense end of grid)", maxF, 100 * dgap))
  if (sgap > 0.10) flag <- c(flag, sprintf("SPARSE-END GAP: data starts at FPR=%.2f (extend sparse end)", minF))
  data.frame(method = method, minFPR = round(minF, 3), maxFPR = round(maxF, 3),
             maxTPR = round(maxT, 3), dense_gap = round(dgap, 3),
             ok = length(flag) == 0, flag = paste(flag, collapse = "; "))
}

# SCORE-SWEEP ROC: pool off-diagonal scores over eval_slices, threshold a grid. (mgm/LUPINE)
roc_from_scores <- function(score_list, true_list, eval_slices, n_thr = 200) {
  P <- nrow(true_list[[1]]); ut <- upper.tri(matrix(0, P, P))
  scores <- numeric(); labels <- logical()
  for (k in eval_slices) {
    scores <- c(scores, abs(score_list[[k]][ut]))
    labels <- c(labels, true_list[[k]][ut] != 0)
  }
  thr <- sort(unique(c(0, quantile(scores, seq(0, 1, length.out = n_thr), na.rm = TRUE),
                       max(scores, na.rm = TRUE) * 1.001)), decreasing = TRUE)
  d <- do.call(rbind, lapply(thr, function(t) {
    pred <- scores >= t
    data.frame(FPR = sum(pred & !labels) / max(1, sum(!labels)),
               TPR = sum(pred &  labels) / max(1, sum( labels)))
  }))
  d <- d[order(d$FPR), ]
  list(roc = d, auc = auc_trap(d))
}

# One LAMBDA-SWEEP operating point: edge = (Omega_hat != 0), pooled over eval_slices.
# Matches metrics_05-20's pred_adj <- (Omega_hat != 0) convention. (TV/CGLasso)
edge_fpr_tpr <- function(Om_list, true_list, eval_slices) {
  P <- nrow(Om_list[[1]]); ut <- upper.tri(matrix(0, P, P))
  tp <- fp <- tpos <- fpos <- 0
  for (k in eval_slices) {
    pe <- Om_list[[k]][ut] != 0
    te <- true_list[[k]][ut] != 0
    tp <- tp + sum(pe & te); fp <- fp + sum(pe & !te)
    tpos <- tpos + sum(te);  fpos <- fpos + sum(!te)
  }
  c(FPR = fp / max(1, fpos), TPR = tp / max(1, tpos))
}

# Build the LAMBDA-SWEEP ROC from a list (ascending lambda) of assembled Omega(t) lists.
roc_from_lambda_path <- function(Om_by_lambda, lambda_grid, true_list, eval_slices) {
  pts <- do.call(rbind, lapply(seq_along(Om_by_lambda), function(j) {
    ft <- edge_fpr_tpr(Om_by_lambda[[j]], true_list, eval_slices)
    data.frame(FPR = ft["FPR"], TPR = ft["TPR"], lambda = lambda_grid[j])
  }))
  pts <- pts[order(pts$FPR), ]; row.names(pts) <- NULL
  list(roc = pts[, c("FPR", "TPR")], auc = auc_trap(pts), detail = pts)
}

# ENTRY-LAMBDA per-edge score (for the symmetric score-sweep): for each slice & edge,
# the largest lambda at which |Omega_hat_ij(t_k)| != 0. Edges never selected -> 0.
# Om_by_lambda must be ASCENDING in lambda. Higher score = survives more shrinkage = stronger.
entry_lambda_scores <- function(Om_by_lambda, lambda_grid, m) {
  P <- nrow(Om_by_lambda[[1]][[1]])
  ord <- order(lambda_grid)                       # ensure ascending
  Om_by_lambda <- Om_by_lambda[ord]; lambda_grid <- lambda_grid[ord]
  lapply(1:m, function(k) {
    S <- matrix(0, P, P)
    for (j in seq_along(lambda_grid)) {
      alive <- Om_by_lambda[[j]][[k]] != 0
      S[alive] <- lambda_grid[j]                  # overwrite -> ends at the LARGEST surviving lambda
    }
    diag(S) <- 0; S
  })
}
