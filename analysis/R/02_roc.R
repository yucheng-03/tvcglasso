# ---------------------------------------------------------------------------
# analysis/R/02_roc.R  --  ONE seed-averaging convention for ROC curves.
#
# The prototype scripts contained five different implementations with four
# different grids and the tie-envelope applied in only two of them. This file
# is the single convention; every figure uses it.
#
# The convention, and why each piece is there:
#
#  (1) TIES IN FPR ARE COLLAPSED TO THE UPPER ENVELOPE. Many lambda values give
#      FPR = 0 with different TPR (the empty graph, plus a vertical
#      zero-false-positive segment). An ROC read as a function of FPR is the
#      upper envelope of those points, so tapply(TPR, FPR, max).
#
#  (2) INTERPOLATION IS approx(rule = 1). Outside a seed's own [min FPR,
#      max FPR] the value is NA, never carried forward. rule = 2 would hold the
#      last TPR flat to FPR = 1 and draw a straight line into a region where
#      that seed has no operating point -- a fabricated corner. This is the
#      project's standing rule and it has published precedent: in Tian et al.'s
#      ROC_dense the curves visibly stop at their real maximum FPR (~0.80) and
#      leave the rest of the axis blank.
#
#  (3) THE MEAN IS DRAWN ONLY WHERE AT LEAST HALF THE SEEDS STILL HAVE REAL
#      POINTS. Beyond that the average would be over a shrinking, self-selected
#      subset (the seeds that happened to reach further), which biases the tail.
#
#  (4) (0,0) IS ANCHORED, (1,1) IS NOT. The empty graph is a REAL operating
#      point at the top of every method's data-adaptive lambda path, so (0,0)
#      is attained. The saturated graph generally is not: glasso on a
#      near-singular per-slice S at n < P stops densifying, and pushing the
#      path deeper does not help (measured: rho_lo_div 2000/5000/10000 all give
#      the same maxFPR at 1x/4x/10x the compute). That gap is a property of the
#      estimator, and it is disclosed in the caption, not filled in.
#
#  NOTE ON SCALARS: no AUC is computed or plotted here. The stored `auc` field
#  (R/roc_utils.R auc_trap) integrates over [0, maxFPR] with no (1,1) anchor,
#  so a method whose sweep stops at FPR 0.7 is not comparable with one that
#  reaches 1.0 -- quoting it would systematically penalise exactly the methods
#  with the honest coverage gap. The quantitative comparison is the as-deployed
#  table (analysis/tables/), not a number on the curve.
# ---------------------------------------------------------------------------

ROC_GRID     <- seq(0, 1, by = 0.0025)   # 401 points
ROC_MIN_FRAC <- 0.5                      # coverage rule (3)

#' Extract one seed's ROC as an FPR-indexed step function.
read_roc <- function(x, method = NULL) {
  r <- x$result
  roc <- r$roc
  if (is.null(roc) || !nrow(roc)) return(NULL)
  o <- order(roc$FPR, roc$TPR)
  list(FPR = roc$FPR[o], TPR = roc$TPR[o], maxFPR = max(roc$FPR),
       minFPR = min(roc$FPR))
}

#' Upper envelope over tied FPR values -- rule (1).
roc_envelope <- function(FPR, TPR) {
  a <- tapply(TPR, FPR, max)
  list(FPR = as.numeric(names(a)), TPR = as.numeric(a))
}

#' Seed-average a list of per-seed ROCs onto the common grid.
#'
#' Returns a data.frame(FPR, TPR, coverage, n) with TPR = NA wherever the
#' coverage rule bites. Plotting code must use type = "l" on the NA-containing
#' vector so the line simply stops -- do not na.omit() and reconnect.
roc_mean <- function(rocs, grid = ROC_GRID, min_frac = ROC_MIN_FRAC,
                     anchor_origin = TRUE) {
  rocs <- Filter(Negate(is.null), rocs)
  if (!length(rocs)) return(NULL)
  M <- vapply(rocs, function(s) {
    e <- roc_envelope(s$FPR, s$TPR)
    if (length(e$FPR) < 2) return(rep(NA_real_, length(grid)))
    stats::approx(e$FPR, e$TPR, xout = grid, rule = 1, ties = "ordered")$y
  }, numeric(length(grid)))
  if (is.null(dim(M))) M <- matrix(M, nrow = length(grid))
  cov <- rowMeans(!is.na(M))
  mu  <- rowMeans(M, na.rm = TRUE)
  mu[cov < min_frac] <- NA_real_
  out <- data.frame(FPR = grid, TPR = mu, coverage = cov, n = length(rocs))
  if (anchor_origin) {
    # (0,0) is a real attained operating point; make sure the curve starts there
    out$TPR[out$FPR == 0] <- 0
  }
  out
}

#' Per-(method, cell) corner coverage -- the honest companion to every ROC.
#'
#' Reports the median over seeds of the smallest and largest REAL FPR reached,
#' so the caption can state how much of the axis each curve actually spans.
roc_coverage <- function(rocs) {
  if (!length(rocs)) return(NULL)
  mn <- vapply(rocs, function(s) s$minFPR, numeric(1))
  mx <- vapply(rocs, function(s) s$maxFPR, numeric(1))
  data.frame(
    n           = length(rocs),
    minFPR_med  = stats::median(mn),
    maxFPR_med  = stats::median(mx),
    maxFPR_min  = min(mx),
    maxFPR_max  = max(mx),
    dense_gap   = 1 - stats::median(mx)
  )
}

#' Stitch a second, denser path segment onto a first one (CGLasso only).
#'
#' The main rho path and the dense-end extension share their seam rho exactly
#' but not their warm-start history, so a small step at the seam is real and is
#' drawn as it is. Nothing is interpolated, smoothed or joined across the seam.
roc_stitch <- function(roc_main, roc_ext) {
  if (is.null(roc_ext)) return(roc_main)
  if (is.null(roc_main)) return(roc_ext)
  FPR <- c(roc_main$FPR, roc_ext$FPR)
  TPR <- c(roc_main$TPR, roc_ext$TPR)
  o <- order(FPR, TPR)
  list(FPR = FPR[o], TPR = TPR[o],
       maxFPR = max(FPR), minFPR = min(FPR))
}
