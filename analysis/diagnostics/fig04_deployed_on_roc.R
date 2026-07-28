# ---------------------------------------------------------------------------
# analysis/scripts/fig04_deployed_on_roc.R
#
# PRODUCES  analysis/figures/fig04_deployed_on_roc.eps  (+ .png preview)
#
# tvcglasso vs CGLasso: the seed-averaged ROC curve of each, with the point
# each method's OWN native selector actually deploys marked on it, and the
# seed-to-seed spread of that point shown as a +/- 1 sd cross.
#
# This is the figure that separates the two things the comparison confounds:
#   the CURVE   = what the estimator can achieve if tuned correctly
#   the POINT   = where its own selector, given no truth, actually lands
#   the CROSS   = how reproducible that landing is across replicates
#
# HONEST CAVEAT ENCODED HERE: only tvcglasso's deployed point is a genuine
# index into its own curve (method_tv.R recomputes it on the pass-1 graph at
# the selected lambda). CGLasso deploys a PER-SLICE AMALGAM -- each time slice
# independently minimises its own BIC -- so its point need NOT lie on the
# pooled single-rho curve, and where it visibly sits off the curve that is a
# real property of the estimator, not a plotting error.
# ---------------------------------------------------------------------------

source(here::here("analysis", "config.R"))
for (f in list.files(here::here("analysis", "R"), full.names = TRUE, pattern = "\\.R$")) source(f)

REFRESH <- as.logical(Sys.getenv("TVCG_REFRESH", "FALSE"))
SHOW    <- c("tvcglasso", "CGLasso")

meta   <- cell_meta()
curves <- readRDS(file.path(CACHE_DIR, "fig01_curves.rds"))
pts    <- readRDS(file.path(CACHE_DIR, "fig02_points.rds"))

dep <- pts[pts$stage == "refit" & pts$method %in% SHOW, ]
summ <- do.call(rbind, lapply(split(dep, list(dep$method, dep$cell), drop = TRUE), function(g)
  data.frame(method = g$method[1], cell = g$cell[1], n = nrow(g),
             FPR = mean(g$FPR, na.rm = TRUE), TPR = mean(g$TPR, na.rm = TRUE),
             FPR_sd = stats::sd(g$FPR, na.rm = TRUE),
             TPR_sd = stats::sd(g$TPR, na.rm = TRUE),
             F1 = mean(g$F1, na.rm = TRUE), stringsAsFactors = FALSE)))

## panel order: rows = (m, depth), cols = (P, n)  -- same as Figure 1
meta$rowkey <- paste0(meta$m, "|", meta$depth)
row_levels <- unique(meta$rowkey[order(meta$m, meta$depth != "low")])
meta$colkey <- paste0(meta$P, "|", meta$n)
col_levels <- unique(meta$colkey[order(meta$P, meta$n)])
row_label <- function(k) { p <- strsplit(k, "\\|")[[1]]; paste0("m = ", p[1], ",  ", p[2], " depth") }
col_label <- function(k) { p <- strsplit(k, "\\|")[[1]]
  bquote(italic(P) == .(p[1]) * "," ~~ italic(n) == .(p[2])) }

draw_fig04 <- function() {
  nr <- length(row_levels); nc <- length(col_levels)
  mat <- rbind(matrix(seq_len(nr * nc), nrow = nr, byrow = TRUE), nr * nc + 1L)
  layout(mat, heights = c(rep(1, nr), 0.20))
  par(mar = c(1.1, 1.1, 1.3, 0.9), oma = c(0.6, 2.6, 0.2, 1.6),
      mgp = c(3, 0.25, 0), tcl = -0.20, xaxs = "i", yaxs = "i", cex = 1, las = 1)

  at <- c(0, 0.5, 1)
  for (ri in seq_len(nr)) for (ci in seq_len(nc)) {
    cc <- meta$cell[meta$rowkey == row_levels[ri] & meta$colkey == col_levels[ci]]
    plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
    segments(0, 0, 1, 1, lty = 3, lwd = LWD_REF, col = "grey55")

    for (mth in SHOW) {
      cv <- curves[[paste(mth, cc)]]$curve
      if (!is.null(cv)) lines(cv$FPR, cv$TPR, col = METHOD_COL[[mth]],
                              lty = METHOD_LTY[[mth]], lwd = LWD_CURVE)
    }
    ## deployed points last, so they sit on top of the curves
    for (mth in SHOW) {
      s <- summ[summ$method == mth & summ$cell == cc, ]
      if (!nrow(s)) next
      col <- METHOD_COL[[mth]]
      ## +/- 1 sd across seeds, in both coordinates
      segments(max(0, s$FPR - s$FPR_sd), s$TPR, min(1, s$FPR + s$FPR_sd), s$TPR,
               col = col, lwd = LWD_AXIS)
      segments(s$FPR, max(0, s$TPR - s$TPR_sd), s$FPR, min(1, s$TPR + s$TPR_sd),
               col = col, lwd = LWD_AXIS)
      points(s$FPR, s$TPR, pch = METHOD_PCH[[mth]], bg = col, col = "white",
             cex = 1.15, lwd = LWD_AXIS)
    }

    panel_axes(at, at, xlab_show = (ri == nr), ylab_show = (ci == 1),
               xfmt = c("0", "0.5", "1"), yfmt = c("0", "0.5", "1"))
    if (ri == 1) mtext(col_label(col_levels[ci]), side = 3, line = 0.25,
                       cex = par("cex") * CEX_STRIP)
    if (ci == nc) mtext(row_label(row_levels[ri]), side = 4, line = 0.35,
                        cex = par("cex") * CEX_STRIP, las = 0)
  }
  mtext("True positive rate", side = 2, outer = TRUE, line = 1.3,
        cex = par("cex") * CEX_LAB, las = 0)

  par(mar = c(0, 0, 0, 0)); plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0.5, 0.88, "False positive rate", adj = c(0.5, 1), cex = par("cex") * CEX_LAB, xpd = NA)
  legend(0.5, 0.30, xjust = 0.5, yjust = 0.5, ncol = 2, bty = "n",
         cex = CEX_LEGEND, seg.len = 2.4,
         legend = paste0(METHOD_LABEL[SHOW], "  (curve, deployed point +/- 1 sd)"),
         col = METHOD_COL[SHOW], lty = METHOD_LTY[SHOW], lwd = LWD_CURVE,
         pch = METHOD_PCH[SHOW], pt.bg = METHOD_COL[SHOW], pt.cex = 1.1, xpd = NA)
}

save_figure("fig04_deployed_on_roc", width = WIDTH_FULL, height = 7.0, draw = draw_fig04)

## the numbers behind the picture
summ <- merge(summ, meta[, c("cell", "P", "n", "m", "depth")], by = "cell")
summ <- summ[order(summ$cell, match(summ$method, METHODS)), ]
write.csv(summ, file.path(TAB_DIR, "tab04_deployed_on_roc.csv"), row.names = FALSE)
message("wrote ", file.path(TAB_DIR, "tab04_deployed_on_roc.csv"))
print(summ[, c("cell", "P", "n", "m", "depth", "method", "FPR", "FPR_sd", "TPR", "TPR_sd", "F1")],
      row.names = FALSE, digits = 3)
