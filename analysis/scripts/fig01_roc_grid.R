# ---------------------------------------------------------------------------
# analysis/scripts/fig01_roc_grid.R
#
# PRODUCES  analysis/figures/fig01_roc_grid.eps   (paper Figure 1)
#           analysis/figures/fig01_roc_grid.png   (preview, not submitted)
#           analysis/tables/tab02_roc_coverage.csv (the honest companion)
#
# The seed-averaged ROC curve of all four methods, one panel per simulation
# cell: rows = (m, sequencing depth), columns = (P, n).
#
# Deliberately NOT on this figure:
#   * no AUC anywhere. The stored `auc` integrates over [0, maxFPR] only, so it
#     penalises exactly the methods whose curves honestly stop short. The
#     quantitative comparison is analysis/tables/tab01_deployed.csv.
#   * no deployed operating-point markers. In published JASA ROC figures the
#     curves carry no point symbols, and three of our four methods deploy a
#     point that does NOT lie on the plotted curve (per-slice, per-node and
#     per-lambda2 amalgams) -- drawing it here would imply otherwise.
#   * no interpolation into the region a method never reached. Each curve stops
#     at its real coverage; the gap is quantified in tab02 and stated in the
#     caption.
# ---------------------------------------------------------------------------

source(here::here("analysis", "config.R"))
for (f in list.files(here::here("analysis", "R"), full.names = TRUE, pattern = "\\.R$")) source(f)

REFRESH <- as.logical(Sys.getenv("TVCG_REFRESH", "FALSE"))

meta <- cell_meta()

## --- gather the seed-averaged curves ---------------------------------------

gather_curves <- function() {
  out <- list()
  for (mth in METHODS) {
    cells_m <- if (mth == "JGL") intersect(CELLS, JGL_CELLS) else CELLS
    fits <- list_fits(mth, cells = cells_m)
    require_complete(fits)
    for (cc in sort(unique(fits$cell))) {
      paths <- fits$path[fits$cell == cc]
      rocs <- lapply(paths, function(p) read_roc(read_fit(p), mth))

      ## CGLasso's rho path was run in two segments (main + dense-end
      ## extension). Both are real fitted points; concatenating them and
      ## sorting by FPR gives one curve with nothing interpolated across the
      ## seam. The small step at the seam is genuine -- the two segments reach
      ## the same rho with different warm-start histories -- and is drawn as is.
      if (mth == "CGLasso" && !is.na(CGLASSO_DENSEEND_DIR) &&
          dir.exists(CGLASSO_DENSEEND_DIR)) {
        ext_roots <- c(CGLasso = CGLASSO_DENSEEND_DIR)
        ext <- tryCatch(list_fits("CGLasso", cells = cc, roots = ext_roots),
                        error = function(e) NULL)
        if (!is.null(ext)) {
          key <- as.integer(sub(".*_seed([0-9]{3})\\.rds$", "\\1", basename(paths)))
          for (i in seq_along(rocs)) {
            j <- which(ext$seed == key[i])
            if (length(j) == 1L) {
              re <- tryCatch(read_roc(read_fit(ext$path[j]), mth),
                             error = function(e) NULL)
              rocs[[i]] <- roc_stitch(rocs[[i]], re)
            }
          }
        }
      }

      rocs <- Filter(Negate(is.null), rocs)
      out[[paste(mth, cc)]] <- list(
        method = mth, cell = cc,
        curve = roc_mean(rocs),
        coverage = roc_coverage(rocs)
      )
    }
  }
  out
}

curves <- with_cache("fig01_curves", gather_curves(), refresh = REFRESH)

## --- panel order ------------------------------------------------------------
## rows: (m, depth) with depth low above high; columns: (P, n).
meta$rowkey <- paste0(meta$m, "|", meta$depth)
row_levels <- unique(meta$rowkey[order(meta$m, meta$depth != "low")])
meta$colkey <- paste0(meta$P, "|", meta$n)
col_levels <- unique(meta$colkey[order(meta$P, meta$n)])

row_label <- function(k) {
  p <- strsplit(k, "\\|")[[1]]
  paste0("m = ", p[1], ",  ", p[2], " depth")
}
col_label <- function(k) {
  p <- strsplit(k, "\\|")[[1]]
  bquote(italic(P) == .(p[1]) * "," ~~ italic(n) == .(p[2]))
}

## --- draw -------------------------------------------------------------------

draw_fig01 <- function() {
  nr <- length(row_levels); nc <- length(col_levels)
  mat <- rbind(matrix(seq_len(nr * nc), nrow = nr, byrow = TRUE), nr * nc + 1L)
  layout(mat, heights = c(rep(1, nr), 0.34))
  ## NOTE: layout()/mfrow silently rescale par("cex") (0.66 at a 4x4 grid),
  ## which would drop 8 pt type to 5.3 pt without any warning. cex = 1 is set
  ## AFTER the layout call for exactly that reason.
  par(mar = c(1.1, 1.1, 1.3, 0.9), oma = c(0.6, 2.6, 0.2, 1.6),
      mgp = c(3, 0.25, 0), tcl = -0.20, xaxs = "i", yaxs = "i",
      cex = 1, las = 1)

  at <- c(0, 0.5, 1)
  for (ri in seq_len(nr)) {
    for (ci in seq_len(nc)) {
      cc <- meta$cell[meta$rowkey == row_levels[ri] & meta$colkey == col_levels[ci]]
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))

      ## chance diagonal: a substantive reference, drawn thin and dotted so it
      ## cannot be confused with a method curve
      segments(0, 0, 1, 1, lty = 3, lwd = LWD_REF, col = "grey55")

      for (mth in METHODS) {
        k <- paste(mth, cc)
        if (is.null(curves[[k]])) next
        cv <- curves[[k]]$curve
        if (is.null(cv)) next
        lines(cv$FPR, cv$TPR, col = METHOD_COL[[mth]], lty = METHOD_LTY[[mth]],
              lwd = LWD_CURVE)
      }

      panel_axes(at, at,
                 xlab_show = (ri == nr), ylab_show = (ci == 1),
                 xfmt = c("0", "0.5", "1"), yfmt = c("0", "0.5", "1"))
      if (ri == 1) mtext(col_label(col_levels[ci]), side = 3, line = 0.25,
                         cex = par("cex") * CEX_STRIP)
      if (ci == nc) mtext(row_label(row_levels[ri]), side = 4, line = 0.35,
                          cex = par("cex") * CEX_STRIP, las = 0)
    }
  }

  ## shared y-axis title, printed once
  mtext("True positive rate", side = 2, outer = TRUE, line = 1.3,
        cex = par("cex") * CEX_LAB, las = 0)

  ## The bottom strip carries the shared x-axis title and then the single
  ## legend for all panels.
  legend_strip(METHODS, type = "line", xlab = "False positive rate")
}

save_figure("fig01_roc_grid", width = WIDTH_FULL, height = 7.0, draw = draw_fig01)

## --- the coverage table that the caption must be able to cite ---------------
cov_rows <- lapply(curves, function(z) {
  if (is.null(z$coverage)) return(NULL)
  cbind(data.frame(method = z$method, cell = z$cell, stringsAsFactors = FALSE),
        z$coverage)
})
cov <- do.call(rbind, Filter(Negate(is.null), cov_rows))
cov <- merge(cov, meta[, c("cell", "P", "n", "m", "depth")], by = "cell")
cov <- cov[order(cov$cell, match(cov$method, METHODS)), ]
write.csv(cov, file.path(TAB_DIR, "tab02_roc_coverage.csv"), row.names = FALSE)
message("wrote ", file.path(TAB_DIR, "tab02_roc_coverage.csv"))
