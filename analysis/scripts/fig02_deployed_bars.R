# ---------------------------------------------------------------------------
# analysis/scripts/fig02_deployed_bars.R
#
# PRODUCES  analysis/figures/fig02_deployed.eps  (paper Figure 2)
#           analysis/figures/fig02_deployed.png  (preview, not submitted)
#
# The AS-DEPLOYED operating point: recall, precision and F1 achieved when each
# method is tuned by its OWN native selector, with no knowledge of the truth.
#   tvcglasso  refit-BIC (joint LNM at the inferred latent Zhat)
#   CGLasso    per-slice BIC
#   JGL        refit-AIC (Danaher et al.'s recommended selector)
#   tvmgm      per-node relaxed EBIC (mgm's paper-native selector)
#
# The selectors are deliberately NOT unified. Uniform-refit fairness in this
# paper means the same relaxed-lasso RECIPE (select -> unpenalised refit ->
# reselect) applied inside each method's own estimation framework with its own
# native selector; it does not mean one shared information criterion. The four
# criteria are not on a common scale and are never compared numerically -- what
# is compared is the deployed GRAPH.
#
# Greyscale safety: T&F print figures in black and white and forbid encoding
# series identity by hue alone. Bars cannot carry a line type, so the four
# fills are separated in LUMINANCE (CIE L* = 0.0 / 29.2 / 49.5 / 72.5) and
# additionally hatched at alternating angles.
# ---------------------------------------------------------------------------

source(here::here("analysis", "config.R"))
for (f in list.files(here::here("analysis", "R"), full.names = TRUE, pattern = "\\.R$")) source(f)

REFRESH <- as.logical(Sys.getenv("TVCG_REFRESH", "FALSE"))
STAGE   <- Sys.getenv("TVCG_STAGE", "refit")   # "refit" (deployed) or "pre"

meta <- cell_meta()

gather_points <- function() {
  rows <- list()
  for (mth in METHODS) {
    cells_m <- if (mth == "JGL") intersect(CELLS, JGL_CELLS) else CELLS
    p <- collect_points(mth, stages = c("refit", "pre"), cells = cells_m)
    if (!is.null(p)) rows[[mth]] <- p
  }
  do.call(rbind, rows)
}

pts  <- with_cache("fig02_points", gather_points(), refresh = REFRESH)
summ <- summarise_points(pts)
write.csv(summ, file.path(TAB_DIR, "tab01_deployed.csv"), row.names = FALSE)
message("wrote ", file.path(TAB_DIR, "tab01_deployed.csv"))

S <- summ[summ$stage == STAGE, , drop = FALSE]

METRICS <- c(TPR = "Recall", precision = "Precision", F1 = "F1")

## Column groups: one bar cluster per (P, n); rows: sequencing depth.
## m = 7 and m = 15 are shown as separate figures because JGL exists only at
## m <= 7 (its ADMM does not scale) -- putting them in one panel would render
## an absent method as a missing bar.
M_TARGET <- as.integer(Sys.getenv("TVCG_M", "7"))
S <- S[S$m == M_TARGET, , drop = FALSE]
meta_m <- meta[meta$m == M_TARGET, , drop = FALSE]

depths  <- c("low", "high")
settings <- unique(meta_m[order(meta_m$P, meta_m$n), c("P", "n")])
set_lab  <- sprintf("%d, %d", settings$P, settings$n)

draw_fig02 <- function() {
  nr <- length(depths); nc <- length(METRICS)
  mat <- rbind(matrix(seq_len(nr * nc), nrow = nr, byrow = TRUE), nr * nc + 1L)
  layout(mat, heights = c(rep(1, nr), 0.34))
  ## cex = 1 AFTER layout(): a panel layout silently rescales par("cex").
  par(mar = c(1.5, 1.4, 1.4, 0.6), oma = c(0.4, 2.2, 0.2, 1.8),
      mgp = c(3, 0.25, 0), tcl = -0.20, cex = 1, las = 1)

  gap_in  <- 0.35   # gap between clusters, in bar widths
  for (ri in seq_along(depths)) {
    for (ci in seq_along(METRICS)) {
      mt <- names(METRICS)[ci]
      plot.new()
      nset <- nrow(settings)
      plot.window(xlim = c(0.4, nset + 0.6), ylim = c(0, 1))

      ## faint horizontal guides only: the one place gridlines are defensible
      ## (reading bar heights), matching the published bar figure in Tian et al.
      abline(h = seq(0.25, 1, by = 0.25), col = "grey88", lwd = LWD_REF)

      meths <- METHODS[METHODS %in% unique(S$method)]
      k  <- length(meths)
      bw <- (1 - gap_in) / k
      for (si in seq_len(nset)) {
        cc <- meta_m$cell[meta_m$P == settings$P[si] & meta_m$n == settings$n[si] &
                            meta_m$depth == depths[ri]]
        for (mi in seq_along(meths)) {
          row <- S[S$cell == cc & S$method == meths[mi], , drop = FALSE]
          if (!nrow(row)) next
          mu <- row[[paste0(mt, "_mean")]]
          ## SPREAD = one standard deviation across seeds, lower end floored at
          ## zero -- the convention of the published JASA bar figure in Tian et
          ## al. It answers "how much does a single replicate vary", which is
          ## what a reader wants from a simulation study. The standard error
          ## (sd/sqrt(100)) is ten times smaller and would say only that the
          ## mean is precisely estimated. Both are in tab01_deployed.csv; the
          ## caption MUST state which one the bars show.
          se <- row[[paste0(mt, "_sd")]]
          xl <- si - 0.5 + gap_in / 2 + (mi - 1) * bw
          rect(xl, 0, xl + bw * 0.92, mu,
               col = METHOD_COL[[meths[mi]]], border = "black", lwd = LWD_AXIS)
          error_bar(x = xl + bw * 0.46, mean = mu, se = se,
                    fill = METHOD_COL[[meths[mi]]], half_width = bw * 0.18)
        }
      }
      axis(2, at = seq(0, 1, 0.25), labels = if (ci == 1) c("0", "0.25", "0.5", "0.75", "1") else FALSE,
           lwd = 0, lwd.ticks = LWD_AXIS, tcl = -0.20, cex.axis = CEX_AXIS,
           mgp = c(3, 0.45, 0))
      axis(1, at = seq_len(nset), labels = if (ri == nr) set_lab else FALSE,
           lwd = 0, lwd.ticks = 0, cex.axis = CEX_AXIS, mgp = c(3, 0.15, 0))
      usr <- par("usr")
      segments(usr[1], 0, usr[2], 0, lwd = LWD_AXIS, xpd = NA)
      segments(usr[1], 0, usr[1], 1, lwd = LWD_AXIS, xpd = NA)

      if (ri == 1) mtext(METRICS[[ci]], side = 3, line = 0.3,
                         cex = par("cex") * CEX_STRIP, font = 2)
      if (ci == nc) mtext(paste0(depths[ri], " depth"), side = 4, line = 0.35,
                          cex = par("cex") * CEX_STRIP, las = 0)
    }
  }
  legend_strip(METHODS[METHODS %in% unique(S$method)], type = "fill",
               xlab = expression(italic(P) * "," ~ italic(n)))
}

save_figure(sprintf("fig02_deployed_m%02d", M_TARGET),
            width = WIDTH_FULL, height = 3.5, draw = draw_fig02)
