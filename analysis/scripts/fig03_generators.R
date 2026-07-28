# ---------------------------------------------------------------------------
# analysis/scripts/fig03_generators.R
#
# PRODUCES  analysis/figures/fig03_edge_shapes.eps  (paper Figure 3)
#           analysis/figures/fig03_edge_shapes.png  (preview, not submitted)
#
# The two true edge-trajectory shapes of the simulation generator
# (mirrorexp_contB, the only generator the publication grid uses):
#
#   (a) EARLY   f_early(t) = A (e^{-rt} - c) / (1 - c)  on [0, 1/2],  0 after
#   (b) LATE    f_late(t)  = f_early(1 - t)             -- the mirror image
#
# with c = e^{-r/2}, A = 0.5, r = 6. Every active edge draws one of the two with
# probability 1/2. Both are continuous at t = 1/2, where they reach EXACTLY
# zero, so each edge has a genuine temporal zero region on half the interval --
# the property the zero-identification results are about. The curve is C0 but
# not C1 there: it arrives at zero with non-zero slope and then stops.
#
# The two lines are defined here because they are local variables inside the
# generator and cannot be reached from outside; they are then CHECKED against
# the generator's actual true_Omega_list, with a hard stop on any mismatch, so
# the figure cannot drift away from the estimator it illustrates.
# ---------------------------------------------------------------------------

source(here::here("analysis", "config.R"))
for (f in list.files(here::here("analysis", "R"), full.names = TRUE, pattern = "\\.R$")) source(f)
source(here::here("simulation", "generators.R"))

## --- the publication parameters (config/paired_grid.R) ----------------------
A      <- 0.5   # edge_strength
RATE   <- 6     # decay rate
M_SHOW <- 7     # observation grid drawn as points

## Mirrors simulation/generators.R:48-50 exactly; verified below.
c0 <- exp(-RATE / 2); Zn <- 1 - c0
f_early <- function(t) ifelse(t <= 0.5, A * (exp(-RATE * t)       - c0) / Zn, 0)
f_late  <- function(t) ifelse(t >= 0.5, A * (exp(-RATE * (1 - t)) - c0) / Zn, 0)

## --- verify against the shipped generator, do not trust the transcription ---
verify_shapes <- function(P = 15, n = 12, m = M_SHOW, seed = 1) {
  cell <- list(generator = "mirrorexp_contB", P = P, n = n, m = m,
               depth_mode = "high", depth_lo = 900, depth_hi = 1800,
               edge_strength = A, rate = RATE, base_diag_mean = 0.1)
  d <- generate_data(cell, seed)
  ts <- seq(0, 1, length.out = m)
  worst <- 0
  for (k in seq_len(m)) {
    Om  <- d$true_Omega_list[[k]]
    off <- Om[upper.tri(Om)]; off <- off[off != 0]
    if (!length(off)) next
    allowed <- c(f_early(ts[k]), f_late(ts[k]))
    worst <- max(worst, max(vapply(off, function(v) min(abs(v - allowed)), numeric(1))))
  }
  worst
}
dev <- verify_shapes()
if (!is.finite(dev) || dev > 1e-10) {
  stop("fig03: f(t) no longer matches simulation/generators.R (max deviation ",
       format(dev), "). Fix the figure, not the check.")
}
message(sprintf("fig03: f(t) matches the generator (max dev %.2e)", dev))

## --- appearance -------------------------------------------------------------
## A solid light fill under the active half makes the support visible at a
## glance and contrasts it with the flat zero half. Solid, never alpha: EPS
## carries no transparency and T&F require flattened artwork.
FILL_ACTIVE <- "#DCE3EA"   # light neutral blue-grey
COL_ZERO    <- "#6E7A88"   # only the "exactly zero" caption, not the curve

tt  <- seq(0, 1, length.out = 2001)
obs <- seq(0, 1, length.out = M_SHOW)

draw_panel <- function(f, strip, active_side = c("left", "right")) {
  active_side <- match.arg(active_side)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(-0.015, 0.56))

  ## shaded support: the half where the edge is active
  xs <- if (active_side == "left") tt[tt <= 0.5] else tt[tt >= 0.5]
  polygon(c(xs[1], xs, xs[length(xs)]), c(0, f(xs), 0),
          col = FILL_ACTIVE, border = NA)

  ## t = 1/2 guide
  segments(0.5, 0, 0.5, 0.56, lty = 3, lwd = LWD_REF, col = "grey60")

  ## The trajectory, drawn as ONE black curve across the whole interval: the
  ## zero half is part of the same function, not a separate object. What marks
  ## it out is the absence of shading plus the label -- colouring it differently
  ## would invite the reader to treat it as a different kind of thing.
  lines(tt, f(tt), lwd = LWD_CURVE * 1.15, col = "black")

  ## the m observation times actually available to the estimator
  points(obs, f(obs), pch = 21, bg = "white", col = "black",
         cex = 0.62, lwd = LWD_REF)

  ## in-panel labels: with one trajectory per panel there is room for these,
  ## and they carry the point of the figure without a legend
  lab_x_active <- if (active_side == "left") 0.20 else 0.80
  lab_x_zero   <- if (active_side == "left") 0.76 else 0.24
  text(lab_x_active, 0.50, "active", cex = CEX_STRIP * 0.95, adj = c(0.5, 0.5))
  text(lab_x_zero,   0.10, "exactly zero", cex = CEX_STRIP * 0.95,
       adj = c(0.5, 0.5), col = COL_ZERO)

  panel_axes(c(0, 0.25, 0.5, 0.75, 1), c(0, 0.25, 0.5),
             xfmt = c("0", "0.25", "0.5", "0.75", "1"),
             yfmt = c("0", "0.25", "0.5"))
  panel_strip(strip)
  mtext("Time t", side = 1, line = 1.25, cex = par("cex") * CEX_LAB)
}

draw_fig03 <- function() {
  par(mfrow = c(1, 2))
  ## cex = 1 AFTER the layout: a panel layout silently rescales par("cex").
  par(mar = c(2.4, 3.7, 1.6, 0.8), oma = c(0.2, 0.2, 0.2, 0.2),
      mgp = c(3, 0.25, 0), tcl = -0.20, cex = 1, las = 1, xaxs = "i", yaxs = "i")

  draw_panel(f_early, "(a)  early", active_side = "left")
  mtext("Edge weight", side = 2, line = 2.6, cex = par("cex") * CEX_LAB, las = 0)
  draw_panel(f_late,  "(b)  late",  active_side = "right")
}

save_figure("fig03_edge_shapes", width = WIDTH_FULL, height = 2.9, draw = draw_fig03)
