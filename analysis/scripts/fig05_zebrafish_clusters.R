# ---------------------------------------------------------------------------
# analysis/scripts/fig05_zebrafish_clusters.R
#
# PRODUCES  analysis/figures/fig05_zeb_clusters_P<P>.eps  (+ .png preview)
#           analysis/tables/tab05_zeb_clusters_P<P>.csv
#
# Zebrafish real data: the estimated edge trajectories Omega_ij(t) grouped by
# k-means on their SHAPE, separately for the infected and uninfected fish.
#
#   rows     the two fish groups
#   columns  the clusters, ordered by when the cluster mean peaks
#   grey     the individual member trajectories
#   black    the cluster mean, with the m = 7 sampling days marked
#
# WHAT IS CLUSTERED, AND WHAT IS DRAWN, ARE DIFFERENT OBJECTS. k-means runs on
# z-scored trajectories, so edges are grouped by the SHAPE of their time course
# and not by how strong they are; the bold line is then the RAW (un-z-scored)
# mean of that cluster's members, so the vertical axis stays interpretable as
# a precision-matrix entry. Keep both facts in the caption -- a reader who
# assumes the clustering used the raw values will misread the panel heights.
#
# ENCODING: grey members + one black mean, and cluster identity carried by
# POSITION rather than colour. With up to six clusters no colour scheme stays
# separable in greyscale, and these are not the four methods, so reusing the
# method palette would be actively misleading.
#
# ---------------------------------------------------------------------------
# ⚠ TWO CAVEATS -- BOTH BELONG IN THE CAPTION.
#
# (1) THE OPERATING POINT. On this data the refit reaches its convergence gate
# at almost no penalty value: across the four fits (P in {15,25} x two groups,
# 50 lambda each) only 21 of 200 lambda converged, and 20 of those are the EMPTY
# graph. Three of the four fits therefore have NO eligible lambda and deploy via
# the all-lambda fallback; the fourth (P=15 infected) has exactly one non-empty
# eligible lambda and deploys a ONE-EDGE graph, which cannot be clustered at all.
# Using each fit's stated deployment would compare a 1-edge network against a
# 51-edge one and read as a dramatic infection effect that is purely an artefact
# of which fits happened to converge: under the uniform fallback the two groups
# have 49 vs 51 edges at P=15 and 91 vs 91 at P=25. This script therefore
# defaults to the UNIFORM path-wide BIC minimum -- which is exactly what three
# of the four fits already deploy -- and reports the convergence status of the
# point it used. Set TVCG_ZEB_POINT=deployed to use each fit's stated rule.
#
# Consequence for reading the panels: the SHAPES are trustworthy, the HEIGHTS
# are not. Solving the same beta problem to stationarity (frozen Z, damped
# Newton, ~2 s) leaves the support identical and the z-scored shapes almost
# unchanged (median correlation 0.988, no sign inversions, 84% of edges keep
# their peak day) but changes max|Omega_ij| by a factor of ~4. The clustering
# survives the non-convergence; the vertical axis does not.
#
# (2) THE BIOLOGY. The phase-composition contrast (infected fish showing early-
# and late-peaking edges with few mid-peaking ones; uninfected fish showing more
# mid-peaking edges) survived a matched-sample-size check, but is NOT
# established as biology. Every mid-transient edge in the uninfected group
# touches a taxon that is absent on some days, whose ALR value is then pinned by
# the pseudocount floor -- the same near-absent-taxon artefact that produced a
# spuriously "stable" edge elsewhere in this data set. A within-group bootstrap
# calls those edges reproducible, which is exactly what a deterministic artefact
# would do. Measured directly: the largest panel with no all-zero (taxon, day)
# cell in either group is P = 12; P = 15 has 2 and 4 such cells and P = 25 has
# 17 and 31. Treat the figure as a description of the fitted networks, not as
# evidence of an infection effect.
# ---------------------------------------------------------------------------

source(here::here("analysis", "config.R"))
for (f in list.files(here::here("analysis", "R"), full.names = TRUE, pattern = "\\.R$")) source(f)
suppressPackageStartupMessages(library(cluster))   # Omega comes straight from the stored fit; no engine needed
`%||%` <- function(a, b) if (is.null(a)) b else a

P_TARGET <- as.integer(Sys.getenv("TVCG_ZEB_P", "15"))
POINT    <- match.arg(Sys.getenv("TVCG_ZEB_POINT", "bicmin"), c("bicmin", "deployed"))
EDGE_TOL <- 0.01     # an edge is "active" if it exceeds this at some time point
SEED     <- 42       # k-means is randomised; fixed so the figure is reproducible
GROUPS   <- c("infected", "not_infected")   # one fit file per group
GROUP_LAB <- c(infected = "Infected", not_infected = "Uninfected")

read_fit <- function(group) {
  f <- file.path(ZEB_FITS_DIR, sprintf("zeb_fit_P%02d_%s.rds", P_TARGET, group))
  if (!file.exists(f)) {
    stop("Zebrafish fit not found: ", f,
         "\n  Produce it with analysis/realdata/fit_zebrafish.R, or set",
         "\n  ZEB_FITS_DIR in analysis/config.R to where the fits live.", call. = FALSE)
  }
  readRDS(f)
}
GROUPS_KEY <- GROUPS
fits <- lapply(GROUPS_KEY, read_fit); names(fits) <- GROUPS_KEY
days <- fits[[1]]$days; m <- length(days); nodes <- fits[[1]]$node_names

#' Active edge trajectories at the deployed operating point.
#'
#' TWO PASSES, TWO ROLES -- and the figure needs the second one. The SUPPORT
#' (which edges are non-zero) and the selected lambda come from pass 1; the
#' MAGNITUDES come from the relaxed refit, which re-estimates the surviving
#' entries without the L1 penalty. A trajectory plot is a plot of magnitudes,
#' so it must read Omega_post, not the penalized pass-1 Omega -- otherwise
#' every curve is shrunk toward zero by exactly the penalty the refit exists to
#' remove, and the vertical axis understates the estimated dependence.
#'
#' The operating point is the method's own native selector (refit-BIC) applied
#' to real data, where there is no truth to tune against. POINT = "deployed"
#' uses each fit's stated rule (BIC minimum among CONVERGED lambda, with an
#' all-lambda fallback); POINT = "bicmin" applies the all-lambda BIC minimum to
#' every fit, which is what three of the four fits already deploy and is the
#' only choice under which all four are drawable and mutually comparable. See
#' caveat (1) in the header.
collect_traj <- function(fit) {
  di <- if (POINT == "deployed") fit$deployed$lambda_index else which.min(fit$detail$BIC)
  rf <- fit$estimate$refit[[di]]
  if (is.null(rf$Omega_post))
    stop("lambda[", di, "] has no refit for ", fit$group,
         " (error: ", rf$error %||% "none", ")", call. = FALSE)
  Om <- rf$Omega_post
  ut <- which(upper.tri(matrix(0, length(nodes), length(nodes))), arr.ind = TRUE)
  keep <- list(); ij <- list()
  for (e in seq_len(nrow(ut))) {
    i <- ut[e, 1]; j <- ut[e, 2]
    tr <- vapply(seq_len(m), function(k) Om[[k]][i, j], numeric(1))
    if (max(abs(tr)) <= EDGE_TOL) next
    keep[[length(keep) + 1L]] <- tr; ij[[length(ij) + 1L]] <- c(i, j)
  }
  if (!length(keep))
    stop("no edge exceeds EDGE_TOL at lambda[", di, "] for ", fit$group, call. = FALSE)
  if (length(keep) < 5L)
    stop(sprintf("only %d edge(s) at lambda[%d] for %s -- too few to cluster. %s",
                 length(keep), di, fit$group,
                 if (POINT == "deployed") "Use TVCG_ZEB_POINT=bicmin (the default)." else ""),
         call. = FALSE)
  list(TR = do.call(rbind, keep), ij = do.call(rbind, ij),
       lambda = fit$detail$lambda[di], lambda_index = di, n_edge = length(keep),
       ## the refit at this lambda: did it reach the gate, and how far off is it?
       refit_converged = isTRUE(fit$detail$converged[di]),
       beta_grad = rf$max_active_grad %||% NA_real_,
       is_deployed = identical(di, fit$deployed$lambda_index))
}

#' k-means on z-scored shapes, k chosen by average silhouette width.
cluster_traj <- function(TR) {
  n <- nrow(TR)
  Z <- t(scale(t(TR)))                       # shape, not magnitude
  kmax <- min(6L, max(2L, floor(n / 5)))
  set.seed(SEED); d <- stats::dist(Z)
  sil <- vapply(2:kmax, function(k)
    mean(cluster::silhouette(stats::kmeans(Z, centers = k, nstart = 50)$cluster, d)[, 3]),
    numeric(1))
  k <- (2:kmax)[which.max(sil)]
  set.seed(SEED); km <- stats::kmeans(Z, centers = k, nstart = 200)
  ## order clusters by WHEN the mean trajectory peaks, so columns read
  ## left-to-right as early -> late and the two rows are comparable
  ord <- order(vapply(seq_len(k), function(kk)
    which.max(abs(colMeans(TR[km$cluster == kk, , drop = FALSE]))), integer(1)))
  list(k = k, cl = km$cluster, order = ord, silhouette = max(sil))
}

dat <- lapply(GROUPS_KEY, function(g) {
  tj <- collect_traj(fits[[g]])
  c(tj, cluster_traj(tj$TR), list(group = g))
})
names(dat) <- GROUPS_KEY

## one y-range for every panel: cluster heights are only comparable on a shared axis
YL <- range(unlist(lapply(dat, function(z) z$TR)))
YL <- YL + c(-1, 1) * 0.04 * diff(YL)
NCOL <- max(vapply(dat, function(z) z$k, integer(1)))

draw_fig05 <- function() {
  nr <- length(dat)
  layout(rbind(matrix(seq_len(nr * NCOL), nrow = nr, byrow = TRUE), nr * NCOL + 1L),
         heights = c(rep(1, nr), 0.16))
  par(mar = c(1.2, 1.1, 1.5, 0.6), oma = c(0.6, 3.0, 0.2, 2.0),
      mgp = c(3, 0.25, 0), tcl = -0.20, cex = 1, las = 1)

  xat <- c(7, 30, 60, 86)
  yat <- pretty(YL, 4)
  for (ri in seq_len(nr)) {
    z <- dat[[ri]]
    for (ci in seq_len(NCOL)) {
      if (ci > z$k) { plot.new(); next }          # ragged: fewer clusters here
      kk <- z$order[ci]; mem <- which(z$cl == kk)
      mu <- colMeans(z$TR[mem, , drop = FALSE])

      plot.new(); plot.window(xlim = range(days), ylim = YL)
      abline(h = 0, col = "grey80", lwd = LWD_REF)
      for (e in mem) lines(days, z$TR[e, ], col = "grey72", lwd = LWD_REF)
      lines(days, mu, col = "black", lwd = LWD_CURVE * 1.4)
      points(days, mu, pch = 16, cex = 0.35, col = "black")

      panel_axes(xat, yat,
                 xlab_show = (ri == nr), ylab_show = (ci == 1),
                 xfmt = as.character(xat), yfmt = format(yat, trim = TRUE))
      mtext(sprintf("cluster %d  (n = %d)", ci, length(mem)), side = 3, line = 0.2,
            cex = par("cex") * CEX_STRIP)
    }
    mtext(GROUP_LAB[[z$group]], side = 4, line = 0.5, cex = par("cex") * CEX_STRIP,
          las = 0, at = 0.5 * diff(YL) + YL[1])
  }
  mtext(expression(hat(Omega)[italic(ij)](italic(t))), side = 2, outer = TRUE,
        line = 1.7, cex = par("cex") * CEX_LAB, las = 0)

  par(mar = c(0, 0, 0, 0)); plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0.5, 0.80, "Days post exposure", adj = c(0.5, 1),
       cex = par("cex") * CEX_LAB, xpd = NA)
  legend(0.5, 0.22, xjust = 0.5, yjust = 0.5, ncol = 2, bty = "n",
         cex = CEX_LEGEND, seg.len = 2.4,
         legend = c("individual edges", "cluster mean"),
         col = c("grey72", "black"), lwd = c(LWD_REF, LWD_CURVE * 1.4), xpd = NA)
}

save_figure(sprintf("fig05_zeb_clusters_P%02d", P_TARGET),
            width = WIDTH_FULL, height = 3.6, draw = draw_fig05)

## --- the membership table behind the picture --------------------------------
rows <- do.call(rbind, lapply(dat, function(z) {
  do.call(rbind, lapply(seq_len(z$k), function(ci) {
    kk <- z$order[ci]; mem <- which(z$cl == kk)
    mu <- colMeans(z$TR[mem, , drop = FALSE])
    data.frame(P = P_TARGET, group = z$group, point = POINT,
               lambda_index = z$lambda_index, lambda = z$lambda,
               refit_converged = z$refit_converged, beta_grad = z$beta_grad,
               cluster = ci, n_edge = length(mem),
               peak_day = days[which.max(abs(mu))],
               mean_at_peak = mu[which.max(abs(mu))],
               net_change = mu[m] - mu[1],
               edges = paste(sprintf("%s-%s", nodes[z$ij[mem, 1]], nodes[z$ij[mem, 2]]),
                             collapse = "; "),
               stringsAsFactors = FALSE)
  }))
}))
f_out <- file.path(TAB_DIR, sprintf("tab05_zeb_clusters_P%02d.csv", P_TARGET))
write.csv(rows, f_out, row.names = FALSE)
message("wrote ", f_out)
message(sprintf("operating point: %s", POINT))
for (g in names(dat)) {
  z <- dat[[g]]
  message(sprintf("  %-13s lambda[%2d]=%.4f  edges=%3d  k=%d (silhouette %.3f)  refit %s (|g|=%.2e)%s",
                  g, z$lambda_index, z$lambda, z$n_edge, z$k, z$silhouette,
                  if (z$refit_converged) "CONVERGED" else "NOT converged", z$beta_grad,
                  if (z$is_deployed) "  [= this fit's deployed point]" else ""))
}
if (any(!vapply(dat, function(z) z$refit_converged, logical(1))))
  message("  NOTE: at least one panel uses a refit that did not reach its convergence gate.\n",
          "        Shapes are robust to this; the vertical scale is not (header caveat 1).")
print(rows[, c("group", "cluster", "n_edge", "peak_day", "mean_at_peak", "net_change")],
      row.names = FALSE, digits = 3)
