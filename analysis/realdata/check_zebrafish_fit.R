# ---------------------------------------------------------------------------
# analysis/realdata/check_zebrafish_fit.R
#
#   Rscript analysis/realdata/check_zebrafish_fit.R            # every fit found
#   TVCG_ZEB_P=15 Rscript analysis/realdata/check_zebrafish_fit.R
#
# Decide whether a zebrafish fit is fit to build a figure on, BEFORE building
# one. On real data there is no truth, so none of the simulation's checks apply.
# What can still be checked is that the estimator behaved: that the penalty path
# spans empty to dense, that the deployed point was actually SELECTED rather than
# arrived at by default or by fallback, and that the fit at that point converged.
#
# The one thing this script must not do is measure a property over a set the
# selector never sees. The deployed point is the BIC-minimiser over the lambda
# whose refit CONVERGED; a BIC surface that looks well separated over the whole
# path says nothing if only one lambda is eligible. Every quantity below is
# therefore computed over the CANDIDATE set, and the candidate set is reported
# first.
#
# Nothing here is a quality measure. A fit can pass every check and still be
# describing noise -- see the caveat in scripts/fig05_zebrafish_clusters.R.
# ---------------------------------------------------------------------------

source(here::here("analysis", "config.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

files <- list.files(ZEB_FITS_DIR, pattern = "^zeb_fit_P[0-9]{2}_.*\\.rds$", full.names = TRUE)
if (nzchar(Sys.getenv("TVCG_ZEB_P")))
  files <- files[grepl(sprintf("_P%02d_", as.integer(Sys.getenv("TVCG_ZEB_P"))), files)]
if (!length(files)) stop("no fits in ", ZEB_FITS_DIR, call. = FALSE)

## edges pooled over slices counts an edge once per slice it is active on; the
## number of distinct (i,j) pairs is what a trajectory figure can actually draw.
n_pairs <- function(Omega_list, P) {
  if (is.null(Omega_list)) return(NA_integer_)
  ut <- upper.tri(matrix(0, P, P))
  M  <- vapply(Omega_list, function(A) A[ut], numeric(sum(ut)))
  sum(rowSums(abs(M) > 0) > 0)
}

report <- function(f) {
  x <- readRDS(f); d <- x$detail; dp <- x$deployed; R <- x$estimate$refit
  n_lam <- nrow(d); di <- dp$lambda_index
  nz <- d$n_edges > 0
  maxpair <- (x$P * (x$P - 1L)) / 2L * x$m
  flag <- function(cond, msg) if (isTRUE(cond)) cat("    !! ", msg, "\n", sep = "")

  cat(sprintf("\n===== %s   P=%d  %s =====\n", basename(f), x$P, x$group))
  cat(sprintf("  n per day %s   N=%d   m=%d   %.0f min\n",
              paste(x$n_per_slice, collapse = ","), x$N, x$m, x$secs / 60))

  ## ---- 1. the path -------------------------------------------------------
  cat(sprintf("\n  [path]      %d lambda in [%.4g, %.4g]\n", n_lam, min(d$lambda), max(d$lambda)))
  cat(sprintf("              %d empty-graph lambda (%.0f%% of the path)\n", sum(!nz), 100 * mean(!nz)))
  cat(sprintf("              densest fit %d edge-slices = %.0f%% of %d possible\n",
              max(d$n_edges), 100 * max(d$n_edges) / maxpair, maxpair))
  flag(sum(!nz) == 0, "path never reaches the empty graph")
  flag(max(d$n_edges) / maxpair < 0.25, "dense end is thin; the path may stop too early")
  flag(mean(!nz) > 0.4, sprintf("%.0f%% of the path is all-empty fits (wasted grid, disclosable)",
                                100 * mean(!nz)))

  ## ---- 2. how the refit exited, and which half of the gate failed --------
  ex <- vapply(R, function(r) r$exit_reason %||% "error", character(1))
  cat(sprintf("\n  [refit]     converged %d | max_outer %d | error %d   (of %d)\n",
              sum(ex == "converged"), sum(ex == "max_outer"), sum(ex == "error"), n_lam))
  bg <- vapply(R, function(r) r$max_active_grad %||% NA_real_, numeric(1))
  zg <- vapply(R, function(r) r$max_z_grad      %||% NA_real_, numeric(1))
  fail <- which(ex == "max_outer")
  if (length(fail)) {
    fb <- sum(bg[fail] > 1e-4, na.rm = TRUE); fz <- sum(zg[fail] > 1e-4, na.rm = TRUE)
    cat(sprintf("              of the %d that hit the cap: %d fail on beta, %d on Z\n",
                length(fail), fb, fz))
    cat(sprintf("              beta-gradient there spans [%.1e, %.1e]  (tolerance 1e-4)\n",
                min(bg[fail], na.rm = TRUE), max(bg[fail], na.rm = TRUE)))
    ## an empty support leaves only the always-free diagonal, whose refit is
    ## strictly convex and richly identified -- it cannot fail for want of a
    ## solution, so a failure there is a solver or setting problem, not the
    ## non-existence that a dense support at n < P genuinely can hit.
    ef <- intersect(fail, which(!nz))
    flag(length(ef) > 0,
         sprintf("%d EMPTY-graph fits did not converge (beta-grad up to %.1e). An empty support fits only the diagonal, whose optimum exists -- this is not non-existence.",
                 length(ef), max(bg[ef], na.rm = TRUE)))
  }

  ## ---- 3. THE CANDIDATE SET -- what the selector could actually choose ---
  cand <- which(d$converged & is.finite(d$BIC))
  cand_ne <- intersect(cand, which(nz))
  cat(sprintf("\n  [candidates] %d lambda are eligible (converged AND finite BIC)\n", length(cand)))
  cat(sprintf("               of these %d are NON-EMPTY\n", length(cand_ne)))
  flag(length(cand) == 0, "NO eligible lambda at all: the deployed point is a FALLBACK, not a selection")
  flag(length(cand) > 0 && length(cand_ne) == 0, "every eligible lambda is the empty graph")
  flag(length(cand_ne) == 1,
       "exactly ONE non-empty eligible lambda: the deployed point is the only alternative to the empty graph, not the winner of a comparison")

  if (length(cand) > 1) {
    b <- d$BIC[cand]; rng <- diff(range(b))
    near <- if (rng > 0) sum(b < min(b) + 0.01 * rng) else length(b)
    cat(sprintf("               BIC over the candidates spans %.1f; %d within 1%% of the minimum\n",
                rng, near))
    flag(near > 0.2 * length(cand), "BIC is nearly flat over the candidates: the argmin is weakly identified")
  }

  ## ---- 4. the deployed point ---------------------------------------------
  fb_used <- !isTRUE(dp$refit_converged)
  np <- n_pairs(R[[di]]$Omega_post, x$P)
  cat(sprintf("\n  [deployed]  lambda[%d/%d] = %.4g   %s\n", di, n_lam, dp$lambda,
              if (fb_used) "<-- FALLBACK (its refit did NOT converge)" else "(refit converged)"))
  cat(sprintf("              %d edge-slices, %s distinct (i,j) pairs, df=%d\n",
              dp$n_edges, format(np), dp$df))
  cat(sprintf("              beta-grad %.2e, z-grad %.2e   (tolerance 1e-4)\n", bg[di], zg[di]))
  bicmin <- which.min(d$BIC)
  cat(sprintf("              path-wide BIC minimum is at lambda[%d] (%d edge-slices, converged=%s)\n",
              bicmin, d$n_edges[bicmin], d$converged[bicmin]))
  flag(fb_used, "deployed via the all-lambda fallback: this point was not selected by the stated rule")
  flag(di == 1L || di == n_lam, "deployed AT a path boundary: the selector wanted to leave the grid")
  flag(isTRUE(dp$n_edges == 0L), "deploys the EMPTY graph")
  flag(!is.na(np) && np < 5, sprintf("only %s distinct edge pair(s): too few for a trajectory-clustering figure", format(np)))
  flag(bicmin != di, "the path-wide BIC minimum is NOT the deployed point (convergence, not BIC, is binding)")
  invisible(NULL)
}
invisible(lapply(sort(files), report))

cat("\n---\nLines marked !! are things to resolve before using a fit in a figure.\n")
cat("Two fits deployed by different mechanisms (one selected, one fallback) are\n")
cat("not comparable to each other, however similar the plots look.\n")
