# ---------------------------------------------------------------------------
# analysis/config.R
#
# The ONLY file in analysis/ that knows machine-specific locations.
# Everything else uses relative paths through here::here() or the values below.
#
# Reproducibility note (JASA ACC Part 3): a reviewer who has the shipped
# per-(method, cell, seed) result files needs to edit nothing but RESULT_ROOTS.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(here))

## --- where the fitted results live -----------------------------------------
## One entry per method. Each directory holds cellCC_seedSSS.rds files written
## by simulation/run_comparison.R. They may sit in different trees (the four
## legs were run as separate campaigns); nothing downstream cares.
##
## Set TVCG_RESULTS_DIR / CGLASSO_RESULTS_DIR / ... in the environment to
## override without editing this file.

.env_or <- function(var, default) {
  v <- Sys.getenv(var, unset = "")
  if (nzchar(v)) v else default
}

## The shipped layout is results/<method>/. During the campaign the four legs
## were run from different trees, so each root is independently overridable.
RESULT_ROOTS <- c(
  tvcglasso = .env_or("TVCG_RESULTS_DIR",    here::here("results", "tvcglasso")),
  CGLasso   = .env_or("CGLASSO_RESULTS_DIR", here::here("results", "CGLasso")),
  tvmgm     = .env_or("TVMGM_RESULTS_DIR",   here::here("results", "tvmgm")),
  JGL       = .env_or("JGL_RESULTS_DIR",     here::here("results", "JGL"))
)

## Optional second segment of the CGLasso rho path (the dense-end extension).
## If present, fig01 stitches it onto the main path -- both segments are REAL
## fitted points, nothing is interpolated across the seam. Set to NA to skip.
CGLASSO_DENSEEND_DIR <- .env_or("CGLASSO_DENSEEND_DIR",
                                here::here("results", "cglasso_denseend", "CGLasso"))

## --- the design grid -------------------------------------------------------
## 16 paired cells, seeds 1:100. Defined in config/paired_grid.R; repeated here
## only as the EXPECTED extent, so the loaders can assert completeness instead
## of silently averaging over whatever happens to be on disk.
## TVCG_CELLS / TVCG_SEEDS accept R range syntax ("1:16", "c(1,3,5)") and exist
## so a reviewer can run a fast subset before committing to the full grid.
.parse_ints <- function(var, default) {
  v <- Sys.getenv(var, unset = "")
  if (!nzchar(v)) return(default)
  as.integer(eval(parse(text = v)))
}
CELLS <- .parse_ints("TVCG_CELLS", 1:16)
SEEDS <- .parse_ints("TVCG_SEEDS", 1:100)

## JGL is feasible only at m <= 7 (its ADMM does not scale); it exists for the
## m = 7 cells only. This is a property of the method, not missing data.
JGL_CELLS <- 1:8

## --- real-data (Zebrafish) fitted objects ----------------------------------
## The per-group TVCGLasso fits used by the real-data figures. These are
## DERIVED objects (fitted Omega paths), not the source counts; the source
## tables belong to Gaulke et al. (2019) and are held back pending
## redistribution clearance (see data/README.md and .gitignore). Point this at
## wherever the fits live, or regenerate them with the real-data driver.
ZEB_FITS_DIR <- .env_or("ZEB_FITS_DIR", here::here("results", "zebrafish"))

## --- output ----------------------------------------------------------------
FIG_DIR   <- here::here("analysis", "figures")
TAB_DIR   <- here::here("analysis", "tables")
CACHE_DIR <- here::here("analysis", "cache")

for (d in c(FIG_DIR, TAB_DIR, CACHE_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## --- strictness ------------------------------------------------------------
## TRUE  : a cell with fewer than SEEDS seeds is a hard error (publication runs).
## FALSE : incomplete cells are allowed and the realised n is reported per curve
##         (use while a campaign leg is still finishing).
STRICT_COMPLETENESS <- as.logical(.env_or("TVCG_STRICT", "FALSE"))
