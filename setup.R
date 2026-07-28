# setup.R — install the R packages the tvcglasso repo needs.
#   Rscript setup.R            # core only (TV engine + refit + generators + run driver)
#   Rscript setup.R baselines  # core + all baseline R-package deps (for CGLasso/tvmgm/JGL runs)
#
# ★ §4-2 NOTE: CGLasso is VENDORED SOURCE (R/baselines/CompoGlasso.R, committed) — NOT an R
# package, so there is nothing to install for CGLasso itself. Its R-package DEPENDENCIES are
# `huge` + `propagate` (used by bigcor) + `glasso`. Those (and mgm/JGL for the other baselines)
# are installed by `Rscript setup.R baselines`.
core <- c("Rcpp", "glasso", "Matrix", "MASS", "here")   # `splines` is base R
baseline_deps <- c("huge", "propagate",                 # CGLasso (vendored source) needs these
                   "mgm",                                # tvmgm baseline
                   "JGL")                                # JGL baseline (small-m)

install_missing <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
  else cat("  (all present)\n")
}

cat("Installing core packages:\n"); install_missing(core)
if (identical(commandArgs(trailingOnly = TRUE)[1], "baselines")) {
  cat("Installing baseline R-package deps (huge/propagate for the vendored CGLasso; mgm; JGL):\n")
  install_missing(baseline_deps)
} else {
  cat("Core ready. Run `Rscript setup.R baselines` to also install CGLasso/tvmgm/JGL deps.\n")
}
