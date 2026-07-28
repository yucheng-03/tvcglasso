# ---------------------------------------------------------------------------
# analysis/make.R  --  regenerate EVERY figure and table in the paper.
#
#   Rscript analysis/make.R            # from the repository root
#
# This is the single documented command referred to by analysis/README.md and
# by ACC Part 3 ("Reproducibility workflow"). It rebuilds all numbered figures
# and tables from the stored per-(method, cell, seed) result files; it does not
# re-run the simulation campaign (see README, Tier 2).
#
# Environment switches:
#   TVCG_REFRESH=TRUE   ignore analysis/cache/ and re-extract from the rds files
#   TVCG_STRICT=TRUE    treat an incomplete cell as an error rather than a warning
# ---------------------------------------------------------------------------

t0 <- Sys.time()
scripts <- c(
  "scripts/fig01_roc_grid.R",       # Figure 1  + tables/tab02_roc_coverage.csv
  "scripts/fig02_deployed_bars.R",  # Figure 2  + tables/tab01_deployed.csv
  "scripts/fig03_generators.R",     # Figure 3  (true edge-trajectory families)
  "scripts/tab03_provenance.R"      # tables/tab03_provenance.csv
)

root <- here::here()
run <- function(s, env = character()) {
  message("\n==== ", s, if (length(env)) paste0("  [", paste(env, collapse = " "), "]") else "", " ====")
  ## names = TRUE is REQUIRED: Sys.getenv() returns an UNNAMED scalar when given
  ## a single variable name (it only names the result when length(x) > 1), so
  ## old[[k]] below is a subscript-out-of-bounds error without it.
  old <- if (length(env)) Sys.getenv(names(env), unset = NA, names = TRUE) else character()
  if (length(env)) do.call(Sys.setenv, as.list(env))
  on.exit({
    for (k in names(env)) {
      if (is.na(old[[k]])) Sys.unsetenv(k)
      else do.call(Sys.setenv, stats::setNames(list(old[[k]]), k))
    }
  }, add = TRUE)
  ## each script is self-contained: it sources config.R and analysis/R/*.R
  ## itself, so it can also be run alone.
  source(file.path(root, "analysis", s), echo = FALSE, local = new.env())
}

for (s in scripts) {
  ## fig02 is drawn once per m: JGL exists only at m <= 7 (its ADMM does not
  ## scale), so putting both m in one panel would render an inapplicable method
  ## as four missing bars. Rendering only the default m silently ships a stale
  ## m = 15 figure, which is what happened before this loop existed.
  if (grepl("fig02_", s)) {
    for (mm in c("7", "15")) run(s, c(TVCG_M = mm))
  } else {
    run(s)
  }
}

message("\nall figures and tables rebuilt in ",
        format(round(difftime(Sys.time(), t0, units = "mins"), 2)))
message("R: ", R.version.string)
writeLines(capture.output(utils::sessionInfo()),
           file.path(here::here("analysis"), "sessionInfo.txt"))
message("wrote analysis/sessionInfo.txt")
