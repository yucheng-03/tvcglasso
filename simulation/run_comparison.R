# ============================================================================
# simulation/run_comparison.R — MASTER driver.
#
#   Rscript simulation/run_comparison.R                        # tiny TV+refit smoke
#   Rscript simulation/run_comparison.R config/publication.R   # config must define CFG
#
# Run and plotting are DECOUPLED: DATA ONLY here. Full per-fit intermediate is saved so
# any figure / selector is recomputable later WITHOUT re-fitting.
#
# ★ P0-3 (2026-07-22): one file PER (METHOD, cell, seed): <out>/<tag>/<method>/cellCC_seedSSS.rds
#   -> a single (cell,seed,method) OR many can be run; adding a method later BACKFILLS (its
#   files are missing -> it runs, the others skip). Data is generated ONCE per (cell,seed)
#   (deterministic in seed -> identical paired data across methods).
# ★ P0-2 (validated-DONE): a method's result file is written ONLY if the fit is VALID (no
#   error). A failed method writes a .FAILED marker (not a result) and KEEPS its checkpoints,
#   so a re-submit re-runs it (never a silent NA baptized as done, never permanently skipped).
# ★ §4-1: a config PATH that does not exist is a hard ERROR (never silently run the smoke).
# ★ §4-3: each result file embeds a reproducibility receipt (config + git hash + data
#   fingerprint + R version + generator/seed).
#
# Discipline: strict generator dispatch (no silent fallback); per-lambda TV + refit
# checkpoints (preempt-safe); each method wrapped in tryCatch; mclapply(mc.preschedule=FALSE);
# SLURM array mode. Baseline packages load ONLY if requested.
# ============================================================================
suppressPackageStartupMessages({
  library(here); library(splines); library(glasso); library(Matrix); library(MASS)
  library(Rcpp); library(parallel)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- CFG: from a config file arg, else a tiny TV+refit smoke ----
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) {
  if (!file.exists(.args[1]))                                         # §4-1: no silent smoke fallback
    stop(sprintf("config file not found: '%s' — refusing to silently run the smoke.", .args[1]))
  source(.args[1]); cat(sprintf("[cfg] loaded %s\n", .args[1]))       # must define CFG
} else {
  CFG <- list(
    tag = "smoke",
    cells = list(list(generator = "mirrorexp_contB", P = 8, n = 15, m = 7,
                      depth_mode = "low", edge_strength = 0.5, rate = 6)),
    seeds = 1:1, methods = c("tvcglasso"),
    lambda_grid_tv = NULL, tv_n_lambda = 8L,
    tv_q = 2L, tv_N_n = 1L, tv_init_mode = "diag", tv_sel_type = "hard",
    tv_weight_mode = "glasso", tv_max_iter = 40L, tv_free_diag = TRUE,
    tv_refit = TRUE, refit_max_outer = 40L, refit_z_align_max = 60L,
    out_dir = here::here("results"), ncores = 1L)
  cat("[cfg] using default tiny TV+refit smoke CFG\n")
}

## ---- source engine + methods (baselines only if requested) ----
source(here::here("R", "tvcglasso.R"))        # TV engine + tv_warm_path + G_beta_Rcpp
source(here::here("R", "refit.R"))            # refit_fixed_support + IC helpers
source(here::here("R", "roc_utils.R"))        # method contract + ROC helpers + %||%
source(here::here("simulation", "generators.R"))  # generate_data (strict dispatch)
source(here::here("R", "methods", "method_tv.R"))
if ("CGLasso" %in% CFG$methods) {
  if (!requireNamespace("huge", quietly = TRUE) || !requireNamespace("propagate", quietly = TRUE))
    stop("CGLasso requires the 'huge' and 'propagate' packages")
  suppressPackageStartupMessages({ library(huge); library(propagate) })
  source(here::here("R", "baselines", "CompoGlasso.R"))     # cg_-namespaced (P0-1): does NOT overwrite the TV engine
  source(here::here("R", "baselines", "cglasso_static.R"))
  source(here::here("R", "baselines", "stars_cglasso.R"))
  source(here::here("R", "baselines", "refit_cglasso_core.R"))   # CG's OWN refit core (cg_-prefixed) — INDEPENDENT of TV's R/refit.R
  source(here::here("R", "baselines", "refit_cglasso.R"))
  source(here::here("R", "methods", "method_cglasso.R"))
}
if ("tvmgm" %in% CFG$methods) {
  if (!requireNamespace("mgm", quietly = TRUE)) stop("tvmgm requires the 'mgm' package")
  suppressPackageStartupMessages(library(mgm))
  source(here::here("R", "baselines", "refit_mgm.R"))
  source(here::here("R", "methods", "method_tvmgm.R"))
}
if ("JGL" %in% CFG$methods) {
  if (!requireNamespace("JGL", quietly = TRUE)) stop("JGL requires the 'JGL' package")
  suppressPackageStartupMessages({ library(JGL); library(glasso) })
  source(here::here("R", "baselines", "refit_jgl.R"))
  source(here::here("R", "methods", "method_jgl.R"))
}

run_one_method <- function(name, dat, cfg, eval_slices) {
  switch(name,
    "tvcglasso" = run_method_tv(dat, cfg, eval_slices),
    "CGLasso" = run_method_cglasso(dat, cfg, eval_slices),
    "tvmgm"   = run_method_tvmgm(dat, cfg, eval_slices),
    "JGL"     = run_method_jgl(dat, cfg, eval_slices),
    stop("unknown method: ", name))
}

## ---- §4-3 provenance: git hash (once) + a cheap deterministic data fingerprint ----
.git_hash <- local({
  env <- Sys.getenv("TVCG_GIT_HASH")                       # injected by the submit script (cluster: not a git repo)
  if (nzchar(env)) return(env)
  tryCatch({ h <- system(sprintf("git -C %s rev-parse --short HEAD 2>/dev/null",
             shQuote(here::here())), intern = TRUE); if (length(h)) h else "(no git)" },
           error = function(e) "(no git)")
})
data_fingerprint <- function(dat) sum(vapply(dat$X, function(M) sum(as.numeric(M)), numeric(1)))

## ---- output layout: <out_dir>/<tag>/{<method>/, _ckpt/, run.log} ----
run_dir  <- file.path(CFG$out_dir, CFG$tag)
ckpt_dir <- file.path(run_dir, "_ckpt"); dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)
for (nm in CFG$methods) dir.create(file.path(run_dir, nm), recursive = TRUE, showWarnings = FALSE)
mfile <- function(nm, ci, sd) file.path(run_dir, nm, sprintf("cell%02d_seed%03d.rds", ci, sd))

## ---- task list = (cell_idx, seed) [data shared across methods] + SLURM array slicing ----
tasks <- do.call(rbind, lapply(seq_along(CFG$cells), function(ci)
  data.frame(cell = ci, seed = CFG$seeds)))
.arr <- Sys.getenv("SLURM_ARRAY_TASK_ID"); array_mode <- nzchar(.arr)
if (array_mode) {
  k <- as.integer(.arr); gran <- CFG$array_by %||% "seed"
  if (gran == "task") { tasks <- tasks[k, , drop = FALSE]
    cat(sprintf("[array] element %s -> cell%d seed%d\n", .arr, tasks$cell[1], tasks$seed[1]))
  } else { this_seed <- CFG$seeds[k]; tasks <- tasks[tasks$seed == this_seed, , drop = FALSE]
    cat(sprintf("[array] element %s -> seed %d (%d tasks)\n", .arr, this_seed, nrow(tasks))) }
}
cat(sprintf("[run] %d cells x %d seeds ; methods: %s%s\n",
            length(CFG$cells), length(CFG$seeds), paste(CFG$methods, collapse = ", "),
            if (array_mode) "  [ARRAY ELEMENT]" else ""))

## ---- one task: for each REQUESTED method missing its per-method file, run + VALIDATE + save ----
run_task <- function(ti) {
  ci <- tasks$cell[ti]; sd <- tasks$seed[ti]; cell <- CFG$cells[[ci]]
  todo <- CFG$methods[!vapply(CFG$methods, function(nm) file.exists(mfile(nm, ci, sd)), logical(1))]
  if (length(todo) == 0) { cat(sprintf("[skip] cell%d seed%d — all methods done\n", ci, sd)); return(invisible(NULL)) }
  dat <- generate_data(cell, sd)                                      # STRICT dispatch; deterministic in seed
  m <- length(dat$X); dat$seed <- sd; dfp <- data_fingerprint(dat)
  dat$tv_q <- cell$tv_q; dat$tv_N_n <- cell$tv_N_n                     # per-cell B-spline basis (scales with m)
  dat$tv_ckpt_file   <- file.path(ckpt_dir, sprintf("tvckpt_c%02d_s%03d.rds", ci, sd))
  dat$refit_ckpt_dir <- file.path(ckpt_dir, sprintf("refit_c%02d_s%03d", ci, sd))
  dir.create(dat$refit_ckpt_dir, showWarnings = FALSE, recursive = TRUE)
  eval_slices <- CFG$eval_slices %||% (1:m)
  # provenance = the REPRODUCIBILITY RECEIPT only (config + git + fingerprint + R + date). cell / seed /
  # generator are NOT duplicated here — they live once at the top level (out$cell, out$seed, cell$generator).
  prov <- list(config = CFG[setdiff(names(CFG), "cells")], git = .git_hash, data_fingerprint = dfp,
               R = R.version.string, date = as.character(Sys.Date()))

  for (nm in todo) {
    mt0 <- Sys.time()
    r <- tryCatch(run_one_method(nm, dat, CFG, eval_slices),
                  error = function(e) structure(list(method = nm, error = conditionMessage(e)), class = "fit_error"))
    r$secs <- as.numeric(difftime(Sys.time(), mt0, units = "secs"))
    valid <- is.null(r$error) && !is.null(r$roc)                      # P0-2: valid = no crash + produced an ROC
    if (valid) {
      out <- list(method = nm, cell_idx = ci, cell = cell, seed = sd, m = m, eval_slices = eval_slices,
                  result = r, true_Omega = dat$true_Omega_list, x_sequence = dat$x_sequence, provenance = prov)
      tmp <- paste0(mfile(nm, ci, sd), ".tmp"); saveRDS(out, tmp); file.rename(tmp, mfile(nm, ci, sd))  # atomic
      cat(sprintf("[ok]   %s cell%d(P=%d,n=%d) seed%d  %.0fs  AUC=%.3f\n",
                  nm, ci, cell$P, cell$n, sd, r$secs, r$auc %||% NA_real_))
    } else {                                                          # P0-2: failure -> marker, keep checkpoints
      writeLines(c(sprintf("FAILED %s cell%d seed%d", nm, ci, sd), r$error %||% "no ROC produced"),
                 paste0(mfile(nm, ci, sd), ".FAILED"))
      cat(sprintf("[FAIL] %s cell%d seed%d: %s (checkpoints kept for resume)\n", nm, ci, sd, r$error %||% "no ROC"))
    }
  }
  # drop scratch checkpoints ONLY if every requested method for this (cell,seed) now has a result file
  if (all(vapply(CFG$methods, function(nm) file.exists(mfile(nm, ci, sd)), logical(1)))) {
    if (file.exists(dat$tv_ckpt_file)) unlink(dat$tv_ckpt_file)
    if (dir.exists(dat$refit_ckpt_dir)) unlink(dat$refit_ckpt_dir, recursive = TRUE)
  }
  invisible(NULL)
}

## ---- run (parallel) ----
ncores <- CFG$ncores %||% { s <- Sys.getenv("SLURM_CPUS_PER_TASK"); if (nzchar(s)) as.integer(s) else max(1, detectCores() - 1) }
t_start <- Sys.time(); cat(sprintf("[run] %d cores, start %s\n", ncores, format(t_start)))
invisible(mclapply(seq_len(nrow(tasks)), run_task, mc.cores = ncores, mc.preschedule = FALSE))

## ---- run.log (SKIPPED in array mode; each array element already wrote its per-method files) ----
if (!array_mode) {
  tryCatch({
    con <- file(file.path(run_dir, "run.log"), "w")
    writeLines(c(
      sprintf("run_comparison — tag=%s", CFG$tag),
      sprintf("date: %s   R: %s   git: %s", Sys.Date(), R.version.string, .git_hash),
      sprintf("TV engine: R/tvcglasso.R (option2 pseudocount [init only, P0-5], free_diag=%s, weight=%s, sel=%s, max_iter=%s, refit=%s)",
              CFG$tv_free_diag %||% TRUE, CFG$tv_weight_mode, CFG$tv_sel_type, CFG$tv_max_iter %||% 150L, isTRUE(CFG$tv_refit %||% TRUE)),
      sprintf("methods: %s", paste(CFG$methods, collapse = ", ")),
      sprintf("cells: %s", paste(sapply(CFG$cells, function(c)
        sprintf("%s/P%d/n%d/m%d/%s", c$generator %||% "?", c$P, c$n, c$m, c$depth_mode %||% "?")), collapse = "  ")),
      sprintf("seeds: %s   wall-clock: %.0f s", paste(range(CFG$seeds), collapse = "-"),
              as.numeric(difftime(Sys.time(), t_start, units = "secs")))), con)
    close(con); cat("[log] run.log saved\n")
  }, error = function(e) cat("run.log FAILED (per-method .rds are safe):", conditionMessage(e), "\n"))
}
cat(sprintf("\n[done] results in %s\n", run_dir))
