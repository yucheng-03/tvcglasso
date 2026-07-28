# ============================================================================
# config/pub_tvcglasso.R — TV (tvcglasso) publication run = OUR METHOD ONLY, 100 seeds (per-method run file; shares paired_grid + tag).
#
#   Rscript simulation/run_comparison.R config/pub_tvcglasso.R
#
# Purpose: produce the DRAFT-paper TV results at PUBLICATION resolution (70 lambda,
# refit at every lambda, free diagonal, raw-count likelihood, n_k-weighted). The
# CGLasso/JGL/mgm baselines are run from their OWN windows on the SAME tag ->
# they BACKFILL results/publication_tv/<method>/ on the IDENTICAL paired data
# (config/paired_grid.R). Topping seeds 30 -> 100 for SUBMISSION re-runs nothing.
# ============================================================================

source(here::here("config", "paired_grid.R"))     # PAIRED_SEEDS, PAIRED_CELLS (shared spec)

CFG <- list(
  tag     = "publication_tv",                     # SAME tag as the full record -> baselines backfill here
  methods = c("tvcglasso"),                       # this window = OUR method only; CGLasso/JGL/mgm backfill later
  seeds   = PAIRED_SEEDS,                          # ALL 100 (submission, 2026-07-24; was [1:30] draft)
  cells   = PAIRED_CELLS,

  # ---- lambda path: 70 points (Yuan/Tian CGLasso convention = same resolution as the baseline),
  #      data-adaptive per (cell,seed): lambda_max = 2*max|off-diag cov(ALR Z_0)| (guarantees the empty
  #      graph / bottom-left corner) down to lambda_max/1000 (dense end), log-spaced. TV's hard-threshold
  #      caps max density so the top-right corner is APPROACHED (maxFPR ~0.9-1.0) with REAL lambda. ----
  lambda_grid_tv = NULL, tv_n_lambda = 50L,   # 50 (2026-07-24, reduced from 70): ALL 4 methods use 50
  tv_q = 2L, tv_N_n = 1L, tv_init_mode = "diag", tv_sel_type = "hard",
  tv_weight_mode = "glasso", tv_max_iter = 150L, tv_sel_threshold = NULL, tv_free_diag = TRUE,

  # ---- refit (pass 2) at EVERY lambda: de-biased estimator + joint-LNM-at-Zhat {AIC,BIC,eBIC}
  #      pre & post, plus base convergence status -> any selector recomputable at analysis, no re-fit. ----
  tv_refit = TRUE,
  tv_refit_max_outer = 500L, refit_inner_max = 10L, refit_z_align_max = 160L,   # TV-OWN cap = 500 (2026-07-24; 150 too low at low-depth n<P, free-diag refit slower)
  refit_conv_tol_Z = 5e-5, refit_conv_tol_z_grad = 1e-4,

  # ---- run/output ----
  array_by = "seed",                              # 1 array element = 1 seed x all 16 cells (mclapply)
  out_dir  = here::here("results"),
  ncores   = NULL                                 # NULL => SLURM_CPUS_PER_TASK
)
