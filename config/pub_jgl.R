# ============================================================================
# config/pub_jgl.R — JGL (Danaher et al. 2014) publication run = THIS BASELINE ONLY.
#
#   Rscript simulation/run_comparison.R config/pub_jgl.R
#
# Per-method run file (the file we actually launch). It sources the SHARED paired-data
# spec but writes to its OWN tag/dir (see below), so it cannot disturb another method's files.
# Running the per-method files separately is bit-identical to running
# config/publication.R once (identical seeds + generator => identical data).
#
# METHOD (see R/methods/method_jgl.R + R/baselines/refit_jgl.R): base = the AUTHORS' own
# JGL::JGL(penalty="fused") on the centered plug-in ALR, over a (lambda1 x lambda2) grid --
# a faithful reproduction, all package defaults respected. Refit = relaxed lasso on the fixed
# pass-1 support: dropping BOTH penalties decouples the slices, so each slice becomes the
# unpenalized constrained Gaussian MLE glasso(S_k, rho=0, zero=non-edges) (Dempster 1972 /
# ESL 2e Alg 17.1), at glasso's OWN default tolerance -- the SAME solver settings as the
# CGLasso refit. That MLE does not exist at the dense end of the lambda1 path (n<P; Grone et al.
# 1984 / Buhl 1993), where the unridged solve does not even TERMINATE, so pass 2 is solved on a
# RIDGED slice covariance (jgl_refit_ridge_eps below) -- ours, disclosed. Any point still not
# finite+PD is flagged refit-unavailable and excluded from the deployment candidates.
#
# ONLY the m=7 cells: JGL's fused ADMM is a small-m baseline (infeasible at m=15/30), and
# run_method_jgl hard-stops m>7. => 8 cells x 100 seeds.
#
# OUTPUT follows the SHARED shell contract (same layout as tvcglasso / CGLasso):
#   $roc $auc $scores + $detail (FLAT per-(lambda1,lambda2) table: FPR/TPR/precision/F1/df/n_edges
#   + POST-refit {neg2loglik,AIC,BIC,eBIC} + PRE-refit {*_pre} + refit_exists/base_pd)
#   + $deployed (the as-deployed point at the native selector, a direct read)
#   + $estimate (payload: base$Theta / refit$Theta at EVERY grid point, row-aligned with $detail).
# So every selector's deployed point and the base-vs-refit comparison are recomputable at
# analysis with NO re-fit. JGL-specific: the path point is 2-D and the payload has Theta only
# (no beta -- no basis; no Z -- no latent layer).
# ============================================================================
source(here::here("config", "paired_grid.R"))   # PAIRED_SEEDS (100), PAIRED_CELLS (16) — shared spec

CFG <- list(
  # ---- OUTPUT ISOLATION (hard rule: each method owns its files and cannot disturb another's).
  #      JGL gets its OWN tag, hence its OWN results dir AND its OWN _ckpt tree:
  #          results/publication_jgl/JGL/cellCC_seedSSS.rds
  #          results/publication_jgl/_ckpt/...
  #      This is NOT cosmetic. simulation/run_comparison.R cleans up scratch checkpoints with
  #          if (all(methods of THIS config have their rds)) unlink(refit_c<CC>_s<SSS>, recursive = TRUE)
  #      and that directory is keyed by (cell,seed) ONLY -- not by method. So any two configs that
  #      SHARE a tag also share that directory, and whichever finishes a (cell,seed) first DELETES
  #      the other's in-progress checkpoints. Separate tags make the interference impossible by
  #      construction rather than by timing. Analysis reads the per-method dirs side by side; the
  #      data are paired regardless (same paired_grid => bit-identical (cell,seed) data).
  tag     = "publication_jgl",
  methods = c("JGL"),             # this window = the JGL baseline only

  # ---- paired data (NEVER redefine locally; identical seeds+generator => bit-identical
  #      data as every other method's window) ----
  seeds = PAIRED_SEEDS,                                        # ALL 100 (submission)
  cells = Filter(function(cc) cc$m == 7L, PAIRED_CELLS),       # the 8 m=7 cells (JGL is m<=7 only)

  # ---- JGL base fit (the authors' JGL::JGL, penalty="fused") ----
  #      lambda1: NULL => DATA-ADAPTIVE per (cell,seed) -- lambda1_max|off-diag S_k|, down to /1e4, jgl_n_lambda1 log-spaced points = 2*max over slices of
  #      max|off-diag S_k|, down to /1000, jgl_n_lambda1 log-spaced points (50 = the shared
  #      path resolution). lambda2 = the fusion grid (the extra knob JGL alone has).
  jgl_lambda1_grid = NULL,
  jgl_n_lambda1    = 50L,
  jgl_lambda2_grid = c(0, 0.02, 0.05, 0.10, 0.20, 0.40),

  # ---- relaxed refit (per-slice unpenalized constrained MLE on the pass-1 support) ----
  jgl_refit        = TRUE,
  jgl_refit_maxit  = 500L,        # glasso backstop, SAME as the CGLasso refit (thr stays glasso's own default)
  jgl_refit_ridge_eps = 1e-3,     # ★ OURS, disclosed — pass-2 ridge on the slice covariance, so the
                                  #   given-graph MLE exists for EVERY support at EVERY n (a strictly
                                  #   PD matrix is its own PD completion). Applied INSIDE the solver
                                  #   only; the likelihood is scored on the unridged S. Off-diagonals
                                  #   stay unpenalized, so the de-biasing is untouched. 0 = exact
                                  #   unridged MLE (and the dense-end stalls). JGL's OWN key — the
                                  #   CGLasso refit has its own (cglasso_refit_ridge_eps); never shared.
                                  #   See the header of R/baselines/refit_jgl.R.
  # ---- lambda2 STRIPING (optional; scheduling only, never changes the estimator) ----
  #      Set JGL_L2_STRIPE="i/n" in the environment to make this task compute only the lambda2
  #      columns of stripe i of n, writing them to the SHARED per-point checkpoint dir. Launch the
  #      n stripes as n separate array jobs; whichever finishes LAST assembles and writes the rds.
  #      Unset => no striping => the whole (lambda1 x lambda2) grid in this task, as before.
  #      The assembled result is bit-identical either way (verified).
  jgl_l2_stripe = { s <- Sys.getenv("JGL_L2_STRIPE")
                    if (!nzchar(s)) NULL else as.integer(strsplit(s, "/", fixed = TRUE)[[1]]) },

  jgl_selector     = "AIC",       # deployed selector (Danaher 2014 §6); BIC/eBIC also recorded for every grid point
  jgl_ebic_gamma   = 0.25,        # eBIC = BIC + 4*gamma*D*log(P) (Foygel-Drton 2010)

  # ---- run/output ----
  # ---- JOB SHAPE (2026-07-25, measured against THIS cluster's scheduler) ----
  #  scontrol show config: bf_interval=60, bf_max_job_user_part=1, bf_max_job_start=5,
  #  bf_window=7200, partition_job_depth=30. MEASURED on a live 800-element 4-core array:
  #  jobs start at EXACTLY 1 per minute, one per minute, with no exception -- i.e. ARRAY TASKS
  #  DO NOT BYPASS bf_max_job_user_part. Throughput is therefore
  #        cores per minute  =  cores per JOB              (job COUNT buys nothing)
  #  so many-small-jobs is the WORST shape here: 4800 single-core tasks would take 80 h merely to
  #  START. Hence: FEW, WIDE jobs. array_by="seed" => 100 elements (1 seed x all 8 cells);
  #  ncores=1 keeps the OUTER loop over cells serial so it cannot fork-storm against the INNER
  #  mclapply in run_method_jgl, which uses SLURM_CPUS_PER_TASK to spread the 300 (l1,l2) grid
  #  points. Submit half to `share` and half to `preempt`: the 1-job/minute quota is PER PARTITION.
  #  lambda2 striping is NOT used for the same reason (it multiplies the job count 6x); the code
  #  stays in place (verified bit-identical) for environments where job count is not the binding
  #  constraint.
  array_by = "seed",              # 1 array element = 1 seed x all 8 cells
  out_dir  = here::here("results"),
  ncores   = 1L                   # OUTER loop serial; the 16 cores go to the INNER grid sweep
)
