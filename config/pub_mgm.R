# ============================================================================
# config/pub_mgm.R — tvmgm (Haslbeck & Waldorp 2020) publication run = THIS BASELINE ONLY.
#
#   Rscript simulation/run_comparison.R config/pub_mgm.R
#
# Per-method run file (the file we actually launch). It sources the SHARED paired-data spec
# and uses the SHARED tag, so it BACKFILLS results/<tag>/tvmgm/ alongside the other methods'
# folders -- running the per-method files separately is bit-identical to running
# config/publication.R once (identical seeds + generator => identical data).
#
# METHOD (see R/methods/method_tvmgm.R + R/baselines/refit_mgm.R): base = the AUTHORS' own
# mgm::tvmgm, a kernel-weighted nodewise lasso, swept over the lambda grid (edge = wadj != 0)
# -- a faithful reproduction with the package defaults respected. Refit = the textbook
# relaxed lasso inside mgm's own nodewise framework: each node is refitted UNPENALISED by
# weighted least squares on ITS OWN directed lasso mask (taken from wadjNodewise, never from
# the OR-symmetrised graph), using tvmgm's own kernel weights and standardised design, and
# re-selected PER NODE by mgm's OWN EBIC (nodeEst's criterion: a unit-variance Gaussian
# deviance + d*log(nadj) + 2*gamma*d*log(P-1), nadj = the summed kernel weights). So pre and
# refit differ ONLY in penalised-vs-unpenalised. Where the exact WLS does not exist the
# candidate is flagged and excluded -- never a silent ridge.
#
# ALL 16 cells: tvmgm handles m = 7 and m = 15 (validated end-to-end on all 16), so this is
# 16 cells x 100 seeds.
#
# OUTPUT follows the SHARED shell contract (same layout as tvcglasso / CGLasso / JGL):
#   $roc $auc $scores + $detail (per-lambda: lambda/n_edges/df/FPR/TPR/status)
#   + $deployed (the as-deployed point at the native selector, a direct read) + $native
#   + $tuning (bandwidth + every mgm knob actually used)
#   + $estimate (payload: base wadj/wadjNodewise at EVERY lambda; refit landings, de-biased
#     networks, and the FULL per-(estpoint,node,lambda) grid of RSS/df/neighbours/coefficients).
# So every selector's deployed point and the base-vs-refit comparison are recomputable at
# analysis with NO re-fit. mgm-specific: the payload holds nodewise weights (no Omega -- mgm
# estimates no joint precision matrix; no Z -- no latent layer), so mgm enters the comparison
# on EDGE SUPPORT, and its weights are not to be read as precision magnitudes.
# ============================================================================
source(here::here("config", "paired_grid.R"))   # PAIRED_SEEDS (100), PAIRED_CELLS (16) — shared spec

CFG <- list(
  tag     = "publication_tv",     # SHARED tag (same as pub_tvcglasso.R / pub_cglasso.R / pub_jgl.R)
  methods = c("tvmgm"),           # this window = the tvmgm baseline only

  # ---- paired data (NEVER redefine locally; identical seeds+generator => bit-identical
  #      data as every other method's window) ----
  seeds = PAIRED_SEEDS,           # ALL 100 (submission)
  cells = PAIRED_CELLS,           # all 16 (tvmgm runs at m = 7 and m = 15)

  # ---- tvmgm base fit (the authors' mgm::tvmgm; package defaults) ----
  #      lambda: NULL => DATA-ADAPTIVE per (cell,seed) at the SELECTED bandwidth, using glmnet's
  #      OWN lambda_max (max over estpoint x node, so the sparse end is the EMPTY graph) down to
  #      lambda_max/1e4 (glmnet's own n>p lambda.min.ratio, so the dense end is the FULL graph)
  #      -- BOTH ROC corners from real fits, no interpolation. 50 points = the shared path
  #      resolution (as CGLasso's rholist and JGL's lambda1); lambda is never compared ACROSS
  #      methods -- each traces its own ROC.
  mgm_lambda_grid = NULL, mgm_n_lambda = 50L, mgm_lambda_min_ratio = 1e-4,
  mgm_lambdaGam   = 0.25,         # EBIC gamma = mgm's package default
  mgm_tune        = TRUE,         # bandwidth by mgm's OWN bwSelect (its out-of-sample CV)
  mgm_bwSeq       = c(0.03, 0.05, 0.1, 0.2, 0.35, 0.6, 1.2),
  mgm_bwFolds     = 1L, mgm_bwFoldsize = 5L,

  # ---- relaxed refit (per-node unpenalised WLS on each node's own directed mask) ----
  mgm_refit = TRUE,               # deployed selector = mgm-native relaxed EBIC; BIC/AIC also recorded

  # ---- run/output ----
  array_by = "task",              # 1 array element = 1 (cell,seed); backfills as small jobs
  out_dir  = here::here("results"),
  ncores   = NULL                 # NULL => SLURM_CPUS_PER_TASK
)
