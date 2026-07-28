# ============================================================================
# config/publication.R — the publication run spec (defines CFG).
#
#   Rscript simulation/run_comparison.R config/publication.R
#
# This file IS the "what did the paper run" record (the reproducibility artifact).
# It ONLY defines parameters; the machinery is simulation/run_comparison.R.
#
# STATUS: all four are DONE — TV (base + refit), the CGLasso baseline (base + glasso-native
# relaxed refit), the JGL baseline (base + relaxed constrained-MLE refit; m=7 cells only) and
# the tvmgm baseline (base + relaxed per-node post-lasso refit; all 16 cells).
# The generator/cells/seeds below are the paired data all methods reuse, so DO NOT change
# them when adding baselines (identical seeds+generator => bit-identical paired data per
# (cell,seed)).
#
# ⚠ COST: refit runs at EVERY lambda (user decision: the whole path, so the
# selected point genuinely differs). PROFILE ONE heaviest fit on the CLUSTER
# before launching the full grid (cluster-submit discipline), then scale.
# ============================================================================

source(here::here("config", "paired_grid.R"))   # PAIRED_SEEDS (100), PAIRED_CELLS (shared, anti-drift)

CFG <- list(
  tag = "publication_tv",

  # ---- METHODS: all four, each at its own refit point and its own native selector.
  #      They reuse the identical paired data (same seeds + generator => bit-identical
  #      per cell), so running the four config/pub_<method>.R files separately gives
  #      results identical to running this one file. ----
  methods = c("tvcglasso", "CGLasso", "JGL", "tvmgm"),   # JGL runs on the m=7 cells only

  # ---- SEEDS: 100 (JASA-grade / SUBMISSION) ; paired across methods (shared PAIRED_SEEDS) ----
  seeds = PAIRED_SEEDS,

  # ---- CELLS (EDIT to the final publication grid; these are a representative
  #      Zebrafish-grounded contB grid — m in {7,15}, P in {15,25}, n in {12,20},
  #      depth low/high). Every cell names its generator EXPLICITLY (strict
  #      dispatch: "mirrorexp_contB" = the continuous exp f; NOT the 0617 jump). ----
  cells = PAIRED_CELLS,   # shared spec (config/paired_grid.R) — identical across all method windows

  # ---- TV base-fit (engine) params ----
  #   lambda_grid_tv = NULL -> DATA-ADAPTIVE lambda_max|off-diag S_k| down to /1e4, 50 log-spaced points per (cell,seed): lambda_max = 2*max|off-diag
  #   cov(ALR Z_0)| down to /1000, so the sparse (empty) AND dense corners are reached with REAL
  #   lambda, no wasted all-empty points. tv_free_diag = TRUE (ADOPTED 2026-07-22, data-validated:
  #   free/estimated diagonal, the glasso/CGLasso/Xue-Shu-Qu standard).
  lambda_grid_tv = NULL, tv_n_lambda = 50L,   # 50 (2026-07-24, reduced from 70): ample ROC/pAUC resolution, ~29% fewer fits; ALL 4 methods use 50
  tv_q = 2L, tv_N_n = 1L, tv_init_mode = "diag", tv_sel_type = "hard",
  tv_weight_mode = "glasso", tv_max_iter = 150L, tv_sel_threshold = NULL, tv_free_diag = TRUE,

  # ---- refit (pass 2) — de-biased estimator + joint-LNM-at-Zhat BIC at every lambda/rho.
  #      SEPARATED CAPS (2026-07-22 engineering discipline): TV and CGLasso refit are INDEPENDENT
  #      (own files R/refit.R vs R/baselines/refit_cglasso_core.R) with OWN caps — never a shared key. ----
  tv_refit = TRUE, cglasso_refit = TRUE,
  tv_refit_max_outer      = 500L,   # TV cap = 500 (2026-07-24 DECISION): valid16 proved 150 too low at low-depth n<P (0/70 refit converge); free-diagonal refit is slower too. refit exits at convergence so only the slow lambdas cost more.
  cglasso_refit_max_outer = 50L,    # = Compo_glasso's own max_iter (Yuan); the refit also uses Yuan's convergence rule (relative-to-magnitude, parameter-change only)
  refit_inner_max = 10L, refit_z_align_max = 160L,
  refit_conv_tol_Z = 5e-5, refit_conv_tol_z_grad = 1e-4,

  # ---- CGLasso base-fit params (used only when "CGLasso" is in methods) ----
  cglasso_length_rholist = 50L, cglasso_option = 2L,

  # ---- JGL params (m<=7 only — run_method_jgl hard-stops m>7, so JGL covers the 8 m=7 cells).
  #      base = the authors' JGL::JGL(penalty="fused") on the centered plug-in ALR (faithful).
  #      lambda1 = NULL -> DATA-ADAPTIVE per (cell,seed): lambda1_max = 2*max over slices of
  #      max|off-diag S_k| down to /1000, 50 log-spaced points (the shared path resolution);
  #      lambda2 = the fusion grid (JGL's extra knob). refit = per-slice unpenalized constrained
  #      MLE glasso(S_k, rho=0, zero=non-edges) at glasso's OWN default tolerance (same solver
  #      settings as the CGLasso refit); a support whose MLE does not exist at this n is flagged
  #      and excluded from deployment. Output follows the SHARED shell: a FLAT `detail` row per
  #      (l1,l2) with POST-refit {neg2loglik,AIC,BIC,eBIC} + PRE-refit {*_pre} + FPR/TPR/precision/
  #      F1/df, a direct-read `deployed`, and an `estimate` payload (base/refit Theta at EVERY grid
  #      point) -- so any selector is recomputable at analysis with no re-fit. ----
  jgl_refit = TRUE, jgl_selector = "AIC", jgl_ebic_gamma = 0.25,
  jgl_lambda1_grid = NULL, jgl_n_lambda1 = 50L,
  jgl_lambda2_grid = c(0, 0.02, 0.05, 0.10, 0.20, 0.40),
  jgl_refit_maxit = 500L,
  jgl_refit_ridge_eps = 1e-3,   # ★ OURS, disclosed (as for CGLasso): pass-2 ridge eps*mean(diag(S_k))*I
                                # inside the refit SOLVER only, so the given-graph MLE exists for every
                                # support at every n; off-diagonals unpenalized, likelihood scored on the
                                # unridged S. JGL's OWN key. See R/baselines/refit_jgl.R.

  # ---- tvmgm params (all 16 cells — tvmgm handles m=7 and m=15).
  #      base = the authors' mgm::tvmgm (kernel-weighted nodewise lasso) swept over 50
  #      DATA-ADAPTIVE lambda (glmnet's OWN lambda_max, max over estpoint x node -> empty graph
  #      at the sparse end, down to lambda_max/1e4 -> full graph at the dense end, so BOTH ROC
  #      corners are real fits; 50 = the shared path resolution), edge = wadj != 0; package defaults
  #      respected (k=2, ruleReg="OR", threshold="none", scale=TRUE, lambdaGam=0.25) and the
  #      bandwidth chosen by mgm's OWN bwSelect. refit = the textbook relaxed lasso inside
  #      mgm's nodewise framework: each node refitted UNPENALISED by weighted least squares on
  #      ITS OWN directed lasso mask (from wadjNodewise, not the OR-symmetrised graph), on the
  #      same kernel weights and standardised design mgm uses, then re-selected PER NODE by
  #      mgm's OWN EBIC (nodeEst's criterion) — so pre vs refit differ ONLY in penalised-vs-
  #      unpenalised; BIC/AIC are recorded too. A support whose exact WLS does not exist is
  #      flagged and excluded (never a silent ridge). Output follows the SHARED shell: per-
  #      lambda `detail`, a direct-read `deployed`, `native` (mgm's penalised EBIC point), and
  #      an `estimate` payload carrying base wadj/wadjNodewise at every lambda plus the FULL
  #      per-(estpoint,node,lambda) grid of RSS/df/neighbours/coefficients — so any selector is
  #      recomputable at analysis with no re-fit. mgm compares on EDGE SUPPORT: it estimates no
  #      joint precision matrix, so its weights are not precision magnitudes. ----
  mgm_lambda_grid = NULL, mgm_n_lambda = 50L, mgm_lambda_min_ratio = 1e-4,
  mgm_lambdaGam = 0.25, mgm_tune = TRUE, mgm_refit = TRUE,
  mgm_bwSeq = c(0.03, 0.05, 0.1, 0.2, 0.35, 0.6, 1.2), mgm_bwFolds = 1L, mgm_bwFoldsize = 5L,

  # ---- run/output ----
  array_by = "task",                      # 1 array element = 1 (cell,seed); backfills as small jobs
  out_dir = here::here("results"),
  ncores = NULL                           # NULL => SLURM_CPUS_PER_TASK (or detectCores()-1 locally)
)
