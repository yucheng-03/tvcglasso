# ============================================================================
# config/pub_cglasso.R — CGLasso (Tian et al. 2023) publication run = THIS BASELINE ONLY.
#
#   Rscript simulation/run_comparison.R config/pub_cglasso.R
#
# Per-method run file (the file we actually launch). It sources the SHARED paired-data
# spec and uses the SHARED tag, so it BACKFILLS results/<tag>/CGLasso/ alongside the
# other methods' folders -- running the per-method files separately is bit-identical to
# running config/publication.R once (identical seeds + generator => identical data).
#
# METHOD (see R/baselines/): base = Yuan's vendored Compo_glasso fit INDEPENDENTLY per
# time slice (the static, no-time-awareness control), on a PER-SLICE data-adaptive rho
# grid (his Functions_ROC.R recipe applied per slice). Refit = relaxed lasso on the fixed
# pass-1 support: the glasso-native constrained MLE glasso(S_k, rho=0, zero=non-edges)
# (Dempster 1972 / ESL 2e Alg 17.1) alternating with Yuan's OWN cg_NR latent-Z update,
# under Compo_glasso's OWN convergence rule. Deployed point = per-slice independent BIC
# (joint-LNM-at-Zhat) over the refit path.
# ============================================================================
source(here::here("config", "paired_grid.R"))   # PAIRED_SEEDS (100), PAIRED_CELLS (16) — shared spec

CFG <- list(
  tag     = "publication_tv",     # SHARED tag (same as pub_tvcglasso.R) -> results/publication_tv/CGLasso/
  methods = c("CGLasso"),         # this window = the CGLasso baseline only

  # ---- paired data (NEVER redefine locally; identical seeds+generator => bit-identical
  #      data as every other method's window) ----
  seeds = PAIRED_SEEDS,           # ALL 100 (submission)
  cells = PAIRED_CELLS,           # 16 cells = P{15,25} x n{12,20} x depth{low,high} x m{7,15}

  # ---- CGLasso base fit (Compo_glasso, per slice) ----
  cglasso_length_rholist = 50L,   # 50-point path (2026-07-24 decision: all 4 methods use 50)
  cglasso_option         = 2L,    # proportional pooled pseudocount (Compo_glasso's default option)

  # ---- relaxed refit (glasso-native constrained MLE + cg_NR latent Z) ----
  #      Convergence follows Compo_glasso's OWN rule (relative to the initial magnitude,
  #      parameter-change only, no gradient test) -- see R/baselines/refit_cglasso_core.R.
  cglasso_refit           = TRUE,
  cglasso_refit_max_outer = 50L,    # = Compo_glasso's own max_iter (Yuan); its convergence rule is Yuan's too
  # ★ OURS, DISCLOSED (2026-07-25) — ridge on the pass-2 slice covariance: the refit solves the
  # given-graph MLE on S_k + eps*mean(diag(S_k))*I instead of S_k, so the constrained MLE exists
  # for every graph and every n (a PD matrix is its own PD completion; Dempster 1972 Biometrics
  # 28:157-175, Uhler 2012 Ann.Statist. 40(1) Thm 2.1). Without it the refit does NOT terminate at
  # either end of the rho path when n<P. Pinned HERE rather than left to the code default so the
  # value lands in every result's provenance. Validated on the 12 previously-stalled cells:
  # 12 x 50 refits, all converged, 0 non-existent, max|Omega_post| in [17,41].
  # See CLAUDE.md, the "CGLasso REFIT STALLS" block.
  cglasso_refit_ridge_eps = 0.01,
  refit_z_align_max       = 160L,   # cap on the pre-point latent-Z alignment sweeps
  refit_inner_max         = 10L,    # (inert for CGLasso: the Omega solve is delegated to glasso)
  # Fallback tolerances, used only if the Yuan rule is disabled (yuan_conv = FALSE):
  refit_conv_tol_Z = 5e-5, refit_conv_tol_z_grad = 1e-4,

  # ---- run/output ----
  array_by = "task",              # 1 array element = 1 (cell,seed); backfills as small jobs
  out_dir  = here::here("results"),
  ncores   = NULL                 # NULL => SLURM_CPUS_PER_TASK
)
