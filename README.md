# TVCGLasso — time-varying microbial interaction networks

Sparse time-varying precision-matrix estimation for longitudinal compositional
count data. Each edge of the precision matrix Ω(t) is a per-edge B-spline function
fit under a logistic-normal-multinomial (LNM) layer, an adaptive sliding-window
group-lasso penalty, and a de-biasing **refit** (relaxed lasso on the selected
support). This repository is the clean, publication version of the method + the
simulation comparison framework.

---

## STATUS — pre-submission snapshot

This is a **working snapshot**, not a finished release. Three things are
deliberately deferred until the manuscript is submitted and the arrangements are
confirmed with all collaborators:

| Deferred | Why | State today |
|---|---|---|
| **LICENSE** | Copyright in the code is shared with co-authors; one author should not license joint work unilaterally. | No `LICENSE` file yet, so default copyright applies. A license will be added at submission. |
| **Zebrafish source data** | The tables belong to the study of Gaulke et al. (2019). | Held back; the cleaning **code** is here and regenerates everything from them. |
| **Real-data results** | Not front-running our own submission. | `analysis/` ships the simulation leg; the zebrafish figure/table outputs are held. |

Everything needed to reproduce the **simulation** results is present.

`R/baselines/CompoGlasso.R` is third-party code; see [`PROVENANCE.md`](PROVENANCE.md)
for its origin, its licensing, and the complete list of our changes to it.

---

## What's here

```
tvcglasso/
├── R/
│   ├── tvcglasso.R          # TV engine: main_function_final() + tv_warm_path()  (+ G_beta_Rcpp)
│   ├── tvcglasso.cpp        # the one compiled helper (G_beta_Rcpp)
│   ├── refit.R              # relaxed refit: refit_fixed_support() + joint-LNM-at-Zhat BIC
│   ├── roc_utils.R          # method output contract + ROC helpers
│   ├── methods/
│   │   ├── method_tv.R      # OUR method: base fit + refit
│   │   ├── method_cglasso.R # static CGLasso control  (base + relaxed refit, ridge-stabilised)
│   │   ├── method_tvmgm.R   # kernel tvmgm baseline   (base + per-node relaxed refit)
│   │   └── method_jgl.R     # fused JGL baseline, m<=7 (base + given-graph relaxed refit)
│   └── baselines/           # vendored CGLasso solver (NOT ours — see its header) + our refits
├── simulation/
│   ├── generators.R         # generate_data() strict dispatch + the two generators
│   └── run_comparison.R     # MASTER driver: run_comparison(CFG)
├── config/
│   ├── paired_grid.R        # THE shared design: 16 cells x seeds 1:100 (all methods)
│   ├── publication.R        # the authoritative 4-method run spec
│   └── pub_<method>.R       # one per method — what we actually launch
├── analysis/                # regenerate every paper figure + table  (see analysis/README.md)
│   ├── make.R, Makefile     # `Rscript analysis/make.R` rebuilds everything
│   ├── R/, scripts/         # shared helpers + one script per paper output
│   ├── cache/               # lets the figures rebuild WITHOUT the ~17 GB results/
│   └── figures/, tables/    # rendered outputs
├── data/
│   ├── prepare_zebrafish.R  # raw -> cleaned; the only script that writes data/
│   ├── README.md            # provenance, citation, cleaning recipe
│   └── (raw/ + derived .rds are NOT committed — see STATUS above)
├── setup.R                  # install dependencies
├── PROVENANCE.md            # how the stored result files map to this code
└── .here                    # here::here() project-root sentinel
```

## Install & run

```r
Rscript setup.R                                             # install deps
Rscript simulation/run_comparison.R                         # tiny TV+refit smoke (results/smoke/)
Rscript simulation/run_comparison.R config/pub_tvcglasso.R  # one method's publication run
Rscript simulation/run_comparison.R config/publication.R    # all four methods
Rscript analysis/make.R                                     # rebuild every figure + table
```

Run from the **repo root** (so `here::here()` anchors on `.here`). The full grid is
a cluster-scale job (16 cells x 100 seeds x 4 methods); `run_comparison.R` reads
`SLURM_ARRAY_TASK_ID` and supports `array_by = "seed"` or `"task"`, and every
`(cell, seed)` is checkpointed and skipped on resubmit. See `PROVENANCE.md` for the
run inventory.

## The method (`method_tv.R`)

The publication estimator is **TV base fit + refit**, controlled by `cfg$tv_refit`:

- **Pass 1** — `tv_warm_path()`: a high→low-λ warm-start continuation that SELECTS
  the support at each λ. Its λ-sweep ROC (edge = `Ω̂ ≠ 0`) is the method's ROC.
- **Pass 2** — `refit_fixed_support()` at **every λ**: re-estimates the active
  off-diagonal magnitudes **and** the latent Z on the FIXED support, unpenalized
  (a genuinely reduced parameterization — not gradient-masking). Refit does not
  change which edges are nonzero, so it does not change the ROC; its role is the
  **de-biased deployed operating point** (refit-BIC) + magnitude accuracy.

Model-selection likelihood = the **joint LNM density at the inferred latent Ẑ**
(multinomial + Gaussian-graphical layers), via `refit_joint_nll_average` →
`refit_information_criteria`. It is **not** the raw-ALR Z₀ Gaussian-only likelihood.

## Saved output (run and plotting are decoupled)

Each `(cell,seed)` writes `results/<tag>/<method>/cellXX_seedYYY.rds` holding the full
intermediate, so **any figure / selector is computable later without re-fitting**:

- `roc`, `auc`, `scores` (per-slice entry-λ), `detail` (per-λ df, edges, min_eig,
  λ-aligned FPR/TPR, and the joint-LNM IC components),
- `deployed` — the as-deployed operating point at the method's own native selector,
- `estimate` — the per-λ fitted objects for both passes (β / Ω / Z as the method has them),
- `provenance` — resolved config, git commit, R version, data fingerprint.

Resume-safe: a completed `(cell,seed)` is skipped on resubmit; TV and refit
checkpoint per-λ (preemption loses ≤ 1 λ).

## Generators (`generators.R`)

- `mirrorexp_contB` — the main-text "exp f": continuous mirror-exp, each active edge
  decays to **exactly 0** at t = 0.5 (kink, no jump) → a genuine temporal zero
  region. Zero-point-identification generator.
- `hetrate_contBfast` — continuous heterogeneous-rate option (fast edges use the
  contB zero-region shape; slow edges stay always-on).

`generate_data(cell, seed)` is a **strict dispatcher**: an unknown `cell$generator`
is a hard error — there is no silent fallback (the discontinuous 0617 mirror-exp is
deliberately unavailable).

## Real data (`data/`)

The application is the OSU zebrafish gut-microbiome time course of
[Gaulke et al. (2019)](https://doi.org/10.1186/s40168-019-0622-9): 207 fish
sampled destructively over **7 irregularly spaced days** (7, 10, 21, 30, 43, 59,
86 post-exposure), split by realized parasite burden into **infected** (n = 81)
and **not infected** (n = 126).

**The data files themselves are not committed** (see STATUS above). The cleaning
code and its documentation are, and reproduce the cleaned dataset exactly once the
source tables are in `data/raw/`:

```r
Rscript data/prepare_zebrafish.R      # from the repo root
```

which writes `data/zebrafish_clean.rds` (genus counts plus per-sample and
per-taxon metadata) and `data/zebrafish_real_depths.rds`. Per-P model input is cut
on demand with `zeb_slices(clean, P, group)` rather than cached. The script ends in
an **assertion gate over 20 landmark quantities** — per-day sample counts, taxon
counts at each prevalence threshold, depth quantiles, the P = 15 node set — so a
drift in the raw data or the recipe is a hard failure rather than a silently
different dataset.

Read **[`data/README.md`](data/README.md)** before touching the raw tables: it
records the provenance and citation, the cleaning recipe step by step, and three
ways these files can be *silently* misparsed.

## Method status

All four methods are complete: base fit + relaxed refit + full intermediate recording.

- **TVCGLasso (ours)** — per-edge B-spline Ω(t) under the LNM layer; refit on the
  fixed support with a genuinely reduced parameterization.
- **CGLasso** — static per-slice control (the vendored Yuan solver); refit is the
  given-graph Gaussian MLE, ridge-stabilised on the slice covariance so it exists at
  n < P. The ridge is **our addition** and is disclosed in `refit_cglasso_core.R`.
- **tvmgm** — kernel nodewise baseline; refit is the textbook per-node OLS-post-lasso
  with a per-node relaxed-EBIC reselection (mgm's paper-native selector).
- **JGL** — fused baseline, m ≤ 7 only (its ADMM does not scale); refit is the
  unpenalized given-graph MLE per slice, likewise ridge-stabilised.

Each method keeps its **own** native selector rather than a single shared criterion:
TVCGLasso → refit-BIC, CGLasso → BIC, tvmgm → per-node EBIC, JGL → AIC. Methods are
compared on **deployed graphs and ROC**, never on IC values across methods.

## Provenance / correctness

- The engine `R/tvcglasso.R` is a dead-code-free consolidation of the project's
  canonical chain (`main_function_final_0624.R` → `main_function_pd_0605.R` →
  `main_function_05-20.R` + `main.cpp`), **verified bit-identical** to it
  (`max|Δβ| = 0` on single-λ and warm-path fits).
- ALR zero-handling = Yuan `z_hat_offset` option 2 (proportional pseudocount),
  pooled across slices. `set.seed(7)` is kept but the caller's RNG is
  saved/restored on exit; the reported Ω̂(t_k) is PD on every slice or the fit
  errors loudly.
- `R/baselines/CompoGlasso.R` is **third-party code** (Tian et al. 2023). Its header
  lists every difference from upstream. See `PROVENANCE.md`.
