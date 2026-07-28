# analysis/ — regenerating every figure and table in the paper

```
Rscript analysis/make.R        # run from the repository root
```

That one command rebuilds every numbered figure and table from the stored
per-`(method, cell, seed)` result files. It needs no cluster, no simulation
re-run, and no packages beyond those in `setup.R`.

---

## What each script produces

| script | paper output |
|---|---|
| `scripts/fig01_roc_grid.R` | **Figure 1** `figures/fig01_roc_grid.eps` — seed-averaged ROC of all four methods, one panel per cell; also `tables/tab02_roc_coverage.csv` |
| `scripts/fig02_deployed_bars.R` | **Figure 2** `figures/fig02_deployed_m07.eps` — recall / precision / F1 at each method's own native selector; also `tables/tab01_deployed.csv` |
| `scripts/fig03_generators.R` | **Figure 3** `figures/fig03_edge_shapes.eps` — the two true edge-trajectory shapes (early / late) of the simulation generator |
| `scripts/tab03_provenance.R` | `tables/tab03_provenance.csv`, `tables/tab03_pairing_check.csv` — commit, R version, unit counts, and the paired-data check |
| `scripts/fig05_zebrafish_clusters.R` | **Figure 5** `figures/fig05_zeb_clusters_P<P>.eps` — real-data edge trajectories grouped by shape; also `tables/tab05_zeb_clusters_P<P>.csv`. Needs the fits from `realdata/` (below), so `make.R` does not run it |

Draft captions for every figure, including the disclosures each one must carry,
are in [`CAPTIONS.md`](CAPTIONS.md).

Each script is self-contained and can be run alone; `make.R` just runs them in
order. `diagnostics/` holds analyses that are **not** paper outputs and that
`make.R` deliberately does not run — currently `fig04_deployed_on_roc.R`, which
overlays each method's deployed operating point on its own ROC curve to separate
estimator quality from tuning-rule quality. Run it by hand if you want it. Set `TVCG_REFRESH=TRUE` to bypass `cache/` and re-extract from the `.rds`
files (deleting `cache/` and re-running must reproduce every figure exactly).

## The real-data application

`realdata/fit_zebrafish.R` fits the publication estimator to the zebrafish time
course, one fit per (P, fish group), and writes `results/zebrafish/`. `make.R`
does **not** run it: it is a multi-hour cluster job, not part of the
regenerate-the-figures path. `scripts/fig05_zebrafish_clusters.R` then reads
those fits and needs nothing else.

```
Rscript analysis/realdata/fit_zebrafish.R          # all four (P, group) fits
TVCG_ZEB_P=15 TVCG_ZEB_GROUP=infected Rscript ...  # one cell
sbatch cluster/submit_zebrafish.sh                 # the four as a SLURM array
```

Two things about this leg are worth knowing before reading its output.

**It is the n < P regime, and that is the point.** There are 15 (or 25) genera
but only 7-14 fish per day in the infected group and 15-23 in the uninfected
one. Every quantity that pools days therefore weights day *k* by *n_k*/*N*, not
by 1/*m* — the objective, the gradient, the information criteria and the
convergence test all do. This is also exactly the regime the simulation flags
as the one where likelihood-based model selection becomes unreliable, so the
simulation's caveat about the deployed operating point applies here directly
rather than by analogy.

**The data is not in this repository.** `data/prepare_zebrafish.R` and
`data/README.md` are committed and reproduce the cleaned dataset exactly from
the source tables, but the tables themselves belong to the original study and
are held back pending redistribution clearance (see `.gitignore`). The fits are
derived objects and are small; where they live is set by `ZEB_FITS_DIR` in
`config.R`.

## Where the inputs come from

`config.R` is the only file that knows machine-specific paths. It expects

```
results/tvcglasso/cellCC_seedSSS.rds
results/CGLasso/  ...
results/tvmgm/    ...
results/JGL/      ...          (m = 7 cells only, see below)
```

and each root can be redirected with an environment variable
(`TVCG_RESULTS_DIR`, `CGLASSO_RESULTS_DIR`, `TVMGM_RESULTS_DIR`,
`JGL_RESULTS_DIR`) without editing any file.

**Two tiers.** Tier 1 is this directory: minutes, on a laptop, from the stored
results. Tier 2 is the simulation campaign that produced those results
(`simulation/run_comparison.R` + `cluster/`), which is a multi-hundred
core-hour SLURM job; the per-`(cell, seed)` design means any single unit can be
recomputed on a laptop to spot-check that Tier 2 reproduces Tier 1's inputs.

**JGL exists only for cells 1–8.** Its fused-penalty ADMM is infeasible at
m = 15, and `method_jgl.R` hard-stops there. That is a property of the method,
not missing data, and the figures show it as an absent curve rather than a gap.

## The conventions, and why they are what they are

These are enforced in code, in one place each, so the two figures cannot drift.

**Visual identity — `R/00_aes.R`.** One method → (colour, line type, fill,
symbol) map. Every series is encoded by colour **and** line type, because
Taylor & Francis print figures in black and white and forbid encoding series
identity by hue alone. The palette is Paul Tol's high-contrast scheme, whose
four inks are separated in luminance (CIE L\* 0.0 / 29.2 / 49.5 / 72.5,
minimum adjacent gap 20.3) as well as in hue, and which is stable under
deuteranopia and protanopia.

**ROC averaging — `R/02_roc.R`.** One convention: ties in FPR collapsed to the
upper envelope; `approx(rule = 1)` so no seed is extrapolated past its own
data; the mean drawn only where at least half the seeds still have real points;
(0,0) anchored because the empty graph is a real operating point at the top of
every method's data-adaptive λ path; **(1,1) never anchored**. Curves therefore
stop where the estimator stops. This is published practice in this literature,
not an omission — the caption states that the difference in curve length is
**coverage, not quality**, and `tables/tab02_roc_coverage.csv` quantifies it.

**No AUC is reported.** The `auc` field stored in the result files integrates
over `[0, maxFPR]` with no (1,1) anchor, so it systematically under-counts
exactly those methods whose curves honestly stop short; quoting it would
penalise them for the disclosure. The quantitative comparison is the
as-deployed table.

**Operating points — `R/03_deployed.R`.** The four methods do not share a row
shape (`deployed` is flat with 10, 11 and 7 fields for tvcglasso, JGL and
tvmgm, and *nested* under `$pre`/`$post` for CGLasso), and "the point before the
refit" is a different object in each. One adapter normalises all of them.
Precision, where a method does not store it, is derived exactly from
(FPR, TPR) and the truth's class sizes rather than inverted from F1.

Each method is tuned by **its own** native selector (tvcglasso refit-BIC,
CGLasso per-slice BIC, JGL refit-AIC, tvmgm per-node relaxed EBIC). The
criteria are on four different scales and are never compared numerically; what
is compared is the deployed graph. For CGLasso, tvmgm and JGL the deployed
point need not lie on the plotted curve (per-slice, per-node and per-λ₂
amalgams respectively); `on_curve` records this and the figures do not imply
otherwise.

**Artwork compliance — `R/04_check.R`.** Every figure is checked against the
Taylor & Francis specification *as produced*, because the common failure modes
are silent: an unresolved font family is replaced without a warning, a panel
layout rescales text below the minimum without a warning, and a stray
transparent object is either dropped or rasterised without a warning. The gate
verifies EPS output, embedded fonts, standard font families only, no Type 3
fonts, no rasterised regions, page size, and a minimum stroke of 0.5 pt (which
satisfies both the 0.3 pt journals rule and the 0.5 pt books rule). Figures are
drawn at final printed size (7.17 in full width) so that nothing is rescaled
afterwards — rescaling shrinks line weights and type proportionally and is the
usual way a compliant figure quietly stops being one.

`.eps` is the submission file. The `.png` beside it is a screen preview and is
never submitted.

## Requirements

R (tested on 4.3.3 and 4.4.3) with `here`; base graphics only — no plotting
package is added on top of `setup.R`. `ghostscript` is optional and used only
by the compliance checker's independent font probe.
