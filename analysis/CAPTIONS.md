# Figure and table captions — drafts

Each caption below carries the things the figure itself does **not** show. Some
of these are Taylor & Francis requirements (a figure that encodes data by
colour, pattern or symbol must have an explanatory key), and some are this
project's own disclosure rules. The sentences marked **[REQUIRED]** should not
be dropped: each one closes a specific way the figure could otherwise be
misread.

Numbers marked `<...>` must be filled from the tables, not from memory.

---

## Figure 1 — `fig01_roc_grid.eps`

> **Figure 1.** Seed-averaged ROC curves for the four methods on the 16
> simulation cells. Columns give the network size and per-slice sample size
> (*P*, *n*); rows give the number of time points *m* and the sequencing-depth
> regime. Curves are averaged over 100 replicates by vertical averaging on a
> common false-positive-rate grid. **[REQUIRED]** Each curve is drawn only over
> the false-positive range its own λ-path actually reaches, and stops there; no
> curve is extrapolated to (1, 1). The differences in curve length are
> therefore differences in **coverage, not in quality** — a method whose path
> stops early has not been penalised for it here. Per-curve coverage is
> reported in `tab02_roc_coverage.csv`. **[REQUIRED]** JGL appears only in the
> *m* = 7 rows: its fused-penalty ADMM does not scale to *m* = 15, so it is
> inapplicable rather than missing. The diagonal marks chance performance.

*Not in the caption, on purpose:* no AUC. The stored `auc` integrates over
`[0, max FPR]` only, so it systematically under-counts exactly the methods
whose curves stop short; quoting it would penalise them for the disclosure.

---

## Figure 2 — `fig02_deployed_m07.eps`, `fig02_deployed_m15.eps`

> **Figure 2.** Recall, precision and F1 at each method's **own** deployed
> operating point, for the *m* = <7 | 15> cells. Bars are means over 100
> simulation replicates; **[REQUIRED]** whiskers span one standard deviation
> across replicates, truncated at zero — they describe how much a single
> replicate varies, not the precision of the mean. **[REQUIRED]** Each method
> is tuned by its own native selector with no access to the truth (tvcglasso:
> refit-BIC; CGLasso: per-slice BIC; JGL: refit-AIC; tvmgm: per-node relaxed
> EBIC); the four criteria are on different scales and are never compared
> numerically — what is compared is the deployed graph. Cells where a method
> contributed fewer than 100 replicates are listed with their realised counts
> in `tab01_deployed.csv`.

*If the m = 15 panel is shown:* add "JGL is absent because it does not scale to
*m* = 15."

---

## Figure 3 — `fig03_edge_shapes.eps`

> **Figure 3.** The two true edge-trajectory shapes of the simulation
> generator. Every active edge takes the early form (a) or its mirror image,
> the late form (b), **[REQUIRED]** each with probability one half. Shading
> marks the half of the interval on which the edge is active; on the other half
> the edge is **exactly** zero, giving each edge a genuine temporal zero
> region. **[REQUIRED]** The trajectory is continuous at *t* = 1/2 but not
> differentiable there: it reaches zero with non-zero slope and then stops
> (C⁰, not C¹). Open circles mark the *m* = 7 sampling times, so the kink is
> observed at exactly one of them. Amplitude *A* = 0.5, decay rate 6.

---

## Figure 4 (diagnostic, `analysis/diagnostics/`) — not for the paper as drawn

> Each method's ROC curve with its own deployed operating point marked.
> **[REQUIRED if ever used]** The plotted point is the mean over replicates of
> points that each lie on their own replicate's curve; because the ROC is
> concave, averaging points that are spread along it places the mean *below*
> the averaged curve. The vertical gap is therefore an averaging artefact, not
> a failure of the method to reach its own curve, and it grows with how
> unstable the selected λ is (correlation between the gap and the standard
> deviation of the selected λ index across the 16 cells: 0.95). **[REQUIRED]**
> Separately, and genuinely: the CGLasso point need not lie on its curve at
> all, because it deploys a per-slice amalgam — each time slice minimises its
> own BIC — while the curve uses one ρ for all slices.

---

## Figure 5 — `fig05_zeb_clusters_P<P>.eps`

> **Figure 5.** Estimated edge trajectories Ω̂*ᵢⱼ*(*t*) in the zebrafish gut
> microbiome, grouped by *k*-means, shown separately for infected and
> uninfected fish (*P* = <15 | 25> genera, <7-14 | 15-23> fish per sampling
> day). Grey lines are individual edges; the black line is the cluster mean,
> with the seven sampling days marked. **[REQUIRED]** *k*-means is run on
> **z-scored** trajectories, so edges are grouped by the shape of their time
> course and not by their magnitude, while the black line is the **raw**
> (un-standardised) mean of the cluster's members — the vertical axis is
> therefore on the precision-matrix scale. The number of clusters is chosen per
> group by average silhouette width; clusters are ordered by when the cluster
> mean peaks. **[REQUIRED]** Magnitudes come from the relaxed refit at the
> BIC-selected penalty; the support and the selected penalty come from the
> penalised first pass.
>
> **[REQUIRED — the honest limitation]** This figure describes the fitted
> networks. It is not evidence of an infection effect: several edges involve
> genera that are absent on some sampling days, whose additive-log-ratio value
> is then determined by the pseudocount rather than by the data, and a
> within-group bootstrap cannot distinguish such a deterministic artefact from
> a reproducible biological signal.

---

## Table 1 — `tab01_deployed.csv`

Column meanings, for the table note: `*_mean` / `*_sd` / `*_se` are the mean,
standard deviation and standard error across replicates; `n_seed` is the number
of replicates actually contributing (not always 100 — see Figure 2's note);
`stage` is `refit` (deployed) or `pre` (the same selector before the relaxed
refit); `selector` names the criterion each method used.

**[REQUIRED]** The information criteria behind `selector` are not comparable
across methods: tvcglasso and CGLasso use the joint logistic-normal-multinomial
deviance at the inferred latent Ẑ on raw counts, JGL a Gaussian deviance on the
centred plug-in ALR, and tvmgm a unit-variance nodewise EBIC. They are used
only to pick each method's own operating point.

---

## Standing rules these captions encode

1. No full-range AUC, and no pAUC, anywhere.
2. Curves stop at real operating points; corner gaps are disclosed, never
   interpolated.
3. Any oracle / truth-aware quantity is labelled as such and is
   simulation-only.
4. mgm is compared on edge support only — never in a precision-magnitude table,
   because it estimates no joint precision matrix.
5. Where a method's deployed point is an amalgam (CGLasso per slice, tvmgm per
   node, JGL across λ₂), say so wherever the point is plotted next to a curve.
