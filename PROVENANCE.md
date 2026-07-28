# Provenance

How the stored result files relate to the code in this repository, and what is
third-party.

## 1. The git hash stamped in the result files

Every `results/<tag>/<method>/cellXX_seedYYY.rds` carries
`provenance$git`, `provenance$R` and `provenance$data_fingerprint`. The stored
runs were produced in a **private development repository**, and their stamp is
`d6a9079`.

**That hash does not exist in this repository, and it should not be used as the
code reference.** Two reasons, both checkable:

1. `d6a9079` predates several files the runs actually required. `git ls-tree -r d6a9079`
   in the development repo contains none of `R/baselines/refit_mgm.R`,
   `config/pub_mgm.R`, `config/pub_jgl.R` or `config/pub_cglasso.R` — yet
   `simulation/run_comparison.R` sources `refit_mgm.R` unconditionally. A checkout of
   `d6a9079` therefore could not have run tvmgm at all.
2. The CGLasso and JGL refits changed substantially after that commit (the disclosed
   pass-2 ridge in `refit_jgl.R` is entirely post-`d6a9079`).

**The initial commit of this repository is the authoritative code reference** for
every reported number. This is stated rather than hidden because a hash that
resolves to a tree which could not have produced the results is worse than no hash.

## 2. Which config produced which result set

| Result directory | Config | Note |
|---|---|---|
| `results/tvcglasso/` | `config/pub_tvcglasso.R` | 16 cells x 100 seeds |
| `results/CGLasso/` | `config/pub_cglasso.R` | 16 cells x 100 seeds |
| `results/tvmgm/` | `config/pub_mgm.R` | 16 cells x 100 seeds |
| `results/JGL/` | `config/pub_jgl.R` | m = 7 cells only (8 x 100); JGL's ADMM does not scale past m = 7 |

All four share `config/paired_grid.R`, so a given `(cell, seed)` yields
**bit-identical data** for every method. This is verified after the fact through
`provenance$data_fingerprint` (the sum of the generated counts), not assumed.

**CGLasso λ range.** The reported CGLasso leg uses the package default
`cglasso_rho_lo_div = 1000`, i.e. `rho_min = rho_max / 1000` — Yuan's own recipe.
Neither `config/publication.R` nor `config/pub_cglasso.R` overrides it, so running
either reproduces the reported leg. A deeper `/10000` variant was explored and
**rejected**: pushing rho below that floor drives Ω large enough that the multinomial
`exp(z)` overflows inside the vendored Newton solver, which lost 107 of 1600 units to
non-random failures concentrated in the hardest cells — a loss that would have
systematically flattered CGLasso exactly where it is weakest.

## 3. Third-party code

`R/baselines/CompoGlasso.R` is **not ours**. It is the solver accompanying

> Tian, Jiang, Hammer, Sharpton & Jiang (2023). *Compositional Graphical Lasso
> Resolves the Impact of Parasitic Infection on Gut Microbial Interaction Networks
> in a Zebrafish Model.* JASA 118(544):1500–1514. doi:10.1080/01621459.2022.2164287

Upstream: <https://github.com/yuanjiang-osu/Comp-gLASSO-JASA>. An earlier release of
the same codebase is published by the same author **under the MIT License** at
<https://github.com/yuanjiang-osu/Comp-gLASSO>.

**Relationship between the two upstream releases** (measured by diff, comments
stripped). The five functions vendored here — `z_hat_offset`, `obj`, `NR`, `NR_para`
and the alternating `Compo_glasso` loop — are all present in the MIT-licensed 2021
release; the JASA release is a refinement of it. Restricted to the vendored
functions, the JASA release differs from the MIT release in exactly four places:

| Delta | JASA release |
|---|---|
| `z_hat_offset` | adds an `option == 0` branch |
| `Compo_glasso` | `max_iter` default 100 → 50 |
| parameter names | `ratio.z` / `ratio.O` → `z_ratio` / `O_ratio` |
| robustness | two `is.na` guards inside the iteration |

The remainder is the MIT-licensed code. We record this because attribution should be
precise about which grant covers what, not because the distinction has ever been in
dispute — the copyright holder is a co-author of the present work.

**Our changes to the vendored file** are listed exhaustively in its header, verified
by diff: three functions renamed with a `cg_` prefix to avoid a namespace collision
with our own engine, one `nr_max_iter` termination guard (the only semantic change),
and `generate_cov` omitted. Nothing else was altered.

Everything else under `R/`, `simulation/`, `config/`, `analysis/` and `data/` is ours.

## 4. Real data

The zebrafish tables of Gaulke et al. (2019) are not committed here — see the STATUS
section of the README. `data/prepare_zebrafish.R` regenerates every derived file from
them and ends in an assertion gate over 18 landmark quantities, so a drift in either
the raw data or the recipe is a hard failure rather than a silently different dataset.

## 5. Reproducibility caveat worth knowing

The data generator is **not** reproducible across R versions: identical cell
parameters and seed give `data_fingerprint = 63718` under R 4.4.3 but `64213` under
R 4.3.3. The paired-data guarantee across methods holds because every method ran on
the same machine and R version. `provenance$R` records it per result file.
