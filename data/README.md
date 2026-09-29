# Data

The real-data application: the OSU zebrafish gut-microbiome time course.

> ### The data files are not distributed in this repository (yet)
>
> The source tables belong to the original study. We have not yet confirmed in
> writing that we may redistribute them, so they and everything derived from
> them are **deliberately excluded from version control** (see `.gitignore`).
> What is committed is our **cleaning code** and this documentation, which
> together reproduce the cleaned dataset exactly once you have the source
> tables. This is a temporary hold, not a permanent policy — see *Provenance*
> below for where the data comes from and how to obtain it.

```
data/
├── prepare_zebrafish.R         raw -> cleaned; the only script that writes here   [committed]
├── README.md                   this file                                          [committed]
├── raw/                        the source tables, exactly as published            [not committed]
│   ├── asv.tab                 237 samples x 841 ASV counts
│   ├── tax.tab                 841 ASVs x 6 taxonomic ranks
│   └── metadata.tab            237 samples x 18 metadata fields
├── zebrafish_clean.rds         the cleaned dataset we analyse                     [not committed]
└── zebrafish_real_depths.rds   207 library sizes; the pool the simulation draws
                                from for depth_mode = "real" (generators.R)
```

Given the three tables in `raw/`, everything else is rebuilt by

```sh
Rscript data/prepare_zebrafish.R      # from the repo root
```

## Provenance

> Gaulke CA, Martins ML, Watral VG, Humphreys IR, Spagnoli ST, Kent ML, Sharpton TJ (2019).
> A longitudinal assessment of host-microbe-parasite interactions resolves the zebrafish gut
> microbiome's link to *Pseudocapillaria tomentosa* infection and pathology.
> *Microbiome* **7**:10. doi:[10.1186/s40168-019-0622-9](https://doi.org/10.1186/s40168-019-0622-9)

Zebrafish (*Danio rerio*) were exposed to the intestinal nematode *P. tomentosa*
and sampled destructively over 86 days. The article is open access under
[CC BY 4.0](http://creativecommons.org/licenses/by/4.0/); the raw sequence reads
are in the NCBI SRA under BioProjects
[PRJNA472413](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA472413) and
[PRJNA472775](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA472775).

The three tables in `raw/` are the ones released with the compositional graphical
lasso paper, which analyses this same cohort and is the static baseline in our
comparison:

> Tian C, Jiang D, Hammer A, Sharpton TJ, Jiang Y (2023). Compositional graphical lasso
> resolves the impact of parasitic infection on gut microbial interaction networks in a
> zebrafish model. *JASA* **118**(543):1500–1514.
> doi:[10.1080/01621459.2022.2164287](https://doi.org/10.1080/01621459.2022.2164287)

Our copies are those tables unmodified. Checksums, so that whoever supplies the
files can confirm they are the ones this code was written against (they are also
recorded in `zebrafish_clean.rds$provenance`, making a mismatch detectable):

| file | bytes | MD5 |
|---|---|---|
| `asv.tab` | 436,438 | `ed5cffdd6ec09adbf501381dde07742b` |
| `tax.tab` | 95,555 | `295782c763b375712b064b5d2ab5013c` |
| `metadata.tab` | 32,529 | `f357dc3d5ae9bb9db217a899bb76d09f` |

## The cleaned dataset

`zebrafish_clean.rds` is a list with `counts` (207 x 260 integer genus counts),
`samples`, `taxa`, `days`, `dropped` and `provenance`. The recipe:

1. **ASVs to genus** — 841 ASVs collapse to 260 columns. The 298 ASVs unassigned
   at genus level are pooled into one column, `NONE`, following the published
   preprocessing. Their reads are kept, so each row still sums to the true
   library size (asserted).
2. **Drop the 30 fish sampled at `DaysPE == 0`** — pre-exposure animals, not an
   independent draw from the post-exposure course. 207 fish over **7 days**
   (7, 10, 21, 30, 43, 59, 86) remain.
3. **Group by realized worm burden** — `infected` = `Total > 0` (n = 81),
   `not_infected` otherwise (n = 126). `Total` is the parasite count, not the
   assigned arm: 24 `Exposed` fish ended with `Total == 0` and 2 `Unexposed`
   fish with `Total > 0`, and the grouping follows the burden.
4. **Depth** = total reads per sample over all genera.
5. **Prevalence / relative abundance** over the 207 retained samples.

`counts` is the single source of truth; per-P analysis slices are deliberately
**not** cached, since they are a re-cut of these same columns. To get model input:

```r
source("data/prepare_zebrafish.R")               # loads the functions, runs nothing
clean <- readRDS("data/zebrafish_clean.rds")
s <- zeb_slices(clean, P = 15, group = "infected")
```

which returns one raw-count matrix per day, `n_k x (P+1)`, with the **ALR
reference taxon in the last column** (so `Z` is `n_k x P` and `Omega` is `P x P`),
plus the day values and their positions on [0,1] at the **real irregular
spacing**. Nodes are ranked by prevalence, ties broken explicitly by relative
abundance then alphabetically — without a stated rule the P = 25 boundary would
be decided by the order genera happen to appear in `tax.tab`.

**`NONE` is never a node; by default it is the ALR denominator.** It pools ASVs
unassigned at genus level across unrelated lineages, so an edge to it would have
no biological reading. It stays in `counts` (the reads are real, and depth must
be the true library size) and, as in the published preprocessing, serves as the
reference taxon: the default `reference = "NONE"` puts it in the last column and
makes the `P` most prevalent named genera the nodes. `reference =
"top_prevalence"` uses the most prevalent named genus (Aeromonas) as the
denominator instead and the next `P` named genera as nodes; every real-data fit
in this project before 2026-09-28 used that setting.

## Guarantees

`prepare_zebrafish.R` ends in an assertion gate over **20 landmarks** — per-day
counts per group (infected 7,10,11,13,14,13,13; not infected 23,20,19,17,15,17,15),
genera above each prevalence threshold (43 / 25 / 16 at 5% / 10% / 20%), depth
min–median–max (7,400 / 22,452 / 53,971), the reads in `NONE`, and the P = 15
reference and node set under both denominators. If the raw data or the recipe drifts it **fails** rather
than writing a different dataset underneath the analyses. Genus aggregation is
also reconciled row-by-row against `asv.tab` for all 237 samples.

`zebrafish_real_depths.rds` is regenerated by the same script and **checked
against the shipped copy**; if it ever differed the script stops rather than
overwriting, because that would silently change every `depth_mode = "real"`
simulation. It currently regenerates identically (MD5
`39eec824d92d27fbbed299194bf57515`).

Three ways these files are silently misread — all closed and asserted against in
the script, and worth knowing before writing any other reader:
`read.table()`'s default `check.names = TRUE` mangles the ASV ids so that **0 of
841** match `tax.tab` and every genus count becomes zero; `metadata.tab` is
comma-separated despite its `.tab` name (`asv.tab` and `tax.tab` genuinely are
tab-separated); and its header has 18 fields against 19 in every data row, so a
reader that does not treat the first field as a row index shifts every column by
one and reads `DaysPE` as a sample-ID string.
