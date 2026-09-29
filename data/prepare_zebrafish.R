# ============================================================================
# data/prepare_zebrafish.R — raw Zebrafish 16S tables -> cleaned data.
#
# Run from the repo root:   Rscript data/prepare_zebrafish.R
#
# Reads ONLY the three unmodified source tables in data/raw/ and
# writes every derived object into data/. Nothing here depends on a fitted
# model, so the cleaned data can be regenerated from scratch at any time.
# Base R only (the published preprocessing uses tidyverse; this repo does not).
#
# Outputs
#   data/zebrafish_clean.rds        the cleaned genus-level dataset (see below)
#   data/zebrafish_real_depths.rds  per-sample total reads — the pool the
#                                   simulation samples for depth_mode = "real"
#
# THREE traps that silently corrupt this dataset, all closed below.
#   (1) read.table()'s default check.names = TRUE rewrites "-" to "." in the ASV
#       ids, after which NONE of the 841 ids matches a row name of tax.tab
#       (measured: 0 of 841) and every genus count aggregates to zero — with no
#       error. Read with check.names = FALSE and assert the ids match.
#   (2) metadata.tab is COMMA-delimited despite the .tab extension; read with
#       read.csv(), not read.table(). Verified: 18 comma fields, 1 tab field.
#   (3) `tax$Genus == g` is NA-valued when g is NA, so `which()` of it returns
#       integer(0) and the unassigned-genus column silently becomes ALL ZEROS —
#       discarding 134,761 of the 4,623,919 reads in the retained samples
#       (2.91%). Branch on is.na() explicitly, as the published preprocessing
#       does. See zeb_aggregate_genus().
# ============================================================================

suppressPackageStartupMessages(library(here))

RAW_DIR  <- here::here("data", "raw")
OUT_DIR  <- here::here("data")
CLEAN_RDS <- file.path(OUT_DIR, "zebrafish_clean.rds")
DEPTH_RDS <- file.path(OUT_DIR, "zebrafish_real_depths.rds")

## Landmark values this script must reproduce. They are asserted, not printed
## for eyeballing: if the raw files or the recipe ever drift, this script FAILS
## rather than quietly writing a different dataset underneath the analyses.
EXPECT <- list(
  n_samples_raw   = 237L,
  n_asv           = 841L,
  n_genus         = 260L,
  n_samples_kept  = 207L,
  days            = c(7L, 10L, 21L, 30L, 43L, 59L, 86L),
  n_infected      = 81L,
  n_not_infected  = 126L,
  per_day_infected     = c(7L, 10L, 11L, 13L, 14L, 13L, 13L),
  per_day_not_infected = c(23L, 20L, 19L, 17L, 15L, 17L, 15L),
  n_prev_gt_05    = 43L,   # prevalence > 5%  (the published filter)
  n_prev_gt_10    = 25L,   # prevalence > 10% -> P = 24 nodes + 1 reference
  n_prev_gt_20    = 16L,   # prevalence > 20% -> P = 15 nodes + 1 reference
  depth_min       = 7400,
  depth_median    = 22452,
  depth_max       = 53971,
  # Guards TRAP (3): if the unassigned-genus column is ever silently zeroed
  # again, this is the check that fails.
  none_total_reads = 134761,
  # The default ALR denominator and the P = 15 node set under it: NONE, as in the
  # published preprocessing, with the 15 most prevalent named genera as nodes.
  reference        = "NONE",
  nodes_P15 = c("Aeromonas", "Cetobacterium", "Pseudomonas", "Plesiomonas",
                "ZOR0006", "Acinetobacter", "Shewanella", "Paucibacter",
                "Chitinibacter", "Crenobacter", "Flavobacterium", "Mycoplasma",
                "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                "Cloacibacterium", "Fluviicola"),
  # The alternative denominator (reference = "top_prevalence") and its P = 15
  # node set. Every real-data fit before 2026-09-28 used this setting, so it is
  # asserted too: those results must stay reproducible.
  reference_top_prevalence = "Aeromonas",
  nodes_P15_top_prevalence = c("Cetobacterium", "Pseudomonas", "Plesiomonas", "ZOR0006",
                "Acinetobacter", "Shewanella", "Paucibacter", "Chitinibacter",
                "Crenobacter", "Flavobacterium", "Mycoplasma",
                "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",
                "Cloacibacterium", "Fluviicola", "Phreatobacter")
)

# ---------------------------------------------------------------------------
# 1. Read the three raw tables exactly as shipped.
# ---------------------------------------------------------------------------

#' Read the raw Zebrafish tables.
#'
#' @return list(asv = samples x ASV counts, tax = ASV x taxonomy,
#'              meta = samples x metadata)
zeb_read_raw <- function(dir = RAW_DIR) {
  need <- c("asv.tab", "tax.tab", "metadata.tab")
  miss <- need[!file.exists(file.path(dir, need))]
  if (length(miss)) {
    stop("missing raw Zebrafish file(s) in ", dir, ": ", paste(miss, collapse = ", "),
         "\n  See data/README.md for where these come from.", call. = FALSE)
  }

  # asv.tab / tax.tab: whitespace-delimited, first field of each data row is the
  # row name (the header has one fewer field), so read.table assigns row names.
  # TRAP (1): check.names = FALSE keeps the "-" in the ASV ids intact.
  asv <- read.table(file.path(dir, "asv.tab"), header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)
  tax <- read.table(file.path(dir, "tax.tab"), header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)

  # The assertion that makes TRAP (1) impossible to reintroduce silently.
  if (!setequal(colnames(asv), rownames(tax))) {
    stop("ASV ids in asv.tab do not match the row names of tax.tab (",
         length(intersect(colnames(asv), rownames(tax))), " of ", ncol(asv),
         " match).\n  If this is 0, the ids were mangled on read ",
         "(check.names must be FALSE).", call. = FALSE)
  }

  # TRAP (2): comma-delimited despite the .tab extension.
  meta <- read.csv(file.path(dir, "metadata.tab"), stringsAsFactors = FALSE)

  list(asv = as.matrix(asv), tax = tax, meta = meta)
}

# ---------------------------------------------------------------------------
# 2. Aggregate ASVs to genus.
# ---------------------------------------------------------------------------

#' Sum ASV counts within each genus.
#'
#' The 298 ASVs whose Genus is NA (unassigned at genus level) are pooled into a
#' single column named "NONE", the same convention as the published
#' preprocessing. Those reads are real and are kept, so the per-sample total is
#' the full library and the composition stays closed.
#'
#' TRAP (3): writing this as `rownames(tax)[which(tax$Genus == g)]` looks
#' equivalent but is not — for g = NA the comparison is NA, `which()` drops it,
#' and "NONE" comes out ALL ZEROS, silently discarding 2.91% of reads. The
#' is.na() branch below is what makes the column real; the rowSums assertion
#' afterwards is what makes the mistake impossible to reintroduce unnoticed.
zeb_aggregate_genus <- function(asv, tax) {
  genera <- unique(as.character(tax$Genus))
  gmat <- vapply(genera, function(g) {
    ids <- if (is.na(g)) rownames(tax)[is.na(tax$Genus)]
           else          rownames(tax)[!is.na(tax$Genus) & tax$Genus == g]
    ids <- intersect(ids, colnames(asv))
    if (!length(ids)) rep(0, nrow(asv)) else rowSums(asv[, ids, drop = FALSE])
  }, numeric(nrow(asv)))

  colnames(gmat) <- ifelse(is.na(genera), "NONE", genera)
  rownames(gmat) <- rownames(asv)

  # If the gsub() trap above were re-broken, every count would be 0 and this
  # is the assertion that catches it.
  stopifnot(all(rowSums(gmat) == rowSums(asv)))
  storage.mode(gmat) <- "integer"
  gmat
}

#' One representative higher-rank lineage per genus.
#'
#' Kingdom..Family are constant within a genus in this table; where they are
#' not, the first non-NA value is kept and the genus is flagged.
zeb_genus_lineage <- function(tax) {
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family")
  key <- ifelse(is.na(tax$Genus), "NONE", as.character(tax$Genus))
  out <- lapply(ranks, function(r) {
    v <- as.character(tax[[r]])
    vapply(split(v, key), function(x) {
      x <- x[!is.na(x)]
      if (!length(x)) NA_character_ else x[1]
    }, character(1))
  })
  names(out) <- ranks
  amb <- vapply(split(seq_len(nrow(tax)), key), function(i) {
    any(vapply(ranks, function(r) length(unique(na.omit(tax[[r]][i]))) > 1, logical(1)))
  }, logical(1))
  data.frame(genus = names(amb), out, lineage_ambiguous = unname(amb),
             row.names = NULL, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# 3. Clean: drop pre-exposure fish, split by parasite burden, describe taxa.
# ---------------------------------------------------------------------------

#' Build the cleaned Zebrafish dataset.
#'
#' Steps, in order:
#'   1. aggregate ASVs to genus (237 x 260)
#'   2. DROP the 30 fish sampled at DaysPE == 0. They are pre-exposure animals,
#'      not an independent draw from the post-exposure time course, and the
#'      published preprocessing drops them too. 207 fish over 7 days remain.
#'   3. group by parasite burden: infected = Total > 0 (n = 81),
#'      not_infected = Total == 0 or NA (n = 126). NA means no worm count was
#'      recorded; those fish are unexposed controls.
#'   4. sequencing depth = total reads per sample, summed over ALL genera
#'      (i.e. the full library, not the analysed subset).
#'   5. prevalence and mean relative abundance are computed over the 207 KEPT
#'      samples, so they describe the analysed cohort.
zeb_clean <- function(raw) {
  gmat_all <- zeb_aggregate_genus(raw$asv, raw$tax)
  meta <- raw$meta
  stopifnot(identical(rownames(gmat_all), meta$Seq_ID))

  keep <- meta$DaysPE != 0
  counts <- gmat_all[keep, , drop = FALSE]
  md <- meta[keep, , drop = FALSE]

  worms <- md$Total
  infected <- !is.na(worms) & worms > 0
  depth <- unname(rowSums(counts))

  samples <- data.frame(
    sample_id  = rownames(counts),
    day        = as.integer(md$DaysPE),
    group      = ifelse(infected, "infected", "not_infected"),
    depth      = depth,
    exposure   = md$Exposure,
    worm_total = worms,
    tank       = md$tank,
    row.names  = NULL, stringsAsFactors = FALSE
  )

  prevalence <- colMeans(counts > 0)
  rel_ab     <- colMeans(counts / rowSums(counts))
  lin        <- zeb_genus_lineage(raw$tax)

  taxa <- data.frame(genus = colnames(counts),
                     prevalence = unname(prevalence),
                     rel_abundance = unname(rel_ab),
                     # "NONE" pools ASVs unassigned at genus level: kept in the
                     # counts, never eligible as a network node (see zeb_slices).
                     is_unassigned = colnames(counts) == "NONE",
                     row.names = NULL, stringsAsFactors = FALSE)
  taxa <- merge(taxa, lin, by = "genus", all.x = TRUE, sort = FALSE)
  taxa <- taxa[match(colnames(counts), taxa$genus), , drop = FALSE]
  taxa$rank <- rank(-taxa$prevalence, ties.method = "min")
  rownames(taxa) <- NULL

  list(
    counts  = counts,
    samples = samples,
    taxa    = taxa,
    days    = sort(unique(samples$day)),
    dropped = list(reason = "DaysPE == 0 (pre-exposure fish)",
                   n = sum(!keep), sample_id = meta$Seq_ID[!keep])
  )
}

# ---------------------------------------------------------------------------
# 4. Taxon ranking and analysis-ready slices.
# ---------------------------------------------------------------------------

#' Rank genera for node selection, deterministically.
#'
#' Primary key is PREVALENCE (fraction of the 207 samples in which the genus is
#' observed) -- the convention used throughout this project, and the same
#' quantity the published preprocessing filters on. Prevalence has exact ties
#' (e.g. Pseudoduganella and Ignatzschineria are both 0.0918, and that tie falls
#' exactly at the P = 25 boundary), so ties are broken explicitly by mean
#' relative abundance and then alphabetically. Without an explicit rule the
#' winner is decided by column order, i.e. by the order genera happen to appear
#' in tax.tab -- reproducible by accident rather than by specification.
#'
#' The explicit rule selects the SAME genus set as plain prevalence ordering at
#' P = 10, 15 and 25; it only fixes the order within tied blocks.
zeb_rank_taxa <- function(clean) {
  t <- clean$taxa
  t$genus[order(-t$prevalence, -t$rel_abundance, t$genus)]
}

#' Cut the cleaned data into the per-day count matrices a model consumes.
#'
#' Returns, for one group, a list of m matrices (one per sampling day), each
#' n_k x (P+1) raw counts. The ALR REFERENCE taxon is the LAST column, which is
#' the layout every estimator in this repo expects (Z is n_k x P, Omega is PxP,
#' and the reference has no row/column in Omega).
#'
#' "NONE" IS NEVER A NODE. It is not a genus — it pools 298 ASVs that could not
#' be assigned at genus level, spanning unrelated lineages (which is why
#' taxa$lineage_ambiguous is TRUE for it) — so an estimated edge to it would
#' have no biological reading. It is kept in the cleaned counts (the reads are
#' real and the depth must be the true library size) but excluded from the node
#' pool; letting it compete on prevalence would rank it 3rd and displace a real
#' genus. By default it is the ALR denominator instead.
#'
#' @param P number of network nodes (the reference is additional).
#' @param reference which taxon becomes the ALR denominator.
#'   "NONE" (default): the pooled unassigned column, which is what the published
#'     preprocessing uses (Tian et al. 2023); the top P named genera are then all
#'     nodes.
#'   "top_prevalence": the most prevalent NAMED genus (Aeromonas, present in
#'     99.5% of samples); the next P named genera are the nodes. This was the
#'     default until 2026-09-28, so every real-data fit before that date used it.
zeb_slices <- function(clean, P = 15L, group = c("infected", "not_infected"),
                       reference = c("NONE", "top_prevalence")) {
  group     <- match.arg(group)
  reference <- match.arg(reference)
  named     <- setdiff(zeb_rank_taxa(clean), "NONE")   # NONE is never a node

  if (reference == "top_prevalence") {
    if (length(named) < P + 1L)
      stop("P too large: only ", length(named), " named genera", call. = FALSE)
    ref_name   <- named[1]
    node_names <- named[2:(P + 1L)]
  } else {
    if (length(named) < P)
      stop("P too large: only ", length(named), " named genera", call. = FALSE)
    ref_name   <- "NONE"
    node_names <- named[seq_len(P)]
  }

  sel <- clean$samples$group == group
  X <- clean$counts[, c(node_names, ref_name), drop = FALSE]
  days <- clean$days
  slices <- lapply(days, function(d) X[sel & clean$samples$day == d, , drop = FALSE])
  names(slices) <- paste0("day", days)

  list(
    X            = slices,
    days         = days,
    t            = (days - min(days)) / (max(days) - min(days)),  # real irregular spacing on [0,1]
    n_per_slice  = vapply(slices, nrow, integer(1)),
    nodes        = node_names,
    reference    = ref_name,
    P            = length(node_names),
    group        = group
  )
}

# ---------------------------------------------------------------------------
# 5. Build, verify, write.
# ---------------------------------------------------------------------------

main <- function() {
  cat("Reading raw tables from ", RAW_DIR, "\n", sep = "")
  raw <- zeb_read_raw()
  stopifnot(nrow(raw$asv) == EXPECT$n_samples_raw,
            ncol(raw$asv) == EXPECT$n_asv,
            nrow(raw$tax) == EXPECT$n_asv)

  clean <- zeb_clean(raw)

  ## --- hard verification gate -------------------------------------------
  s <- clean$samples
  inf <- s$group == "infected"
  pd_inf <- vapply(clean$days, function(d) sum(inf & s$day == d), integer(1))
  pd_not <- vapply(clean$days, function(d) sum(!inf & s$day == d), integer(1))
  prev <- clean$taxa$prevalence

  stopifnot(
    ncol(clean$counts)              == EXPECT$n_genus,
    nrow(clean$counts)              == EXPECT$n_samples_kept,
    identical(as.integer(clean$days), EXPECT$days),
    sum(inf)                        == EXPECT$n_infected,
    sum(!inf)                       == EXPECT$n_not_infected,
    identical(pd_inf, EXPECT$per_day_infected),
    identical(pd_not, EXPECT$per_day_not_infected),
    sum(prev > 0.05)                == EXPECT$n_prev_gt_05,
    sum(prev > 0.10)                == EXPECT$n_prev_gt_10,
    sum(prev > 0.20)                == EXPECT$n_prev_gt_20,
    min(s$depth)                    == EXPECT$depth_min,
    median(s$depth)                 == EXPECT$depth_median,
    max(s$depth)                    == EXPECT$depth_max,
    sum(clean$counts[, "NONE"])     == EXPECT$none_total_reads
  )
  p15 <- zeb_slices(clean, 15L, "infected")
  p15_tp <- zeb_slices(clean, 15L, "infected", reference = "top_prevalence")
  stopifnot(identical(p15$reference, EXPECT$reference),
            setequal(p15$nodes, EXPECT$nodes_P15),
            length(p15$nodes) == 15L,
            identical(p15_tp$reference, EXPECT$reference_top_prevalence),
            setequal(p15_tp$nodes, EXPECT$nodes_P15_top_prevalence),
            length(p15_tp$nodes) == 15L)
  cat("verification gate: PASSED (", length(EXPECT), " landmark checks incl. the ",
      "P=15 node set under both denominators and the unassigned-read total)\n", sep = "")

  ## --- provenance --------------------------------------------------------
  src <- file.path(RAW_DIR, c("asv.tab", "tax.tab", "metadata.tab"))
  clean$provenance <- list(
    script       = "data/prepare_zebrafish.R",
    generated_on = as.character(Sys.Date()),
    R            = R.version.string,
    source_files = basename(src),
    source_md5   = unname(tools::md5sum(src)),
    source_bytes = unname(file.size(src)),
    recipe       = paste(
      "ASVs summed to genus (unassigned -> 'NONE');",
      "DaysPE == 0 dropped;",
      "infected = Total > 0;",
      "depth = total reads over all genera;",
      "prevalence and relative abundance over the 207 retained samples.")
  )

  ## --- depth pool used by the simulation ---------------------------------
  ## Plain unnamed numeric, in the original sample order. generators.R does
  ## sample(rd, n, replace = TRUE) on it for depth_mode = "real".
  depths <- as.numeric(clean$samples$depth)
  if (file.exists(DEPTH_RDS)) {
    old <- as.numeric(readRDS(DEPTH_RDS))
    if (!identical(old, depths)) {
      stop("REGENERATED DEPTH POOL DIFFERS FROM THE SHIPPED ONE.\n",
           "  shipped: n=", length(old), " min=", min(old), " median=", median(old), "\n",
           "  new:     n=", length(depths), " min=", min(depths), " median=", median(depths), "\n",
           "  Overwriting would silently change every depth_mode='real' simulation.\n",
           "  Resolve before proceeding.", call. = FALSE)
    }
    cat("depth pool: regenerated value is IDENTICAL to the shipped file (n = ",
        length(depths), ")\n", sep = "")
  }
  saveRDS(depths, DEPTH_RDS)

  ## `counts` is the single source of truth. Per-P analysis slices are NOT
  ## stored: they are a pure re-cut of these same columns, so caching them would
  ## duplicate the data and create a second thing to keep in sync. Get them with
  ##     source("data/prepare_zebrafish.R"); zeb_slices(clean, P = 15, "infected")
  saveRDS(clean, CLEAN_RDS)

  ## --- report ------------------------------------------------------------
  cat("\n", strrep("-", 70), "\n", sep = "")
  cat("cleaned dataset -> ", CLEAN_RDS, "\n", sep = "")
  cat("  genus counts   : ", nrow(clean$counts), " samples x ", ncol(clean$counts), " genera\n", sep = "")
  cat("  days           : ", paste(clean$days, collapse = ", "), "\n", sep = "")
  cat("  infected       : n = ", sum(inf), "  per day ", paste(pd_inf, collapse = ","), "\n", sep = "")
  cat("  not_infected   : n = ", sum(!inf), "  per day ", paste(pd_not, collapse = ","), "\n", sep = "")
  cat("  depth (reads)  : min ", min(s$depth), "  median ", median(s$depth),
      "  max ", max(s$depth), "\n", sep = "")
  cat("  genera by prevalence: >5% ", sum(prev > 0.05), "   >10% ", sum(prev > 0.10),
      "   >20% ", sum(prev > 0.20), "\n", sep = "")
  cat("depth pool      -> ", DEPTH_RDS, " (n = ", length(depths), ")\n", sep = "")
  cat("\nexample slice set  P = 15, infected  (zeb_slices, not stored):\n")
  cat("  reference taxon: ", p15$reference, "   (last column of every matrix)\n", sep = "")
  cat("  nodes          : ", paste(p15$nodes, collapse = ", "), "\n", sep = "")
  cat("  n per day      : ", paste(p15$n_per_slice, collapse = ", "), "\n", sep = "")
  cat(strrep("-", 70), "\n", sep = "")
  invisible(clean)
}

if (sys.nframe() == 0L) main()
