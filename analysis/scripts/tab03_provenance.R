# ---------------------------------------------------------------------------
# analysis/scripts/tab03_provenance.R
#
# PRODUCES  analysis/tables/tab03_provenance.csv
#
# The reproducibility receipt for the code-availability section: which commit,
# which R, how many (cell, seed) units each method actually contributed, and --
# the load-bearing one -- evidence that the four methods were run on
# BIT-IDENTICAL DATA.
#
# Pairing is checked, not assumed. Every result file stores
# provenance$data_fingerprint (the sum of all counts in that unit's data set);
# generate_data() is deterministic in (cell, seed), so equal fingerprints for a
# given (cell, seed) across method folders means the methods saw the same data
# and per-seed paired inference is legitimate rather than merely marginal.
# ---------------------------------------------------------------------------

source(here::here("analysis", "config.R"))
for (f in list.files(here::here("analysis", "R"), full.names = TRUE, pattern = "\\.R$")) source(f)

REFRESH <- as.logical(Sys.getenv("TVCG_REFRESH", "FALSE"))

gather_prov <- function() {
  rows <- list()
  for (mth in METHODS) {
    cells_m <- if (mth == "JGL") intersect(CELLS, JGL_CELLS) else CELLS
    fits <- tryCatch(list_fits(mth, cells = cells_m), error = function(e) NULL)
    if (is.null(fits)) { message("skipping ", mth, " (no results)"); next }
    for (i in seq_len(nrow(fits))) {
      x <- read_fit(fits$path[i])
      p <- x$provenance
      rows[[length(rows) + 1L]] <- data.frame(
        method      = mth,
        cell        = fits$cell[i],
        seed        = fits$seed[i],
        git         = p$git %||% NA_character_,
        R           = p$R %||% NA_character_,
        fingerprint = p$data_fingerprint %||% NA_real_,
        secs        = x$result$secs %||% NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

prov <- with_cache("tab03_provenance_raw", gather_prov(), refresh = REFRESH)

## --- per-method summary -----------------------------------------------------
summ <- do.call(rbind, lapply(split(prov, prov$method), function(g) {
  fm <- count_failed_markers(g$method[1])
  data.frame(
    method       = g$method[1],
    n_units      = nrow(g),
    n_cells      = length(unique(g$cell)),
    seeds_min    = min(g$seed), seeds_max = max(g$seed),
    git          = paste(unique(g$git), collapse = ";"),
    R            = paste(unique(g$R), collapse = ";"),
    core_hours   = round(sum(g$secs, na.rm = TRUE) / 3600, 1),
    failed_total = fm[["total"]], failed_stale = fm[["stale"]],
    failed_live  = fm[["live"]],
    stringsAsFactors = FALSE
  )
}))
summ <- summ[match(METHODS, summ$method), , drop = FALSE]
summ <- summ[!is.na(summ$method), , drop = FALSE]

## --- pairing check ----------------------------------------------------------
## For each (cell, seed) present in >1 method, do all methods report the same
## fingerprint?
key <- paste(prov$cell, prov$seed, sep = "_")
sp  <- split(prov$fingerprint, key)
sp  <- sp[vapply(sp, length, integer(1)) > 1L]
agree <- vapply(sp, function(v) length(unique(v)) == 1L, logical(1))
pairing <- data.frame(
  units_compared = length(sp),
  units_agreeing = sum(agree),
  units_disagree = sum(!agree),
  stringsAsFactors = FALSE
)

write.csv(summ, file.path(TAB_DIR, "tab03_provenance.csv"), row.names = FALSE)
write.csv(pairing, file.path(TAB_DIR, "tab03_pairing_check.csv"), row.names = FALSE)
message("wrote ", file.path(TAB_DIR, "tab03_provenance.csv"))
print(summ)
print(pairing)
if (pairing$units_disagree > 0) {
  warning("PAIRING VIOLATION: ", pairing$units_disagree,
          " (cell,seed) units have different data fingerprints across methods. ",
          "Paired per-seed comparison is NOT valid for those units.", call. = FALSE)
}
