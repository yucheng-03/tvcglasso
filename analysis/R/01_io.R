# ---------------------------------------------------------------------------
# analysis/R/01_io.R  --  locating, reading and validating fitted results.
#
# One glob, one reader, one completeness gate. Every figure/table script goes
# through here so that a partially-synced results tree cannot silently produce
# a thin curve averaged over a different number of seeds than its neighbours.
# ---------------------------------------------------------------------------

#' Enumerate the result files for one method.
#'
#' Globs *.rds only. A failed fit writes a sibling `<file>.rds.FAILED` text
#' marker instead of a result (simulation/run_comparison.R), and stale markers
#' can sit next to a later successful .rds, so the marker count is reported but
#' never used to exclude a unit that has a real file.
list_fits <- function(method, cells = CELLS, seeds = SEEDS,
                      roots = RESULT_ROOTS) {
  dir <- roots[[method]]
  if (is.null(dir) || !dir.exists(dir)) {
    stop("results directory for '", method, "' not found: ",
         if (is.null(dir)) "<unset>" else dir,
         "\n  Set it in analysis/config.R or via the environment override.",
         call. = FALSE)
  }
  files <- list.files(dir, pattern = "^cell[0-9]{2}_seed[0-9]{3}\\.rds$",
                      full.names = TRUE)
  if (!length(files)) {
    stop("no cellCC_seedSSS.rds files in ", dir, call. = FALSE)
  }
  base <- basename(files)
  out <- data.frame(
    method = method,
    cell   = as.integer(substr(base, 5, 6)),
    seed   = as.integer(substr(base, 12, 14)),
    path   = files,
    stringsAsFactors = FALSE
  )
  out <- out[out$cell %in% cells & out$seed %in% seeds, , drop = FALSE]
  out[order(out$cell, out$seed), , drop = FALSE]
}

#' Read one result file.
read_fit <- function(path) readRDS(path)

#' Per-cell design metadata, read from the files themselves.
#'
#' The cell definition travels inside every result (`x$cell`), so the figure
#' layout is derived from the data rather than from a hardcoded label vector
#' that can silently go stale.
cell_meta <- function(method = "tvcglasso", cells = CELLS, roots = RESULT_ROOTS) {
  fits <- list_fits(method, cells = cells, roots = roots)
  first <- fits[!duplicated(fits$cell), , drop = FALSE]
  rows <- lapply(seq_len(nrow(first)), function(i) {
    x <- read_fit(first$path[i])
    data.frame(cell  = first$cell[i],
               P     = x$cell$P,
               n     = x$cell$n,
               m     = x$cell$m,
               depth = x$cell$depth_mode,
               stringsAsFactors = FALSE)
  })
  meta <- do.call(rbind, rows)
  meta[order(meta$cell), , drop = FALSE]
}

#' Assert that a method x cell has the expected number of seeds.
#'
#' STRICT_COMPLETENESS = TRUE turns a short cell into an error; FALSE reports
#' it and lets the caller record the realised n per curve. Either way the
#' shortfall is never silent.
require_complete <- function(fits, expected = length(SEEDS),
                             strict = STRICT_COMPLETENESS) {
  tab <- table(fits$cell)
  short <- tab[tab < expected]
  if (length(short)) {
    msg <- paste0(
      "incomplete cells for method '", fits$method[1], "': ",
      paste(sprintf("cell%02d=%d/%d", as.integer(names(short)),
                    as.integer(short), expected), collapse = ", ")
    )
    if (isTRUE(strict)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(tab)
}

#' Count stale/live .FAILED markers (reported in the provenance table).
count_failed_markers <- function(method, roots = RESULT_ROOTS) {
  dir <- roots[[method]]
  if (is.null(dir) || !dir.exists(dir)) return(c(total = NA, stale = NA, live = NA))
  mk <- list.files(dir, pattern = "\\.rds\\.FAILED$", full.names = TRUE)
  if (!length(mk)) return(c(total = 0L, stale = 0L, live = 0L))
  rds <- sub("\\.FAILED$", "", mk)
  stale <- sum(file.exists(rds))
  c(total = length(mk), stale = stale, live = length(mk) - stale)
}

#' Memoise an expensive extraction to analysis/cache/.
#'
#' The cache is a convenience only: deleting analysis/cache/ and re-running
#' must reproduce every figure bit-for-bit.
with_cache <- function(name, expr, dir = CACHE_DIR, refresh = FALSE) {
  f <- file.path(dir, paste0(name, ".rds"))
  if (!refresh && file.exists(f)) return(readRDS(f))
  val <- force(expr)
  saveRDS(val, f)
  val
}

#' Number of true positives / negatives available in a fit's truth.
#'
#' Micro-averaged over eval_slices and the strict upper triangle -- exactly the
#' convention edge_fpr_tpr() uses (R/roc_utils.R), so counts derived here are
#' consistent with the stored FPR/TPR.
truth_counts <- function(x) {
  sl <- x$eval_slices
  npos <- 0; nneg <- 0
  for (k in sl) {
    T <- x$true_Omega[[k]]
    ut <- upper.tri(T)
    npos <- npos + sum(T[ut] != 0)
    nneg <- nneg + sum(T[ut] == 0)
  }
  c(npos = npos, nneg = nneg)
}

#' Precision and F1 implied EXACTLY by (FPR, TPR) and the truth's class sizes.
#'
#' Several methods store FPR/TPR but not precision (CGLasso stores no precision
#' at all; tvcglasso's per-lambda `detail` carries neither precision nor F1).
#' Because the confusion table is fully determined by (FPR, TPR, npos, nneg),
#' these are exact identities, NOT the F1-inversion fallback the prototype
#' scripts used.
pr_from_rates <- function(FPR, TPR, npos, nneg) {
  TP <- TPR * npos
  FP <- FPR * nneg
  prec <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  f1   <- ifelse((2 * TP + FP + (npos - TP)) > 0,
                 2 * TP / (2 * TP + FP + (npos - TP)), NA_real_)
  list(precision = prec, F1 = f1)
}
