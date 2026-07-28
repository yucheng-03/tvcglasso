# ---------------------------------------------------------------------------
# analysis/R/03_deployed.R  --  ONE normalised operating-point row from four
# structurally different `deployed` objects.
#
# WHY THIS FILE EXISTS. The four methods do NOT share a row shape, so the
# as-deployed table cannot be built by rbind():
#
#   tvcglasso  flat, 10 fields   selector lambda_index lambda FPR TPR precision
#                                F1 n_edges df refit_converged
#   JGL        flat, 11 fields   ... lambda1 lambda2 ... refit_exists
#   tvmgm      flat,  7 fields   selector bw FPR TPR precision F1 n_edges
#                                (no lambda at all -- the refit re-selects PER
#                                NODE, so there is no single grid index)
#   CGLasso    NESTED            selector n_valid_* n_rho_total pre post
#                                incomplete; the numbers live one level down in
#                                $post, which additionally has NO precision and
#                                calls the edge count `edges`, not `n_edges`
#
# WHAT "pre" MEANS IS ALSO METHOD-SPECIFIC (it is each method's own penalised
# native selector, before the relaxed refit):
#   tvcglasso  argmin of detail$BIC_pre (the Z-aligned pre criterion)
#   CGLasso    deployed$pre              (per-slice BIC on the pass-1 path)
#   JGL        refit$by_selector$AIC$pre (a ready-made detail row)
#   tvmgm      result$native             (mgm's published penalised per-node EBIC)
#
# PRECISION IS DERIVED, NOT INVERTED. Where a method stores no precision we
# recover it exactly from (FPR, TPR) and the truth's class sizes -- the
# confusion table is fully determined by those four numbers. The prototype
# scripts inverted F1 instead, which is algebraically equivalent only when F1
# and recall are both finite and 2*rec - F1 > 0.
#
# HARD RULES ENCODED HERE:
#  * tvcglasso: `detail$converged` is the BASE fit's flag while `detail$BIC` in
#    the same row is refit-POST. Re-selecting with
#    which.min(detail$BIC[detail$converged]) silently filters the refit
#    criterion by base convergence and does NOT reproduce the shipped
#    deployed$lambda_index. We read estimate$refit[[i]]$converged instead.
#  * The deployed point lies on the plotted ROC for tvcglasso only. CGLasso
#    deploys a per-SLICE amalgam, tvmgm a per-NODE amalgam, and JGL's deployed
#    lambda2 can differ from the lambda2 whose curve is drawn. `on_curve` records
#    this so a figure never implies otherwise.
# ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

.num <- function(v) if (is.null(v) || !length(v)) NA_real_ else as.numeric(v)[1]

.row <- function(method, x, stage, selector, FPR, TPR, precision, F1,
                 n_edges, df, index, on_curve, status = NA_character_) {
  data.frame(
    method    = method,
    cell      = x$cell_idx,
    seed      = x$seed,
    P         = x$cell$P,
    n         = x$cell$n,
    m         = x$cell$m,
    depth     = x$cell$depth_mode,
    stage     = stage,
    selector  = selector %||% NA_character_,
    FPR       = .num(FPR),
    TPR       = .num(TPR),
    precision = .num(precision),
    F1        = .num(F1),
    n_edges   = .num(n_edges),
    df        = .num(df),
    index     = .num(index),
    on_curve  = on_curve,
    status    = status,
    stringsAsFactors = FALSE
  )
}

#' The eligible candidate set for tvcglasso, matching method_tv.R exactly.
.tvcg_eligible <- function(r) {
  rf <- r$estimate$refit
  if (is.null(rf)) return(rep(FALSE, nrow(r$detail)))
  vapply(seq_along(rf), function(i) {
    e <- rf[[i]]
    isTRUE(e$converged) && !is.null(e$post)
  }, logical(1))
}

#' Normalised operating point for one fit.
#'
#' @param stage "refit" (the deployed point, each method's native selector
#'   applied after the relaxed refit), "pre" (the same selector before the
#'   refit), or "oracle" (the truth-maximising F1 point -- TRUTH-AWARE, valid
#'   only as a labelled simulation ceiling, never as a deployable result).
operating_point <- function(x, method = x$method, stage = c("refit", "pre", "oracle")) {
  stage <- match.arg(stage)
  r  <- x$result
  tc <- truth_counts(x)

  derive <- function(FPR, TPR) pr_from_rates(FPR, TPR, tc[["npos"]], tc[["nneg"]])

  ## ---- tvcglasso ----------------------------------------------------------
  if (method == "tvcglasso") {
    d <- r$detail
    if (stage == "refit") {
      dp <- r$deployed
      if (is.null(dp)) return(NULL)
      return(.row(method, x, "refit", dp$selector, dp$FPR, dp$TPR,
                  dp$precision, dp$F1, dp$n_edges, dp$df, dp$lambda_index,
                  on_curve = TRUE,
                  status = if (isTRUE(dp$refit_converged)) "converged" else "not_converged"))
    }
    if (stage == "pre") {
      ok <- .tvcg_eligible(r) & is.finite(d$BIC_pre)
      if (!any(ok)) return(NULL)
      i <- which(ok)[which.min(d$BIC_pre[ok])]
      pr <- derive(d$FPR[i], d$TPR[i])
      return(.row(method, x, "pre", "BIC_pre", d$FPR[i], d$TPR[i],
                  pr$precision, pr$F1, d$n_edges[i], d$df[i], i,
                  on_curve = TRUE))
    }
    pr <- derive(d$FPR, d$TPR)
    i <- which.max(pr$F1)
    if (!length(i)) return(NULL)
    return(.row(method, x, "oracle", "max-F1 (truth-aware)", d$FPR[i], d$TPR[i],
                pr$precision[i], pr$F1[i], d$n_edges[i], d$df[i], i,
                on_curve = TRUE))
  }

  ## ---- CGLasso ------------------------------------------------------------
  if (method == "CGLasso") {
    dp <- r$deployed
    if (stage %in% c("refit", "pre")) {
      if (is.null(dp)) return(NULL)
      s <- if (stage == "refit") dp$post else dp$pre
      if (is.null(s)) return(NULL)
      pr <- derive(s$FPR, s$TPR)
      return(.row(method, x, stage, dp$selector, s$FPR, s$TPR,
                  pr$precision, s$F1, s$edges, NA, NA,
                  on_curve = FALSE,          # per-slice amalgam
                  status = if (isTRUE(all(s$sel_converged))) "converged" else "mixed"))
    }
    d <- r$detail
    pr <- derive(d$FPR, d$TPR)
    i <- which.max(pr$F1)
    if (!length(i)) return(NULL)
    return(.row(method, x, "oracle", "max-F1 (truth-aware)", d$FPR[i], d$TPR[i],
                pr$precision[i], pr$F1[i], d$n_edges[i], d$df[i], d$index[i],
                on_curve = TRUE))
  }

  ## ---- JGL ----------------------------------------------------------------
  if (method == "JGL") {
    if (stage == "refit") {
      dp <- r$deployed
      if (is.null(dp)) return(NULL)
      return(.row(method, x, "refit", dp$selector, dp$FPR, dp$TPR,
                  dp$precision, dp$F1, dp$n_edges, dp$df, dp$lambda_index,
                  on_curve = FALSE,          # deployed lambda2 may differ from
                                             # the lambda2 whose curve is drawn
                  status = if (isTRUE(dp$refit_exists)) "refit_exists" else "no_refit"))
    }
    if (stage == "pre") {
      p <- r$refit$by_selector[["AIC"]]$pre
      if (is.null(p)) return(NULL)
      return(.row(method, x, "pre", "AIC_pre", p$FPR, p$TPR, p$precision,
                  p$F1, p$n_edges, p$df, NA, on_curve = FALSE))
    }
    o <- r$refit$oracle
    if (is.null(o)) return(NULL)
    return(.row(method, x, "oracle", "max-F1 (truth-aware)", o$FPR, o$TPR,
                o$precision, o$F1, o$n_edges, o$df, r$refit$oracle_index,
                on_curve = FALSE))
  }

  ## ---- tvmgm --------------------------------------------------------------
  if (method == "tvmgm") {
    if (stage == "refit") {
      dp <- r$deployed
      if (is.null(dp)) return(NULL)
      return(.row(method, x, "refit", dp$selector, dp$FPR, dp$TPR,
                  dp$precision, dp$F1, dp$n_edges, NA, NA,
                  on_curve = FALSE))         # per-node amalgam
    }
    if (stage == "pre") {
      nv <- r$native
      if (is.null(nv)) return(NULL)
      return(.row(method, x, "pre", nv$selector, nv$FPR, nv$TPR, nv$precision,
                  nv$F1, nv$edges, NA, NA, on_curve = FALSE))
    }
    o <- r$estimate$refit$oracle
    if (is.null(o)) return(NULL)
    return(.row(method, x, "oracle", "max-F1 (truth-aware)", o$FPR, o$TPR,
                o$precision, o$F1, o$edges, NA, NA, on_curve = TRUE))
  }

  stop("unknown method: ", method, call. = FALSE)
}

#' Collect operating points across a method's fits.
collect_points <- function(method, stages = c("refit", "pre"),
                           cells = CELLS, seeds = SEEDS) {
  fits <- list_fits(method, cells = cells, seeds = seeds)
  require_complete(fits)
  rows <- list()
  for (i in seq_len(nrow(fits))) {
    x <- read_fit(fits$path[i])
    for (st in stages) {
      p <- tryCatch(operating_point(x, method, st), error = function(e) NULL)
      if (!is.null(p)) rows[[length(rows) + 1L]] <- p
    }
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

#' Aggregate operating points over seeds: mean, standard deviation, standard error.
#'
#' BOTH spread measures are stored because they answer different questions and
#' the figure must say which one it draws:
#'   sd  -- how much a SINGLE replicate varies (what a practitioner would see)
#'   se  -- how precisely the MEAN is estimated (sd / sqrt(n); at n = 100 seeds
#'          this is ten times smaller and says almost nothing about spread)
#' The published JASA bar figure in Tian et al. draws mean +/- 1 sd with the
#' lower end floored at zero, and analysis/scripts/fig02 follows it.
#'
#' n_seed is reported so a cell whose campaign leg is short is visible in the
#' table rather than hidden inside the mean.
summarise_points <- function(pts, metrics = c("FPR", "TPR", "precision", "F1",
                                              "n_edges", "df")) {
  key <- interaction(pts$method, pts$stage, pts$cell, drop = TRUE, sep = "\r")
  out <- lapply(split(pts, key), function(g) {
    base <- g[1, c("method", "stage", "cell", "P", "n", "m", "depth", "selector")]
    for (mt in metrics) {
      v <- g[[mt]]
      v <- v[is.finite(v)]
      base[[paste0(mt, "_mean")]] <- if (length(v)) mean(v) else NA_real_
      base[[paste0(mt, "_sd")]]   <- if (length(v) > 1) stats::sd(v) else NA_real_
      base[[paste0(mt, "_se")]]   <- if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA_real_
    }
    base$n_seed <- nrow(g)
    base
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res[order(res$cell, match(res$method, METHODS), res$stage), , drop = FALSE]
}
