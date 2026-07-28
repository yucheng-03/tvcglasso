# ============================================================================
# methods/method_tvmgm.R — tvmgm (Haslbeck & Waldorp 2020) BASELINE: base fit + relaxed refit.
#
# The baseline is called FAITHFULLY -- we reproduce mgm as published and never modify its
# internals. Both passes stay inside mgm's own nodewise framework:
#   pass 1  mgm::tvmgm, kernel-weighted nodewise lasso, swept over a lambda grid -> the ROC;
#   pass 2  OUR relaxed refit (R/baselines/refit_mgm.R): each node refitted UNPENALISED on
#           ITS OWN directed lasso mask, re-selected per node by mgm's OWN EBIC.
#
# INPUT = ALR (PI decision 2026-06-24): per-slice CENTERED ALR (dat$Z_0, last-taxon
#   reference), unified with tvcglasso/CGLasso/JGL so every method targets the SAME
#   ALR-precision the truth is defined in (this removed the old CLR-space metric mismatch).
#   P-dimensional -- no P+1 CLR subsetting.
#
# ROC = LAMBDA-SWEEP (method-symmetric with the other three): tvmgm is fit at each grid
#   lambda (forced through a single-value lambdaSeq) and an edge is wadj != 0. This replaces
#   a score-sweep of one EBIC fit's weights, which is degenerate -- at the EBIC fit ~87% of
#   edges are lasso-killed, so thresholding the survivors gives a straight line at high FPR
#   (AUC 0.573 vs the lambda-sweep's 0.931 on P=15/n=20/real seed 1; diagnosed 2026-06-20).
#
# TIME = mgm's kernel: at each estimation point every observation is weighted by
#   exp(-0.5*((t_obs - t*)/bw)^2) and ONE weighted nodewise lasso is fit (tvmgm does not
#   average per-timepoint graphs). The refit reuses these exact weights.
#
# BANDWIDTH: mgm-native bwSelect (its own out-of-sample CV) when cfg$mgm_tune. NB bwSelect
#   was designed for a single dense time series (one row per timepoint) and indexes its test
#   models by observation ORDER, which is approximate for our replicated cross-section (n
#   samples share each slice time). Disclosed, not patched: it is mgm's own tuning routine and
#   selecting the bandwidth by anything else (e.g. the truth-maximising AUC) would be oracle
#   tuning of a competing method. The bandwidth MATTERS -- holding a cell fixed and varying
#   only bw over bwSeq, AUC moved 0.70-0.83 and the deployed F1 0.065-0.415 -- so the selected
#   value is recorded (bw, bw_sel, bwSeq, bw_source in $tuning); a bw_sel at an endpoint of
#   bwSeq means the CV wanted to leave the grid and should be read with that in mind.
#   `tuned` is TRUE only if bwSelect actually returned (a failure falls back to a fixed
#   bandwidth and is reported honestly via bw_source).
#
# SCALE: mgm standardises every Gaussian column internally (scale=TRUE). The refit is given
#   that SAME standardised design, so (i) its EBIC -- which is a unit-variance Gaussian
#   deviance, see refit_mgm.R -- is calibrated exactly as mgm's is, and (ii) the de-biased
#   coefficients live on the same scale as mgm's own wadj.
#
# LAMBDA GRID: DATA-ADAPTIVE by default (cfg$mgm_lambda_grid = NULL), using glmnet's OWN
#   lambda_max -- i.e. exactly the path mgm would generate if we did not force one, and the
#   same "anchor lambda_max on the data" convention tvcglasso and CGLasso follow. For each
#   estimation point and node, glmnet's lambda_max is the smallest penalty that zeroes every
#   coefficient of that weighted, standardised node regression; taking the MAX over all
#   (estpoint, node) guarantees the EMPTY graph, so the ROC starts at a REAL (0,0). The dense
#   end is lambda_max/1e4 (glmnet's own lambda.min.ratio for n > p), which reaches the FULL
#   graph, so the ROC also ENDS at a REAL (1,1): both corners come from real fits, never from
#   interpolation. (A hardcoded grid was reaching only FPR ~0.88-0.95, leaving the dense corner
#   uncovered, and its lambda_max = 1.2 was ~2x larger than the data needs -- wasted all-empty
#   points. Verified: data-adaptive lambda_max ~0.66 gives 0 edges at the sparse end and 100%
#   of edges at the dense end, on both a high- and a low-depth cell.)
#
# A lambda whose tvmgm fit ERRORS is recorded as status="fit_error" and EXCLUDED; it is never
# replaced by an all-zero graph (which would counterfeit a perfect sparse ROC corner).
#
# OUTPUT = the shared SHELL + a method-specific PAYLOAD (same contract as method_tv.R):
#   SHELL   roc / auc / scores / detail / deployed / native / tuning
#   PAYLOAD estimate$base  {wadj, wadjNodewise, lambda}
#           estimate$refit {landings, de-biased nets, FULL per-(estpoint,node,lambda) grid,
#                           diagnostics} -> every selector recomputable with NO re-fit.
# ============================================================================

# glmnet's OWN lambda_max for the kernel-weighted nodewise problem, maximised over every
# (estimation point, node) so that ALL node models are empty at lambda_max -> the empty graph.
# Zs = the standardised design mgm's node models see; w = that estpoint's kernel weights.
# For a weighted Gaussian lasso with standardisation, glmnet's lambda_max is
#   max_j | sum_i wn_i * xs_ij * (y_i - ybar_w) |,  wn = w/sum(w), xs = weighted-standardised X
# (verified against glmnet's own fit$lambda[1]: agreement to 4/4 exact, ratio 1.0000).
mgm_lambda_max <- function(Zs, tp, ep, bw) {
  P <- ncol(Zs)
  max(vapply(ep, function(e) {
    w <- exp(-0.5 * ((tp - e) / bw)^2); wn <- w / sum(w)
    max(vapply(seq_len(P), function(p) {
      y <- Zs[, p]; X <- Zs[, -p, drop = FALSE]
      mu <- colSums(wn * X); sdv <- sqrt(colSums(wn * sweep(X, 2, mu)^2))
      xs <- sweep(sweep(X, 2, mu), 2, pmax(sdv, 1e-12), "/")
      max(abs(colSums(wn * xs * (y - sum(wn * y)))))
    }, numeric(1)))
  }, numeric(1)))
}

run_method_tvmgm <- function(dat, cfg, eval_slices) {
  m  <- length(dat$Z_0); P <- ncol(dat$Z_0[[1]])
  tp <- rep(seq(0, 1, length.out = m), times = sapply(dat$Z_0, nrow))  # each observation's slice time
  ep <- seq(0, 1, length.out = m)                                      # estimation points = the slice times
  Za <- do.call(rbind, lapply(dat$Z_0, function(z) scale(z, center = TRUE, scale = FALSE)))
  gamma <- cfg$mgm_lambdaGam %||% 0.25                                 # mgm's EBIC gamma (package default)
  bwSeq <- cfg$mgm_bwSeq     %||% c(0.03, 0.05, 0.1, 0.2, 0.35, 0.6, 1.2)

  ## ---- bandwidth: mgm-native bwSelect (faithful) when requested ----
  bw <- cfg$mgm_bw %||% 0.25; bw_sel <- NA_real_; bw_source <- "fixed"
  if (isTRUE(cfg$mgm_tune)) {
    bwo <- tryCatch(bwSelect(data = Za, type = rep("g", P), level = rep(1, P), bwSeq = bwSeq,
                             bwFolds = cfg$mgm_bwFolds %||% 1L, bwFoldsize = cfg$mgm_bwFoldsize %||% 5L,
                             modeltype = "mgm", timepoints = tp, k = 2, lambdaSel = "EBIC",
                             lambdaGam = gamma, ruleReg = "OR", pbar = FALSE, signInfo = FALSE),
                    error = function(e) NULL)
    if (!is.null(bwo)) { bw_sel <- bwSeq[which.min(bwo$meanError)]; bw <- bw_sel; bw_source <- "bwSelect" }
    else                 bw_source <- "bwSelect_failed_fixed_fallback"
  }
  tuned <- identical(bw_source, "bwSelect")

  ## ---- lambda grid: DATA-ADAPTIVE (glmnet's own lambda_max at the SELECTED bandwidth, since
  ##      the kernel weights enter it) unless the config pins an explicit grid ----
  Zs  <- scale(Za)                                                     # == the design mgm's node models see
  lam_max <- NA_real_
  LAM <- cfg$mgm_lambda_grid
  if (is.null(LAM)) {
    lam_max <- mgm_lambda_max(Zs, tp, ep, bw)
    LAM <- exp(seq(log(lam_max * (cfg$mgm_lambda_min_ratio %||% 1e-4)),  # glmnet's own n>p default ratio
                   log(lam_max), length.out = cfg$mgm_n_lambda %||% 50L))
  }

  ## ---- pass 1: one tvmgm fit per lambda -> symmetric wadj (ROC) + directed wadjNodewise (refit) ----
  .fit_at <- function(lam) {
    mg <- tryCatch(tvmgm(data = Za, type = rep("g", P), level = rep(1, P),
                         timepoints = tp, estpoints = ep, bandwidth = bw, k = 2,
                         lambdaSeq = lam, lambdaSel = "EBIC", lambdaGam = gamma, ruleReg = "OR",
                         threshold = "none", scale = TRUE, pbar = FALSE, signInfo = FALSE),
                   error = function(e) NULL)
    if (is.null(mg)) return(NULL)                                      # excluded, never a zero graph
    zap <- function(W) { W[abs(W) <= 1e-10] <- 0; W }
    list(sym = lapply(1:m, function(k) zap(mg$tvmodels[[k]]$pairwise$wadj[1:P, 1:P])),
         dir = lapply(1:m, function(k) zap(mg$tvmodels[[k]]$pairwise$wadjNodewise[1:P, 1:P])))
  }
  # embarrassingly parallel over lambda; with array_by="task" the outer mclapply holds one
  # task per array element, so this inner parallelism does not oversubscribe.
  ncore <- { s <- Sys.getenv("SLURM_CPUS_PER_TASK"); if (nzchar(s)) as.integer(s) else max(1L, parallel::detectCores() - 1L) }
  raw <- if (ncore > 1L) parallel::mclapply(LAM, .fit_at, mc.cores = ncore, mc.preschedule = FALSE)
         else lapply(LAM, .fit_at)

  ok <- !vapply(raw, is.null, logical(1))
  if (!any(ok)) return(list(method = "tvmgm", error = "every lambda fit failed"))   # -> master .FAILED marker
  LAMok         <- LAM[ok]
  Om_by_lambda  <- lapply(raw[ok], `[[`, "sym")
  Dir_by_lambda <- lapply(raw[ok], `[[`, "dir")

  rp     <- roc_from_lambda_path(Om_by_lambda, LAMok, dat$true_Omega_list, eval_slices)
  scores <- entry_lambda_scores(Om_by_lambda, LAMok, m)                # per-edge entry-lambda score

  ## ---- per-lambda detail, ascending lambda, lambda-aligned FPR/TPR (deployed point = direct read) ----
  detail <- do.call(rbind, lapply(seq_along(LAMok), function(j) {
    Oml <- Om_by_lambda[[j]]; ft <- edge_fpr_tpr(Oml, dat$true_Omega_list, eval_slices)
    ne  <- sum(vapply(Oml, function(W) { d <- W; diag(d) <- 0; sum(d[upper.tri(d)] != 0) }, integer(1)))
    data.frame(lambda = LAMok[j], n_edges = ne, df = ne,
               FPR = unname(ft["FPR"]), TPR = unname(ft["TPR"]), status = "ok", row.names = NULL)
  }))
  if (any(!ok)) detail <- rbind(detail, data.frame(lambda = LAM[!ok], n_edges = NA_integer_, df = NA_integer_,
                                 FPR = NA_real_, TPR = NA_real_, status = "fit_error", row.names = NULL))
  detail <- detail[order(detail$lambda), ]; row.names(detail) <- NULL

  ## ---- mgm's NATIVE penalised per-node EBIC point ("pre" = the method as published) ----
  ## one full-grid EBIC fit, computed once and reused for `native`, the refit's `pre` baseline,
  ## and the deployed point when the refit is switched off.
  native <- NULL
  mgE <- tryCatch(tvmgm(data = Za, type = rep("g", P), level = rep(1, P),
                        timepoints = tp, estpoints = ep, bandwidth = bw, k = 2,
                        lambdaSeq = LAM, lambdaSel = "EBIC", lambdaGam = gamma, ruleReg = "OR",
                        threshold = "none", scale = TRUE, pbar = FALSE, signInfo = FALSE),
                  error = function(e) NULL)
  if (!is.null(mgE)) {
    OmE <- lapply(1:m, function(k) { W <- mgE$pairwise$wadj[1:P, 1:P, k]; W[abs(W) <= 1e-10] <- 0; W })
    cm  <- .mgm_confusion(OmE, dat$true_Omega_list, eval_slices)
    native <- list(selector = "EBIC_penalised", FPR = unname(cm["FPR"]), TPR = unname(cm["TPR"]),
                   precision = unname(cm["PR"]), F1 = unname(cm["F1"]), edges = unname(cm["edges"]),
                   wadj = OmE)
  }

  ## ---- pass 2: relaxed refit on the standardised design mgm itself uses (scale=TRUE) ----
  refit <- NULL
  if (isTRUE(cfg$mgm_refit)) {
    Wk <- lapply(ep, function(e) exp(-0.5 * ((tp - e) / bw)^2))        # tvmgm's Gaussian kernel weights (Zs built above)
    refit <- tryCatch(
      mgm_refit_pernode(Zs, Wk, Dir_by_lambda, Om_by_lambda, LAMok,
                        dat$true_Omega_list, eval_slices, gamma),
      error = function(e) list(error = conditionMessage(e)))
    if (is.null(refit$error)) {
      refit$pre <- native                                              # penalised baseline for the comparison
      refit$scale_sd <- as.numeric(attr(Zs, "scaled:scale"))           # provenance of the standardisation
    }
  }

  ## ---- as-deployed operating point = mgm's native selector, relaxed (refit-EBIC) ----
  ##      the "EBIC selects too sparse" fix applied to mgm; falls back to the penalised
  ##      EBIC point if the refit is off or failed.
  deployed <- if (!is.null(refit) && is.null(refit$error)) {
    d <- refit$refit
    list(selector = "refit-EBIC", bw = bw, FPR = d$FPR, TPR = d$TPR,
         precision = d$precision, F1 = d$F1, n_edges = d$edges)
  } else if (!is.null(native)) {
    list(selector = "EBIC_penalised", bw = bw, FPR = native$FPR, TPR = native$TPR,
         precision = native$precision, F1 = native$F1, n_edges = native$edges)
  } else NULL

  list(method = "tvmgm", convention = "lambda-sweep",
       # ---- SHELL (shared across the four methods) ----
       roc = rp$roc, auc = rp$auc, scores = scores, detail = detail,
       deployed = deployed, native = native,
       n_per_slice = sapply(dat$Z_0, nrow), N = nrow(Za),
       tuning = list(bw = bw, bw_sel = bw_sel, bw_source = bw_source, tuned = tuned,
                     bwSeq = bwSeq, lambdaGam = gamma, ruleReg = "OR", threshold = "none",
                     scale = TRUE, k = 2L,
                     lambda_max = lam_max,                             # NA if the grid was pinned by config
                     lambda_adaptive = is.null(cfg$mgm_lambda_grid),
                     lambda_min_ratio = cfg$mgm_lambda_min_ratio %||% 1e-4, n_lambda = length(LAM)),
       # ---- PAYLOAD: mgm's native fitted objects, both passes, every successful lambda ----
       estimate = list(
         base  = list(wadj = Om_by_lambda, wadjNodewise = Dir_by_lambda, lambda = LAMok),
         refit = refit))
}
