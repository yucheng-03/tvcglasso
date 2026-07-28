# ============================================================================
# methods/method_jgl.R — Joint Graphical Lasso (Danaher et al. 2014), m<=7 SECONDARY baseline.
#   input = CENTERED ALR (dat$Z_0; the plug-in Gaussian data — JGL has NO LNM layer / no Zhat)
#   base  = JGL(penalty="fused", lambda1 x lambda2, return.whole.theta=TRUE); edge = |theta| > tol
#
# refit (pass 2, cfg$jgl_refit) — the STANDARD relaxed graphical lasso / refitted-MLE: drop BOTH
# lambda1 + lambda2 => per-slice UNPENALIZED constrained Gaussian MLE (Dempster covariance
# selection) on the fused-selected support, glasso(S_k, rho=0, zero=off-support), FREE diagonal.
# Existence is governed by whether the graph-restricted S_k admits a PD completion (Grone et al.
# 1984; equivalently n_k - 1 >= MLT(E_k) — n_k MINUS ONE, since Y is centered — Buhl 1993): sparse (deployed) supports exist even at n<P,
# but the DENSE END of the lambda1 sweep does not, and there the unridged solve does not merely
# fail — it does NOT TERMINATE (measured: > 2 h 34 min on ONE 15x15 slice). We therefore solve
# pass 2 on a RIDGED slice covariance, S_k + eps*mean(diag(S_k))*I with eps = 1e-3
# (cfg$jgl_refit_ridge_eps) — ★ OURS, DISCLOSED, not part of Danaher's JGL; it is a penalty on the
# DIAGONAL of Omega only, so the selected off-diagonals stay unpenalized and the de-biasing is
# intact, and it is confined to the solver (the likelihood is scored on the UNRIDGED S_k). Any
# remaining point whose solve is not finite+PD is FLAGGED refit_exists = FALSE and excluded from
# the deployment candidates — never silently substituted. Full derivation, literature and the
# eps calibration table: R/baselines/refit_jgl.R. Refit does NOT move the support -> ROC
# unchanged; its role is the de-biased DEPLOYED point.
#
# DEPLOYED SELECTOR = AIC (cfg$jgl_selector). AIC is the approximate selector Danaher 2014 §6 /
# Eq.(6.21) evaluates on the PENALIZED JGL fit; the JGL R package ships NO auto-selector, and the
# REFIT-AIC below is OUR extension (not "JGL-native"). The FULL (l1,l2) grid records base & refit
# {neg2loglik, AIC, BIC, eBIC} + FPR/TPR/F1/df + refit-existence, so ANY selector's deployed point
# and the base-vs-refit comparison are recomputable at analysis with NO re-fit. JGL is a JOINT
# method: ONE (l1,l2) fits all m slices jointly, so the deployed point is ONE grid point selected
# by the JOINT (summed-over-slice) IC — NOT per-slice-independent (that is CGLasso).
#
# LAMBDA1 — DATA-ADAPTIVE (matches TV/CGLasso): lambda1_max = 2 * max over slices of the max
# off-diagonal |S_k| (>= the empty-graph threshold => the sparse/left-bottom corner is reached
# with a REAL lambda1), lambda1_min = lambda1_max/1e4 (2026-07-25: /1e3 did not reach FPR=1), jgl_n_lambda1 (default 50) log-spaced
# points — the SAME 50-point path resolution all four methods use (2026-07-24 decision), and the
# same recipe as CGLasso's exp(seq(log(2 maxSigma), log(maxSigma/1000), n)). NB JGL is a JOINT
# method: ONE lambda1 governs ALL m slices, so the grid is anchored on the max ACROSS slices
# (CGLasso, being per-slice independent, anchors each slice's own grid).
# cfg$jgl_lambda1_grid overrides. LAMBDA2 = the fusion grid (cfg$jgl_lambda2_grid).
#
# REPORTED ROC — lambda2 is JGL's EXTRA (fused) knob; to trace one curve we give JGL its BEST
# fusion level (a deliberate steelman: maximise the competitor, then show it still loses). The
# lambda2 is chosen by AGGREGATING SEEDS then taking the single best lambda2 — an ANALYSIS-time
# step, from $per_l2 (this seed's per-lambda2 AUC) pooled across seeds. The worker's roc/auc here
# is a WITHIN-SEED reference = this seed's best-lambda2 curve; the paper's JGL curve is the
# across-seed best-lambda2. $envelope = the best-over-(l1,l2) oracle upper bound (truth-aware,
# LABELED, sim-only). scores = per-edge ENTRY-lambda1 at this seed's best-lambda2.
#
# ── OUTPUT: the SHARED SHELL (identical layout across all four methods) ───────
#   $method $convention $roc $auc $scores $detail $deployed $n_per_slice $N $estimate
#   $detail   FLAT data.frame, one row per (lambda1,lambda2) path point, path-ALIGNED with
#             $estimate: lambda1 lambda2 df n_edges FPR TPR precision F1 |
#             neg2loglik AIC BIC eBIC            <- POST-refit (de-biased; the deployed criterion)
#             neg2loglik_pre AIC_pre BIC_pre eBIC_pre  <- PRE-refit (penalized base)
#             refit_exists base_pd               <- JGL-specific validity flags
#   $deployed the as-deployed point at the NATIVE selector, a DIRECT READ (no recompute):
#             selector lambda_index lambda1 lambda2 FPR TPR precision F1 n_edges df refit_exists
#   $estimate the PAYLOAD slot (shared NAME, method-specific CONTENT): base$Theta / refit$Theta,
#             each a length-nrow(grid) list of m precision matrices, row-aligned with $detail.
#             JGL has Theta ONLY — no beta (no basis) and no Z (no latent layer).
#   JGL-specific siblings: $per_l2 $best_l2_this_seed $envelope $env_auc $lambda1_grid
#             $lambda2_grid $refit (per-selector pre/post deployed indices + oracle) $tuning.
#
# References: Danaher-Wang-Witten 2014 (JRSS-B 76:373) JGL; Dempster 1972 covariance selection;
# Buhl 1993 MLE existence; Meinshausen 2007 relaxed lasso. Requires: JGL, glasso, parallel,
# R/roc_utils.R, R/baselines/refit_jgl.R.
# ============================================================================

run_method_jgl <- function(dat, cfg, eval_slices) {
  m <- length(dat$X); P <- ncol(dat$Z_0[[1]])
  # ── m<=7 HARD STOP: JGL's fused ADMM is a small-m secondary baseline; INFEASIBLE at m=15/30.
  #    ★ MEASURED 2026-07-26 — the quantitative basis for that claim (the earlier basis was only
  #    "a 2026-06-07 m=15 attempt was killed after 6.5 h", whose cause was never diagnosed; this
  #    replaces it). SAME cell (P=15, n=12, low depth), SAME 50x6 grid, SAME 16 cores, the ONLY
  #    variable is m:
  #        m =  7 :  71.00 grid points / min  ->    4.2 min per (cell,seed)
  #        m = 15 :   1.11 grid points / min  ->  ~300   min per (cell,seed)     [64x SLOWER]
  #    and the m=15 rate is FLAT across the run (per-point timestamps stay 1-3 min apart — it does
  #    not warm up). Decomposition: per-ITERATION cost grows only ~2.1x (one eigendecomposition per
  #    slice, O(m P^3)); the remaining ~30x is the ADMM ITERATION COUNT blowing up, because the
  #    fused penalty must reach consensus across all m slices and doubling the slice count makes
  #    that far harder. Extrapolated to the full design: ~64,000 core-hours for 8 cells x 100 seeds
  #    (vs 1,804 at m=7) = ~64 h wall-clock on 1,000 cores, with every task exceeding a 4 h
  #    walltime. Hence the 8 m=15 cells are covered by tvcglasso/CGLasso/tvmgm only.
  #    `jgl_allow_large_m=TRUE` exists solely to reproduce this diagnostic — NOT for a campaign.
  #    Refuse m>7 BEFORE any fit unless that explicit, named override is set. ──
  if (m > 7L && !isTRUE(cfg$jgl_allow_large_m))
    stop(sprintf("JGL is an m<=7 secondary baseline (m=%d); set cfg$jgl_allow_large_m=TRUE to override.", m))

  Y <- lapply(dat$Z_0, function(z) scale(z, center = TRUE, scale = FALSE))   # centered ALR per slice (plug-in Gaussian)
  nk <- vapply(Y, nrow, integer(1))                                          # per-slice n
  S_list <- lapply(seq_len(m), function(k) crossprod(Y[[k]]) / nk[k])        # 1/n MLE cov (Y centered)
  true <- dat$true_Omega_list

  # LAMBDA1 — data-adaptive (unless supplied): lambda1_max = 2 * max over slices of max|off-diag S_k|
  # (>= empty-graph threshold), down to /1e4, jgl_n_lambda1 log-spaced points.
  LAM1 <- cfg$jgl_lambda1_grid
  if (is.null(LAM1)) {
    max_offdiag <- max(vapply(S_list, function(S) max(abs(S[upper.tri(S)])), numeric(1)))
    l1_max <- 2 * max_offdiag
    # DENSE END = l1_max / 1e4 (2026-07-25, user decision; was /1e3). The /1e3 end did NOT reach
    # the top-right ROC corner with REAL operating points -- measured at the smallest lambda1 on
    # the finished lambda2 blocks: FPR 0.943 / 0.951 / 0.980 / 0.983 / 0.990, i.e. 1-6 points short
    # of FPR = 1, and maxFPR was attained EXACTLY at the smallest lambda1 (so the path was monotone
    # there: the shortfall is the grid's lower bound, not a missing point). The project rule is to
    # close a corner gap with MORE REAL lambda, never by interpolating to (1,1).
    LAM1 <- exp(seq(log(l1_max), log(l1_max / 1e4), length.out = cfg$jgl_n_lambda1 %||% 50L))
  }
  LAM2 <- cfg$jgl_lambda2_grid %||% cfg$jgl_lambda2 %||% c(0, 0.02, 0.05, 0.10, 0.20, 0.40)
  EBIC_G  <- cfg$jgl_ebic_gamma %||% 0.25
  EDGE_TOL <- cfg$jgl_edge_tol %||% 1e-5                                      # JGL's OWN truncate default = the honest support tolerance
  REFIT_MAXIT <- cfg$jgl_refit_maxit %||% 500L                               # glasso(rho=0) backstop, SAME as the CGLasso refit (thr stays glasso's own default 1e-4 — never hand-tightened)
  # ★ OURS, disclosed: pass-2 ridge on the slice covariance, guaranteeing the given-graph MLE
  #   exists at every grid point (see the header block of R/baselines/refit_jgl.R). JGL's OWN key
  #   — never shares cglasso_refit_ridge_eps. 0 = exact unridged MLE (and the dense-end stalls).
  REFIT_RIDGE <- cfg$jgl_refit_ridge_eps %||% 1e-3
  grid <- expand.grid(l1 = LAM1, l2 = LAM2, KEEP.OUT.ATTRS = FALSE)
  do_refit <- isTRUE(cfg$jgl_refit %||% TRUE)
  sel <- cfg$jgl_selector %||% "AIC"
  if (!sel %in% c("AIC", "BIC", "eBIC"))    # fail NOW, not after 300 fits (JASA review #11)
    stop(sprintf("cfg$jgl_selector must be one of AIC/BIC/eBIC; got \"%s\".", sel))

  # EXPLICIT + recorded JGL package params (nothing left to a silent default; provenance).
  # NB base maxiter = 500 = the author's JGL default (faithful); it is SEPARATE from REFIT_MAXIT (glasso).
  jgl_pkg <- list(rho = cfg$jgl_rho %||% 1, weights = cfg$jgl_weights %||% "equal",
                  penalize.diagonal = FALSE, maxiter = cfg$jgl_maxiter %||% 500L,
                  tol = cfg$jgl_tol %||% 1e-5, truncate = cfg$jgl_truncate %||% 1e-5,
                  screening = cfg$jgl_screening %||% "fast")
  jgl_call <- function(l1, l2) JGL::JGL(Y = Y, penalty = "fused", lambda1 = l1, lambda2 = l2,
                                        rho = jgl_pkg$rho, weights = jgl_pkg$weights,
                                        penalize.diagonal = jgl_pkg$penalize.diagonal,
                                        maxiter = jgl_pkg$maxiter, tol = jgl_pkg$tol,
                                        truncate = jgl_pkg$truncate, screening = jgl_pkg$screening,
                                        return.whole.theta = TRUE)
  supp_of <- function(fit) lapply(fit$theta, function(O) { O[abs(O) <= EDGE_TOL] <- 0; O })  # edge = |theta|>tol

  # per-(l1,l2) checkpoint (discipline #0 / preempt-safe): each point atomically saved + skipped on resume.
  NC <- { s <- Sys.getenv("SLURM_CPUS_PER_TASK"); if (nzchar(s)) max(1L, as.integer(s)) else 1L }  # 1 locally: no nested fork-storm (matches CGLasso)
  ckdir <- if (!is.null(dat$refit_ckpt_dir)) file.path(dat$refit_ckpt_dir, "jgl") else NULL
  if (!is.null(ckdir)) dir.create(ckdir, showWarnings = FALSE, recursive = TRUE)

  # ── CHECKPOINT IDENTITY: everything that DEFINES the stored point. A cached point is reused only
  #    if its (lambda1,lambda2) AND this stamp match bit-for-bit; otherwise it is RECOMPUTED.
  #    Matching on (l1,l2) alone would silently resurrect points computed under a DIFFERENT
  #    estimator or a different detail schema — this repo has already changed both (the pass-2
  #    ridge changes the refit estimator; the detail columns were renamed), and mixing them would
  #    produce a table of mixed provenance that still rbind()s cleanly. Same guard, same reason as
  #    the CGLasso refit's `config` check (refit_cglasso_core.R:467-473). ──
  CK_STAMP <- list(schema = "jgl-detail-v2", edge_tol = EDGE_TOL, ebic_gamma = EBIC_G,
                   do_refit = do_refit, refit_maxit = REFIT_MAXIT, refit_ridge_eps = REFIT_RIDGE,
                   jgl_pkg = jgl_pkg, eval_slices = eval_slices, P = P, m = m, n_per_slice = nk)

  # ── detail ROW = the SHARED per-path-point contract, JGL-flavoured. Column names follow the
  #    tvcglasso/CGLasso `detail` layout so a client reads any method's table the same way:
  #    UNSUFFIXED {neg2loglik,AIC,BIC,eBIC} = the POST-refit (de-biased, deployed) criterion;
  #    *_pre = the PRE-refit penalized base. JGL-SPECIFIC: the path point is 2-D (lambda1 x
  #    lambda2, its extra fusion knob) and it carries validity flags (refit_exists = the
  #    constrained MLE exists at this n [Buhl 1993]; base_pd = the ADMM Theta is PD). ──
  na_row <- function(l1, l2) data.frame(
      lambda1 = l1, lambda2 = l2, df = NA_integer_, n_edges = NA_integer_,
      FPR = NA_real_, TPR = NA_real_, precision = NA_real_, F1 = NA_real_,
      neg2loglik = NA_real_, AIC = NA_real_, BIC = NA_real_, eBIC = NA_real_,
      neg2loglik_pre = NA_real_, AIC_pre = NA_real_, BIC_pre = NA_real_, eBIC_pre = NA_real_,
      refit_exists = NA, base_pd = NA)

  fit_one <- function(r) {
    l1 <- grid$l1[r]; l2 <- grid$l2[r]
    ckf <- if (!is.null(ckdir)) file.path(ckdir, sprintf("pt_%04d.rds", r)) else NULL
    if (!is.null(ckf) && file.exists(ckf)) {                                  # resume: reuse ONLY a bit-compatible checkpoint
      cached <- tryCatch(readRDS(ckf), error = function(e) NULL)
      if (!is.null(cached) && identical(cached$stamp, CK_STAMP) &&
          isTRUE(all.equal(cached$row$lambda1, l1)) && isTRUE(all.equal(cached$row$lambda2, l2)))
        return(cached)
    }
    fit <- tryCatch(jgl_call(l1, l2), error = function(e) NULL)
    if (is.null(fit)) { out <- list(row = na_row(l1, l2), Theta = NULL, Theta_refit = NULL, stamp = CK_STAMP)
                        if (!is.null(ckf)) .atomic_save(out, ckf); return(out) }
    th <- supp_of(fit)
    cm <- .jgl_confusion(th, true, eval_slices)
    d_per_slice <- vapply(th, function(O) sum(O[upper.tri(O)] != 0), integer(1)); D <- sum(d_per_slice)
    base_ll <- jgl_gaussian_neg2ll(S_list, th, nk)                           # PD-safe (NA if any base slice non-PD)
    b <- jgl_slice_ic(base_ll, nk, d_per_slice, P, EBIC_G)
    row <- na_row(l1, l2)
    row$df <- D; row$n_edges <- D                                            # df = #free off-diag params = #edges (diagonal is a constant offset)
    row$FPR <- cm["FPR"]; row$TPR <- cm["TPR"]; row$precision <- cm["precision"]; row$F1 <- cm["F1"]
    row$neg2loglik_pre <- b["neg2loglik"]; row$AIC_pre <- b["AIC"]
    row$BIC_pre <- b["BIC_slice"]; row$eBIC_pre <- b["eBIC_slice"]
    row$base_pd <- !is.na(base_ll)
    th_ref <- NULL
    if (do_refit) {                                                          # rho=0 constrained MLE, used where it exists (Buhl 1993)
      th_ref <- jgl_refit_supports(S_list, th, maxit = REFIT_MAXIT, ridge_eps = REFIT_RIDGE)
      exists_all <- all(!vapply(th_ref, is.null, logical(1)))
      row$refit_exists <- exists_all
      if (exists_all) {
        rf <- jgl_slice_ic(jgl_gaussian_neg2ll(S_list, th_ref, nk), nk, d_per_slice, P, EBIC_G)
        row$neg2loglik <- rf["neg2loglik"]; row$AIC <- rf["AIC"]
        row$BIC <- rf["BIC_slice"]; row$eBIC <- rf["eBIC_slice"]
      } else th_ref <- NULL                                                  # partial refit is not a usable estimate -> store nothing
    }
    rownames(row) <- NULL
    # PAYLOAD (per grid point): BOTH passes' fitted precision matrices, one list of m per point.
    # Completeness over size (project rule) — and it makes the deployed points a LOOKUP into
    # `estimate` rather than 3 extra JGL re-fits at assembly time.
    out <- list(row = row, Theta = th, Theta_refit = th_ref, stamp = CK_STAMP)
    if (!is.null(ckf)) .atomic_save(out, ckf)
    out
  }

  # ── LAMBDA2 STRIPING (cfg$jgl_l2_stripe = c(i, n)) — OPTIONAL, PURELY A SCHEDULING SPLIT ──
  #    lambda2 values are INDEPENDENT fits, so one (cell,seed)'s 50 x |LAM2| grid can be spread
  #    over n concurrent SLURM tasks: task i computes only the lambda2 columns with
  #    ((j - 1) %% n) + 1 == i and writes them to the SHARED per-point checkpoint dir. Whichever
  #    striped task finishes LAST finds all points cached, assembles, and returns the full result
  #    (so no separate "assembler" job is needed and no harness change is required); the earlier
  #    ones return a partial marker, which the harness records without blocking a resubmit and
  #    WITHOUT deleting the checkpoints. Striping changes NOTHING about the estimator or the
  #    output: every point is computed by the same fit_one and keyed by the same CK_STAMP, and
  #    the assembled result is bit-identical to an unstriped run.
  #    NULL (the default) = no striping = the whole grid in this task, exactly as before.
  STRIPE <- cfg$jgl_l2_stripe
  l2_index_of_row <- ((seq_len(nrow(grid)) - 1L) %/% length(LAM1)) + 1L      # expand.grid varies l1 fastest
  if (!is.null(STRIPE)) {
    if (is.null(ckdir))
      stop("cfg$jgl_l2_stripe requires a checkpoint dir (dat$refit_ckpt_dir): striped tasks share state through it.")
    if (length(STRIPE) != 2L || STRIPE[1] < 1L || STRIPE[1] > STRIPE[2])
      stop(sprintf("cfg$jgl_l2_stripe must be c(i, n) with 1 <= i <= n; got %s", paste(STRIPE, collapse = ",")))
    rows_todo <- which(((l2_index_of_row - 1L) %% STRIPE[2]) + 1L == STRIPE[1])
    cat(sprintf("[JGL] lambda2 stripe %d/%d: %d of %d grid points in this task\n",
                STRIPE[1], STRIPE[2], length(rows_todo), nrow(grid)))
  } else rows_todo <- seq_len(nrow(grid))

  computed <- parallel::mclapply(rows_todo, fit_one, mc.cores = NC, mc.preschedule = FALSE)
  names(computed) <- as.character(rows_todo)

  # Assemble the FULL grid: this task's own points, plus any point another stripe has already
  # checkpointed. A point is usable only if its checkpoint carries a matching CK_STAMP.
  res <- lapply(seq_len(nrow(grid)), function(r) {
    x <- computed[[as.character(r)]]
    if (is.null(x) && !is.null(ckdir)) {
      f <- file.path(ckdir, sprintf("pt_%04d.rds", r))
      if (file.exists(f)) {
        cached <- tryCatch(readRDS(f), error = function(e) NULL)
        if (!is.null(cached) && identical(cached$stamp, CK_STAMP)) x <- cached
      }
    }
    if (inherits(x, "try-error") || !is.list(x) || is.null(x$row)) NULL else x
  })

  missing <- which(vapply(res, is.null, logical(1)))
  if (!is.null(STRIPE) && length(missing) > 0L)                             # other stripes not done yet
    return(list(method = "JGL", stripe = STRIPE, n_points_done = nrow(grid) - length(missing),
                n_points_total = nrow(grid),
                error = sprintf("lambda2 stripe %d/%d complete (%d/%d grid points present); awaiting the other stripes — this is BY DESIGN, not a failure. The last stripe to finish assembles the result. Checkpoints are kept.",
                                STRIPE[1], STRIPE[2], nrow(grid) - length(missing), nrow(grid))))
  res[missing] <- lapply(missing, function(r)                                # unstriped: a crashed point -> NA row
    list(row = na_row(grid$l1[r], grid$l2[r]), Theta = NULL, Theta_refit = NULL))
  gdf <- do.call(rbind, lapply(res, `[[`, "row")); rownames(gdf) <- NULL   # = the FLAT `detail` table

  # ── per-lambda2 AUC (this seed). The reported curve chooses lambda2 by AGGREGATING SEEDS at
  #    analysis (from per_l2); the worker's roc/auc below is a WITHIN-SEED reference = this seed's
  #    best-lambda2 curve (NOT a per-seed truth-tuned deployable estimator — it is labeled as such). ──
  per_l2 <- do.call(rbind, lapply(sort(unique(gdf$lambda2)), function(L) {
    d <- gdf[abs(gdf$lambda2 - L) < 1e-12 & is.finite(gdf$FPR), c("FPR", "TPR")]; d <- d[order(d$FPR), ]
    data.frame(lambda2 = L, auc = if (nrow(d) >= 2) auc_trap(d) else NA_real_, n_pts = nrow(d)) }))
  best_l2 <- if (all(is.na(per_l2$auc))) per_l2$lambda2[1] else per_l2$lambda2[which.max(per_l2$auc)]

  mm <- abs(gdf$lambda2 - best_l2) < 1e-12 & is.finite(gdf$FPR)
  main <- gdf[mm, c("FPR", "TPR", "lambda1")]; main <- main[order(main$FPR), ]; rownames(main) <- NULL
  # A 0-/1-row `main` would still pass the harness's validity test (`!is.null(r$roc)` is TRUE for a
  # 0-row data.frame), so a TOTAL grid failure would be written out as a VALID result with AUC=NA,
  # its checkpoints deleted and that (cell,seed) skipped forever. Fail LOUDLY instead (JASA review
  # #9): the harness then writes a .FAILED marker and KEEPS the checkpoints for a resume.
  if (nrow(main) < 2L)
    stop(sprintf("JGL: no usable ROC — %d of %d grid points have a finite FPR, and only %d row(s) at the selected lambda2 = %g.",
                 sum(is.finite(gdf$FPR)), nrow(gdf), nrow(main), best_l2))
  main_auc <- auc_trap(main[, c("FPR", "TPR")])
  cov_rep  <- roc_sanity_report(main[, c("FPR", "TPR")], "JGL")

  # entry-lambda1 scores at this seed's best-lambda2 (method-symmetric score-sweep input).
  im <- which(abs(grid$l2 - best_l2) < 1e-12)
  supp_bl2 <- lapply(im, function(r) { Th <- res[[r]]$Theta
                                       if (is.null(Th)) NULL else lapply(Th, function(O) { diag(O) <- 0; O }) })
  keep <- !vapply(supp_bl2, is.null, logical(1))
  scores <- if (any(keep)) tryCatch(entry_lambda_scores(supp_bl2[keep], grid$l1[im][keep], m), error = function(e) NULL) else NULL

  # ── DEPLOYED = argmin IC over the (l1,l2) grid. pre over PD points; post(refit) over points where
  #    the constrained MLE exists (Buhl 1993). pre & post for AIC/BIC/eBIC ALL recorded, each as a
  #    ROW INDEX into `detail` (and hence into `estimate`) => 免re-fit, no duplicated matrices. ──
  pick_min <- function(v) { ok <- is.finite(v); if (!any(ok)) NA_integer_ else which(ok)[which.min(v[ok])] }
  dep_for <- function(ic) {
    ipre <- pick_min(gdf[[paste0(ic, "_pre")]]); ipost <- pick_min(gdf[[ic]])
    list(pre_index = ipre, post_index = ipost,
         pre  = if (is.na(ipre))  NULL else gdf[ipre, ],
         post = if (is.na(ipost)) NULL else gdf[ipost, ],
         n_valid_pre = sum(is.finite(gdf[[paste0(ic, "_pre")]])), n_valid_post = sum(is.finite(gdf[[ic]]))) }
  D_by <- list(AIC = dep_for("AIC"), BIC = dep_for("BIC"), eBIC = dep_for("eBIC"))
  io <- pick_min(-gdf$F1)                                                   # oracle = max-F1 (LABELED, truth-aware)

  # ── `deployed` = the as-deployed operating point at JGL's NATIVE selector, a DIRECT READ
  #    (shared shell contract). JGL-specific: the point is a (lambda1,lambda2) PAIR; lambda_index
  #    indexes BOTH `detail` and `estimate`. refit_exists = the constrained MLE existed there. ──
  di <- D_by[[sel]]$post_index; if (is.na(di)) di <- D_by[[sel]]$pre_index
  deployed <- if (is.na(di)) NULL else { rr <- gdf[di, ]
    list(selector = paste0(if (do_refit && !is.na(D_by[[sel]]$post_index)) "refit-" else "base-", sel),
         lambda_index = di, lambda1 = rr$lambda1, lambda2 = rr$lambda2,
         FPR = rr$FPR, TPR = rr$TPR, precision = rr$precision, F1 = rr$F1,
         n_edges = rr$n_edges, df = rr$df, refit_exists = rr$refit_exists) }

  # oracle best-over-(l1,l2) envelope (LABELED oracle upper bound — truth-aware, sim-only).
  fpr_grid <- sort(unique(c(0, gdf$FPR[is.finite(gdf$FPR)])))
  env_tpr  <- cummax(vapply(fpr_grid, function(fp) max(gdf$TPR[is.finite(gdf$FPR) & gdf$FPR <= fp + 1e-9], 0), numeric(1)))
  env_auc  <- if (length(fpr_grid) >= 2) auc_trap(data.frame(FPR = fpr_grid, TPR = env_tpr)) else NA_real_

  # ── `estimate` = the PAYLOAD slot (shared NAME, method-specific CONTENT). JGL is a Gaussian
  #    method on the plug-in ALR: it has Theta only — NO beta (no basis) and NO Z (no latent
  #    layer), unlike tvcglasso {beta,Omega,Z} / CGLasso {Omega,Z}. BOTH passes at EVERY one of the
  #    nrow(grid) points, row-aligned with `detail` (refit NULL where the MLE does not exist). ──
  estimate <- list(base  = list(Theta = lapply(res, `[[`, "Theta")),
                   refit = list(Theta = lapply(res, `[[`, "Theta_refit")))

  list(method = "JGL", convention = "lambda-sweep",
       roc = main[, c("FPR", "TPR")], auc = main_auc, scores = scores,
       detail = gdf, deployed = deployed,                                   # SHELL: flat table + direct-read point
       n_per_slice = nk, N = sum(nk),
       # ---- JGL-specific siblings (the extra fusion knob + its steelman bookkeeping) ----
       per_l2 = per_l2, best_l2_this_seed = best_l2,
       envelope = data.frame(FPR = fpr_grid, TPR = env_tpr), env_auc = env_auc,
       lambda1_grid = LAM1, lambda2_grid = LAM2,
       refit = list(selector = sel, did_refit = do_refit, by_selector = D_by,
                    oracle_index = io, oracle = if (is.na(io)) NULL else gdf[io, ],
                    n_refit_exist = sum(gdf$refit_exists %in% TRUE), n_grid = nrow(gdf)),
       tuning = list(ebic_gamma = EBIC_G, edge_tol = EDGE_TOL, refit_maxit = REFIT_MAXIT,
                     refit_ridge_eps = REFIT_RIDGE,     # ★ OURS, disclosed — recorded per result
                     jgl_pkg = jgl_pkg, coverage = cov_rep),
       estimate = estimate)
}

# atomic checkpoint save (tmp -> rename) so a preempt mid-write never leaves a corrupt point.
.atomic_save <- function(obj, path) { tmp <- paste0(path, ".tmp"); saveRDS(obj, tmp); file.rename(tmp, path) }
