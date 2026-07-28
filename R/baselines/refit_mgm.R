# ============================================================================
# baselines/refit_mgm.R — tvmgm (Haslbeck & Waldorp 2020) RELAXED nodewise refit.
#
# mgm is a NODEWISE neighbourhood-selection estimator (one L1 GLM per node): no joint
# precision matrix, no latent Z. Like JGL it therefore does NOT route through the LNM
# refit core (R/refit.R) used by tvcglasso/CGLasso -- its relaxed lasso stays inside
# mgm's OWN nodewise framework (the textbook OLS-post-lasso; Belloni-Chernozhukov 2013,
# Meinshausen 2007):
#
#   pass 1 (mgm, penalised)  each node's L1 regression selects its DIRECTED neighbourhood;
#   pass 2 (here, relaxed)   each node is refitted UNPENALISED on ITS OWN neighbourhood by
#                            weighted least squares with the SAME per-estpoint Gaussian
#                            kernel weights tvmgm uses -> de-biased magnitudes; the
#                            operating point is re-selected PER NODE by the relaxed
#                            (post-lasso) nodewise EBIC.
#
# FAITHFULNESS TO mgm (the design rule for this baseline: stay as close to the published
# method as possible; our only addition is removing the penalty in pass 2):
#
#  * SUPPORT = each node's OWN directed lasso mask, taken from the DIRECTED
#    `wadjNodewise` (verified: wadjNodewise[q,p] = |node p's coefficient on q|, i.e.
#    COLUMN p is node p's own regression), NOT the OR-symmetrised `wadj`. A one-sided
#    OR edge is therefore never forced into both endpoints' regressions. The OR combine
#    happens only at the END, on logical masks (mgm's Reg2Graph order).
#
#  * INFORMATION CRITERION = mgm's OWN. mgm::nodeEst computes, per node and lambda,
#        EBIC = -2*LL + d*log(nadj) + 2*gamma*d*log(ncol(X))
#    and for a Gaussian node calcLL/glmnet give  -2*LL = deviance + const, where the
#    Gaussian deviance IS the weighted RSS (verified numerically) and the additive
#    anchor (LL_sat = nulldev/2 + LL_null) does NOT depend on lambda. So mgm's selection
#    is argmin over lambda of
#        RSS_w + d*log(nadj) + 2*gamma*d*log(P-1),
#    i.e. a UNIT-VARIANCE Gaussian deviance (sigma^2 = 1), which is calibrated because
#    mgm standardises every Gaussian column internally (scale=TRUE). We reproduce that
#    criterion EXACTLY -- same functional form, same nadj = sum of kernel weights
#    (tvmgm's Ne), same gamma, same log(P-1) -- on the SAME standardised scale (the
#    caller passes the standardised design; see method_tvmgm.R). The lambda-independent
#    LL_sat anchor is omitted: it shifts every candidate of a node equally and so cannot
#    change that node's argmin. BIC/AIC are recorded in the same -2LL convention.
#    => pre (penalised) vs refit (relaxed) differ ONLY in penalised-vs-unpenalised.
#
#  * AGGREGATION = mgm's Reg2Graph OR rule: the undirected weight is the mean of the two
#    directions' ABSOLUTE strengths, wadj[p,q] = (|b_pq| + |b_qp|)/2 (a one-sided edge
#    keeps half weight), with the signs stored separately. Support is the OR of the two
#    directed logical masks -- never inferred from a floating coefficient being non-zero.
#
#  * SINGULARITY: pass 2's weighted OLS is well-posed almost everywhere here -- all
#    n_total = n*m observations carry positive kernel weight, so the weighted design is
#    generically full column rank and the exact WLS exists. (mgm itself never meets this
#    question: it is always L1-penalised via glmnet.) We therefore compute the exact WLS
#    by Cholesky and add NO ridge: a ridge would silently make the estimator something
#    other than the unpenalised post-lasso it claims to be, and would invalidate the
#    integer-df IC. A genuine Cholesky failure (exact collinearity) marks that candidate
#    NA so it cannot be selected. Diagnostics (d, n_eff, chol_ok, RSS) are recorded, and
#    the 16-cell validation found chol_fail = 0 and no selected candidate with d >= n_eff
#    in ANY cell (including n < P), so no additional filtering rule is needed.
#
#  * COMPLETENESS: the FULL per-(estpoint, node, lambda) grid {RSS, df, neighbours,
#    coefficients, n_eff} is returned, so ANY selector (a different gamma, BIC, AIC, a
#    different IC form) is recomputable at analysis time WITHOUT re-fitting.
#
# Called by R/methods/method_tvmgm.R when cfg$mgm_refit, reusing the base lambda-sweep's
# stored supports (no re-fit of pass 1).
# ============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

# --- weighted OLS of y on X (+ intercept) by Cholesky; NO ridge ---------------
# Returns the coefficients on the NEIGHBOUR columns (intercept dropped), the weighted
# RSS (= the Gaussian deviance mgm's EBIC uses), and the rank diagnostics.
.mgm_wls <- function(y, X, w) {
  n_eff <- sum(w)
  if (is.null(X) || ncol(X) == 0L) {                    # empty neighbourhood -> intercept only
    mu <- sum(w * y) / n_eff
    return(list(coef = numeric(0), rss = sum(w * (y - mu)^2), d = 0L, chol_ok = TRUE))
  }
  d  <- ncol(X); Xd <- cbind(1, X); sw <- sqrt(w); Xw <- Xd * sw
  rc <- tryCatch(chol(crossprod(Xw)), error = function(e) NULL)
  if (is.null(rc))                                      # exact collinearity: unavailable, never ridged
    return(list(coef = rep(NA_real_, d), rss = NA_real_, d = d, chol_ok = FALSE))
  b   <- backsolve(rc, backsolve(rc, crossprod(Xw, y * sw), transpose = TRUE))
  list(coef = b[-1], rss = sum(w * (y - as.numeric(Xd %*% b))^2), d = d, chol_ok = TRUE)
}

# --- per-estpoint post-lasso grid --------------------------------------------
# dir_list: list (over lambda) of the P x P DIRECTED wadjNodewise at THIS estpoint.
# Node p's own neighbourhood at lambda j = COLUMN p. Adjacent lambdas usually share a
# node's neighbourhood, so the weighted OLS is cached across the path.
mgm_refit_slice_grid <- function(Za, w, dir_list, P) {
  n_eff <- sum(w); n_lam <- length(dir_list); chol_fail <- 0L
  node <- lapply(seq_len(P), function(p) {
    rss <- numeric(n_lam); df <- integer(n_lam); ok <- logical(n_lam)
    coefs <- vector("list", n_lam); Nb <- vector("list", n_lam)
    prev_key <- NULL; cache <- NULL
    for (j in seq_len(n_lam)) {
      nb  <- which(dir_list[[j]][, p] != 0)             # COLUMN p = node p's OWN directed mask
      key <- paste(nb, collapse = ",")
      if (!identical(key, prev_key)) {
        cache <- .mgm_wls(Za[, p], if (length(nb)) Za[, nb, drop = FALSE] else NULL, w)
        prev_key <- key
        if (!isTRUE(cache$chol_ok)) chol_fail <<- chol_fail + 1L
      }
      rss[j] <- cache$rss; df[j] <- length(nb); coefs[[j]] <- cache$coef
      Nb[[j]] <- nb; ok[j] <- cache$chol_ok
    }
    list(rss = rss, df = df, coefs = coefs, Nb = Nb, chol_ok = ok)
  })
  list(node = node, n_eff = n_eff, chol_fail = chol_fail)
}

# --- nodewise IC over the lambda grid, in mgm's OWN -2LL convention -----------
# -2*LL = the weighted Gaussian deviance = RSS (unit variance, standardised design),
# matching mgm::nodeEst; the lambda-independent LL_sat anchor is omitted (it cannot
# change a node's argmin). EBIC is mgm's native selector; BIC/AIC use the same -2LL.
.mgm_node_ic <- function(rss, df, n_eff, P, gamma, penalty) {
  switch(penalty,
         EBIC = rss + df * log(n_eff) + 2 * gamma * df * log(max(P - 1, 1)),
         BIC  = rss + df * log(n_eff),
         AIC  = rss + 2 * df)
}

# --- per-node selection under one penalty -> de-biased symmetric graph --------
# Each node minimises its OWN relaxed IC (ties -> the larger/sparser lambda). The
# directed refit coefficients are assembled COLUMN-wise (column p = node p's fit), then
# aggregated by mgm's OR rule: mean of the two ABSOLUTE directed strengths, signs apart.
mgm_select_pernode <- function(grid, LAM, P, penalty, gamma) {
  B <- matrix(0, P, P); selj <- rep(NA_integer_, P); sel_d <- rep(NA_integer_, P)
  for (p in seq_len(P)) {
    ic <- .mgm_node_ic(grid$node[[p]]$rss, grid$node[[p]]$df, grid$n_eff, P, gamma, penalty)
    ok <- is.finite(ic); if (!any(ok)) next
    best <- min(ic[ok]); idx <- which(ic == best & ok); j <- idx[which.max(LAM[idx])]
    selj[p] <- j; nb <- grid$node[[p]]$Nb[[j]]; sel_d[p] <- length(nb)
    if (length(nb)) B[nb, p] <- grid$node[[p]]$coefs[[j]]     # COLUMN p = node p's coefficients
  }
  A <- abs(B); W <- (A + t(A)) / 2; diag(W) <- 0               # mgm OR rule: mean absolute strength
  S <- (B + t(B)) / 2; diag(S) <- 0                            # signed companion (magnitudes only)
  list(wadj = W, signed = S, sel_lambda_index = selj, sel_d = sel_d, n_eff = grid$n_eff)
}

# confusion of a per-slice network list against truth, pooled over eval_slices (edge = != 0)
.mgm_confusion <- function(Om_list, true_list, eval_slices) {
  P <- nrow(Om_list[[1]]); ut <- upper.tri(matrix(0, P, P))
  TP <- FP <- FN <- TN <- 0
  for (k in eval_slices) {
    pe <- Om_list[[k]][ut] != 0; te <- true_list[[k]][ut] != 0
    TP <- TP + sum(pe & te);  FP <- FP + sum(pe & !te)
    FN <- FN + sum(!pe & te); TN <- TN + sum(!pe & !te)
  }
  TPR <- TP / max(TP + FN, 1); FPR <- FP / max(FP + TN, 1); PR <- TP / max(TP + FP, 1)
  F1  <- if (TPR + PR > 0) 2 * TPR * PR / (TPR + PR) else 0
  c(TPR = TPR, FPR = FPR, PR = PR, F1 = F1, edges = TP + FP)
}

.mgm_pick <- function(v, lam, maximize = FALSE) {
  ok <- is.finite(v); if (!any(ok)) return(NA_integer_)
  best <- if (maximize) max(v[ok]) else min(v[ok])
  idx <- which(v == best & ok); idx[which.max(lam[idx])]
}

# ============================================================================
# mgm_refit_pernode(Zs, Wk, Dir_by_lambda, Om_by_lambda, LAM, true_Omega_list,
#                   eval_slices, gamma)
#   Zs            = the STANDARDISED pooled ALR design (n_total x P) -- the same scale
#                   mgm's own node models see (scale=TRUE); supplied by the caller.
#   Wk            = list (length m) of per-estpoint Gaussian kernel weights (tvmgm's).
#   Dir_by_lambda = list (ascending lambda) of per-slice P x P DIRECTED wadjNodewise.
#   Om_by_lambda  = list (ascending lambda) of per-slice symmetric base supports (used
#                   only for the labelled oracle ceiling).
#   LAM           = the ascending lambda grid of the SUCCESSFUL base fits.
# Returns the deployed landings (refit = mgm-native relaxed EBIC, plus BIC/AIC/oracle),
# the de-biased networks, the FULL per-(k,node,lambda) grid, and the rank diagnostics.
# ============================================================================
mgm_refit_pernode <- function(Zs, Wk, Dir_by_lambda, Om_by_lambda, LAM,
                              true_Omega_list, eval_slices, gamma = 0.25) {
  P <- ncol(Zs); m <- length(Wk)
  grid_by_k  <- lapply(seq_len(m), function(k)
    mgm_refit_slice_grid(Zs, Wk[[k]], lapply(Dir_by_lambda, `[[`, k), P))
  chol_fail  <- sum(vapply(grid_by_k, function(g) g$chol_fail, integer(1)))
  n_eff_by_k <- vapply(grid_by_k, function(g) g$n_eff, numeric(1))
  # how often the RSS-collapse regime (support as large as the effective sample) is even
  # REACHABLE on the grid -- the quantity that would justify an extra filtering rule
  grid_overfit <- sum(vapply(grid_by_k, function(g)
    sum(vapply(g$node, function(nd) sum(nd$df >= g$n_eff, na.rm = TRUE), integer(1))), integer(1)))

  strip  <- function(d) list(FPR = d$FPR, TPR = d$TPR, precision = d$PR, F1 = d$F1, edges = d$edges)
  deploy <- function(penalty) {
    picks <- lapply(seq_len(m), function(k) mgm_select_pernode(grid_by_k[[k]], LAM, P, penalty, gamma))
    nets  <- lapply(picks, `[[`, "wadj")
    cm <- .mgm_confusion(nets, true_Omega_list, eval_slices)
    dg <- do.call(rbind, lapply(seq_len(m), function(k) data.frame(
      penalty = penalty, k = k, node = seq_len(P), sel_lambda_index = picks[[k]]$sel_lambda_index,
      sel_d = picks[[k]]$sel_d, n_eff = picks[[k]]$n_eff,
      d_ge_neff = picks[[k]]$sel_d >= picks[[k]]$n_eff, row.names = NULL)))
    list(nets = nets, signed = lapply(picks, `[[`, "signed"),
         sel_lambda_index = lapply(picks, `[[`, "sel_lambda_index"), diag = dg,
         FPR = unname(cm["FPR"]), TPR = unname(cm["TPR"]), PR = unname(cm["PR"]),
         F1 = unname(cm["F1"]), edges = unname(cm["edges"]))
  }
  dEBIC <- deploy("EBIC"); dBIC <- deploy("BIC"); dAIC <- deploy("AIC")

  # oracle ceiling = argmax F1 over the base single-lambda sweep (uses truth -> labelled oracle)
  f1  <- vapply(Om_by_lambda, function(Oml)
    unname(.mgm_confusion(Oml, true_Omega_list, eval_slices)["F1"]), numeric(1))
  oi  <- .mgm_pick(f1, LAM, maximize = TRUE)
  cmO <- if (!is.na(oi)) .mgm_confusion(Om_by_lambda[[oi]], true_Omega_list, eval_slices)
         else c(FPR = NA, TPR = NA, PR = NA, F1 = NA, edges = NA)

  list(
    # ---- deployed landings (the as-deployed table) ----
    refit     = strip(dEBIC),                 # mgm-native relaxed EBIC = the deployed selector
    refit_BIC = strip(dBIC),
    refit_AIC = strip(dAIC),
    oracle    = list(lambda = if (!is.na(oi)) LAM[oi] else NA_real_,
                     FPR = unname(cmO["FPR"]), TPR = unname(cmO["TPR"]),
                     precision = unname(cmO["PR"]), F1 = unname(cmO["F1"]), edges = unname(cmO["edges"])),
    # ---- the de-biased networks at the deployed (EBIC) point ----
    wadj_refit       = dEBIC$nets,
    signed_refit     = dEBIC$signed,
    sel_lambda_index = dEBIC$sel_lambda_index,
    # ---- FULL per-(estpoint, node, lambda) grid: ANY selector recomputable, no re-fit ----
    grid = lapply(seq_len(m), function(k) list(
      n_eff  = grid_by_k[[k]]$n_eff,
      rss    = do.call(rbind, lapply(grid_by_k[[k]]$node, `[[`, "rss")),   # P x n_lambda
      df     = do.call(rbind, lapply(grid_by_k[[k]]$node, `[[`, "df")),    # P x n_lambda
      chol_ok = do.call(rbind, lapply(grid_by_k[[k]]$node, `[[`, "chol_ok")),
      neighbours = lapply(grid_by_k[[k]]$node, `[[`, "Nb"),                # [[p]][[j]] = indices
      coefs      = lapply(grid_by_k[[k]]$node, `[[`, "coefs"))),           # [[p]][[j]] = de-biased coefs
    # ---- rank / over-fit diagnostics ----
    diagnostics = list(n_eff_by_slice = n_eff_by_k, chol_fail = chol_fail,
                       grid_overfit = grid_overfit,
                       deployed = rbind(dEBIC$diag, dBIC$diag, dAIC$diag)),
    gamma = gamma)
}
