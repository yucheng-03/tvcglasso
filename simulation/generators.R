# ============================================================================
# simulation/generators.R — publication simulation-data generators.
#
# TWO generators only (LUPINE / oracle / the discontinuous-0617 mirror-exp are
# deliberately NOT here):
#   generate_simulation_data_mirrorexp_contB  — the MAIN-TEXT "exp f" (contB):
#       continuous mirror-exp, each active edge decays to EXACT 0 at t=0.5 (kink,
#       no jump) => a genuine temporal zero region.  (Zero-point-ID generator.)
#   generate_hetrate_contBfast                — the CONTINUOUS heterogeneous-rate
#       option: 2 FAST edges use the contB truncated-to-zero shape, 2 SLOW edges
#       (const / gentle exp) stay always-on.  (Heterogeneity generator.)
#
# generate_data(cell, seed) is a STRICT dispatcher: an unknown cell$generator is a
# hard error (NO silent fallback to any default generator — this is what caused the
# batch8 "discontinuous 0617" mix-up). Both generators are copied verbatim from the
# project; the ONLY change is the real_depths_path default -> data/.
#
# Requires: MASS, Matrix, here.  Run from the repo root so here::here() anchors.
# ============================================================================
suppressPackageStartupMessages({ library(MASS); library(Matrix) })

`%||%` <- function(a, b) if (is.null(a)) b else a


# --- MAIN-TEXT generator: continuous-B mirror-exp ("exp f") --------------------
# f_early(t) = es*(exp(-rt)-c)/(1-c) on [0,0.5], else 0 ; c = exp(-rate/2); f_late = mirror.
# Each active edge is EXACTLY 0 on its off-half (clean binary zero-region truth).
# Z_0 = +0.5 pseudocount ALR init (unchanged from the project).
generate_simulation_data_mirrorexp_contB <- function(
    n = 20, P = 40, m = 7,
    structure = "random",
    seed = 1,
    edge_strength = 0.5,
    rate = 6,
    random_edge_prob = NULL,                 # default 3/P (Yuan JASA)
    depth_mode = c("low", "high", "real"),
    depth_lo = NULL, depth_hi = NULL,          # custom uniform band U(depth_lo*P, depth_hi*P); overrides depth_mode when both given
    real_depths_path = here::here("data", "zebrafish_real_depths.rds"),
    base_diag_mean = 0.1
) {
  depth_mode <- match.arg(depth_mode)
  if (is.null(random_edge_prob)) random_edge_prob <- 3 / P
  set.seed(seed)
  parameters <- as.list(environment())
  x_sequence <- seq(0, 1, length.out = m)          # EQUAL spacing

  # --- CONTINUOUS-B mirror-exp time functions (continuous at t=0.5, half-zero support) ---
  c0 <- exp(-rate / 2); Zn <- 1 - c0
  f_early <- function(t) ifelse(t <= 0.5, edge_strength * (exp(-rate * t)       - c0) / Zn, 0)
  f_late  <- function(t) ifelse(t >= 0.5, edge_strength * (exp(-rate * (1 - t)) - c0) / Zn, 0)

  # --- masks: each pair -> early / late with prob random_edge_prob, ~50/50 ---
  Mask_early <- matrix(0, P, P); Mask_late <- matrix(0, P, P)
  for (i in 1:(P - 1)) for (j in (i + 1):P) {
    if (runif(1) < random_edge_prob) {
      if (runif(1) < 0.5) { Mask_early[i, j] <- Mask_early[j, i] <- 1 }
      else                { Mask_late[i, j]  <- Mask_late[j, i]  <- 1 }
    }
  }

  # --- true Omega(t): diagonal-dominant (PD) ---
  true_Omega_list <- lapply(x_sequence, function(t) {
    Off <- Mask_early * f_early(t) + Mask_late * f_late(t)
    Omega <- Off; diag(Omega) <- rowSums(abs(Off)) + base_diag_mean; Omega
  })

  # --- depths: low=U(20P,40P), high=U(500P,600P), real=sampled Zebrafish ---
  draw_depths <- function(n) {
    if (!is.null(depth_lo) && !is.null(depth_hi)) runif(n, depth_lo * P, depth_hi * P)
    else if (depth_mode == "low")  runif(n, 20 * P, 40 * P)
    else if (depth_mode == "high") runif(n, 500 * P, 600 * P)
    else { rd <- readRDS(real_depths_path); sample(rd, n, replace = TRUE) }
  }

  # --- Z ~ N(0, Omega^{-1}) per slice -> LNM multinomial counts ---
  mu_z <- rep(0, P)
  Z_t <- lapply(true_Omega_list, function(omega) {
    omega <- (omega + t(omega)) / 2
    sigma <- tryCatch(chol2inv(chol(omega)), error = function(e) solve(as.matrix(nearPD(omega)$mat)))
    mvrnorm(n = n, mu = mu_z, Sigma = sigma)
  })
  p_total <- lapply(Z_t, function(Z) { eZ <- exp(Z); d <- rowSums(eZ) + 1; cbind(eZ / d, 1 / d) })
  X_total <- lapply(p_total, function(pm) {
    Md <- draw_depths(nrow(pm))
    t(sapply(1:nrow(pm), function(i) rmultinom(1, size = round(Md[i]), prob = pm[i, ])))
  })

  # Z_0 (ALR init) — +0.5 pseudocount
  Z_0 <- lapply(X_total, function(Xc) { Xc <- Xc + 0.5; ph <- Xc / rowSums(Xc); log(ph[, -(P + 1)] / ph[, P + 1]) })

  list(X = X_total, Z_0 = Z_0, true_Omega_list = true_Omega_list,
       x_sequence = x_sequence, n_edges_early = sum(Mask_early) / 2,
       n_edges_late = sum(Mask_late) / 2, parameters = parameters)
}


# --- CONTINUOUS heterogeneous-rate option: hetrate + contB fast shapes ---------
# 2 FAST edges (fast_e/fast_l) use the contB truncated-to-EXACT-0 shape; 2 SLOW
# edges (const, gentle exp) stay always-on. Rate heterogeneity WITH a real zero
# region on the fast edges. Z_0 = +1 pseudocount, CENTERED (the hetrate convention).
.contBfast_traj <- function(typ, xs, A, fast_rate, slow_rate) {
  c0 <- exp(-fast_rate / 2)                                   # contB continuity constant (zero at t=0.5)
  switch(typ,
    fast_e = ifelse(xs <= 0.5, A * (exp(-fast_rate * xs)       - c0) / (1 - c0), 0),   # sharp decay -> EXACT 0 on (0.5,1]
    fast_l = ifelse(xs >= 0.5, A * (exp(-fast_rate * (1 - xs)) - c0) / (1 - c0), 0),   # EXACT 0 on [0,0.5) -> sharp rise
    const  = rep(A, length(xs)),                               # slow group: constant (always on)
    slow   = A * exp(-slow_rate * xs),                         # slow group: gentle decay (always on, never 0)
    stop("unknown contBfast type: ", typ))
}

generate_hetrate_contBfast <- function(n, P, m, seed, A = 0.5, edge_prob = 3 / P,
                                       fast_rate = 6, slow_rate = 1, frac_fast = 0.5,
                                       depth_mode = "real",
                                       real_depths_path = here::here("data", "zebrafish_real_depths.rds")) {
  stopifnot(requireNamespace("MASS", quietly = TRUE), requireNamespace("Matrix", quietly = TRUE))
  set.seed(seed)
  xs <- seq(0, 1, length.out = m)
  edges <- list()
  for (i in 1:(P - 1)) for (j in (i + 1):P) if (runif(1) < edge_prob) {
    typ <- if (runif(1) < frac_fast) sample(c("fast_e", "fast_l"), 1) else sample(c("const", "slow"), 1)
    edges[[length(edges) + 1]] <- list(i = i, j = j,
      vals = .contBfast_traj(typ, xs, A, fast_rate, slow_rate), typ = typ)
  }
  tol <- lapply(seq_len(m), function(k) {
    Off <- matrix(0, P, P)
    for (e in edges) Off[e$i, e$j] <- Off[e$j, e$i] <- e$vals[k]
    Om <- Off; diag(Om) <- rowSums(abs(Off)) + 0.1; Om
  })
  depths <- if (depth_mode == "real" && file.exists(real_depths_path)) readRDS(real_depths_path)
            else if (depth_mode == "high") round(runif(max(n, 100), 500 * P, 600 * P))   # high = U(500P,600P) (match contB)
            else round(runif(max(n, 100), 20 * P, 40 * P))                                # low  = U(20P,40P)
  Z_t <- lapply(tol, function(om) {
    om <- (om + t(om)) / 2
    sg <- tryCatch(chol2inv(chol(om)), error = function(e) solve(as.matrix(Matrix::nearPD(om)$mat)))
    MASS::mvrnorm(n, rep(0, P), sg)
  })
  X <- lapply(Z_t, function(Z) {
    eZ <- exp(Z); d <- rowSums(eZ) + 1; pm <- cbind(eZ / d, 1 / d)
    Md <- sample(depths, nrow(pm), TRUE)
    t(sapply(1:nrow(pm), function(i) rmultinom(1, round(Md[i]), pm[i, ])))
  })
  Z_0 <- lapply(X, function(Xi) { Xo <- Xi + 1; z <- log(Xo[, -ncol(Xo), drop = FALSE] / Xo[, ncol(Xo)]); scale(z, center = TRUE, scale = FALSE) })
  list(X = X, true_Omega_list = tol, x_sequence = xs, Z_0 = Z_0,
       edge_types = sapply(edges, function(e) e$typ), edge_ij = t(sapply(edges, function(e) c(e$i, e$j))),
       n_edge_types = table(sapply(edges, function(e) e$typ)), n_edges = length(edges))
}


# --- STRICT dispatcher --------------------------------------------------------
# Maps cell$generator -> the matching generator. An unknown/absent generator is a
# HARD ERROR: there is NO silent fallback (the batch8 "0617 mirror-exp" mix-up came
# from an else-branch default). Allowed: "mirrorexp_contB", "hetrate_contBfast".
generate_data <- function(cell, seed) {
  g <- cell$generator
  if (is.null(g) || !nzchar(g))
    stop("cell$generator is required — no silent default. Use 'mirrorexp_contB' or 'hetrate_contBfast'.")
  rdp <- cell$real_depths_path %||% here::here("data", "zebrafish_real_depths.rds")
  if (identical(g, "mirrorexp_contB")) {
    generate_simulation_data_mirrorexp_contB(
      n = cell$n, P = cell$P, m = cell$m, seed = seed,
      structure        = cell$structure        %||% "random",
      edge_strength    = cell$edge_strength    %||% 0.5,
      rate             = cell$rate             %||% 6,
      random_edge_prob = cell$random_edge_prob,          # NULL => 3/P
      depth_mode       = cell$depth_mode       %||% "real",
      depth_lo = cell$depth_lo, depth_hi = cell$depth_hi,
      real_depths_path = rdp,
      base_diag_mean   = cell$base_diag_mean   %||% 0.1)
  } else if (identical(g, "hetrate_contBfast")) {
    generate_hetrate_contBfast(
      n = cell$n, P = cell$P, m = cell$m, seed = seed,
      A          = cell$edge_strength %||% 0.5,
      edge_prob  = cell$edge_prob     %||% (3 / cell$P),
      fast_rate  = cell$fast_rate     %||% 6,
      slow_rate  = cell$slow_rate     %||% 1,
      frac_fast  = cell$frac_fast     %||% 0.5,
      depth_mode = cell$depth_mode    %||% "real",
      real_depths_path = rdp)
  } else {
    stop("unknown generator: '", g, "'. Allowed: 'mirrorexp_contB', 'hetrate_contBfast'. ",
         "(The discontinuous 0617 mirror-exp is deliberately unavailable.)")
  }
}
