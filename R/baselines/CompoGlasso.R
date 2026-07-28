# ============================================================================
# R/baselines/CompoGlasso.R — VENDORED third-party solver. NOT OUR CODE.
#
# ORIGIN
#   Compositional Graphical Lasso ("CGLasso"), the solver accompanying
#     Tian, Jiang, Hammer, Sharpton & Jiang (2023), "Compositional Graphical
#     Lasso Resolves the Impact of Parasitic Infection on Gut Microbial
#     Interaction Networks in a Zebrafish Model", JASA 118(544):1500-1514.
#     doi:10.1080/01621459.2022.2164287
#   Upstream file: https://github.com/yuanjiang-osu/Comp-gLASSO-JASA
#                  CompoGlasso.R (278 lines)
#   An earlier release of the same codebase is published by the same author
#   under the MIT License: https://github.com/yuanjiang-osu/Comp-gLASSO
#   Copyright of this file remains with its original authors.
#
# WHY VENDORED
#   So the baseline depends on a fixed, version-controlled solver rather than on
#   a large external reference tree that is not part of this repository.
#
# EXACTLY WHAT WE CHANGED — the complete list, verified by diff against upstream:
#   1. RENAMED three functions to avoid a global-namespace collision:
#        obj -> cg_obj, NR -> cg_NR, NR_para -> cg_NR_para.
#      Our TV engine defines its own obj()/NR(); without the prefix, sourcing
#      this file would silently redefine them and change the TV estimator.
#      `z_hat_offset` and `Compo_glasso` keep their upstream names.
#   2. ADDED `nr_max_iter = 100L`, a termination guard on the per-sample Newton
#      while-loop in cg_NR / cg_NR_para. It is the ONLY semantic change in this
#      file, it is inert whenever the loop converges normally, and it is
#      documented in full at the head of cg_NR.
#   3. OMITTED `generate_cov` (the upstream data generator); this project uses
#      its own generators in simulation/generators.R.
#   Upstream inline comments are otherwise preserved. One commented-out line
#   ("# Omegas.1[, , i.rho.1] <- Omega.1") was dropped.
#   No other function body, default, threshold or convergence rule was touched.
#
# WHAT THIS FILE PROVIDES
#   Compo_glasso(x, rho.list, ...) is a STATIC single-network LNM solver: it takes
#   ONE time-slice count matrix x (n x (K+1), reference taxon = last column) and
#   returns the precision path over rho.list. The multi-timepoint "static control"
#   wrapper is ours and lives in R/baselines/cglasso_static.R.
# ============================================================================

# Function for additive log-ratio transformation
# Option give different options to handle zero's in the numerator/denominator in the ratio
# Only option = 0 or 2 is used in our work
z_hat_offset <- function(x, offset, option)
{
  K <- dim(x)[2] - 1
  M <- apply(x, 1, sum)
  N <- dim(x)[1]

  # no offset
  if(option == 0){
    z.hat <- log(x[, -(K + 1)]/x[, (K + 1)])
  }
  # x offset by a constant
  if(option == 1){
    p.hat <- rep(1/(K + 1), K + 1)
    x.adj <- t(t(x) + p.hat * offset)
    z.hat <- log(x.adj[, -(K + 1)]/x.adj[, (K + 1)])
  }
  # x offset proportionally
  else if(option == 2){
    p.hat <- colMeans(x/M)
    x.adj <- t(t(x) + p.hat * offset)
    z.hat <- log(x.adj[, -(K + 1)]/x.adj[, (K + 1)])
  }
  # ratio offset proportionally only for zero x_j or x_{K+1}
  else if(option == 3){
    p.hat <- colMeans(x/M)
    z.hat <- matrix(NA, nrow = nrow(x), ncol = K)
    for(j in 1 : K){
      zero.ind <- (x[, j] == 0) | (x[, K + 1] == 0)
      z.hat[!zero.ind, j] <- log(x[!zero.ind, j]/x[!zero.ind, (K + 1)])
      z.hat[zero.ind, j] <- log(p.hat[j]/p.hat[K + 1])
    }
  }
  # ratio offset proportionally for all x's
  else if(option == 4){
    p.hat <- colMeans(x/M)
    z.hat <- matrix(NA, nrow = nrow(x), ncol = K)
    for(j in 1 : K){
      zero.ind <- (x[, j] == 0) | (x[, K + 1] == 0)
      z.hat[!zero.ind, j] <- log(x[!zero.ind, j]/x[!zero.ind, (K + 1)] + p.hat[j]/p.hat[K + 1])
      z.hat[zero.ind, j] <- log(p.hat[j]/p.hat[K + 1])
    }
  }
  # ratio offset proportionally only for zero x_j or x_{K+1}
  else if(option == 5){
    p.hat <- matrix(0, N, K + 1)
    z.hat <- matrix(NA, nrow = nrow(x), ncol = K)
    for (i in 1:N) {
        zeros <- which(x[i, ] == 0)
        nzeros <- which(x[i, ] != 0)
        p.hat[i, zeros] <- (x[i, zeros] + offset)/(M[i] + offset * length(zeros))
        p.hat[i, nzeros] <- (x[i, nzeros])/(M[i] + offset * length(zeros))
        z.hat[i, ] <- log(p.hat[i, -(K + 1)]/p.hat[i, K + 1])
      }
  }
  return(z.hat)
}

# Function to calculate the objective function in equation (7) of the paper
# ★ P0-1 (2026-07-22): renamed obj/NR/NR_para -> cg_obj/cg_NR/cg_NR_para so this vendored
# CGLasso solver does NOT overwrite the TV engine's global obj()/NR() when both are sourced
# (the TV estimator must be identical regardless of which baselines are loaded).
cg_obj <- function(x, z, Omega, K) {
  M <- sum(x)
  mu = mean(z)
  f = 1 / 2 * t(z - mu) %*% Omega %*% (z - mu) - (t(x) %*% z - M * log(as.numeric(t(rep(1, K)) %*% exp(z) + 1)))
  return(as.numeric(f))
}

# Newton-Raphson procedure to optimize the objective function in equation (7) of the paper (with parallelization)
 # Same SAFETY CAP as cg_NR (see its header): bounds ONLY the per-sample Newton while-loop, so a
 # diverged Omega cannot spin it forever. Not used by our pipeline (we call with para_NR = FALSE),
 # capped for consistency.
 cg_NR_para <- function(x, z.0, Omega.0, alpha_0 = 1, delta = 5 , epsilon = 0.01, threshold = 0.0001, num_cores = 1,
                        nr_max_iter = 100L)
{
  # Initialization
  n = dim(z.0)[1]
  K = dim(z.0)[2]
  M <- as.numeric(apply(x, 1, sum))
  mu.0 = apply(z.0, 2, mean)
  z.1 <- matrix(0, n, K)

  # Parallelization begins
  library(foreach)
  library(doParallel)
  cl<-makeCluster(num_cores)
  registerDoParallel(cl)
  z.new <- foreach (j = 1:n, .combine = rbind, .export = "cg_obj") %dopar% {
    z.iter <- 0
    alpha <- alpha_0
    # Loop to update z
    while (mean((z.0[j,] - z.1[j,]) ^ 2) > threshold && z.iter < nr_max_iter) {
      if (z.iter != 0 && (cg_obj(x[j, 1:K], z.1[j,], Omega.0, K) <= (cg_obj(x[j, 1:K], z.0[j,], Omega.0, K) + epsilon * alpha * h_0))) {
        z.0[j,] <- z.1[j,]
      }
      z.iter <- z.iter + 1
      # Gradient of the objective function in equation (7) with respect to z
      dipi <- M[j] * exp(z.0[j,]) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1) - x[j, 1:K] +  as.vector(Omega.0 %*% (z.0[j,] - mu.0))
      # Hessian matrix of the objective function in equation (7) with respect to z
      tripi <- M[j] * diag(exp(z.0[j,])) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1) -
        M[j] * (exp(z.0[j,])) %*% t(exp(z.0[j,])) / (as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1)) ^ 2 + Omega.0
      # Newton-Raphson's updating rule
      z.1[j,] <- z.0[j,] - alpha * solve(tripi) %*% dipi
      # Shrink the step size using Armijo's Rule
      dk = (-1) * solve(tripi) %*% dipi
      h_0 = t(z.0[j,] - mu.0) %*% Omega.0 %*% dk - x[j, 1:K] %*% dk +
        M[j] * as.numeric(t(dk) %*% exp(z.1[j,])) / as.numeric(t(rep(1, K)) %*% exp(z.1[j,]) + 1)
      if (cg_obj(x[j, 1:K], z.1[j,], Omega.0, K) > (cg_obj(x[j, 1:K], z.0[j,], Omega.0, K) + epsilon * alpha * h_0)) {
        alpha <- alpha / delta
      }
    }
    z.1[j,]
  }
  stopCluster(cl)

  return(unname(z.new, force = TRUE))
}

# Newton-Raphson procedure to optimize the objective function in equation (7) of the paper
 # ★ SAFETY CAP ONLY (2026-07-25, added by us; NOT a change to the author's method).
 # The author's per-sample Newton loop `while (mean((z.0-z.1)^2) > threshold)` has NO iteration
 # bound. In the author's own setting (one static network, well-conditioned Omega) it converges in
 # a few steps, so no bound is needed. In OUR relaxed-refit reuse it can be fed a DIVERGED Omega
 # (ill-conditioned dense supports at n<P make the Z<->Omega alternation blow up: observed
 # dOmega ~ 5e9); then exp(z) overflows, the while condition can never be satisfied, and the loop
 # spins forever (measured: 20+ min at 99.7% CPU, no output, the job never returns).
 # `nr_max_iter` bounds THAT LOOP AND NOTHING ELSE: the update rule, the Armijo line search, the
 # threshold and the returned value are byte-identical to the author's code whenever the loop
 # converges normally (i.e. every well-conditioned case). A rho that hits the cap returns its
 # current iterate; the refit then flags it non-converged and its garbage likelihood makes its BIC
 # huge, so it is never the selected operating point.
 cg_NR <- function(x, z.0, Omega.0, alpha_0 = 1, delta = 5 , epsilon = 0.01, threshold = 0.0001,
                   nr_max_iter = 100L)
{
  # Initialization
  n = dim(z.0)[1]
  K = dim(z.0)[2]
  M <- as.numeric(apply(x, 1, sum))
  mu.0 = apply(z.0, 2, mean)
  z.1 <- matrix(0, n, K)

  for (j in 1:n)
  {
    z.iter <- 0
    alpha <- alpha_0
    # Loop to update z
    while (mean((z.0[j,] - z.1[j,]) ^ 2) > threshold && z.iter < nr_max_iter) {
      if (z.iter != 0 && (cg_obj(x[j, 1:K], z.1[j,], Omega.0, K) <= (cg_obj(x[j, 1:K], z.0[j,], Omega.0, K) + epsilon * alpha * h_0))) {
        z.0[j,] <- z.1[j,]
      }
      z.iter <- z.iter + 1
      # Gradient of the objective function in equation (7) with respect to z
      dipi <- M[j] * exp(z.0[j,]) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1) - x[j, 1:K] +  as.vector(Omega.0 %*% (z.0[j,] - mu.0))
      # Hessian matrix of the objective function in equation (7) with respect to z
      tripi <- M[j] * diag(exp(z.0[j,])) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1) -
        M[j] * (exp(z.0[j,])) %*% t(exp(z.0[j,])) / (as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1)) ^ 2 + Omega.0
      # Newton-Raphson's updating rule
      z.1[j,] <- z.0[j,] - alpha * solve(tripi) %*% dipi
      # Shrink the step size using Armijo's Rule
      dk = (-1) * solve(tripi) %*% dipi
      h_0 = t(z.0[j,] - mu.0) %*% Omega.0 %*% dk - x[j, 1:K] %*% dk +
        M[j] * as.numeric(t(dk) %*% exp(z.0[j,])) / as.numeric(t(rep(1, K)) %*% exp(z.0[j,]) + 1)
      if (cg_obj(x[j, 1:K], z.1[j,], Omega.0, K) > (cg_obj(x[j, 1:K], z.0[j,], Omega.0, K) + epsilon * alpha * h_0)) {
        alpha <- alpha / delta
      }
    }
  }
  return(z.1)
}

# Main function for compositional graphical lasso
Compo_glasso <- function(x, rho.list, option = 2, offset = (K + 1), para_NR = FALSE, num_cores_NR = 1, z_ratio = 1000, O_ratio = 1000, max_iter = 50)
{
  gc()
  library(MASS)
  library(glasso)
  library(huge)
  library(propagate)

  n <- dim(x)[1] # number of samples
  K <- dim(x)[2] - 1 # number of OTUs (reference OTU not counted)
  M <- as.numeric(apply(x, 1, sum)) # Sequencing depths

  results <- list() # List to save results

  # Additive log-ratio transformation
  z.hat <- z_hat_offset(x, offset, option)
  colnames(z.hat) = NULL
  rownames(z.hat) = NULL
  cat("z_hat comupted \n")

  rho.list.1  <- rho.list
  n1.rho <- length(rho.list.1)
  Sigma.2 <- bigcor(z.hat, fun = "cov", verbose = FALSE)
  Sigma.2 <- Sigma.2[1:nrow(Sigma.2), 1:ncol(Sigma.2)]

  # Run graphical lasso on the additive log-ratio transformed data
  path.2 <- huge(Sigma.2, method = "glasso", lambda = max(rho.list))
  Omegas.2 <- path.2$icov
  Sigma.0 <- Sigma.2
  Omegas.1 <- array(0, dim = c(K, K, n1.rho))

  for(i.rho.1 in 1 : n1.rho)
  {
    if (n1.rho != 1) cat("i.rho.1 =", i.rho.1, "\n")
    rho <- rho.list.1[i.rho.1]

    # Initialization
    if (i.rho.1 == 1) {
      Omega.0 <- as.matrix(Omegas.2[[1]])
    }
    else if (i.rho.1 > 1) {
      Omega.0 <- Omegas.1[, , (i.rho.1 - 1)]
    }

    # Iteration between Newton-Raphson and graphical lasso
    Omega.1 <- matrix(0, K, K)
    iter <- 0
    z.start <- z.hat
    z.end <- matrix(0, n, K)
    O_thr <- mean((Omega.0 - Omega.1) ^ 2) / O_ratio # convergence threshold
    z_thr <- mean((z.start - z.end) ^ 2) / z_ratio # convergence threshold
    while ((mean((Omega.0 - Omega.1) ^ 2) > O_thr || (mean((z.start - z.end) ^ 2) > z_thr)) && iter <= max_iter)
    {
      cat("iter = ", iter + 1, "mean((Omega.0 - Omega.1) ^ 2) = ", mean((Omega.0 - Omega.1) ^ 2), "\n")
      if (iter != 0) {
        Omega.0 <- Omega.1
        z.start <- z.end
      }
      iter <- iter + 1

      # Update z: Newton-Raphson
      if (para_NR == TRUE) z.end <- cg_NR_para(x, z.start, Omega.0, num_cores = num_cores_NR) # with parallelization
      else z.end <- cg_NR(x, z.start, Omega.0) # without parallelization
      cat("mean square z.end - z.start =", mean((z.start - z.end) ^ 2), "\n")

      # Update Omega: Graphical Lasso
      Sigma.1 <- bigcor(z.end, fun = "cov", verbose = FALSE)
      Sigma.1 <- Sigma.1[1:nrow(Sigma.1), 1:ncol(Sigma.1)]
      mod <- huge(x = Sigma.1, lambda = rho, method = "glasso")
      Omega.1 = as.matrix(mod$icov[[1]])

      # Added 04/08/2022 by Yuan Jiang to avoid ERROR
      if (is.na(mean((Omega.0 - Omega.1) ^ 2)) || is.na(mean((z.start - z.end) ^ 2)) )
        break
    }

    # Updated 04/08/2022 by Yuan Jiang to avoid ERROR
    if (any(is.na(Omega.1))) {
      Omegas.1[, , i.rho.1] <- Omega.0
    } else {
      Omegas.1[, , i.rho.1] <- Omega.1
    }

  }

  results <- list(Omega.1, Omegas.1, rho.list.1)
  names(results) <- c("Omega.1", "Omegas.1", "rho.list.1")
  return(results)
}
