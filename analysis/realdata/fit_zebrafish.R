# ---------------------------------------------------------------------------
# analysis/realdata/fit_zebrafish.R
#
# Fit the PUBLICATION TVCGLasso estimator to the zebrafish time course, one fit
# per (P, fish group), and store the same intermediate the simulation stores so
# any selector or figure is recomputable later without re-fitting.
#
#   Rscript analysis/realdata/fit_zebrafish.R            # P = 15 and 25, both groups
#   TVCG_ZEB_P=15 Rscript analysis/realdata/fit_zebrafish.R
#
# WRITES  results/zebrafish/zeb_fit_P<P>_<group>.rds
#
# WHAT MAKES THIS THE PUBLICATION ESTIMATOR AND NOT THE EARLIER EXPLORATORY FIT:
#   free_diag = TRUE      the diagonal is estimated, not frozen at the lambda_max
#                         glasso fit (adopted 2026-07-22; it changes the selected
#                         graph, so an old frozen-diagonal fit is not comparable)
#   50 lambda             the locked path resolution, data-adaptive endpoints
#   refit at EVERY lambda relaxed (unpenalised) re-estimation on the fixed pass-1
#                         support -- the paper's methodological addition
#   deployed = refit-BIC  on the JOINT LNM likelihood at the inferred latent Zhat
#                         (multinomial + Gaussian), NOT the Gaussian-only
#                         likelihood of the raw plug-in ALR Z_0
#
# NO TRUTH EXISTS HERE. There is no ROC, no FPR/TPR, no oracle. The deployed
# point is whatever the method's own selector picks, which is the whole reason
# the selector had to be settled on simulated data first.
#
# UNEQUAL n PER DAY. The zebrafish days carry 7-14 (infected) and 15-23
# (uninfected) fish, so every quantity that pools slices must weight slice k by
# n_k/N rather than 1/m. The objective, the gradient and the information
# criteria in R/refit.R already do (the P0-6 rewrite); see the note at the end
# of this file for the one place that did not.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here); library(splines); library(glasso); library(Matrix)
  library(MASS); library(Rcpp)
})
source(here::here("analysis", "config.R"))
source(here::here("R", "tvcglasso.R"))     # engine + G_beta_Rcpp
source(here::here("R", "refit.R"))         # refit_fixed_support + IC helpers
source(here::here("R", "methods", "method_tv.R"))   # tv_lambda_grid
source(here::here("data", "prepare_zebrafish.R"))   # zeb_slices

`%||%` <- function(a, b) if (is.null(a)) b else a

## --- reproducibility --------------------------------------------------------
## The estimator itself is deterministic: neither R/tvcglasso.R nor the refit
## path draws a random number (the only rnorm() in R/refit.R is inside
## refit_fd_check, a finite-difference gradient checker that nothing calls and
## that seeds itself). The seed is set anyway, for two reasons:
##   1. it makes the claim checkable rather than a matter of trust, and it will
##      keep the run reproducible if any component ever gains a stochastic step;
##   2. mclapply gives each forked child its OWN seed by default, so a future
##      random step would silently make the parallel path disagree with the
##      serial one. RNGkind("L'Ecuyer-CMRG") is the stream generator
##      parallel::mc* uses to make forked streams reproducible.
RNGkind("L'Ecuyer-CMRG")
SEED <- 20260727L
set.seed(SEED)

OUT_DIR <- here::here("results", "zebrafish")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
CKPT    <- file.path(OUT_DIR, "_ckpt"); dir.create(CKPT, showWarnings = FALSE)

## TVCG_ZEB_P / TVCG_ZEB_GROUP restrict the run to one cell, which is how the
## SLURM array splits the four (P, group) units across elements. Unset = all four
## in sequence, which is what a local run does.
P_LIST <- if (nzchar(Sys.getenv("TVCG_ZEB_P"))) as.integer(Sys.getenv("TVCG_ZEB_P")) else c(15L, 25L)
GROUPS <- { g <- Sys.getenv("TVCG_ZEB_GROUP")
            if (nzchar(g)) g else c("infected", "not_infected") }
stopifnot(all(GROUPS %in% c("infected", "not_infected")))

## Publication settings, mirroring config/pub_tvcglasso.R.
CFG <- list(
  ## TVCG_ZEB_NLAMBDA exists only so the pipeline can be smoke-tested in
  ## seconds; the publication path length is 50 and nothing else should change it.
  tv_n_lambda      = as.integer(Sys.getenv("TVCG_ZEB_NLAMBDA", "50")),
  tv_init_mode     = "diag",
  tv_sel_type      = "hard",
  tv_weight_mode   = "glasso",
  tv_free_diag     = TRUE,
  tv_max_iter      = 150L,
  tv_refit_max_outer = 500L,
  refit_inner_max    = 10L,
  refit_z_align_max  = 160L,
  refit_conv_tol_Z      = 5e-5,
  refit_conv_tol_z_grad = 1e-4
)
## m = 7 for both groups -> the m=7 basis of the publication grid.
Q <- 2L; N_KNOT <- 1L

clean <- readRDS(here::here("data", "zebrafish_clean.rds"))

fit_one <- function(P, group) {
  tag  <- sprintf("P%02d_%s", P, group)
  fout <- file.path(OUT_DIR, sprintf("zeb_fit_%s.rds", tag))
  if (file.exists(fout)) { message("skip (exists): ", basename(fout)); return(invisible(NULL)) }
  t0 <- Sys.time()

  ## zeb_slices supplies everything the estimator needs, INCLUDING the time
  ## axis. Do not recompute it here: `t` is the real sampling days rescaled to
  ## [0,1] with the irregular spacing preserved (days 7,10,21,30,43,59,86 are
  ## not equally spaced, and pretending they are would misplace the knots).
  sl <- zeb_slices(clean, P = P, group = group, reference = "top_prevalence")
  X  <- sl$X; m <- length(X)
  days <- sl$days; xseq <- as.numeric(sl$t); n_per_slice <- sl$n_per_slice
  message(sprintf("\n=== %s : P=%d  m=%d  n per day = %s  (N=%d) ===",
                  tag, P, m, paste(n_per_slice, collapse = ","), sum(n_per_slice)))

  mid   <- seq(xseq[1], xseq[m], length.out = N_KNOT + 2)[-c(1, N_KNOT + 2)]
  basis <- bs(xseq, degree = Q, knots = mid,
              Boundary.knots = c(xseq[1], xseq[m]), intercept = TRUE)

  ## lambda grid: the SAME data-adaptive rule the simulation uses
  ## (tv_lambda_grid: lambda_max = 2 * max over slices of max |off-diagonal
  ## cov(ALR Z_0)|, down to lambda_max/1000, log-spaced). It is part of the
  ## method, so the real-data analysis uses it unchanged -- inventing a
  ## different anchoring rule here would make the application inconsistent
  ## with the simulation and would be an undisclosed change to the estimator.
  lam_grid <- tv_lambda_grid(X, CFG$tv_n_lambda)
  message(sprintf("  lambda grid [%.5f, %.4f]", min(lam_grid), max(lam_grid)))

  ## ---- pass 1 ----------------------------------------------------------
  fits <- tv_warm_path(X, lam_grid, xseq, q = Q, N_n = N_KNOT,
                       init_mode = CFG$tv_init_mode, sel_type = CFG$tv_sel_type,
                       weight_mode = CFG$tv_weight_mode, free_diag = CFG$tv_free_diag,
                       Max_iterations = CFG$tv_max_iter,
                       ckpt_file = file.path(CKPT, sprintf("base_%s.rds", tag)))
  lam <- as.numeric(names(fits))

  Om_by_lambda <- lapply(fits, function(f) G_beta_Rcpp(f$beta, basis, m))
  n_edges <- vapply(Om_by_lambda, function(Ol)
    sum(vapply(Ol, function(O) sum(O[upper.tri(O)] != 0), integer(1))), integer(1))
  df_by_lambda <- vapply(fits, function(f)
    sum(vapply(f$beta, function(b) sum(b[upper.tri(b)] != 0), integer(1))), integer(1))

  ## ---- pass 2: refit at every lambda, IN PARALLEL -----------------------
  ## The refits at different lambda are INDEPENDENT: each starts from its own
  ## pass-1 fit (there is no cross-lambda warm start), uses no RNG, and writes
  ## its own checkpoint. Evaluating them concurrently therefore returns the
  ## same list, in the same order, with the same values -- this mirrors the
  ## TVCG_REFIT_CORES path in R/methods/method_tv.R, where the equivalence is
  ## documented. Pass 1 stays serial because it IS a continuation: each lambda
  ## warm-starts from the previous one.
  X_work <- refit_prepare_counts(X)
  n_core <- max(1L, min(as.integer(Sys.getenv("TVCG_REFIT_CORES",
                                              as.character(parallel::detectCores() - 1L))),
                        length(fits)))
  message(sprintf("  refit: %d lambda on %d core(s)", length(fits), n_core))
  .refit_one <- function(i) {
    f <- fits[[i]]
    r <- tryCatch(
      refit_fixed_support(X_work = X_work, Z_start = f$Z, beta_start = f$beta, basis = basis,
                          checkpoint_file = file.path(CKPT, sprintf("refit_%s_lam%02d.rds", tag, i)),
                          progress_log    = file.path(CKPT, sprintf("refit_%s_lam%02d.log", tag, i)),
                          max_outer = CFG$tv_refit_max_outer, inner_max = CFG$refit_inner_max,
                          z_align_max = CFG$refit_z_align_max,
                          conv_tol_Z = CFG$refit_conv_tol_Z,
                          conv_tol_z_grad = CFG$refit_conv_tol_z_grad),
      error = function(e) structure(list(msg = conditionMessage(e)), class = "refit_error"))
    if (inherits(r, "refit_error"))
      return(list(lambda = lam[i], lambda_index = i, df = NA_integer_, error = r$msg))
    ic <- function(cr) {
      v <- refit_information_criteria(cr$total, n_per_slice, r$df, P)
      list(mult = cr$multinomial, gauss = cr$neg_logdet + cr$trace, nll = cr$total,
           neg2loglik = unname(v["neg2loglik"]), AIC = unname(v["AIC"]),
           BIC = unname(v["BIC"]), eBIC = unname(v["eBIC"]))
    }
    list(lambda = lam[i], lambda_index = i, df = r$df, support = r$support,
         pre = ic(r$pre$criterion), post = ic(r$post$criterion),
         beta_post = r$post$beta, Omega_post = r$post$criterion$Omega, Z_post = r$post$Z,
         converged = isTRUE(r$converged), exit_reason = r$exit_reason,
         z_align_converged = isTRUE(r$z_align_converged),
         max_active_grad = r$max_active_grad, max_z_grad = r$max_z_grad)
  }
  refit <- if (n_core > 1L)
    parallel::mclapply(seq_along(fits), .refit_one, mc.cores = n_core, mc.preschedule = FALSE)
  else lapply(seq_along(fits), .refit_one)
  ## A forked worker that dies outside R's condition system (OOM/segfault) returns a
  ## try-error rather than a list; normalise it to the SAME error element the serial
  ## tryCatch produces so everything downstream sees one shape.
  refit <- lapply(seq_along(refit), function(i) {
    x <- refit[[i]]
    if (is.list(x) && !is.null(x$lambda_index)) return(x)
    list(lambda = lam[i], lambda_index = i, df = NA_integer_,
         error = if (inherits(x, "try-error")) as.character(x) else "parallel refit worker failed")
  })

  ic_col <- function(pass, key) vapply(refit, function(x) {
    p <- x[[pass]]; if (is.null(p) || is.null(p[[key]])) NA_real_ else p[[key]] }, numeric(1))
  detail <- data.frame(
    lambda = lam, lambda_index = seq_along(lam),
    df = df_by_lambda, n_edges = n_edges,
    converged = vapply(refit, function(x) isTRUE(x$converged), logical(1)),
    AIC = ic_col("post","AIC"), BIC = ic_col("post","BIC"), eBIC = ic_col("post","eBIC"),
    AIC_pre = ic_col("pre","AIC"), BIC_pre = ic_col("pre","BIC"), eBIC_pre = ic_col("pre","eBIC"),
    row.names = NULL)

  ## ---- deployed point: refit-BIC argmin over CONVERGED lambdas ----------
  bic_conv <- ifelse(detail$converged, detail$BIC, NA_real_)
  di <- if (any(is.finite(bic_conv))) which.min(bic_conv) else which.min(detail$BIC)
  deployed <- list(selector = "refit-BIC", lambda_index = di, lambda = lam[di],
                   df = detail$df[di], n_edges = detail$n_edges[di],
                   refit_converged = isTRUE(refit[[di]]$converged))

  res <- list(
    method = "tvcglasso", data = "zebrafish", P = P, group = group,
    ## NOTE the node ORDER here is prepare_zebrafish.R's explicit rule
    ## (prevalence, then mean relative abundance, then alphabetical). It is the
    ## same SET as earlier exploratory fits but not the same order, so edge
    ## labels must be read from `node_names` and never from a stored index.
    node_names = sl$nodes, ref_name = sl$reference,
    days = days, x_sequence = xseq, basis = basis,
    n_per_slice = n_per_slice, N = sum(n_per_slice), m = m, q = Q, N_n = N_KNOT,
    detail = detail, deployed = deployed,
    estimate = list(base = list(beta = lapply(fits, `[[`, "beta"),
                                Omega = Om_by_lambda,
                                Z = lapply(fits, `[[`, "Z")),
                    refit = refit),
    provenance = list(config = CFG, git = Sys.getenv("TVCG_GIT_HASH", unset = NA),
                      R = R.version.string, date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                      data_provenance = clean$provenance),
    secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  saveRDS(res, fout)
  message(sprintf("  -> %s   deployed lambda[%d]=%.4f  df=%d  edges=%d  refit_converged=%s  (%.1f min)",
                  basename(fout), di, lam[di], detail$df[di], detail$n_edges[di],
                  deployed$refit_converged, res$secs / 60))
  invisible(res)
}

for (P in P_LIST) for (g in GROUPS) fit_one(P, g)
message("\nall zebrafish fits done -> ", OUT_DIR)
