# ============================================================================
# config/paired_grid.R — the SHARED paired-data spec (anti-drift guarantee).
#
# BOTH the tvcglasso window AND the CGLasso/JGL/mgm windows MUST source this file so
# every method runs on BIT-IDENTICAL data per (cell,seed). generate_data() is
#   (a) deterministic in seed (set.seed(seed) is its first line), and
#   (b) independent of `methods` (data is built ONCE per (cell,seed), before the
#       method loop, and handed to every method identically),
# so identical cells + seeds => identical data. VERIFY post-hoc via the
# `data_fingerprint` (= sum of all counts) stored in every result's provenance:
# equal fingerprint per (cell,seed) across windows == provably the same data.
#
# Change this file in ONE place and every method's run inherits it.
# ============================================================================

# SUBMISSION = 100 seeds (JASA-grade). A run may take a PREFIX (e.g. PAIRED_SEEDS[1:30]
# for the draft); seeds are additive/paired, so topping 30 -> 100 re-runs nothing.
PAIRED_SEEDS <- 1:100

# 16 cells: P in {15,25} x n in {12,20} x depth in {low,high} x m in {7,15}, contB
# (continuous exp f, NOT the 0617 jump), edge_strength 0.5, rate 6. P=15 is the small-P
# publication cell (matches the real-Zebrafish P=15 analysis, per the P-DESIGN-LOCK).
# Per-cell B-spline basis scales with m (discover convention): m=7 -> q=2/N_n=1;
# m=15 -> q=3 (cubic)/N_n=3. tv_q/tv_N_n travel WITH the cell.
#
# DEPTH bands are UNIFORM U(lo*P, hi*P) (Yuan/Tian convention). NB (recorded 2026-07-22):
# the REAL Zebrafish depth distribution is NOT uniform (bell-shaped, median ~22.5k reads);
# we deliberately keep a UNIFORM band for the SIMULATION. `high` is P-DEPENDENT, grounded
# on real reads/P (median x/P falls with P: P=25->898, P=15->1497, P=10->2245):
#   P=25 -> U(500P,1000P) ; P=15 -> U(900P,1800P) ; P=10 -> U(1500P,2800P) [P=10 retained
#   only as an extra smaller-P trend point]. `low` = U(20P,40P) for BOTH P.
#
# EVERY cell stores its depth band (depth_lo/depth_hi) EXPLICITLY (not implicit in the
# generator) and a human-readable f-expression string, so a client reads the DGP off the file.
.low_band  <- function(P) c(20, 40)
.high_band <- function(P) if (P >= 25) c(500, 1000) else if (P >= 15) c(900, 1800) else c(1500, 2800)

# Human-readable true-f(t) description (mirrorexp_contB): each off-diagonal edge draws one of
# two continuous-exp shapes, active on half the interval and EXACTLY 0 on the other half;
# the diagonal is set for positive-definiteness (diagonal dominance).
.f_expr <- function(edge_strength, rate) sprintf(
  paste0("mirrorexp_contB: c0=exp(-%g/2); ",
         "f_early(t)=%g*(exp(-%g*t)-c0)/(1-c0) for t<=0.5 else 0; ",
         "f_late(t)=f_early mirrored about t=0.5 (i.e. f_early(1-t)); ",
         "each active edge (i,j) draws early- or late-type (prob 1/2); ",
         "diag(Omega(t))=rowSums(|off-diag(t)|)+base_diag (diagonal-dominant PD)."),
  rate, edge_strength, rate)

PAIRED_CELLS <- local({
  grid <- expand.grid(P = c(15, 25), n = c(12, 20), depth = c("low", "high"), m = c(7L, 15L),
                      stringsAsFactors = FALSE)
  lapply(seq_len(nrow(grid)), function(i) {
    b  <- if (grid$m[i] <= 7L) list(q = 2L, N_n = 1L) else list(q = 3L, N_n = 3L)
    hb <- if (grid$depth[i] == "high") .high_band(grid$P[i]) else .low_band(grid$P[i])
    list(generator = "mirrorexp_contB", P = grid$P[i], n = grid$n[i], m = grid$m[i],
         depth_mode = grid$depth[i], depth_lo = hb[1], depth_hi = hb[2],   # explicit band, every cell
         edge_strength = 0.5, rate = 6, base_diag_mean = 0.1,
         f_expr = .f_expr(0.5, 6),                                          # human-readable true f(t)
         tv_q = b$q, tv_N_n = b$N_n)
  })
})
