// ============================================================================
// tvcglasso.cpp — Rcpp helper for the TVCGLasso engine (publication version).
//
// Only G_beta_Rcpp is kept: it assembles the precision slices
//   Omega(t_k) = sum_h B_h(t_k) * beta_h
// from the per-edge B-spline coefficients. It is the only compiled function the
// public API (main_function_final / tv_warm_path / the refit) actually calls.
// The five other exports in the original main.cpp (generate_alpha_tau_Rcpp,
// derivative_first_penalty_Rcpp, S_Z_t_Rcpp, the two NaN-replacers) were dead in
// the final engine (the penalty is reimplemented in R) and are dropped.
//
// Copied verbatim from code/main.cpp (lines 261-292); behavior is unchanged.
// ============================================================================
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List G_beta_Rcpp(const List& beta, const NumericMatrix& x_b_spline_base, int m, double tol = 1e-3) {
  int P = as<NumericMatrix>(beta[0]).nrow();
  int J_n = beta.size();
  List G_temp(m);

  for (int i = 0; i < m; ++i) {
    NumericMatrix G_t(P, P);
    for (int j = 0; j < P; ++j) {
      for (int q = j; q < P; ++q) {
        NumericVector beta_temp(J_n);
        for (int z = 0; z < J_n; ++z) {
          NumericMatrix beta_z = as<NumericMatrix>(beta[z]);
          beta_temp[z] = beta_z(j, q);
        }
        double sum_val = 0;
        for (int k = 0; k < J_n; ++k) {
          sum_val += x_b_spline_base(i, k) * beta_temp[k];
        }
        G_t(j, q) = sum_val;
      }
    }
    // mirror the upper triangle into the lower triangle -> symmetric Omega(t_i)
    for (int j = 0; j < P; ++j) {
      for (int q = 0; q < j; ++q) {
        G_t(j, q) = G_t(q, j);
      }
    }
    G_temp[i] = G_t;
  }
  return G_temp;
}
