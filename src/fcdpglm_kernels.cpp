#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// C++ versions of the inner loops of R/fcdpglm.R, used when
// fit_rat_dpglm(..., use_cpp_kernels = TRUE).

// b_i(theta_i; w) = log sum_k w_k exp(theta_i k), via log-sum-exp per row.
// [[Rcpp::export]]
NumericVector rat_btheta_cpp(const NumericVector& theta, const NumericVector& atoms,
                             const NumericVector& logw) {
  const int n = theta.size(), K = atoms.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double mx = R_NegInf;
    for (int k = 0; k < K; ++k) {
      const double v = logw[k] + theta[i] * atoms[k];
      if (v > mx) mx = v;
    }
    double s = 0.0;
    for (int k = 0; k < K; ++k) s += std::exp(logw[k] + theta[i] * atoms[k] - mx);
    out[i] = mx + std::log(s);
  }
  return out;
}

// log psi(k) = log sum_i u_i exp(theta_i k), via log-sum-exp per support point.
// [[Rcpp::export]]
NumericVector rat_log_psi_cpp(const NumericVector& log_u, const NumericVector& theta,
                              const NumericVector& atoms) {
  const int n = theta.size(), K = atoms.size();
  NumericVector out(K);
  for (int k = 0; k < K; ++k) {
    double mx = R_NegInf;
    for (int i = 0; i < n; ++i) {
      const double v = log_u[i] + theta[i] * atoms[k];
      if (v > mx) mx = v;
    }
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += std::exp(log_u[i] + theta[i] * atoms[k] - mx);
    out[k] = mx + std::log(s);
  }
  return out;
}

// Variance of each tilted distribution, V_i = sum_k k^2 p_ik - (sum_k k p_ik)^2,
// with p_ik = exp(logw_k + theta_i k - b_i).
// [[Rcpp::export]]
NumericVector rat_tilted_var_cpp(const NumericVector& theta, const NumericVector& atoms,
                                 const NumericVector& logw, const NumericVector& btheta) {
  const int n = theta.size(), K = atoms.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double m1 = 0.0, m2 = 0.0;
    for (int k = 0; k < K; ++k) {
      const double p = std::exp(logw[k] + theta[i] * atoms[k] - btheta[i]);
      m1 += p * atoms[k];
      m2 += p * atoms[k] * atoms[k];
    }
    out[i] = m2 - m1 * m1;
  }
  return out;
}
