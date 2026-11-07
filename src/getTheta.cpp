#include <Rcpp.h>
#include <limits>
#include <vector>
using namespace Rcpp;

namespace {

inline double g_transform(const double mu, const double m, const double M) {
  return std::log(mu - m) - std::log(M - mu);
}

inline double clamp_interior(const double x, const double lower, const double upper) {
  if (x < lower) return lower;
  if (x > upper) return upper;
  return x;
}

void compute_bprime_stats(const NumericVector& spt,
                          const NumericVector& f0,
                          const NumericVector& theta,
                          const IntegerVector& indices,
                          NumericVector& bPrime,
                          NumericVector& bPrime2,
                          NumericVector& bPrimeErr,
                          const NumericVector& mu) {
  const int sptN = spt.size();

  // log(f0) and the per-observation weights are computed once and reused.
  std::vector<double> log_f0(sptN), weight(sptN);
  for (int j = 0; j < sptN; ++j) log_f0[j] = std::log(f0[j]);

  for (int k = 0; k < indices.size(); ++k) {
    const int idx = indices[k];
    const double theta_i = theta[idx];

    double max_log = R_NegInf;
    for (int j = 0; j < sptN; ++j) {
      const double log_term = log_f0[j] + spt[j] * theta_i;
      if (log_term > max_log) {
        max_log = log_term;
      }
    }

    double denom = 0.0;
    double numer = 0.0;
    for (int j = 0; j < sptN; ++j) {
      weight[j] = std::exp(log_f0[j] + spt[j] * theta_i - max_log);
      denom += weight[j];
      numer += spt[j] * weight[j];
    }

    const double mean_val = numer / denom;

    double var_numer = 0.0;
    for (int j = 0; j < sptN; ++j) {
      const double diff = spt[j] - mean_val;
      var_numer += diff * diff * weight[j];
    }

    bPrime[idx] = mean_val;
    bPrime2[idx] = var_numer / denom;
    bPrimeErr[idx] = mean_val - mu[idx];
  }
}

IntegerVector which_false(const LogicalVector& x) {
  std::vector<int> out;
  out.reserve(x.size());
  for (int i = 0; i < x.size(); ++i) {
    if (!x[i]) {
      out.push_back(i);
    }
  }
  return wrap(out);
}

}  // namespace

// [[Rcpp::export]]
List getTheta_cpp(const NumericVector& spt,
                  const NumericVector& f0,
                  const NumericVector& mu,
                  Nullable<NumericVector> thetaStart = R_NilValue,
                  const double eps = 1e-10,
                  const int maxiter = 100,
                  const int maxhalf = 20,
                  const double maxtheta = 500.0,
                  const bool logit = true) {
  const int sptN = spt.size();
  const int n = mu.size();
  const double boundary_eps = std::max(1e-12 * (max(spt) - min(spt)), std::numeric_limits<double>::epsilon());
  const double var_eps = 1e-12;

  if (f0.size() != sptN) {
    stop("spt and f0 must be vectors of equal length.");
  }
  for (int j = 0; j < sptN; ++j) {
    if (f0[j] < 0) {
      stop("f0 values cannot be negative.");
    }
  }

  const double m = min(spt);
  const double M = max(spt);
  const double safe_lower = m + boundary_eps;
  const double safe_upper = M - boundary_eps;
  if (min(mu) < m || max(mu) > M) {
    stop("mu starting values must lie within the range of spt.");
  }

  NumericVector theta(n);
  if (thetaStart.isNotNull()) {
    NumericVector theta_start_vec(thetaStart);
    if (theta_start_vec.size() != n) {
      stop("thetaStart must be a vector with length equal length(mu)");
    }
    theta = clone(theta_start_vec);
  } else {
    std::fill(theta.begin(), theta.end(), 0.0);
  }

  NumericVector gMu(n);
  for (int i = 0; i < n; ++i) {
    const double mu_safe = clamp_interior(mu[i], safe_lower, safe_upper);
    gMu[i] = g_transform(mu_safe, m, M);
  }

  NumericVector thetaOld(n, NA_REAL);
  NumericVector bPrimeErrOld(n, NA_REAL);
  NumericVector bPrime(n), bPrime2(n), bPrimeErr(n);
  LogicalVector conv(n, false);
  LogicalVector maxedOut(n, false);

  IntegerVector all_idx(n);
  for (int i = 0; i < n; ++i) {
    all_idx[i] = i;
  }
  compute_bprime_stats(spt, f0, theta, all_idx, bPrime, bPrime2, bPrimeErr, mu);

  for (int i = 0; i < n; ++i) {
    conv[i] = (std::abs(bPrimeErr[i]) < eps) ||
      (theta[i] == maxtheta && bPrimeErr[i] < 0) ||
      (theta[i] == -maxtheta && bPrimeErr[i] > 0);
  }

  IntegerVector s = which_false(conv);
  int iter = 0;

  while (s.size() > 0 && iter < maxiter) {
    ++iter;

    for (int k = 0; k < s.size(); ++k) {
      const int idx = s[k];
      bPrimeErrOld[idx] = bPrimeErr[idx];
      thetaOld[idx] = theta[idx];

      double theta_new = theta[idx];
      if (logit) {
        const double bPrimeSafe = clamp_interior(bPrime[idx], safe_lower, safe_upper);
        const double bPrime2Safe = std::max(bPrime2[idx], var_eps);
        const double tPrime = (M - m) / ((bPrimeSafe - m) * (M - bPrimeSafe)) * bPrime2Safe;
        if (R_finite(tPrime) && std::abs(tPrime) > var_eps) {
          theta_new = theta[idx] - (g_transform(bPrimeSafe, m, M) - gMu[idx]) / tPrime;
        } else {
          theta_new = theta[idx] - bPrimeErr[idx] / std::max(bPrime2[idx], var_eps);
        }
      } else {
        theta_new = theta[idx] - bPrimeErr[idx] / std::max(bPrime2[idx], var_eps);
      }

      if (theta_new > maxtheta) theta_new = maxtheta;
      if (theta_new < -maxtheta) theta_new = -maxtheta;
      theta[idx] = theta_new;
    }

    IntegerVector ss = clone(s);
    int nhalf = 0;

    while (ss.size() > 0 && nhalf < maxhalf) {
      compute_bprime_stats(spt, f0, theta, ss, bPrime, bPrime2, bPrimeErr, mu);

      std::vector<int> bad_idx;
      bad_idx.reserve(ss.size());
      for (int k = 0; k < ss.size(); ++k) {
        const int idx = ss[k];
        if (std::abs(bPrimeErr[idx]) > std::abs(bPrimeErrOld[idx])) {
          bad_idx.push_back(idx);
        }
      }

      if (!bad_idx.empty()) {
        ++nhalf;
        for (size_t k = 0; k < bad_idx.size(); ++k) {
          const int idx = bad_idx[k];
          theta[idx] = (theta[idx] + thetaOld[idx]) / 2.0;
        }
      }

      ss = wrap(bad_idx);
    }

    for (int k = 0; k < ss.size(); ++k) {
      const int idx = ss[k];
      maxedOut[idx] = true;
      theta[idx] = thetaOld[idx];
    }

    for (int k = 0; k < s.size(); ++k) {
      const int idx = s[k];
      conv[idx] = (std::abs(bPrimeErr[idx]) < eps) ||
        (theta[idx] == maxtheta && bPrimeErr[idx] < 0) ||
        (theta[idx] == -maxtheta && bPrimeErr[idx] > 0);
    }

    std::vector<int> next_s;
    next_s.reserve(s.size());
    for (int k = 0; k < s.size(); ++k) {
      const int idx = s[k];
      if (!conv[idx] && !maxedOut[idx]) {
        next_s.push_back(idx);
      }
    }
    s = wrap(next_s);
  }

  return List::create(
    Named("theta") = theta,
    Named("bPrime") = bPrime,
    Named("bPrime2") = bPrime2,
    Named("conv") = conv,
    Named("iter") = iter
  );
}
