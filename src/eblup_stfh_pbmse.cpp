// src/eblup_stfh_pbmse.cpp
// Parametric Bootstrap MSE for Spatio-Temporal Fay-Herriot Model
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "eblup_stfh_core.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// ============================================================================
// Helper: Set R random seed
// ============================================================================
static void set_r_seed(int seed) {
  if (seed >= 0) {
    Rcpp::Environment base("package:base");
    Rcpp::Function setSeed = base["set.seed"];
    setSeed(seed);
  }
}

// ============================================================================
// Helper: Generate SAR(1) random effects via Cholesky
// ============================================================================
static void generate_u1_sar(int D, double sigma21, double rho1,
                           const mat& W, const mat& L_omega1,
                           vec& u1) {
  u1.zeros(D);
  if (sigma21 <= 0) return;

  vec z(D);
  for (int i = 0; i < D; ++i) {
    z(i) = R::rnorm(0.0, 1.0);
  }
  u1 = std::sqrt(sigma21) * (L_omega1 * z);
}

// ============================================================================
// Helper: Generate AR(1) random effects for one domain
// ============================================================================
static void generate_u2_ar1(int Tt, double sigma22, double rho2, vec& u2) {
  u2.zeros(Tt);
  if (sigma22 <= 0) return;

  double var_u2 = sigma22 / (1.0 - rho2 * rho2);
  double sd_init = std::sqrt(var_u2);

  u2(0) = R::rnorm(0.0, sd_init);

  for (int t = 1; t < Tt; ++t) {
    u2(t) = rho2 * u2(t - 1) + R::rnorm(0.0, std::sqrt(sigma22));
  }
}

// ============================================================================
// Parametric Bootstrap MSE for ST-FH
// ============================================================================
// [[Rcpp::export(.pbmse_stfh)]]
List pbmse_stfh(
    const arma::mat& Xall,
    const arma::vec& yall,
    const arma::vec& vardirall,
    const arma::mat& proxmat,
    int D,
    int Tt,
    std::string model = "ST",
    int maxiter = 100,
    double precision = 1e-4,
    int B = 100,
    int n_threads = 0,
    int seed = -1
) {
  if (seed >= 0) set_r_seed(seed);

  const int M = D * Tt;
  const double med_vardir = median(vardirall);

  // ============================================================================
  // 1) Initial fit using the robust .eblup_stfh_core function
  // ============================================================================
  bool isST = (model == "ST");

  List init_fit = eblup_stfh_core(
    Xall, yall, vardirall, proxmat, D, Tt,
    model, maxiter, precision,
    0.5 * med_vardir, 0.5,
    0.5 * med_vardir, 0.5
  );

  if (!init_fit.containsElementNamed("estcoef") || Rf_isNull(init_fit["estcoef"])) {
    Rcpp::stop("Initial fit failed");
  }

  // Extract estimates from initial fit
  DataFrame estvarcomp = as<DataFrame>(init_fit["estvarcomp"]);
  NumericVector est = estvarcomp["estimate"];
  const double sigma21_est = est(0);
  const double rho1_est = est(1);
  const double sigma22_est = est(2);
  const double rho2_est = isST ? est(3) : 0.0;

  // Extract beta and theta
  DataFrame estcoef = as<DataFrame>(init_fit["estcoef"]);
  NumericVector beta_ = estcoef["beta"];
  vec beta_est = as<vec>(beta_);
  vec theta_est = Xall * beta_est;
  DataFrame df_eblup = as<DataFrame>(init_fit["df_eblup"]);
  vec eblup_est = as<vec>(df_eblup["eblup"]);

  bool converged = as<bool>(init_fit["convergence"]);
  int n_iter_init = as<int>(init_fit["n_iter"]);

  // ============================================================================
  // 2) Cholesky of Omega1 for efficient u1 generation
  // ============================================================================
  mat ImrW_init = eye<mat>(D, D) - rho1_est * proxmat;
  mat Omega1_init;
  if (!inv_sympd(Omega1_init, ImrW_init.t() * ImrW_init)) {
    Rcpp::stop("Failed to compute Omega1");
  }
  // Lower factor: Omega1 = L L'  =>  Var(L z) = Omega1. (arma::chol's default
  // is the UPPER factor R with Omega1 = R'R, for which Var(R z) = R R' != Omega1.)
  mat L_omega1;
  if (!chol(L_omega1, Omega1_init, "lower")) {
    Rcpp::stop("Failed to compute Cholesky of Omega1");
  }

  // ============================================================================
  // 3)-6) Rounds of replicates until exactly B valid ones (as in the SFH
  // bootstrap): a replicate whose REML fit does not converge is replaced by a
  // fresh draw, up to 5 * B attempts. Random numbers are drawn SEQUENTIALLY
  // (R RNG) before each parallel round, so results do not depend on n_threads.
  // MSE = mean over replicates of (eblup* - mu*)^2, mu* = X beta^ + u1* + u2*
  // (as in sae::pbmseSTFH).
  // ============================================================================
#ifdef _OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
#endif

  vec mse_pb(M, fill::zeros);
  int n_valid = 0, n_attempted = 0;
  const int max_attempts = 5 * B;

  while (n_valid < B && n_attempted < max_attempts) {
    const int need = std::min(B - n_valid, max_attempts - n_attempted);

    mat U1_boot(D, need);
    mat U2_boot(M, need);
    mat Eps_boot(M, need);
    for (int b = 0; b < need; ++b) {
      vec u1_d(D);
      generate_u1_sar(D, sigma21_est, rho1_est, proxmat, L_omega1, u1_d);
      U1_boot.col(b) = u1_d;

      if (isST) {
        for (int d = 0; d < D; ++d) {
          vec u2_t(Tt);
          generate_u2_ar1(Tt, sigma22_est, rho2_est, u2_t);
          U2_boot.col(b).rows(d * Tt, (d + 1) * Tt - 1) = u2_t;
        }
      } else {
        // Model "S": u2dt iid N(0, sigma22) per domain and time.
        for (int d = 0; d < D; ++d) {
          for (int t = 0; t < Tt; ++t) {
            U2_boot(d * Tt + t, b) = R::rnorm(0.0, std::sqrt(sigma22_est));
          }
        }
      }

      for (int i = 0; i < M; ++i) {
        Eps_boot(i, b) = R::rnorm(0.0, std::sqrt(vardirall(i)));
      }
    }

    mat SqDiff(M, need, fill::zeros);
    std::vector<int> valid(need, 0);

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < need; ++b) {
      vec u1_expanded(M);
      for (int d = 0; d < D; ++d) {
        u1_expanded.subvec(d * Tt, (d + 1) * Tt - 1).fill(U1_boot(d, b));
      }
      const vec mudt_b = theta_est + u1_expanded + U2_boot.col(b);
      const vec y_boot = mudt_b + Eps_boot.col(b);

      // Re-estimate theta by REML in every replicate (as in Marhuenda et al.
      // 2013 / sae::pbmseSTFH); with theta held fixed the bootstrap would only
      // capture g1 + g2 and miss the variance-component uncertainty (g3).
      // Same default starting values as the original fit.
      STFitArma fit_b = stfh_fit_arma(
        Xall, y_boot, vardirall, proxmat, D, Tt, isST,
        maxiter, precision,
        0.5 * med_vardir, 0.5, 0.5 * med_vardir, 0.5
      );

      if (fit_b.ok && fit_b.converged) {
        SqDiff.col(b) = square(fit_b.eblup - mudt_b);
        valid[b] = 1;
      }
    }

    for (int b = 0; b < need; ++b) {
      if (valid[b]) {
        mse_pb += SqDiff.col(b);
        ++n_valid;
      }
    }
    n_attempted += need;
  }

  if (n_valid < B) {
    Rcpp::warning("Only %d of %d bootstrap replicates converged (after %d attempts).",
                  n_valid, B, n_attempted);
  }
  if (n_valid > 0) {
    mse_pb /= n_valid;
  }

  // ============================================================================
  // 7) Assemble output
  // ============================================================================
  List out;
  out["eblup"] = eblup_est;
  out["mse_pb"] = mse_pb;
  out["B"] = n_valid;
  out["B_total"] = n_attempted;
  out["sigma21"] = sigma21_est;
  out["rho1"] = rho1_est;
  out["sigma22"] = sigma22_est;
  out["rho2"] = rho2_est;
  out["beta"] = beta_est;
  out["n_iter"] = n_iter_init;
  out["convergence"] = converged;

  return out;
}
