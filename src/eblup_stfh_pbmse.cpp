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
// Helper: Build Omega2 (AR(1) covariance matrix) and its derivative
// ============================================================================
static void build_omega2_with_deriv(int Tt, double rho2,
                                    arma::mat& Omega2,
                                    arma::mat& dOmega2) {
  Omega2.zeros(Tt, Tt);
  dOmega2.zeros(Tt, Tt);
  const double one_m_rho2sq = 1.0 - rho2 * rho2;
  const double denom = one_m_rho2sq;
  const double denom_sq = denom * denom;

  for (int i = 0; i < Tt; ++i) {
    for (int j = 0; j < i; ++j) {
      const int lag = std::abs(i - j);
      const double rho_lag = std::pow(rho2, (double)lag);
      const double val = rho_lag / denom;
      Omega2(i, j) = val;
      Omega2(j, i) = val;

      double raw_deriv = (lag > 0) ? lag * std::pow(rho2, (double)(lag - 1)) : 0.0;
      double dval = (raw_deriv * denom - rho_lag * (-2.0 * rho2)) / denom_sq;
      dOmega2(i, j) = dval;
      dOmega2(j, i) = dval;
    }
  }
  Omega2.diag().fill(1.0 / denom);
  dOmega2.diag().fill(2.0 * rho2 / denom_sq);
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
// Pure Armadillo ST-FH EBLUP computation for bootstrap
// Uses FIXED theta from initial fit - no re-estimation
// ============================================================================
static STFitArma fit_stfh_eblup_only(
    const arma::mat& X,
    const arma::vec& y,
    const arma::vec& vardir,
    const arma::mat& W,
    int D,
    int Tt,
    bool isST,
    double sigma21_est,
    double rho1_est,
    double sigma22_est,
    double rho2_est
) {
  STFitArma res;
  res.ok = false;
  res.n_iter = 0;
  res.converged = true;
  res.sigma21 = sigma21_est;
  res.rho1 = rho1_est;
  res.sigma22 = sigma22_est;
  res.rho2 = isST ? rho2_est : 0.0;

  const int M = D * Tt;

  if ((int)X.n_rows != M || (int)y.n_elem != M || (int)vardir.n_elem != M) {
    return res;
  }

  // Build matrices
  mat Id = eye<mat>(D, D);
  mat tX = X.t();

  // Build Omega1 = (I - rho1*W)^{-1} (I - rho1*W')^{-1}
  mat ImrW = Id - rho1_est * W;
  mat Omega1;
  if (!inv_sympd(Omega1, ImrW.t() * ImrW)) return res;

  // Build Omega2 for ST model
  mat Omega2, dOmega2;
  if (isST) {
    build_omega2_with_deriv(Tt, rho2_est, Omega2, dOmega2);
  } else {
    Omega2 = eye<mat>(Tt, Tt);
  }

  // Build Vu1
  mat Vu1 = sigma21_est * Omega1;
  mat invVu1;
  if (sigma21_est > 1e-10) {
    if (!inv_sympd(invVu1, Vu1)) return res;
  } else {
    invVu1.zeros(D, D);
  }

  // Build invA (block diagonal M x M)
  mat invA(M, M, fill::zeros);
  for (int d = 0; d < D; ++d) {
    const int first = d * Tt;
    const int last = first + Tt - 1;
    mat Ved = diagmat(vardir.subvec(first, last));
    mat Ad = sigma22_est * Omega2 + Ved;
    mat invAd;
    if (!inv_sympd(invAd, Ad)) return res;
    invA.submat(first, first, last, last) = invAd;
  }

  // Build invV
  mat invV;
  if (sigma21_est <= 1e-10) {
    // V = A when sigma21 = 0
    invV = invA;
  } else {
    // Build Z1' invA
    mat invAZ1(M, D, fill::zeros);
    for (int d = 0; d < D; ++d) {
      const int first = d * Tt, last = first + Tt - 1;
      invAZ1.submat(first, d, last, d) = sum(invA.submat(first, first, last, last), 1);
    }

    // C = Vu1^{-1} + Z1' A^{-1} Z1
    vec diagC(D);
    for (int d = 0; d < D; ++d) {
      diagC(d) = dot(invAZ1.col(d), invAZ1.col(d));
    }
    mat Cmat = invVu1 + diagmat(diagC);
    mat invC;
    if (!inv_sympd(invC, Cmat)) return res;

    // Woodbury identity
    invV = invA - invAZ1 * invC * invAZ1.t();
  }

  // Fixed effects beta
  mat tXinvV = tX * invV;
  mat tXinvVX = tXinvV * X;
  mat Qxx;
  if (!inv_sympd(Qxx, tXinvVX)) return res;
  res.beta = Qxx * (tXinvV * y);

  // Residuals
  vec resid = y - X * res.beta;
  vec invVresid = invV * resid;

  // tZ1' invV resid
  vec tZ1invVresid(D);
  for (int d = 0; d < D; ++d) {
    const int first = d * Tt, last = first + Tt - 1;
    tZ1invVresid(d) = sum(invVresid.subvec(first, last));
  }

  // Random effects
  vec u1 = Vu1 * tZ1invVresid;
  vec u2(M);
  if (isST) {
    for (int d = 0; d < D; ++d) {
      const int first = d * Tt, last = first + Tt - 1;
      u2.subvec(first, last) = sigma22_est * Omega2 * invVresid.subvec(first, last);
    }
  } else {
    u2 = sigma22_est * invVresid;
  }

  // Expand u1 to domain-major order
  vec u1_expanded(M);
  for (int d = 0; d < D; ++d) {
    u1_expanded.subvec(d * Tt, (d + 1) * Tt - 1).fill(u1(d));
  }

  res.eblup = X * res.beta + u1_expanded + u2;
  res.theta = X * res.beta;
  res.ok = true;

  return res;
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
    stop("Initial fit failed");
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
    stop("Failed to compute Omega1");
  }
  mat L_omega1;
  if (!chol(L_omega1, Omega1_init)) {
    stop("Failed to compute Cholesky of Omega1");
  }

  // ============================================================================
  // 3) Generate bootstrap random effects and errors SEQUENTIALLY
  // ============================================================================
  mat U1_boot(D, B);
  mat U2_boot(M, B);
  mat Eps_boot(M, B);

  for (int b = 0; b < B; ++b) {
    vec u1_d(D);
    generate_u1_sar(D, sigma21_est, rho1_est, proxmat, L_omega1, u1_d);
    U1_boot.col(b) = u1_d;

    if (isST) {
      cube U2_st(Tt, D, 1);
      for (int d = 0; d < D; ++d) {
        vec u2_t(Tt);
        generate_u2_ar1(Tt, sigma22_est, rho2_est, u2_t);
        U2_st.slice(0).col(d) = u2_t;
      }
      for (int d = 0; d < D; ++d) {
        U2_boot.rows(d * Tt, (d + 1) * Tt - 1).col(b) = U2_st.slice(0).col(d);
      }
    } else {
      vec u2_d(Tt);
      for (int t = 0; t < Tt; ++t) {
        u2_d(t) = R::rnorm(0.0, std::sqrt(sigma22_est));
      }
      for (int d = 0; d < D; ++d) {
        U2_boot.rows(d * Tt, (d + 1) * Tt - 1).col(b) = u2_d;
      }
    }

    for (int i = 0; i < M; ++i) {
      Eps_boot(i, b) = R::rnorm(0.0, std::sqrt(vardirall(i)));
    }
  }

  // ============================================================================
  // 4) Storage for bootstrap EBLUP replicates
  // ============================================================================
  mat Eblup_boot(M, B, fill::zeros);
  std::vector<int> valid_count(B, 0);

#ifdef _OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
#endif

  // ============================================================================
  // 5) Parallel bootstrap loop
  // ============================================================================
#pragma omp parallel for schedule(dynamic)
  for (int b = 0; b < B; ++b) {
    // Construct bootstrap y
    vec u1_b = U1_boot.col(b);
    vec u2_b = U2_boot.col(b);

    vec u1_expanded(M);
    for (int d = 0; d < D; ++d) {
      u1_expanded.subvec(d * Tt, (d + 1) * Tt - 1).fill(u1_b(d));
    }

    vec eps_boot = Eps_boot.col(b);

    // y_boot = theta + u1_expanded + u2 + epsilon
    vec y_boot = theta_est + u1_expanded + u2_b + eps_boot;

    // Compute EBLUP using fixed theta (no re-estimation)
    STFitArma fit_b = fit_stfh_eblup_only(
      Xall, y_boot, vardirall, proxmat, D, Tt,
      isST,
      sigma21_est, rho1_est,
      sigma22_est, rho2_est
    );

    if (fit_b.ok) {
      Eblup_boot.col(b) = fit_b.eblup;
      valid_count[b] = 1;
    } else {
      Eblup_boot.col(b) = theta_est;
      valid_count[b] = 0;
    }
  }

  // ============================================================================
  // 6) Compute MSE (matching sae::pbmseSTFH formula)
  // ---------------------------------------------------------------------------
  // In sae: MSE = mean((theta_boot - mudt_b)^2)
  // where mudt_b = X*beta + u1_expanded + u2  (using ORIGINAL beta)
  // ============================================================================
  vec mse_pb(M, fill::zeros);
  int n_valid = 0;
  for (int b = 0; b < B; ++b) {
    if (valid_count[b] == 1) {
      // mudt_b = theta_est + u1_expanded + u2_b
      vec u1_b = U1_boot.col(b);
      vec u2_b = U2_boot.col(b);

      vec u1_expanded(M);
      for (int d = 0; d < D; ++d) {
        u1_expanded.subvec(d * Tt, (d + 1) * Tt - 1).fill(u1_b(d));
      }

      vec mudt_b = theta_est + u1_expanded + u2_b;
      mse_pb += square(Eblup_boot.col(b) - mudt_b);
      n_valid++;
    }
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
  out["B_total"] = B;
  out["sigma21"] = sigma21_est;
  out["rho1"] = rho1_est;
  out["sigma22"] = sigma22_est;
  out["rho2"] = rho2_est;
  out["beta"] = beta_est;
  out["n_iter"] = n_iter_init;
  out["convergence"] = converged;

  return out;
}
