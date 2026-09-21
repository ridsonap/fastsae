// src/eblup_sfh_pbmse.cpp
// Parametric Bootstrap MSE for Spatial Fay-Herriot Model
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "eblup_sfh.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export(.seblup_pbmse)]]
List seblup_pbmse(
    const arma::mat& X,
    const arma::vec& y,
    const arma::vec& vardir,
    const arma::mat& W,
    std::string method = "REML",
    int maxiter = 100,
    double precision = 1e-4,
    int B = 100,
    int n_threads = 0,
    int max_attempts_factor = 5,
    int seed = -1
) {
  if (seed >= 0) set_r_seed(seed);

  // ================================================================
  // Main fit
  // ================================================================
  List result_awal = seblup_core(
    X, y, vardir, W,
    method, maxiter, precision,
    true
  );

  if (!(bool) result_awal["convergence"]) {
    Rcpp::warning("Initial fit did not converge; bootstrap MSE results may be unreliable.");
  }

  const int m = X.n_rows;

  double sigma2_boot = result_awal["sigma2_u"];
  double rho_boot = result_awal["rho"];
  mat beta_boot = result_awal["beta"];
  mat Xbeta = result_awal["Xbeta"];
  mat GVi = result_awal["GVi"];
  mat Q_mat = result_awal["Q"];
  mat XtVi_mat = result_awal["XtVi"];
  vec g1 = result_awal["g1d"];
  vec g2 = result_awal["g2d"];

  mat I = eye<mat>(m, m);
  mat Irhoproxmat = I - rho_boot * W;
  mat QXtVi = Q_mat * XtVi_mat;

  // ================================================================
  // Batch scheme: repeat until exactly B valid replicates
  // ================================================================
  vec sum_mse(m, fill::zeros);
  vec sum_g3(m, fill::zeros);
  vec sum_g1(m, fill::zeros);
  vec sum_g2(m, fill::zeros);

  int n_valid_total = 0;
  int n_attempted_total = 0;
  int round_no = 0;
  const int max_total_attempts = B * max_attempts_factor;

#ifdef _OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
#endif

  while (n_valid_total < B) {
    ++round_no;
    int need = B - n_valid_total;

    if (n_attempted_total + need > max_total_attempts) {
      int remaining_budget = max_total_attempts - n_attempted_total;
      if (remaining_budget <= 0) {
        Rcpp::stop("Failed to collect %d valid replicates after %d total attempts "
             "(too many replicates failed to converge).", B, n_attempted_total);
      }
      need = remaining_budget;
    }

    // Generate random numbers SEQUENTIALLY for this batch
    mat U_boot(m, need);
    mat E_boot(m, need);
    for (int b = 0; b < need; ++b) {
      U_boot.col(b) = as<arma::vec>(rnorm(m, 0.0, std::sqrt(sigma2_boot)));
      for (int i = 0; i < m; ++i) {
        E_boot(i, b) = R::rnorm(0.0, std::sqrt(vardir(i)));
      }
    }

    mat V_boot = arma::solve(Irhoproxmat, U_boot);

    mat SqDiff_r(m, need, fill::zeros);
    mat SqDiffPb_r(m, need, fill::zeros);
    mat G1_r(m, need, fill::zeros);
    mat G2_r(m, need, fill::zeros);
    vec valid_flag_r(need, fill::zeros);

#pragma omp parallel for schedule(dynamic)
    for (int b = 0; b < need; ++b) {
      vec v_boot = V_boot.col(b);
      vec theta_boot = X * beta_boot + v_boot;
      vec direct_boot = theta_boot + E_boot.col(b);

      SeblupFitArma res = seblup_fit_core_arma(
        X, direct_boot, vardir, W, method, maxiter, precision
      );

      bool valid = res.converged &&
        res.sigma2 >= 0.0 &&
        res.theta.is_finite() &&
        res.g1d.is_finite() &&
        res.g2d.is_finite();

      if (!valid) continue;

      // sblup estimate using initial fit weights
      vec Bstim_sblup = QXtVi * direct_boot;
      vec theta_sblup = X * Bstim_sblup + GVi * (direct_boot - X * Bstim_sblup);

      SqDiff_r.col(b) = arma::square(res.theta - theta_boot);
      SqDiffPb_r.col(b) = arma::square(res.theta - theta_sblup);
      G1_r.col(b) = res.g1d;
      G2_r.col(b) = res.g2d;
      valid_flag_r(b) = 1.0;
    }

    int n_valid_round = (int) arma::accu(valid_flag_r);

    sum_mse += arma::sum(SqDiff_r, 1);
    sum_g3 += arma::sum(SqDiffPb_r, 1);
    sum_g1 += arma::sum(G1_r, 1);
    sum_g2 += arma::sum(G2_r, 1);

    n_valid_total += n_valid_round;
    n_attempted_total += need;

    if (n_valid_round == 0) {
      Rcpp::warning(
        "Round %d: 0 out of %d replicates valid (all failed to converge).",
        round_no, need
      );
    }
  }

  // ================================================================
  // Sequential aggregation
  // ================================================================
  vec mse_pb = sum_mse / (double) B;
  vec g3_pb = sum_g3 / (double) B;
  vec g1_pb = sum_g1 / (double) B;
  vec g2_pb = sum_g2 / (double) B;

  // bias corrected
  vec mse_pb2 = 2.0 * (g1 + g2) - g1_pb - g2_pb + g3_pb;

  result_awal = seblup_core(
    X, y, vardir, W,
    method, maxiter, precision,
    false
  );

  DataFrame df_eblup = Rcpp::as<Rcpp::DataFrame>(result_awal["df_eblup"]);
  df_eblup.push_back(mse_pb, "mse_pb");
  df_eblup.push_back(mse_pb2, "mse_pbbc");
  result_awal["df_eblup"] = df_eblup;

  result_awal.push_back(round_no, "n_rounds");
  result_awal.push_back(n_attempted_total, "n_attempted_total");

  return result_awal;
}
