// src/seblup_core.h
// Forward declarations for seblup functions

#ifndef SEBLUP_CORE_H
#define SEBLUP_CORE_H

#include <RcppArmadillo.h>

struct SeblupFitArma {
  arma::vec theta;
  arma::vec g1d;
  arma::vec g2d;
  double sigma2;
  double rho_fix;
  int n_iter;
  bool converged;
};

void set_r_seed(int seed);

// Marginal SAR covariance of the sampled domains from the FULL W
// (see eblup_sfh_core.cpp). Pure Armadillo.
bool sar_marginal(double rho, const arma::mat& Wf, const arma::mat& WtWf,
                  const arma::uvec& s, bool all,
                  arma::mat& Sig, arma::mat& SCS,
                  arma::mat* SCSCS, arma::mat* SWWS, arma::mat* Sfull,
                  arma::mat* SCS_full = nullptr, arma::mat* D2_full = nullptr);

// X, y, vardir: sampled domains only; W: FULL weight matrix; idx_s: rows of
// W that are sampled (0-based).
SeblupFitArma seblup_fit_core_arma(
    const arma::mat& X,
    const arma::vec& y,
    const arma::vec& vardir,
    const arma::mat& W,
    const arma::uvec& idx_s,
    const std::string& method,
    int maxiter,
    double precision
);

// Forward declare seblup_core (Rcpp export)
Rcpp::List seblup_core(
    const arma::mat& Xall,
    const arma::vec& yall,
    const arma::vec& vardirall,
    const arma::mat& Wall,
    std::string method,
    int maxiter,
    double precision,
    bool only_core
);

#endif // SEBLUP_CORE_H
