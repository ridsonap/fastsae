// src/eblup_stfh_core.h
// Header for ST-FH core functions

#ifndef EBLUP_STFH_CORE_H
#define EBLUP_STFH_CORE_H

#include <RcppArmadillo.h>

struct STFitArma {
  bool ok;
  int n_iter;
  bool converged;
  arma::vec beta;
  double sigma21;
  double rho1;
  double sigma22;
  double rho2;
  arma::vec eblup;
  arma::vec theta;
};

// Rcpp wrap specialization for STFitArma (inline to avoid duplicate symbols)
namespace Rcpp {
  template <> inline SEXP wrap(const STFitArma& s) {
    List out;
    out["ok"] = s.ok;
    out["n_iter"] = s.n_iter;
    out["converged"] = s.converged;
    out["beta"] = s.beta;
    out["sigma21"] = s.sigma21;
    out["rho1"] = s.rho1;
    out["sigma22"] = s.sigma22;
    out["rho2"] = s.rho2;
    out["eblup"] = s.eblup;
    out["theta"] = s.theta;
    return out;
  }
}

// Exported function declarations
Rcpp::List eblup_stfh_core(
    const arma::mat& X,
    const arma::vec& y,
    const arma::vec& vardir,
    const arma::mat& proxmat,
    int D,
    int Tt,
    std::string model,
    int maxiter,
    double precision,
    double sigma21_start,
    double rho1_start,
    double sigma22_start,
    double rho2_start
);

Rcpp::List pbmse_stfh(
    const arma::mat& X,
    const arma::vec& y,
    const arma::vec& vardir,
    const arma::mat& proxmat,
    int D,
    int Tt,
    std::string model,
    int maxiter,
    double precision,
    int B,
    int n_threads,
    int seed
);

// Internal function (not exported) - declaration in eblup_stfh_pbmse.cpp
// STFitArma fit_stfh_arma(...)

#endif // EBLUP_STFH_CORE_H
