// src/eblup_stfh_core.h
// Header for ST-FH core functions

#ifndef EBLUP_STFH_CORE_H
#define EBLUP_STFH_CORE_H

#include <RcppArmadillo.h>
#include <vector>

// Observed cells of the D x T grid (domain-major). The model is defined on the
// full grid; only the cells in `obs` have a direct estimate. Rows of the
// "compact" data (X, y, vardir of the observed cells) are domain-major too:
// domain d owns compact rows off[d] .. off[d+1]-1, at periods tix[d].
struct StIndex {
  int D = 0, T = 0, n = 0;
  bool full = true;                 // every cell observed (balanced panel)
  std::vector<int> off;             // D + 1 offsets into the compact rows
  std::vector<arma::uvec> tix;      // observed periods of each domain
  arma::uvec obs;                   // grid positions (d*T + t) of the compact rows
};

// Build the index from the observed-cell mask of the D x T grid (domain-major).
StIndex st_make_index(int D, int T, const arma::uvec& observed);

struct STFitArma {
  bool ok;
  int n_iter;
  bool converged;
  arma::vec beta;
  double sigma21;
  double rho1;
  double sigma22;
  double rho2;
  arma::vec eblup;   // all D x T cells (BLUP for the unobserved ones)
  arma::vec theta;   // X_full beta
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

// Exported function declarations (full-grid inputs; NA y = unobserved cell)
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

// Pure Armadillo REML fit + EBLUP (thread-safe; defined in eblup_stfh_core.cpp).
// Used to re-estimate theta in every bootstrap replicate. X, y, vardir are the
// compact (observed) data; Xfull holds the covariates of all D x T cells, and
// the returned eblup covers all of them.
STFitArma stfh_fit_arma(
    const arma::mat& X, const arma::vec& y, const arma::vec& vardir,
    const arma::mat& W, const StIndex& ix, bool isST,
    int maxiter, double precision,
    double sigma21_start, double rho1_start,
    double sigma22_start, double rho2_start,
    const arma::mat& Xfull
);

#endif // EBLUP_STFH_CORE_H
