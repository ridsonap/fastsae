// src/sar_engine.h
// Factorisations for the SAR Fay-Herriot model on D domains (m sampled):
//   A(rho) = (I - rho W')(I - rho W)          (precision of the SAR effects)
//   M      = A + s2 E' Psi^-1 E              (posterior precision of v / s2)
// with E the m x D selection of the sampled domains. Everything the SFH fit
// and MSE need follows from solves with A and M (see eblup_sfh_core.cpp),
// e.g. V^-1 = Psi^-1 - s2 Psi^-1 E M^-1 E' Psi^-1 and V^-1 E Sigma = Psi^-1 E M^-1,
// so no D x D inverse or D^3 product is formed.
//
// Sparse W: Eigen simplicial LDL' with AMD ordering on the fixed pattern of
// I + W + W' + W'W (symbolic analysis once). Dense W (more than 25% non-zero
// in that pattern): Eigen dense LLT / LDLT. Pure C++ (thread-safe per object).
#ifndef SAR_ENGINE_H
#define SAR_ENGINE_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <cmath>
#include <vector>

class SarEngine {
 public:
  SarEngine(const arma::mat& W, const arma::uvec& s) : s_(s), D_(W.n_rows) {
    Wsp_ = arma::sp_mat(W);
    Wssp_ = arma::sp_mat(W + W.t());
    WtWsp_ = Wsp_.t() * Wsp_;
    const arma::sp_mat pat = arma::speye(D_, D_) + arma::abs(Wssp_) + arma::abs(WtWsp_);
    const double nnz_lower = 0.5 * (pat.n_nonzero + D_);
    sparse_ = nnz_lower <= 0.25 * 0.5 * D_ * (D_ + 1.0);
    if (sparse_) build_pattern();
    else { Wd_ = W; WtWd_ = W.t() * W; }
  }

  bool is_sparse() const { return sparse_; }
  arma::uword D() const { return D_; }
  const arma::uvec& s() const { return s_; }

  // Factorise A(rho) and M = A + s2 E' Psi^-1 E (psi: sampled domains)
  bool factor(double s2, double rho, const arma::vec& psi) {
    rho_ = rho;
    arma::vec dpsi(D_, arma::fill::zeros);
    dpsi.elem(s_) = s2 / psi;
    if (sparse_) {
      double* a = Asp_.valuePtr();
      double* mv = Msp_.valuePtr();
      const double r2 = rho * rho;
      for (size_t k = 0; k < vI_.size(); ++k) {
        a[k] = vI_[k] - rho * vWs_[k] + r2 * vWtW_[k];
        mv[k] = a[k];
      }
      for (arma::uword i = 0; i < D_; ++i) mv[diag_pos_[i]] += dpsi(i);
      ldA_.factorize(Asp_);
      if (ldA_.info() != Eigen::Success) return false;
      ldM_.factorize(Msp_);
      if (ldM_.info() != Eigen::Success) return false;
      logdetA_ = sum_log(ldA_.vectorD(), &A_pd_);
      logdetM_ = sum_log(ldM_.vectorD(), &M_pd_);
      return A_pd_ && std::isfinite(logdetM_);
    }
    arma::mat A = arma::symmatu(arma::eye(D_, D_) - rho * (Wd_ + Wd_.t()) + (rho * rho) * WtWd_);
    arma::mat M = A;
    M.diag() += dpsi;
    Eigen::Map<Eigen::MatrixXd> Ae(A.memptr(), D_, D_), Me(M.memptr(), D_, D_);
    dA_.compute(Ae);
    if (dA_.info() != Eigen::Success || !dA_.isPositive()) return false;
    dM_.compute(Me);
    if (dM_.info() != Eigen::Success) return false;
    A_pd_ = true;
    logdetA_ = sum_log(dA_.vectorD(), &A_pd_);
    logdetM_ = sum_log(dM_.vectorD(), &M_pd_);
    return A_pd_ && std::isfinite(logdetM_);
  }

  double logdetA() const { return logdetA_; }
  double logdetM() const { return logdetM_; }
  bool M_pd() const { return M_pd_; }

  arma::mat solveA(const arma::mat& B) const { return sparse_ ? solve_(ldA_, B) : solve_(dA_, B); }
  arma::mat solveM(const arma::mat& B) const { return sparse_ ? solve_(ldM_, B) : solve_(dM_, B); }

  // C(rho) B, C = dA/drho = 2 rho W'W - (W + W')
  arma::mat mulC(const arma::mat& B) const {
    if (sparse_) return arma::mat(2.0 * rho_ * (WtWsp_ * B) - Wssp_ * B);
    return 2.0 * rho_ * (WtWd_ * B) - (Wd_ * B + Wd_.t() * B);
  }
  // A(rho) B
  arma::mat mulA(const arma::mat& B) const {
    if (sparse_) return arma::mat(B - rho_ * (Wssp_ * B) + (rho_ * rho_) * (WtWsp_ * B));
    return B - rho_ * (Wd_ * B + Wd_.t() * B) + (rho_ * rho_) * (WtWd_ * B);
  }
  // W B and W'W B
  arma::mat mulW(const arma::mat& B) const { return sparse_ ? arma::mat(Wsp_ * B) : arma::mat(Wd_ * B); }
  arma::mat mulWtW(const arma::mat& B) const { return sparse_ ? arma::mat(WtWsp_ * B) : arma::mat(WtWd_ * B); }

  // D x k matrix with rows idx = B (others zero)
  arma::mat scatter(const arma::uvec& idx, const arma::mat& B) const {
    arma::mat out(D_, B.n_cols, arma::fill::zeros);
    out.rows(idx) = B;
    return out;
  }
  arma::mat scatter(const arma::mat& B) const { return scatter(s_, B); }

 private:
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::SimplicialLDLT<SpMat, Eigen::Lower, Eigen::AMDOrdering<int> > SpLDLT;

  template <class Solver>
  static arma::mat solve_(const Solver& f, const arma::mat& B) {
    arma::mat X(B.n_rows, B.n_cols);
    Eigen::Map<const Eigen::MatrixXd> Be(B.memptr(), B.n_rows, B.n_cols);
    Eigen::Map<Eigen::MatrixXd> Xe(X.memptr(), X.n_rows, X.n_cols);
    Xe = f.solve(Be);
    return X;
  }

  template <class Vec>
  static double sum_log(const Vec& d, bool* pd) {
    double out = 0.0;
    *pd = true;
    for (Eigen::Index i = 0; i < d.size(); ++i) {
      if (!(d(i) > 0.0)) *pd = false;
      out += std::log(std::fabs(d(i)));
    }
    return out;
  }

  // Lower-triangular union pattern of I, W + W' and W'W with the three value
  // arrays aligned to it, so A(rho) and M are filled in O(nnz).
  void build_pattern() {
    std::vector<Eigen::Triplet<double> > tp;
    auto add_lower = [&](const arma::sp_mat& S, std::vector<Eigen::Triplet<double> >& t, bool vals) {
      for (arma::sp_mat::const_iterator it = S.begin(); it != S.end(); ++it)
        if (it.row() >= it.col()) t.emplace_back(it.row(), it.col(), vals ? double(*it) : 0.0);
    };
    const arma::sp_mat I = arma::speye(D_, D_);
    add_lower(I, tp, false);
    add_lower(Wssp_, tp, false);
    add_lower(WtWsp_, tp, false);
    auto build = [&](const arma::sp_mat& S) {
      std::vector<Eigen::Triplet<double> > t = tp;
      add_lower(S, t, true);
      SpMat Mx(D_, D_);
      Mx.setFromTriplets(t.begin(), t.end());
      Mx.makeCompressed();
      return Mx;
    };
    SpMat mI = build(I), mWs = build(Wssp_), mWtW = build(WtWsp_);
    vI_.assign(mI.valuePtr(), mI.valuePtr() + mI.nonZeros());
    vWs_.assign(mWs.valuePtr(), mWs.valuePtr() + mWs.nonZeros());
    vWtW_.assign(mWtW.valuePtr(), mWtW.valuePtr() + mWtW.nonZeros());
    Asp_ = mI;
    Msp_ = mI;
    diag_pos_.resize(D_);
    for (int j = 0; j < (int)D_; ++j)
      for (SpMat::InnerIterator it(mI, j); it; ++it)
        if (it.row() == j) diag_pos_[j] = &it.valueRef() - mI.valuePtr();
    ldA_.analyzePattern(Asp_);
    ldM_.analyzePattern(Msp_);
  }

  arma::uvec s_;
  arma::uword D_;
  bool sparse_ = true;
  double rho_ = 0.0, logdetA_ = 0.0, logdetM_ = 0.0;
  bool A_pd_ = false, M_pd_ = false;
  arma::sp_mat Wsp_, Wssp_, WtWsp_;
  SpMat Asp_, Msp_;
  std::vector<double> vI_, vWs_, vWtW_;
  std::vector<std::ptrdiff_t> diag_pos_;
  SpLDLT ldA_, ldM_;
  arma::mat Wd_, WtWd_;
  Eigen::LDLT<Eigen::MatrixXd> dA_, dM_;
};

#endif  // SAR_ENGINE_H
