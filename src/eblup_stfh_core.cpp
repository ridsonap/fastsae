// src/eblup_stfh_core.cpp
// Core Spatio-Temporal Fay-Herriot EBLUP estimation (Fisher-scoring, REML-style)
//
// The fit, the Fisher information, the SEs and the analytical MSE all use the
// structured representation of V (block-diagonal A plus the rank-D term
// Z Vu1 Z'; see "Structured evaluation" below): no M x M matrix is formed.
//
// The model is defined on the full D x T grid. Cells without a direct
// estimate (NA y: intermittent gaps or whole unsampled domains) are left out of
// the likelihood, which is then the marginal likelihood of the observed cells
// (the SAR effect u1 still lives on all D domains, the AR(1) effect u2 on all
// T periods), and they are predicted by the BLUP (see "All D x T cells" in
// .eblup_stfh_core()). With every cell observed, the arithmetic is the
// balanced-panel one.
//
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "eblup_stfh_core.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// ============================================================================
// build_omega2()
//
// Membangun matriks kovarians AR(1) "Omega2(rho2)" berukuran T x T beserta
// turunannya terhadap rho2.
// ============================================================================
static void build_omega2(int Tt, double rho2, arma::mat& Omega2, arma::mat& dOmega2) {
  Omega2.zeros(Tt, Tt);
  dOmega2.zeros(Tt, Tt);
  const double one_m_rho2sq = 1.0 - rho2 * rho2;

  for (int i = 0; i < Tt; ++i) {
    for (int j = 0; j < i; ++j) {
      const int lag = i - j;
      const double rho_lag = std::pow(rho2, (double) lag);
      const double val = rho_lag / one_m_rho2sq;
      Omega2(i, j) = val;
      Omega2(j, i) = val;

      const double raw_deriv = (double) lag * std::pow(rho2, (double)(lag - 1));
      const double dval = raw_deriv / one_m_rho2sq + (2.0 * rho2 / one_m_rho2sq) * val;
      dOmega2(i, j) = dval;
      dOmega2(j, i) = dval;
    }
  }
  // Diagonal (lag = 0): rho2^0 / (1-rho2^2) = 1/(1-rho2^2)
  Omega2.diag().fill(1.0 / one_m_rho2sq);
  // Turunan diagonal: d/drho2 [1/(1-rho2^2)] = 2*rho2 / (1-rho2^2)^2.
  dOmega2.diag().fill(2.0 * rho2 / (one_m_rho2sq * one_m_rho2sq));
}

// ============================================================================
// build_omega2_d2()
//
// Second derivative of the AR(1) matrix Omega2(rho2)_ij = rho2^l / (1 - rho2^2),
// l = |i - j|. With f = rho^l and h = 1 / (1 - rho^2) (h' = 2 rho h^2,
// h'' = 2 h^2 + 8 rho^2 h^3):  Omega2'' = f'' h + 2 f' h' + f h''.
// ============================================================================
static arma::mat build_omega2_d2(int Tt, double rho2) {
  arma::mat d2(Tt, Tt);
  const double h = 1.0 / (1.0 - rho2 * rho2);
  const double h1 = 2.0 * rho2 * h * h;
  const double h2 = 2.0 * h * h + 8.0 * rho2 * rho2 * h * h * h;
  for (int i = 0; i < Tt; ++i) {
    for (int j = 0; j < Tt; ++j) {
      const int l = std::abs(i - j);
      const double f = std::pow(rho2, (double) l);
      const double f1 = (l >= 1) ? l * std::pow(rho2, (double)(l - 1)) : 0.0;
      const double f2 = (l >= 2) ? l * (l - 1) * std::pow(rho2, (double)(l - 2)) : 0.0;
      d2(i, j) = f2 * h + 2.0 * f1 * h1 + f * h2;
    }
  }
  return d2;
}

// ============================================================================
// st_make_index(): observed cells of the D x T grid -> per-domain layout of
// the compact (observed-only) rows.
// ============================================================================
StIndex st_make_index(int D, int T, const arma::uvec& observed) {
  StIndex ix;
  ix.D = D; ix.T = T;
  ix.off.assign(D + 1, 0);
  ix.tix.resize(D);
  std::vector<arma::uword> pos;
  for (int d = 0; d < D; ++d) {
    std::vector<arma::uword> ts;
    for (int t = 0; t < T; ++t) {
      if (observed(d * T + t)) { ts.push_back(t); pos.push_back(d * T + t); }
    }
    ix.tix[d] = arma::conv_to<arma::uvec>::from(ts);
    ix.off[d + 1] = ix.off[d] + (int) ts.size();
  }
  ix.n = ix.off[D];
  ix.full = (ix.n == D * T);
  ix.obs = arma::conv_to<arma::uvec>::from(pos);
  return ix;
}

// ============================================================================
// Structured evaluation (no M x M matrices)
//
// V = A + Z Vu1 Z', A = blockdiag(A_d) (n_d x n_d per domain: the observed
// periods; empty for an unsampled domain), Z = the domain indicator of each row.
// With Woodbury, V^{-1} = A^{-1} - (A^{-1}Z) C^{-1} (A^{-1}Z)',
// C = Vu1^{-1} + Z'A^{-1}Z, and
//   P = V^{-1} - V^{-1}X Q X'V^{-1} = A^{-1} - U Wm U',
//   U = [A^{-1}Z, V^{-1}X] (M x (D + p)), Wm = blockdiag(C^{-1}, Q).
// The derivatives are of two kinds: V_a = Z B_a Z' (sigma21, rho1) and
// V_a = I_D (x) E_a (sigma22, rho2). Scores and Fisher information then only
// need D x D, (D+p) x (D+p) and T x T blocks:
//   tr(P Z Ba Z' P Z Bb Z') = tr(Ba Z'PZ Bb Z'PZ)
//   tr(P Z Ba Z' P (I(x)E)) = tr(Ba (PZ)'(I(x)E)(PZ))
//   tr(P Ej P Ek) = sum_d tr(A_d^{-1} Ej A_d^{-1} Ek)
//                   - 2 tr(Wm (Ek U)'(A^{-1}Ej U)) + tr(Wm U'EjU Wm U'EkU)
// i.e. O(D^3 + M T (D+p)) instead of O(M^3) per iteration. Algebraically
// identical to the dense computation. Pure Armadillo (thread-safe).
// ============================================================================
namespace {
struct StPieces {
  int D = 0, T = 0, M = 0;
  bool isST = true, hasVu1 = true, full = true;
  std::vector<int> off;          // compact rows of domain d: off[d] .. off[d+1]-1
  std::vector<arma::uvec> tix;   // observed periods of domain d
  std::vector<arma::mat> Ainv;   // D blocks, n_d x n_d (empty: unsampled domain)
  arma::mat AZ;                  // M x D, A^{-1} Z
  arma::vec c;                   // D, 1'A_d^{-1}1
  arma::mat Cmat, Cinv;          // D x D (when sigma21 != 0)
  arma::mat Omega1, Omega2, dOmega2;
  double logdetA = 0.0;
};
// Compact rows f..l of domain d; false for a domain without observed cells.
inline bool st_rng(const StPieces& p, int d, int& f, int& l) {
  f = p.off[d];
  l = p.off[d + 1] - 1;
  return l >= f;
}
// The block of a T x T matrix E for the observed periods of domain d.
inline arma::mat st_sub(const StPieces& p, const arma::mat& E, int d) {
  return p.full ? E : arma::mat(E.submat(p.tix[d], p.tix[d]));
}
}  // namespace

static bool st_build(const arma::vec& th, bool isST, const StIndex& ix,
                     const arma::vec& vardir, const arma::mat& W, StPieces& p) {
  const int D = ix.D, T = ix.T;
  p.D = D; p.T = T; p.M = ix.n; p.isST = isST;
  p.full = ix.full; p.off = ix.off; p.tix = ix.tix;
  const arma::mat I = arma::eye<arma::mat>(D, D);
  const arma::mat ImrW = I - th(1) * W;
  const arma::mat A1 = arma::symmatu(ImrW.t() * ImrW);
  if (!arma::inv_sympd(p.Omega1, A1) && !arma::inv(p.Omega1, A1)) return false;
  if (isST) build_omega2(T, th(3), p.Omega2, p.dOmega2);
  p.Ainv.resize(D);
  p.AZ.zeros(p.M, D);
  p.c.set_size(D);
  p.logdetA = 0.0;
  for (int d = 0; d < D; ++d) {
    int f, l;
    if (!st_rng(p, d, f, l)) { p.c(d) = 0.0; p.Ainv[d].reset(); continue; }
    const arma::vec psi_d = vardir.subvec(f, l);
    arma::mat Ad;
    if (isST) {
      if (p.full) Ad = arma::symmatu(th(2) * p.Omega2 + arma::diagmat(psi_d));
      else Ad = arma::symmatu(th(2) * st_sub(p, p.Omega2, d) + arma::diagmat(psi_d));
    } else {
      const arma::vec bd = th(2) + psi_d;
      if (arma::any(bd <= 0.0)) return false;
      Ad = arma::diagmat(bd);
    }
    if (!arma::inv_sympd(p.Ainv[d], Ad) && !arma::inv(p.Ainv[d], Ad)) return false;
    double ld, sg;
    if (!arma::log_det(ld, sg, Ad) || sg <= 0) return false;
    p.logdetA += ld;
    const arma::vec a = arma::sum(p.Ainv[d], 1);
    p.AZ.submat(f, d, l, d) = a;
    p.c(d) = arma::accu(a);
  }
  p.hasVu1 = (th(0) != 0.0);
  if (p.hasVu1) {
    // V_u1^{-1} = (sigma21 Omega1)^{-1} = A1 / sigma21 (no second inversion)
    p.Cmat = arma::symmatu(A1 / th(0) + arma::diagmat(p.c));
    if (!arma::inv_sympd(p.Cinv, p.Cmat) && !arma::inv(p.Cinv, p.Cmat)) return false;
  }
  return true;
}

// A^{-1} Y, (I (x) E) Y, blockdiag(A_d^{-1} E) Y, Z'Y, V^{-1} Y
static arma::mat st_Ainv(const StPieces& p, const arma::mat& Y) {
  arma::mat R(Y.n_rows, Y.n_cols);
  int f, l;
  for (int d = 0; d < p.D; ++d) {
    if (!st_rng(p, d, f, l)) continue;
    R.rows(f, l) = p.Ainv[d] * Y.rows(f, l);
  }
  return R;
}
// blockdiag(E[tix_d, tix_d]) Y for a T x T matrix E
static arma::mat st_blk(const StPieces& p, const arma::mat& E, const arma::mat& Y) {
  arma::mat R(Y.n_rows, Y.n_cols);
  int f, l;
  for (int d = 0; d < p.D; ++d) {
    if (!st_rng(p, d, f, l)) continue;
    if (p.full) R.rows(f, l) = E * Y.rows(f, l);
    else R.rows(f, l) = st_sub(p, E, d) * Y.rows(f, l);
  }
  return R;
}
static arma::mat st_Zt(const StPieces& p, const arma::mat& Y) {
  arma::mat R(p.D, Y.n_cols, fill::zeros);
  int f, l;
  for (int d = 0; d < p.D; ++d) {
    if (!st_rng(p, d, f, l)) continue;
    R.row(d) = arma::sum(Y.rows(f, l), 0);
  }
  return R;
}
// Z S: row i of the result is row d(i) of S (S has D rows)
static arma::mat st_Z(const StPieces& p, const arma::mat& S) {
  arma::mat R(p.M, S.n_cols);
  int f, l;
  for (int d = 0; d < p.D; ++d) {
    if (!st_rng(p, d, f, l)) continue;
    R.rows(f, l) = arma::repmat(S.row(d), l - f + 1, 1);
  }
  return R;
}
static arma::mat st_Vinv(const StPieces& p, const arma::mat& Y) {
  arma::mat R = st_Ainv(p, Y);
  if (p.hasVu1) R -= p.AZ * (p.Cinv * (p.AZ.t() * Y));
  return R;
}

// Low-rank factor U = [A^-1 Z, V^-1 X] of P = A^-1 - U Wm U' (Wm = diag(C^-1, Q)).
// The A^-1 Z part has a single non-zero per row (a_i in column d(i)), so
//   U'Y = [Z' A^-1 Y ; X' V^-1 Y]  and  U S  cost O(M k) instead of O(M r k).
namespace {
struct StLR {
  const StPieces& p;
  const arma::mat& ViX;
  arma::vec avec;
  int off, px, r;
  StLR(const StPieces& p_, const arma::mat& ViX_) : p(p_), ViX(ViX_) {
    px = ViX.n_cols;
    off = p.hasVu1 ? p.D : 0;
    r = off + px;
    avec.zeros(p.M);
    if (p.hasVu1)
      for (int d = 0; d < p.D; ++d) {
        int f, l;
        if (!st_rng(p, d, f, l)) continue;
        avec.subvec(f, l) = p.AZ.submat(f, d, l, d);
      }
  }
  // U' Y (r x k)
  arma::mat Ut(const arma::mat& Y) const {
    arma::mat R(r, Y.n_cols);
    if (p.hasVu1) R.rows(0, p.D - 1) = st_Zt(p, st_Ainv(p, Y));
    R.rows(off, off + px - 1) = ViX.t() * Y;
    return R;
  }
  // U S (M x k)
  arma::mat Umul(const arma::mat& Sm) const {
    arma::mat R = ViX * Sm.rows(off, off + px - 1);
    if (p.hasVu1)
      for (int d = 0; d < p.D; ++d) {
        int f, l;
        if (!st_rng(p, d, f, l)) continue;
        R.rows(f, l) += avec.subvec(f, l) * Sm.row(d);
      }
    return R;
  }
  // rowsum(Y % U) for Y (M x r)
  arma::vec Urow(const arma::mat& Y) const {
    arma::vec out = arma::sum(Y.cols(off, off + px - 1) % ViX, 1);
    if (p.hasVu1)
      for (int d = 0; d < p.D; ++d) {
        int f, l;
        if (!st_rng(p, d, f, l)) continue;
        out.subvec(f, l) += avec.subvec(f, l) % Y.submat(f, d, l, d);
      }
    return out;
  }
  // U as a dense M x r matrix
  arma::mat dense() const {
    arma::mat Ud(p.M, r, fill::zeros);
    Ud.cols(off, off + px - 1) = ViX;
    if (p.hasVu1)
      for (int d = 0; d < p.D; ++d) {
        int f, l;
        if (!st_rng(p, d, f, l)) continue;
        Ud.submat(f, d, l, d) = avec.subvec(f, l);
      }
    return Ud;
  }
  arma::mat Wm(const arma::mat& Q) const {
    arma::mat W_(r, r, fill::zeros);
    if (p.hasVu1) W_.submat(0, 0, p.D - 1, p.D - 1) = p.Cinv;
    W_.submat(off, off, off + px - 1, off + px - 1) = Q;
    return W_;
  }
};
}  // namespace

// An M x r factor F = [a-part | Xp] with the same layout as U: column d of the
// a-part is non-zero only in block d (entries a_i), Xp is dense M x px.
// Block-diagonal transforms keep this layout, so P E U, A^-1 X U, ... stay cheap.
namespace {
struct BSF {
  int D = 0, T = 0, M = 0, off = 0, px = 0, r = 0;
  bool hasA = false;
  std::vector<int> roff;          // compact rows of each domain (as StPieces::off)
  arma::vec a;
  arma::mat Xp;
  // F S  (M x k)
  arma::mat mul(const arma::mat& Sm) const {
    arma::mat R = Xp * Sm.rows(off, off + px - 1);
    if (hasA)
      for (int d = 0; d < D; ++d) {
        const int f = roff[d], l = roff[d + 1] - 1;
        if (l < f) continue;
        R.rows(f, l) += a.subvec(f, l) * Sm.row(d);
      }
    return R;
  }
  // rowsum(Y % F) for Y (M x r)
  arma::vec row(const arma::mat& Y) const {
    arma::vec out = arma::sum(Y.cols(off, off + px - 1) % Xp, 1);
    if (hasA)
      for (int d = 0; d < D; ++d) {
        const int f = roff[d], l = roff[d + 1] - 1;
        if (l < f) continue;
        out.subvec(f, l) += a.subvec(f, l) % Y.submat(f, d, l, d);
      }
    return out;
  }
};
BSF bsf_of(const StLR& L) {
  BSF F;
  F.D = L.p.D; F.T = L.p.T; F.M = L.p.M; F.off = L.off; F.px = L.px; F.r = L.r;
  F.hasA = L.p.hasVu1; F.a = L.avec; F.Xp = L.ViX; F.roff = L.p.off;
  return F;
}
// blockdiag(Xb) F
BSF bsf_transform(const BSF& F, const std::vector<arma::mat>& Xb) {
  BSF G = F;
  for (int d = 0; d < F.D; ++d) {
    const int f = F.roff[d], l = F.roff[d + 1] - 1;
    if (l < f) continue;
    if (F.hasA) G.a.subvec(f, l) = Xb[d] * F.a.subvec(f, l);
    G.Xp.rows(f, l) = Xb[d] * F.Xp.rows(f, l);
  }
  return G;
}
}  // namespace

// REML score S and Fisher information F at th (pieces p already built);
// also returns ViX, Q and Py. th may have negative variances (classical
// iteration): no logarithms are taken here. Cost O(D^3 + M T r), no M x M
// or M x r x r products.
static bool st_score_info(const StPieces& p, const arma::vec& th, const arma::mat& X,
                          const arma::vec& y, const arma::mat& W, const arma::mat& WtW,
                          arma::vec& S, arma::mat& F,
                          arma::mat& ViX, arma::mat& Q, arma::vec& Py) {
  const int D = p.D, T = p.T, nparam = p.isST ? 4 : 3;
  ViX = st_Vinv(p, X);
  const arma::mat XtViX = arma::symmatu(X.t() * ViX);
  if (!arma::inv_sympd(Q, XtViX) && !arma::inv(Q, XtViX)) return false;
  const arma::vec Viy = st_Vinv(p, arma::mat(y));
  Py = Viy - ViX * (Q * (X.t() * Viy));

  const StLR L(p, ViX);
  const arma::mat Wm = L.Wm(Q);
  const arma::mat Ud = L.dense();

  // P Z = A^-1 Z - U G with G = Wm U'Z; U'Z = [diag(c); X'V^-1 Z]. The top
  // (A^-1 Z) blocks of U'Z, U'E U and U'E A^-1 Z are diagonal, so the products
  // with them below are column scalings (rmul / tmul), not D^3 products.
  const int off = L.off, px = L.px;
  const arma::mat ZtViX = st_Zt(p, ViX);                      // D x px
  arma::mat UtZ(L.r, D);
  if (p.hasVu1) UtZ.rows(0, D - 1) = arma::diagmat(p.c);
  UtZ.rows(off, off + px - 1) = ZtViX.t();
  // M_ (r x r) times Y (r x k), where M_ = [diag(dg) R12 ; R21 R22] (hasVu1)
  auto tmul = [&](const arma::vec& dg, const arma::mat& Mfull, const arma::mat& Y) -> arma::mat {
    if (!p.hasVu1) return Mfull * Y;
    arma::mat R(L.r, Y.n_cols);
    const arma::mat Yt = Y.rows(0, D - 1), Yb = Y.rows(off, off + px - 1);
    R.rows(0, D - 1) = (Yt.each_col() % dg) + Mfull.submat(0, off, D - 1, off + px - 1) * Yb;
    R.rows(off, off + px - 1) = Mfull.submat(off, 0, off + px - 1, D - 1) * Yt +
      Mfull.submat(off, off, off + px - 1, off + px - 1) * Yb;
    return R;
  };
  // Wm Y with Wm = diag(C^-1, Q)
  auto wmul = [&](const arma::mat& Y) -> arma::mat {
    arma::mat R(L.r, Y.n_cols);
    if (p.hasVu1) R.rows(0, D - 1) = p.Cinv * Y.rows(0, D - 1);
    R.rows(off, off + px - 1) = Q * Y.rows(off, off + px - 1);
    return R;
  };
  arma::mat G(L.r, D);                                        // Wm U'Z
  if (p.hasVu1) G.rows(0, D - 1) = p.Cinv.each_row() % p.c.t();
  G.rows(off, off + px - 1) = Q * ZtViX.t();
  // Z'PZ = diag(c) - (U'Z)' G
  arma::mat ZPZ = -(ZtViX * G.rows(off, off + px - 1));
  if (p.hasVu1) ZPZ -= G.rows(0, D - 1).each_col() % p.c;
  ZPZ.diag() += p.c;
  ZPZ = 0.5 * (ZPZ + ZPZ.t());
  const arma::vec zPy = st_Zt(p, arma::mat(Py));

  // B_rho1 = -sigma21 Omega1 C Omega1 with the sparse C = 2 rho W'W - W - W'
  const arma::sp_mat Cw(arma::mat(2.0 * th(1) * WtW - (W + W.t())));
  const arma::mat B[2] = { p.Omega1, -th(0) * (p.Omega1 * arma::mat(Cw * p.Omega1)) };
  const arma::mat BZ[2] = { B[0] * ZPZ, B[1] * ZPZ };
  std::vector<arma::mat> E;
  if (p.isST) { E.push_back(p.Omega2); E.push_back(th(2) * p.dOmega2); }
  else E.push_back(arma::eye<arma::mat>(T, T));
  const int nE = (int) E.size();

  S.set_size(nparam);
  F.set_size(nparam, nparam);
  for (int a = 0; a < 2; ++a) {
    S(a) = -0.5 * arma::trace(BZ[a]) + 0.5 * arma::dot(zPy, B[a] * zPy);
    for (int b = a; b < 2; ++b) {
      F(a, b) = 0.5 * arma::accu(BZ[a] % BZ[b].t());
      F(b, a) = F(a, b);
    }
  }
  std::vector<arma::mat> EU(nE), UtEU(nE);
  std::vector<arma::vec> tdgE(nE);
  for (int j = 0; j < nE; ++j) {
    EU[j] = st_blk(p, E[j], Ud);
    UtEU[j] = L.Ut(EU[j]);                                     // U' E_j U
    double trAE = 0.0;
    int f, l;
    for (int d = 0; d < D; ++d)
      if (!st_rng(p, d, f, l)) continue;
      else if (p.full) trAE += arma::accu(p.Ainv[d] % E[j].t());
      else trAE += arma::accu(p.Ainv[d] % st_sub(p, E[j], d).t());
    const arma::vec EPy = st_blk(p, E[j], arma::mat(Py));
    S(2 + j) = -0.5 * (trAE - arma::accu(Wm % UtEU[j].t())) + 0.5 * arma::dot(Py, EPy);
    // Z'P E_j P Z = (AZ - U G)' E_j (AZ - U G)
    //             = Z'A^-1 E A^-1 Z - H'G - G'H + G' (U'E U) G,  H = U' E A^-1 Z
    const arma::mat EAZ = st_blk(p, E[j], p.AZ);
    const arma::vec tdg = arma::diagvec(st_Zt(p, st_Ainv(p, EAZ)));   // diagonal
    const arma::mat Hb = ViX.t() * EAZ;                          // bottom of H (px x D)
    arma::mat HtG = Hb.t() * G.rows(off, off + px - 1);
    if (p.hasVu1) HtG += G.rows(0, D - 1).each_col() % tdg;
    arma::mat PZtEPZ = -(HtG + HtG.t()) + G.t() * tmul(tdg, UtEU[j], G);
    PZtEPZ.diag() += tdg;
    tdgE[j] = tdg;
    for (int a = 0; a < 2; ++a) {
      F(a, 2 + j) = 0.5 * arma::accu(B[a] % PZtEPZ.t());
      F(2 + j, a) = F(a, 2 + j);
    }
  }
  // Wm (U'E_j U): the top-left block of U'E_j U is diagonal (tdgE[j])
  auto wmulE = [&](int j) -> arma::mat {
    if (!p.hasVu1) return wmul(UtEU[j]);
    const arma::mat& M_ = UtEU[j];
    arma::mat R(L.r, L.r);
    R.submat(0, 0, D - 1, D - 1) = p.Cinv.each_row() % tdgE[j].t();
    R.submat(0, off, D - 1, off + px - 1) = p.Cinv * M_.submat(0, off, D - 1, off + px - 1);
    R.rows(off, off + px - 1) = Q * M_.rows(off, off + px - 1);
    return R;
  };
  for (int j = 0; j < nE; ++j) {
    for (int k = j; k < nE; ++k) {
      double T1 = 0.0;
      int f, l;
      for (int d = 0; d < D; ++d) {
        if (!st_rng(p, d, f, l)) continue;
        if (p.full) {
          const arma::mat AEj = p.Ainv[d] * E[j], AEk = p.Ainv[d] * E[k];
          T1 += arma::accu(AEj % AEk.t());
        } else {
          const arma::mat Ejd = st_sub(p, E[j], d), Ekd = st_sub(p, E[k], d);
          const arma::mat AEj = p.Ainv[d] * Ejd, AEk = p.Ainv[d] * Ekd;
          T1 += arma::accu(AEj % AEk.t());
        }
      }
      // U' E_k A^-1 E_j U
      const arma::mat UEAEU = L.Ut(st_blk(p, E[k], st_Ainv(p, EU[j])));
      const double T2 = arma::accu(Wm % UEAEU.t());
      const arma::mat WEj = wmulE(j), WEk = wmulE(k);
      const double T3 = arma::accu(WEj % WEk.t());
      F(2 + j, 2 + k) = 0.5 * (T1 - 2.0 * T2 + T3);
      F(2 + k, 2 + j) = F(2 + j, 2 + k);
    }
  }
  return S.is_finite() && F.is_finite();
}

// ============================================================================
// stfh_reml_eval()
//
// Restricted log-likelihood (up to a constant)
//   ll = -1/2 [ log|V| + log|X'V^{-1}X| + y'Py ],
// with log|V| = log|A| + log|Vu1| + log|Vu1^{-1} + Z1'A^{-1}Z1| (determinant
// lemma), and optionally the REML score S and Fisher information F.
// Pure Armadillo (thread-safe).
// ============================================================================
static bool stfh_reml_eval(
    const arma::vec& th, bool isST, const StIndex& ix,
    const arma::mat& X, const arma::vec& y, const arma::vec& vardir,
    const arma::mat& W, bool want_SF,
    double& ll, arma::vec& S, arma::mat& F
) {
  StPieces p;
  if (!st_build(th, isST, ix, vardir, W, p)) return false;
  arma::mat ViX, Q;
  arma::vec Py;
  arma::vec S0;
  arma::mat F0;
  if (want_SF) {
    if (!st_score_info(p, th, X, y, W, W.t() * W, S0, F0, ViX, Q, Py)) return false;
  } else {
    // log-likelihood only (line search): no score / information
    ViX = st_Vinv(p, X);
    const arma::mat XtViX = arma::symmatu(X.t() * ViX);
    if (!arma::inv_sympd(Q, XtViX)) return false;
    const arma::vec Viy = st_Vinv(p, arma::mat(y));
    Py = Viy - ViX * (Q * (X.t() * Viy));
  }

  double logdetV = p.logdetA, ld, sg;
  if (p.hasVu1) {
    if (th(0) <= 0.0) return false;
    if (!arma::log_det(ld, sg, p.Omega1) || sg <= 0) return false;
    logdetV += ix.D * std::log(th(0)) + ld;
    if (!arma::log_det(ld, sg, p.Cmat) || sg <= 0) return false;
    logdetV += ld;
  }
  if (!arma::log_det(ld, sg, arma::mat(arma::symmatu(X.t() * ViX))) || sg <= 0) return false;
  ll = -0.5 * (logdetV + ld + arma::dot(y, Py));
  if (!std::isfinite(ll)) return false;
  if (want_SF) { S = S0; F = F0; }
  return true;
}

// ============================================================================
// stfh_fit_safeguarded()
//
// Fallback when the classical (sae-compatible) Fisher scoring fails, which
// happens when a variance iterate goes negative: with rho1 -> 1 the matrix V
// then becomes strongly indefinite and the iteration wanders. Here instead:
//   * projected Fisher scoring on sigma21, sigma22 >= 0, |rho| <= 0.999;
//   * parameters at an active bound are held fixed (a variance at 0 whose
//     score points outwards, and the correlation of a variance that is 0,
//     which is then unidentified);
//   * step halving until the REML log-likelihood does not decrease.
// Returns the constrained REML maximiser in `th`.
// ============================================================================
static bool stfh_fit_safeguarded(
    const arma::mat& X, const arma::vec& y, const arma::vec& vardir,
    const arma::mat& W, const StIndex& ix, bool isST,
    int maxiter, double precision, arma::vec& th, int& iters
) {
  const int nparam = isST ? 4 : 3;
  auto project = [&](arma::vec t) {
    t(0) = std::max(t(0), 0.0);
    t(2) = std::max(t(2), 0.0);
    t(1) = std::min(std::max(t(1), -0.999), 0.999);
    if (isST) t(3) = std::min(std::max(t(3), -0.999), 0.999);
    return t;
  };
  th = project(th);
  double ll;
  arma::vec S;
  arma::mat F;
  if (!stfh_reml_eval(th, isST, ix, X, y, vardir, W, true, ll, S, F)) return false;

  // Near a nearly unidentified ridge (e.g. sigma22 -> 0 with rho2 -> 1, where
  // only sigma22 / (1 - rho2^2) is well determined) progress per iteration is
  // slow, so the fallback gets its own, larger budget.
  const int max_it = 5 * maxiter;
  for (iters = 1; iters <= max_it; ++iters) {
    // Active set: a variance at 0 whose score points outwards; the
    // correlation of a zero variance (unidentified); a correlation at its
    // bound whose score points outwards.
    std::vector<arma::uword> fr;
    for (int a = 0; a < nparam; ++a) {
      bool fixed = false;
      if ((a == 0 || a == 2) && th(a) == 0.0 && S(a) <= 0.0) fixed = true;
      if ((a == 1 && th(0) == 0.0) || (a == 3 && th(2) == 0.0)) fixed = true;
      if ((a == 1 || a == 3) && ((th(a) >= 0.999 && S(a) >= 0.0) ||
                                 (th(a) <= -0.999 && S(a) <= 0.0))) fixed = true;
      if (!fixed) fr.push_back(a);
    }
    if (fr.empty()) return true;
    const arma::uvec f = arma::conv_to<arma::uvec>::from(fr);

    // Converged when every free score is negligible in standard-error units:
    // |S_a| / sqrt(F_aa) is the size of a scoring step measured in SEs (the
    // remaining log-likelihood gain is ~ z^2 / 2). A relative change in the
    // parameter is a poor criterion for a variance close to 0, where Fisher
    // scoring also zig-zags because the observed curvature differs from F.
    // Default precision 1e-4 => steps below 0.1% of an SE.
    {
      double zmax = 0.0;
      for (arma::uword a : fr) {
        if (F(a, a) > 0.0) zmax = std::max(zmax, std::fabs(S(a)) / std::sqrt(F(a, a)));
      }
      if (zmax < 10.0 * precision) return true;
    }

    // Candidate directions: Fisher-scoring step on the free set, then (if the
    // line search fails on it) the diagonally scaled gradient.
    std::vector<arma::vec> dirs;
    {
      arma::vec sf, step(nparam, fill::zeros);
      // no_approx: a singular F block yields no Newton direction (the scaled
      // gradient below is used) instead of a printed warning from a thread.
      if (arma::solve(sf, arma::mat(F.submat(f, f)), arma::vec(S.elem(f)),
                      arma::solve_opts::no_approx)) {
        step.elem(f) = sf;
        dirs.push_back(step);
      }
      arma::vec g(nparam, fill::zeros);
      const arma::vec dF = F.diag();
      for (arma::uword a : fr) g(a) = S(a) / std::max(dF(a), 1e-12);
      dirs.push_back(g);
    }

    double t = 1.0, llc = 0.0;
    arma::vec cand, Sc;
    arma::mat Fc;
    bool accepted = false;
    // line search on the log-likelihood only; score / information are
    // computed once, for the accepted point
    for (const arma::vec& step : dirs) {
      t = 1.0;
      for (int h = 0; h < 20; ++h) {
        cand = project(th + t * step);
        if (stfh_reml_eval(cand, isST, ix, X, y, vardir, W, false, llc, Sc, Fc) &&
            llc > ll) {
          accepted = true;
          break;
        }
        t *= 0.5;
      }
      if (accepted) break;
    }
    if (!accepted) return true;  // no ascent direction left: constrained maximum
    if (!stfh_reml_eval(cand, isST, ix, X, y, vardir, W, true, llc, Sc, Fc)) return false;
    const bool full_step = (t == 1.0);

    const double rel = arma::max(arma::abs(cand - th) / arma::max(arma::abs(th), 1e-4 * arma::ones<arma::vec>(nparam)));
    const double gain = llc - ll;
    th = cand; ll = llc; S = Sc; F = Fc;
    if (rel < precision) return true;
    // On a (nearly) flat ridge -- e.g. rho1 when sigma21 ~ 0 -- the parameter
    // can keep moving while the likelihood no longer changes. Only trusted
    // for a full (un-halved) scoring step, so that a step shortened by the
    // line search is never mistaken for convergence.
    if (full_step && iters > 1 && gain <= 1e-10 * (1.0 + std::fabs(ll))) return true;
  }
  return false;
}

// ============================================================================
// st_full_u(): random effects on the full D x T grid from V^-1 (y - X beta)
//   u1 = sigma21 Omega1 Z'V^-1 r  (all D domains: kriging for unsampled ones)
//   u2_d = sigma22 Omega2[., tix_d] (V^-1 r)_d  (all T periods of domain d;
//          model "S": sigma22 (V^-1 r)_i at the observed cells, 0 elsewhere)
// Returned as two D*T vectors (domain-major). At observed cells these equal the
// EBLUP's random effects; elsewhere they are the BLUP of the unobserved effects.
// ============================================================================
static void st_full_u(const StPieces& p, double s21, double s22,
                      const arma::vec& invVresid, arma::vec& u1_full, arma::vec& u2_full) {
  const int D = p.D, T = p.T;
  arma::vec zr(D, fill::zeros);
  {
    int f, l;
    for (int d = 0; d < D; ++d) if (st_rng(p, d, f, l)) zr(d) = arma::sum(invVresid.subvec(f, l));
  }
  const arma::mat Vu1 = s21 * p.Omega1;                  // stored first, as before
  const arma::vec u1 = Vu1 * zr;
  u1_full.set_size(D * T);
  u2_full.zeros(D * T);
  int f, l;
  for (int d = 0; d < D; ++d) {
    u1_full.subvec(d * T, d * T + T - 1).fill(u1(d));
    if (!st_rng(p, d, f, l)) continue;
    const arma::vec r_d = invVresid.subvec(f, l);
    if (p.isST) {
      if (p.full) u2_full.subvec(d * T, d * T + T - 1) = s22 * (p.Omega2 * r_d);
      else u2_full.subvec(d * T, d * T + T - 1) = s22 * (arma::mat(p.Omega2.cols(p.tix[d])) * r_d);
    } else {
      for (arma::uword k = 0; k < p.tix[d].n_elem; ++k)
        u2_full(d * T + p.tix[d](k)) = s22 * r_d(k);
    }
  }
}

// ============================================================================
// stfh_fit_arma()
//
// Pure Armadillo REML Fisher-scoring fit + EBLUP (no R objects, no R RNG),
// safe to call from an OpenMP parallel region. Same algorithm as the
// iteration in .eblup_stfh_core(): same update, clamping, convergence rule
// and final truncation. Used to RE-ESTIMATE theta in every bootstrap
// replicate (otherwise the bootstrap MSE misses the g3 contribution).
// ============================================================================
STFitArma stfh_fit_arma(
    const arma::mat& X, const arma::vec& y, const arma::vec& vardir,
    const arma::mat& W, const StIndex& ix, bool isST,
    int maxiter, double precision,
    double sigma21_start, double rho1_start,
    double sigma22_start, double rho2_start,
    const arma::mat& Xfull
) {
  STFitArma res;
  res.ok = false;
  res.converged = false;
  res.n_iter = 0;
  res.sigma21 = res.rho1 = res.sigma22 = res.rho2 = 0.0;

  const int nparam = isST ? 4 : 3;
  const arma::mat WtW = W.t() * W;

  arma::vec thetak(nparam), thetak1(nparam);
  thetak1(0) = sigma21_start;
  thetak1(1) = rho1_start;
  thetak1(2) = sigma22_start;
  if (isST) thetak1(3) = rho2_start;

  int k = 0;
  double diff = precision + 1.0;
  bool primary_failed = false;
  while (diff > precision && k < maxiter) {
    ++k;
    thetak = thetak1;

    StPieces sp;
    if (!st_build(thetak, isST, ix, vardir, W, sp)) { primary_failed = true; break; }
    arma::vec S;
    arma::mat F, ViX_k, Q_k;
    arma::vec Py_k;
    if (!st_score_info(sp, thetak, X, y, W, WtW, S, F, ViX_k, Q_k, Py_k)) { primary_failed = true; break; }
    arma::mat Finv;
    if (!arma::inv_sympd(Finv, F)) {
      if (!arma::inv(Finv, F)) { primary_failed = true; break; }
    }
    thetak1 = thetak + Finv * S;
    if (!thetak1.is_finite()) { primary_failed = true; break; }

    if (thetak1(1) <= -1) thetak1(1) = -0.999;
    else if (thetak1(1) >= 1) thetak1(1) = 0.999;
    if (isST) {
      if (thetak1(3) <= -1) thetak1(3) = -0.999;
      else if (thetak1(3) >= 1) thetak1(3) = 0.999;
    }

    arma::vec thetak_safe = thetak;
    for (int i = 0; i < nparam; ++i) if (thetak_safe(i) == 0.0) thetak_safe(i) = 1e-4;
    diff = arma::max(arma::abs((thetak_safe - thetak1) / thetak_safe));
  }
  if (primary_failed || diff > precision) {
    // classical iteration failed: constrained REML with step halving
    arma::vec th(nparam);
    th(0) = sigma21_start; th(1) = rho1_start; th(2) = sigma22_start;
    if (isST) th(3) = rho2_start;
    int it2 = 0;
    if (!stfh_fit_safeguarded(X, y, vardir, W, ix, isST, maxiter, precision, th, it2)) {
      res.n_iter = k + it2;
      return res;
    }
    thetak1 = th;
    k += it2;
  } else if (thetak1(0) < 0.0 || thetak1(2) < 0.0) {
    // same constrained re-maximisation as in .eblup_stfh_core()
    arma::vec th = thetak1;
    int it2 = 0;
    if (!stfh_fit_safeguarded(X, y, vardir, W, ix, isST, maxiter, precision, th, it2)) {
      res.n_iter = k + it2;
      return res;
    }
    thetak1 = th;
    k += it2;
  }
  res.n_iter = k;

  const double s21 = std::max(thetak1(0), 0.0);
  const double r1 = thetak1(1);
  const double s22 = std::max(thetak1(2), 0.0);
  const double r2 = isST ? thetak1(3) : 0.0;

  arma::vec thf(nparam);
  thf(0) = s21; thf(1) = r1; thf(2) = s22;
  if (isST) thf(3) = r2;
  StPieces spf;
  if (!st_build(thf, isST, ix, vardir, W, spf)) return res;
  const arma::mat ViXf = st_Vinv(spf, X);
  arma::mat Q;
  if (!arma::inv_sympd(Q, arma::symmatu(X.t() * ViXf))) return res;
  res.beta = Q * (ViXf.t() * y);
  const arma::vec invVresid = st_Vinv(spf, arma::mat(y - X * res.beta));

  // EBLUP of every D x T cell: x' beta + u1_d + u2_dt (BLUP where unobserved)
  arma::vec u1_full, u2_full;
  st_full_u(spf, s21, s22, invVresid, u1_full, u2_full);
  res.theta = Xfull * res.beta;
  // same association as the balanced-panel code: X beta + (u1 + u2)
  const arma::vec eblup = res.theta + (u1_full + u2_full);

  res.eblup = eblup;
  res.sigma21 = s21;
  res.rho1 = r1;
  res.sigma22 = s22;
  res.rho2 = r2;
  res.ok = eblup.is_finite();
  res.converged = res.ok;
  return res;
}

// ============================================================================
// .eblup_stfh_core()
//
// Output structure matching seblup_area + estvarcomp for spatio-temporal
// ============================================================================
// [[Rcpp::export(.eblup_stfh_core)]]
List eblup_stfh_core(
    const arma::mat& Xall,
    const arma::vec& yall,
    const arma::vec& vardirall,
    const arma::mat& proxmat,
    int D,
    int Tt,
    std::string model = "ST",
    int maxiter = 100,
    double precision = 1e-4,
    double sigma21_start = -1.0,
    double rho1_start = 0.5,
    double sigma22_start = -1.0,
    double rho2_start = 0.5
) {
  if (model != "S" && model != "ST") Rcpp::stop("Argument model must be \"S\" or \"ST\".");

  const int Mfull = D * Tt;
  const int p = Xall.n_cols;
  if ((int) Xall.n_rows != Mfull) Rcpp::stop("nrow(Xall) must equal D*T.");
  if ((int) yall.n_elem != Mfull) Rcpp::stop("length(yall) must equal D*T.");
  if ((int) vardirall.n_elem != Mfull) Rcpp::stop("length(vardirall) must equal D*T.");
  if ((int) proxmat.n_rows != D || (int) proxmat.n_cols != D)
    Rcpp::stop("proxmat must be a square D x D matrix.");
  if (!Xall.is_finite()) Rcpp::stop("Xall must be finite in every cell (also where y is NA).");

  // =====================================================================
  // Observed cells: y not NA. The others (intermittent gaps, unsampled
  // domains) are predicted from the model fitted to the observed cells.
  // X, y, vardir below are the compact (observed-only) data, domain-major.
  // =====================================================================
  arma::uvec observed(Mfull);
  for (int i = 0; i < Mfull; ++i) observed(i) = std::isfinite(yall(i)) ? 1u : 0u;
  const StIndex ix = st_make_index(D, Tt, observed);
  const int M = ix.n;
  if (M == 0) Rcpp::stop("All cells are unsampled (y is NA everywhere) -- cannot fit the model.");
  // balanced panel: use the inputs themselves (same memory, same arithmetic)
  const arma::mat Xc = ix.full ? arma::mat() : arma::mat(Xall.rows(ix.obs));
  const arma::vec yc = ix.full ? arma::vec() : arma::vec(yall.elem(ix.obs));
  const arma::vec vc = ix.full ? arma::vec() : arma::vec(vardirall.elem(ix.obs));
  const arma::mat& X = ix.full ? Xall : Xc;
  const arma::vec& y = ix.full ? yall : yc;
  const arma::vec& vardir = ix.full ? vardirall : vc;
  if (!vardir.is_finite() || arma::any(vardir <= 0.0))
    Rcpp::stop("vardir must be finite and positive in every sampled cell.");
  if (M <= p + (model == "ST" ? 4 : 3))
    Rcpp::stop("Too few sampled cells (%d) for %d fixed effects and the variance components.", M, p);

  const double med_vardir = arma::median(vardir);
  if (sigma21_start < 0) sigma21_start = 0.5 * med_vardir;
  if (sigma22_start < 0) sigma22_start = 0.5 * med_vardir;
  if (rho1_start <= -1 || rho1_start >= 1) Rcpp::stop("rho1_start must be in (-1,1).");
  if (rho2_start <= -1 || rho2_start >= 1) Rcpp::stop("rho2_start must be in (-1,1).");

  const bool isST = (model == "ST");
  const int nparam = isST ? 4 : 3;

  const arma::mat Id = arma::eye<arma::mat>(D, D);
  const arma::mat& W = proxmat;
  const arma::mat Wt = W.t();
  const arma::mat WtW = Wt * W;
  const arma::mat WpWt = W + Wt;
  const arma::mat tX = X.t();

  CharacterVector thetanames = isST
  ? CharacterVector::create("sigma21", "rho1", "sigma22", "rho2")
    : CharacterVector::create("sigma21", "rho1", "sigma22");

  arma::vec thetak(nparam), thetak1(nparam);
  thetak1(0) = sigma21_start;
  thetak1(1) = rho1_start;
  thetak1(2) = sigma22_start;
  if (isST) thetak1(3) = rho2_start;

  arma::vec S(nparam, fill::zeros);
  arma::mat F(nparam, nparam, fill::zeros);
  arma::mat Finv(nparam, nparam, fill::zeros);

  int k = 0;
  double diff = precision + 1.0;

  auto make_fail_result = [&](bool conv) -> List {
    return List::create(
      _["estcoef"] = R_NilValue,
      _["estvarcomp"] = R_NilValue,
      _["goodness"] = R_NilValue,
      _["df_eblup"] = R_NilValue,
      _["model"] = model,
      _["convergence"] = conv,
      _["n_iter"] = k
    );
  };

  bool used_fallback = false;
  bool primary_failed = false;
  while (diff > precision && k < maxiter) {
    ++k;
    thetak = thetak1;

    // Score and Fisher information, structured (no M x M matrices; see
    // st_score_info). Algebraically identical to the dense computation.
    StPieces sp;
    if (!st_build(thetak, isST, ix, vardir, W, sp)) { primary_failed = true; break; }
    arma::mat ViX_k, Q_k;
    arma::vec Py_k;
    if (!st_score_info(sp, thetak, X, y, W, WtW, S, F, ViX_k, Q_k, Py_k)) { primary_failed = true; break; }

    // Update theta
    arma::mat symF = arma::symmatu(F);
    bool okF = arma::inv_sympd(Finv, symF);
    if (!okF) okF = arma::inv(Finv, symF);
    if (!okF) { primary_failed = true; break; }

    thetak1 = thetak + Finv * S;

    // Clamp parameters
    if (thetak1(1) <= -1) thetak1(1) = -0.999;
    else if (thetak1(1) >= 1) thetak1(1) = 0.999;
    if (isST) {
      if (thetak1(3) <= -1) thetak1(3) = -0.999;
      else if (thetak1(3) >= 1) thetak1(3) = 0.999;
    }

    // Convergence check
    arma::vec thetak_safe = thetak;
    for (int i = 0; i < nparam; ++i) if (thetak_safe(i) == 0.0) thetak_safe(i) = 1e-4;
    diff = arma::max(arma::abs((thetak_safe - thetak1) / thetak_safe));
  }

  if (primary_failed || (k >= maxiter && diff >= precision)) {
    // The classical (sae-compatible) Fisher scoring failed -- typically a
    // variance iterate went negative and V became indefinite. Fall back to
    // constrained REML (projected scoring + step halving) from the starting
    // values. Fits on which the classical iteration converges are unchanged.
    arma::vec th(nparam);
    th(0) = sigma21_start; th(1) = rho1_start; th(2) = sigma22_start;
    if (isST) th(3) = rho2_start;
    int it2 = 0;
    const bool ok2 = stfh_fit_safeguarded(X, y, vardir, W, ix, isST,
                                          maxiter, precision, th, it2);
    k += it2;
    if (!ok2) return make_fail_result(false);
    thetak1 = th;
    used_fallback = true;
  } else if (thetak1(0) < 0.0 || thetak1(2) < 0.0) {
    // The classical iteration converged to a point with a NEGATIVE variance;
    // truncating it to 0 (as sae::eblupSTFH does) leaves the other parameters
    // at values that are not optimal given that variance = 0. Re-maximise the
    // REML likelihood under the constraints, starting from the truncated point.
    arma::vec th = thetak1;
    int it2 = 0;
    const bool ok2 = stfh_fit_safeguarded(X, y, vardir, W, ix, isST,
                                          maxiter, precision, th, it2);
    k += it2;
    if (!ok2) return make_fail_result(false);
    thetak1 = th;
    used_fallback = true;
  }

  // ==========================================================================
  // Finalization
  // ==========================================================================
  thetak1(0) = std::max(thetak1(0), 0.0);
  thetak1(2) = std::max(thetak1(2), 0.0);

  const double sigma21_f = thetak1(0);
  const double rho1_f = thetak1(1);
  const double sigma22_f = thetak1(2);
  const double rho2_f = isST ? thetak1(3) : 0.0;

  const bool param_invalid = (sigma21_f < 0) || (rho1_f < -1) || (rho1_f > 1) ||
    (sigma22_f < 0) || (isST && (rho2_f < -1 || rho2_f > 1));

  NumericVector est_vec(thetak1.begin(), thetak1.end());
  est_vec.names() = thetanames;

  if (param_invalid) {
    DataFrame estvarcomp = DataFrame::create(
      _["estimate"] = est_vec,
      _["std.error"] = NumericVector(nparam, 0.0)
    );
    return List::create(
      _["estcoef"] = R_NilValue,
      _["estvarcomp"] = estvarcomp,
      _["goodness"] = R_NilValue,
      _["df_eblup"] = R_NilValue,
      _["model"] = model,
      _["convergence"] = false,
      _["n_iter"] = k
    );
  }

  auto fail_final = [&]() -> List {
    DataFrame estvarcomp = DataFrame::create(
      _["estimate"] = est_vec,
      _["std.error"] = NumericVector(nparam, 0.0)
    );
    return List::create(
      _["estcoef"] = R_NilValue,
      _["estvarcomp"] = estvarcomp,
      _["goodness"] = R_NilValue,
      _["df_eblup"] = R_NilValue,
      _["model"] = model,
      _["convergence"] = false,
      _["n_iter"] = k
    );
  };

  // ==========================================================================
  // Everything below uses the structured representation of V (StPieces):
  //   V^-1 = A^-1 - A^-1 Z C^-1 Z' A^-1 and P = A^-1 - U Wm U', with A^-1
  //   block diagonal (T x T blocks), U = [A^-1 Z, V^-1 X], Wm = diag(C^-1, Q).
  // No M x M matrix is formed (M = number of observed cells).
  // ==========================================================================
  StPieces spf;
  if (!st_build(thetak1, isST, ix, vardir, W, spf)) return fail_final();
  const bool haveVu1 = spf.hasVu1;

  // REML score/information at the final estimates, V^-1 X and Q
  arma::vec Sf, Pyf;
  arma::mat Ff, ViXf, Q;
  if (!st_score_info(spf, thetak1, X, y, W, WtW, Sf, Ff, ViXf, Q, Pyf)) return fail_final();
  Ff = arma::symmatu(Ff);

  // Beta and residuals
  arma::vec beta = Q * (ViXf.t() * y);
  arma::vec resid = y - X * beta;
  arma::vec invVresid = st_Vinv(spf, arma::mat(resid));

  // Random effects
  const arma::vec tZ1invVresid = st_Zt(spf, arma::mat(invVresid));
  arma::vec u1est = (sigma21_f * spf.Omega1) * tZ1invVresid;
  arma::vec u1dt(M);
  {
    int first, last;
    for (int d = 0; d < D; ++d)
      if (st_rng(spf, d, first, last)) u1dt.subvec(first, last).fill(u1est(d));
  }

  arma::vec u2dt(M);
  if (!isST) {
    u2dt = sigma22_f * invVresid;
  } else {
    arma::mat sigma22Omega2 = sigma22_f * spf.Omega2;
    int first, last;
    for (int d = 0; d < D; ++d) {
      if (!st_rng(spf, d, first, last)) continue;
      if (spf.full) u2dt.subvec(first, last) = sigma22Omega2 * invVresid.subvec(first, last);
      else u2dt.subvec(first, last) = st_sub(spf, sigma22Omega2, d) * invVresid.subvec(first, last);
    }
  }

  arma::vec eblup = X * beta + u1dt + u2dt;

  // Log-likelihood: log|V| = log|A| + log|V_u1| + log|C|
  double logdetV;
  const double quadform = arma::dot(resid, invVresid);
  if (haveVu1) {
    arma::mat ImrWf = Id - rho1_f * W;
    arma::mat A1f = ImrWf.t() * ImrWf;
    double logdetA1f, sign_;
    arma::log_det(logdetA1f, sign_, A1f);
    const double logdetVu1 = D * std::log(sigma21_f) - logdetA1f;

    double logdetC, signC;
    arma::log_det(logdetC, signC, spf.Cmat);

    logdetV = spf.logdetA + logdetVu1 + logdetC;
  } else {
    logdetV = spf.logdetA;
  }

  // Full (ML-form) log-likelihood at the REML estimates, as in sae::eblupSTFH,
  // so that AIC/BIC compare models with different covariates; the restricted
  // log-likelihood (the REML objective) is reported next to it.
  const double loglike = -0.5 * (M * std::log(2.0 * M_PI) + logdetV + quadform);
  const double AIC = -2.0 * loglike + 2.0 * (p + nparam);
  const double BIC = -2.0 * loglike + std::log((double) M) * (p + nparam);
  double ldQ = 0.0, sgQ = 0.0, ldXX = 0.0, sgXX = 0.0;
  arma::log_det(ldQ, sgQ, Q);
  arma::log_det(ldXX, sgXX, arma::mat(X.t() * X));
  const double loglike_reml = loglike + 0.5 * p * std::log(2.0 * M_PI) + 0.5 * ldQ + 0.5 * ldXX;

  NumericVector goodness = NumericVector::create(
    _["loglikelihood"] = loglike, _["AIC"] = AIC, _["BIC"] = BIC,
    _["loglikelihood_reml"] = loglike_reml
  );

  // SE for beta
  arma::vec stderr_beta = arma::sqrt(Q.diag());
  arma::vec tvalue = beta / stderr_beta;
  arma::vec pvalue(p);
  for (int i = 0; i < p; ++i)
    pvalue(i) = 2.0 * R::pnorm(std::fabs(tvalue(i)), 0.0, 1.0, 0, 0);

  DataFrame estcoef = DataFrame::create(
    _["beta"] = NumericVector(beta.begin(), beta.end()),
    _["std.error"] = NumericVector(stderr_beta.begin(), stderr_beta.end()),
    _["tvalue"] = NumericVector(tvalue.begin(), tvalue.end()),
    _["pvalue"] = NumericVector(pvalue.begin(), pvalue.end())
  );

  arma::uvec unidentified;
  // Inverse information of the last classical iteration (what
  // sae::eblupSTFH reports as the SEs of theta); kept for the SEs.
  const arma::mat Finv_iter = Finv;
  {
    // At a boundary a parameter may be unidentified: with sigma21 = 0,
    // dV/drho1 = 0, so its row/column of F is zero. Invert F on the
    // identified parameters only; an unidentified one gets zero variance
    // (which is also its zero contribution to g3) and an NA standard error.
    // Correlations estimated on the boundary (|rho| >= 0.999) are treated as
    // fixed as well: the Prasad-Rao expansion in rho is not valid there.
    const double tol_id = 1e-10 * std::max(1.0, Ff.diag().max());
    arma::uvec keep(nparam, fill::ones);
    for (int a = 0; a < nparam; ++a) if (Ff(a, a) <= tol_id) keep(a) = 0;
    if (std::fabs(rho1_f) >= 0.999) keep(1) = 0;
    if (isST && std::fabs(rho2_f) >= 0.999) keep(3) = 0;
    arma::uvec id_par = arma::find(keep == 1);
    Finv.zeros(nparam, nparam);
    if (id_par.n_elem > 0) {
      arma::mat Fsub = Ff.submat(id_par, id_par), Fsub_inv;
      if (!arma::inv_sympd(Fsub_inv, Fsub)) Fsub_inv = arma::pinv(Fsub);
      Finv.submat(id_par, id_par) = Fsub_inv;
    }
    unidentified.zeros(nparam);
    for (int a = 0; a < nparam; ++a) if (keep(a) == 0) unidentified(a) = 1;
  }

  // SE for theta: as sae::eblupSTFH (last classical iteration) when that
  // iteration converged; from the information at the final estimates when
  // the constrained fallback was needed (unidentified parameters: NA).
  arma::vec diagFinv = used_fallback ? arma::vec(Finv.diag()) : arma::vec(Finv_iter.diag());
  if (arma::any(diagFinv < 0) || diagFinv.has_nan()) return fail_final();
  arma::vec stderr_theta = arma::sqrt(diagFinv);
  if (used_fallback) stderr_theta.elem(arma::find(unidentified == 1)).fill(NA_REAL);

  DataFrame estvarcomp = DataFrame::create(
    _["estimate"] = est_vec,
    _["std.error"] = NumericVector(stderr_theta.begin(), stderr_theta.end())
  );

  // ================================================================
  // Analytical MSE (Prasad-Rao / Marhuenda et al. 2013) in structured form
  // ================================================================
  // Helpers for U = [A^-1 Z, V^-1 X]: the A^-1 Z part has one non-zero per
  // row (a_i in column d(i)), so U'Y = [Z'A^-1 Y; X'V^-1 Y] and U S cost O(M).
  const StLR LR(spf, ViXf);
  const BSF Uf = bsf_of(LR);
  const arma::mat Wm = LR.Wm(Q);
  auto Ut = [&](const arma::mat& Y) { return LR.Ut(Y); };          // U' Y
  const arma::vec& avec = LR.avec;
  const arma::mat Udense = LR.dense();
  auto blk_apply = [&](const std::vector<arma::mat>& Xb, const arma::mat& Y) -> arma::mat {
    arma::mat R(Y.n_rows, Y.n_cols);
    int f, l;
    for (int d = 0; d < D; ++d) {
      if (!st_rng(spf, d, f, l)) continue;
      R.rows(f, l) = Xb[d] * Y.rows(f, l);
    }
    return R;
  };
  // per-domain blocks E[tix_d, tix_d] of a T x T matrix E
  auto blocks_of = [&](const arma::mat& E) -> std::vector<arma::mat> {
    std::vector<arma::mat> B(D);
    int f, l;
    for (int d = 0; d < D; ++d) if (st_rng(spf, d, f, l)) B[d] = st_sub(spf, E, d);
    return B;
  };
  auto Papply = [&](const arma::mat& Y) -> arma::mat {    // P Y
    return st_Ainv(spf, Y) - Uf.mul(Wm * Ut(Y));
  };
  // diag(P X P) for a block-diagonal X (blocks Xb), with P = A^-1 - U Wm U':
  //   diag(A^-1 X A^-1) - diag(A^-1 X U Wm U') - diag(U Wm U' X A^-1)
  //   + diag(U Wm U'X U Wm U')
  auto diag_PXP = [&](const std::vector<arma::mat>& Xb) -> arma::vec {
    arma::vec out(M);
    std::vector<arma::mat> AX(D), AXt(D);
    int f, l;
    for (int d = 0; d < D; ++d) {
      if (!st_rng(spf, d, f, l)) continue;
      AX[d] = spf.Ainv[d] * Xb[d];
      AXt[d] = spf.Ainv[d] * Xb[d].t();
      out.subvec(f, l) = arma::diagvec(AX[d] * spf.Ainv[d]);
    }
    const BSF AXU = bsf_transform(Uf, AX), AXtU = bsf_transform(Uf, AXt);
    out -= Uf.row(AXU.mul(Wm));
    out -= AXtU.row(Uf.mul(Wm));
    out += Uf.row(Uf.mul(Wm * Ut(blk_apply(Xb, Udense)) * Wm));
    return out;
  };

  // g1 = psi - psi^2 diag(V^-1), diag(V^-1) = diag(A^-1) - a_i^2 C^-1_dd
  arma::vec invV_diag(M);
  for (int d = 0; d < D; ++d) {
    int f, l;
    if (!st_rng(spf, d, f, l)) continue;
    invV_diag.subvec(f, l) = arma::diagvec(spf.Ainv[d]);
    if (haveVu1) invV_diag.subvec(f, l) -= arma::square(avec.subvec(f, l)) * spf.Cinv(d, d);
  }
  arma::vec g1 = vardir - (arma::square(vardir) % invV_diag);
  for (uword i = 0; i < g1.n_elem; ++i) {
    if (g1(i) < 0.0) g1(i) = 0.0;
  }

  // g2 component: uncertainty from beta estimation
  arma::mat R_mat = arma::diagmat(vardir) * ViXf;
  arma::vec g2 = arma::sum((R_mat * Q) % R_mat, 1);
  for (uword i = 0; i < g2.n_elem; ++i) {
    if (g2(i) < 0.0) g2(i) = 0.0;
  }

  // g3 component: uncertainty from variance component estimation.
  // The EBLUP can be written as theta~(th) = y - Psi P(th) y, so its exact
  // derivative is d theta~/d th_a = Psi P V_a P y, with
  //   Cov(Psi P V_a P y, Psi P V_b P y)_ii = psi_i^2 (P V_a P V_b P)_ii
  // (using P V P = P). Hence
  //   g3_i = psi_i^2 sum_ab Finv_ab (P V_a P V_b P)_ii.
  // This agrees with the textbook Prasad-Rao form
  //   psi_i^2 sum_ab Finv_ab (V^-1 V_a V^-1 V V^-1 V_b V^-1)_ii
  // up to o(1/m), but unlike it, it is unaffected by components of V_a
  // lying in col(X), which are absorbed by beta~ (e.g. at rho1 -> 1 with a
  // row-standardised W, Omega1 is dominated by the constant vector and the
  // textbook form explodes).
  // V_a is either Z B_a Z' (sigma21, rho1) or block diagonal I (x) E_j
  // (sigma22, rho2), so with PZ = P Z:
  //   Z-Z: (P V_a P V_b P)_ii = (PZ B_a Z'PZ B_b PZ')_ii
  //   Z-E: (PZ B_a PZ' E_j P)_ii = rowsum(PZ B_a % P E_j PZ)
  //   E-E: (P E_j P E_k P)_ii = (P E_j A^-1 E_k P)_ii - ((P E_j U) Wm (P E_k U)')_ii
  arma::mat derivrho1 = -WpWt + 2.0 * rho1_f * WtW;
  const arma::mat OCO = spf.Omega1 * derivrho1 * spf.Omega1;
  const arma::mat Bz[2] = { spf.Omega1, -sigma21_f * OCO };
  std::vector<arma::mat> Eb;
  if (isST) { Eb.push_back(spf.Omega2); Eb.push_back(sigma22_f * spf.dOmega2); }
  else Eb.push_back(arma::eye<arma::mat>(Tt, Tt));
  const int nE = (int) Eb.size();
  auto wgt = [&](int a, int b) { return (a == b) ? Finv(a, b) : 2.0 * Finv(a, b); };

  // P Z = A^-1 Z - U G, G = Wm U'Z (U'Z = [diag(c); X'V^-1 Z]). With
  // sigma21 > 0, A^-1 Z is the a-part of U, so P Z = U K with K = [I; 0] - G
  // and diag(PZ S PZ') = rowsum(U (K S K') % U): O(D^3 + M r), no M x D^2.
  const arma::mat ZtViX = st_Zt(spf, ViXf);
  arma::mat UtZ(LR.r, D);
  if (haveVu1) UtZ.rows(0, D - 1) = arma::diagmat(spf.c);
  UtZ.rows(LR.off, LR.off + p - 1) = ZtViX.t();
  const arma::mat G = Wm * UtZ;
  arma::mat Kz;
  arma::mat PZd;                                   // dense P Z (sigma21 = 0 only)
  if (haveVu1) {
    Kz = -G;
    Kz.submat(0, 0, D - 1, D - 1).diag() += 1.0;
  } else {
    PZd = spf.AZ - Uf.mul(G);
  }
  auto PZmul = [&](const arma::mat& Sm) -> arma::mat {        // (P Z) S
    return haveVu1 ? Uf.mul(Kz * Sm) : arma::mat(PZd * Sm);
  };
  auto diagPZ = [&](const arma::mat& Sm) -> arma::vec {       // diag(PZ S PZ')
    return haveVu1 ? Uf.row(Uf.mul(Kz * Sm * Kz.t())) : arma::vec(arma::sum((PZd * Sm) % PZd, 1));
  };
  arma::mat ZPZ = arma::diagmat(spf.c) - UtZ.t() * G;
  ZPZ = 0.5 * (ZPZ + ZPZ.t());

  arma::vec g3s(M, fill::zeros);
  {
    arma::mat Szz(D, D, fill::zeros);
    for (int a = 0; a < 2; ++a)
      for (int b = a; b < 2; ++b)
        if (Finv(a, b) != 0.0) Szz += wgt(a, b) * (Bz[a] * ZPZ * Bz[b]);
    if (arma::any(arma::vectorise(Szz) != 0.0)) g3s += diagPZ(Szz);

    const arma::mat PZ = PZmul(arma::eye<arma::mat>(D, D));
    std::vector<arma::mat> UtEU(nE);
    for (int j = 0; j < nE; ++j) {
      const int aj = 2 + j;
      bool needZE = false;
      for (int a = 0; a < 2; ++a) needZE = needZE || (Finv(a, aj) != 0.0);
      if (needZE) {
        const arma::mat Yj = Papply(st_blk(spf, Eb[j], PZ));        // P E_j P Z
        arma::mat Bsum(D, D, fill::zeros);
        for (int a = 0; a < 2; ++a) if (Finv(a, aj) != 0.0) Bsum += wgt(a, aj) * Bz[a];
        g3s += arma::sum(PZmul(Bsum) % Yj, 1);
      }
      UtEU[j] = Ut(st_blk(spf, Eb[j], Udense));                     // U' E_j U
    }
    // E-E: (P E_j P E_k P)_ii = diag(P E_j A^-1 E_k P) - diag((P E_j U) Wm (P E_k U)'),
    // P E_j U = A^-1 E_j U - U N_j with N_j = Wm U'E_j U
    for (int j = 0; j < nE; ++j) {
      for (int kk = j; kk < nE; ++kk) {
        const double w = wgt(2 + j, 2 + kk);
        if (w == 0.0) continue;
        std::vector<arma::mat> Xblk(D), AEj(D), AEk(D);
        int f, l;
        for (int d = 0; d < D; ++d) {
          if (!st_rng(spf, d, f, l)) continue;
          if (spf.full) {
            Xblk[d] = Eb[j] * spf.Ainv[d] * Eb[kk];
            AEj[d] = spf.Ainv[d] * Eb[j];
            AEk[d] = spf.Ainv[d] * Eb[kk];
          } else {
            const arma::mat Ejd = st_sub(spf, Eb[j], d), Ekd = st_sub(spf, Eb[kk], d);
            Xblk[d] = Ejd * spf.Ainv[d] * Ekd;
            AEj[d] = spf.Ainv[d] * Ejd;
            AEk[d] = spf.Ainv[d] * Ekd;
          }
        }
        const BSF Bj = bsf_transform(Uf, AEj), Bk = bsf_transform(Uf, AEk);
        const arma::mat Nj = Wm * UtEU[j], Nk = Wm * UtEU[kk];
        const arma::vec cross = Bk.row(Bj.mul(Wm)) - Uf.row(Bj.mul(Wm * Nk.t())) -
          Bk.row(Uf.mul(Nj * Wm)) + Uf.row(Uf.mul(Nj * Wm * Nk.t()));
        g3s += w * (diag_PXP(Xblk) - cross);
      }
    }
  }
  const arma::vec psi2 = arma::square(vardir);
  arma::vec g3 = psi2 % g3s;
  for (uword i = 0; i < g3.n_elem; ++i) {
    if (g3(i) < 0.0) g3(i) = 0.0;
  }

  // g4 component: V is nonlinear in rho1 and rho2, so the plug-in g1(th^)
  // has an extra bias 1/2 sum_jk Finv_jk psi_i^2 (V^-1 V_jk V^-1)_ii on top of
  // the -g3 part (as the g4 term of the SFH model, Singh et al. 2005). In the
  // P form used for g3 (equal up to o(1/m), robust to col(X) components):
  //   g4_i = 1/2 sum_jk Finv_jk psi_i^2 (P V_jk P)_ii,
  // with the non-zero second derivatives V_jk = d^2 V / d th_j d th_k:
  //   (sigma21, rho1): (d Omega1 / d rho1) (x) J = -(Omega1 C Omega1) (x) J
  //   (rho1, rho1)   : sigma21 (2 Omega1 C Omega1 C Omega1 - 2 Omega1 W'W Omega1) (x) J
  //   (sigma22, rho2): I (x) dOmega2               [ST only]
  //   (rho2, rho2)   : sigma22 I (x) d2Omega2      [ST only]
  arma::vec g4(M, fill::zeros);
  {
    auto add_z = [&](int a, int b, const arma::mat& K) {      // V_ab = Z K Z'
      const double w = wgt(a, b);
      if (w == 0.0) return;
      g4 += 0.5 * w * (psi2 % diagPZ(K));
    };
    auto add_e = [&](int a, int b, const arma::mat& E) {      // V_ab = I (x) E
      const double w = wgt(a, b);
      if (w == 0.0) return;
      g4 += 0.5 * w * (psi2 % diag_PXP(blocks_of(E)));
    };
    add_z(0, 1, -OCO);
    const arma::mat d2Om1 = 2.0 * OCO * derivrho1 * spf.Omega1 -
      2.0 * spf.Omega1 * WtW * spf.Omega1;
    add_z(1, 1, sigma21_f * d2Om1);
    if (isST) {
      add_e(2, 3, spf.dOmega2);
      add_e(3, 3, sigma22_f * build_omega2_d2(Tt, rho2_f));
    }
  }

  arma::vec mse_anal = g1 + g2 + 2.0 * g3 - g4;

  // ================================================================
  // All D x T cells. Observed cells keep the EBLUP and MSE above. An
  // unobserved cell i = (d, t) (gap, or a domain without any sample) gets the
  // BLUP of theta_i = x_i' beta + u1_d + u2_dt from the observed data
  // (Saei & Chambers 2005; Benedetti et al. 2024 for the SAR effect),
  //   theta~_i = x_i' beta~ + a_i' V^-1 (y - X beta~),
  //   a_i = Cov(theta_i, y) : sigma21 Omega1[d, d(k)] + 1{d(k) = d} sigma22 Omega2[t, t(k)]
  // (model "S": no second term, the u2 of an unobserved cell is independent
  // of the data). Its MSE is
  //   g1 = G_ii - a_i' V^-1 a_i,  g2 = R_i Q R_i',  R_i = x_i' - a_i' V^-1 X,
  //   g3 = sum_ab Finv_ab l_a' P l_b,  l_a' = A_a - (a'V^-1 + R_i Q X'V^-1) V_a
  //   (exact derivative of the predictor in the variance parameters),
  //   MSE = g1 + g2 + g3.
  // The second-order correction of the plug-in g1, 1/2 sum_jk Finv_jk d2 g1 /
  // d th_j d th_k, is not applied: it is numerically unstable when rho1 is
  // close to +-1 (nearly singular Omega1) or a variance component is weakly
  // identified, and in Monte Carlo (benchmarks/sim_stfh_unsampled.R) g1 + g2 +
  // g3 has a relative bias within a few percent and ~95% coverage while the
  // corrected form is dominated by a few replicates with an exploding term.
  // ================================================================
  const arma::uvec unobs = arma::find(observed == 0);
  const arma::uword r = unobs.n_elem;
  arma::vec u1_full, u2_full;
  st_full_u(spf, sigma21_f, sigma22_f, invVresid, u1_full, u2_full);
  for (int d = 0; d < D; ++d) u1_full.subvec(d * Tt, d * Tt + Tt - 1).fill(u1est(d));
  arma::vec eblup_full = Xall * beta + (u1_full + u2_full);
  eblup_full.elem(ix.obs) = eblup;
  u2_full.elem(ix.obs) = u2dt;
  arma::vec mse_full(Mfull), g1_full(Mfull), g2_full(Mfull), g3_full(Mfull), g4_full(Mfull);
  mse_full.elem(ix.obs) = mse_anal;
  g1_full.elem(ix.obs) = g1;
  g2_full.elem(ix.obs) = g2;
  g3_full.elem(ix.obs) = g3;
  g4_full.elem(ix.obs) = g4;

  if (r > 0) {
    const arma::uvec dr = unobs / Tt, tr = unobs - dr * Tt;      // domain, period
    const arma::uvec dk = ix.obs / Tt, tk = ix.obs - dk * Tt;     // of the observed rows
    // r x M blocks of a and its derivatives:
    //   ZZ(K)_ik = K[d_i, d(k)],   EE(E)_ik = 1{d(k) = d_i} E[t_i, t(k)]
    auto ZZ = [&](const arma::mat& K) -> arma::mat { return K.submat(dr, dk); };
    auto EE = [&](const arma::mat& E) -> arma::mat {
      arma::mat R_(r, M, fill::zeros);
      int f, l;
      for (arma::uword i = 0; i < r; ++i) {
        if (!st_rng(spf, (int) dr(i), f, l)) continue;
        R_.submat(i, f, i, l) = E.submat(arma::uvec{tr(i)}, spf.tix[dr(i)]);
      }
      return R_;
    };
    auto dZZ = [&](const arma::mat& K) -> arma::vec {             // K[d_i, d_i]
      arma::vec v(r);
      for (arma::uword i = 0; i < r; ++i) v(i) = K(dr(i), dr(i));
      return v;
    };
    auto dEE = [&](const arma::mat& E) -> arma::vec {             // E[t_i, t_i]
      arma::vec v(r);
      for (arma::uword i = 0; i < r; ++i) v(i) = E(tr(i), tr(i));
      return v;
    };

    const arma::mat d1Om1 = -OCO;                                  // d Omega1 / d rho1
    const arma::mat Ir = arma::eye<arma::mat>(Tt, Tt);
    const arma::mat zero_rM(r, M, fill::zeros);

    // a (r x M) and Var(u_i)
    arma::mat A = sigma21_f * ZZ(spf.Omega1);
    if (isST && sigma22_f != 0.0) A += sigma22_f * EE(spf.Omega2);
    const arma::vec Grr = sigma21_f * dZZ(spf.Omega1) +
      (isST ? arma::vec(sigma22_f * dEE(spf.Omega2)) : arma::vec(r, fill::value(sigma22_f)));

    // first derivatives of a (parameter order sigma21, rho1, sigma22, rho2)
    std::vector<arma::mat> Aj = { ZZ(spf.Omega1), sigma21_f * ZZ(d1Om1),
                                  isST ? EE(spf.Omega2) : zero_rM };
    if (isST) Aj.push_back(sigma22_f * EE(spf.dOmega2));

    // V_j Y for the compact data (Y: M x k)
    auto Vz = [&](const arma::mat& K, const arma::mat& Y) -> arma::mat { return st_Z(spf, K * st_Zt(spf, Y)); };
    auto Vj = [&](int j, const arma::mat& Y) -> arma::mat {
      switch (j) {
      case 0: return Vz(spf.Omega1, Y);
      case 1: return Vz(sigma21_f * d1Om1, Y);
      case 2: return isST ? st_blk(spf, spf.Omega2, Y) : st_blk(spf, Ir, Y);
      default: return st_blk(spf, sigma22_f * spf.dOmega2, Y);
      }
    };
    auto Kmul = [&](const arma::mat& Y) -> arma::mat { return st_Vinv(spf, Y); };            // V^-1 Y
    auto Pmul = [&](const arma::mat& Y) -> arma::mat {                                        // P Y
      return st_Vinv(spf, Y) - ViXf * (Q * (ViXf.t() * Y));
    };
    auto dg = [](const arma::mat& M1, const arma::mat& M2) -> arma::vec { return arma::sum(M1 % M2, 1); };

    const arma::mat KAt = Kmul(A.t());                    // V^-1 a      (M x r)
    const arma::mat AK = KAt.t();                         // a' V^-1     (r x M)
    const arma::mat Xr = Xall.rows(unobs);
    const arma::mat Rr = Xr - AK * X;                     // x_i' - a'V^-1 X
    arma::vec g1o = Grr - dg(A, AK);
    arma::vec g2o = dg(Rr * Q, Rr);

    // g3: exact derivative of the BLUP, d theta~_i / d th_a = l_a' P y
    const arma::mat Bmat = AK + Rr * Q * ViXf.t();        // r x M
    std::vector<arma::mat> L(nparam), PLt(nparam);
    for (int a = 0; a < nparam; ++a) {
      L[a] = Aj[a] - Vj(a, Bmat.t()).t();
      PLt[a] = Pmul(L[a].t());
    }
    arma::vec g3o(r, fill::zeros);
    for (int a = 0; a < nparam; ++a)
      for (int b = 0; b < nparam; ++b)
        if (Finv(a, b) != 0.0) g3o += Finv(a, b) * dg(L[a], PLt[b].t());

    g1o.transform([](double v) { return v < 0.0 ? 0.0 : v; });
    g2o.transform([](double v) { return v < 0.0 ? 0.0 : v; });
    g3o.transform([](double v) { return v < 0.0 ? 0.0 : v; });

    mse_full.elem(unobs) = g1o + g2o + g3o;
    g1_full.elem(unobs) = g1o;
    g2_full.elem(unobs) = g2o;
    g3_full.elem(unobs) = g3o;
    g4_full.elem(unobs).zeros();
  }

  arma::vec rse_full(Mfull);
  for (int i = 0; i < Mfull; ++i) {
    rse_full(i) = (std::abs(eblup_full(i)) < 1e-12 || mse_full(i) < 0.0) ? NA_REAL
      : (std::sqrt(mse_full(i)) / std::abs(eblup_full(i)) * 100.0);
  }
  LogicalVector sampled(Mfull);
  for (int i = 0; i < Mfull; ++i) sampled[i] = observed(i) == 1u;

  // df_eblup with random effects, for all D x T cells (domain-major)
  DataFrame df_eblup = DataFrame::create(
    _["eblup"] = NumericVector(eblup_full.begin(), eblup_full.end()),
    _["mse"] = NumericVector(mse_full.begin(), mse_full.end()),
    _["rse"] = NumericVector(rse_full.begin(), rse_full.end()),
    _["random_effect_u1"] = NumericVector(u1_full.begin(), u1_full.end()),
    _["random_effect_u2"] = NumericVector(u2_full.begin(), u2_full.end()),
    _["sampled"] = sampled
  );

  // Return structure matching seblup_area + estvarcomp + analytical MSE
  return List::create(
    _["estcoef"] = estcoef,
    _["estvarcomp"] = estvarcomp,
    _["goodness"] = goodness,
    _["df_eblup"] = df_eblup,
    _["mse_analytical"] = NumericVector(mse_full.begin(), mse_full.end()),
    _["rse_analytical"] = NumericVector(rse_full.begin(), rse_full.end()),
    _["g1"] = NumericVector(g1_full.begin(), g1_full.end()),
    _["g2"] = NumericVector(g2_full.begin(), g2_full.end()),
    _["g3"] = NumericVector(g3_full.begin(), g3_full.end()),
    _["g4"] = NumericVector(g4_full.begin(), g4_full.end()),
    _["n_sampled"] = M,
    _["model"] = model,
    _["convergence"] = true,
    _["n_iter"] = k
  );
}
