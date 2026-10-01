// src/eblup_tfh.cpp
// Two-fold Fay-Herriot model (area-level with subareas), Torabi & Rao (2014).
//
//   Level 1 (sampling): y_ij = theta_ij + e_ij,  e_ij ~ N(0, psi_ij)
//   Level 2 (linking):  theta_ij = x_ij'beta + v_i + u_ij,
//                       v_i ~ N(0, s2v), u_ij ~ N(0, s2u)
//
// Paper notation: s2v = AREA variance, s2u = SUBAREA variance.
//
// Variance components (s2v, s2u) by ML/REML Fisher-scoring exploiting the
// block-diagonal covariance (Woodbury per area block). EBLUP from the BLUP
// equations. MSE = g1 + g2 + g3 (Prasad-Rao type), with g3 from an exact
// matrix derivation:
//
//   d mu~/d s2u = (V^{-1} - S V^{-2}) r  =: Mu r
//   d mu~/d s2v = (11'V^{-1} - S V^{-1} 11' V^{-1}) r =: Mv r
//   E[(d mu~/d s2u)^2] = diag(Mu V Mu'), etc.
//
// NOTE (2026-10-01): Torabi & Rao (2014) eq. (3.4) was investigated and found
// to contain a spurious [A][B] term and a structurally incorrect var_u
// coefficient K (wrong shape, not just scale; paper s4 also has wrong sign).
// The matrix-based g3 above was validated against finite-difference
// derivatives (cor 0.97/0.96/0.74). A proposed sigma_u^-2 "correction" was
// rejected as dimensionally wrong.
//
// REVISION NOTES (this file):
//  [bug]  beta/Q returned by tfh_fit now correspond to the FINAL (s2v, s2u)
//         (before: they came from the last iterate's variances).
//  [bug]  MSE of non-sampled subareas in sampled areas: g1*/g2* were wrong
//         (missing -s2v^2 * 1'V^{-1}1 and +s2v^2 (X'w)'Q(X'w)); fixed and
//         checked numerically against the exact variance.
//  [bug]  mse_all was uninitialised for non-sampled subareas in areas with NO
//         sampled subarea; now x'Qx + s2v + s2u (synthetic predictor).
//  [bug]  Bootstrap predictor for non-sampled subareas in sampled areas now
//         includes the area-effect BLUP (consistent with the analytical path).
//  [bug]  log-likelihood now includes the 2*pi constant (ML) / REML terms.
//  [stat] Inverse information for g3 now uses the method-specific (ML or
//         REML) Fisher information at the final estimates (before: always
//         ML-type). Effect is second order; see comment in tfh_fit.
//  [misc] area index validation; v/u variable names aligned with the paper.
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

struct AreaScratch {
  vec dd;      // 1 / (s2u + psi_d)
  double s = 0.0, c = 0.0;
  double dd2 = 0.0, dd3 = 0.0;
  mat A;       // V_d^{-1} X_d
  vec viny;    // V_d^{-1} y_d
};

struct TfhFit {
  vec beta;
  double s2v = 0.0, s2u = 0.0;  // paper notation: v=area, u=subarea
  mat Q;      // (X'V^{-1}X)^{-1}
  mat InfoInv;  // inverse Fisher information for (s2v, s2u)
  int niter = 0;
  bool conv = false;
  double loglik = 0.0;
};

// Everything that depends on (s2v, s2u) in one pass over the areas:
// beta, Q, score vector and Fisher information (ML or REML).
struct TfhPass {
  vec beta;
  mat Q;
  vec sc;     // (score_v, score_u)
  mat Info;   // 2x2, order (v, u)
};

TfhPass tfh_pass(const mat& Xs, const vec& ys, const vec& psis,
                 const std::vector<uvec>& sidx, int m_all, bool is_ml,
                 double s2v, double s2u) {
  const int p = Xs.n_cols;
  TfhPass o;

  mat XtVX(p, p, fill::zeros);
  vec Xty(p, fill::zeros);
  std::vector<AreaScratch> scr(m_all);

  for (int d = 0; d < m_all; ++d) {
    const uvec& id = sidx[d];
    if (id.n_elem == 0) continue;
    AreaScratch& a = scr[d];
    const mat Xd = Xs.rows(id);
    const vec yd = ys.elem(id);
    const vec psid = psis.elem(id);

    a.dd = 1.0 / (s2u + psid);
    a.s = accu(a.dd);
    a.c = (s2v > 0.0) ? s2v / (1.0 + s2v * a.s) : 0.0;
    a.dd2 = accu(square(a.dd));
    a.dd3 = accu(a.dd % square(a.dd));

    const double ddy = dot(a.dd, yd);
    a.viny = a.dd % yd - a.c * a.dd * ddy;

    a.A = Xd.each_col() % a.dd;
    a.A -= a.c * a.dd * (a.dd.t() * Xd);

    XtVX += Xd.t() * a.A;
    Xty += Xd.t() * a.viny;
  }

  mat Q;
  if (!inv_sympd(Q, XtVX)) Q = solve(XtVX, eye<mat>(p, p));
  const vec tvec = Q * Xty;

  // Score and information for (s2v, s2u). Order: 0=v (area), 1=u (subarea).
  double sv = 0.0, su = 0.0, Ivv = 0.0, Iuu = 0.0, Ivu = 0.0;
  for (int d = 0; d < m_all; ++d) {
    const uvec& id = sidx[d];
    if (id.n_elem == 0) continue;
    const AreaScratch& a = scr[d];

    const vec r = a.viny - a.A * tvec;        // V^{-1}(y - X beta) = P y (REML)
    const double sumr = accu(r);
    const double trV = accu(a.dd) - a.c * a.dd2;
    const double s1 = a.s * (1.0 - a.c * a.s); // 1'V^{-1}1
    const vec z = a.dd * (1.0 - a.c * a.s);    // V^{-1}1

    if (!is_ml) {
      const vec avec = sum(a.A, 0).t();        // X'V^{-1}1
      const mat M = a.A.t() * a.A;             // X'V^{-2}X
      const double qd = s1 - as_scalar(avec.t() * Q * avec);
      const double trP = trV - trace(Q * M);
      sv += -0.5 * qd + 0.5 * sumr * sumr;
      su += -0.5 * trP + 0.5 * dot(r, r);
      Ivv += 0.5 * qd * qd;

      mat VinvA = a.A.each_col() % a.dd;
      VinvA -= a.c * a.dd * (a.dd.t() * a.A);
      const mat G = a.A.t() * VinvA;
      const mat MQ = M * Q;
      const double trV2 = a.dd2 - 2.0 * a.c * a.dd3 + a.c * a.c * a.dd2 * a.dd2;
      const double trP2 = trV2 - 2.0 * trace(Q * G) + trace(MQ * MQ);
      Iuu += 0.5 * trP2;

      const rowvec zA = z.t() * a.A;
      const double term = dot(z, z) - 2.0 * as_scalar(zA * Q * avec) +
        as_scalar(avec.t() * Q * M * Q * avec);
      Ivu += 0.5 * term;
    } else {
      sv += -0.5 * s1 + 0.5 * sumr * sumr;
      su += -0.5 * trV + 0.5 * dot(r, r);
      Ivv += 0.5 * s1 * s1;
      Iuu += 0.5 * (a.dd2 - 2.0 * a.c * a.dd3 + a.c * a.c * a.dd2 * a.dd2);
      Ivu += 0.5 * dot(z, z);
    }
  }

  o.beta = tvec;
  o.Q = Q;
  o.sc = vec(2);
  o.sc(0) = sv; o.sc(1) = su;
  o.Info = mat(2, 2);
  o.Info(0, 0) = Ivv; o.Info(0, 1) = Ivu;
  o.Info(1, 0) = Ivu; o.Info(1, 1) = Iuu;
  return o;
}

// Fisher-scoring for (s2v, s2u). sidx[d] = row positions of area d.
// want_loglik = false skips the log-likelihood (used inside the bootstrap).
TfhFit tfh_fit(const mat& Xs, const vec& ys, const vec& psis,
               const std::vector<uvec>& sidx, int m_all,
               const std::string& method, int maxiter, double precision,
               bool want_loglik = true) {
  const int p = Xs.n_cols;
  const bool is_ml = (method == "ML");
  TfhFit f;

  const double med = median(psis);
  double s2v = 0.5 * med, s2u = 0.5 * med;

  double diff = precision + 1.0;
  int k = 0;
  while ((diff > precision) && (k < maxiter)) {
    const TfhPass o = tfh_pass(Xs, ys, psis, sidx, m_all, is_ml, s2v, s2u);
    const double Ivv = o.Info(0, 0), Iuu = o.Info(1, 1), Ivu = o.Info(0, 1);
    const double sv = o.sc(0), su = o.sc(1);

    vec delta(2, fill::zeros);
    bool solved = false;
    const double det = Ivv * Iuu - Ivu * Ivu;
    if (std::fabs(det) > 1e-300 && std::isfinite(det)) {
      mat II;
      if (inv_sympd(II, o.Info)) {
        delta = II * o.sc;
        solved = true;
      }
    }
    if (!solved) {
      delta(0) = (Ivv > 1e-300) ? sv / Ivv : 0.0;
      delta(1) = (Iuu > 1e-300) ? su / Iuu : 0.0;
    }

    double s2v_new = s2v + delta(0);
    double s2u_new = s2u + delta(1);
    if (!std::isfinite(s2v_new)) s2v_new = s2v;
    if (!std::isfinite(s2u_new)) s2u_new = s2u;
    if (s2v_new < 0.0) s2v_new = 0.0;
    if (s2u_new < 0.0) s2u_new = 0.0;

    const double d1 = std::fabs((s2v_new - s2v) / std::max(s2v, 1e-12));
    const double d2 = std::fabs((s2u_new - s2u) / std::max(s2u, 1e-12));
    diff = std::max(d1, d2);
    s2v = s2v_new;
    s2u = s2u_new;
    ++k;
  }

  // Final pass at the FINAL (s2v, s2u): beta, Q and information are now
  // consistent with the reported variance components.
  const TfhPass fin = tfh_pass(Xs, ys, psis, sidx, m_all, is_ml, s2v, s2u);

  f.s2v = s2v;
  f.s2u = s2u;
  f.beta = fin.beta;
  f.Q = fin.Q;
  f.niter = k;
  f.conv = (diff <= precision);

  // Inverse information (ML or REML, matching `method`) for the g3 variances.
  {
    mat II;
    if (inv_sympd(II, fin.Info)) f.InfoInv = II;
    else f.InfoInv = zeros<mat>(2, 2);
  }

  if (want_loglik) {
    double core = 0.0;  // sum_d [ log|V_d| + r'V_d^{-1}r ]
    for (int d = 0; d < m_all; ++d) {
      const uvec& id = sidx[d];
      if (id.n_elem == 0) continue;
      const mat Xd = Xs.rows(id);
      const vec yd = ys.elem(id);
      const vec psid = psis.elem(id);
      const vec dd = 1.0 / (s2u + psid);
      const double s = accu(dd);
      const double c = (s2v > 0.0) ? s2v / (1.0 + s2v * s) : 0.0;
      const vec res = yd - Xd * f.beta;
      const vec vr = dd % res - c * dd * dot(dd, res);
      core += accu(log(s2u + psid)) + std::log(1.0 + s2v * s) + dot(res, vr);
    }
    const double n = (double)Xs.n_rows;
    const double l2pi = std::log(2.0 * datum::pi);
    if (is_ml) {
      f.loglik = -0.5 * (n * l2pi + core);
    } else {
      // Restricted log-likelihood:
      // -1/2 [ (n-p) log 2pi + log|V| + log|X'V^{-1}X| - log|X'X| + r'V^{-1}r ]
      double ldQ = 0.0, sgQ = 0.0, ldXX = 0.0, sgXX = 0.0;
      log_det(ldQ, sgQ, f.Q);                       // log|Q| = -log|X'V^{-1}X|
      const mat XtX = Xs.t() * Xs;
      log_det(ldXX, sgXX, XtX);
      f.loglik = -0.5 * ((n - p) * l2pi + core - ldQ - ldXX);
    }
  }
  return f;
}

}  // namespace

// [[Rcpp::export(.eblup_tfh_core)]]
List eblup_tfh_core(const arma::mat& Xall, const arma::vec& yall,
                    const arma::vec& vardirall, const arma::ivec& area,
                    std::string method = "REML", std::string mse_type = "analytical",
                    int B = 200, int maxiter = 100, double precision = 1e-4) {
  const int N = Xall.n_rows;
  const int p = Xall.n_cols;

  if (!(method == "ML" || method == "REML")) Rcpp::stop("method must be 'ML' or 'REML'.");
  if (!(mse_type == "analytical" || mse_type == "bootstrap"))
    Rcpp::stop("mse_type must be 'analytical' or 'bootstrap'.");
  if (mse_type == "bootstrap" && B < 1) Rcpp::stop("B must be >= 1 for bootstrap MSE.");
  if (area.n_elem != (uword)N) Rcpp::stop("area length must match data rows.");
  if (vardirall.n_elem != (uword)N || yall.n_elem != (uword)N)
    Rcpp::stop("y, vardir and X must have the same number of rows.");
  if (arma::any(area < 0)) Rcpp::stop("area must be a 0-based index (>= 0).");

  uvec idx_s = find_finite(yall);
  uvec idx_ns = find_nonfinite(yall);
  const bool adaNA = idx_ns.n_elem > 0;

  mat Xs = Xall.rows(idx_s);
  vec ys = yall.elem(idx_s);
  vec psis = vardirall.elem(idx_s);
  ivec areas_s = area.elem(idx_s);

  const int Ns = Xs.n_rows;
  if (Ns == 0) Rcpp::stop("No sampled subareas found.");
  if (Ns <= p) Rcpp::stop("Not enough sampled subareas for the model.");
  if (arma::any(psis <= 0)) Rcpp::stop("vardir must be strictly positive for sampled subareas.");

  const int m_all = max(area) + 1;
  std::vector<uvec> sidx(m_all);   // row positions (within the sampled data) per area
  int m_fit = 0;
  for (int d = 0; d < m_all; ++d) {
    sidx[d] = find(areas_s == d);
    if (sidx[d].n_elem > 0) ++m_fit;
  }
  if (m_fit < 2) Rcpp::stop("Two-fold model needs at least 2 areas with sampled subareas.");

  // global row indices of NON-sampled subareas, per area
  std::vector<uvec> nsg(m_all);
  if (adaNA) {
    const ivec areas_ns = area.elem(idx_ns);
    for (int d = 0; d < m_all; ++d) nsg[d] = idx_ns.elem(find(areas_ns == d));
  }

  // ---- main fit (paper notation: s2v=area, s2u=subarea)
  TfhFit f = tfh_fit(Xs, ys, psis, sidx, m_all, method, maxiter, precision);
  const vec& beta = f.beta;
  const double s2v = f.s2v, s2u = f.s2u;
  const double var_v = f.InfoInv(0,0), var_u = f.InfoInv(1,1), cov_vu = f.InfoInv(0,1);

  // ---- EBLUP for sampled; non-sampled: x'beta (+ area effect if area is sampled)
  vec eblup_all(N), mse_all(N);
  vec v_area(N), u_sub(N);          // v = area effect, u = subarea effect (paper)
  eblup_all.fill(datum::nan);
  mse_all.fill(datum::nan);
  v_area.fill(NA_REAL);
  u_sub.fill(NA_REAL);
  vec Xbeta_all = Xall * beta;
  if (adaNA) eblup_all.elem(idx_ns) = Xbeta_all.elem(idx_ns);   // synthetic default

  for (int d = 0; d < m_all; ++d) {
    const uvec& id = sidx[d];
    if (id.n_elem == 0) continue;
    uvec gid = idx_s.elem(id);
    const mat Xd = Xs.rows(id);
    const vec yd = ys.elem(id);
    const vec psid = psis.elem(id);

    // Woodbury quantities
    const vec dd = 1.0 / (s2u + psid);
    const double s = accu(dd);
    const double c = (s2v > 0.0) ? s2v / (1.0 + s2v * s) : 0.0;
    const vec res = yd - Xd * beta;
    const double ddr = dot(dd, res);
    const vec vr = dd % res - c * dd * ddr;   // V_d^{-1} r
    const double vhat = s2v * accu(vr);       // area effect BLUP
    const vec uhat = s2u * vr;                // subarea effect BLUP
    eblup_all.elem(gid) = Xd * beta + vhat + uhat;
    v_area.elem(gid).fill(vhat);
    u_sub.elem(gid) = uhat;

    // non-sampled subareas of a sampled area: x'beta + vhat_i
    for (uword k = 0; k < nsg[d].n_elem; ++k) {
      const uword gi = nsg[d](k);
      eblup_all(gi) = Xbeta_all(gi) + vhat;
      v_area(gi) = vhat;
      u_sub(gi) = 0.0;
    }
  }

  // ---- MSE: analytical (g1+g2+g3) or parametric bootstrap
  if (mse_type == "analytical") {
    // Sampled: g1+g2+g3 via explicit per-area matrices
    for (int d = 0; d < m_all; ++d) {
      const uvec& id = sidx[d];
      if (id.n_elem == 0) continue;
      uvec gid = idx_s.elem(id);
      const int nd = id.n_elem;
      const mat Xd = Xs.rows(id);
      const vec psid = psis.elem(id);

      mat Vd = s2v * ones<mat>(nd, nd) + diagmat(s2u + psid);
      mat Sd = s2v * ones<mat>(nd, nd) + s2u * eye<mat>(nd, nd);
      mat Vinv;
      if (!inv_sympd(Vinv, Vd)) Vinv = pinv(Vd);
      mat Vinv2 = Vinv * Vinv;

      // g1 = diag(Sd - Sd*Vinv*Sd)
      mat SVinv = Sd * Vinv;
      mat G1m = Sd - SVinv * Sd;
      vec g1 = G1m.diag();

      // g2 = diag( D Q D' ), D = (I - Sd*Vinv) Xd
      mat Dd = (eye<mat>(nd, nd) - SVinv) * Xd;
      mat DQ = Dd * f.Q;
      vec g2(nd);
      for (int j = 0; j < nd; ++j) g2(j) = dot(DQ.row(j), Dd.row(j));

      // g3: Mu = Vinv - Sd*Vinv2 ; Mv = 11'Vinv - Sd*Vinv*11'Vinv
      mat Mu = Vinv - Sd * Vinv2;
      vec onev = ones<vec>(nd);
      vec w = Vinv * onev;
      mat Mv = onev * w.t() - Sd * Vinv * onev * w.t();
      mat MuV = Mu * Vd;
      mat MvV = Mv * Vd;
      vec E_uu(nd), E_vv(nd), E_vu(nd);
      for (int j = 0; j < nd; ++j) {
        E_uu(j) = dot(MuV.row(j), Mu.row(j));
        E_vv(j) = dot(MvV.row(j), Mv.row(j));
        E_vu(j) = dot(MvV.row(j), Mu.row(j));
      }
      vec g3 = E_vv * var_v + E_uu * var_u + 2.0 * E_vu * cov_vu;
      vec mse_d = g1 + g2 + g3;
      for (int j = 0; j < nd; ++j) if (mse_d(j) < 0 && mse_d(j) > -1e-8) mse_d(j) = 0;
      mse_all.elem(gid) = mse_d;
    }

    // Non-sampled subareas. Predictor: x'beta + vhat_i (sampled area) or
    // x'beta (area without any sampled subarea).
    if (adaNA) {
      for (int d = 0; d < m_all; ++d) {
        const uvec& ng = nsg[d];
        if (ng.n_elem == 0) continue;

        if (sidx[d].n_elem == 0) {
          // synthetic: error = x'(b^-b) - v_i - u_ij  =>  x'Qx + s2v + s2u
          for (uword k = 0; k < ng.n_elem; ++k) {
            const uword gi = ng(k);
            const vec xj = Xall.row(gi).t();
            double m0 = as_scalar(xj.t() * f.Q * xj) + s2v + s2u;
            mse_all(gi) = m0;
          }
          continue;
        }

        const uvec& id = sidx[d];
        const mat Xd = Xs.rows(id);
        const vec psid = psis.elem(id);
        const int nd = id.n_elem;
        const vec dd = 1.0 / (s2u + psid);
        const double s = accu(dd);
        const double c = (s2v > 0.0) ? s2v / (1.0 + s2v * s) : 0.0;

        const vec w = dd * (1.0 - c * s);            // V^{-1} 1
        const double w1 = accu(w);                   // 1'V^{-1}1
        const vec Xtw = Xd.t() * w;                  // X'V^{-1}1
        const vec Vinv_w = dd % w - c * dd * dot(dd, w);   // V^{-2} 1

        // g1*: Var(v_i^ - v_i) + Var(u_ij) ;  v_i^ = s2v 1'V^{-1}(y - X b)
        const double g1s = s2v - s2v * s2v * w1 + s2u;

        // g3*: derivatives of b = s2v V^{-1}1 w.r.t. s2v and s2u
        const vec a_v = w - s2v * w1 * w;
        const vec a_u = -s2v * Vinv_w;
        const vec ones_nd = ones<vec>(nd);
        const vec Va_v = a_v / dd + s2v * accu(a_v) * ones_nd;   // V a_v
        const vec Va_u = a_u / dd + s2v * accu(a_u) * ones_nd;   // V a_u
        const double Evv = dot(a_v, Va_v);
        const double Euu = dot(a_u, Va_u);
        const double Evu = dot(a_v, Va_u);
        const double g3s = Evv * var_v + Euu * var_u + 2.0 * Evu * cov_vu;

        for (uword k = 0; k < ng.n_elem; ++k) {
          const uword gi = ng(k);
          const vec xj = Xall.row(gi).t();
          // g2*: (x - s2v X'V^{-1}1)' Q (x - s2v X'V^{-1}1)
          const vec av = xj - s2v * Xtw;
          const double g2s = as_scalar(av.t() * f.Q * av);
          double mse_star = g1s + g2s + g3s;
          if (mse_star < 0 && mse_star > -1e-8) mse_star = 0;
          mse_all(gi) = mse_star;
        }
      }
    }
  } else {
    // ---- parametric bootstrap MSE (paper notation: s2v=area, s2u=subarea)
    // OpenMP parallel: pre-generate all random draws (RNG not thread-safe),
    // each thread writes its own column of SqDiff.
    const double sd_v = std::sqrt(s2v), sd_u = std::sqrt(s2u);
    mat V_star = randn<mat>(m_all, B) * sd_v;  // area effects
    mat U_star = randn<mat>(N, B) * sd_u;      // subarea effects
    mat E_star = randn<mat>(Ns, B);
    E_star.each_col() %= sqrt(psis);           // sampling errors
    mat SqDiff(N, B, fill::zeros);

#pragma omp parallel for schedule(static)
    for (int b = 0; b < B; ++b) {
      vec theta_star(N);
      vec v_star = V_star.col(b);
      vec u_star = U_star.col(b);
      for (int i = 0; i < N; ++i)
        theta_star(i) = Xbeta_all(i) + v_star(area(i)) + u_star(i);
      vec y_star = theta_star.elem(idx_s) + E_star.col(b);

      TfhFit fb = tfh_fit(Xs, y_star, psis, sidx, m_all, method, maxiter, precision,
                          false);
      vec eb_star(N);
      vec Xb = Xall * fb.beta;
      const double bs2v = fb.s2v, bs2u = fb.s2u;
      if (adaNA) eb_star.elem(idx_ns) = Xb.elem(idx_ns);   // default: synthetic
      for (int d = 0; d < m_all; ++d) {
        const uvec& id = sidx[d];
        if (id.n_elem == 0) continue;
        uvec gid = idx_s.elem(id);
        const mat Xd = Xs.rows(id);
        const vec psid = psis.elem(id);
        const vec dd = 1.0 / (bs2u + psid);
        const double s = accu(dd);
        const double c = (bs2v > 0.0) ? bs2v / (1.0 + bs2v * s) : 0.0;
        const vec res = y_star.elem(id) - Xd * fb.beta;
        const vec vr = dd % res - c * dd * dot(dd, res);
        const double vh = bs2v * accu(vr);
        eb_star.elem(gid) = Xd * fb.beta + vh + bs2u * vr;
        // same predictor as the point estimate for non-sampled subareas
        for (uword k = 0; k < nsg[d].n_elem; ++k) {
          const uword gi = nsg[d](k);
          eb_star(gi) = Xb(gi) + vh;
        }
      }
      SqDiff.col(b) = square(eb_star - theta_star);
    }
    mse_all = mean(SqDiff, 1);
  }

  // ---- RSE (%)
  vec rse = sqrt(mse_all);
  vec abse = arma::abs(eblup_all);
  for (uword i = 0; i < rse.n_elem; ++i) {
    rse(i) = (abse(i) > 0) ? 100.0 * rse(i) / abse(i) : datum::nan;
  }

  // ---- coef table
  vec stderr_beta = sqrt(f.Q.diag());
  vec zvalue = beta / stderr_beta;
  vec pvalue(zvalue.n_elem);
  for (uword i = 0; i < pvalue.n_elem; ++i)
    pvalue(i) = 2.0 * R::pnorm5(-std::fabs(zvalue(i)), 0.0, 1.0, 1, 0);

  DataFrame df_coef = DataFrame::create(
    _["beta"] = beta, _["std.error"] = stderr_beta,
    _["stderr_beta"] = stderr_beta, _["zvalue"] = zvalue, _["pvalue"] = pvalue);

  DataFrame df_eblup = DataFrame::create(
    _["y"] = yall, _["eblup"] = eblup_all, _["vardir"] = vardirall,
    _["random_effect_area"] = v_area, _["random_effect_subarea"] = u_sub,
    _["mse"] = mse_all, _["rse"] = rse);

  // Paper notation: sigma2_v = area, sigma2_u = subarea
  NumericVector s2 = NumericVector::create(_["sigma2_v"] = s2v, _["sigma2_u"] = s2u);
  NumericVector goodness = NumericVector::create(
    _["loglikelihood"] = f.loglik,
    _["AIC"] = -2.0 * f.loglik + 2.0 * (p + 2),
    _["BIC"] = -2.0 * f.loglik + (p + 2) * std::log((double)m_fit));

  return List::create(
    _["random_effect_var"] = s2, _["estcoef"] = df_coef, _["df_eblup"] = df_eblup,
    _["goodness"] = goodness, _["n_iter"] = f.niter,
    _["convergence"] = f.conv, _["method"] = "eblup", _["level"] = "subarea");
}
