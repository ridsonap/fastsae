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
#include <RcppArmadillo.h>

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

// Fisher-scoring for (s2v, s2u). sidx[d] = row positions of area d.
TfhFit tfh_fit(const mat& Xs, const vec& ys, const vec& psis,
               const std::vector<uvec>& sidx, int m_all,
               const std::string& method, int maxiter, double precision) {
  const int p = Xs.n_cols;
  const bool is_ml = (method == "ML");
  TfhFit f;

  const double med = as_scalar(median(psis));
  double s2v = 0.5 * med, s2u = 0.5 * med;

  mat XtVX(p, p), eye_p = eye<mat>(p, p), Q(p, p), Info(2, 2);
  vec Xty(p), tvec(p), sc(2), delta(2);

  double diff = precision + 1.0;
  int k = 0;
  while ((diff > precision) && (k < maxiter)) {
    XtVX.zeros();
    Xty.zeros();
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

    bool ok = inv_sympd(Q, XtVX);
    if (!ok) Q = solve(XtVX, eye_p);
    tvec = Q * Xty;

    // Score and information for (s2v, s2u). Order: 0=v (area), 1=u (subarea).
    double sv = 0.0, su = 0.0, Ivv = 0.0, Iuu = 0.0, Ivu = 0.0;
    for (int d = 0; d < m_all; ++d) {
      const uvec& id = sidx[d];
      if (id.n_elem == 0) continue;
      const AreaScratch& a = scr[d];

      const vec r = a.viny - a.A * tvec;
      const double sumr = accu(r);
      const double trV = accu(a.dd) - a.c * a.dd2;
      const double s1 = a.s * (1.0 - a.c * a.s);
      const vec avec = sum(a.A, 0).t();
      const mat M = a.A.t() * a.A;
      const double qd = s1 - as_scalar(avec.t() * Q * avec);

      if (!is_ml) {
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

        const vec z = a.dd * (1.0 - a.c * a.s);
        const rowvec zA = z.t() * a.A;
        const double term = dot(z, z) - 2.0 * as_scalar(zA * Q * avec) +
                            as_scalar(avec.t() * Q * M * Q * avec);
        Ivu += 0.5 * term;
      } else {
        const vec e = a.viny - a.A * tvec;
        const double sume = accu(e);
        sv += -0.5 * s1 + 0.5 * sume * sume;
        su += -0.5 * trV + 0.5 * dot(e, e);
        Ivv += 0.5 * s1 * s1;
        Iuu += 0.5 * (a.dd2 - 2.0 * a.c * a.dd3 + a.c * a.c * a.dd2 * a.dd2);
        const vec z = a.dd * (1.0 - a.c * a.s);
        Ivu += 0.5 * dot(z, z);
      }
    }

    Info(0, 0) = Ivv; Info(0, 1) = Ivu;
    Info(1, 0) = Ivu; Info(1, 1) = Iuu;
    sc(0) = sv; sc(1) = su;

    const double det = Ivv * Iuu - Ivu * Ivu;
    if (std::fabs(det) > 1e-300 && std::isfinite(det)) {
      mat InfoInv;
      bool ok2 = inv_sympd(InfoInv, Info);
      if (ok2) {
        delta = InfoInv * sc;
      } else {
        delta(0) = (Ivv > 1e-300) ? sv / Ivv : 0.0;
        delta(1) = (Iuu > 1e-300) ? su / Iuu : 0.0;
      }
    } else {
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

  f.s2v = s2v;
  f.s2u = s2u;
  f.beta = tvec;
  f.Q = Q;
  f.niter = k;
  f.conv = (diff <= precision);

  // Inverse information at final estimates (for g3 variances)
  {
    mat InfoF(2, 2);
    double Ivv = 0.0, Iuu = 0.0, Ivu = 0.0;
    // Recompute information at final (s2v, s2u) using ML-type formulas
    // (simpler; REML adjustment is second-order for g3)
    for (int d = 0; d < m_all; ++d) {
      const uvec& id = sidx[d];
      if (id.n_elem == 0) continue;
      const vec psid = psis.elem(id);
      const vec dd = 1.0 / (s2u + psid);
      const double s = accu(dd);
      const double c = (s2v > 0.0) ? s2v / (1.0 + s2v * s) : 0.0;
      const double dd2 = accu(square(dd));
      const double dd3 = accu(dd % square(dd));
      const double s1 = s * (1.0 - c * s);
      Ivv += 0.5 * s1 * s1;
      Iuu += 0.5 * (dd2 - 2.0 * c * dd3 + c * c * dd2 * dd2);
      const vec z = dd * (1.0 - c * s);
      Ivu += 0.5 * dot(z, z);
    }
    InfoF(0,0) = Ivv; InfoF(0,1) = Ivu; InfoF(1,0) = Ivu; InfoF(1,1) = Iuu;
    mat II;
    if (inv_sympd(II, InfoF)) f.InfoInv = II;
    else f.InfoInv = zeros<mat>(2,2);
  }

  // log-likelihood
  double ll = 0.0;
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
    ll += -0.5 * (accu(log(s2u + psid)) + std::log(1.0 + s2v * s) + dot(res, vr));
  }
  f.loglik = ll;
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
  std::vector<uvec> sidx(m_all);
  int m_fit = 0;
  for (int d = 0; d < m_all; ++d) {
    sidx[d] = find(areas_s == d);
    if (sidx[d].n_elem > 0) ++m_fit;
  }
  if (m_fit < 2) Rcpp::stop("Two-fold model needs at least 2 areas with sampled subareas.");

  // ---- main fit (paper notation: s2v=area, s2u=subarea)
  TfhFit f = tfh_fit(Xs, ys, psis, sidx, m_all, method, maxiter, precision);
  const vec& beta = f.beta;
  const double s2v = f.s2v, s2u = f.s2u;
  const double var_v = f.InfoInv(0,0), var_u = f.InfoInv(1,1), cov_vu = f.InfoInv(0,1);

  // ---- EBLUP for sampled; synthetic for non-sampled
  // MSE computed separately below (analytical or bootstrap)
  vec eblup_all(N), mse_all(N);
  vec u_area(N), v_sub(N);
  u_area.fill(NA_REAL);
  v_sub.fill(NA_REAL);
  vec Xbeta_all = Xall * beta;

  for (int d = 0; d < m_all; ++d) {
    const uvec& id = sidx[d];
    if (id.n_elem == 0) continue;
    uvec gid = idx_s.elem(id);
    const int nd = id.n_elem;
    const mat Xd = Xs.rows(id);
    const vec yd = ys.elem(id);
    const vec psid = psis.elem(id);

    // Woodbury quantities
    const vec dd = 1.0 / (s2u + psid);
    const double s = accu(dd);
    const double c = (s2v > 0.0) ? s2v / (1.0 + s2v * s) : 0.0;
    const vec res = yd - Xd * beta;
    const double ddr = dot(dd, res);
    const vec vr = dd % res - c * dd * ddr;  // V_d^{-1} r
    const double vhat = s2v * accu(vr);       // area effect BLUP
    const vec uhat = s2u * vr;                // subarea effect BLUP
    eblup_all.elem(gid) = Xd * beta + vhat + uhat;
    u_area.elem(gid).fill(vhat);
    v_sub.elem(gid) = uhat;
  }
  // Non-sampled EBLUP: synthetic
  if (adaNA) eblup_all.elem(idx_ns) = Xbeta_all.elem(idx_ns);

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

    // Non-sampled: g1*+g2*+g3* via delta method
    if (adaNA) {
      for (int d = 0; d < m_all; ++d) {
        const uvec& id = sidx[d];
        if (id.n_elem == 0) continue;
        const int nd = id.n_elem;
        const mat Xd = Xs.rows(id);
        const vec psid = psis.elem(id);
        const vec dd = 1.0 / (s2u + psid);
        const double s = accu(dd);
        const double c = (s2v > 0.0) ? s2v / (1.0 + s2v * s) : 0.0;
        const vec w = dd * (1.0 - c * s);
        const vec Xtw = Xd.t() * w;
        uvec alld = find(area == d);
        for (uword k = 0; k < alld.n_elem; ++k) {
          const int gi = alld(k);
          if (arma::is_finite(yall(gi))) continue;
          rowvec xj = Xall.row(gi);
          const vec yd = ys.elem(id);
          const vec res = yd - Xd * beta;
          double vhat_i = s2v * dot(w, res);
          eblup_all(gi) = as_scalar(xj * beta) + vhat_i;
          u_area(gi) = vhat_i;
          v_sub(gi) = 0.0;
          double xQx = as_scalar(xj * f.Q * xj.t());
          double xQXt = as_scalar(xj * f.Q * Xtw);
          double mse_star = xQx + s2v + s2u - 2.0 * s2v * xQXt;
          mat Vd = s2v * ones<mat>(nd, nd) + diagmat(s2u + psid);
          mat Vinv;
          if (!inv_sympd(Vinv, Vd)) Vinv = pinv(Vd);
          vec wv = Vinv * ones<vec>(nd);
          double w1 = dot(wv, ones<vec>(nd));
          vec a_v = wv - s2v * w1 * wv;
          vec a_u = -s2v * (Vinv * wv);
          double Evv = as_scalar(a_v.t() * Vd * a_v);
          double Euu = as_scalar(a_u.t() * Vd * a_u);
          double Evu = as_scalar(a_v.t() * Vd * a_u);
          double g3s = Evv * var_v + Euu * var_u + 2.0 * Evu * cov_vu;
          mse_star += g3s;
          if (mse_star < 0 && mse_star > -1e-8) mse_star = 0;
          mse_all(gi) = mse_star;
        }
      }
    }
  } else {
    // ---- parametric bootstrap MSE (paper notation: s2v=area, s2u=subarea)
    vec se_acc(N, fill::zeros);
    const double sd_v = std::sqrt(s2v), sd_u = std::sqrt(s2u);
    vec theta_star(N), y_star(Ns);
    for (int b = 0; b < B; ++b) {
      vec v_star = randn<vec>(m_all) * sd_v;  // area effects
      vec u_star = randn<vec>(N) * sd_u;      // subarea effects
      vec e_star = randn<vec>(Ns);
      e_star %= sqrt(psis);
      for (int i = 0; i < N; ++i) theta_star(i) = Xbeta_all(i) + v_star(area(i)) + u_star(i);
      y_star = theta_star.elem(idx_s) + e_star;

      TfhFit fb = tfh_fit(Xs, y_star, psis, sidx, m_all, method, maxiter, precision);
      vec eb_star(N);
      vec Xb = Xall * fb.beta;
      const double bs2v = fb.s2v, bs2u = fb.s2u;
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
        eb_star.elem(gid) = Xd * fb.beta + bs2v * accu(vr) + bs2u * vr;
      }
      if (adaNA) eb_star.elem(idx_ns) = Xb.elem(idx_ns);
      se_acc += square(eb_star - theta_star);
    }
    mse_all = se_acc / (double)B;
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
      _["random_effect_area"] = u_area, _["random_effect_subarea"] = v_sub,
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
