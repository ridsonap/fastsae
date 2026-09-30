// src/eblup_tfh.cpp
// Two-fold Fay-Herriot model (area-level with subareas), Torabi & Rao (2014).
//
//   Level 1 (sampling): y_di = theta_di + e_di,  e_di ~ N(0, psi_di)
//   Level 2 (linking):  theta_di = x_di'beta + u_d + v_di,
//                       u_d ~ N(0, s2u), v_di ~ N(0, s2v)
//
// Variance components (s2u, s2v) by ML/REML Fisher-scoring exploiting the
// block-diagonal covariance (Woodbury per area block). EBLUP from the BLUP
// equations. MSE by parametric bootstrap (same pattern as seblup/stfh pbmse);
// the closed-form Prasad-Rao g3 term of Torabi & Rao (2014) is left for a
// future release.
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

// ---------------------------------------------------------------- per-area
// scratch for one Fisher-scoring iteration
struct AreaScratch {
  vec dd;      // 1 / (s2v + psi_d)
  double s = 0.0, c = 0.0;
  double dd2 = 0.0, dd3 = 0.0;  // sums of dd^2, dd^3
  mat A;       // V_d^{-1} X_d
  vec viny;    // V_d^{-1} y_d
};

struct TfhFit {
  vec beta;
  double s2u = 0.0, s2v = 0.0;
  mat Q;  // (X'V^{-1}X)^{-1}
  int niter = 0;
  bool conv = false;
  double loglik = 0.0;
};

// Fisher-scoring for (s2u, s2v). sidx[d] = row positions (into the sampled
// arrays) of area d; areas without sampled rows are skipped.
TfhFit tfh_fit(const mat& Xs, const vec& ys, const vec& psis,
               const std::vector<uvec>& sidx, int m_all,
               const std::string& method, int maxiter, double precision) {
  const int p = Xs.n_cols;
  const bool is_ml = (method == "ML");
  TfhFit f;

  const double med = as_scalar(median(psis));
  double s2u = 0.5 * med, s2v = 0.5 * med;

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

      a.dd = 1.0 / (s2v + psid);
      a.s = accu(a.dd);
      a.c = (s2u > 0.0) ? s2u / (1.0 + s2u * a.s) : 0.0;
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

    double su = 0.0, sv = 0.0, Iuu = 0.0, Ivv = 0.0, Iuv = 0.0;
    for (int d = 0; d < m_all; ++d) {
      const uvec& id = sidx[d];
      if (id.n_elem == 0) continue;
      const AreaScratch& a = scr[d];

      const vec r = a.viny - a.A * tvec;
      const double sumr = accu(r);
      const double trV = accu(a.dd) - a.c * a.dd2;  // tr(V_d^{-1})
      const double s1 = a.s * (1.0 - a.c * a.s);                    // 1'V_d^{-1}1
      const vec avec = sum(a.A, 0).t();
      const mat M = a.A.t() * a.A;
      const double qd = s1 - as_scalar(avec.t() * Q * avec);  // 1'P_d 1

      if (!is_ml) {
        const double trP = trV - trace(Q * M);
        su += -0.5 * qd + 0.5 * sumr * sumr;
        sv += -0.5 * trP + 0.5 * dot(r, r);
        Iuu += 0.5 * qd * qd;

        mat VinvA = a.A.each_col() % a.dd;
        VinvA -= a.c * a.dd * (a.dd.t() * a.A);
        const mat G = a.A.t() * VinvA;
        const mat MQ = M * Q;
        const double trV2 = a.dd2 - 2.0 * a.c * a.dd3 + a.c * a.c * a.dd2 * a.dd2;
        const double trP2 = trV2 - 2.0 * trace(Q * G) + trace(MQ * MQ);
        Ivv += 0.5 * trP2;

        const vec z = a.dd * (1.0 - a.c * a.s);  // V_d^{-1} 1
        const rowvec zA = z.t() * a.A;
        const double term = dot(z, z) - 2.0 * as_scalar(zA * Q * avec) +
                            as_scalar(avec.t() * Q * M * Q * avec);
        Iuv += 0.5 * term;
      } else {
        const vec e = a.viny - a.A * tvec;
        const double sume = accu(e);
        su += -0.5 * s1 + 0.5 * sume * sume;
        sv += -0.5 * trV + 0.5 * dot(e, e);
        Iuu += 0.5 * s1 * s1;
        Ivv += 0.5 * (a.dd2 - 2.0 * a.c * a.dd3 + a.c * a.c * a.dd2 * a.dd2);
        const vec z = a.dd * (1.0 - a.c * a.s);
        Iuv += 0.5 * dot(z, z);
      }
    }

    Info(0, 0) = Iuu; Info(0, 1) = Iuv;
    Info(1, 0) = Iuv; Info(1, 1) = Ivv;
    sc(0) = su; sc(1) = sv;

    const double det = Iuu * Ivv - Iuv * Iuv;
    if (std::fabs(det) > 1e-300 && std::isfinite(det)) {
      mat InfoInv;
      bool ok2 = inv_sympd(InfoInv, Info);
      if (ok2) {
        delta = InfoInv * sc;
      } else {
        delta(0) = (Iuu > 1e-300) ? su / Iuu : 0.0;
        delta(1) = (Ivv > 1e-300) ? sv / Ivv : 0.0;
      }
    } else {
      delta(0) = (Iuu > 1e-300) ? su / Iuu : 0.0;
      delta(1) = (Ivv > 1e-300) ? sv / Ivv : 0.0;
    }

    double s2u_new = s2u + delta(0);
    double s2v_new = s2v + delta(1);
    if (!std::isfinite(s2u_new)) s2u_new = s2u;
    if (!std::isfinite(s2v_new)) s2v_new = s2v;
    if (s2u_new < 0.0) s2u_new = 0.0;
    if (s2v_new < 0.0) s2v_new = 0.0;

    const double d1 = std::fabs((s2u_new - s2u) / std::max(s2u, 1e-12));
    const double d2 = std::fabs((s2v_new - s2v) / std::max(s2v, 1e-12));
    diff = std::max(d1, d2);
    s2u = s2u_new;
    s2v = s2v_new;
    ++k;
  }

  f.s2u = s2u;
  f.s2v = s2v;
  f.beta = tvec;
  f.Q = Q;
  f.niter = k;
  f.conv = (diff <= precision);

  // log-likelihood at the final estimates
  double ll = 0.0;
  for (int d = 0; d < m_all; ++d) {
    const uvec& id = sidx[d];
    if (id.n_elem == 0) continue;
    const mat Xd = Xs.rows(id);
    const vec yd = ys.elem(id);
    const vec psid = psis.elem(id);
    const vec dd = 1.0 / (s2v + psid);
    const double s = accu(dd);
    const double c = (s2u > 0.0) ? s2u / (1.0 + s2u * s) : 0.0;
    const vec res = yd - Xd * f.beta;
    const vec vr = dd % res - c * dd * dot(dd, res);
    ll += -0.5 * (accu(log(s2v + psid)) + std::log(1.0 + s2u * s) + dot(res, vr));
  }
  f.loglik = ll;
  return f;
}

}  // namespace

// [[Rcpp::export(.eblup_tfh_core)]]
List eblup_tfh_core(const arma::mat& Xall, const arma::vec& yall,
                    const arma::vec& vardirall, const arma::ivec& area,
                    std::string method = "REML", int B = 200,
                    int maxiter = 100, double precision = 1e-4) {
  const int N = Xall.n_rows;
  const int p = Xall.n_cols;

  if (!(method == "ML" || method == "REML")) Rcpp::stop("method must be 'ML' or 'REML'.");
  if (B < 1) Rcpp::stop("B must be >= 1.");
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

  // ---- main fit
  TfhFit f = tfh_fit(Xs, ys, psis, sidx, m_all, method, maxiter, precision);
  const vec& beta = f.beta;
  const double s2u = f.s2u, s2v = f.s2v;

  // ---- EBLUP for all rows
  vec eblup_all(N), u_area(N), v_sub(N);
  u_area.fill(NA_REAL);
  v_sub.fill(NA_REAL);
  vec Xbeta_all = Xall * beta;
  for (int d = 0; d < m_all; ++d) {
    const uvec& id = sidx[d];
    if (id.n_elem == 0) continue;
    // map sampled positions back to global rows
    uvec gid = idx_s.elem(id);
    const mat Xd = Xs.rows(id);
    const vec yd = ys.elem(id);
    const vec psid = psis.elem(id);
    const vec dd = 1.0 / (s2v + psid);
    const double s = accu(dd);
    const double c = (s2u > 0.0) ? s2u / (1.0 + s2u * s) : 0.0;
    const vec res = yd - Xd * beta;
    const vec vr = dd % res - c * dd * dot(dd, res);
    const double uhat = s2u * accu(vr);
    const vec vhat = s2v * vr;
    eblup_all.elem(gid) = Xd * beta + uhat + vhat;
    u_area.elem(gid).fill(uhat);
    v_sub.elem(gid) = vhat;
  }
  if (adaNA) eblup_all.elem(idx_ns) = Xbeta_all.elem(idx_ns);

  // ---- parametric bootstrap MSE
  vec se_acc(N, fill::zeros);
  const double sd_u = std::sqrt(s2u), sd_v = std::sqrt(s2v);
  vec theta_star(N), y_star(Ns);
  for (int b = 0; b < B; ++b) {
    vec u_star = randn<vec>(m_all) * sd_u;
    vec v_star = randn<vec>(N) * sd_v;
    vec e_star = randn<vec>(Ns);
    e_star %= sqrt(psis);
    for (int i = 0; i < N; ++i) theta_star(i) = Xbeta_all(i) + u_star(area(i)) + v_star(i);
    y_star = theta_star.elem(idx_s) + e_star;

    TfhFit fb = tfh_fit(Xs, y_star, psis, sidx, m_all, method, maxiter, precision);
    // EBLUP under bootstrap fit
    vec eb_star(N);
    vec Xb = Xall * fb.beta;
    const double bs2u = fb.s2u, bs2v = fb.s2v;
    for (int d = 0; d < m_all; ++d) {
      const uvec& id = sidx[d];
      if (id.n_elem == 0) continue;
      uvec gid = idx_s.elem(id);
      const mat Xd = Xs.rows(id);
      const vec psid = psis.elem(id);
      const vec dd = 1.0 / (bs2v + psid);
      const double s = accu(dd);
      const double c = (bs2u > 0.0) ? bs2u / (1.0 + bs2u * s) : 0.0;
      const vec res = y_star.elem(id) - Xd * fb.beta;
      const vec vr = dd % res - c * dd * dot(dd, res);
      eb_star.elem(gid) = Xd * fb.beta + bs2u * accu(vr) + bs2v * vr;
    }
    if (adaNA) eb_star.elem(idx_ns) = Xb.elem(idx_ns);
    se_acc += square(eb_star - theta_star);
  }
  vec mse_all = se_acc / (double)B;

  // ---- RSE (%)
  vec rse = sqrt(mse_all);
  vec relerr = rse / arma::abs(eblup_all);
  relerr.elem(arma::find(eblup_all == 0)).fill(datum::nan);
  rse = 100.0 * relerr;

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

  NumericVector s2 = NumericVector::create(_["sigma2_u"] = s2u, _["sigma2_v"] = s2v);
  NumericVector goodness = NumericVector::create(
      _["loglikelihood"] = f.loglik,
      _["AIC"] = -2.0 * f.loglik + 2.0 * (p + 2),
      _["BIC"] = -2.0 * f.loglik + (p + 2) * std::log((double)m_fit));

  return List::create(
      _["random_effect_var"] = s2, _["estcoef"] = df_coef, _["df_eblup"] = df_eblup,
      _["goodness"] = goodness, _["n_iter"] = f.niter,
      _["convergence"] = f.conv, _["method"] = "eblup", _["level"] = "subarea",
      _["B"] = B);
}
