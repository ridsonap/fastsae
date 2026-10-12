// src/eblup_sfh_core.cpp
// Core Spatial Fay-Herriot EBLUP estimation with unified ML/REML loops
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "eblup_sfh.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// ============================================================================
// sar_marginal()
//
// SAR random effects are defined on ALL D domains of the (full) weight
// matrix Wf: v = (I - rho Wf)^{-1} u, Var(v) = s2 * Sf, Sf = [(I-rho Wf')(I-rho Wf)]^{-1}.
// The sampled domains s therefore have covariance s2 * Sf[s,s] -- which is
// NOT the SAR covariance built from the sub-matrix Wf[s,s]. All derivatives
// are taken on the full matrix and then restricted to s:
//   Sig   = Sf[s,s]
//   SCS   = (Sf C Sf)[s,s],          C = dA/drho = 2 rho Wf'Wf - Wf - Wf'
//   SCSCS = (Sf C Sf C Sf)[s,s]      (for g4)
//   SWWS  = (Sf Wf'Wf Sf)[s,s]       (for g4)
// With all domains sampled (all = true) this is exactly the classical SFH.
// Pure Armadillo (thread-safe).
// ============================================================================
bool sar_marginal(double rho, const arma::mat& Wf, const arma::mat& WtWf,
                  const arma::uvec& s, bool all,
                  arma::mat& Sig, arma::mat& SCS,
                  arma::mat* SCSCS, arma::mat* SWWS, arma::mat* Sfull) {
  const arma::uword D = Wf.n_rows;
  const arma::mat I = arma::eye<arma::mat>(D, D);
  // symmetric in exact arithmetic; symmetrise the round-off
  const arma::mat A = arma::symmatu((I - rho * Wf.t()) * (I - rho * Wf));
  arma::mat Sf;
  bool ok = arma::inv_sympd(Sf, A);
  if (!ok) Sf = arma::pinv(A);
  const arma::mat C = 2 * rho * WtWf - (Wf + Wf.t());
  const arma::mat CS = C * Sf;
  const arma::mat SCSf = Sf * CS;
  Sig = all ? Sf : arma::mat(Sf.submat(s, s));
  SCS = all ? SCSf : arma::mat(SCSf.submat(s, s));
  if (SCSCS) {
    arma::mat t = SCSf * CS;
    *SCSCS = all ? t : arma::mat(t.submat(s, s));
  }
  if (SWWS) {
    arma::mat t = Sf * WtWf * Sf;
    *SWWS = all ? t : arma::mat(t.submat(s, s));
  }
  if (Sfull) *Sfull = Sf;
  return ok;
}



// ============================================================================
// sfh_lower_rho()
//
// sae::eblupSFH maps an estimate at the lower bound -0.999 to rho = -1. That
// is only valid when (I + W')(I + W) is non-singular, i.e. -1 is not an
// eigenvalue of W -- which fails for bipartite graphs such as a rook grid.
// Map to -1 only when it is valid; otherwise keep -0.999.
// ============================================================================
static double sfh_lower_rho(double rho, const arma::mat& W) {
  if (rho > -0.999) return rho;
  const arma::uword D = W.n_rows;
  const arma::mat I = arma::eye<arma::mat>(D, D);
  const arma::mat A = (I + W.t()) * (I + W);
  return (arma::rcond(A) > 1e-10) ? -1.0 : -0.999;
}

// ============================================================================
// sfh_eval()
//
// (Restricted) log-likelihood of the SFH model at th = (sigma2, rho), up to a
// constant, with optionally the score and Fisher information:
//   REML: -1/2 [ log|V| + log|X'V^{-1}X| + y'Py ]
//   ML  : -1/2 [ log|V| + (y - X b)'V^{-1}(y - X b) ]   (b = GLS estimate)
// V = sigma2 * Sigma_ss(rho) + Psi (Sigma_ss from sar_marginal). Pure Armadillo.
// ============================================================================
static bool sfh_eval(const arma::vec& th, const arma::mat& X, const arma::vec& y,
                     const arma::vec& vardir, const arma::mat& W, const arma::mat& WtW,
                     const arma::uvec& idx_s, bool all_s, bool isML, bool want_SF,
                     double& ll, arma::vec& S, arma::mat& F) {
  arma::mat Sig, SCS;
  if (!sar_marginal(th(1), W, WtW, idx_s, all_s, Sig, SCS, nullptr, nullptr, nullptr)) return false;
  const arma::mat V = arma::symmatu(th(0) * Sig + arma::diagmat(vardir));
  arma::mat R;
  if (!arma::chol(R, V)) return false;
  const arma::mat Rinv = arma::inv(arma::trimatu(R));
  const arma::mat Vi = Rinv * Rinv.t();
  const double logdetV = 2.0 * arma::accu(arma::log(R.diag()));
  const arma::mat XtVi = X.t() * Vi;
  const arma::mat XtViX = arma::symmatu(XtVi * X);
  arma::mat Q;
  if (!arma::inv_sympd(Q, XtViX)) return false;
  const arma::mat P = Vi - XtVi.t() * Q * XtVi;
  const arma::vec Py = P * y;
  if (isML) {
    ll = -0.5 * (logdetV + arma::dot(y, Py));
  } else {
    double ld, sg;
    if (!arma::log_det(ld, sg, XtViX) || sg <= 0) return false;
    ll = -0.5 * (logdetV + ld + arma::dot(y, Py));
  }
  if (!std::isfinite(ll)) return false;
  if (!want_SF) return true;

  const arma::mat Vs = Sig, Vr = -th(0) * SCS;   // dV/dsigma2, dV/drho
  const arma::mat& M = isML ? Vi : P;             // ML: V^{-1}; REML: P
  const arma::mat MS = M * Vs, MR = M * Vr;
  S.set_size(2);
  F.set_size(2, 2);
  S(0) = -0.5 * arma::trace(MS) + 0.5 * arma::as_scalar(Py.t() * Vs * Py);
  S(1) = -0.5 * arma::trace(MR) + 0.5 * arma::as_scalar(Py.t() * Vr * Py);
  F(0, 0) = 0.5 * arma::accu(MS % MS.t());
  F(0, 1) = F(1, 0) = 0.5 * arma::accu(MS % MR.t());
  F(1, 1) = 0.5 * arma::accu(MR % MR.t());
  return S.is_finite() && F.is_finite();
}

// ============================================================================
// sfh_fit_safeguarded()
//
// Constrained (RE)ML for the SFH model: sigma2 >= 0, |rho| <= 0.999, by
// projected Fisher scoring with an active set and step halving (same scheme
// as stfh_fit_safeguarded). Used only when the classical, sae-compatible
// iteration fails or ends with sigma2 < 0 (which sae simply truncates to 0,
// leaving rho and beta at values that are not optimal for sigma2 = 0).
// ============================================================================
static bool sfh_fit_safeguarded(const arma::mat& X, const arma::vec& y, const arma::vec& vardir,
                         const arma::mat& W, const arma::mat& WtW,
                         const arma::uvec& idx_s, bool all_s, bool isML,
                         int maxiter, double precision, arma::vec& th, int& iters) {
  auto project = [](arma::vec t) {
    t(0) = std::max(t(0), 0.0);
    t(1) = std::min(std::max(t(1), -0.999), 0.999);
    return t;
  };
  th = project(th);
  double ll;
  arma::vec S;
  arma::mat F;
  if (!sfh_eval(th, X, y, vardir, W, WtW, idx_s, all_s, isML, true, ll, S, F)) return false;

  const int max_it = 5 * maxiter;
  for (iters = 1; iters <= max_it; ++iters) {
    // active set: sigma2 at 0 with outward score; rho unidentified when
    // sigma2 = 0; rho at a bound with outward score
    std::vector<arma::uword> fr;
    if (!(th(0) == 0.0 && S(0) <= 0.0)) fr.push_back(0);
    if (!(th(0) == 0.0) && !((th(1) >= 0.999 && S(1) >= 0.0) || (th(1) <= -0.999 && S(1) <= 0.0)))
      fr.push_back(1);
    if (fr.empty()) return true;
    const arma::uvec f = arma::conv_to<arma::uvec>::from(fr);

    // converged when the free scores are negligible in SE units
    double zmax = 0.0;
    for (arma::uword a : fr) if (F(a, a) > 0.0) zmax = std::max(zmax, std::fabs(S(a)) / std::sqrt(F(a, a)));
    if (zmax < 10.0 * precision) return true;

    std::vector<arma::vec> dirs;
    {
      arma::vec sf, step(2, fill::zeros);
      if (arma::solve(sf, arma::mat(F.submat(f, f)), arma::vec(S.elem(f)),
                      arma::solve_opts::no_approx)) {
        step.elem(f) = sf;
        dirs.push_back(step);
      }
      arma::vec g(2, fill::zeros);
      for (arma::uword a : fr) g(a) = S(a) / std::max(F(a, a), 1e-12);
      dirs.push_back(g);
    }

    double t = 1.0, llc = 0.0;
    arma::vec cand, Sc;
    arma::mat Fc;
    bool accepted = false;
    for (const arma::vec& step : dirs) {
      t = 1.0;
      for (int h = 0; h < 30; ++h) {
        cand = project(th + t * step);
        if (sfh_eval(cand, X, y, vardir, W, WtW, idx_s, all_s, isML, true, llc, Sc, Fc) && llc > ll) {
          accepted = true;
          break;
        }
        t *= 0.5;
      }
      if (accepted) break;
    }
    if (!accepted) return true;  // no ascent direction left: constrained maximum
    th = cand; ll = llc; S = Sc; F = Fc;
  }
  return false;
}

// ============================================================================
// seblup_fit_core_arma()
//
// Murni Armadillo: TIDAK membuat objek R dan TIDAK memanggil R RNG/Math.
// Aman dipanggil dari OpenMP parallel region.
// ============================================================================
SeblupFitArma seblup_fit_core_arma(
    const arma::mat& X,
    const arma::vec& y,
    const arma::vec& vardir,
    const arma::mat& W,
    const arma::uvec& idx_s,
    const std::string& method,
    int maxiter,
    double precision
) {
  const int m = X.n_rows;
  const int p = X.n_cols;

  // W is the FULL weight matrix; idx_s are the sampled rows (see sar_marginal)
  const bool all_s = (idx_s.n_elem == W.n_rows);
  const mat WtW = W.t() * W;
  mat eye_p = eye<mat>(p, p);

  // initial values
  vec sigma2_u(maxiter + 1, fill::zeros);
  vec rho(maxiter + 1, fill::zeros);
  sigma2_u(0) = median(vardir);
  rho(0) = 0.5;

  // pre-allocate temporaries
  vec s(2, fill::zeros);
  mat Idev(2, 2, fill::zeros);
  vec stime_fin(2, fill::zeros);
  mat Vi(m, m), derSigma(m, m);
  mat derVRho(m, m), V(m, m);
  mat XtVi(p, m), XtViX(p, p), Q(p, p);
  vec Py(m);
  mat P(m, m);

  double diff = precision + 1.0;
  int k = 0;
  const bool isML = (method == "ML");

  // ================================================================
  // Fisher-Scoring Iteration (unified ML/REML)
  // ================================================================
  while ((diff > precision) && (k < maxiter)) {
    ++k;

    // ----- spatial covariance structure -----
    mat SCS;
    sar_marginal(rho(k - 1), W, WtW, idx_s, all_s, derSigma, SCS,
                 nullptr, nullptr, nullptr);
    derVRho = -sigma2_u(k - 1) * SCS;
    bool ok_inv;

    V = symmatu(sigma2_u(k - 1) * derSigma + diagmat(vardir));
    ok_inv = inv_sympd(Vi, V);
    if (!ok_inv) Vi = pinv(V);

    // ----- weighted least squares -----
    XtVi = X.t() * Vi;
    XtViX = XtVi * X;

    ok_inv = inv_sympd(Q, XtViX);
    if (!ok_inv) Q = solve(XtViX, eye_p);

    mat QXtVi = Q * XtVi;  // p x m
    P = Vi - trans(XtVi) * QXtVi;
    Py = P * y;

    // ----- derivatives -----
    mat PD = P * derSigma;
    mat PR = P * derVRho;

    // ----- score and information -----
    if (isML) {
      mat VID = Vi * derSigma;
      mat VIR = Vi * derVRho;

      s(0) = -0.5 * trace(VID) + 0.5 * as_scalar(y.t() * PD * Py);
      s(1) = -0.5 * trace(VIR) + 0.5 * as_scalar(y.t() * PR * Py);

      Idev(0, 0) = 0.5 * trace(VID * VID);
      Idev(0, 1) = 0.5 * trace(VID * VIR);
      Idev(1, 0) = Idev(0, 1);
      Idev(1, 1) = 0.5 * trace(VIR * VIR);
    } else {
      s(0) = -0.5 * trace(PD) + 0.5 * as_scalar(y.t() * PD * Py);
      s(1) = -0.5 * trace(PR) + 0.5 * as_scalar(y.t() * PR * Py);

      Idev(0, 0) = 0.5 * trace(PD * PD);
      Idev(0, 1) = 0.5 * trace(PD * PR);
      Idev(1, 0) = Idev(0, 1);
      Idev(1, 1) = 0.5 * trace(PR * PR);
    }

    // ----- parameter update -----
    vec par_stim(2);
    par_stim(0) = sigma2_u(k - 1);
    par_stim(1) = rho(k - 1);

    vec step;
    bool ok_solve = solve(step, Idev, s);
    if (!ok_solve || !step.is_finite()) {
      k = maxiter;
      break;
    }

    stime_fin = par_stim + step;

    // bound rho ke (-0.999, 0.999)
    if (stime_fin(1) <= -0.999) stime_fin(1) = -0.999;
    if (stime_fin(1) >= 0.999)  stime_fin(1) = 0.999;

    sigma2_u(k) = stime_fin(0);
    rho(k) = stime_fin(1);

    diff = max(abs(stime_fin - par_stim) / abs(par_stim + 1e-12));
    if (!std::isfinite(diff)) {
      k = maxiter;
      break;
    }
  }

  // ================================================================
  // Final estimates
  // ================================================================
  // Constrained fallback, as in seblup_core().
  bool conv = (k < maxiter) && (diff <= precision);
  double s2_hat = sigma2_u(k), rho_hat = rho(k);
  if (!conv || s2_hat < 0.0) {
    arma::vec th(2);
    th(0) = conv ? s2_hat : median(vardir);
    th(1) = conv ? rho_hat : 0.5;
    int it2 = 0;
    conv = sfh_fit_safeguarded(X, y, vardir, W, WtW, idx_s, all_s, isML,
                               maxiter, precision, th, it2);
    k += it2;
    s2_hat = th(0);
    rho_hat = th(1);
  }

  // Same boundary handling as seblup_core().
  double rho_fix = sfh_lower_rho(rho_hat, W);
  if (rho_fix >= 0.999)  rho_fix = 0.999;

  double sigma2 = s2_hat;
  if (sigma2 < 0.0) sigma2 = 0.0;

  // spatial covariance
  mat SCS_f;
  bool ok_inv = sar_marginal(rho_fix, W, WtW, idx_s, all_s, derSigma, SCS_f,
                             nullptr, nullptr, nullptr);

  mat G = sigma2 * derSigma;
  V = symmatu(G + diagmat(vardir));
  ok_inv = inv_sympd(Vi, V) && ok_inv;
  if (!ok_inv) Vi = pinv(V);

  // beta estimate
  XtVi = X.t() * Vi;
  XtViX = XtVi * X;
  mat XtViy = XtVi * y;

  ok_inv = inv_sympd(Q, XtViX) && ok_inv;
  if (!ok_inv) Q = solve(XtViX, eye_p);

  mat beta = Q * XtViy;

  // EBLUP
  vec Xbeta = X * beta;
  vec resid = y - Xbeta;
  mat GVi = G * Vi;
  vec eblup = Xbeta + GVi * resid;

  // MSE components
  mat Ga = G - GVi * G;
  mat Gb = GVi * X;
  mat R = X - Gb;
  vec g1d = Ga.diag();
  vec g2d = sum((R * Q) % R, 1);

  // return
  // Validity must not depend on the scale of y (no absolute thresholds).
  bool ok_all = ok_inv && beta.is_finite() && eblup.is_finite() && g1d.is_finite() && g2d.is_finite() &&
                arma::all(g1d >= 0.0);

  SeblupFitArma out;
  out.theta     = eblup;
  out.g1d       = g1d;
  out.g2d       = g2d;
  out.sigma2    = sigma2;
  out.rho_fix   = rho_fix;
  out.n_iter    = k;
  out.converged = conv && ok_all;
  return out;
}

// ============================================================================
// set_r_seed()
//
// Mengatur seed RNG R dari dalam C++.
// Selalu dipanggil SEBELUM tahap pembangkitan angka acak (sekuensial).
// ============================================================================
void set_r_seed(int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

// [[Rcpp::export(.seblup_core)]]
List seblup_core(
    const arma::mat& Xall,
    const arma::vec& yall,
    const arma::vec& vardirall,
    const arma::mat& Wall,
    std::string method = "REML",
    int maxiter = 100,
    double precision = 1e-4,
    bool only_core = false
) {
  // ================================================================
  // Validasi input
  // ================================================================
  if (Wall.n_rows != Xall.n_rows || Wall.n_cols != Xall.n_rows) {
    Rcpp::stop("'Wall' must be a square matrix with dimension equal to nrow(Xall).");
  }

  // ================================================================
  // hasil
  // ================================================================
  arma::vec eblup_all(Xall.n_rows);
  arma::vec mse_all(Xall.n_rows);

  // ================================================================
  // Handle NA values (unsampled areas)
  // ================================================================
  arma::mat Xns, X;
  arma::vec y, vardir;
  arma::uvec idx_ns, idx_s;
  arma::vec u_ns; // kriged random effects for unsampled areas (filled when adaNA)
  bool adaNA = yall.has_nan();

  if (adaNA) {
    idx_ns = arma::find_nonfinite(yall);
    idx_s  = arma::find_finite(yall);

    Xns    = Xall.rows(idx_ns);
    X      = Xall.rows(idx_s);
    y      = yall.elem(idx_s);
    vardir = vardirall.elem(idx_s);
  } else {
    X      = Xall;
    y      = yall;
    vardir = vardirall;
  }

  // The SAR model is defined on ALL domains of Wall; the sampled domains
  // have the marginal covariance s2 * Sigma_full[s,s] (see sar_marginal).
  // (Previously the model was fitted with the sub-matrix Wall[s,s], i.e. a
  // different SAR model than the one used to krige the unsampled domains.)
  const mat& W = Wall;
  if (!adaNA) idx_s = arma::regspace<arma::uvec>(0, Xall.n_rows - 1);
  const bool all_s = !adaNA;

  const int m = X.n_rows;
  const int p = X.n_cols;

  // safety checks
  if (m == 0) {
    Rcpp::stop("All areas are unsampled -- cannot fit the model.");
  }
  if ((int) vardir.n_elem != m) {
    Rcpp::stop("Length of 'vardir' must equal nrow(X).");
  }
  if ((int) y.n_elem != m) {
    Rcpp::stop("Length of 'y' must equal nrow(X).");
  }
  if (!(method == "ML" || method == "REML")) {
    Rcpp::stop("method must be 'ML' or 'REML'.");
  }

  // ================================================================
  // Precomputations & preallocations
  // ================================================================
  mat Xt = X.t();
  mat eye_p = eye<mat>(p, p);
  const mat WtW = W.t() * W;

  vec sigma2_u(maxiter + 1, fill::zeros);
  vec rho(maxiter + 1, fill::zeros);
  sigma2_u(0) = median(vardir);
  rho(0) = 0.5;

  vec s(2, fill::zeros);
  mat Idev(2, 2, fill::zeros);
  vec stime_fin(2, fill::zeros);
  mat Vi(m, m), derSigma(m, m);
  mat derVRho(m, m), V(m, m);
  mat XtVi(p, m), XtViX(p, p), Q(p, p);
  vec Py(m);
  mat P(m, m);

  double diff = precision + 1.0;
  int k = 0;
  const bool isML = (method == "ML");

  // ================================================================
  // Fisher-Scoring Iteration (unified ML/REML)
  // ================================================================
  while ((diff > precision) && (k < maxiter)) {
    ++k;

    // spatial covariance structure
    mat SCS;
    sar_marginal(rho(k - 1), W, WtW, idx_s, all_s, derSigma, SCS,
                 nullptr, nullptr, nullptr);
    derVRho = -sigma2_u(k - 1) * SCS;
    bool ok_inv;

    V = symmatu(sigma2_u(k - 1) * derSigma + diagmat(vardir));
    ok_inv = inv_sympd(Vi, V);
    if (!ok_inv) Vi = pinv(V);

    // weighted least squares
    XtVi = X.t() * Vi;
    XtViX = XtVi * X;

    ok_inv = inv_sympd(Q, XtViX);
    if (!ok_inv) Q = solve(XtViX, eye_p);

    mat QXtVi = Q * XtVi;
    P = Vi - trans(XtVi) * QXtVi;
    Py = P * y;

    // derivatives
    mat PD = P * derSigma;
    mat PR = P * derVRho;

    // score and information
    if (isML) {
      mat VID = Vi * derSigma;
      mat VIR = Vi * derVRho;

      s(0) = -0.5 * trace(VID) + 0.5 * as_scalar(y.t() * PD * Py);
      s(1) = -0.5 * trace(VIR) + 0.5 * as_scalar(y.t() * PR * Py);

      Idev(0, 0) = 0.5 * trace(VID * VID);
      Idev(0, 1) = 0.5 * trace(VID * VIR);
      Idev(1, 0) = Idev(0, 1);
      Idev(1, 1) = 0.5 * trace(VIR * VIR);
    } else {
      s(0) = -0.5 * trace(PD) + 0.5 * as_scalar(y.t() * PD * Py);
      s(1) = -0.5 * trace(PR) + 0.5 * as_scalar(y.t() * PR * Py);

      Idev(0, 0) = 0.5 * trace(PD * PD);
      Idev(0, 1) = 0.5 * trace(PD * PR);
      Idev(1, 0) = Idev(0, 1);
      Idev(1, 1) = 0.5 * trace(PR * PR);
    }

    // parameter update
    vec par_stim(2);
    par_stim(0) = sigma2_u(k - 1);
    par_stim(1) = rho(k - 1);

    vec step;
    bool ok_solve = solve(step, Idev, s);
    if (!ok_solve || !step.is_finite()) {
      k = maxiter;
      break;
    }

    stime_fin = par_stim + step;

    if (stime_fin(1) <= -0.999) stime_fin(1) = -0.999;
    if (stime_fin(1) >= 0.999)  stime_fin(1) = 0.999;

    sigma2_u(k) = stime_fin(0);
    rho(k) = stime_fin(1);

    diff = max(abs(stime_fin - par_stim) / abs(par_stim + 1e-12));
    if (!std::isfinite(diff)) {
      k = maxiter;
      break;
    }
  }

  // ================================================================
  // Final estimates
  // ================================================================
  // If the classical (sae-compatible) iteration failed, or converged to a
  // NEGATIVE sigma2 (which sae truncates to 0 without re-optimising rho),
  // maximise the (RE)ML likelihood under the constraints instead. Interior
  // solutions of the classical iteration are unchanged.
  bool is_converged = (k < maxiter) && (diff <= precision);
  double s2_hat = sigma2_u(k), rho_hat = rho(k);
  if (!is_converged || s2_hat < 0.0) {
    arma::vec th(2);
    th(0) = is_converged ? s2_hat : median(vardir);
    th(1) = is_converged ? rho_hat : 0.5;
    int it2 = 0;
    is_converged = sfh_fit_safeguarded(X, y, vardir, W, WtW, idx_s, all_s, isML,
                                       maxiter, precision, th, it2);
    k += it2;
    s2_hat = th(0);
    rho_hat = th(1);
  }

  // Lower bound: -0.999 -> -1 as sae::eblupSFH, when rho = -1 is valid
  // (see sfh_lower_rho). Upper bound: never mapped to +1 (I - W singular).
  double rho_fix = sfh_lower_rho(rho_hat, W);

  double sigma2 = s2_hat;
  if (sigma2 < 0.0) sigma2 = 0.0;

  // spatial covariance (+ rho-derivative pieces for g3/g4 and the full
  // Sigma for kriging the unsampled domains)
  mat der3, SCSCS_f, SWWS_f, Sigma_full;
  bool ok_inv = sar_marginal(rho_fix, W, WtW, idx_s, all_s, derSigma, der3,
                             &SCSCS_f, &SWWS_f, &Sigma_full);

  mat G = sigma2 * derSigma;
  V = symmatu(G + diagmat(vardir));
  ok_inv = inv_sympd(Vi, V) && ok_inv;
  if (!ok_inv) Vi = pinv(V);

  // beta
  XtVi = X.t() * Vi;
  XtViX = XtVi * X;
  mat XtViy = XtVi * y;

  ok_inv = inv_sympd(Q, XtViX) && ok_inv;
  if (!ok_inv) Q = solve(XtViX, eye_p);

  mat beta = Q * XtViy;

  // only_core mode: return lightweight results
  if (only_core) {
    vec Xbeta = X * beta;
    vec resid = y - Xbeta;
    mat GVi = G * Vi;
    vec eblup = Xbeta + GVi * resid;

    mat Ga = G - GVi * G;
    mat Gb = GVi * X;
    mat R = X - Gb;
    vec g2d = sum((R * Q) % R, 1);
    vec g1d = Ga.diag();

    return List::create(
      _["convergence"] = is_converged,
      _["sigma2_u"] = sigma2,
      _["rho"]      = rho_fix,
      _["beta"]     = beta,
      _["Xbeta"]    = Xbeta,
      _["theta"]    = eblup,
      _["GVi"]      = GVi,
      _["g1d"]      = g1d,
      _["g2d"]      = g2d,
      _["Q"]        = Q,
      _["XtVi"]     = XtVi,
      _["Vi"]       = Vi,
      _["G"]        = G
    );
  }

  // ================================================================
  // Full output mode
  // ================================================================
  vec stderr_beta = sqrt(Q.diag());
  vec zvalue = beta / stderr_beta;
  vec abs_z = arma::abs(zvalue);
  vec pvalue(zvalue.n_elem);
  for (uword i = 0; i < pvalue.n_elem; ++i) {
    pvalue(i) = 2.0 * R::pnorm5(-abs_z(i), 0.0, 1.0, 1, 0);
  }

  vec Xbeta = X * beta;
  vec resid = y - Xbeta;
  mat GVi = G * Vi;
  vec eblup = Xbeta + GVi * resid;

  // loglikelihood
  double sign, logdetV;
  log_det(logdetV, sign, V);
  double quadform = dot(resid, Vi * resid);
  const double loglike = -0.5 * (m * std::log(2.0 * M_PI) + logdetV + quadform);
  const double AIC = -2.0 * loglike + 2.0 * (p + 2);
  const double BIC = -2.0 * loglike + (p + 2) * std::log((double) m);

  Rcpp::NumericVector goodness = Rcpp::NumericVector::create(
    Rcpp::Named("loglikelihood") = loglike,
    Rcpp::Named("AIC") = AIC,
    Rcpp::Named("BIC") = BIC
  );

  // ================================================================
  // MSE computations
  // ================================================================
  mat Ga = G - GVi * G;
  mat Gb = GVi * X;
  mat R = X - Gb;
  vec g1d = Ga.diag();
  vec g2d = sum((R * Q) % R, 1);

  vec g3d(m, fill::zeros);
  vec g4d(m, fill::zeros);

  // turunan rho
  // der3 = (Sigma C Sigma)[s,s] = -dSigma/drho  (from sar_marginal)
  mat Amat = -sigma2 * der3;

  mat QXtVi = Q * XtVi;
  P = Vi - XtVi.t() * QXtVi;
  mat PD = P * derSigma;
  mat PR = P * Amat;

  // informasi deviance
  Idev.zeros();
  Idev(0, 0) = 0.5 * trace(PD * PD);
  Idev(0, 1) = 0.5 * trace(PD * PR);
  Idev(1, 0) = Idev(0, 1);
  Idev(1, 1) = 0.5 * trace(PR * PR);

  // For ML, the asymptotic bias b_ML = 1/2 I_ML^{-1} h and Var(th^) = I_ML^{-1}
  // use the ML information 1/2 tr(V^-1 V_a V^-1 V_b). (sae::mseSFH uses the
  // REML information also for ML; the two agree asymptotically, but when V_a
  // has a large component in col(X) -- e.g. sigma2 ~ 0 with rho ~ 1 and a
  // row-standardised W, where Sigma is dominated by the constant vector -- h
  // and grad g1 become huge and only I_ML cancels them. In simulations the ML
  // information also gave a smaller MSE bias with few domains.)
  if (isML) {
    const mat VD = Vi * derSigma, VR = Vi * Amat;
    Idev(0, 0) = 0.5 * accu(VD % VD.t());
    Idev(0, 1) = Idev(1, 0) = 0.5 * accu(VD % VR.t());
    Idev(1, 1) = 0.5 * accu(VR % VR.t());
  }

  mat Idevi;
  bool ok_idev = inv_sympd(Idevi, Idev);
  if (!ok_idev) Idevi = pinv(Idev);

  // rho estimated on the boundary of the parameter space (|rho| = 0.999, or
  // mapped to -1): the Prasad-Rao expansion in rho is not valid there (and can
  // explode, e.g. for rho -> -1 on a bipartite W), so rho is treated as fixed:
  // only the sigma2 block of the information is inverted. This drops the rho
  // terms from g3, g4, b_ML and the kriging MSE (and gives an NA SE for rho).
  const bool rho_at_bound = (rho_fix <= -0.999) || (rho_fix >= 0.999);
  if (rho_at_bound) {
    Idevi.zeros();
    if (Idev(0, 0) > 0.0) Idevi(0, 0) = 1.0 / Idev(0, 0);
  }

  // g3d
  mat ViD = Vi * derSigma;
  mat ViR = Vi * Amat;

  mat l1 = ViD - sigma2 * ViD * ViD;
  mat l2 = ViR - sigma2 * ViR * ViD;
  mat l1t = l1.t();
  mat l2t = l2.t();

  for (int i = 0; i < m; ++i) {
    mat L(2, m, fill::zeros);
    L.row(0) = l1t.row(i);
    L.row(1) = l2t.row(i);
    g3d(i) = trace(L * V * L.t() * Idevi);
  }

  // g4d
  mat psi = diagmat(vardir);
  mat D12aux = -der3;
  mat D22aux = 2 * sigma2 * SCSCS_f
  - 2 * sigma2 * SWWS_f;

  mat psiVi = psi * Vi;
  mat Vipsi = Vi * psi;
  mat D = (psiVi * D12aux * Vipsi) * (Idevi(0, 1) + Idevi(1, 0))
    + (psiVi * D22aux * Vipsi) * Idevi(1, 1);

  for (int i = 0; i < m; ++i) {
    g4d(i) = 0.5 * D(i, i);
  }

  vec mse2d = g1d + g2d + 2.0 * g3d - g4d;

  // ML bias correction
  vec bML(2, fill::zeros);
  if (isML) {
    mat ViX = Vi * X;
    double h1 = -trace(QXtVi * derSigma * ViX);
    double h2 = -trace(QXtVi * Amat * ViX);

    vec h(2);
    h(0) = h1;
    h(1) = h2;
    bML = Idevi * h / 2.0;

    mat GViCi = GVi * derSigma;
    mat GViAmat = GVi * Amat;

    mat dg1_dA = derSigma - 2.0 * GViCi + sigma2 * GViCi * ViD;
    mat dg1_dp = Amat - 2.0 * GViAmat + sigma2 * GViAmat * ViD;

    vec bMLgradg1(m, fill::zeros);
    for (int i = 0; i < m; ++i) {
      vec gradg1d(2);
      gradg1d(0) = dg1_dA(i, i);
      gradg1d(1) = dg1_dp(i, i);
      bMLgradg1(i) = as_scalar(bML.t() * gradg1d);
    }

    mse2d -= bMLgradg1;
  }

  // ================================================================
  // Synthetic estimators untuk area tidak tersampel
  // ================================================================
  if (adaNA) {
    // same full-domain SAR covariance as used in the fit
    mat G_full = sigma2 * Sigma_full;
    mat G_rs = G_full.submat(idx_ns, idx_s);
    mat G_rr = G_full.submat(idx_ns, idx_ns);

    mat KrigW = G_rs * Vi;
    u_ns = KrigW * resid;
    vec eblup_ns = Xns * beta + u_ns;

    mat Gb_ns = KrigW * X;
    mat R_ns = Xns - Gb_ns;

    mat Ga_ns = G_rr - KrigW * G_rs.t();
    vec g1d_ns = Ga_ns.diag();
    vec g2d_ns = sum((R_ns * Q) % R_ns, 1);

    // ---- second-order terms for the kriged (unsampled) domains ----
    // Parameters th = (sigma2, rho); G = sigma2 * Sigma_full(rho) on ALL
    // domains, with Sigma' = -Sigma C Sigma, Sigma'' = 2 Sigma C Sigma C Sigma
    // - 2 Sigma W'W Sigma (C = 2 rho W'W - W - W'). For a block B of G:
    //   B_0 = Sigma[B], B_1 = sigma2 Sigma'[B], B_01 = Sigma'[B],
    //   B_11 = sigma2 Sigma''[B], B_00 = 0.
    const mat C_full = 2.0 * rho_fix * WtW - (W + W.t());
    const mat Sd1 = -(Sigma_full * C_full * Sigma_full);                      // Sigma'
    const mat Sd2 = 2.0 * (-Sd1) * C_full * Sigma_full - 2.0 * Sigma_full * WtW * Sigma_full;  // Sigma''
    auto blk = [&](const mat& Mf, const uvec& r, const uvec& c) { return mat(Mf.submat(r, c)); };
    // A = G_rs and its derivatives; V_j = G_j[s,s]
    const mat A = G_rs;
    std::vector<mat> Aj = { blk(Sigma_full, idx_ns, idx_s), sigma2 * blk(Sd1, idx_ns, idx_s) };
    std::vector<mat> Vj = { derSigma, Amat };        // = Sigma[s,s], sigma2 Sigma'[s,s]
    // second derivatives (j,k): (0,0) = 0, (0,1) = Sigma', (1,1) = sigma2 Sigma''
    auto Ajk = [&](int j, int k) -> mat {
      if (j == 0 && k == 0) return mat(A.n_rows, A.n_cols, fill::zeros);
      if (j + k == 1) return blk(Sd1, idx_ns, idx_s);
      return sigma2 * blk(Sd2, idx_ns, idx_s);
    };
    auto Vjk = [&](int j, int k) -> mat {
      if (j == 0 && k == 0) return mat(m, m, fill::zeros);
      if (j + k == 1) return blk(Sd1, idx_s, idx_s);
      return sigma2 * blk(Sd2, idx_s, idx_s);
    };
    auto Grr_jk = [&](int j, int k) -> vec {
      if (j == 0 && k == 0) return vec(idx_ns.n_elem, fill::zeros);
      if (j + k == 1) return vec(blk(Sd1, idx_ns, idx_ns).diag());
      return sigma2 * vec(blk(Sd2, idx_ns, idx_ns).diag());
    };
    const mat& K = Vi;
    std::vector<mat> Kj = { -K * Vj[0] * K, -K * Vj[1] * K };
    auto dg = [](const mat& M1, const mat& M2) -> vec { return sum(M1 % M2, 1); };  // diag(M1 M2')

    // g3 (exact derivative of the kriging predictor):
    //   d theta~_r / d th_a = l_a' P y,
    //   l_a' = A_a - (A V^{-1} + R_r Q X'V^{-1}) V_a ,
    // so Cov = l_a' P l_b and g3_r = sum_ab Idevi_ab l_a' P l_b.
    const mat Bmat = KrigW + R_ns * Q * XtVi;   // (A V^{-1} + R Q X' V^{-1})
    std::vector<mat> L = { Aj[0] - Bmat * Vj[0], Aj[1] - Bmat * Vj[1] };
    vec g3d_ns(idx_ns.n_elem, fill::zeros);
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b)
        if (Idevi(a, b) != 0.0) g3d_ns += Idevi(a, b) * dg(L[a] * P, L[b]);

    // Bias of the plug-in g1_r = G_rr - A K A' (K = V^{-1}):
    //   E g1(th^) - g1 ~ b' grad g1 + 1/2 sum_jk Idevi_jk d2 g1 / d th_j d th_k
    //   d_j g1  = G_rr,j - 2 diag(A_j K A') - diag(A K_j A')
    //   d_jk g1 = G_rr,jk - 2 diag(A_jk K A') - 2 diag(A_j K_k A') - 2 diag(A_k K_j A')
    //            - 2 diag(A_j K A_k') - diag(A K_jk A'),
    //   K_j = -K V_j K, K_jk = K V_j K V_k K + K V_k K V_j K - K V_jk K.
    // For sampled domains this reproduces g1 + g2 + 2 g3 - g4 up to o(1/m).
    vec half_tr(idx_ns.n_elem, fill::zeros);
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        if (Idevi(j, k) == 0.0) continue;
        const mat Kjk = K * Vj[j] * K * Vj[k] * K + K * Vj[k] * K * Vj[j] * K - K * Vjk(j, k) * K;
        const vec d2 = Grr_jk(j, k) - 2.0 * dg(Ajk(j, k) * K, A) - 2.0 * dg(Aj[j] * Kj[k], A)
          - 2.0 * dg(Aj[k] * Kj[j], A) - 2.0 * dg(Aj[j] * K, Aj[k]) - dg(A * Kjk, A);
        half_tr += 0.5 * Idevi(j, k) * d2;
      }
    }
    vec mse_ns = g1d_ns + g2d_ns + g3d_ns - half_tr;
    if (isML) {
      const vec Grr_0 = vec(blk(Sigma_full, idx_ns, idx_ns).diag());
      const vec Grr_1 = sigma2 * vec(blk(Sd1, idx_ns, idx_ns).diag());
      const vec dg1_0 = Grr_0 - 2.0 * dg(Aj[0] * K, A) - dg(A * Kj[0], A);
      const vec dg1_1 = Grr_1 - 2.0 * dg(Aj[1] * K, A) - dg(A * Kj[1], A);
      mse_ns -= bML(0) * dg1_0 + bML(1) * dg1_1;
    }

    eblup_all.elem(idx_s)  = eblup;
    eblup_all.elem(idx_ns)  = eblup_ns;
    mse_all.elem(idx_s)     = mse2d;
    mse_all.elem(idx_ns)    = mse_ns;
  } else {
    eblup_all = eblup;
    mse_all = mse2d;
  }

  vec u_all(Xall.n_rows);
  if (!adaNA) {
    u_all = GVi * resid;
  } else {
    u_all.elem(idx_s) = GVi * resid;
    u_all.elem(idx_ns) = u_ns;
  }

  // ================================================================
  // Relative Standard Error (%)
  // ================================================================
  vec rse = sqrt(mse_all);
  vec relerr = rse / arma::abs(eblup_all);
  relerr.elem(arma::find(eblup_all == 0)).fill(datum::nan);
  rse = 100.0 * relerr;

  // ================================================================
  // Output
  // ================================================================
  DataFrame df_coef = DataFrame::create(
    _["beta"]        = beta,
    _["std.error"]   = stderr_beta,
    _["stderr_beta"] = stderr_beta,
    _["zvalue"]      = zvalue,
    _["pvalue"]      = pvalue
  );

  DataFrame df_eblup = DataFrame::create(
    _["y"]             = yall,
    _["vardir"]        = vardirall,
    _["eblup"]         = eblup_all,
    _["random_effect"] = u_all,
    _["mse"]           = mse_all,
    _["rse"]           = rse
  );

  // Standardized output structure
  double se_sigma2 = (Idevi(0, 0) > 0.0) ? std::sqrt(Idevi(0, 0)) : NA_REAL;
  double se_rho = (Idevi(1, 1) > 0.0) ? std::sqrt(Idevi(1, 1)) : NA_REAL;
  DataFrame estvarcomp = DataFrame::create(
    _["parameter"] = CharacterVector::create("sigma2_u", "rho"),
    _["estimate"]  = NumericVector::create(sigma2, rho_fix),
    _["std.error"] = NumericVector::create(se_sigma2, se_rho)
  );

  List out = List::create(
    _["estcoef"]          = df_coef,
    _["random_effect_var"] = sigma2,
    _["rho"]              = rho_fix,
    _["estvarcomp"]       = estvarcomp,
    _["goodness"]         = goodness,
    _["df_eblup"]         = df_eblup,
    _["model"]            = "SFH",
    _["level"]            = "area",
    _["n_iter"]           = k,
    _["convergence"]      = is_converged,
    _["method"]           = method
  );

  return out;
}
