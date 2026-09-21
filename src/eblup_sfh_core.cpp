// src/eblup_sfh_core.cpp
// Core Spatial Fay-Herriot EBLUP estimation with unified ML/REML loops
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "eblup_sfh.h"
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

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
    const std::string& method,
    int maxiter,
    double precision
) {
  const int m = X.n_rows;
  const int p = X.n_cols;

  // precompute konstanta
  mat Wt = W.t();
  mat eye_p = eye<mat>(p, p);
  mat I = eye<mat>(m, m);
  mat WtW = Wt * W;
  mat WpWt = W + Wt;

  // initial values
  vec sigma2_u(maxiter + 1, fill::zeros);
  vec rho(maxiter + 1, fill::zeros);
  sigma2_u(0) = median(vardir);
  rho(0) = 0.5;

  // pre-allocate temporaries
  vec s(2, fill::zeros);
  mat Idev(2, 2, fill::zeros);
  vec stime_fin(2, fill::zeros);
  mat A(m, m), Vi(m, m), derSigma(m, m);
  mat derRho(m, m), derVRho(m, m), V(m, m);
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
    A = (I - rho(k - 1) * Wt) * (I - rho(k - 1) * W);
    bool ok_inv = inv_sympd(derSigma, A);
    if (!ok_inv) derSigma = pinv(A);

    derRho = 2 * rho(k - 1) * WtW - WpWt;
    derVRho = -sigma2_u(k - 1) * (derSigma * derRho * derSigma);

    V = sigma2_u(k - 1) * derSigma + diagmat(vardir);
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
  double rho_fix = rho(k);
  if (rho_fix <= -0.999) rho_fix = -0.999;
  if (rho_fix >= 0.999)  rho_fix = 0.999;

  double sigma2 = sigma2_u(k);
  if (sigma2 < 0.0) sigma2 = 0.0;

  // spatial covariance
  A = (I - rho_fix * Wt) * (I - rho_fix * W);
  bool ok_inv = inv_sympd(derSigma, A);
  if (!ok_inv) derSigma = pinv(A);

  mat G = sigma2 * derSigma;
  V = G + diagmat(vardir);
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
  bool ok_all = ok_inv && beta.is_finite() && eblup.is_finite() && g1d.is_finite() && g2d.is_finite() &&
                arma::all(g1d >= 0.0) && (arma::max(g1d) < 1e6) && (arma::max(arma::abs(eblup)) < 1e6);

  SeblupFitArma out;
  out.theta     = eblup;
  out.g1d       = g1d;
  out.g2d       = g2d;
  out.sigma2    = sigma2;
  out.rho_fix   = rho_fix;
  out.n_iter    = k;
  out.converged = (k < maxiter) && (diff <= precision) && ok_all;
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

  // W untuk area tersampel saja
  mat W = adaNA ? Wall.submat(idx_s, idx_s) : Wall;

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
  mat Wt = W.t();
  mat eye_p = eye<mat>(p, p);
  mat I = eye<mat>(m, m);
  mat WtW = Wt * W;
  mat WpWt = W + Wt;

  vec sigma2_u(maxiter + 1, fill::zeros);
  vec rho(maxiter + 1, fill::zeros);
  sigma2_u(0) = median(vardir);
  rho(0) = 0.5;

  vec s(2, fill::zeros);
  mat Idev(2, 2, fill::zeros);
  vec stime_fin(2, fill::zeros);
  mat A(m, m), Vi(m, m), derSigma(m, m);
  mat derRho(m, m), derVRho(m, m), V(m, m);
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
    A = (I - rho(k - 1) * Wt) * (I - rho(k - 1) * W);
    bool ok_inv = inv_sympd(derSigma, A);
    if (!ok_inv) derSigma = pinv(A);

    derRho = 2 * rho(k - 1) * WtW - WpWt;
    derVRho = -sigma2_u(k - 1) * (derSigma * derRho * derSigma);

    V = sigma2_u(k - 1) * derSigma + diagmat(vardir);
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
  double rho_fix = rho(k);
  if (rho_fix == -0.999) rho_fix = -1.0;
  else if (rho_fix == 0.999) rho_fix = 1.0;

  double sigma2 = sigma2_u(k);
  if (sigma2 < 0.0) sigma2 = 0.0;

  // spatial covariance
  A = (I - rho_fix * Wt) * (I - rho_fix * W);
  bool ok_inv = inv_sympd(derSigma, A);
  if (!ok_inv) derSigma = pinv(A);

  mat G = sigma2 * derSigma;
  V = G + diagmat(vardir);
  ok_inv = inv_sympd(Vi, V) && ok_inv;
  if (!ok_inv) Vi = pinv(V);

  // beta
  XtVi = X.t() * Vi;
  XtViX = XtVi * X;
  mat XtViy = XtVi * y;

  ok_inv = inv_sympd(Q, XtViX) && ok_inv;
  if (!ok_inv) Q = solve(XtViX, eye_p);

  mat beta = Q * XtViy;

  bool is_converged = (k < maxiter) && (diff <= precision);

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
  mat derRho_final = 2 * rho_fix * WtW - WpWt;
  mat rhosigma = derRho_final * derSigma;
  mat der3 = derSigma * rhosigma;
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

  mat Idevi;
  bool ok_idev = inv_sympd(Idevi, Idev);
  if (!ok_idev) Idevi = pinv(Idev);

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
  mat D22aux = 2 * sigma2 * der3 * rhosigma
  - 2 * sigma2 * derSigma * WtW * derSigma;

  mat psiVi = psi * Vi;
  mat Vipsi = Vi * psi;
  mat D = (psiVi * D12aux * Vipsi) * (Idevi(0, 1) + Idevi(1, 0))
    + (psiVi * D22aux * Vipsi) * Idevi(1, 1);

  for (int i = 0; i < m; ++i) {
    g4d(i) = 0.5 * D(i, i);
  }

  vec mse2d = g1d + g2d + 2.0 * g3d - g4d;

  // ML bias correction
  if (isML) {
    mat ViX = Vi * X;
    double h1 = -trace(QXtVi * derSigma * ViX);
    double h2 = -trace(QXtVi * Amat * ViX);

    vec h(2);
    h(0) = h1;
    h(1) = h2;
    vec bML = Idevi * h / 2.0;

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
    const int m_total = Xall.n_rows;
    mat I_full = eye<mat>(m_total, m_total);
    mat Wallt = Wall.t();
    mat A_full_target = (I_full - rho_fix * Wallt) * (I_full - rho_fix * Wall);
    mat A_full;
    bool ok_full = inv_sympd(A_full, A_full_target);
    if (!ok_full) A_full = pinv(A_full_target);

    mat G_full = sigma2 * A_full;
    mat G_rs = G_full.submat(idx_ns, idx_s);
    mat G_rr = G_full.submat(idx_ns, idx_ns);

    mat KrigW = G_rs * Vi;
    mat u_ns = KrigW * resid;
    vec eblup_ns = Xns * beta + u_ns;

    mat Gb_ns = KrigW * X;
    mat R_ns = Xns - Gb_ns;

    mat Ga_ns = G_rr - KrigW * G_rs.t();
    vec g1d_ns = Ga_ns.diag();
    vec g2d_ns = sum((R_ns * Q) % R_ns, 1);
    vec mse_ns = g1d_ns + g2d_ns;

    vec u_sampled = GVi * resid;
    vec u_all(Xall.n_rows);
    u_all.elem(idx_s) = u_sampled;
    u_all.elem(idx_ns) = u_ns;

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
    // already populated inside if (adaNA)
    u_all.elem(idx_s) = GVi * resid;
    // idx_ns already set
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
