// src/eblup_sfh_core.cpp
// Core Spatial Fay-Herriot EBLUP estimation with unified ML/REML loops
#include <RcppArmadillo.h>
#include "sar_engine.h"
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
                  arma::mat* SCSCS, arma::mat* SWWS, arma::mat* Sfull,
                  arma::mat* SCS_full, arma::mat* D2_full) {
  const arma::uword D = Wf.n_rows;
  const arma::mat I = arma::eye<arma::mat>(D, D);
  // A = (I - rho W')(I - rho W) = I - rho (W + W') + rho^2 W'W, using the
  // precomputed W'W: O(D^2) instead of a D^3 product per call.
  // (symmetric in exact arithmetic; symmetrise the round-off)
  const arma::mat A = arma::symmatu(I - rho * (Wf + Wf.t()) + (rho * rho) * WtWf);
  arma::mat Sf;
  bool ok = arma::inv_sympd(Sf, A);
  if (!ok) Sf = arma::pinv(A);
  const arma::mat C = 2 * rho * WtWf - (Wf + Wf.t());
  const arma::mat CS = C * Sf;
  const arma::mat SCSf = Sf * CS;
  Sig = all ? Sf : arma::mat(Sf.submat(s, s));
  SCS = all ? SCSf : arma::mat(SCSf.submat(s, s));
  arma::mat t1, t2;
  if (SCSCS || D2_full) t1 = SCSf * CS;            // Sigma C Sigma C Sigma
  if (SWWS || D2_full) t2 = Sf * WtWf * Sf;        // Sigma W'W Sigma
  if (SCSCS) *SCSCS = all ? t1 : arma::mat(t1.submat(s, s));
  if (SWWS) *SWWS = all ? t2 : arma::mat(t2.submat(s, s));
  if (Sfull) *Sfull = Sf;
  if (SCS_full) *SCS_full = SCSf;                  // = -Sigma'  (all domains)
  if (D2_full) *D2_full = 2.0 * t1 - 2.0 * t2;     // =  Sigma'' (all domains)
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
// SFH evaluation through SarEngine (see sar_engine.h)
//
// With A = (I - rho W')(I - rho W), Sigma = A^-1, E the m x D selection of the
// sampled domains, Psi = diag(vardir) and M = A + s2 E' Psi^-1 E:
//   V^-1         = Psi^-1 - s2 Psi^-1 E M^-1 E' Psi^-1          (Woodbury)
//   V^-1 E Sigma = Psi^-1 E M^-1                                 (push-through)
//   log|V|       = log|Psi| + log|M| - log|A|
//   G - G V^-1 G = s2 (M^-1)_ss,   G V^-1 = s2 (M^-1)_ss Psi^-1
// so the fit and the MSE only need solves with the sparse A and M (m or D
// right-hand sides) and no D x D inverse or D^3 product.
// ============================================================================
struct SfhState {
  mat Z;      // M^-1 E'      (D x m)
  mat Vi;     // V^-1         (m x m)
  mat XtVi;   // X' V^-1
  mat XtViX;  // X' V^-1 X
  mat Q;      // (X' V^-1 X)^-1
  vec Py;     // P y
};

static bool sfh_state(SarEngine& S, double s2, double rho, const mat& X, const vec& y,
                      const vec& psi, SfhState& st) {
  if (!S.factor(s2, rho, psi)) return false;
  const uword m = psi.n_elem;
  st.Z = S.solveM(S.scatter(eye<mat>(m, m)));
  const vec ipsi = 1.0 / psi;
  mat Vi = st.Z.rows(S.s());
  Vi.each_col() %= ipsi;
  Vi.each_row() %= ipsi.t();
  Vi *= -s2;
  Vi.diag() += ipsi;
  st.Vi = symmatu(Vi);
  st.XtVi = X.t() * st.Vi;
  st.XtViX = symmatu(st.XtVi * X);
  if (!inv_sympd(st.Q, st.XtViX)) {
    if (!solve(st.Q, st.XtViX, eye<mat>(X.n_cols, X.n_cols))) return false;
  }
  st.Py = st.Vi * y - st.XtVi.t() * (st.Q * (st.XtVi * y));
  return st.Z.is_finite() && st.Vi.is_finite() && st.Py.is_finite();
}

// V^-1 E Sigma = Psi^-1 Z'  (m x D)
static mat sfh_left_V(const SfhState& st, const vec& psi) {
  mat L = st.Z.t();
  L.each_col() /= psi;
  return L;
}

// P E Sigma = V^-1 E Sigma - V^-1 X Q X' V^-1 E Sigma
static mat sfh_left_P(const SfhState& st, const mat& X, const mat& Lv) {
  return Lv - st.XtVi.t() * (st.Q * (X.t() * Lv));
}

// For L = K E Sigma (K = V^-1 or P): KD = K Sigma_ss = L E' and
// KR = K V_rho = -s2 K E Sigma C Sigma E' = -s2 (rows_s(Sigma C L'))'
static void sfh_dmats(const SarEngine& S, double s2, const mat& L, mat& KD, mat& KR) {
  KD = L.cols(S.s());
  KR = -s2 * trans(mat(S.solveA(S.mulC(L.t())).rows(S.s())));
}

// Score quadratic forms: Py' Sigma_ss Py and Py' V_rho Py
static void sfh_quad(const SarEngine& S, double s2, const vec& Py, double& qs, double& qr) {
  const vec a = S.solveA(S.scatter(mat(Py)));  // Sigma E' P y
  qs = dot(Py, a.elem(S.s()));
  qr = -s2 * dot(a, vec(S.mulC(a)));
}

// Score and Fisher information of (s2, rho): REML uses P, ML uses V^-1 in the
// traces (the score quadratic form always uses P y)
static bool sfh_score_info(const SarEngine& S, const SfhState& st, double s2, const mat& X,
                           const vec& psi, bool isML, vec& Sc, mat& F) {
  const mat Lv = sfh_left_V(st, psi);
  mat KD, KR;
  sfh_dmats(S, s2, isML ? Lv : sfh_left_P(st, X, Lv), KD, KR);
  double qs, qr;
  sfh_quad(S, s2, st.Py, qs, qr);
  Sc.set_size(2);
  F.set_size(2, 2);
  Sc(0) = -0.5 * trace(KD) + 0.5 * qs;
  Sc(1) = -0.5 * trace(KR) + 0.5 * qr;
  F(0, 0) = 0.5 * accu(KD % KD.t());
  F(0, 1) = F(1, 0) = 0.5 * accu(KD % KR.t());
  F(1, 1) = 0.5 * accu(KR % KR.t());
  return Sc.is_finite() && F.is_finite();
}

// ============================================================================
// sfh_eval()
//
// (Restricted) log-likelihood of the SFH model at th = (sigma2, rho), up to a
// constant, with optionally the score and Fisher information:
//   REML: -1/2 [ log|V| + log|X'V^{-1}X| + y'Py ]
//   ML  : -1/2 [ log|V| + (y - X b)'V^{-1}(y - X b) ]   (b = GLS estimate)
// ============================================================================
static bool sfh_eval(SarEngine& S, const vec& th, const mat& X, const vec& y, const vec& psi,
                     bool isML, bool want_SF, double& ll, vec& Sc, mat& F) {
  SfhState st;
  if (!sfh_state(S, th(0), th(1), X, y, psi, st) || !S.M_pd()) return false;
  const double logdetV = accu(log(psi)) + S.logdetM() - S.logdetA();
  const double yPy = dot(y, st.Py);
  if (isML) {
    ll = -0.5 * (logdetV + yPy);
  } else {
    double ld, sg;
    if (!log_det(ld, sg, st.XtViX) || sg <= 0) return false;
    ll = -0.5 * (logdetV + ld + yPy);
  }
  if (!std::isfinite(ll)) return false;
  if (!want_SF) return true;
  return sfh_score_info(S, st, th(0), X, psi, isML, Sc, F);
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
static bool sfh_fit_safeguarded(SarEngine& S, const mat& X, const vec& y, const vec& psi,
                                bool isML, int maxiter, double precision, vec& th, int& iters) {
  auto project = [](vec t) {
    t(0) = std::max(t(0), 0.0);
    t(1) = std::min(std::max(t(1), -0.999), 0.999);
    return t;
  };
  th = project(th);
  double ll;
  vec Sc;
  mat F;
  if (!sfh_eval(S, th, X, y, psi, isML, true, ll, Sc, F)) return false;

  const int max_it = 5 * maxiter;
  for (iters = 1; iters <= max_it; ++iters) {
    // active set: sigma2 at 0 with outward score; rho unidentified when
    // sigma2 = 0; rho at a bound with outward score
    std::vector<uword> fr;
    if (!(th(0) == 0.0 && Sc(0) <= 0.0)) fr.push_back(0);
    if (!(th(0) == 0.0) && !((th(1) >= 0.999 && Sc(1) >= 0.0) || (th(1) <= -0.999 && Sc(1) <= 0.0)))
      fr.push_back(1);
    if (fr.empty()) return true;
    const uvec f = conv_to<uvec>::from(fr);

    // converged when the free scores are negligible in SE units
    double zmax = 0.0;
    for (uword a : fr) if (F(a, a) > 0.0) zmax = std::max(zmax, std::fabs(Sc(a)) / std::sqrt(F(a, a)));
    if (zmax < 10.0 * precision) return true;

    std::vector<vec> dirs;
    {
      vec sf, step(2, fill::zeros);
      if (solve(sf, mat(F.submat(f, f)), vec(Sc.elem(f)), solve_opts::no_approx)) {
        step.elem(f) = sf;
        dirs.push_back(step);
      }
      vec g(2, fill::zeros);
      for (uword a : fr) g(a) = Sc(a) / std::max(F(a, a), 1e-12);
      dirs.push_back(g);
    }

    double llc = 0.0;
    vec cand, Scc;
    mat Fc;
    bool accepted = false;
    for (const vec& step : dirs) {
      double t = 1.0;
      for (int h = 0; h < 30; ++h) {
        cand = project(th + t * step);
        if (sfh_eval(S, cand, X, y, psi, isML, true, llc, Scc, Fc) && llc > ll) {
          accepted = true;
          break;
        }
        t *= 0.5;
      }
      if (accepted) break;
    }
    if (!accepted) return true;  // no ascent direction left: constrained maximum
    th = cand; ll = llc; Sc = Scc; F = Fc;
  }
  return false;
}

// ============================================================================
// sfh_fit()
//
// Classical Fisher scoring as sae::eblupSFH (start sigma2 = median(vardir),
// rho = 0.5; rho clipped to [-0.999, 0.999]; relative-change stopping rule),
// followed by the constrained fallback when it fails or ends with sigma2 < 0.
// Interior solutions of the classical iteration are unchanged.
// ============================================================================
struct SfhFitRes {
  double s2, rho;
  int k;
  bool conv;
};

static SfhFitRes sfh_fit(SarEngine& S, const mat& X, const vec& y, const vec& psi, bool isML,
                         int maxiter, double precision) {
  vec sig(maxiter + 1, fill::zeros), rho(maxiter + 1, fill::zeros);
  sig(0) = median(psi);
  rho(0) = 0.5;
  double diff = precision + 1.0;
  int k = 0;
  SfhState st;
  vec Sc;
  mat F;
  while ((diff > precision) && (k < maxiter)) {
    ++k;
    if (!sfh_state(S, sig(k - 1), rho(k - 1), X, y, psi, st) ||
        !sfh_score_info(S, st, sig(k - 1), X, psi, isML, Sc, F)) {
      k = maxiter;
      break;
    }
    vec par(2);
    par(0) = sig(k - 1);
    par(1) = rho(k - 1);
    vec step;
    if (!solve(step, F, Sc) || !step.is_finite()) {
      k = maxiter;
      break;
    }
    vec nw = par + step;
    if (nw(1) <= -0.999) nw(1) = -0.999;
    if (nw(1) >= 0.999) nw(1) = 0.999;
    sig(k) = nw(0);
    rho(k) = nw(1);
    diff = max(abs(nw - par) / abs(par + 1e-12));
    if (!std::isfinite(diff)) {
      k = maxiter;
      break;
    }
  }

  SfhFitRes out;
  out.conv = (k < maxiter) && (diff <= precision);
  out.s2 = sig(k);
  out.rho = rho(k);
  out.k = k;
  if (!out.conv || out.s2 < 0.0) {
    vec th(2);
    th(0) = out.conv ? out.s2 : median(psi);
    th(1) = out.conv ? out.rho : 0.5;
    int it2 = 0;
    out.conv = sfh_fit_safeguarded(S, X, y, psi, isML, maxiter, precision, th, it2);
    out.k += it2;
    out.s2 = th(0);
    out.rho = th(1);
  }
  if (out.s2 < 0.0) out.s2 = 0.0;
  return out;
}

// ============================================================================
// seblup_fit_core_arma()
//
// Murni C++: TIDAK membuat objek R dan TIDAK memanggil R RNG/Math.
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
  SeblupFitArma out;
  out.converged = false;
  SarEngine S(W, idx_s);
  const SfhFitRes fr = sfh_fit(S, X, y, vardir, method == "ML", maxiter, precision);

  // Same boundary handling as seblup_core().
  double rho_fix = sfh_lower_rho(fr.rho, W);
  if (rho_fix >= 0.999) rho_fix = 0.999;
  const double sigma2 = fr.s2;

  SfhState st;
  out.sigma2 = sigma2;
  out.rho_fix = rho_fix;
  out.n_iter = fr.k;
  if (!sfh_state(S, sigma2, rho_fix, X, y, vardir, st)) return out;

  const mat Zs = st.Z.rows(idx_s);
  const mat beta = st.Q * (st.XtVi * y);
  const vec Xbeta = X * beta;
  mat GVi = sigma2 * Zs;                 // G V^-1 = s2 (M^-1)_ss Psi^-1
  GVi.each_row() /= vardir.t();
  const vec eblup = Xbeta + GVi * (y - Xbeta);
  const vec g1d = sigma2 * vec(Zs.diag());
  const mat R = X - GVi * X;
  const vec g2d = sum((R * st.Q) % R, 1);

  // Validity must not depend on the scale of y (no absolute thresholds).
  const bool ok_all = beta.is_finite() && eblup.is_finite() && g1d.is_finite() &&
                      g2d.is_finite() && arma::all(g1d >= 0.0);
  out.theta = eblup;
  out.g1d = g1d;
  out.g2d = g2d;
  out.converged = fr.conv && ok_all;
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

  arma::vec eblup_all(Xall.n_rows);
  arma::vec mse_all(Xall.n_rows);

  // ================================================================
  // Handle NA values (unsampled areas)
  // ================================================================
  arma::mat Xns, X;
  arma::vec y, vardir;
  arma::uvec idx_ns, idx_s;
  arma::vec u_ns;  // kriged random effects for unsampled areas (filled when adaNA)
  const bool adaNA = yall.has_nan();

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
    idx_s  = arma::regspace<arma::uvec>(0, Xall.n_rows - 1);
  }

  // The SAR model is defined on ALL domains of Wall; the sampled domains
  // have the marginal covariance s2 * Sigma_full[s,s].
  const mat& W = Wall;
  const int m = X.n_rows;
  const int p = X.n_cols;

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
  const bool isML = (method == "ML");

  // ================================================================
  // Estimation (classical Fisher scoring + constrained fallback)
  // ================================================================
  SarEngine S(W, idx_s);
  const SfhFitRes fr = sfh_fit(S, X, y, vardir, isML, maxiter, precision);
  const bool is_converged = fr.conv;
  const int k = fr.k;

  // Lower bound: -0.999 -> -1 as sae::eblupSFH, when rho = -1 is valid
  // (see sfh_lower_rho). Upper bound: never mapped to +1 (I - W singular).
  const double rho_fix = sfh_lower_rho(fr.rho, W);
  const double sigma2 = fr.s2;

  // ================================================================
  // Final estimates
  // ================================================================
  SfhState st;
  if (!sfh_state(S, sigma2, rho_fix, X, y, vardir, st)) {
    Rcpp::stop("The SFH covariance matrix is singular at the final estimates (sigma2 = %g, rho = %g).",
               sigma2, rho_fix);
  }
  const vec ipsi = 1.0 / vardir;
  const mat Zs = symmatu(mat(st.Z.rows(idx_s)));          // (M^-1)_ss
  const mat SE = S.solveA(S.scatter(eye<mat>(m, m)));      // Sigma[:, s]
  const mat derSigma = symmatu(mat(SE.rows(idx_s)));       // Sigma_ss
  const mat G = sigma2 * derSigma;
  const mat& Vi = st.Vi;
  const mat& XtVi = st.XtVi;
  const mat& Q = st.Q;

  const mat beta = Q * (XtVi * y);
  const vec Xbeta = X * beta;
  const vec resid = y - Xbeta;
  mat GVi = sigma2 * Zs;                                    // G V^-1 = s2 (M^-1)_ss Psi^-1
  GVi.each_row() %= ipsi.t();
  const vec eblup = Xbeta + GVi * resid;

  const vec g1d = sigma2 * vec(Zs.diag());                  // diag(G - G V^-1 G)
  const mat R = X - GVi * X;
  const vec g2d = sum((R * Q) % R, 1);

  if (only_core) {
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

  // loglikelihood: log|V| = log|Psi| + log|M| - log|A|
  const double logdetV = accu(log(vardir)) + S.logdetM() - S.logdetA();
  const double quadform = dot(resid, Vi * resid);
  // Full (ML-form) log-likelihood at the estimates for ML and REML, as in
  // sae::eblupSFH, so that AIC/BIC compare models with different covariates.
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
  // Fisher information at the estimates: REML 1/2 tr(P V_a P V_b). For ML,
  // the asymptotic bias b_ML = 1/2 I_ML^{-1} h and Var(th^) = I_ML^{-1} use
  // the ML information 1/2 tr(V^-1 V_a V^-1 V_b). (sae::mseSFH uses the
  // REML information also for ML; the two agree asymptotically, but when V_a
  // has a large component in col(X) -- e.g. sigma2 ~ 0 with rho ~ 1 and a
  // row-standardised W, where Sigma is dominated by the constant vector -- h
  // and grad g1 become huge and only I_ML cancels them. In simulations the ML
  // information also gave a smaller MSE bias with few domains.)
  const mat Lv = sfh_left_V(st, vardir);
  mat VD, VR;                                    // V^-1 Sigma_ss, V^-1 V_rho
  sfh_dmats(S, sigma2, Lv, VD, VR);
  mat Idev(2, 2);
  {
    mat KD = VD, KR = VR;
    if (!isML) sfh_dmats(S, sigma2, sfh_left_P(st, X, Lv), KD, KR);
    Idev(0, 0) = 0.5 * accu(KD % KD.t());
    Idev(0, 1) = Idev(1, 0) = 0.5 * accu(KD % KR.t());
    Idev(1, 1) = 0.5 * accu(KR % KR.t());
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

  // g3_i = psi_i^2 sum_ab I^ab (V^-1 V_a V^-1 V_b V^-1)_ii. With
  // V_a = E Sigma N_a Sigma E' (N_s2 = A, N_rho = -s2 C) and V^-1 E Sigma =
  // Psi^-1 E M^-1:
  //   g3_i = sum_ab I^ab [ s2 (Dp Y_a)_i' Sigma (Dp Y_b)_i + (Y_a)_i' Dp (Y_b)_i ],
  // Y_a = M^-1 N_a Z, Dp = E' Psi^-1 E (columns i = sampled domains).
  const uword D = W.n_rows;
  vec dpsi(D, fill::zeros);
  dpsi.elem(idx_s) = ipsi;
  vec g3d(m, fill::zeros);
  {
    std::vector<mat> Y(2), DY(2), SDY(2);
    const bool need_r = (Idevi(0, 1) != 0.0 || Idevi(1, 0) != 0.0 || Idevi(1, 1) != 0.0);
    Y[0] = S.solveM(S.mulA(st.Z));
    if (need_r) Y[1] = -sigma2 * S.solveM(S.mulC(st.Z));
    for (int a = 0; a < 2; ++a) {
      if (a == 1 && !need_r) break;
      DY[a] = Y[a];
      DY[a].each_col() %= dpsi;
      SDY[a] = S.solveA(DY[a]);
    }
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b)
        if (Idevi(a, b) != 0.0)
          g3d += Idevi(a, b) * (sigma2 * sum(DY[a] % SDY[b], 0).t() + sum(Y[a] % DY[b], 0).t());
  }

  // g4 = 1/2 sum_jk I^jk psi_i^2 (V^-1 V_jk V^-1)_ii with V_{s2 rho} = Sigma'_ss
  // = -(Sigma C Sigma)_ss and V_{rho rho} = s2 Sigma''_ss = 2 s2 (Sigma C Sigma C
  // Sigma - Sigma W'W Sigma)_ss; (V^-1 E Sigma K Sigma E' V^-1)_ii =
  // psi_i^-2 (Z' K Z)_ii.
  const mat CZ = S.mulC(st.Z);
  const vec zCz = sum(st.Z % CZ, 0).t();          // diag(Z' C Z)
  vec g4d(m, fill::zeros);
  g4d += -0.5 * (Idevi(0, 1) + Idevi(1, 0)) * zCz;
  if (Idevi(1, 1) != 0.0) {
    const mat WZ = S.mulW(st.Z);
    g4d += Idevi(1, 1) * sigma2 * (sum(CZ % S.solveA(CZ), 0).t() - sum(WZ % WZ, 0).t());
  }

  vec mse2d = g1d + g2d + 2.0 * g3d - g4d;

  // ML bias correction: - b_ML' grad g1, b_ML = 1/2 I^-1 h,
  // h_a = -tr(Q X'V^-1 V_a V^-1 X); g1 = s2 (M^-1)_ss, so
  // d g1 / d s2 = (M^-1)_ii - s2 sum_j (M^-1)_ij^2 / psi_j and
  // d g1 / d rho = -s2 (Z' C Z)_ii.
  vec bML(2, fill::zeros);
  if (isML) {
    const mat ViX = XtVi.t();
    const mat SViX = S.solveA(S.scatter(ViX));                  // Sigma E' V^-1 X
    const mat VsViX = SViX.rows(idx_s);
    const mat VrViX = -sigma2 * mat(S.solveA(S.mulC(SViX)).rows(idx_s));
    vec h(2);
    h(0) = -trace(Q * ViX.t() * VsViX);
    h(1) = -trace(Q * ViX.t() * VrViX);
    bML = Idevi * h / 2.0;

    mat Zp = Zs;
    Zp.each_row() %= ipsi.t();
    const vec dg_s = vec(Zs.diag()) - sigma2 * sum(Zs % Zp, 1);
    const vec dg_r = -sigma2 * zCz;
    mse2d -= bML(0) * dg_s + bML(1) * dg_r;
  }

  // ================================================================
  // Synthetic (kriging) estimators for unsampled areas
  // ================================================================
  if (adaNA) {
    const uword r = idx_ns.n_elem;
    // same full-domain SAR covariance as used in the fit
    const mat SEn = S.solveA(S.scatter(idx_ns, eye<mat>(r, r)));   // Sigma[:, ns]
    const mat Sig_rs = SE.rows(idx_ns);
    const vec Sig_rr = vec(mat(SEn.rows(idx_ns)).diag());
    const mat G_rs = sigma2 * Sig_rs;

    const mat KrigW = G_rs * Vi;
    u_ns = KrigW * resid;
    vec eblup_ns = Xns * beta + u_ns;

    const mat R_ns = Xns - KrigW * X;
    const vec g1d_ns = sigma2 * Sig_rr - sum(KrigW % G_rs, 1);
    const vec g2d_ns = sum((R_ns * Q) % R_ns, 1);

    // ---- second-order terms for the kriged (unsampled) domains ----
    // Parameters th = (sigma2, rho); G = sigma2 * Sigma(rho) on ALL
    // domains, with Sigma' = -Sigma C Sigma, Sigma'' = 2 Sigma C Sigma C Sigma
    // - 2 Sigma W'W Sigma (C = 2 rho W'W - W - W'). For a block B of G:
    //   B_0 = Sigma[B], B_1 = sigma2 Sigma'[B], B_01 = Sigma'[B],
    //   B_11 = sigma2 Sigma''[B], B_00 = 0.
    // Blocks via solves with A: Sigma'[:, s] = -A^-1 C Sigma[:, s],
    // Sigma''[:, s] = A^-1 (2 C (Sigma C Sigma)[:, s] - 2 W'W Sigma[:, s]),
    // diagonals of the [ns, ns] blocks from Sigma[:, ns].
    const mat CS1 = S.solveA(S.mulC(SE));                     // (Sigma C Sigma)[:, s]
    const mat Sd1_rs = -CS1.rows(idx_ns);
    const mat Sd1_ss = -symmatu(mat(CS1.rows(idx_s)));
    const mat T2 = S.solveA(2.0 * S.mulC(CS1) - 2.0 * S.mulWtW(SE));
    const mat Sd2_rs = T2.rows(idx_ns);
    const mat Sd2_ss = symmatu(mat(T2.rows(idx_s)));
    const mat CSEn = S.mulC(SEn), WSEn = S.mulW(SEn);
    const vec Sd1_rr = -sum(SEn % CSEn, 0).t();
    const vec Sd2_rr = 2.0 * sum(CSEn % S.solveA(CSEn), 0).t() - 2.0 * sum(WSEn % WSEn, 0).t();
    const mat Amat = sigma2 * Sd1_ss;                         // V_rho

    // A = G_rs and its derivatives; V_j = G_j[s,s]
    const mat& A = G_rs;
    std::vector<mat> Aj = { Sig_rs, sigma2 * Sd1_rs };
    std::vector<mat> Vj = { derSigma, Amat };
    // second derivatives (j,k): (0,0) = 0, (0,1) = Sigma', (1,1) = sigma2 Sigma''
    auto Ajk = [&](int j, int kk) -> mat {
      if (j == 0 && kk == 0) return mat(A.n_rows, A.n_cols, fill::zeros);
      if (j + kk == 1) return Sd1_rs;
      return sigma2 * Sd2_rs;
    };
    auto Vjk = [&](int j, int kk) -> mat {
      if (j == 0 && kk == 0) return mat(m, m, fill::zeros);
      if (j + kk == 1) return Sd1_ss;
      return sigma2 * Sd2_ss;
    };
    auto Grr_jk = [&](int j, int kk) -> vec {
      if (j == 0 && kk == 0) return vec(r, fill::zeros);
      if (j + kk == 1) return Sd1_rr;
      return sigma2 * Sd2_rr;
    };
    const mat& K = Vi;
    auto dg = [](const mat& M1, const mat& M2) -> vec { return sum(M1 % M2, 1); };  // diag(M1 M2')
    // r x m products only (r = #unsampled): AK = A K, AjK = A_j K
    const mat AK = A * K;
    std::vector<mat> AjK = { Aj[0] * K, Aj[1] * K };
    std::vector<mat> AKVj = { AK * Vj[0], AK * Vj[1] };   // A K V_j

    // g3 (exact derivative of the kriging predictor):
    //   d theta~_r / d th_a = l_a' P y,
    //   l_a' = A_a - (A V^{-1} + R_r Q X'V^{-1}) V_a ,
    // so Cov = l_a' P l_b and g3_r = sum_ab Idevi_ab l_a' P l_b.
    const mat P = Vi - XtVi.t() * Q * XtVi;
    const mat Bmat = KrigW + R_ns * Q * XtVi;   // (A V^{-1} + R Q X' V^{-1})
    std::vector<mat> L = { Aj[0] - Bmat * Vj[0], Aj[1] - Bmat * Vj[1] };
    vec g3d_ns(r, fill::zeros);
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
    vec half_tr(r, fill::zeros);
    for (int j = 0; j < 2; ++j) {
      for (int kk = 0; kk < 2; ++kk) {
        if (Idevi(j, kk) == 0.0) continue;
        // diag(A K_jk A') = diag(A K V_j K V_k K A') + diag(A K V_k K V_j K A')
        //                  - diag(A K V_jk K A'),  with (K A')' = A K
        const vec AKjkA = dg(AKVj[j] * K * Vj[kk], AK) + dg(AKVj[kk] * K * Vj[j], AK)
                        - dg(AK * Vjk(j, kk), AK);
        // diag(A_j K_k A') = -diag(A_j K V_k K A') = -diag((A_j K V_k)(A K)')
        const vec AjKkA = -dg(AjK[j] * Vj[kk], AK);
        const vec AkKjA = -dg(AjK[kk] * Vj[j], AK);
        const vec d2 = Grr_jk(j, kk) - 2.0 * dg(Ajk(j, kk) * K, A) - 2.0 * AjKkA
          - 2.0 * AkKjA - 2.0 * dg(AjK[j], Aj[kk]) - AKjkA;
        half_tr += 0.5 * Idevi(j, kk) * d2;
      }
    }
    vec mse_ns = g1d_ns + g2d_ns + g3d_ns - half_tr;
    if (isML) {
      const vec Grr_1 = sigma2 * Sd1_rr;
      // diag(A K_j A') = -diag((A K V_j)(A K)')
      const vec dg1_0 = Sig_rr - 2.0 * dg(AjK[0], A) + dg(AKVj[0], AK);
      const vec dg1_1 = Grr_1 - 2.0 * dg(AjK[1], A) + dg(AKVj[1], AK);
      mse_ns -= bML(0) * dg1_0 + bML(1) * dg1_1;
    }

    eblup_all.elem(idx_s)  = eblup;
    eblup_all.elem(idx_ns) = eblup_ns;
    mse_all.elem(idx_s)    = mse2d;
    mse_all.elem(idx_ns)   = mse_ns;
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
