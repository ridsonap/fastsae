// src/eblup_bhf_core.cpp
// REML/ML fit and parametric bootstrap MSE for the Battese-Harter-Fuller
// nested error model  y_dj = x_dj' beta + u_d + e_dj.
//
// With gamma = s2u / s2e the marginal covariance of domain d is
// s2e * (I + gamma J); beta and s2e are profiled out analytically, so the fit
// is a 1-D maximisation in gamma >= 0. Every evaluation only needs the
// within-domain cross-products (fixed) and the domain means:
//   X'V0^-1 X = Wxx + sum_d w_d xbar_d xbar_d',  w_d = n_d / (1 + n_d gamma)
// so its cost is O(D p^2 + p^3), independent of the number of units.
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

namespace {

// Sufficient statistics of the sampled domains (S of them)
struct BhfStats {
  double n = 0;
  int p = 0;
  mat Wxx;   // within-domain cross-products of X (p x p)
  vec Wxy;   // within-domain cross-products of X and y (p)
  double Wyy = 0;
  mat xbar;  // S x p domain means of X
  vec ybar;  // S domain means of y
  vec nd;    // S sample sizes
};

struct BhfEval {
  bool ok = false;
  double ll = 0, grad = 0, Q = 0, logdetA = 0;
  vec beta, rbar, w;
  mat Rchol;  // upper Cholesky factor of X'V0^-1 X
};

// Profiled (restricted) log-likelihood in gamma and its derivative
BhfEval bhf_eval(const BhfStats& s, double g, bool reml) {
  BhfEval e;
  e.w = s.nd / (1.0 + s.nd * g);
  mat A = s.Wxx + s.xbar.t() * (s.xbar.each_col() % e.w);
  const vec b = s.Wxy + s.xbar.t() * (e.w % s.ybar);
  if (!chol(e.Rchol, symmatu(A))) return e;
  const mat Rt = trimatl(e.Rchol.t());
  e.beta = solve(trimatu(e.Rchol), solve(Rt, b, solve_opts::fast), solve_opts::fast);
  e.logdetA = 2.0 * accu(log(e.Rchol.diag()));
  e.rbar = s.ybar - s.xbar * e.beta;
  const double Wrr = std::max(0.0, s.Wyy - 2.0 * dot(e.beta, s.Wxy) +
                                       as_scalar(e.beta.t() * s.Wxx * e.beta));
  e.Q = Wrr + accu(e.w % square(e.rbar));
  if (!(e.Q > 0) || !std::isfinite(e.Q)) return e;
  const double ldV = accu(log1p(s.nd * g));
  const vec w2 = square(e.w);
  const double dQ = -accu(w2 % square(e.rbar));  // envelope: beta at its optimum
  if (reml) {
    const double df = s.n - s.p;
    const mat T = solve(Rt, s.xbar.t(), solve_opts::fast);  // p x S
    const rowvec q = sum(square(T), 0);                      // xbar_d' A^-1 xbar_d
    e.ll = -0.5 * (df * std::log(e.Q) + ldV + e.logdetA);
    e.grad = -0.5 * (df * dQ / e.Q + accu(e.w) - dot(w2, q.t()));
  } else {
    e.ll = -0.5 * (s.n * std::log(e.Q) + ldV);
    e.grad = -0.5 * (s.n * dQ / e.Q + accu(e.w));
  }
  e.ok = std::isfinite(e.ll) && std::isfinite(e.grad);
  return e;
}

// Brent's root finder on the score in [a, b] with grad(a) > 0 >= grad(b)
double bhf_zeroin(const BhfStats& s, bool reml, double a, double fa, double b, double fb) {
  const double EPS = std::numeric_limits<double>::epsilon();
  double c = a, fc = fa;
  for (int it = 0; it < 200; ++it) {
    const double prev = b - a;
    if (std::fabs(fc) < std::fabs(fb)) {
      a = b; b = c; c = a;
      fa = fb; fb = fc; fc = fa;
    }
    const double tol = 2.0 * EPS * std::fabs(b) + 1e-300;
    const double m = 0.5 * (c - b);
    if (std::fabs(m) <= tol || fb == 0.0) return b;
    double step = m;
    if (std::fabs(prev) >= tol && std::fabs(fa) > std::fabs(fb)) {
      double p, q, t = (c - b);
      const double s1 = fb / fa;
      if (a == c) {
        p = t * s1;
        q = 1.0 - s1;
      } else {
        const double q1 = fa / fc, r = fb / fc;
        p = s1 * (t * q1 * (q1 - r) - (b - a) * (r - 1.0));
        q = (q1 - 1.0) * (r - 1.0) * (s1 - 1.0);
      }
      if (p > 0) q = -q; else p = -p;
      if (p < (0.75 * t * q - std::fabs(tol * q) / 2) && p < std::fabs(prev * q / 2)) step = p / q;
    }
    if (std::fabs(step) < tol) step = (step > 0) ? tol : -tol;
    a = b; fa = fb;
    b += step;
    const BhfEval eb = bhf_eval(s, b, reml);
    fb = eb.ok ? eb.grad : -1.0;  // failure: treat as beyond the root
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
      c = a; fc = fa;
    }
  }
  return b;
}

struct BhfFit {
  bool ok = false, converged = false;
  double gamma = 0, s2e = 0, s2u = 0;
  BhfEval e;
};

// Global maximiser of the profiled criterion over gamma >= 0: scan a log grid
// for sign changes of the score, refine each with Brent, keep the best; the
// boundary gamma = 0 is a candidate when the score there is non-positive.
BhfFit bhf_fit(const BhfStats& s, bool reml) {
  BhfFit f;
  const BhfEval e0 = bhf_eval(s, 0.0, reml);
  if (!e0.ok) return f;
  double best_g = -1, best_ll = -INFINITY;
  if (e0.grad <= 0) { best_g = 0; best_ll = e0.ll; }
  double ga = 0, fa = e0.grad;
  bool top_open = false;
  for (int k = 0; k <= 32; ++k) {
    const double gb = std::pow(10.0, -8.0 + 0.5 * k);
    const BhfEval eb = bhf_eval(s, gb, reml);
    if (!eb.ok) { ga = gb; fa = -1; continue; }
    if (fa > 0 && eb.grad <= 0) {
      const double r = bhf_zeroin(s, reml, ga, fa, gb, eb.grad);
      const BhfEval er = bhf_eval(s, r, reml);
      if (er.ok && er.ll > best_ll) { best_ll = er.ll; best_g = r; }
    }
    ga = gb; fa = eb.grad;
    if (k == 32 && eb.grad > 0) {
      top_open = true;
      if (eb.ll > best_ll) { best_ll = eb.ll; best_g = gb; }
    }
  }
  if (best_g < 0) return f;
  f.e = bhf_eval(s, best_g, reml);
  if (!f.e.ok) return f;
  f.ok = true;
  f.converged = !top_open || best_g < std::pow(10.0, 8.0);
  f.gamma = best_g;
  f.s2e = f.e.Q / (reml ? (s.n - s.p) : s.n);
  f.s2u = f.gamma * f.s2e;
  return f;
}

// Domain bookkeeping shared by the fit and the bootstrap
struct BhfDesign {
  int n = 0, p = 0, D = 0;
  std::vector<int> sidx;  // sampled domain ids (0-based, ascending)
  ivec pos;               // domain id -> position among sampled (-1 if none)
  ivec spos;              // unit -> position among sampled domains
  vec nd_all;             // D sample sizes
  mat Xc;                 // within-domain centred X
  BhfStats base;          // y-independent part of the statistics
};

BhfDesign bhf_design(const mat& Xs, const ivec& dom, int D) {
  BhfDesign d;
  d.n = Xs.n_rows; d.p = Xs.n_cols; d.D = D;
  d.nd_all.zeros(D);
  for (int j = 0; j < d.n; ++j) d.nd_all(dom(j)) += 1;
  d.pos.set_size(D); d.pos.fill(-1);
  for (int k = 0; k < D; ++k) if (d.nd_all(k) > 0) { d.pos(k) = d.sidx.size(); d.sidx.push_back(k); }
  const int S = d.sidx.size();
  d.spos.set_size(d.n);
  d.base.xbar.zeros(S, d.p);
  d.base.nd.set_size(S);
  for (int i = 0; i < S; ++i) d.base.nd(i) = d.nd_all(d.sidx[i]);
  for (int j = 0; j < d.n; ++j) {
    d.spos(j) = d.pos(dom(j));
    d.base.xbar.row(d.spos(j)) += Xs.row(j);
  }
  d.base.xbar.each_col() /= d.base.nd;
  d.Xc = Xs;
  for (int j = 0; j < d.n; ++j) d.Xc.row(j) -= d.base.xbar.row(d.spos(j));
  d.base.Wxx = d.Xc.t() * d.Xc;
  d.base.n = d.n;
  d.base.p = d.p;
  return d;
}

// Fill the y-dependent statistics
void bhf_ystats(const BhfDesign& d, const vec& y, BhfStats& s) {
  const int S = d.sidx.size();
  s.ybar.zeros(S);
  for (int j = 0; j < d.n; ++j) s.ybar(d.spos(j)) += y(j);
  s.ybar /= s.nd;
  s.Wxy = d.Xc.t() * y;
  double w = 0;
  for (int j = 0; j < d.n; ++j) { const double r = y(j) - s.ybar(d.spos(j)); w += r * r; }
  s.Wyy = w;
}

// EBLUP of the population means (sampled: f ybar + (Xbar - f xbar)' b + (1 - f) u;
// unsampled: Xbar' b)
vec bhf_eblup(const BhfDesign& d, const BhfStats& s, const BhfFit& f,
              const mat& meanxpop, const vec& popnsize) {
  vec out = meanxpop * f.e.beta;
  for (size_t i = 0; i < d.sidx.size(); ++i) {
    const int k = d.sidx[i];
    const double fd = s.nd(i) / popnsize(k);
    const double u = f.gamma * f.e.w(i) * f.e.rbar(i);
    out(k) = fd * s.ybar(i) - fd * dot(s.xbar.row(i), f.e.beta) + out(k) + (1.0 - fd) * u;
  }
  return out;
}

}  // namespace

// [[Rcpp::export(.bhf_fit_cpp)]]
Rcpp::List bhf_fit_cpp(const arma::mat& Xs, const arma::vec& ys, const arma::ivec& dom,
                       int D, bool reml) {
  BhfDesign d = bhf_design(Xs, dom, D);
  BhfStats s = d.base;
  bhf_ystats(d, ys, s);
  const BhfFit f = bhf_fit(s, reml);
  if (!f.ok) Rcpp::stop("BHF variance component estimation failed.");
  vec u(D); u.fill(NA_REAL);
  for (size_t i = 0; i < d.sidx.size(); ++i) u(d.sidx[i]) = f.gamma * f.e.w(i) * f.e.rbar(i);
  const mat Rinv = inv(trimatu(f.e.Rchol));
  const mat vcov = f.s2e * (Rinv * Rinv.t());
  return Rcpp::List::create(
    Rcpp::_["beta"] = f.e.beta,
    Rcpp::_["vcov"] = vcov,
    Rcpp::_["sigma2_u"] = f.s2u,
    Rcpp::_["sigma2_e"] = f.s2e,
    Rcpp::_["u"] = u,
    Rcpp::_["convergence"] = f.converged);
}

// Parametric bootstrap MSE (Gonzalez-Manteiga et al. 2008 style) of the BHF
// EBLUP. Random draws come from R's RNG in a fixed order per replicate
// (u: D, e: n, non-sampled mean error: D), so results depend on the seed only,
// not on the number of threads.
// [[Rcpp::export(.pbmse_bhf_cpp)]]
Rcpp::List pbmse_bhf_cpp(const arma::mat& Xs, const arma::ivec& dom, int D,
                         const arma::mat& meanxpop, const arma::vec& popnsize,
                         const arma::vec& beta, double sigma2_u, double sigma2_e,
                         int B, bool reml, int n_threads) {
  const BhfDesign d = bhf_design(Xs, dom, D);
  const int n = d.n;
  const vec Xb = Xs * beta;
  const vec mu = meanxpop * beta;
  const double su = std::sqrt(sigma2_u), se = std::sqrt(sigma2_e);
  // population mean error: ebar_d = (n_d ebar_{d,s} + r_d ebar_{d,r}) / N_d
  vec w_s(D), w_r(D), sd_r(D);
  for (int k = 0; k < D; ++k) {
    const double nk = d.nd_all(k), Nk = popnsize(k), rk = std::max(Nk - nk, 0.0);
    if (nk > 0 && rk > 0) { w_s(k) = nk / Nk; w_r(k) = rk / Nk; sd_r(k) = std::sqrt(sigma2_e / rk); }
    else if (nk > 0) { w_s(k) = 1; w_r(k) = 0; sd_r(k) = 0; }
    else { w_s(k) = 0; w_r(k) = 1; sd_r(k) = std::sqrt(sigma2_e / Nk); }
  }

  int nt = 1;
#ifdef _OPENMP
  nt = std::max(1, n_threads);
#endif
  const int chunk = std::max(1, std::min(B, 32 * nt));
  mat U(D, chunk), E(n, chunk), ER(D, chunk), SQ(D, chunk);
  std::vector<int> valid(chunk);
  vec mse(D, fill::zeros);
  int valid_B = 0;

  for (int b0 = 0; b0 < B; b0 += chunk) {
    const int nb = std::min(chunk, B - b0);
    for (int c = 0; c < nb; ++c) {  // R's RNG: sequential, master thread only
      for (int k = 0; k < D; ++k) U(k, c) = R::rnorm(0.0, su);
      for (int j = 0; j < n; ++j) E(j, c) = R::rnorm(0.0, se);
      for (int k = 0; k < D; ++k) ER(k, c) = R::rnorm(0.0, sd_r(k));
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nt)
#endif
    for (int c = 0; c < nb; ++c) {
      vec y(n);
      for (int j = 0; j < n; ++j) y(j) = Xb(j) + U(dom(j), c) + E(j, c);
      BhfStats s = d.base;
      bhf_ystats(d, y, s);
      const BhfFit f = bhf_fit(s, reml);
      valid[c] = f.ok;
      if (!f.ok) continue;
      const vec est = bhf_eblup(d, s, f, meanxpop, popnsize);
      vec esd(D, fill::zeros);
      for (int j = 0; j < n; ++j) esd(dom(j)) += E(j, c);
      for (int k = 0; k < D; ++k) if (d.nd_all(k) > 0) esd(k) /= d.nd_all(k);
      const vec truth = mu + U.col(c) + w_s % esd + w_r % ER.col(c);
      SQ.col(c) = square(est - truth);
    }
    for (int c = 0; c < nb; ++c) if (valid[c]) { mse += SQ.col(c); ++valid_B; }
    Rcpp::checkUserInterrupt();
  }
  if (valid_B > 0) mse /= valid_B; else mse.fill(NA_REAL);
  return Rcpp::List::create(Rcpp::_["mse"] = mse, Rcpp::_["B"] = valid_B);
}
