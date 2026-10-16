// Reference model for fastsae::hb_area(): the same likelihood, latent effect and
// priors as the INLA fit, sampled with NUTS.
data {
  int<lower=1> D;
  int<lower=1> p;
  matrix[D, p] X;                    // first column is the intercept
  int<lower=0> Ns;
  array[Ns] int<lower=1, upper=D> idx;
  int<lower=1, upper=6> fam;         // 1 gaussian (FH), 2 binomial, 3 poisson, 4 beta, 5 gamma, 6 nbinomial
  vector[Ns] y;
  array[Ns] int yi;
  vector[Ns] psi;                    // gaussian: sampling variances
  array[Ns] int ntr;                 // binomial: trials
  vector[Ns] logE;                   // poisson, nbinomial: log exposure
  vector[Ns] sc;                     // beta: phi_d; gamma: shape s_d
  int<lower=1, upper=5> re;          // 1 iid, 2 scaled besag, 3 bym2, 4 leroux, 5 slm
  int<lower=0> Ne;
  array[Ne] int n1;
  array[Ne] int n2;
  real<lower=0> gscale;              // geometric mean of the ICAR marginal variances
  matrix[D, D] R;                    // ICAR structure matrix (leroux)
  matrix[D, D] Wst;                  // row-standardised W (slm)
  real rho_min;
  real rho_max;
  real<lower=0> lam;                 // PC prior rate on sigma
  int<lower=0, upper=1> leroux_unif; // 1: beta ~ U(0, 1); 0: logit(beta) ~ N(0, sd 1/sqrt(0.1))
  real<lower=0> nb_a;                // nbinomial: size ~ gamma(nb_a, nb_b)
  real<lower=0> nb_b;
}
parameters {
  real b0;
  vector[p - 1] b;
  real<lower=0> sigma;
  vector[D] z;
  vector[D] u;
  real<lower=0, upper=1> phi;        // bym2 mixing, leroux beta or slm rho*
  real<lower=0> size;
}
transformed parameters {
  vector[D] eta;
  {
    vector[D] fixed = b0 + X[, 2:p] * b;
    if (re == 1) eta = fixed + sigma * z;
    else if (re == 2) eta = fixed + sigma * u / sqrt(gscale);
    else if (re == 3) eta = fixed + sigma * (sqrt(1 - phi) * z + sqrt(phi) * u / sqrt(gscale));
    else if (re == 4) {
      matrix[D, D] L = cholesky_decompose(add_diag(phi * R, 1 - phi));
      eta = fixed + sigma * mdivide_right_tri_low(z', L)';
    } else {
      real rho = rho_min + phi * (rho_max - rho_min);
      eta = mdivide_left(add_diag(-rho * Wst, 1.0), fixed + sigma * z);
    }
  }
}
model {
  sigma ~ exponential(lam);
  if (re == 5) {             // slm: all coefficients inside the latent field, Q.beta = 1e-4
    b0 ~ normal(0, 100);
    b ~ normal(0, 100);
  } else {                   // INLA defaults: flat intercept, precision 0.001
    b ~ normal(0, sqrt(1000));
  }
  z ~ std_normal();
  if (re == 2 || re == 3) {
    target += -0.5 * dot_self(u[n1] - u[n2]);
    sum(u) ~ normal(0, 0.001 * D);
  } else {
    u ~ std_normal();
  }
  if (re == 4 && leroux_unif == 0) {
    target += normal_lpdf(logit(phi) | 0, 1 / sqrt(0.1)) - log(phi) - log1m(phi);
  }
  size ~ gamma(nb_a, nb_b);
  vector[Ns] es = eta[idx];
  if (fam == 1) y ~ normal(es, sqrt(psi));
  else if (fam == 2) yi ~ binomial_logit(ntr, es);
  else if (fam == 3) yi ~ poisson_log(logE + es);
  else if (fam == 4) { vector[Ns] mu = inv_logit(es); y ~ beta(mu .* sc, (1 - mu) .* sc); }
  else if (fam == 5) y ~ gamma(sc, sc ./ exp(es));
  else yi ~ neg_binomial_2_log(logE + es, size);
}
generated quantities {
  vector[D] theta;
  if (fam == 1) theta = eta;
  else if (fam == 2 || fam == 4) theta = inv_logit(eta);
  else theta = exp(eta);
}
