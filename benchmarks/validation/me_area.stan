// Reference model for fastsae::hb_area(var_x = ): an area-level model whose
// covariates W are measured with error (W = X + e, e ~ N(0, C), C known), with
// the likelihood and priors of the INLA fit, sampled with NUTS.
// structural = 1: X ~ N(mu, sx^2), mu ~ N(mean w, (100 sd w)^2), sx ~ half-normal(sd w);
// structural = 0 (functional): a flat prior on X.
data {
  int<lower=1> D;
  int<lower=0> p;                    // error-free covariates besides the intercept
  matrix[D, p] Z;
  int<lower=1> K;                    // covariates measured with error
  matrix[D, K] W;                    // observed values
  matrix<lower=0>[D, K] C;           // error variances (> 0)
  int<lower=0, upper=1> structural;
  int<lower=0> Ns;
  array[Ns] int<lower=1, upper=D> idx;
  int<lower=1, upper=2> fam;         // 1 gaussian (Fay-Herriot), 2 binomial
  vector[Ns] y;
  array[Ns] int yi;
  vector[Ns] psi;                    // gaussian: sampling variances
  array[Ns] int ntr;                 // binomial: trials
  real<lower=0> sd_scale;            // half-normal scale of the area SD
  vector[K] w_mean;
  vector<lower=0>[K] w_sd;
}
parameters {
  real b0;
  vector[p] b;
  vector[K] bx;
  real<lower=0> sigma;
  vector[D] z;
  matrix[D, K] zx;
  vector[K] mu;
  vector<lower=0>[K] sx;
}
transformed parameters {
  matrix[D, K] X;
  for (k in 1:K) {
    if (structural == 1) {
      // X | W: the conditional of the true values given the observed ones (as INLA's mec)
      vector[D] prec = 1 / square(sx[k]) + 1 ./ col(C, k);
      vector[D] m = (mu[k] / square(sx[k]) + col(W, k) ./ col(C, k)) ./ prec;
      X[, k] = m + col(zx, k) ./ sqrt(prec);
    } else {
      X[, k] = col(W, k) - sqrt(col(C, k)) .* col(zx, k);
    }
  }
  vector[D] eta = b0 + Z * b + X * bx + sigma * z;
}
model {
  // INLA defaults: flat intercept, precision 0.001; the same prior for bx
  b ~ normal(0, sqrt(1000));
  bx ~ normal(0, sqrt(1000));
  sigma ~ normal(0, sd_scale);
  z ~ std_normal();
  to_vector(zx) ~ std_normal();
  if (structural == 1) {
    mu ~ normal(w_mean, 100 * w_sd);
    sx ~ normal(0, w_sd);
    // marginal of the observed values given mu and sx
    for (k in 1:K) col(W, k) ~ normal(mu[k], sqrt(square(sx[k]) + col(C, k)));
  } else {
    // not in the functional model; proper priors keep them sampled
    mu ~ std_normal();
    sx ~ std_normal();
  }
  if (fam == 1) y ~ normal(eta[idx], sqrt(psi));
  else yi ~ binomial_logit(ntr, eta[idx]);
}
generated quantities {
  vector[D] theta = fam == 1 ? eta : inv_logit(eta);
  vector[K] sigma2_x = square(sx);
}
