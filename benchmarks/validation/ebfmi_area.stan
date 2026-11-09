// The non-centered iid area model of the NumPyro backend (fit_numpyro_hb), on its
// standardised scale: same parameters, same priors, same unconstraining transforms.
data {
  int<lower=1> D;
  int<lower=1> P;
  matrix[D, P] X;
  int<lower=1> Ns;
  array[Ns] int<lower=1, upper=D> idx;
  int<lower=1, upper=2> fam;           // 1 gaussian, 2 binomial
  int<lower=0, upper=1> has_vardir;
  vector[Ns] y_g;
  vector[Ns] psi;
  array[Ns] int y_b;
  array[Ns] int tr;
  real<lower=0> beta_scale;
  real<lower=0> sd_scale;
}
parameters {
  vector[P] beta;
  real<lower=0> sigma_s;
  vector[D] z_s;
  array[(fam == 1 && has_vardir == 0) ? 1 : 0] real<lower=0> sigma_e;
}
model {
  vector[D] eta = X * beta + sigma_s * z_s;
  beta ~ normal(0, beta_scale);
  sigma_s ~ normal(0, sd_scale);
  z_s ~ std_normal();
  if (fam == 1) {
    if (has_vardir == 1) {
      y_g ~ normal(eta[idx], sqrt(psi));
    } else {
      sigma_e[1] ~ normal(0, 1);
      y_g ~ normal(eta[idx], sigma_e[1]);
    }
  } else {
    y_b ~ binomial_logit(tr, eta[idx]);
  }
}
