// Reference for fastsae::hb_area() with a space-time interaction. Each intrinsic
// effect is written as B z with z ~ N(0, I) and B = V+ Lambda+^(-1/2) from the
// eigen-decomposition of its structure matrix: this is the Gaussian with that
// precision under hard constraints spanning its null space (sum-to-zero for the
// main effects, the Knorr-Held constraints for the interaction).
data {
  int<lower=1> N;
  int<lower=1> p;
  matrix[N, p] X;
  int<lower=1> D;
  int<lower=2> Tt;
  array[N] int<lower=1, upper=D> dom;
  array[N] int<lower=1, upper=Tt> tim;
  int<lower=0> Ns;
  array[Ns] int idx;
  int<lower=1, upper=2> fam;            // 1 gaussian (FH), 2 poisson
  vector[Ns] y;
  array[Ns] int yi;
  vector[Ns] psi;
  vector[Ns] logE;
  int<lower=1> rs; matrix[D, rs] Bs;    // spatial main effect basis
  int<lower=1> rt; matrix[Tt, rt] Bt;   // temporal main effect basis
  int<lower=0, upper=1> trend;          // rw2 main effect: unpenalised linear trend
  vector[Tt] tc;
  int<lower=1> Nint;
  array[N] int<lower=1, upper=Nint> iint;
  int<lower=1> ri; matrix[Nint, ri] Bi; // interaction basis
  real<lower=0> lam;
}
parameters {
  real b0;
  vector[p - 1] b;
  real<lower=0> sig_s;
  real<lower=0> sig_t;
  real<lower=0> sig_i;
  vector[rs] zs;
  vector[rt] zt;
  vector[ri] zi;
  real a_trend;
}
transformed parameters {
  vector[Tt] ut = sig_t * (Bt * zt) + (trend == 1 ? a_trend * tc : rep_vector(0, Tt));
  vector[N] eta = b0 + X[, 2:p] * b + sig_s * (Bs * zs)[dom] + ut[tim] + sig_i * (Bi * zi)[iint];
}
model {
  sig_s ~ exponential(lam);
  sig_t ~ exponential(lam);
  sig_i ~ exponential(lam);
  b ~ normal(0, sqrt(1000));
  zs ~ std_normal();
  zt ~ std_normal();
  zi ~ std_normal();
  if (trend == 0) a_trend ~ std_normal();   // unused
  if (fam == 1) y ~ normal(eta[idx], sqrt(psi));
  else yi ~ poisson_log(logE + eta[idx]);
}
generated quantities {
  vector[N] theta = fam == 1 ? eta : exp(eta);
}
