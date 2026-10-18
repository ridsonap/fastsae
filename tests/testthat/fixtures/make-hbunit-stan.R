# Builds hbunit-stan.rds: Stan posteriors of the unit-level models that
# hb_unit(engine = jax_engine()) fits, with the same priors and parameterisation
# as the NumPyro backend (inst/python/hb_area_numpyro.py), the reference for
# test_hb_unit_engine.R. Agreement with it checks the implementation; the priors
# themselves differ from INLA's and hbsae's.
# Run from the package root: Rscript tests/testthat/fixtures/make-hbunit-stan.R
# (a few minutes; needs cmdstanr and CmdStan).

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> D;
  int<lower=1, upper=3> fam;          // 1 gaussian, 2 bernoulli, 3 poisson
  matrix[N, P] X;
  array[N] int<lower=1, upper=D> dom;
  vector[N] y;
  array[N] int y_int;
  matrix[D, P] Xp;                    // prediction rows: population covariate means
  real<lower=0> beta_scale;
  int<lower=0, upper=1> bym2;
  int<lower=0> K;
  matrix[D, K] basis;                 // unit-variance scaled ICAR basis
  real y_c;
  real<lower=0> y_s;
  vector[P] int_mask;
}
parameters {
  vector[P] beta;
  real<lower=0> sigma_s;
  vector[D] z_s;
  vector[K] z_icar;
  real<lower=0, upper=1> phi;
  real<lower=0> sigma_e;
}
transformed parameters {
  vector[D] u;
  if (bym2 == 1) {
    u = sigma_s * (sqrt(1 - phi) * z_s + sqrt(phi) * (basis * z_icar));
  } else {
    u = sigma_s * z_s;
  }
}
model {
  beta ~ normal(0, beta_scale);
  sigma_s ~ normal(0, 1);
  z_s ~ std_normal();
  z_icar ~ std_normal();
  phi ~ beta(1, 1);
  sigma_e ~ normal(0, 1);
  {
    vector[N] eta = X * beta + u[dom];
    if (fam == 1) {
      (y - y_c) / y_s ~ normal(eta, sigma_e);
    } else if (fam == 2) {
      y_int ~ bernoulli_logit(eta);
    } else {
      y_int ~ poisson_log(eta);
    }
  }
}
generated quantities {
  vector[D] mu;
  real sigma2_u;
  real sigma2_e;
  vector[P] beta_y;
  {
    vector[D] eta_p = Xp * beta + u;
    if (fam == 1) {
      mu = y_c + y_s * eta_p;
      sigma2_u = square(y_s * sigma_s);
      sigma2_e = square(y_s * sigma_e);
      beta_y = y_s * beta + y_c * int_mask;
    } else {
      mu = fam == 2 ? inv_logit(eta_p) : exp(eta_p);
      sigma2_u = square(sigma_s);
      sigma2_e = not_a_number();
      beta_y = beta;
    }
  }
}
"

devtools::load_all(quiet = TRUE)

# 30 domains in a 5 x 6 grid, three of them without a sample
make_data <- function(seed = 11) {
  set.seed(seed)
  D <- 30
  mx <- rnorm(D)
  u <- rnorm(D)
  pop <- data.frame(dom = seq_len(D), x1 = mx, N = sample(200:800, D))
  n_d <- sample(3:12, D, replace = TRUE)
  n_d[c(4, 17, 29)] <- 0
  s <- do.call(rbind, lapply(seq_len(D), function(d) {
    if (n_d[d] == 0) {
      return(NULL)
    }
    data.frame(dom = d, x1 = rnorm(n_d[d], mx[d], 1))
  }))
  s$y_gauss <- 2 + 1.5 * s$x1 + u[s$dom] + rnorm(nrow(s), sd = 2)
  s$y_bin <- rbinom(nrow(s), 1, plogis(-0.5 + 0.8 * s$x1 + 0.6 * u[s$dom]))
  s$y_pois <- rpois(nrow(s), exp(0.3 + 0.4 * s$x1 + 0.4 * u[s$dom]))
  grid <- expand.grid(r = 1:5, c = 1:6)
  W <- 1 * (as.matrix(dist(grid, method = "manhattan")) == 1)
  dimnames(W) <- NULL
  list(sample = s, pop = pop, W = W)
}

dat <- make_data()
s <- dat$sample
X <- cbind(1, s$x1)
Xp <- cbind(1, dat$pop$x1)
spec <- .icar_spectral(dat$W)
pos <- spec$eig_values > 1e-6
basis <- spec$eig_vectors[, pos] / rep(sqrt(spec$eig_values[pos] * spec$scale_factor), each = nrow(dat$W))

model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stan_code))

fit_ref <- function(response, fam, bym2 = FALSE) {
  y <- s[[response]]
  gauss <- fam == 1
  # The backend's standardisation: centred (model with intercept) and divided by
  # the population SD of the observed responses
  sdat <- list(
    N = nrow(s), P = 2, D = nrow(dat$pop), fam = fam, X = X, dom = s$dom,
    y = as.numeric(y), y_int = if (gauss) rep(0L, nrow(s)) else as.integer(y), Xp = Xp,
    beta_scale = if (fam == 2) 2.5 else 10, bym2 = as.integer(bym2),
    K = if (bym2) ncol(basis) else 0L,
    basis = if (bym2) basis else matrix(0, nrow(dat$pop), 0),
    y_c = if (gauss) mean(y) else 0, y_s = if (gauss) sqrt(mean((y - mean(y))^2)) else 1,
    int_mask = c(1, 0)
  )
  fit <- model$sample(
    data = sdat, chains = 4, parallel_chains = 4, iter_warmup = 1000,
    iter_sampling = 5000, seed = 7, refresh = 0, adapt_delta = 0.95,
    show_messages = FALSE
  )
  pars <- c("mu", "sigma2_u", "sigma2_e", "beta_y")
  smr <- fit$summary(pars, "rhat", "ess_bulk")
  info <- list(
    max_rhat = max(smr$rhat, na.rm = TRUE),
    min_ess = min(smr$ess_bulk, na.rm = TRUE),
    divergent = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
  )
  stopifnot(info$max_rhat < 1.01, info$min_ess > 1000, info$divergent == 0)
  dr <- fit$draws(pars, format = "draws_matrix")
  col <- function(p) dr[, grep(paste0("^", p, "(\\[|$)"), colnames(dr)), drop = FALSE]
  mu <- unclass(col("mu"))
  list(
    mean = colMeans(mu), sd = apply(mu, 2, sd),
    sigma2_u = mean(col("sigma2_u")), sigma2_e = mean(col("sigma2_e")),
    beta = colMeans(unclass(col("beta_y"))), info = info
  )
}

ref <- list(
  gaussian = fit_ref("y_gauss", 1),
  gaussian_bym2 = fit_ref("y_gauss", 1, bym2 = TRUE),
  binomial = fit_ref("y_bin", 2),
  poisson = fit_ref("y_pois", 3)
)
str(lapply(ref, `[[`, "info"))
saveRDS(list(data = dat, ref = ref), "tests/testthat/fixtures/hbunit-stan.rds", version = 2)
