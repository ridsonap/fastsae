# Validation of the E-BFMI of jax_engine() fits against CmdStan: hb_area() and
# ebfmi_area.stan fit the same non-centered iid area model, on the same standardised
# scale and with the same priors and unconstraining transforms, so the energies of
# the two samplers are comparable (the E-BFMI depends on the parameterization).
# Same NUTS settings: 4 chains, 1000 + 1000 iterations, adapt_delta = target_accept
# = 0.85, maximum tree depth 10. Three models over 20 seeds:
#   vardir:    Fay-Herriot on mys, known sampling variances (well behaved),
#   no_vardir: mys with sigma2_e estimated (weakly identified; neither sampler converges),
#   binomial:  iid logit model on simulated counts. The backend's clipping of eta at
#              +-8 is not reproduced; it has no effect on these data.
# The posterior means are compared too, as a check that both target the same model.
# Run from the package root: Rscript benchmarks/validation/ebfmi_stan.R
Sys.setenv(JAX_PLATFORMS = "cpu")
suppressMessages(devtools::load_all(".", quiet = TRUE))
vdir <- "benchmarks/validation"
mod <- cmdstanr::cmdstan_model(file.path(vdir, "ebfmi_area.stan"), quiet = TRUE)
data("mys", package = "fastsae")

reps <- 20L
chains <- 4L
warmup <- 1000L
samples <- 1000L
delta <- 0.85

set.seed(101)
Db <- 60L
bin <- data.frame(area = seq_len(Db), x1 = rnorm(Db), n = sample(20:200, Db, TRUE))
bin$cases <- rbinom(Db, bin$n, plogis(-0.8 + 0.5 * bin$x1 + rnorm(Db, 0, 0.5)))

# The data of the backend's Gaussian models are standardised: y_s is the SD with
# denominator n, the prior scale of the area SD uses the SD with denominator n - 1
stan_data <- function(y, X, vardir = NULL, tr = NULL) {
  obs <- which(!is.na(y))
  gauss <- is.null(tr)
  y_c <- if (gauss) mean(y[obs]) else 0
  y_s <- if (gauss) sqrt(mean((y[obs] - mean(y[obs]))^2)) else 1
  sd_prior <- if (gauss) sd(y[obs]) else 1
  list(
    stan = list(
      D = nrow(X), P = ncol(X), X = X, Ns = length(obs), idx = as.array(obs),
      fam = if (gauss) 1L else 2L,
      has_vardir = as.integer(!is.null(vardir)),
      y_g = as.array(if (gauss) (y[obs] - y_c) / y_s else rep(0, length(obs))),
      psi = as.array(if (!is.null(vardir)) vardir[obs] / y_s^2 else rep(1, length(obs))),
      y_b = as.array(if (gauss) rep(0L, length(obs)) else as.integer(y[obs])),
      tr = as.array(if (gauss) rep(1L, length(obs)) else as.integer(tr[obs])),
      beta_scale = sqrt(1000) / y_s, sd_scale = sd_prior / y_s
    ),
    y_c = y_c, y_s = y_s, coef = colnames(X)
  )
}

Xm <- model.matrix(~ x1 + x2 + x3, mys)
Xb <- model.matrix(~x1, bin)
cases <- list(
  vardir = list(
    fit = function(e) hb_area(y ~ x1 + x2 + x3, vardir = "vardir", data = mys, engine = e, print_result = FALSE),
    data = stan_data(mys$y, Xm, vardir = mys$vardir)
  ),
  no_vardir = list(
    fit = function(e) hb_area(y ~ x1 + x2 + x3, data = mys, engine = e, print_result = FALSE),
    data = stan_data(mys$y, Xm)
  ),
  binomial = list(
    fit = function(e) {
      hb_area(cases ~ x1, data = bin, domain = "area", family = "binomial", trials = "n", engine = e, print_result = FALSE)
    },
    data = stan_data(bin$cases, Xb, tr = bin$n)
  )
)

fit_stan <- function(dat, seed) {
  sf <- mod$sample(
    data = dat$stan, seed = seed, chains = chains, parallel_chains = chains,
    iter_warmup = warmup, iter_sampling = samples, adapt_delta = delta, max_treedepth = 10,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE, diagnostics = NULL
  )
  diag <- sf$sampler_diagnostics(inc_warmup = FALSE)
  energy <- posterior::extract_variable_matrix(diag, "energy__")
  has_e <- dat$stan$fam == 1L && dat$stan$has_vardir == 0L
  draws <- posterior::as_draws_df(sf$draws(c("beta", "sigma_s", if (has_e) "sigma_e"), format = "draws_matrix"))
  beta <- as.matrix(draws[, grep("^beta\\[", names(draws))]) * dat$y_s
  beta[, 1] <- beta[, 1] + dat$y_c
  colnames(beta) <- dat$coef
  post <- cbind(beta, sigma2_u = (draws$sigma_s * dat$y_s)^2)
  if (has_e) post <- cbind(post, sigma2_e = (draws$`sigma_e[1]` * dat$y_s)^2)
  list(
    energy = unclass(energy),
    ebfmi = as.numeric(sf$diagnostic_summary(quiet = TRUE)$ebfmi),
    divergences = sum(posterior::extract_variable_matrix(diag, "divergent__")),
    post = post
  )
}

ebfmi_rows <- post_rows <- list()
for (nm in names(cases)) {
  for (r in seq_len(reps)) {
    seed <- 1000L + r
    eng <- jax_engine(
      warmup = warmup, samples = samples, chains = chains, device = "cpu", seed = seed,
      progress = FALSE, target_accept = delta
    )
    fj <- suppressWarnings(suppressMessages(cases[[nm]]$fit(eng)))
    st <- fit_stan(cases[[nm]]$data, seed)
    if (r == 1L) {
      # the backend's function on Stan's energies must give cmdstanr's values
      backend <- .get_numpyro_backend(device = "cpu")
      np <- reticulate::import("numpy", convert = FALSE)
      py <- as.numeric(reticulate::py_to_r(backend$ebfmi(np$asarray(t(st$energy)))))
      cat(sprintf("%s: backend ebfmi() vs cmdstanr on the same energies, max |diff| = %.1e\n", nm, max(abs(py - st$ebfmi))))
    }
    ebfmi_rows[[length(ebfmi_rows) + 1L]] <- data.frame(
      model = nm, sampler = c(rep("numpyro", chains), rep("stan", chains)), rep = r,
      chain = rep(seq_len(chains), 2), ebfmi = c(fj$ebfmi, st$ebfmi),
      divergences = c(rep(fj$divergences, chains), rep(st$divergences, chains))
    )
    tab <- fj$mcmc_diagnostics
    pars <- colnames(st$post)
    post_rows[[length(post_rows) + 1L]] <- data.frame(
      model = nm, rep = r, par = pars, stan = colMeans(st$post),
      numpyro = tab$mean[match(pars, tab$parameter)], sd_stan = apply(st$post, 2, sd), row.names = NULL
    )
    cat(nm, r, "numpyro", round(mean(fj$ebfmi), 3), "stan", round(mean(st$ebfmi), 3), "\n")
  }
}
ebfmi <- do.call(rbind, ebfmi_rows)
post <- do.call(rbind, post_rows)

# Summary: chains for the distribution, the fit (mean of its chains) as the replicate
# of the comparison, since chains of one fit share the data and the seed
summary_tab <- do.call(rbind, lapply(names(cases), function(nm) {
  d <- ebfmi[ebfmi$model == nm, ]
  fit <- aggregate(ebfmi ~ rep + sampler, d, mean)
  tt <- t.test(fit$ebfmi[fit$sampler == "numpyro"], fit$ebfmi[fit$sampler == "stan"])
  low <- aggregate(ebfmi ~ rep + sampler, d, function(v) any(v < 0.3))
  data.frame(
    model = nm,
    mean_numpyro = mean(d$ebfmi[d$sampler == "numpyro"]), sd_numpyro = sd(d$ebfmi[d$sampler == "numpyro"]),
    mean_stan = mean(d$ebfmi[d$sampler == "stan"]), sd_stan = sd(d$ebfmi[d$sampler == "stan"]),
    diff = unname(tt$estimate[1] - tt$estimate[2]), ci_low = tt$conf.int[1], ci_high = tt$conf.int[2], p = tt$p.value,
    fits_low_numpyro = sum(low$ebfmi[low$sampler == "numpyro"]), fits_low_stan = sum(low$ebfmi[low$sampler == "stan"]),
    divergences_numpyro = mean(d$divergences[d$sampler == "numpyro"]), divergences_stan = mean(d$divergences[d$sampler == "stan"])
  )
}))
post$z <- (post$numpyro - post$stan) / post$sd_stan
post_tab <- aggregate(z ~ model + par, post, function(v) c(mean = mean(v), max_abs = max(abs(v))))

saveRDS(list(ebfmi = ebfmi, posterior = post, summary = summary_tab), file.path(vdir, "results_ebfmi.rds"))
options(width = 200)
print(summary_tab, digits = 3, row.names = FALSE)
print(post_tab, digits = 3, row.names = FALSE)
