# Validation of hb_area() against Stan: Leroux, SLM and negative binomial.
# Run from the package root: Rscript benchmarks/validation/run_spatial.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
vdir <- "benchmarks/validation"
mod <- cmdstanr::cmdstan_model(file.path(vdir, "area.stan"), quiet = TRUE)
dir.create(file.path(vdir, "out"), showWarnings = FALSE)

D <- 49
W <- as.matrix(unclass(sim_spatial_weights(D, type = "grid", seed = 1)))
A <- (W > 0 | t(W) > 0) * 1
diag(A) <- 0
ed <- which(upper.tri(A) & A > 0, arr.ind = TRUE)
R <- diag(rowSums(A)) - A
gscale <- exp(mean(log(diag(MASS::ginv(R)))))
Wst <- W / rowSums(W)
ev <- Re(eigen(Wst, only.values = TRUE)$values)
rho_min <- max(1 / min(ev), -0.999)
rho_max <- min(1 / max(ev), 0.999)

simulate <- function(fam, seed) {
  set.seed(seed)
  x1 <- rnorm(D)
  eff <- 0.5 * rnorm(D)
  d <- data.frame(area = seq_len(D), x1 = x1)
  if (fam == "gaussian") {
    d$vardir <- runif(D, 0.5, 2)
    d$y <- 10 + 2 * x1 + eff + rnorm(D, sd = sqrt(d$vardir))
  }
  if (fam == "poisson") {
    d$E <- runif(D, 5, 50)
    d$y <- rpois(D, d$E * exp(-1 + 0.4 * x1 + eff))
  }
  if (fam == "nbinomial") {
    d$E <- runif(D, 5, 50)
    d$y <- rnbinom(D, size = 5, mu = d$E * exp(-1 + 0.4 * x1 + eff))
  }
  d$y[sort(sample(D, 12))] <- NA
  d
}

compare <- function(label, fam, spatial, seed, leroux_unif = 0L, adapt_delta = 0.99,
                    iter_sampling = 4000, inla_extra = list()) {
  d <- simulate(fam, seed)
  s <- !is.na(d$y)
  # the PC prior of the Stan model (area.stan), not the half-normal default of fastsae >= 1.0.0
  extra <- list(prior_prec = list(prior = "pc.prec", param = c(1, 0.01)))
  if (spatial == "leroux") {
    # uniform (the default since fastsae 1.0.0) or INLA's own default, as in area.stan
    extra$prior_rho <- if (leroux_unif == 1L) list(prior = "logitbeta", param = c(1, 1)) else list(prior = "gaussian", param = c(0, 0.1))
  }
  if (spatial == "bym2") extra$prior_phi <- list(prior = "logitbeta", param = c(1, 1))
  if (fam == "nbinomial") extra$control.family <- list(hyper = list(theta = list(prior = "loggamma", param = c(1, 0.01))))
  extra <- c(extra, inla_extra)
  fits <- list()
  times <- c()
  for (mode in c("default", "compact", "classic")) {
    eng_args <- if (mode == "default") extra else c(extra, list(inla.mode = mode))
    args <- list(
      formula = y ~ x1, data = d, domain = "area", family = fam, spatial = spatial,
      engine = do.call(inla_engine, eng_args), print_result = FALSE
    )
    if (spatial != "none") args$W <- W
    if (fam == "gaussian") args$vardir <- "vardir" else args$exposure <- "E"
    t0 <- Sys.time()
    fits[[mode]] <- suppressWarnings(do.call(hb_area, args))
    times[mode] <- as.numeric(Sys.time() - t0, units = "secs")
  }
  sd_ <- list(
    D = D, p = 2, X = cbind(1, d$x1), Ns = sum(s), idx = which(s),
    fam = match(fam, c("gaussian", "binomial", "poisson", "beta", "gamma", "nbinomial")),
    y = if (fam == "gaussian") d$y[s] else rep(0, sum(s)),
    yi = if (fam == "gaussian") rep(0L, sum(s)) else as.integer(d$y[s]),
    psi = if (fam == "gaussian") d$vardir[s] else rep(1, sum(s)), ntr = rep(1L, sum(s)),
    logE = if (fam == "gaussian") rep(0, sum(s)) else log(d$E[s]), sc = rep(1, sum(s)),
    re = match(spatial, c("none", "besag", "bym2", "leroux", "slm")), Ne = nrow(ed), n1 = ed[, 1], n2 = ed[, 2],
    gscale = gscale, R = R, Wst = Wst, rho_min = rho_min, rho_max = rho_max, lam = -log(0.01),
    leroux_unif = leroux_unif, nb_a = 1, nb_b = 0.01
  )
  fit <- mod$sample(
    data = sd_, chains = 4, parallel_chains = 4, iter_warmup = 2000, iter_sampling = iter_sampling,
    refresh = 0, show_messages = FALSE, seed = 1, adapt_delta = adapt_delta,
    output_dir = file.path(vdir, "out")
  )
  dr <- fit$draws("theta", format = "matrix")
  sm <- colMeans(dr)
  ss <- apply(dr, 2, sd)
  div <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
  ess <- min(posterior::summarise_draws(dr, "ess_bulk")$ess_bulk)
  out <- lapply(names(fits), function(mode) {
    f <- fits[[mode]]
    z <- (f$df_hb$hb - sm) / ss
    r <- f$df_hb$sd / ss
    data.frame(
      config = label, seed = seed, mode = mode, max_dmean_sd = max(abs(z)), mean_dmean_sd = mean(abs(z)),
      sd_ratio_min = min(r), sd_ratio_max = max(r), time_s = times[[mode]], divergent = div, ess_min = ess
    )
  })
  res <- do.call(rbind, out)
  print(res, digits = 3, row.names = FALSE)
  res
}

res <- rbind(
  compare("leroux gaussian", "gaussian", "leroux", 1),
  compare("leroux gaussian", "gaussian", "leroux", 2),
  compare("leroux gaussian U(0,1)", "gaussian", "leroux", 1, leroux_unif = 1L),
  compare("leroux poisson", "poisson", "leroux", 1),
  compare("slm gaussian", "gaussian", "slm", 1),
  compare("slm gaussian", "gaussian", "slm", 2),
  compare("slm poisson", "poisson", "slm", 1),
  compare("nbinomial iid", "nbinomial", "none", 1, iter_sampling = 15000),
  compare("nbinomial bym2", "nbinomial", "bym2", 1, iter_sampling = 15000),
  compare("nbinomial iid, INLA grid", "nbinomial", "none", 1,
    iter_sampling = 15000,
    inla_extra = list(control.inla = list(int.strategy = "grid", dz = 0.5, diff.logdens = 8))
  )
)
saveRDS(res, file.path(vdir, "results_spatial.rds"))
