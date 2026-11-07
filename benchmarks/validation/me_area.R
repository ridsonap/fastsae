# Validation of hb_area(var_x = ) (covariates measured with error) against Stan:
# the same likelihood and priors (me_area.stan, where the mean and SD of the true
# values of the structural model are integrated out), against
#   default:    hb_area() as shipped (INLA's classic mode; structural: mean and SD
#               of the true values at their posterior mode given the observed values),
#   compact:    the same in INLA's compact mode (non-Gaussian models only),
#   full_bayes: the structural model with the mean and SD as two more INLA
#               hyperparameters (INLA integrates over four, with CCD).
# Run from the package root: Rscript benchmarks/validation/me_area.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
vdir <- "benchmarks/validation"
mod <- cmdstanr::cmdstan_model(file.path(vdir, "me_area.stan"), quiet = TRUE)
dir.create(file.path(vdir, "out"), showWarnings = FALSE)

# full_bayes: free the mean and precision of the true values in the mec terms
me_terms_shipped <- get(".me_inla_terms", asNamespace("fastsae"))
me_terms_full_bayes <- function(me, ...) {
  out <- me_terms_shipped(me, ...)
  for (k in seq_along(out$terms)) {
    wg <- out$terms[[k]]$values
    vw <- stats::var(wg)
    cg <- 1 / out$terms[[k]]$scale
    out$terms[[k]]$hyper$mean.x <- list(prior = "gaussian", param = c(mean(wg), 1e-4), initial = mean(wg), fixed = FALSE)
    out$terms[[k]]$hyper$prec.x <- c(.half_normal_prec(1), list(initial = -log(max(vw - mean(cg), 0.05 * vw)), fixed = FALSE))
  }
  out
}
with_me_terms <- function(f, expr) {
  utils::assignInNamespace(".me_inla_terms", f, "fastsae")
  on.exit(utils::assignInNamespace(".me_inla_terms", me_terms_shipped, "fastsae"))
  expr
}

simulate <- function(fam, D, K, seed, n_out = 5) {
  set.seed(seed)
  X <- matrix(rnorm(D * K, 5, 2), D, K)
  C <- matrix(runif(D * K, 0.2, 1.5), D, K)
  W_obs <- X + matrix(rnorm(D * K, 0, sqrt(C)), D, K)
  z <- rnorm(D)
  bx <- c(0.8, -0.5)[seq_len(K)]
  eff <- rnorm(D, 0, 0.7)
  d <- data.frame(area = seq_len(D), z = z)
  for (k in seq_len(K)) {
    d[[paste0("w", k)]] <- W_obs[, k]
    d[[paste0("c", k)]] <- C[, k]
  }
  if (fam == "gaussian") {
    d$vardir <- runif(D, 0.5, 1.5)
    d$y <- 1 + 0.5 * z + drop(X %*% bx) + eff + rnorm(D, 0, sqrt(d$vardir))
  } else {
    d$n <- rpois(D, 30) + 5
    d$y <- rbinom(D, d$n, plogis(-1 + 0.3 * z + drop((X - 5) %*% (bx / 4)) + eff / 2))
  }
  d$y[sort(sample(D, n_out))] <- NA
  d
}

compare <- function(fam, D, K, model, seed) {
  label <- sprintf("%s D=%d K=%d", fam, D, K)
  d <- simulate(fam, D, K, seed)
  s <- !is.na(d$y)
  wn <- paste0("w", seq_len(K))
  cn <- paste0("c", seq_len(K))
  variants <- list(default = list(engine = list(), terms = me_terms_shipped))
  if (fam != "gaussian") variants$compact <- list(engine = list(inla.mode = "compact"), terms = me_terms_shipped)
  if (model == "structural") variants$full_bayes <- list(engine = list(), terms = me_terms_full_bayes)
  fits <- list()
  times <- c()
  for (v in names(variants)) {
    args <- list(
      formula = stats::reformulate(c("z", wn), response = "y"), data = d, domain = "area",
      family = fam, var_x = stats::setNames(as.list(cn), wn), me_model = model,
      engine = do.call(inla_engine, variants[[v]]$engine), print_result = FALSE
    )
    if (fam == "gaussian") args$vardir <- "vardir" else args$trials <- "n"
    t0 <- Sys.time()
    fits[[v]] <- with_me_terms(variants[[v]]$terms, suppressWarnings(do.call(hb_area, args)))
    times[v] <- as.numeric(Sys.time() - t0, units = "secs")
  }
  W <- as.matrix(d[wn])
  sd_ <- list(
    D = D, p = 1, Z = cbind(d$z), K = K, W = W, C = as.matrix(d[cn]),
    structural = as.integer(model == "structural"), Ns = sum(s), idx = which(s),
    fam = if (fam == "gaussian") 1L else 2L,
    y = if (fam == "gaussian") d$y[s] else rep(0, sum(s)),
    yi = if (fam == "gaussian") rep(0L, sum(s)) else as.integer(d$y[s]),
    psi = if (fam == "gaussian") d$vardir[s] else rep(1, sum(s)),
    ntr = if (fam == "gaussian") rep(1L, sum(s)) else as.integer(d$n[s]),
    sd_scale = if (fam == "gaussian") stats::sd(d$y, na.rm = TRUE) else 1,
    w_mean = colMeans(W), w_sd = apply(W, 2, stats::sd)
  )
  fit <- mod$sample(
    data = sd_, chains = 4, parallel_chains = 4, iter_warmup = 2000, iter_sampling = 4000,
    refresh = 0, show_messages = FALSE, seed = 1, adapt_delta = 0.95,
    output_dir = file.path(vdir, "out")
  )
  dr <- fit$draws("theta", format = "matrix")
  sm <- colMeans(dr)
  ss <- apply(dr, 2, stats::sd)
  bx <- fit$draws("bx", format = "matrix")
  div <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
  ess <- min(posterior::summarise_draws(dr, "ess_bulk")$ess_bulk)
  out <- lapply(names(fits), function(v) {
    f <- fits[[v]]
    z <- (f$df_hb$hb - sm) / ss
    r <- f$df_hb$sd / ss
    zb <- (f$estcoef[wn, "beta"] - colMeans(bx)) / apply(bx, 2, stats::sd)
    rb <- f$estcoef[wn, "std.error"] / apply(bx, 2, stats::sd)
    data.frame(
      config = label, model = model, seed = seed, inla = v,
      max_dmean_sd = max(abs(z)), mean_dmean_sd = mean(abs(z)),
      sd_ratio_min = min(r), sd_ratio_max = max(r),
      beta_max_dmean_sd = max(abs(zb)), beta_sd_ratio_min = min(rb), beta_sd_ratio_max = max(rb),
      time_s = times[[v]], divergent = div, ess_min = ess
    )
  })
  res <- do.call(rbind, out)
  print(res, digits = 3, row.names = FALSE)
  res
}

res <- rbind(
  compare("gaussian", 30, 1, "structural", 1),
  compare("gaussian", 30, 1, "functional", 1),
  compare("gaussian", 50, 1, "structural", 1),
  compare("gaussian", 50, 1, "structural", 2),
  compare("gaussian", 50, 1, "functional", 1),
  compare("gaussian", 100, 1, "structural", 1),
  compare("gaussian", 100, 1, "functional", 1),
  compare("gaussian", 50, 2, "structural", 1),
  compare("gaussian", 50, 2, "functional", 1),
  compare("binomial", 50, 1, "structural", 1),
  compare("binomial", 50, 1, "structural", 2),
  compare("binomial", 50, 1, "functional", 1),
  compare("binomial", 100, 1, "structural", 1)
)
saveRDS(res, file.path(vdir, "results_me.rds"))
options(width = 200)
print(res, digits = 3, row.names = FALSE)
