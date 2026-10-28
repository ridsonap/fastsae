# Binomial models with effective sample sizes (hb_area(), hb_twofold(); family = "binomial",
# trials = n / deff): the effective numbers of cases n_eff * p_hat are not integers, and INLA
# fits them with its "xbinomial" likelihood (Chen, Wakefield and Lumley, 2014; Franco and
# Bell, 2015).
#
# A. INLA against long MCMC runs (jax_engine(), NumPyro NUTS, 4 x 5000 draws, same model and
#    priors) on one simulated data set each for hb_area() (40 areas) and hb_twofold()
#    (12 provinces, 4-8 districts each), with non-integer effective sizes and some zero
#    direct estimates. Differences of the posterior means in posterior SD and ratios of the
#    posterior SDs.
# B. Design-based coverage with a two-stage cluster sample from each district (4-12 PSUs of
#    10 units; PSU prevalence logit-normal around the district, SD 0.6, so deff 1.1-1.7),
#    12 provinces x 6 districts, 200 replicates. hb_twofold() with
#      n        trials = sample size (design ignored),
#      true     trials = n / deff with the true design effect,
#      pooled   trials = n / deff with deff estimated from the PSU means and pooled per
#               province (the design effects of single districts are too noisy, and
#               undefined when p_hat is 0).
#    Truth: the district mean of the PSU prevalences, and the population-weighted province
#    mean.
#
# Run from the package root: Rscript benchmarks/validation/binomial_neff.R (~15 minutes; part
# A needs Python with jax and numpyro, see setup_numpyro_env())
#
# Result (2026-10-29):
# A. |difference of posterior means| / posterior SD, median (max); SD ratio MCMC / INLA
#      hb_area               0.004 (0.078)   0.95-1.02
#      hb_twofold district   0.025 (0.063)   0.96-1.04
#      hb_twofold province   0.041 (0.078)   0.96-1.05
# B. coverage of 95% intervals (mean width; RMSE), 0.3 zero direct estimates per replicate
#                 district                  province
#      n          0.911 (0.141; 0.042)      0.897 (0.073; 0.023)
#      true       0.929 (0.147; 0.040)      0.946 (0.086; 0.023)
#      pooled     0.925 (0.148; 0.040)      0.937 (0.086; 0.023)
#    Ignoring the design (trials = n) gives intervals that are too short, more so for the
#    province aggregates; effective sample sizes with the true or the pooled design effects
#    bring the coverage close to 95% with the same accuracy.
devtools::load_all(quiet = TRUE)
q <- function(e) suppressWarnings(suppressMessages(e))

# ---- A. INLA vs MCMC ---------------------------------------------------------------
part_a <- function() {
  if (!q(check_numpyro_available())) {
    message("NumPyro not available: part A skipped")
    return(NULL)
  }
  eng <- jax_engine(warmup = 2000, samples = 5000, chains = 4, device = "cpu", progress = FALSE)
  cmp <- function(fi, fj, est, sd) {
    z <- abs(fi[[est]] - fj[[est]]) / fi[[sd]]
    c(median_diff_sd = stats::median(z), max_diff_sd = max(z),
      sd_ratio_min = min(fj[[sd]] / fi[[sd]]), sd_ratio_max = max(fj[[sd]] / fi[[sd]]))
  }
  set.seed(101)
  D <- 40
  n <- sample(25:200, D, replace = TRUE)
  da <- data.frame(area = seq_len(D), x = rnorm(D), n_eff = n / runif(D, 1.2, 2.5))
  da$p_hat <- rbinom(D, n, plogis(-1.5 + 0.6 * da$x + rnorm(D, sd = 0.5))) / n
  da$p_hat[c(3, 17)] <- 0
  args <- list(p_hat ~ x, data = da, domain = "area", family = "binomial", trials = "n_eff", print_result = FALSE)
  fi <- q(do.call(hb_area, args))
  fj <- q(do.call(hb_area, c(args, list(engine = eng))))
  out <- rbind(hb_area = cmp(fi$df_hb, fj$df_hb, "hb", "sd"))

  set.seed(102)
  prov <- rep(1:12, times = sample(4:8, 12, replace = TRUE))
  n <- sample(25:200, length(prov), replace = TRUE)
  dt <- data.frame(prov = prov, kab = seq_along(prov), x = rnorm(length(prov)),
                   pop = runif(length(prov), 1e4, 1e5), n_eff = n / runif(length(prov), 1.2, 2.5))
  dt$p_hat <- rbinom(length(prov), n, plogis(-1.5 + 0.6 * dt$x + rnorm(12, sd = 0.4)[prov] +
                                               rnorm(length(prov), sd = 0.3))) / n
  dt$p_hat[c(2, 9)] <- 0
  args <- list(p_hat ~ x, data = dt, domain = "prov", subarea = "kab", weight = "pop",
               family = "binomial", trials = "n_eff", print_result = FALSE)
  fi <- q(do.call(hb_twofold, args))
  fj <- q(do.call(hb_twofold, c(args, list(engine = eng))))
  rbind(out,
        hb_twofold_district = cmp(fi$df_hb, fj$df_hb, "hb", "sd"),
        hb_twofold_province = cmp(fi$df_area, fj$df_area, "hb_area", "sd_area"))
}

# ---- B. design-based coverage ------------------------------------------------------
P <- 12
K <- 6
b <- 10 # units per PSU
s_c <- 0.6 # SD of the PSU prevalences on the logit scale
psu_moments <- function(mu) {
  # E[p_c] and E[p_c (1 - p_c)] for logit(p_c) ~ N(mu, s_c^2), by quadrature
  z <- stats::qnorm(stats::ppoints(2000))
  pc <- stats::plogis(mu + s_c * z)
  c(mean = mean(pc), var_between = stats::var(pc), within = mean(pc * (1 - pc)))
}
set.seed(2026)
design <- data.frame(prov = rep(seq_len(P), each = K), kab = seq_len(P * K))
design$x <- rnorm(P * K)
design$a <- sample(4:12, P * K, replace = TRUE) # PSUs sampled
design$n <- design$a * b
design$pop <- round(runif(P * K, 2e4, 2e5))
design$mu <- -1.5 + 0.5 * design$x + rnorm(P, sd = 0.4)[design$prov] + rnorm(P * K, sd = 0.3)
mom <- t(vapply(design$mu, psu_moments, numeric(3)))
design$theta <- mom[, "mean"]
# Var(p_hat) = (Var(p_c) + E[p_c (1 - p_c)] / b) / a for a sample with replacement of PSUs
v_true <- (mom[, "var_between"] + mom[, "within"] / b) / design$a
design$deff <- v_true / (design$theta * (1 - design$theta) / design$n)
theta_prov <- tapply(design$theta * design$pop, design$prov, sum) / tapply(design$pop, design$prov, sum)

one_rep <- function(r) {
  set.seed(5e5 + r)
  d <- design
  p_psu <- lapply(seq_len(nrow(d)), function(i) {
    pc <- stats::plogis(d$mu[i] + s_c * rnorm(d$a[i]))
    rbinom(d$a[i], b, pc) / b
  })
  d$p_hat <- vapply(p_psu, mean, 0)
  v_hat <- vapply(p_psu, function(m) stats::var(m) / length(m), 0)
  deff_hat <- v_hat / (d$p_hat * (1 - d$p_hat) / d$n)
  deff_hat[!is.finite(deff_hat)] <- NA
  # pooled per province: ratio of the summed design variances to the summed SRS variances
  srs <- d$p_hat * (1 - d$p_hat) / d$n
  pool <- tapply(v_hat, d$prov, sum) / tapply(srs, d$prov, sum)
  d$deff_pooled <- pmax(pool[as.character(d$prov)], 1)
  out <- NULL
  for (cfg in c("n", "true", "pooled")) {
    d$n_eff <- switch(cfg, n = d$n, true = d$n / d$deff, pooled = d$n / d$deff_pooled)
    f <- tryCatch(q(hb_twofold(p_hat ~ x,
      domain = "prov", subarea = "kab", weight = "pop", data = d,
      family = "binomial", trials = "n_eff", print_result = FALSE,
      engine = inla_engine(num.threads = "1:1")
    )), error = function(e) NULL)
    if (is.null(f)) next
    h <- f$df_hb
    a <- f$df_area
    out <- rbind(out,
      data.frame(cfg = cfg, level = "district", rep = r,
                 cover = d$theta >= h$ci_lower & d$theta <= h$ci_upper,
                 width = h$ci_upper - h$ci_lower, err = h$hb - d$theta),
      data.frame(cfg = cfg, level = "province", rep = r,
                 cover = theta_prov >= a$ci_lower_area & theta_prov <= a$ci_upper_area,
                 width = a$ci_upper_area - a$ci_lower_area, err = a$hb_area - theta_prov)
    )
  }
  if (!is.null(out)) out$zeros <- sum(d$p_hat == 0)
  out
}

res_a <- part_a()
res_b <- do.call(rbind, parallel::mclapply(1:200, one_rep, mc.cores = max(1, parallel::detectCores() - 2)))
summ_b <- do.call(rbind, lapply(split(res_b, list(res_b$cfg, res_b$level), drop = TRUE), function(s) {
  data.frame(cfg = s$cfg[1], level = s$level[1], coverage = mean(s$cover),
             mean_width = mean(s$width), rmse = sqrt(mean(s$err^2)), bias = mean(s$err))
}))
cat("design effects (min, median, max):", sprintf("%.2f", stats::quantile(design$deff, c(0, 0.5, 1))), "\n")
cat("zero direct estimates per replicate:", mean(res_b$zeros[res_b$cfg == "n" & res_b$level == "province"]), "\n")
print(res_a, digits = 3)
print(summ_b[order(summ_b$level, summ_b$cfg), ], digits = 3, row.names = FALSE)
