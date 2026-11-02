#!/usr/bin/env Rscript

# ==============================================================================
# PIT of diagnose() for count likelihoods (binomial and Poisson hb_area() fits).
#
# 1. INLA's CPO against refitting without the area (exact leave-one-out up to
#    INLA's approximation of the refit), for a few areas.
# 2. Replicated correct models: rejection rate of the KS test and share of PIT
#    values in the central 90% band for INLA's unrandomized PIT (fastsae 1.0.0),
#    INLA's randomized PIT and the randomized mixed predictive PIT (default).
#
# Usage (from the package root, needs INLA and sn; PIT_REPS=0 runs part 1 only):
#   R CMD INSTALL . --preclean && Rscript benchmarks/validation/pit_discrete.R
#
# Result (2026-11-01, INLA 26.8.7, ~25 minutes):
# 1. Binomial (300 areas, 20 trials, 68% zeros): INLA CPO 0.027, 0.103, 0.133 for three
#    zero counts where refitting gives 0.711, 0.583, 0.634 (mixed predictive 0.706, 0.578,
#    0.637). Poisson: the three areas agree (0.768 / 0.770 / 0.768).
# 2. 20 replicates, KS rejection at 5% and mean central-90% share [range]:
#      binomial  INLA PIT            1.00   0.774 [0.497, 0.873]
#                INLA, randomized    0.45   0.777 [0.370, 0.930]
#                mixed, randomized   0.00   0.902 [0.863, 0.927]
#      poisson   INLA PIT            0.95   0.831 [0.257, 0.940]
#                INLA, randomized    0.90   0.800 [0.230, 0.923]
#                mixed, randomized   0.00   0.900 [0.873, 0.930]
#    INLA's CPO fails in some replicates of both families (central share down to 0.23),
#    with no CPO failure flagged; the mixed predictive is calibrated in all of them.
# ==============================================================================

suppressPackageStartupMessages({
  library(fastsae)
  library(INLA)
})
INLA::inla.setOption(num.threads = "1:1")
R <- as.integer(Sys.getenv("PIT_REPS", "20"))

sim_binomial <- function(seed, D = 300, n = 20) {
  set.seed(seed)
  x <- rnorm(D)
  p <- plogis(-4 + 0.5 * x + rnorm(D, 0, 0.3))
  data.frame(area = seq_len(D), y = rbinom(D, n, p) / n, n = n, x = x)
}
sim_poisson <- function(seed, D = 300) {
  set.seed(seed)
  x <- rnorm(D)
  E <- runif(D, 2, 10)
  data.frame(area = seq_len(D), y = rpois(D, E * exp(-2 + 0.4 * x + rnorm(D, 0, 0.4))), E = E, x = x)
}
fit_binomial <- function(d) suppressMessages(hb_area(y ~ x, data = d, family = "binomial", trials = "n", domain = "area", print_result = FALSE))
fit_poisson <- function(d) suppressMessages(hb_area(y ~ x, data = d, family = "poisson", exposure = "E", domain = "area", print_result = FALSE))

# ---- 1. CPO against refitting ---------------------------------------------------
loo <- function(d, fitter, dens, ids) {
  f <- fitter(d)
  mp <- diagnose(f, pit_method = "mixed")$df_diag$cpo
  do.call(rbind, lapply(ids, function(i) {
    d2 <- d
    d2$y[i] <- NA
    s <- inla.posterior.sample(4000, fitter(d2)$fit, selection = list(Predictor = i), seed = 1)
    eta <- vapply(s, function(z) z$latent[1], numeric(1))
    data.frame(area = i, y = d$y[i], inla_cpo = f$fit$cpo$cpo[i], mixed = mp[i], refit = mean(dens(i, eta)))
  }))
}
db <- sim_binomial(112)
cat("== Binomial, D = 300, n = 20: predictive probability of the observed count\n")
print(loo(db, fit_binomial, function(i, eta) dbinom(round(db$y[i] * 20), 20, plogis(eta)), c(5, 7, 8, 1, 2)), digits = 3)
dp <- sim_poisson(3)
cat("== Poisson, D = 300\n")
print(loo(dp, fit_poisson, function(i, eta) dpois(dp$y[i], dp$E[i] * exp(eta)), which(dp$y <= 1)[1:3]), digits = 3)

# ---- 2. replicated correct models ------------------------------------------------
summ <- function(u) c(ks_p = suppressWarnings(ks.test(u, "punif")$p.value), central90 = mean(u >= 0.05 & u <= 0.95))
rep_one <- function(r, sim, fitter) {
  f <- fitter(sim(1000 + r))
  cp <- f$fit$cpo
  set.seed(r)
  c(
    inla = summ(cp$pit[cp$pit > 0 & cp$pit < 1]),
    inla_rand = summ(cp$pit - runif(length(cp$pit)) * cp$cpo),
    mixed = summ(diagnose(f, pit_method = "mixed", seed = r)$df_diag$pit)
  )
}
for (fam in if (R > 0) c("binomial", "poisson")) {
  res <- t(vapply(seq_len(R), rep_one, numeric(6),
    sim = if (fam == "binomial") sim_binomial else sim_poisson,
    fitter = if (fam == "binomial") fit_binomial else fit_poisson
  ))
  cat(sprintf("\n== %s, %d replicates: KS rejection at 5%% and central-90%% share (mean, range)\n", fam, R))
  for (m in c("inla", "inla_rand", "mixed")) {
    cat(sprintf("  %-10s reject %.2f  central90 %.3f [%.3f, %.3f]\n", m,
      mean(res[, paste0(m, ".ks_p")] < 0.05), mean(res[, paste0(m, ".central90")]),
      min(res[, paste0(m, ".central90")]), max(res[, paste0(m, ".central90")])))
  }
}
