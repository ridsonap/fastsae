# Validation of the area-mean MSE of eblup_twofold() (df_area).
#
# Design: m = 40 areas with 2-6 sub-areas, about 20% of the sub-areas without a
# direct estimate and 3 areas without any; X, vardir and the population weights
# are fixed, the random effects and sampling errors are redrawn in every
# replicate. Target: theta_i = sum_j w_ij theta_ij over the sub-areas in the data.
# Reported per estimator: relative bias of the mean MSE estimate against the
# empirical MSE, and the coverage of eblup +/- 1.96 sqrt(mse).
#
# Run from the package root: Rscript benchmarks/validation/twofold_area_mse.R
#
# Result (2026-10-18; 1000 replicates, bootstrap 300 replicates with B = 200):
#   estimator          relative bias (sampled / no sample)  coverage     sum w^2 mse
#   REML analytical    -1.2% / -2.9%                        94.6 / 93.5  -47% / -63%
#   ML analytical      +0.2% / +2.4%                        94.9 / 93.6  -47% / -61%
#   REML bootstrap     +1.8% / -6.9%                        94.9 / 93.8  -46% / -65%
# The last column is the relative bias of sum_j w_ij^2 mse_ij, which ignores the
# covariance of the sub-area errors within an area.

devtools::load_all(quiet = TRUE)

set.seed(2026)
m <- 40
s2v <- 0.6
s2u <- 0.3
design <- do.call(rbind, lapply(seq_len(m), function(d) {
  nd <- sample(2:6, 1)
  data.frame(area = d, subarea = paste0(d, "-", seq_len(nd)),
             x1 = rnorm(nd), vardir = runif(nd, 0.2, 1.2), pop = rpois(nd, 200) + 10)
}))
n <- nrow(design)
miss <- runif(n) < 0.2 | design$area %in% c(7, 21, 33)
w <- design$pop / ave(design$pop, design$area, FUN = sum)
xb <- 1 + 0.8 * design$x1

run <- function(R, method, mse_type, B = 100) {
  est <- mse <- naive <- truth <- matrix(NA_real_, m, R)
  for (r in seq_len(R)) {
    theta <- xb + rnorm(m, sd = sqrt(s2v))[design$area] + rnorm(n, sd = sqrt(s2u))
    dat <- design
    dat$y <- theta + rnorm(n, sd = sqrt(design$vardir))
    dat$y[miss] <- NA
    fit <- suppressWarnings(eblup_twofold(
      y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea", data = dat,
      weight = "pop", method = method, mse_type = mse_type, B = B, print_result = FALSE
    ))
    est[, r] <- fit$df_area$eblup_area
    mse[, r] <- fit$df_area$mse_area
    # sum_j w_ij^2 mse_ij ignores the covariance of the sub-area errors
    naive[, r] <- as.numeric(rowsum(w^2 * fit$df_eblup$mse, design$area))
    truth[, r] <- as.numeric(rowsum(w * theta, design$area))
  }
  emp <- rowMeans((est - truth)^2)
  cover <- rowMeans(abs(est - truth) <= 1.96 * sqrt(mse), na.rm = TRUE)
  sampled <- !(seq_len(m) %in% c(7, 21, 33))
  summ <- function(v) c(sampled = mean(v[sampled]), unsampled = mean(v[!sampled]))
  rb <- rowMeans(mse, na.rm = TRUE) / emp - 1
  rb_naive <- rowMeans(naive, na.rm = TRUE) / emp - 1
  data.frame(
    method = method, mse_type = mse_type,
    group = c("sampled areas", "areas without sample"),
    rel_bias_mean = summ(rb), rel_bias_median = c(median(rb[sampled]), median(rb[!sampled])),
    coverage = summ(cover), rel_bias_naive = summ(rb_naive),
    row.names = NULL
  )
}

t0 <- Sys.time()
res <- rbind(
  run(1000, "REML", "analytical"),
  run(1000, "ML", "analytical"),
  run(300, "REML", "bootstrap", B = 200)
)
print(res, digits = 3)
cat("elapsed:", format(Sys.time() - t0), "\n")
