library(fastsae)

sim_ratio_data <- function(D = 40, seed = 7) {
  set.seed(seed)
  d <- data.frame(area = sprintf("a%02d", D:1), x = rnorm(D), pop = rpois(D, 1e4), prov = rep(1:4, length.out = D))
  d$p0 <- plogis(-2 + 0.3 * rnorm(D))
  d$p0_se <- 0.15 * d$p0
  d$n_poor <- rpois(D, 30) + 1
  d$s_dir <- rbinom(D, d$n_poor, plogis(-1.5 + 0.5 * d$x + rnorm(D, 0, 0.3))) / d$n_poor
  d
}

test_that("hb_ratio estimates stay below the parent and match the sub-proportion model", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("sn")
  d <- sim_ratio_data()
  d$s_dir[3] <- NA # unsampled area
  fit <- hb_ratio(s_dir ~ x, data = d, trials = "n_poor", parent = "p0", domain = "area",
                  n_draws = 500, print_result = FALSE)
  df <- fit$df_hb
  expect_s3_class(fit, c("fastsae_hb_ratio", "fastsae_hb", "fastsae"))
  expect_identical(as.character(df$domain), d$area)
  expect_true(all(df$hb <= d$p0) && all(df$ci_upper <= d$p0))
  expect_equal(df$hb, df$ratio * d$p0)
  expect_equal(df$rse, df$rse_cond)
  expect_equal(df$rse_cond, 100 * df$ratio_sd / df$ratio)
  expect_true(is.na(df$y[3]) && !is.na(df$hb[3]))
  expect_equal(df$y[-3], d$s_dir[-3] * d$p0[-3])
  # the posterior mean of s agrees with the hb_area fit up to Monte Carlo error
  expect_lt(max(abs(df$ratio - fit$ratio_fit$df_hb$hb) / fit$ratio_fit$df_hb$sd), 0.25)
  expect_equal(dim(fit$draws), c(nrow(d), 500L))
  # reproducible up to the run-to-run noise of INLA (~1e-8)
  fit2 <- hb_ratio(s_dir ~ x, data = d, trials = "n_poor", parent = "p0", domain = "area",
                   n_draws = 500, print_result = FALSE)
  expect_equal(fit2$df_hb$hb, df$hb, tolerance = 1e-5)
  expect_message(print(fit), "Ratio model")
})

test_that("hb_ratio benchmarks on the logit scale and adds the parent's uncertainty", {
  skip_on_cran()
  skip_if_not_installed("INLA")
  skip_if_not_installed("sn")
  d <- sim_ratio_data()
  fit <- hb_ratio(s_dir ~ x, data = d, trials = "n_poor", parent = "p0", parent_se = "p0_se",
                  domain = "area", benchmark_target = 0.03, benchmark_weight = "pop",
                  n_draws = 300, print_result = FALSE)
  df <- fit$df_hb
  expect_equal(sum(d$pop * df$hb) / sum(d$pop), 0.03, tolerance = 1e-8)
  expect_equal(unname(fit$benchmark$achieved), 0.03, tolerance = 1e-8)
  expect_true(all(df$hb <= d$p0))
  expect_true(all(df$rse > df$rse_cond))
  # group targets
  tg <- c("1" = 0.02, "2" = 0.03, "3" = 0.025, "4" = 0.035)
  fg <- hb_ratio(s_dir ~ x, data = d, trials = "n_poor", parent = "p0", domain = "area",
                 benchmark_target = tg, benchmark_weight = "pop", benchmark_group = "prov",
                 n_draws = 200, print_result = FALSE)
  ach <- tapply(d$pop * fg$df_hb$hb, d$prov, sum) / tapply(d$pop, d$prov, sum)
  expect_equal(as.numeric(ach), as.numeric(tg), tolerance = 1e-8)
})

test_that("hb_ratio checks its input", {
  d <- sim_ratio_data(D = 10)
  args <- list(formula = s_dir ~ x, data = d, trials = "n_poor", parent = "p0", print_result = FALSE)
  run <- function(...) do.call(hb_ratio, utils::modifyList(args, list(...)))
  expect_error(run(family = "poisson"), "binomial")
  expect_error(run(self_benchmark = TRUE), "benchmark_target")
  expect_error(run(engine = "jax"), "INLA engine")
  expect_error(run(n_draws = 5), "n_draws")
  d_bad <- d
  d_bad$p0[1] <- -0.1
  expect_error(run(data = d_bad), "non-negative")
  expect_error(run(parent_se = rep(NA, 10)), "parent_se")
  expect_error(run(benchmark_target = 0.02), "benchmark_weight")
  expect_error(run(benchmark_target = 0.5, benchmark_weight = "pop"), "below the weighted mean")
  expect_error(run(benchmark_target = c("1" = 0.02), benchmark_weight = "pop", benchmark_group = "prov"), "named by every group")
  d_w <- d
  d_w$pop[] <- 0
  expect_error(run(data = d_w, benchmark_target = 0.02, benchmark_weight = "pop"), "not all zero")
})
