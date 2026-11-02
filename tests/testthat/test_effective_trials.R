library(fastsae)

et_data <- function() {
  data.frame(
    prov = c(1, 1, 1, 2, 2, 3),
    y = c(0.10, 0.20, 0, 0.05, 0.15, 0),
    n = c(100, 50, 80, 200, 40, 30),
    fb = c(NA, NA, 2.5, NA, NA, NA)
  ) |>
    transform(vardir = c(1.5, 2, 0, 1.2, 1.8, 0) * y * (1 - y) / n)
}

test_that("effective_trials pools design effects per group", {
  d <- et_data()
  e <- effective_trials(~y, vardir = ~vardir, n = ~n, group = ~prov, data = d)
  expect_named(e, c("n", "deff", "deff_source", "trials"))
  srs <- d$y * (1 - d$y) / d$n
  ok <- d$y > 0
  p1 <- sum(d$vardir[1:2]) / sum(srs[1:2])
  p2 <- sum(d$vardir[4:5]) / sum(srs[4:5])
  overall <- sum(d$vardir[ok]) / sum(srs[ok])
  expect_equal(e$deff, c(p1, p1, p1, p2, p2, overall))
  expect_equal(e$deff_source, c(rep("group", 5), "overall"))
  expect_equal(e$trials, d$n / e$deff)
})

test_that("effective_trials uses own design effects, then the fallback, then the pooled ones", {
  d <- et_data()
  e <- effective_trials("y", "vardir", "n", group = "prov", data = d, method = "area", deff_fallback = "fb")
  expect_equal(e$deff[c(1, 2, 4, 5)], c(1.5, 2, 1.2, 1.8))
  expect_equal(e$deff[3], 2.5)
  expect_equal(e$deff_source, c("area", "area", "fallback", "area", "area", "overall"))
  # vectors without data, floor at deff_min
  e2 <- effective_trials(d$y, d$vardir, d$n, method = "area", deff_min = 1.6)
  expect_equal(e2$deff[c(1, 2, 4, 5)], c(1.6, 2, 1.6, 1.8))
  expect_equal(e2$deff_source[c(3, 6)], c("group", "group"))
})

test_that("effective_trials handles missing areas and checks its input", {
  d <- et_data()
  d$y[2] <- NA
  e <- effective_trials(~y, ~vardir, ~n, data = d)
  expect_true(is.na(e$trials[2]) && is.na(e$deff[2]) && is.na(e$deff_source[2]))
  expect_false(anyNA(e$trials[-2]))
  expect_error(effective_trials(~y, ~vardir, ~n), "data")
  expect_error(effective_trials(d$y, d$vardir, d$n[-1]), "length")
  expect_error(effective_trials(c(0.1, 1.2), c(1, 1), c(5, 5)), "proportions")
  expect_error(effective_trials(c(0.1, 0.2), c(1, 1), c(5, 0)), "positive")
  expect_error(effective_trials(c(0.1, 0.2), c(1, -1), c(5, 5)), "negative")
  expect_error(effective_trials(c(0, 1), c(0, 0), c(5, 5)), "cannot be estimated")
  expect_error(effective_trials(c(0.1, 0.2), c(1, 1), c(5, 5), deff_min = 0), "deff_min")
})

test_that("effective_trials output feeds hb_area", {
  skip_if_not_installed("INLA")
  set.seed(4)
  d <- data.frame(area = 1:30, prov = rep(1:5, each = 6), n = rpois(30, 60), x = rnorm(30))
  d$y <- rbinom(30, d$n, plogis(-2.5 + 0.4 * d$x)) / d$n
  d$vardir <- 1.4 * d$y * (1 - d$y) / d$n
  d$trials <- effective_trials(~y, ~vardir, ~n, group = ~prov, data = d)$trials
  fit <- suppressMessages(hb_area(y ~ x, data = d, family = "binomial", trials = "trials", domain = "area"))
  expect_equal(fit$df_hb$trials, d$trials)
})
