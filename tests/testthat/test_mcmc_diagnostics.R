# MCMC convergence diagnostics of jax_engine() fits: the R-hat and ESS of the
# NumPyro backend against the R package posterior, the summary, warnings and
# plots on the R side (fake backend), and an end-to-end fit.

skip_if_no_numpyro_diag <- function() {
  skip_on_cran()
  skip_if_not_installed("reticulate")
  if (!suppressMessages(check_numpyro_available())) skip("NumPyro / JAX not available")
}

# Backend output with the diagnostics of 2 chains x 50 draws of three parameters
fake_diag_fit <- function(rhat = c(1.002, 1.001, 1.003), ess = c(400, 450, 500), divergent = 0L) {
  set.seed(5)
  list(
    rhat_max = max(rhat),
    diagnostics = list(
      parameter = list("sigma2_u", "beta[2]", "hb[1]"), mean = c(0.3, 0.5, 1), sd = c(0.1, 0.1, 0.2),
      rhat = rhat, ess_bulk = ess, ess_tail = ess + 10
    ),
    trace = list(names = list("sigma2_u", "beta[2]", "hb[1]"), draws = array(rnorm(2 * 50 * 3), c(2, 50, 3))),
    num_divergent = divergent,
    treedepth_hits = 0L
  )
}

test_that("rank-normalized R-hat and bulk/tail ESS match the posterior package", {
  skip_if_no_numpyro_diag()
  skip_if_not_installed("posterior")
  backend <- suppressMessages(.get_numpyro_backend(device = "cpu"))
  set.seed(1)
  C <- 4
  N <- 251 # odd: split chains drop the middle draw
  cases <- list(
    good = matrix(rnorm(C * N), N, C),
    autocorrelated = apply(matrix(rnorm(C * N), N, C), 2, function(e) as.numeric(stats::filter(e, 0.9, "recursive"))),
    shifted = sweep(matrix(rnorm(C * N), N, C), 2, c(0, 0, 0, 1.5), "+"),
    heavy_tails = matrix(rt(C * N, df = 1.5), N, C),
    scale = sweep(matrix(rnorm(C * N), N, C), 2, c(1, 1, 1, 4), "*")
  )
  arr <- array(NA_real_, c(C, N, length(cases)))
  for (k in seq_along(cases)) arr[, , k] <- t(cases[[k]])
  np <- reticulate::import("numpy", convert = FALSE)
  x <- np$asarray(arr)
  py <- cbind(
    rhat = reticulate::py_to_r(backend$rank_rhat(x)),
    bulk = reticulate::py_to_r(backend$ess_bulk(x)),
    tail = reticulate::py_to_r(backend$ess_tail(x))
  )
  ref <- t(vapply(cases, function(m) c(posterior::rhat(m), posterior::ess_bulk(m), posterior::ess_tail(m)), numeric(3)))
  expect_equal(unname(py), unname(ref), tolerance = 1e-10)
  # the shifted chain is caught by R-hat, the inflated scale by the tail ESS
  expect_gt(py[3, "rhat"], 1.1)
  expect_lt(py[5, "tail"], 100)
})

test_that(".mcmc_summary names the coefficients and orders the trace by iteration and chain", {
  mc <- .mcmc_summary(fake_diag_fit(), c("(Intercept)", "x1"), chains = 2)
  expect_equal(mc$table$parameter, c("sigma2_u", "x1", "hb[1]"))
  expect_equal(mc$ess_min, 400)
  expect_equal(dim(mc$trace), c(50, 2, 3))
  expect_equal(dimnames(mc$trace)$parameter, c("sigma2_u", "x1", "hb[1]"))
  draws <- fake_diag_fit()$trace$draws
  expect_equal(mc$trace[, 2, "x1"], draws[2, , 2])
  # a backend without diagnostics leaves them NULL
  expect_null(.mcmc_summary(list(), "x1", 2)$table)
})

test_that("jax fits warn about R-hat, low ESS and divergent transitions", {
  ok <- .mcmc_summary(fake_diag_fit(), c("(Intercept)", "x1"), 2)
  expect_no_warning(expect_no_message(.warn_mcmc(1.003, ok, chains = 2)))
  expect_message(.warn_mcmc(1.02, ok, chains = 2), "1.01")
  expect_warning(.warn_mcmc(1.08, ok, chains = 2), "not converged.*sigma2_u|not converged")
  low <- .mcmc_summary(fake_diag_fit(ess = c(150, 450, 500)), c("(Intercept)", "x1"), 2)
  expect_warning(.warn_mcmc(1.003, low, chains = 2), "ESS = 150 \\(sigma2_u\\)")
  # 100 per chain: 150 is enough for one chain
  expect_no_warning(.warn_mcmc(1.003, low, chains = 1))
  div <- .mcmc_summary(fake_diag_fit(divergent = 3L), c("(Intercept)", "x1"), 2)
  expect_warning(.warn_mcmc(1.003, div, chains = 2), "3 divergent transitions")
})

test_that("autoplot() draws trace and rank plots of a jax fit", {
  mc <- .mcmc_summary(fake_diag_fit(), c("(Intercept)", "x1"), 2)
  fit <- structure(c(list(engine = "jax", draws = matrix(rnorm(4 * 100), 4, 100)), .mcmc_fields(list(mcmc = mc))),
    class = c("fastsae_hb_area", "fastsae_hb", "fastsae")
  )
  p <- autoplot(fit, type = "trace")
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 50 * 2 * 3)
  expect_match(levels(p$data$parameter)[1], "sigma2_u\nR-hat 1.002, ESS 400", fixed = TRUE)

  # hb[j] outside the stored trace comes from the chain-major draws
  p <- autoplot(fit, type = "trace", pars = c("x1", "hb[4]"))
  hb4 <- p$data[startsWith(as.character(p$data$parameter), "hb[4]"), ]
  expect_equal(hb4$value[hb4$chain == "2"], fit$draws[4, 51:100])

  r <- autoplot(fit, type = "rank", bins = 10)
  expect_s3_class(r, "ggplot")
  # each chain's counts over the bins add up to its draws (the last bin is repeated
  # at the right edge for the step)
  cnt <- r$data[r$data$bin < 10, ]
  expect_equal(unname(tapply(cnt$count, list(cnt$parameter, cnt$chain), sum)), matrix(50, 3, 2))

  expect_error(autoplot(fit, type = "trace", pars = "nope"), "Unknown parameter")
  inla_fit <- structure(list(engine = "inla"), class = c("fastsae_hb_area", "fastsae_hb", "fastsae"))
  expect_error(autoplot(inla_fit, type = "rank"), "jax_engine")
})

test_that("print() reports the ESS and divergent transitions of a jax fit", {
  mc <- .mcmc_summary(fake_diag_fit(divergent = 2L), c("(Intercept)", "x1"), 2)
  fit <- c(list(convergence = TRUE, rhat_max = 1.003, engine = "jax"), .mcmc_fields(list(mcmc = mc)))
  expect_equal(.convergence_detail(fit), "max R-hat = 1.003, min ESS = 400, 2 divergent")
})

test_that("jax fits diagnose every parameter, the variances included", {
  skip_if_no_numpyro_diag()
  data("mys", package = "fastsae")
  # Without vardir the area and residual variances of an area-level model are weakly
  # identified; the R-hat used to cover only the coefficients and the area means
  eng <- jax_engine(warmup = 300, samples = 300, chains = 2, device = "cpu", progress = FALSE)
  fit <- suppressWarnings(suppressMessages(
    hb_area(y ~ x1 + x2 + x3, data = mys, engine = eng, print_result = FALSE)
  ))
  tab <- fit$mcmc_diagnostics
  expect_true(all(c("sigma2_u", "sigma2_e", "(Intercept)", "x3", "hb[1]") %in% tab$parameter))
  expect_equal(sum(startsWith(tab$parameter, "hb[")), nrow(mys))
  expect_equal(fit$rhat_max, max(tab$rhat, na.rm = TRUE))
  expect_equal(fit$ess_min, min(tab$ess_bulk, tab$ess_tail, na.rm = TRUE))
  expect_true(is.numeric(fit$divergences) && is.numeric(fit$treedepth_hits))
  expect_equal(dim(fit$trace)[1:2], c(300, 2))
  expect_true(all(c("sigma2_u", "sigma2_e") %in% dimnames(fit$trace)$parameter))
})
