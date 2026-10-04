library(fastsae)

test_that("eblup_fh supports transform = 'log' with Slud-Maiti bias correction", {
  # Subset to sampled positive areas
  mys_pos <- mys[!is.na(mys$y) & mys$y > 0, ]

  # Standard linear FH
  fit_linear <- eblup_fh(
    y ~ x1 + x2 + x3,
    data = mys_pos,
    vardir = "vardir"
  )

  # Log FH with Slud-Maiti bias correction
  fit_log <- eblup_fh(
    y ~ x1 + x2 + x3,
    data = mys_pos,
    vardir = "vardir",
    transform = "log"
  )

  expect_s3_class(fit_log, "fastsae")
  expect_equal(fit_log$model, "Log-FH")
  expect_equal(fit_log$transform, "log")
  expect_true(fit_log$convergence)

  # EBLUPs and MSE on original scale
  expect_true(all(fit_log$df_eblup$eblup > 0))
  expect_true(all(fit_log$df_eblup$mse > 0))
  expect_true(all(fit_log$df_eblup$rse > 0))

  # Log scale columns preserved
  expect_true("eblup_log" %in% names(fit_log$df_eblup))
  expect_true("mse_log" %in% names(fit_log$df_eblup))

  # Slud-Maiti bias correction: eblup >= exp(eblup_log)
  expect_true(all(fit_log$df_eblup$eblup >= exp(fit_log$df_eblup$eblup_log) - 1e-10))

  # Response matches original y
  expect_equal(fit_log$df_eblup$y, mys_pos$y)
  expect_equal(fit_log$df_eblup$vardir, mys_pos$vardir)
})

test_that("Log-FH integrates with self-benchmarking", {
  mys_pos <- mys[!is.na(mys$y) & mys$y > 0, ]

  fit_log_bench <- eblup_fh(
    y ~ x1 + x2 + x3,
    data = mys_pos,
    vardir = "vardir",
    transform = "log",
    self_benchmark = TRUE,
    benchmark_weight = "n"
  )

  expect_true(fit_log_bench$self_benchmark)
  expect_s3_class(fit_log_bench$benchmark_summary, "data.frame")

  # Weighted benchmark target on log scale is exact (Wang, Fuller, and Qu, 2008)
  w <- mys_pos$n / sum(mys_pos$n)
  direct_log_bench <- sum(w * log(mys_pos$y))
  est_log_bench <- sum(w * fit_log_bench$df_eblup$eblup_log)
  expect_equal(est_log_bench, direct_log_bench, tolerance = 1e-4)

  # On original scale, it is close to direct benchmark target
  direct_bench <- sum(w * mys_pos$y)
  est_bench <- sum(w * fit_log_bench$df_eblup$eblup)
  expect_equal(est_bench, direct_bench, tolerance = 0.1)
})

test_that("Log-FH validation rejects non-positive response values", {
  df_neg <- mys[!is.na(mys$y), ]
  df_neg$y[1] <- -1.0

  expect_error(
    eblup_fh(
      y ~ x1 + x2 + x3,
      data = df_neg,
      vardir = "vardir",
      transform = "log"
    ),
    "strictly positive"
  )
})
