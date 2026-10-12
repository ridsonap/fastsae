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

  # The benchmark holds exactly on the original scale (ratio benchmarking of
  # the back-transformed estimates; a WFQ constraint would only hold on the
  # log scale)
  w <- mys_pos$n / sum(mys_pos$n)
  direct_bench <- sum(w * mys_pos$y)
  est_bench <- sum(w * fit_log_bench$df_eblup$eblup)
  expect_equal(est_bench, direct_bench, tolerance = 1e-10)
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

test_that("Log-FH unsampled areas use the synthetic correction exp(x'b + s2/2) even if vardir is given", {
  data(mys, package = "fastsae")
  d <- mys[!is.na(mys$y), ]
  d$y[1:3] <- NA # unsampled, but keep their (original-scale) vardir
  expect_false(anyNA(d$vardir[1:3]))
  fit <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = d, transform = "log", print_result = FALSE)
  X <- model.matrix(~ x1 + x2, d)
  xb <- as.vector(X[1:3, ] %*% fit$estcoef$beta)
  expect_equal(fit$df_eblup$eblup[1:3], exp(xb + 0.5 * fit$random_effect_var), tolerance = 1e-10)
})
