library(fastsae)

test_that("eblup_stfh computes fast analytical MSE", {
  library(dplyr)
  mys_panel_nona <- mys_panel |>
    filter(!is.na(y) & year >= 2024)

  # Measure time for analytical MSE
  t0 <- Sys.time()
  fit_anal <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat[-c(21, 25), -c(21, 25)],
    model = "ST",
    compute_mse = TRUE,
    mse_type = "analytical",
    print_result = FALSE
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Must be fast (< 1.0 second)
  expect_lt(elapsed, 1.5)
  expect_true(fit_anal$convergence)
  expect_equal(fit_anal$mse_type, "analytical")

  # MSE and RSE checks
  expect_true("mse" %in% names(fit_anal$df_eblup))
  expect_true("rse" %in% names(fit_anal$df_eblup))
  expect_true(all(fit_anal$df_eblup$mse > 0))
  expect_true(all(fit_anal$df_eblup$rse > 0))
  expect_false(any(is.na(fit_anal$df_eblup$mse)))

  # MSE components checks
  expect_true(!is.null(fit_anal$g1))
  expect_true(!is.null(fit_anal$g2))
  expect_true(!is.null(fit_anal$g3))
  expect_true(all(fit_anal$g1 >= 0))
  expect_true(all(fit_anal$g2 >= 0))
  expect_true(all(fit_anal$g3 >= 0))

  # MSE is reasonable compared to direct sampling variance
  expect_true(mean(fit_anal$df_eblup$mse) < mean(fit_anal$df_eblup$vardir) * 1.5)
})

test_that("eblup_stfh analytical MSE works for Spatial-only model S", {
  mys_panel_nona <- mys_panel |>
    filter(!is.na(y))

  fit_s_anal <- eblup_stfh(
    y ~ x1 + x2,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat[-c(21, 25), -c(21, 25)],
    model = "S",
    compute_mse = TRUE,
    mse_type = "analytical",
    print_result = FALSE
  )

  expect_true(fit_s_anal$convergence)
  expect_equal(fit_s_anal$model, "S")
  expect_true(all(fit_s_anal$df_eblup$mse > 0))
  expect_true(all(fit_s_anal$df_eblup$rse > 0))
})
