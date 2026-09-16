library(testthat)
library(fastsae)

skip_if_not_installed("sae")

# ------------------------------------------------------------------
# Shared objects
# ------------------------------------------------------------------

library(dplyr)
mys_panel_nona <- mys_panel |>
  filter(!is.na(y))
mys_proxmat_nona <- mys_proxmat[-c(21, 25), -c(21, 25)]


# -----------------------------------------------------------------
# Test: compute_mse = FALSE (default)
# ------------------------------------------------------------------

test_that("eblup_stfh without MSE returns NA for mse and rse", {
  fit <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    compute_mse = FALSE,
    print_result = FALSE
  )

  expect_true("mse" %in% names(fit$df_eblup))
  expect_true("rse" %in% names(fit$df_eblup))
  expect_true(all(is.na(fit$df_eblup$mse)))
  expect_true(all(is.na(fit$df_eblup$rse)))
  expect_true(is.na(fit$B))
})

# -----------------------------------------------------------------
# Test: compute_mse = TRUE basic structure
# ------------------------------------------------------------------

test_that("eblup_stfh with MSE returns valid structure", {
  fit <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    compute_mse = TRUE,
    B = 50,
    seed = 42,
    print_result = FALSE
  )

  expect_true("mse" %in% names(fit$df_eblup))
  expect_true("rse" %in% names(fit$df_eblup))
  expect_false(any(is.na(fit$df_eblup$mse)))
  expect_false(any(is.na(fit$df_eblup$rse)))
  expect_equal(fit$B, 50)
})

# -----------------------------------------------------------------
# Test: MSE values are positive
# ------------------------------------------------------------------

test_that("MSE values are positive", {
  fit <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    compute_mse = TRUE,
    B = 50,
    seed = 42,
    print_result = FALSE
  )

  expect_true(all(fit$df_eblup$mse > 0))
  expect_true(all(fit$df_eblup$rse > 0))
})

# -----------------------------------------------------------------
# Test: MSE is reasonable compared to direct variance
# ------------------------------------------------------------------

test_that("MSE is in reasonable range compared to vardir", {
  fit <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    compute_mse = TRUE,
    B = 100,
    seed = 123,
    print_result = FALSE
  )

  # MSE should generally be less than or comparable to vardir
  # (SAE borrows strength, so MSE should be smaller than raw variance)
  median_mse <- median(fit$df_eblup$mse)
  median_vardir <- median(mys_panel_nona$vardir)

  expect_true(median_mse < median_vardir * 2,
              info = paste("median MSE", median_mse, "vs median vardir", median_vardir))
})

# -----------------------------------------------------------------
# Test: compare with sae::pbmseSTFH (informational)
# Note: MSE comparison with bootstrap methods has high variability.
# This test compares patterns and order, not exact values.
# ------------------------------------------------------------------

test_that("MSE pattern correlates with sae::pbmseSTFH", {
  skip("Long-running test - run manually when needed")
  # This test is informational only. Bootstrap MSE has high variance
  # and exact values will differ between implementations due to:
  # 1. Different random number generators/sequences
  # 2. Different re-fitting convergence paths
  # 3. Different bootstrap sample generation algorithms

  # Reference from sae package
  fit_sae <- sae::pbmseSTFH(
    mys_panel_nona$y ~ mys_panel_nona$x1 + mys_panel_nona$x2 + mys_panel_nona$x3,
    vardir = mys_panel_nona$vardir,
    D = length(unique(mys_panel_nona$area)),
    T = length(unique(mys_panel_nona$year)),
    proxmat = mys_proxmat_nona,
    B = 50
  )

  # Our implementation
  fit_fast <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    compute_mse = TRUE,
    B = 50,
    seed = 42,
    print_result = FALSE
  )

  # EBLUP should match perfectly
  expect_equal(
    fit_fast$df_eblup$eblup,
    as.numeric(fit_sae$est$eblup),
    tolerance = 1e-10
  )

  # MSE should have positive correlation (same ranking of areas)
  cor_mse <- cor(fit_fast$df_eblup$mse, fit_sae$mse, method = "spearman")
  expect_gt(cor_mse, 0.5,
             info = paste("Spearman correlation of MSE:", cor_mse))

  # Mean MSE should be in same order of magnitude
  ratio <- mean(fit_fast$df_eblup$mse) / mean(fit_sae$mse)
  expect_gt(ratio, 0.5,
           info = paste("Ratio of mean MSE:", ratio))
  expect_lt(ratio, 2.0,
           info = paste("Ratio of mean MSE:", ratio))
})

# -----------------------------------------------------------------
# Test: EBLUP estimates are unchanged when compute_mse = TRUE
# ------------------------------------------------------------------

test_that("EBLUP estimates are unchanged when compute_mse = TRUE", {
  fit_no_mse <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    compute_mse = FALSE,
    print_result = FALSE
  )

  fit_with_mse <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    compute_mse = TRUE,
    B = 50,
    seed = 42,
    print_result = FALSE
  )

  # EBLUP should be identical (MSE computation doesn't change estimates)
  expect_equal(
    fit_no_mse$df_eblup$eblup,
    fit_with_mse$df_eblup$eblup,
    tolerance = 1e-10
  )
})

# -----------------------------------------------------------------
# Test: parallel computation (n_threads > 1)
# ------------------------------------------------------------------

test_that("parallel computation gives same results as sequential", {
  fit_seq <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    compute_mse = TRUE,
    B = 50,
    n_threads = 1,
    seed = 42,
    print_result = FALSE
  )

  fit_par <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    compute_mse = TRUE,
    B = 50,
    n_threads = 2,
    seed = 42,
    print_result = FALSE
  )

  # EBLUP should be identical with same seed
  expect_equal(
    fit_seq$df_eblup$eblup,
    fit_par$df_eblup$eblup,
    tolerance = 1e-10
  )

  # MSE should be very close (within bootstrap variance)
  expect_equal(
    fit_seq$df_eblup$mse,
    fit_par$df_eblup$mse,
    tolerance = 1e-6
  )
})

# -----------------------------------------------------------------
# Test: spatial-only model (model = "S")
# ------------------------------------------------------------------

test_that("spatial-only model (model = 'S') computes MSE correctly", {
  fit <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "S",
    compute_mse = TRUE,
    B = 50,
    seed = 42,
    print_result = FALSE
  )

  expect_true("mse" %in% names(fit$df_eblup))
  expect_false(any(is.na(fit$df_eblup$mse)))
  expect_true(all(fit$df_eblup$mse > 0))
})

# -----------------------------------------------------------------
# Test: B parameter affects precision
# ------------------------------------------------------------------

test_that("larger B gives more stable MSE estimates", {
  fit_small_B <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    compute_mse = TRUE,
    B = 20,
    seed = 42,
    print_result = FALSE
  )

  fit_large_B <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    compute_mse = TRUE,
    B = 200,
    seed = 42,
    print_result = FALSE
  )

  # Both should have positive MSE
  expect_true(all(fit_small_B$df_eblup$mse > 0))
  expect_true(all(fit_large_B$df_eblup$mse > 0))

  # Variance of MSE estimates should be lower with larger B
  var_small_B <- var(fit_small_B$df_eblup$mse)
  var_large_B <- var(fit_large_B$df_eblup$mse)

  expect_true(var_large_B < var_small_B * 1.5,
              info = paste("Variance small B:", var_small_B, ", large B:", var_large_B))
})
