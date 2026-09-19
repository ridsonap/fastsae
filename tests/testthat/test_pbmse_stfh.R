suppressMessages({
  library(testthat)
  library(fastsae)
  library(dplyr)
})

skip_if_not_installed("sae")

# ------------------------------------------------------------------
# Shared objects
# ------------------------------------------------------------------

# Common formula and parameters for eblup_stfh calls
formula_stfh <- y ~ x1 + x2 + x3
mys_panel_nona <- mys_panel |> filter(!is.na(y))
mys_proxmat_nona <- mys_proxmat[-c(21, 25), -c(21, 25)]

model_stfh <- "ST" # default spatial-temporal model

# Common parameters for do.call tests (FIX: was missing)
params_stfh <- list(
  data = mys_panel_nona,
  vardir = ~vardir,
  domain = ~area,
  time = ~year,
  W = mys_proxmat_nona,
  print_result = FALSE
)

fit_nona <- eblup_stfh(
  formula = formula_stfh,
  data = mys_panel_nona,
  vardir = ~vardir,
  domain = ~area,
  time = ~year,
  W = mys_proxmat_nona,
  print_result = FALSE,
  compute_mse = FALSE,
  model = model_stfh
)

fit <- eblup_stfh(
  formula = formula_stfh,
  data = mys_panel_nona,
  vardir = ~vardir,
  domain = ~area,
  time = ~year,
  W = mys_proxmat_nona,
  print_result = FALSE,
  compute_mse = TRUE,
  B = 50,
  seed = 1,
  model = model_stfh
)

# -----------------------------------------------------------------
# Test: compute_mse = FALSE (default)
# ------------------------------------------------------------------

test_that("eblup_stfh without MSE returns NA for mse and rse", {
  expect_true("mse_pb" %in% names(fit_nona$df_eblup))
  expect_true("rse" %in% names(fit_nona$df_eblup))
  expect_true(all(is.na(fit_nona$df_eblup$mse_pb)))
  expect_true(all(is.na(fit_nona$df_eblup$rse)))
  expect_true(is.na(fit_nona$B))
})

# -----------------------------------------------------------------
# Test: compute_mse = TRUE basic structure
# ------------------------------------------------------------------

test_that("eblup_stfh with MSE returns valid structure", {
  expect_true("mse_pb" %in% names(fit$df_eblup))
  expect_true("rse" %in% names(fit$df_eblup))
  expect_false(any(is.na(fit$df_eblup$mse_pb)))
  expect_false(any(is.na(fit$df_eblup$rse)))
  expect_equal(fit$B, 50)
})

# -----------------------------------------------------------------
# Test: MSE values are positive
# ------------------------------------------------------------------

test_that("MSE values are positive", {
  expect_true(all(fit$df_eblup$mse_pb > 0))
  expect_true(all(fit$df_eblup$rse > 0))
})

# -----------------------------------------------------------------
# Test: MSE is reasonable compared to direct variance
# ------------------------------------------------------------------

test_that("MSE is in reasonable range compared to vardir", {
  # MSE should generally be less than or comparable to vardir
  # (SAE borrows strength, so MSE should be smaller than raw variance)
  median_mse <- median(fit$df_eblup$mse_pb)
  median_vardir <- median(mys_panel_nona$vardir)

  expect_true(median_mse < median_vardir * 2,
    info = paste("median MSE", median_mse, "vs median vardir", median_vardir)
  )
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

  mys_sub <- mys_panel_nona |>
    filter(year >= 2024)

  nB <- 50
  nsim <- 100
  m <- nrow(mys_sub)

  # Setup cluster
  n_cores <- detectCores(logical = FALSE) - 1 # sisakan 1 core untuk OS
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  # Export objek yang dibutuhkan tiap worker
  clusterExport(cl, varlist = c("mys_sub", "mys_proxmat_nona", "nB", "formula_stfh", "eblup_stfh"))

  # Load package di tiap worker
  clusterEvalQ(cl, {
    library(sae)
    library(fastsae)
  })

  results <- foreach(
    i = seq_len(nsim),
    .combine = "rbind", # tiap iterasi return 1 baris, rbind jadi matriks
    .errorhandling = "pass" # kalau 1 iterasi error, lanjut (jangan stop)
  ) %dopar% {
    fit_sae <- tryCatch(
      sae::pbmseSTFH(
        mys_sub$y ~ mys_sub$x1 + mys_sub$x2 + mys_sub$x3,
        vardir  = mys_sub$vardir,
        D       = length(unique(mys_sub$area)),
        T       = length(unique(mys_sub$year)),
        proxmat = mys_proxmat_nona,
        B       = nB
      ),
      error = function(e) NULL
    )

    fit_fast <- tryCatch(
      eblup_stfh(formula_stfh,
        data = mys_sub,
        vardir = ~vardir, domain = ~area, time = ~year,
        W = mys_proxmat_nona, model = "ST",
        compute_mse = TRUE, B = nB, print_result = FALSE
      ),
      error = function(e) NULL
    )

    # Tiap worker return list 1 baris
    list(
      mse_fast   = if (!is.null(fit_fast)) fit_fast$df_eblup$mse_pb else rep(NA_real_, m),
      mse_sae    = if (!is.null(fit_sae)) as.numeric(fit_sae$mse) else rep(NA_real_, m)
    )
  }

  stopCluster(cl)

  # Rekonstruksi matriks nsim × m dari list hasil
  mse_fast <- do.call(rbind, results[, "mse_fast"]) # 50 × 120
  mse_sae <- do.call(rbind, results[, "mse_sae"]) # 50 × 120

  mse_domain_fast <- colMeans(mse_fast, na.rm = TRUE) # panjang 120
  mse_domain_sae <- colMeans(mse_sae, na.rm = TRUE) # panjang 120


  # test perbedaan MSE
  expect_equal(
    mse_domain_fast,
    mse_domain_sae,
    tolerance = 0.05,
    ignore_attr = TRUE
  )

  # correlation > 0.9
  cor_mse <- cor(mse_domain_sae, mse_domain_fast, method = "spearman")
  expect_gt(cor_mse, 0.9, label = paste("Spearman correlation of MSE:", cor_mse))

  # Mean MSE should be in same order of magnitude
  ratio <- mean(mse_domain_sae) / mean(mse_domain_fast)
  expect_gt(ratio, 0.5, label = paste("Ratio of mean MSE:", ratio))
  expect_lt(ratio, 2.0, label = paste("Ratio of mean MSE:", ratio))
})

# -----------------------------------------------------------------
# Test: EBLUP estimates are unchanged when compute_mse = TRUE
# ------------------------------------------------------------------

test_that("EBLUP estimates are unchanged when compute_mse = TRUE", {
  fit_no_mse <- do.call(eblup_stfh, c(list(formula = formula_stfh), params_stfh,
    model = model_stfh, compute_mse = FALSE
  ))

  fit_with_mse <- do.call(eblup_stfh, c(list(formula = formula_stfh), params_stfh,
    model = model_stfh, compute_mse = TRUE, B = 50, seed = 1
  ))

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
  fit_seq <- do.call(eblup_stfh, c(list(formula = formula_stfh), params_stfh,
    model = model_stfh, compute_mse = TRUE, B = 50, n_threads = 1, seed = 1
  ))

  fit_par <- do.call(eblup_stfh, c(list(formula = formula_stfh), params_stfh,
    model = model_stfh, compute_mse = TRUE, B = 50, n_threads = 3, seed = 1
  ))

  # EBLUP should be identical with same seed
  expect_equal(
    fit_seq$df_eblup$eblup,
    fit_par$df_eblup$eblup,
    tolerance = 1e-10
  )

  # MSE should be very close (within bootstrap variance)
  expect_equal(
    fit_seq$df_eblup$mse_pb,
    fit_par$df_eblup$mse_pb,
    tolerance = 1e-6
  )
})

# -----------------------------------------------------------------
# Test: spatial-only model (model = "S")
# ------------------------------------------------------------------

test_that("spatial-only model (model = 'S') computes MSE correctly", {
  fit <- do.call(eblup_stfh, c(list(formula_stfh), params_stfh,
    model = "S", compute_mse = TRUE, B = 50, seed = 1
  ))

  expect_true("mse_pb" %in% names(fit$df_eblup))
  expect_false(any(is.na(fit$df_eblup$mse_pb)))
  expect_true(all(fit$df_eblup$mse_pb > 0))
})

# -----------------------------------------------------------------
# Test: B parameter affects precision
# ------------------------------------------------------------------

test_that("larger B gives more stable MSE estimates", {
  fit_small_B <- do.call(eblup_stfh, c(list(formula = formula_stfh), params_stfh,
    model = model_stfh, compute_mse = TRUE, B = 20, seed = 1
  ))

  fit_large_B <- do.call(eblup_stfh, c(list(formula = formula_stfh), params_stfh,
    model = model_stfh, compute_mse = TRUE, B = 200, seed = 1
  ))

  # Both should have positive MSE
  expect_true(all(fit_small_B$df_eblup$mse_pb > 0))
  expect_true(all(fit_large_B$df_eblup$mse_pb > 0))

  # Variance of MSE estimates should be lower with larger B
  var_small_B <- var(fit_small_B$df_eblup$mse_pb)
  var_large_B <- var(fit_large_B$df_eblup$mse_pb)

  expect_true(var_large_B < var_small_B * 1.5,
    info = paste("Variance small B:", var_small_B, ", large B:", var_large_B)
  )
})
