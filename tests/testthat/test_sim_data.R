test_that("sim_spatial_weights produces valid matrices across types and styles", {
  # 1. KNN with Binary style
  W_knn <- sim_spatial_weights(D = 30, type = "knn", k = 4, style = "B", seed = 123)
  expect_equal(dim(W_knn), c(30, 30))
  expect_equal(unname(diag(W_knn)), rep(0, 30))
  expect_true(isSymmetric(W_knn))
  expect_true(all(W_knn %in% c(0, 1)))
  expect_equal(attr(W_knn, "type"), "knn")
  expect_equal(dim(attr(W_knn, "coords")), c(30, 2))

  # 2. KNN with Row-standardized style
  W_knn_w <- sim_spatial_weights(D = 25, type = "knn", k = 3, style = "W", seed = 456)
  expect_equal(dim(W_knn_w), c(25, 25))
  expect_equal(unname(rowSums(W_knn_w)), rep(1, 25), tolerance = 1e-6)

  # 3. Grid Lattice
  W_grid <- sim_spatial_weights(D = 36, type = "grid", style = "B")
  expect_equal(dim(W_grid), c(36, 36))
  expect_true(isSymmetric(W_grid))
  expect_true(all(rowSums(W_grid) >= 2)) # in a grid, corners have at least 2 neighbors

  # 4. Ring Lattice
  W_ring <- sim_spatial_weights(D = 20, type = "ring", style = "B")
  expect_equal(dim(W_ring), c(20, 20))
  expect_true(isSymmetric(W_ring))
  expect_true(all(rowSums(W_ring) == 2)) # exactly 2 neighbors in 1D ring

  # 5. Reproducibility
  W1 <- sim_spatial_weights(D = 15, seed = 999)
  W2 <- sim_spatial_weights(D = 15, seed = 999)
  expect_identical(W1, W2)

  # 6. Error handling
  expect_error(sim_spatial_weights(D = 1), "D.*must be an integer >= 2")
})

test_that("sim_area_data generates complete multi-distribution dataset", {
  sim <- sim_area_data(D = 40, spatial_type = "knn", n_unsampled = 5, seed = 101)

  expect_s3_class(sim, "fastsae_sim_data")
  expect_named(sim, c("data", "W", "W_std", "coords", "u_spatial", "u_iid", "u_total", "phi", "rho"))

  df <- sim$data
  expect_equal(nrow(df), 40)
  expect_true(all(c("domain", "x1", "x2", "y_gaussian", "vardir", "y_poisson",
                    "exposure", "y_binomial", "trials", "y_beta", "y_nbinomial",
                    "y_gamma", "x_coord", "y_coord") %in% names(df)))

  # Check unsampled counts
  expect_equal(sum(is.na(df$y_gaussian)), 5)
  expect_equal(sum(is.na(df$y_poisson)), 5)
  expect_equal(sum(is.na(df$y_binomial)), 5)
  expect_equal(sum(is.na(df$y_beta)), 5)

  # Check auxiliary values are not NA
  expect_false(any(is.na(df$vardir)))
  expect_false(any(is.na(df$exposure)))
  expect_false(any(is.na(df$trials)))

  # Distribution constraints
  sampled_beta <- stats::na.omit(df$y_beta)
  expect_true(all(sampled_beta > 0 & sampled_beta < 1))

  sampled_pois <- stats::na.omit(df$y_poisson)
  expect_true(all(sampled_pois >= 0))

  sampled_bin <- stats::na.omit(df$y_binomial)
  sampled_trials <- df$trials[!is.na(df$y_binomial)]
  expect_true(all(sampled_bin >= 0 & sampled_bin <= sampled_trials))

  # Print method runs cleanly
  expect_output(print(sim), "fastsae Simulated Multi-Distribution Area Data")
})

test_that("sim_area_data works with pre-specified spatial weights W", {
  data("mys_proxmat", package = "fastsae")
  sim_custom <- sim_area_data(W = mys_proxmat, n_unsampled = 6, seed = 202)

  expect_equal(nrow(sim_custom$data), nrow(mys_proxmat))
  expect_equal(dim(sim_custom$W), dim(mys_proxmat))
  expect_equal(sum(is.na(sim_custom$data$y_gaussian)), 6)
})

test_that("ebp_area fits successfully on simulated data across families", {
  skip_if_not(requireNamespace("INLA", quietly = TRUE), "INLA not installed")

  sim <- sim_area_data(D = 35, spatial_type = "knn", n_unsampled = 5, seed = 789)
  df <- sim$data
  W <- sim$W

  # 1. Gaussian Fay-Herriot with simulated data
  fit_gauss <- ebp_area(y_gaussian ~ x1 + x2, data = df, vardir = "vardir")
  expect_s3_class(fit_gauss, "fastsae")
  expect_equal(nrow(fit_gauss$df_ebp), 35)

  # 2. Spatial Poisson with BYM2 and simulated W
  fit_pois <- ebp_area(y_poisson ~ x1 + x2, data = df, exposure = "exposure",
                       family = "poisson", spatial = "bym2", W = W)
  expect_s3_class(fit_pois, "fastsae")
  expect_equal(nrow(fit_pois$df_ebp), 35)
  expect_false(any(is.na(fit_pois$df_ebp$ebp)))

  # 3. Spatial Binomial with BYM2 and simulated W
  fit_bin <- ebp_area(y_binomial ~ x1 + x2, data = df, trials = "trials",
                      family = "binomial", spatial = "bym2", W = W)
  expect_s3_class(fit_bin, "fastsae")
  expect_equal(nrow(fit_bin$df_ebp), 35)

  # 4. Beta regression
  fit_beta <- ebp_area(y_beta ~ x1 + x2, data = df, family = "beta")
  expect_s3_class(fit_beta, "fastsae")
  expect_equal(nrow(fit_beta$df_ebp), 35)
})

test_that("built-in sae_area_multi dataset loads and integrates with mys_proxmat", {
  data("sae_area_multi", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  expect_equal(nrow(sae_area_multi), 42)
  expect_equal(sae_area_multi$area, 1:42)
  expect_true("y_gaussian" %in% names(sae_area_multi))
  expect_true("y_poisson" %in% names(sae_area_multi))
  expect_true("y_binomial" %in% names(sae_area_multi))
  expect_true("y_beta" %in% names(sae_area_multi))

  skip_if_not(requireNamespace("INLA", quietly = TRUE), "INLA not installed")
  # Test spatial fit with mys_proxmat
  fit_multi_spatial <- ebp_area(y_poisson ~ x1 + x2, data = sae_area_multi,
                                exposure = "exposure", family = "poisson",
                                spatial = "bym2", W = mys_proxmat)
  expect_s3_class(fit_multi_spatial, "fastsae")
  expect_equal(nrow(fit_multi_spatial$df_ebp), 42)
})
