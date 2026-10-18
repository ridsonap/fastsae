test_that("spatial_test works with standard data and returns expected structure", {
  data(mys)
  data(mys_proxmat)

  st <- spatial_test(
    formula = y ~ x1 + x2,
    vardir = "vardir",
    data = mys,
    W = mys_proxmat
  )

  expect_s3_class(st, "fastsae_spatial_test")
  expect_true("response_test" %in% names(st))
  expect_true("residual_test" %in% names(st))
  expect_true("model_comparison" %in% names(st))
  expect_true("spatial_recommended" %in% names(st))
  expect_type(st$spatial_recommended, "logical")
  expect_type(st$reason, "character")

  # 1. Direct response Moran's I
  expect_true(!is.na(st$response_test$I))
  expect_true(!is.na(st$response_test$p_value))
  expect_equal(st$response_test$expected, -1 / (st$n - 1))

  # 2. Residual Moran's I and LM tests
  expect_true(!is.na(st$residual_test$moran$I))
  expect_true(!is.na(st$residual_test$lm_error$statistic))
  expect_true(!is.na(st$residual_test$lm_error$p_value))
  expect_true(!is.na(st$residual_test$lm_lag$statistic))

  # 3. Model comparison
  expect_true(!is.null(st$model_comparison))
  expect_true(!is.na(st$model_comparison$lrt_stat))
  expect_true(!is.na(st$model_comparison$delta_aic))
  expect_true(!is.na(st$model_comparison$rho))
  expect_identical(st$model_comparison$method, "ML")
  # the likelihood ratio uses the ML fits of FH and SFH
  fh_ml <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys, method = "ML", print_result = FALSE)
  sfh_ml <- eblup_sfh(y ~ x1 + x2,
    vardir = "vardir", data = mys, W = mys_proxmat, method = "ML",
    print_result = FALSE
  )
  expect_equal(st$model_comparison$ll_fh, unname(fh_ml$goodness["loglikelihood"]))
  expect_equal(
    st$model_comparison$lrt_stat,
    max(0, 2 * unname(sfh_ml$goodness["loglikelihood"] - fh_ml$goodness["loglikelihood"]))
  )
  # one decision test: the LRT when vardir is given
  expect_identical(st$decision_test, "lrt")
  expect_identical(st$spatial_recommended, st$model_comparison$lrt_pval < st$alpha)

  # 4. Print method
  expect_message(print(st), "Spatial")
  expect_message(print(st), "Direct Response Variable")
})

test_that("spatial_test handles data without vardir (ESDA only)", {
  data(mys)
  data(mys_proxmat)

  # Run spatial test without vardir
  st_no_var <- spatial_test(
    formula = y ~ x1 + x2,
    data = mys,
    W = mys_proxmat
  )

  expect_s3_class(st_no_var, "fastsae_spatial_test")
  expect_true(!is.na(st_no_var$response_test$I))
  expect_true(!is.na(st_no_var$residual_test$moran$I))
  expect_null(st_no_var$model_comparison)
  expect_identical(st_no_var$decision_test, "lm_error")
  expect_identical(st_no_var$spatial_recommended, st_no_var$residual_test$lm_error$p_value < st_no_var$alpha)
  expect_message(print(st_no_var), "Spatial")
})

test_that("spatial_test auto-constructs W when data is sf and W is NULL", {
  skip_if_not_installed("sf")

  # Create small 9-polygon grid
  grid_geom <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0, 0, 3, 0, 3, 3, 0, 3, 0, 0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1)
  )
  set.seed(42)
  df_sf <- sf::st_sf(
    domain = 1:9,
    y = rnorm(9, 10, 2),
    x1 = runif(9, 1, 5),
    vardir = rep(0.5, 9),
    geometry = grid_geom
  )

  # Call spatial_test with W = NULL
  st_sf <- spatial_test(
    formula = y ~ x1,
    vardir = "vardir",
    data = df_sf
  )

  expect_s3_class(st_sf, "fastsae_spatial_test")
  expect_true(!is.null(st_sf$W))
  expect_message(print(st_sf), "Spatial")
})

test_that("spatial_test handles collinear covariates without crashing in LM-lag", {
  set.seed(42)
  n <- 25
  x1 <- runif(n, 1, 5)
  x2 <- runif(n, 2, 8)
  x3 <- x1 * 2 # exactly collinear
  y <- 3 + 1.2 * x1 - 0.5 * x2 + rnorm(n, 0, 0.5)
  coords <- matrix(runif(n * 2), ncol = 2)
  W <- create_weights(coords, method = "knn", k = 3, style = "W")
  df <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)

  # Should not error or crash with 'missing value where TRUE/FALSE needed'
  expect_no_error({
    st_collin <- spatial_test(y ~ x1 + x2 + x3, data = df, W = W)
  })
  expect_true(is.finite(st_collin$residual_test$lm_lag$statistic))
  expect_true(is.finite(st_collin$residual_test$lm_lag$p_value))
})

test_that("spatial_test throws informative error when auxiliary variables contain NA", {
  set.seed(42)
  n <- 20
  x1 <- runif(n)
  x1[3] <- NA # NA in sampled area
  y <- rnorm(n)
  coords <- matrix(runif(n * 2), ncol = 2)
  W <- create_weights(coords, method = "knn", k = 3, style = "W")
  df <- data.frame(y = y, x1 = x1)

  expect_error(
    spatial_test(y ~ x1, data = df, W = W),
    "Auxiliary variables in .* contain NA values"
  )
})


test_that("spatial_test agrees with spdep and keeps its size under no spatial effect", {
  skip_if_not_installed("spdep")
  set.seed(5)
  D <- 80
  coords <- cbind(runif(D), runif(D))
  nb <- spdep::make.sym.nb(spdep::knn2nb(spdep::knearneigh(coords, k = 4)))
  W <- spdep::nb2mat(nb, style = "W")
  lw <- spdep::nb2listw(nb, style = "W")

  x1 <- rnorm(D)
  psi <- runif(D, 0.3, 1.5)
  dat <- data.frame(y = 1 + x1 + rnorm(D, sd = 0.8) + rnorm(D, sd = sqrt(psi)), x1 = x1, psi = psi)
  st <- spatial_test(y ~ x1, data = dat, W = W)
  ols <- lm(y ~ x1, data = dat)
  expect_equal(st$residual_test$moran$I, unname(spdep::lm.morantest(ols, lw)$estimate[1]))
  expect_equal(
    st$residual_test$lm_error$statistic,
    unname(spdep::lm.RStests(ols, lw, test = "LMerr")[[1]]$statistic)
  )
  expect_equal(
    st$residual_test$lm_lag$statistic,
    unname(spdep::lm.RStests(ols, lw, test = "LMlag")[[1]]$statistic)
  )

  # false recommendation rate of the single decision test stays near alpha
  rec <- vapply(1:40, function(r) {
    x1 <- rnorm(D)
    psi <- runif(D, 0.3, 1.5)
    dat <- data.frame(y = 1 + x1 + rnorm(D, sd = 0.8) + rnorm(D, sd = sqrt(psi)), x1 = x1, psi = psi)
    suppressMessages(spatial_test(y ~ x1, vardir = "psi", data = dat, W = W))$spatial_recommended
  }, logical(1))
  expect_lte(mean(rec), 0.15)
})
