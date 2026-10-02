library(testthat)
library(fastsae)

test_that("S3 methods work for eblup_fh", {
  fit <- eblup_fh(y ~ x1 + x2, vardir = ~vardir, data = mys, method = "REML", print_result = FALSE)

  # print
  expect_no_error(print(fit))

  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.fastsae")
  expect_no_error(print(s))

  # coef, fitted, residuals
  cf <- coef(fit)
  expect_named(cf)
  expect_equal(length(cf), 3)

  ft <- fitted(fit)
  expect_equal(length(ft), nrow(mys))
  expect_equal(ft, fit$df_eblup$eblup)

  res <- residuals(fit)
  expect_equal(length(res), nrow(mys))
  expect_equal(res, mys$y - ft)

  # plot method routes to autoplot
  p <- plot(fit)
  expect_s3_class(p, "ggplot")
})

test_that("S3 methods work for eblup_sfh", {
  fit <- eblup_sfh(y ~ x1 + x2, vardir = ~vardir, data = mys, W = mys_proxmat, method = "REML", print_result = FALSE)

  expect_output(print(fit))

  s <- summary(fit)
  expect_output(print(s))

  cf <- coef(fit)
  expect_equal(length(cf), 3)

  ft <- fitted(fit)
  expect_equal(length(ft), nrow(mys))

  res <- residuals(fit)
  expect_equal(length(res), nrow(mys))
})

test_that("S3 methods work for eblup_stfh", {
  mys_panel_nona <- mys_panel[!is.na(mys_panel$y), ]
  mys_proxmat_nona <- mys_proxmat[-c(21, 25), -c(21, 25)]

  fit <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    print_result = FALSE
  )

  expect_output(print(fit))

  s <- summary(fit)
  expect_output(print(s))

  cf <- coef(fit)
  expect_true(length(cf) >= 3)

  ft <- fitted(fit)
  expect_equal(length(ft), nrow(mys_panel_nona))

  res <- residuals(fit)
  expect_equal(length(res), nrow(mys_panel_nona))
})

test_that("S3 methods work for eblup_bhf", {
  skip_if_not_installed("lme4")
  set.seed(42)
  N_d <- 10
  m <- 5
  d_id <- rep(seq_len(m), each = N_d)
  x_ij <- rnorm(m * N_d)
  u_i <- rnorm(m, 0, 1)
  y_ij <- 1 + 2 * x_ij + u_i[d_id] + rnorm(m * N_d, 0, 0.5)
  df_sample <- data.frame(y = y_ij, x = x_ij, domain = factor(d_id))
  meanxpop <- data.frame(domain = factor(seq_len(m)), x = tapply(x_ij, d_id, mean), N = rep(100, m))

  fit <- eblup_bhf(
    formula = y ~ x,
    unit_data = df_sample,
    Xpop = meanxpop,
    domain_var = "domain",
    popsize_var = "N",
    print_result = FALSE
  )

  expect_output(print(fit))

  s <- summary(fit)
  expect_output(print(s))

  cf <- coef(fit)
  expect_named(cf)
  expect_equal(length(cf), 2)

  ft <- fitted(fit)
  expect_equal(length(ft), m)

  res <- residuals(fit)
  expect_equal(length(res), nrow(df_sample))
})

test_that("S3 methods cover all print, summary, coef, fitted, residuals, and plot branches", {
  # 1. Non-converged model with two-fold variances, temporal & st_interaction, df_area
  rev_vec <- c(area = 0.5, subarea = 0.25)
  dummy_mod <- structure(
    list(
      call = quote(dummy_call()),
      convergence = FALSE,
      model = "TFH",
      temporal = "ar1",
      st_interaction = "type1",
      method = "REML",
      random_effect_var = rev_vec,
      random_effect_var_time = 0.15,
      rho = 0.35,
      rho_time = 0.45,
      phi = 0.65,
      estcoef = data.frame(beta = c(2, 3), std.error = c(0.1, 0.2), zvalue = c(20, 15), pvalue = c(0.001, 0.001), row.names = c("(Intercept)", "z")),
      hyperpar = data.frame(param = "theta", value = 1.2),
      goodness = c(AIC = 50.5),
      df_hb = data.frame(domain = 1:8, y = 1:8, hb = (1:8) * 1.1, linear_pred = 1:8, sd = rep(0.2, 8), mse = rep(0.04, 8), rse = rep(5, 8)),
      df_area = data.frame(area = 1:8, estimate = 1:8)
    ),
    class = "fastsae"
  )

  expect_no_error(print(dummy_mod))
  s_dummy <- summary(dummy_mod)
  expect_no_error(print(s_dummy))

  # 2. Model with EBP estimates and estvarcomp
  dummy_ebp <- structure(
    list(
      model = "TWOFOLD",
      convergence = TRUE,
      n_iter = 5,
      estvarcomp = data.frame(component = c("area", "subarea"), variance = c(0.4, 0.2)),
      df_ebp = data.frame(domain = 1:3, y = 10:12, ebp = 10:12, mse = rep(1, 3), rse = rep(10, 3)),
      fit = list(beta = c(1, 2))
    ),
    class = "fastsae"
  )
  expect_no_error(print(dummy_ebp))
  s_ebp <- summary(dummy_ebp)
  expect_no_error(print(s_ebp))

  # 3. Model with S model type
  dummy_s <- structure(list(model = "S", convergence = TRUE), class = "fastsae")
  expect_no_error(print(dummy_s))

  # 4. coef extraction fallbacks
  expect_equal(coef(dummy_ebp), c(1, 2))
  expect_null(coef(structure(list(), class = "fastsae")))

  # 5. fitted extraction fallbacks
  expect_equal(fitted(dummy_mod), dummy_mod$df_hb$hb)
  expect_equal(fitted(dummy_ebp), dummy_ebp$df_ebp$ebp)
  expect_null(fitted(structure(list(), class = "fastsae")))

  # 6. residuals extraction fallbacks
  expect_equal(residuals(dummy_mod), dummy_mod$df_hb$y - dummy_mod$df_hb$hb)
  expect_equal(residuals(dummy_ebp), dummy_ebp$df_ebp$y - dummy_ebp$df_ebp$ebp)
  expect_null(residuals(structure(list(), class = "fastsae")))

  # 7. plot method with two models comparison and spatial map
  fit1 <- eblup_fh(y ~ x1 + x2, vardir = ~vardir, data = mys, method = "REML", print_result = FALSE)
  fit2 <- eblup_fh(y ~ x1, vardir = ~vardir, data = mys, method = "REML", print_result = FALSE)
  p2 <- plot(fit1, fit2)
  expect_s3_class(p2, "ggplot")

  if (requireNamespace("sf", quietly = TRUE)) {
    grid_sf <- sf::st_make_grid(
      sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
      cellsize = c(1, 1),
      what = "polygons"
    )[seq_len(nrow(mys))]
    mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)
    p_map1 <- plot(fit1, sf_geom = mys_sf)
    expect_s3_class(p_map1, "ggplot")
    p_map2 <- plot(fit1, fit2, sf_geom = mys_sf)
    expect_s3_class(p_map2, "ggplot")
  }
})

test_that("print and summary print handle non-standard coefficient columns and eblup-only df", {
  mod_custom <- structure(
    list(
      model = "Custom Model",
      estcoef = data.frame(gamma = c(1.2, 0.5)),
      df_eblup = data.frame(domain = 1:3, eblup = 1:3, mse = 0.1, rse = 5)
    ),
    class = "fastsae"
  )
  expect_output(print(mod_custom), "Fixed Effects Coefficients")

  sum_mod <- structure(
    list(
      model = "Custom Summary",
      coefficients = data.frame(gamma = c(1.2, 0.5)),
      df_eblup = data.frame(domain = 1:3, eblup = 1:3, mse = 0.1, rse = 5)
    ),
    class = "summary.fastsae"
  )
  expect_output(print(sum_mod), "EBLUP Summary Statistics")
})
