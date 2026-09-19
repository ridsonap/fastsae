test_that("S3 methods work for eblup_fh", {
  fit <- eblup_fh(y ~ x1 + x2, vardir = ~vardir, data = mys, method = "REML")
  
  # print
  expect_output(print(fit), "Fay-Herriot")
  
  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.fastsae")
  expect_output(print(s), "sigma2_u")
  
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
})

test_that("S3 methods work for eblup_sfh", {
  fit <- eblup_sfh(y ~ x1 + x2, vardir = ~vardir, data = mys, W = mys_proxmat, method = "REML")
  
  expect_output(print(fit), "Spatial Fay-Herriot")
  s <- summary(fit)
  expect_output(print(s), "rho")
  
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
  
  expect_output(print(fit), "Spatio-Temporal Fay-Herriot")
  s <- summary(fit)
  expect_output(print(s), "Variance & Correlation Components")
  
  cf <- coef(fit)
  expect_true(length(cf) >= 3)
  
  ft <- fitted(fit)
  expect_equal(length(ft), nrow(mys_panel_nona))
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

  expect_output(print(fit), "Battese-Harter-Fuller")
  s <- summary(fit)
  expect_output(print(s), "Variance Components")
  
  cf <- coef(fit)
  expect_named(cf)
  expect_equal(length(cf), 2)
  
  ft <- fitted(fit)
  expect_equal(length(ft), m)
  
  res <- residuals(fit)
  expect_equal(length(res), nrow(df_sample))
})
