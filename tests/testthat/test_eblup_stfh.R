library(testthat)
library(fastsae)

skip_if_not_installed("sae")

# ------------------------------------------------------------------
# Shared objects
# ------------------------------------------------------------------

library(dplyr)
# mys_panel is stored time-major; sae::eblupSTFH() requires domain-major
# order (all years of area 1, then area 2, ...), so sort before comparing.
mys_panel_nona <- mys_panel |>
  filter(!is.na(y)) |>
  arrange(area, year)
mys_proxmat_nona <- mys_proxmat[-c(21, 25), -c(21, 25)]

# With the correct (domain-major) order, sae::eblupSTFH() does not converge on
# the observed mys_panel response, so the agreement tests use a response
# simulated from the ST model on the same covariates, vardir and W.
local({
  D <- 40
  Tt <- 11
  set.seed(2)
  W <- mys_proxmat_nona
  Om1 <- solve(t(diag(D) - 0.5 * W) %*% (diag(D) - 0.5 * W))
  Om2 <- outer(1:Tt, 1:Tt, function(i, j) 0.5^abs(i - j)) / (1 - 0.5^2)
  u1 <- sqrt(0.5) * t(chol(Om1)) %*% rnorm(D)
  u2 <- sqrt(0.3) * as.vector(t(chol(Om2)) %*% matrix(rnorm(D * Tt), Tt))
  X <- model.matrix(~ x1 + x2 + x3, mys_panel_nona)
  mys_panel_nona$y <<- as.vector(X %*% c(4, 0, 0.05, 0.01)) + u1[rep(1:D, each = Tt)] + u2 +
    rnorm(D * Tt, sd = sqrt(mys_panel_nona$vardir))
})


fit_fast <- eblup_stfh(
  y ~ x1 + x2 + x3,
  data = mys_panel_nona,
  vardir = ~vardir,
  domain = ~area,
  time = ~year,
  W = mys_proxmat_nona,
  model = "ST",
  print_result = FALSE
)

fit_sae <- sae::eblupSTFH(
  mys_panel_nona$y ~ mys_panel_nona$x1 + mys_panel_nona$x2 + mys_panel_nona$x3,
  vardir = mys_panel_nona$vardir,
  D = length(unique(mys_panel_nona$area)),
  T = length(unique(mys_panel_nona$year)),
  proxmat = mys_proxmat_nona
)


tol <- 1e-6


# -----------------------------------------------------------------
# Structure
# ------------------------------------------------------------------

test_that("eblup_stfh returns valid structure matching seblup_area", {
  expect_true(is.list(fit_fast))

  # Top-level fields (matching seblup_area structure)
  expect_true("estcoef" %in% names(fit_fast))
  expect_true("goodness" %in% names(fit_fast))
  expect_true("df_eblup" %in% names(fit_fast))
  expect_true("convergence" %in% names(fit_fast))
  expect_true("n_iter" %in% names(fit_fast))

  # Additional for spatio-temporal
  expect_true("estvarcomp" %in% names(fit_fast))
  expect_true("model" %in% names(fit_fast))
  expect_true("formula" %in% names(fit_fast))

  # df_eblup has eblup column
  expect_true("eblup" %in% names(fit_fast$df_eblup))
  expect_true("random_effect_u1" %in% names(fit_fast$df_eblup))
  expect_true("random_effect_u2" %in% names(fit_fast$df_eblup))

  expect_length(fit_fast$df_eblup$eblup, nrow(mys_panel_nona))
})

# ------------------------------------------------------------------
# EBLUP
# ------------------------------------------------------------------

test_that("EBLUP agrees with sae::eblupSTFH", {
  expect_equal(
    fit_fast$df_eblup$eblup,
    as.numeric(fit_sae$eblup),
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# Variance component (estvarcomp - spatio-temporal specific)
# ------------------------------------------------------------------

test_that("estvarcomp agrees with sae::eblupSTFH", {
  expect_equal(
    fit_fast$estvarcomp$estimate,
    as.numeric(fit_sae$fit$estvarcomp$estimate),
    tolerance = tol
  )
  expect_equal(
    fit_fast$estvarcomp$std.error,
    as.numeric(fit_sae$fit$estvarcomp$std.error),
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# Goodness of fit
# ------------------------------------------------------------------

test_that("goodness statistics agree with sae::eblupSTFH", {
  expect_named(fit_fast$goodness, c("loglikelihood", "AIC", "BIC"))
  expect_equal(
    as.numeric(fit_fast$goodness),
    as.numeric(fit_sae$fit$goodness),
    tolerance = tol
  )
})

test_that("the log-likelihood matches the dense Gaussian likelihood", {
  D <- length(unique(mys_panel_nona$area))
  Tt <- length(unique(mys_panel_nona$year))
  X <- model.matrix(~ x1 + x2 + x3, mys_panel_nona)
  th <- fit_fast$estvarcomp$estimate
  pr <- .stfh_dense_parts(th, as.matrix(mys_proxmat_nona), D, Tt)
  ll <- .loglik_dense(mys_panel_nona$y, X, pr$G + diag(mys_panel_nona$vardir))
  expect_equal(unname(fit_fast$goodness["loglikelihood"]), ll, tolerance = 1e-6)
})

# ------------------------------------------------------------------
# Regression coefficients
# ------------------------------------------------------------------

test_that("beta estimates agree with sae::eblupSTFH", {
  expect_identical(
    names(fit_fast$estcoef$beta),
    names(fit_sae$fit$estcoef$beta)
  )

  expect_equal(
    fit_fast$estcoef$beta,
    fit_sae$fit$estcoef$beta,
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# Beta Standard errors
# ------------------------------------------------------------------

test_that("beta standard errors agree with sae::eblupSTFH", {
  expect_equal(
    fit_fast$estcoef$std.error,
    fit_sae$fit$estcoef$std.error,
    tolerance = tol
  )
})


# ------------------------------------------------------------------
# Error handling
# ------------------------------------------------------------------

test_that("non-square proximity matrix throws error", {
  W_bad <- mys_proxmat_nona[-1, ]

  expect_error(
    eblup_stfh(
      y ~ x1 + x2 + x3,
      data = mys_panel_nona,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = W_bad,
      model = "ST"
    )
  )
})

test_that("negative sampling variance throws error", {
  dat_bad <- mys_panel_nona
  dat_bad$vardir[1] <- -1

  expect_error(
    eblup_stfh(
      y ~ x1 + x2 + x3,
      data = dat_bad,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = mys_proxmat_nona,
      model = "ST"
    )
  )
})

test_that("eblup_stfh validation and edge cases", {
  # Unsampled areas (NA in response) are predicted, not rejected
  # (see test-eblup_stfh_unsampled.R)
  fit_na <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat,
    model = "ST",
    print_result = FALSE
  )
  expect_equal(sum(!fit_na$df_eblup$sampled), sum(is.na(mys_panel$y)))

  # W contains NA
  W_na <- mys_proxmat_nona
  W_na[1, 2] <- NA
  expect_error(
    eblup_stfh(
      y ~ x1 + x2 + x3,
      data = mys_panel_nona,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = W_na,
      model = "ST",
      print_result = FALSE
    ),
    "Argument W contains NA"
  )

  # Dimension mismatch
  bad_panel <- mys_panel_nona[-1, ]
  expect_error(
    eblup_stfh(
      y ~ x1 + x2 + x3,
      data = bad_panel,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = mys_proxmat_nona,
      model = "ST",
      print_result = FALSE
    ),
    "Dimensions mismatch"
  )

  # Custom starting values
  fit_custom_start <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    sigma21_start = 1.0,
    sigma22_start = 0.5,
    rho1_start = 0.3,
    rho2_start = 0.4,
    print_result = FALSE
  )
  expect_s3_class(fit_custom_start, "fastsae")

  # Non-convergence: maxiter = 1 and an unattainable precision (the
  # constrained fallback, with its 5 * maxiter budget, cannot reach it either)
  fit_st_nc <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat_nona,
    model = "ST",
    maxiter = 1,
    precision = 1e-14,
    print_result = FALSE
  )
  expect_false(fit_st_nc$convergence)

  # print_result = TRUE with model = "S"
  expect_output(
    eblup_stfh(
      y ~ x1 + x2 + x3,
      data = mys_panel_nona,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = mys_proxmat_nona,
      model = "S",
      print_result = TRUE
    )
  )

  # Parameter bounds errors
  expect_error(
    eblup_stfh(y ~ x1 + x2,
      data = mys_panel_nona, vardir = ~vardir, domain = ~area, time = ~year,
      W = mys_proxmat_nona, sigma21_start = -1
    ),
    "sigma21_start must be >= 0"
  )
  expect_error(
    eblup_stfh(y ~ x1 + x2,
      data = mys_panel_nona, vardir = ~vardir, domain = ~area, time = ~year,
      W = mys_proxmat_nona, sigma22_start = -1
    ),
    "sigma22_start must be >= 0"
  )
  expect_error(
    eblup_stfh(y ~ x1 + x2,
      data = mys_panel_nona, vardir = ~vardir, domain = ~area, time = ~year,
      W = mys_proxmat_nona, rho1_start = -1.5
    ),
    "rho1_start must be in the interval"
  )
  expect_error(
    eblup_stfh(y ~ x1 + x2,
      data = mys_panel_nona, vardir = ~vardir, domain = ~area, time = ~year,
      W = mys_proxmat_nona, model = "ST", rho2_start = 1.5
    ),
    "rho2_start must be in the interval"
  )

  # Parametric bootstrap MSE
  expect_output(
    fit_pb <- eblup_stfh(
      y ~ x1 + x2,
      data = mys_panel_nona,
      vardir = ~vardir,
      domain = ~area,
      time = ~year,
      W = mys_proxmat_nona,
      model = "ST",
      mse_type = "bootstrap",
      B = 5,
      seed = 42,
      print_result = TRUE
    ),
    "Fixed Effects Coefficients"
  )
  expect_true("mse_pb" %in% names(fit_pb$df_eblup))
  expect_true(all(fit_pb$df_eblup$mse_pb >= 0))
  expect_equal(fit_pb$B, 5)

  # Direct nrow(mf) != length(vardir) check
  testthat::with_mocked_bindings(
    .get_variable = function(data, variable) c(1, 2),
    .package = "fastsae",
    {
      expect_error(
        eblup_stfh(
          y ~ x1 + x2,
          data = mys_panel_nona,
          vardir = ~vardir,
          domain = ~area,
          time = ~year,
          W = mys_proxmat_nona,
          print_result = FALSE
        ),
        "Length of 'vardir' must equal number of observations"
      )
    }
  )
})


test_that("eblup_stfh is invariant to the row order of the panel", {
  # time-major input (as stored in mys_panel) must give the same fit as the
  # domain-major panel, with results returned in the caller's row order
  tm <- mys_panel_nona[order(mys_panel_nona$year, mys_panel_nona$area), ]
  fit_tm <- eblup_stfh(y ~ x1 + x2 + x3,
    data = tm, vardir = ~vardir, domain = ~area,
    time = ~year, W = mys_proxmat_nona, model = "ST", print_result = FALSE
  )
  expect_equal(fit_tm$estvarcomp$estimate, fit_fast$estvarcomp$estimate, tolerance = 1e-8)
  expect_equal(fit_tm$estcoef$beta, fit_fast$estcoef$beta, tolerance = 1e-8)
  expect_equal(fit_tm$df_eblup$domain, tm$area)
  expect_equal(fit_tm$df_eblup$time, tm$year)
  key_tm <- paste(fit_tm$df_eblup$domain, fit_tm$df_eblup$time)
  key_dm <- paste(fit_fast$df_eblup$domain, fit_fast$df_eblup$time)
  expect_equal(fit_tm$df_eblup$eblup, fit_fast$df_eblup$eblup[match(key_tm, key_dm)], tolerance = 1e-8)
})

test_that("eblup_stfh does not stop at a zero variance when its correlation elsewhere is better", {
  # A variance at 0 leaves its correlation unidentified; the fit used to keep
  # the correlation wherever it was when the variance reached 0, although a
  # positive variance at another correlation had a larger REML likelihood.
  stfh_small <- function(seed) {
    sim <- sim_series_data(
      D = 12, T = 5, n_unsampled = 0, prop_intermittent = 0,
      sigma_s = 0.05, sigma_t = 0.05, rho_s = 0.4, rho_t = 0.5, seed = seed
    )
    d <- sim$data[order(sim$data$area, sim$data$year), ]
    W <- unclass(sim$W)
    W <- W / rowSums(W)
    attributes(W)[c("coords", "type", "style")] <- NULL
    list(d = d, W = W, X = cbind(1, d$x1, d$x2))
  }
  # Largest REML log-likelihood from a few L-BFGS-B starts (dense reference)
  reml_max <- function(z, model) {
    np <- if (model == "ST") 4 else 3
    f <- function(t) {
      v <- tryCatch(
        as.numeric(.stfh_dense_reml(c(t, 0)[1:4], z$X, z$d$y_gaussian, z$d$vardir, z$W, 12, 5, model)),
        error = function(e) -1e10
      )
      if (is.finite(v)) v else -1e10
    }
    best <- -Inf
    for (r in c(-0.8, 0.8)) for (s in c(0.02, 0.3)) {
      o <- optim(c(s, r, s, -r)[1:np], f,
        method = "L-BFGS-B", lower = c(0, -0.999, 0, -0.999)[1:np],
        upper = c(Inf, 0.999, Inf, 0.999)[1:np], control = list(fnscale = -1, maxit = 500)
      )
      best <- max(best, o$value)
    }
    list(max = best, f = f)
  }
  # seed 11, "ST": sigma22 stopped at 0 (rho2 0.29), REML maximum at rho2 near -0.98;
  # seed 2, "S": sigma21 stopped at 0 (rho1 0.95), REML maximum at rho1 = -0.999
  for (case in list(list(seed = 11, model = "ST"), list(seed = 2, model = "S"))) {
    z <- stfh_small(case$seed)
    fit <- suppressMessages(eblup_stfh(y_gaussian ~ x1 + x2,
      data = z$d, vardir = ~vardir, domain = ~area, time = ~year,
      W = z$W, model = case$model, mse_type = "none", print_result = FALSE
    ))
    expect_true(fit$convergence)
    r <- reml_max(z, case$model)
    expect_gt(r$f(fit$estvarcomp$estimate), r$max - 1e-3)
  }
})
