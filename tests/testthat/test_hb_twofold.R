library(testthat)
library(fastsae)

test_that("hb_twofold works for Gaussian two-fold model and matches eblup_tfh", {
  skip_if_not_installed("INLA")

  set.seed(42)
  m <- 12
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:4, 1)
    data.frame(
      area = paste0("Area_", d),
      subarea = paste0("Area_", d, "_sub", seq_len(nd)),
      x1 = stats::rnorm(nd),
      vardir = stats::runif(nd, 0.2, 0.8)
    )
  }))
  dat$y <- 1 + 1.5 * dat$x1 + stats::rnorm(m)[as.integer(factor(dat$area))] +
    stats::rnorm(nrow(dat), sd = 0.4) + stats::rnorm(nrow(dat), sd = sqrt(dat$vardir))

  # Fit EBLUP Twofold
  fit_eblup <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                         data = dat, print_result = FALSE)

  # Fit HB Twofold
  fit_hb <- hb_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                       data = dat, print_result = FALSE)

  expect_s3_class(fit_hb, "fastsae_hb_twofold")
  expect_s3_class(fit_hb, "fastsae_hb")
  expect_s3_class(fit_hb, "fastsae")

  # Sub-area columns
  expected_cols_sub <- c(
    "domain", "subarea", "y", "hb", "linear_pred", "vardir", "sd", "mse",
    "rse", "ci_lower", "ci_upper", "random_effect_area", "random_effect_subarea"
  )
  expect_true(all(expected_cols_sub %in% names(fit_hb$df_hb)))
  expect_equal(nrow(fit_hb$df_hb), nrow(dat))

  # Area columns
  expected_cols_area <- c(
    "domain", "hb_area", "sd_area", "mse_area", "rse_area",
    "ci_lower_area", "ci_upper_area", "n_subareas"
  )
  expect_true(all(expected_cols_area %in% names(fit_hb$df_area)))
  expect_equal(nrow(fit_hb$df_area), m)

  # High correlation with frequentist EBLUP
  r_val <- stats::cor(fit_eblup$df_eblup$eblup, fit_hb$df_hb$hb)
  expect_gt(r_val, 0.90)

  # Fixed effect coefficients consistent
  expect_equal(nrow(fit_hb$estcoef), 2)
  expect_lt(abs(fit_hb$estcoef["x1", "beta"] - fit_eblup$estcoef["x1", "beta"]), 0.15)

  # Variance components
  s2 <- fit_hb$random_effect_var
  expect_named(s2, c("sigma2_v", "sigma2_u"))
  expect_true(!is.na(s2[["sigma2_v"]]))
  expect_gt(s2[["sigma2_v"]], 0)

  # Methods
  expect_length(coef(fit_hb), 2)
  expect_equal(fitted(fit_hb), fit_hb$df_hb$hb)
  expect_equal(residuals(fit_hb), fit_hb$df_hb$y - fit_hb$df_hb$hb)
  expect_s3_class(summary(fit_hb), "summary.fastsae")
})

test_that("hb_twofold handles unsampled subareas with synthetic prediction", {
  skip_if_not_installed("INLA")

  set.seed(11)
  m <- 8
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:4, 1)
    data.frame(
      area = d,
      subarea = paste0(d, "-", seq_len(nd)),
      x1 = stats::rnorm(nd),
      vardir = stats::runif(nd, 0.2, 0.8)
    )
  }))
  dat$y <- 1 + 1.2 * dat$x1 + stats::rnorm(m)[dat$area] +
    stats::rnorm(nrow(dat), sd = sqrt(dat$vardir))
  dat$y[c(3, 8)] <- NA_real_

  fit_sub <- hb_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                        data = dat, print_result = FALSE)

  expect_true(all(is.finite(fit_sub$df_hb$hb)))
  expect_true(all(is.finite(fit_sub$df_hb$sd)))
  expect_false(is.na(fit_sub$df_hb$hb[3]))
  expect_false(is.na(fit_sub$df_hb$hb[8]))
  expect_true(is.na(fit_sub$df_hb$y[3]))
  expect_true(is.na(fit_sub$df_hb$y[8]))
})

test_that("hb_twofold computes area-level aggregates correctly", {
  skip_if_not_installed("INLA")

  set.seed(42)
  m <- 5
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- 3
    data.frame(
      area = d,
      subarea = seq_len(nd),
      x1 = stats::rnorm(nd),
      vardir = rep(0.5, nd),
      w = c(0.2, 0.3, 0.5)
    )
  }))
  dat$y <- 2 + dat$x1 + stats::rnorm(nrow(dat), sd = sqrt(dat$vardir))

  fit_agg <- hb_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                        weight = "w", data = dat, print_result = FALSE)

  expect_equal(nrow(fit_agg$df_area), m)
  for (d in 1:m) {
    sub_hb <- fit_agg$df_hb$hb[fit_agg$df_hb$domain == d]
    sub_w <- dat$w[dat$area == d]
    expected_area_hb <- sum(sub_w * sub_hb) / sum(sub_w)
    expect_equal(fit_agg$df_area$hb_area[d], expected_area_hb, tolerance = 1e-6)
  }
})

test_that("hb_twofold works with spatial BYM2 structure on area level", {
  skip_if_not_installed("INLA")

  set.seed(42)
  m <- 6
  W <- matrix(0, m, m)
  for (i in 1:(m-1)) { W[i, i+1] <- 1; W[i+1, i] <- 1 }
  rownames(W) <- colnames(W) <- paste0("Area_", 1:m)

  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:4, 1)
    data.frame(
      area = paste0("Area_", d),
      subarea = seq_len(nd),
      x1 = stats::rnorm(nd),
      vardir = stats::runif(nd, 0.2, 0.8)
    )
  }))
  dat$y <- 1 + dat$x1 + stats::rnorm(nrow(dat), sd = sqrt(dat$vardir))

  fit_sp <- hb_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                       spatial = "bym2", W = W, data = dat, print_result = FALSE)

  expect_s3_class(fit_sp, "fastsae_hb_twofold")
  expect_equal(fit_sp$spatial, "bym2")
  expect_false(is.na(fit_sp$phi))
  expect_true(fit_sp$phi >= 0 && fit_sp$phi <= 1)
})

test_that("hb_twofold works for Binomial and Poisson response models", {
  skip_if_not_installed("INLA")

  set.seed(123)
  m <- 6
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- 3
    data.frame(
      area = d,
      subarea = seq_len(nd),
      x1 = stats::rnorm(nd),
      n_trials = sample(50:100, nd, replace = TRUE),
      exposure = sample(100:200, nd, replace = TRUE)
    )
  }))

  # 1. Binomial
  p <- 1 / (1 + exp(-(-0.5 + 0.8 * dat$x1)))
  dat$y_bin <- stats::rbinom(nrow(dat), size = dat$n_trials, prob = p)

  fit_bin <- hb_twofold(y_bin ~ x1, domain = "area", subarea = "subarea",
                        trials = "n_trials", family = "binomial", data = dat,
                        print_result = FALSE)

  expect_s3_class(fit_bin, "fastsae_hb_twofold")
  expect_equal(fit_bin$family, "binomial")
  expect_true(all(fit_bin$df_hb$hb >= 0 & fit_bin$df_hb$hb <= 1))

  # 2. Poisson
  lambda <- exp(0.2 + 0.5 * dat$x1)
  dat$y_pois <- stats::rpois(nrow(dat), lambda = dat$exposure * lambda)

  fit_pois <- hb_twofold(y_pois ~ x1, domain = "area", subarea = "subarea",
                         exposure = "exposure", family = "poisson", data = dat,
                         print_result = FALSE)

  expect_s3_class(fit_pois, "fastsae_hb_twofold")
  expect_equal(fit_pois$family, "poisson")
  expect_true(all(fit_pois$df_hb$hb >= 0))
})

test_that("hb_twofold input validation throws informative errors", {
  skip_if_not_installed("INLA")

  dat <- data.frame(area = c(1, 1, 2), subarea = 1:3, y = 1:3, vardir = c(0.1, 0.2, 0.3))

  # Missing domain
  expect_error(hb_twofold(y ~ 1, vardir = "vardir", data = dat), "domain")

  # Missing vardir for Gaussian
  expect_error(hb_twofold(y ~ 1, domain = "area", data = dat), "vardir")

  # Spatial without W
  expect_error(hb_twofold(y ~ 1, vardir = "vardir", domain = "area", spatial = "bym2", data = dat), "Spatial weight")
})

test_that("hb_twofold integrates with compare_sae() and diagnose()", {
  skip_if_not_installed("INLA")

  set.seed(42)
  m <- 8
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:4, 1)
    data.frame(
      area = d,
      subarea = paste0(d, "-", seq_len(nd)),
      x1 = stats::rnorm(nd),
      vardir = stats::runif(nd, 0.2, 0.8)
    )
  }))
  dat$y <- 1 + 1.2 * dat$x1 + stats::rnorm(m)[dat$area] +
    stats::rnorm(nrow(dat), sd = sqrt(dat$vardir))

  fit_eblup <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                         data = dat, print_result = FALSE)
  fit_hb <- hb_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                       data = dat, print_result = FALSE)

  # diagnose
  diag <- diagnose(fit_hb)
  expect_s3_class(diag, "fastsae_diagnose")

  # compare_sae
  comp <- compare_sae(fit_eblup, fit_hb)
  expect_s3_class(comp, "fastsae_comparison")
  cor_val <- as.numeric(comp$metrics[comp$metrics$Metric == "Pearson Correlation (r)", "Value"])
  expect_gt(cor_val, 0.90)
})

test_that("hb_tfh alias is identical to hb_twofold", {
  expect_identical(hb_tfh, hb_twofold)
})

test_that("hb_twofold handles besag spatial and subarea = NULL", {
  skip_if_not_installed("INLA")
  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  # Create synthetic subarea data
  set.seed(42)
  dat <- do.call(rbind, lapply(seq_len(nrow(mys)), function(d) {
    nd <- 2
    data.frame(
      area = mys$area[d],
      x1 = mys$x1[d] + stats::rnorm(nd, 0, 0.1),
      vardir = stats::runif(nd, 0.2, 0.8),
      y = mys$y[d] + stats::rnorm(nd, 0, 0.2)
    )
  }))

  # 1. subarea = NULL and spatial = besag
  fit_besag <- hb_twofold(
    y ~ x1,
    vardir = "vardir",
    domain = "area",
    subarea = NULL,
    spatial = "besag",
    W = mys_proxmat,
    data = dat,
    print_result = FALSE
  )
  expect_s3_class(fit_besag, "fastsae_hb_twofold")
  expect_equal(fit_besag$spatial, "besag")
  expect_equal(nrow(fit_besag$df_hb), nrow(dat))
})

test_that("hb_twofold supports binomial, poisson, and fallback aggregation", {
  skip_if_not_installed("INLA")
  set.seed(42)
  dat_bin <- data.frame(
    area = rep(c("A1", "A2"), each = 3),
    subarea = paste0("sub_", 1:6),
    x = rnorm(6),
    n = c(50, 60, 40, 55, 65, 45),
    y = c(10, 15, 8, 12, 18, 9)
  )

  # Binomial with compute_area = TRUE
  fit_bin <- hb_twofold(
    y ~ x,
    domain = "area",
    subarea = "subarea",
    trials = "n",
    family = "binomial",
    data = dat_bin,
    compute_area = TRUE,
    n_samples = 50,
    print_result = FALSE
  )
  expect_s3_class(fit_bin, "fastsae_hb_twofold")
  expect_true(!is.null(fit_bin$df_area))

  # Poisson with compute_area = TRUE
  dat_pois <- data.frame(
    area = rep(c("A1", "A2"), each = 3),
    subarea = paste0("sub_", 1:6),
    x = rnorm(6),
    y = c(3, 5, 2, 4, 6, 3)
  )
  fit_pois <- hb_twofold(
    y ~ x,
    domain = "area",
    subarea = "subarea",
    family = "poisson",
    data = dat_pois,
    compute_area = TRUE,
    n_samples = 50,
    print_result = FALSE
  )
  expect_s3_class(fit_pois, "fastsae_hb_twofold")

  # Fallback aggregation when n_samples = 0
  dat_gauss <- data.frame(
    area = rep(c("A1", "A2"), each = 3),
    subarea = paste0("sub_", 1:6),
    x = rnorm(6),
    vardir = rep(0.5, 6),
    y = rnorm(6, 5, 1)
  )
  expect_output(
    fit_fallback <- hb_twofold(
      y ~ x,
      vardir = "vardir",
      domain = "area",
      subarea = "subarea",
      data = dat_gauss,
      compute_area = TRUE,
      n_samples = 0,
      print_result = TRUE
    ),
    "Area-Level Aggregates"
  )
  expect_s3_class(fit_fallback, "fastsae_hb_twofold")
  expect_true(!is.null(fit_fallback$df_area))
})

test_that("hb_twofold handles INLA failure and fallback aggregation", {
  skip_if_not_installed("INLA")

  dat <- data.frame(
    area = rep(c("A1", "A2"), each = 3),
    subarea = paste0("sub_", 1:6),
    x = rnorm(6),
    vardir = rep(0.5, 6),
    y = exp(rnorm(6))
  )

  # 1. INLA error -> cli_abort (L296-299)
  testthat::with_mocked_bindings(
    inla = function(...) stop("Simulated INLA crash"),
    .package = "INLA",
    {
      expect_error(
        hb_twofold(
          y ~ x,
          vardir = "vardir",
          domain = "area",
          subarea = "subarea",
          data = dat,
          print_result = FALSE
        ),
        "Fitting two-fold model with INLA failed"
      )
    }
  )

  # 2. Fallback aggregation when n_samples=0 (L394 ps_ok=FALSE path)
  fit_fb <- hb_twofold(
    y ~ x,
    vardir = "vardir",
    domain = "area",
    subarea = "subarea",
    data = dat,
    n_samples = 0L,
    print_result = FALSE
  )
  expect_s3_class(fit_fb, "fastsae")
  expect_false(any(is.na(fit_fb$df_area$hb)))
})


