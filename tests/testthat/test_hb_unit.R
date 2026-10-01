test_that("hb_unit works for Gaussian Battese-Harter-Fuller model and matches eblup_bhf", {
  skip_if_not_installed("INLA")

  data("cornsoybean", package = "fastsae")
  data("cornsoybeanmeans", package = "fastsae")

  df_pop <- cornsoybeanmeans
  names(df_pop)[names(df_pop) == "MeanCornPixPerSeg"] <- "CornPix"
  names(df_pop)[names(df_pop) == "MeanSoyBeansPixPerSeg"] <- "SoyBeansPix"

  df_sample <- cornsoybean
  names(df_sample)[names(df_sample) == "County"] <- "CountyIndex"

  # Fit EBLUP BHF
  fit_eblup <- eblup_bhf(
    formula = CornHec ~ CornPix + SoyBeansPix,
    unit_data = df_sample,
    Xpop = df_pop,
    domain_var = "CountyIndex",
    popsize_var = "PopnSegments",
    print_result = FALSE
  )

  # Fit HB Unit
  fit_hb <- hb_unit(
    formula = CornHec ~ CornPix + SoyBeansPix,
    unit_data = df_sample,
    Xpop = df_pop,
    domain_var = "CountyIndex",
    popsize_var = "PopnSegments",
    family = "gaussian",
    print_result = FALSE
  )

  expect_s3_class(fit_hb, "fastsae_hb_unit")
  expect_s3_class(fit_hb, "fastsae_hb")
  expect_s3_class(fit_hb, "fastsae")

  # Verify df_hb columns
  expected_cols <- c(
    "domain", "y", "hb", "linear_pred", "vardir", "sd", "mse", "rse",
    "ci_lower", "ci_upper", "samp_size", "sample_mean", "random_effect",
    "pop_size", "estimated_total"
  )
  expect_true(all(expected_cols %in% names(fit_hb$df_hb)))
  expect_equal(nrow(fit_hb$df_hb), nrow(df_pop))

  # Correlation between frequentist BHF and Bayesian HB should be high (> 0.90)
  r_val <- stats::cor(fit_eblup$df_eblup$eblup, fit_hb$df_hb$hb)
  expect_gt(r_val, 0.90)

  # Fixed effect coefficients should be consistent with EBLUP
  expect_equal(nrow(fit_hb$estcoef), 3)
  expect_true(all(c("beta", "std.error", "zvalue", "pvalue") %in% names(fit_hb$estcoef)))
  expect_lt(abs(fit_hb$estcoef["CornPix", "beta"] - fit_eblup$estcoef["CornPix", "beta"]), 0.05)

  # Variance components
  expect_true(!is.na(fit_hb$random_effect_var))
  expect_true(!is.na(fit_hb$residual_var))
  expect_gt(fit_hb$random_effect_var, 0)
  expect_gt(fit_hb$residual_var, 0)

  # Check methods
  cf <- coef(fit_hb)
  expect_length(cf, 3)
  ft <- fitted(fit_hb)
  expect_equal(ft, fit_hb$df_hb$hb)
  sm <- summary(fit_hb)
  expect_s3_class(sm, "summary.fastsae")
})

test_that("hb_unit handles unsampled domains with synthetic prediction", {
  skip_if_not_installed("INLA")

  data("cornsoybean", package = "fastsae")
  data("cornsoybeanmeans", package = "fastsae")

  df_pop <- cornsoybeanmeans
  names(df_pop)[names(df_pop) == "MeanCornPixPerSeg"] <- "CornPix"
  names(df_pop)[names(df_pop) == "MeanSoyBeansPixPerSeg"] <- "SoyBeansPix"

  df_sample <- cornsoybean
  names(df_sample)[names(df_sample) == "County"] <- "CountyIndex"

  # Remove CountyIndex 1 and 2 from sample data
  df_sample_sub <- df_sample[!df_sample$CountyIndex %in% c(1, 2), ]

  fit_sub <- hb_unit(
    formula = CornHec ~ CornPix + SoyBeansPix,
    unit_data = df_sample_sub,
    Xpop = df_pop,
    domain_var = "CountyIndex",
    popsize_var = "PopnSegments",
    print_result = FALSE
  )

  # Unsampled domains identified
  expect_setequal(fit_sub$unsampled_domains, c("1", "2"))

  # Verify predictions exist for unsampled domains
  row1 <- fit_sub$df_hb[fit_sub$df_hb$domain == 1, ]
  row2 <- fit_sub$df_hb[fit_sub$df_hb$domain == 2, ]

  expect_equal(row1$samp_size, 0)
  expect_equal(row2$samp_size, 0)
  expect_true(is.na(row1$sample_mean))
  expect_true(is.na(row2$sample_mean))
  expect_false(is.na(row1$hb))
  expect_false(is.na(row2$hb))
  expect_false(is.na(row1$sd))
  expect_false(is.na(row2$sd))
})

test_that("hb_unit works without popsize_var (superpopulation mean)", {
  skip_if_not_installed("INLA")

  data("cornsoybean", package = "fastsae")
  data("cornsoybeanmeans", package = "fastsae")

  df_pop <- cornsoybeanmeans
  names(df_pop)[names(df_pop) == "MeanCornPixPerSeg"] <- "CornPix"
  names(df_pop)[names(df_pop) == "MeanSoyBeansPixPerSeg"] <- "SoyBeansPix"

  df_sample <- cornsoybean
  names(df_sample)[names(df_sample) == "County"] <- "CountyIndex"

  fit_nopop <- hb_unit(
    formula = CornHec ~ CornPix + SoyBeansPix,
    unit_data = df_sample,
    Xpop = df_pop,
    domain_var = "CountyIndex",
    popsize_var = NULL,
    print_result = FALSE
  )

  expect_s3_class(fit_nopop, "fastsae_hb_unit")
  expect_false("pop_size" %in% names(fit_nopop$df_hb))
  expect_false("estimated_total" %in% names(fit_nopop$df_hb))
  expect_equal(nrow(fit_nopop$df_hb), nrow(df_pop))
  expect_false(any(is.na(fit_nopop$df_hb$hb)))
})

test_that("hb_unit works for Binomial (logistic) and Poisson unit models", {
  skip_if_not_installed("INLA")

  set.seed(123)
  D <- 8
  n_per_d <- 20
  domain <- rep(1:D, each = n_per_d)
  x1 <- stats::rnorm(D * n_per_d)
  u <- stats::rnorm(D, 0, 0.4)
  eta <- 0.2 + 0.6 * x1 + u[domain]

  # 1. Binomial
  p <- 1 / (1 + exp(-eta))
  y_bin <- stats::rbinom(D * n_per_d, 1, p)
  df_sample_bin <- data.frame(y = y_bin, x1 = x1, domain = domain)
  Xpop_bin <- data.frame(
    domain = 1:D,
    x1 = tapply(x1, domain, mean),
    N = rep(500, D)
  )

  fit_bin <- hb_unit(
    formula = y ~ x1,
    unit_data = df_sample_bin,
    Xpop = Xpop_bin,
    domain_var = "domain",
    popsize_var = "N",
    family = "binomial",
    print_result = FALSE
  )

  expect_s3_class(fit_bin, "fastsae_hb_unit")
  expect_equal(fit_bin$family, "binomial")
  expect_true(all(fit_bin$df_hb$hb >= 0 & fit_bin$df_hb$hb <= 1))
  expect_true(all(fit_bin$df_hb$ci_lower >= 0 & fit_bin$df_hb$ci_upper <= 1))

  # 2. Poisson
  lambda <- exp(0.5 + 0.3 * x1 + u[domain])
  y_pois <- stats::rpois(D * n_per_d, lambda)
  df_sample_pois <- data.frame(y = y_pois, x1 = x1, domain = domain)

  fit_pois <- hb_unit(
    formula = y ~ x1,
    unit_data = df_sample_pois,
    Xpop = Xpop_bin,
    domain_var = "domain",
    popsize_var = "N",
    family = "poisson",
    print_result = FALSE
  )

  expect_s3_class(fit_pois, "fastsae_hb_unit")
  expect_equal(fit_pois$family, "poisson")
  expect_true(all(fit_pois$df_hb$hb >= 0))
  expect_true(all(fit_pois$df_hb$ci_lower >= 0))
})

test_that("hb_unit works with spatial BYM2 structure", {
  skip_if_not_installed("INLA")

  data("mys", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  # Create synthetic unit records from mys
  set.seed(42)
  n_units <- 4
  sim_units <- do.call(rbind, lapply(seq_len(nrow(mys)), function(i) {
    data.frame(
      area = mys$area[i],
      x1 = mys$x1[i] + stats::rnorm(n_units, 0, 0.1),
      x2 = mys$x2[i] + stats::rnorm(n_units, 0, 0.1),
      y = mys$y[i] + stats::rnorm(n_units, 0, 0.3)
    )
  }))

  Xpop_spatial <- mys[, c("area", "x1", "x2")]
  Xpop_spatial$N <- 1000

  fit_sp <- hb_unit(
    formula = y ~ x1 + x2,
    unit_data = sim_units,
    Xpop = Xpop_spatial,
    domain_var = "area",
    popsize_var = "N",
    spatial = "bym2",
    W = mys_proxmat,
    print_result = FALSE
  )

  expect_s3_class(fit_sp, "fastsae_hb_unit")
  expect_equal(fit_sp$spatial, "bym2")
  expect_false(is.na(fit_sp$phi))
  expect_equal(nrow(fit_sp$df_hb), nrow(mys))
})

test_that("hb_unit input validation throws informative errors", {
  skip_if_not_installed("INLA")

  data("cornsoybean", package = "fastsae")
  data("cornsoybeanmeans", package = "fastsae")

  df_pop <- cornsoybeanmeans
  names(df_pop)[names(df_pop) == "MeanCornPixPerSeg"] <- "CornPix"
  names(df_pop)[names(df_pop) == "MeanSoyBeansPixPerSeg"] <- "SoyBeansPix"

  df_sample <- cornsoybean
  names(df_sample)[names(df_sample) == "County"] <- "CountyIndex"

  # Missing domain_var in unit_data
  expect_error(
    hb_unit(CornHec ~ CornPix, unit_data = df_sample, Xpop = df_pop, domain_var = "NonExistent"),
    "Domain identifier"
  )

  # Missing domain_var in Xpop
  expect_error(
    hb_unit(CornHec ~ CornPix, unit_data = df_sample, Xpop = df_sample, domain_var = "CountyName"),
    "Domain identifier"
  )

  # Spatial without W
  expect_error(
    hb_unit(CornHec ~ CornPix, unit_data = df_sample, Xpop = df_pop, domain_var = "CountyIndex", spatial = "bym2"),
    "Spatial weight"
  )

  # Missing domain in Xpop
  df_pop_incomplete <- df_pop[-1, ]
  expect_error(
    hb_unit(CornHec ~ CornPix + SoyBeansPix, unit_data = df_sample, Xpop = df_pop_incomplete, domain_var = "CountyIndex"),
    "missing from"
  )
})

test_that("hb_unit integrates with diagnose() and compare_sae()", {
  skip_if_not_installed("INLA")

  data("cornsoybean", package = "fastsae")
  data("cornsoybeanmeans", package = "fastsae")

  df_pop <- cornsoybeanmeans
  names(df_pop)[names(df_pop) == "MeanCornPixPerSeg"] <- "CornPix"
  names(df_pop)[names(df_pop) == "MeanSoyBeansPixPerSeg"] <- "SoyBeansPix"

  df_sample <- cornsoybean
  names(df_sample)[names(df_sample) == "County"] <- "CountyIndex"

  fit_hb <- hb_unit(
    formula = CornHec ~ CornPix + SoyBeansPix,
    unit_data = df_sample,
    Xpop = df_pop,
    domain_var = "CountyIndex",
    popsize_var = "PopnSegments",
    print_result = FALSE
  )

  # diagnose() works
  diag <- diagnose(fit_hb)
  expect_s3_class(diag, "fastsae_diagnose")
  expect_true(is.list(diag$precision))
  expect_true(is.list(diag$brown_test))

  # compare_sae() works
  fit_eblup <- eblup_bhf(
    formula = CornHec ~ CornPix + SoyBeansPix,
    unit_data = df_sample,
    Xpop = df_pop,
    domain_var = "CountyIndex",
    popsize_var = "PopnSegments",
    print_result = FALSE
  )

  comp <- compare_sae(fit_eblup, fit_hb)
  expect_s3_class(comp, "fastsae_comparison")
  cor_val <- as.numeric(comp$metrics[comp$metrics$Metric == "Pearson Correlation (r)", "Value"])
  expect_gt(cor_val, 0.80)
})
