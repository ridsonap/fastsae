test_that("sim_series_data produces valid spatio-temporal panel dataset", {
  sim <- sim_series_data(D = 25, T = 4, time_start = 2020, seed = 123)

  expect_s3_class(sim, "fastsae_sim_series")
  expect_named(sim, c("data", "W", "W_std", "coords", "u_spatial", "u_temporal", "parameters"))

  df <- sim$data
  expect_equal(nrow(df), 25 * 4)
  expect_true(all(c("area", "year", "x1", "x2", "y_gaussian", "vardir",
                    "y_poisson", "exposure", "y_binomial", "trials",
                    "y_beta", "y_nbinomial", "y_gamma", "x_coord", "y_coord") %in% names(df)))

  # Check dimensions of random effects
  expect_equal(length(sim$u_spatial), 25)
  expect_equal(dim(sim$u_temporal), c(25, 4))
  expect_equal(dim(sim$W), c(25, 25))
  expect_equal(dim(sim$W_std), c(25, 25))

  # Check domain-major sorting
  expect_equal(df$area[1:4], rep(1, 4))
  expect_equal(df$year[1:4], 2020:2023)
  expect_equal(df$area[5:8], rep(2, 4))
  expect_equal(df$year[5:8], 2020:2023)

  # Check time-major sorting option
  sim_tm <- sim_series_data(D = 20, T = 3, time_start = 2021, sort_order = "time-major", seed = 456)
  df_tm <- sim_tm$data
  expect_equal(df_tm$year[1:20], rep(2021, 20))
  expect_equal(df_tm$area[1:20], 1:20)
})

test_that("sim_series_data handles unsampled domains and distribution bounds", {
  sim <- sim_series_data(D = 30, T = 5, n_unsampled = 4, prop_intermittent = 0.05, seed = 789)
  df <- sim$data

  # Check that at least 4*5 = 20 entries are NA due to persistent unsampled domains
  na_count <- sum(is.na(df$y_gaussian))
  expect_true(na_count >= 20)

  # Auxiliary variables should not be NA
  expect_false(any(is.na(df$vardir)))
  expect_false(any(is.na(df$exposure)))
  expect_false(any(is.na(df$trials)))

  # Beta bounds: strictly between 0 and 1
  sampled_beta <- stats::na.omit(df$y_beta)
  expect_true(all(sampled_beta > 0 & sampled_beta < 1))

  # Poisson bounds: non-negative integers
  sampled_pois <- stats::na.omit(df$y_poisson)
  expect_true(all(sampled_pois >= 0))

  # Binomial bounds: between 0 and trials
  sampled_bin <- stats::na.omit(df$y_binomial)
  sampled_trials <- df$trials[!is.na(df$y_binomial)]
  expect_true(all(sampled_bin >= 0 & sampled_bin <= sampled_trials))

  # Print method
  expect_output(print(sim), "fastsae Simulated Spatio-Temporal Multi-Distribution Panel")
})

test_that("sim_series_data integrates directly with eblup_stfh", {
  # Simulate panel with no missing data for numerical stability in Fisher scoring
  sim <- sim_series_data(D = 25, T = 4, rho_s = 0.4, rho_t = 0.4,
                         n_unsampled = 0, prop_intermittent = 0, seed = 101)

  fit_stfh <- eblup_stfh(
    y_gaussian ~ x1 + x2,
    data = sim$data,
    domain = ~area,
    time = ~year,
    vardir = ~vardir,
    W = sim$W_std,
    print_result = FALSE
  )

  expect_s3_class(fit_stfh, "fastsae")
  expect_true(fit_stfh$convergence)
  expect_equal(nrow(fit_stfh$df_eblup), 25 * 4)
  expect_true(!is.null(fit_stfh$estvarcomp))
})

test_that("built-in sim_panel dataset loads and integrates with eblup_stfh", {
  data("sim_panel", package = "fastsae")
  data("mys_proxmat", package = "fastsae")

  expect_equal(nrow(sim_panel), 42 * 5)
  expect_true("area" %in% names(sim_panel))
  expect_true("year" %in% names(sim_panel))
  expect_true("y_gaussian" %in% names(sim_panel))
  expect_true("y_poisson" %in% names(sim_panel))
  expect_true("y_binomial" %in% names(sim_panel))
  expect_true("y_beta" %in% names(sim_panel))

  # Test eblup_stfh fit on complete domains of sim_panel
  na_areas <- unique(sim_panel$area[is.na(sim_panel$y_gaussian)])
  df_complete <- sim_panel[!sim_panel$area %in% na_areas, ]
  W_complete <- mys_proxmat[!1:42 %in% na_areas, !1:42 %in% na_areas]
  rs <- rowSums(W_complete)
  rs[rs == 0] <- 1
  W_complete <- W_complete / rs

  fit_panel <- eblup_stfh(
    y_gaussian ~ x1 + x2,
    data = df_complete,
    domain = ~area,
    time = ~year,
    vardir = ~vardir,
    W = W_complete,
    print_result = FALSE
  )

  expect_s3_class(fit_panel, "fastsae")
  expect_equal(nrow(fit_panel$df_eblup), nrow(df_complete))
})
