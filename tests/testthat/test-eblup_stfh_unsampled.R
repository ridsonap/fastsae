# eblup_stfh() with cells without a direct estimate (NA y): intermittent gaps
# and domains that are never sampled. The C++ results are compared with a
# dense reference (helper-stfh-dense.R) that forms V, the BLUP and the MSE
# pieces by direct algebra.

stfh_case <- function(seed, D = 12, T = 5, n_unsampled = 2, prop = 0.15) {
  sim <- sim_series_data(
    D = D, T = T, n_unsampled = n_unsampled, prop_intermittent = prop,
    sigma_s = 0.8, sigma_t = 0.6, rho_s = 0.4, rho_t = 0.5, seed = seed
  )
  d <- sim$data
  d <- d[order(d$area, d$year), ]
  W <- unclass(sim$W)
  W <- W / rowSums(W)
  attributes(W)[c("coords", "type", "style")] <- NULL
  list(d = d, W = W, X = cbind(1, d$x1, d$x2), D = D, T = T)
}

fit_case <- function(z, y = z$d$y_gaussian, ...) {
  d <- z$d
  d$y_gaussian <- y
  eblup_stfh(y_gaussian ~ x1 + x2,
    data = d, vardir = ~vardir, domain = ~area, time = ~year,
    W = z$W, print_result = FALSE, ...
  )
}

test_that("predictions, beta, g1 and g2 of every cell match the dense BLUP", {
  for (model in c("ST", "S")) {
    z <- stfh_case(5)
    m <- fit_case(z, model = model)
    expect_true(m$convergence)
    expect_true(any(!m$df_eblup$sampled))
    ref <- .stfh_dense_predict(
      m$estvarcomp$estimate, z$X, z$d$y_gaussian, z$d$vardir,
      z$W, z$D, z$T, model
    )
    expect_equal(m$df_eblup$eblup, ref$pred, tolerance = 1e-8)
    expect_equal(m$estcoef$beta, ref$beta, tolerance = 1e-8)
    expect_equal(m$g1, ref$g1, tolerance = 1e-8)
    expect_equal(m$g2, ref$g2, tolerance = 1e-8)
  }
})

test_that("the estimates maximise the restricted likelihood of the sampled cells", {
  z <- stfh_case(13)
  m <- fit_case(z)
  th <- m$estvarcomp$estimate
  f <- function(t) .stfh_dense_reml(t, z$X, z$d$y_gaussian, z$d$vardir, z$W, z$D, z$T)
  h <- 1e-5
  grad <- vapply(seq_along(th), function(a) {
    e <- replace(numeric(length(th)), a, h)
    (f(th + e) - f(th - e)) / (2 * h)
  }, numeric(1))
  expect_lt(max(abs(grad)), 1e-3)
})

test_that("g3 (all cells) and g4 (sampled cells) match the dense second-order terms", {
  for (seed in c(5, 13)) {
    z <- stfh_case(seed)
    m <- fit_case(z)
    th <- m$estvarcomp$estimate
    expect_true(all(abs(th[c(2, 4)]) < 0.99) && all(th[c(1, 3)] > 1e-3)) # interior estimates
    ref <- .stfh_dense_mse(th, z$X, z$d$y_gaussian, z$d$vardir, z$W, z$D, z$T)
    s <- m$df_eblup$sampled
    expect_equal(m$g3, ref$g3, tolerance = 1e-6)
    expect_equal(m$g4[s], ref$g4[s], tolerance = 1e-2) # finite-difference second derivatives
    expect_true(all(m$g4[!s] == 0)) # no plug-in correction off-sample
    expect_equal(m$df_eblup$mse[s], (m$g1 + m$g2 + 2 * m$g3 - m$g4)[s])
    expect_equal(m$df_eblup$mse[!s], (m$g1 + m$g2 + m$g3)[!s])
  }
})

test_that("unsampled domains, single observed periods and gaps at both ends are predicted", {
  z <- stfh_case(5)
  y <- z$d$y_gaussian
  doms <- unique(z$d$area)
  y[z$d$area == doms[3] & z$d$year != min(z$d$year)] <- NA # one observed period
  y[z$d$area == doms[4] & z$d$year == min(z$d$year)] <- NA # gap at t = 1
  y[z$d$area == doms[5] & z$d$year == max(z$d$year)] <- NA # gap at t = T
  m <- fit_case(z, y = y)
  df <- m$df_eblup
  expect_equal(df$sampled, !is.na(y))
  expect_true(all(is.finite(df$eblup)) && all(df$mse > 0))
  ref <- .stfh_dense_predict(m$estvarcomp$estimate, z$X, y, z$d$vardir, z$W, z$D, z$T)
  expect_equal(df$eblup, ref$pred, tolerance = 1e-8)
  # a domain without any sample has no period effect: u2 = 0 there
  never <- tapply(is.na(y), z$d$area, all)
  expect_true(all(df$random_effect_u2[z$d$area %in% names(never)[never]] == 0))
  expect_true(all(is.na(df$y[!df$sampled])) && all(is.na(df$vardir[!df$sampled])))
})

test_that("results follow the input row order", {
  z <- stfh_case(8)
  m <- fit_case(z)
  o <- order(z$d$year, match(z$d$area, unique(z$d$area))) # time-major input
  d2 <- z$d[o, ]
  m2 <- eblup_stfh(y_gaussian ~ x1 + x2,
    data = d2, vardir = ~vardir, domain = ~area,
    time = ~year, W = z$W, mse_type = "analytical", print_result = FALSE
  )
  expect_equal(m2$df_eblup$eblup, m$df_eblup$eblup[o])
  expect_equal(m2$df_eblup$mse, m$df_eblup$mse[o])
  expect_equal(m2$df_eblup$sampled, m$df_eblup$sampled[o])
})

test_that("the bootstrap MSE covers every cell and does not depend on the threads", {
  z <- stfh_case(5)
  b1 <- fit_case(z, mse_type = "bootstrap", B = 40, seed = 3, n_threads = 1)
  b2 <- fit_case(z, mse_type = "bootstrap", B = 40, seed = 3, n_threads = 2)
  expect_identical(b1$df_eblup$mse, b2$df_eblup$mse)
  expect_true(all(is.finite(b1$df_eblup$mse)) && all(b1$df_eblup$mse > 0))
  expect_equal(b1$B, 40L)
  an <- fit_case(z)
  s <- an$df_eblup$sampled
  # same order of magnitude as the analytical MSE in both kinds of cells
  expect_true(abs(log(median(b1$df_eblup$mse[!s] / an$df_eblup$mse[!s]))) < log(2))
  expect_true(abs(log(median(b1$df_eblup$mse[s] / an$df_eblup$mse[s]))) < log(2))
})

test_that("mys_panel is fitted with its unsampled domains", {
  m <- eblup_stfh(y ~ x1 + x2 + x3,
    data = mys_panel, vardir = ~vardir, domain = ~area,
    time = ~year, W = mys_proxmat, mse_type = "analytical", print_result = FALSE
  )
  df <- m$df_eblup
  expect_equal(nrow(df), nrow(mys_panel))
  expect_equal(sum(!df$sampled), sum(is.na(mys_panel$y)))
  expect_true(all(is.finite(df$mse)))
  expect_gt(min(df$mse[!df$sampled]), median(df$mse[df$sampled]))
  expect_message(print(m), "without a direct estimate")
  expect_s3_class(diagnose(m), "fastsae_diagnose")
  flags <- add_reliability_flags(m)
  expect_equal(nrow(flags), nrow(mys_panel))
})

test_that("inputs that cannot be predicted are rejected", {
  z <- stfh_case(5)
  d <- z$d
  d$x1[1] <- NA
  expect_error(
    eblup_stfh(y_gaussian ~ x1 + x2,
      data = d, vardir = ~vardir, domain = ~area, time = ~year,
      W = z$W, print_result = FALSE
    ),
    "auxiliary variables"
  )
  d <- z$d
  d$vardir[which(!is.na(d$y_gaussian))[1]] <- NA
  expect_error(
    eblup_stfh(y_gaussian ~ x1 + x2,
      data = d, vardir = ~vardir, domain = ~area, time = ~year,
      W = z$W, print_result = FALSE
    ),
    "vardir"
  )
  d <- z$d
  d$y_gaussian <- NA_real_
  expect_error(
    eblup_stfh(y_gaussian ~ x1 + x2,
      data = d, vardir = ~vardir, domain = ~area, time = ~year,
      W = z$W, print_result = FALSE
    ),
    "every cell"
  )
  expect_error(
    eblup_stfh(y_gaussian ~ x1 + x2,
      data = z$d, vardir = ~vardir, domain = ~area, time = ~year,
      W = z$W[-1, -1], print_result = FALSE
    ),
    "square matrix"
  )
})
