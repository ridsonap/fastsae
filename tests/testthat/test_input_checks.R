# Input checks, convergence warnings and mse_type shared by the EBLUP models

sampled_rows <- which(!is.na(mys$y))

twofold_data <- function(m = 15, seed = 1) {
  set.seed(seed)
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:4, 1)
    data.frame(area = d, subarea = paste0(d, "-", seq_len(nd)), x1 = rnorm(nd),
               vardir = runif(nd, 0.2, 0.8))
  }))
  dat$y <- 1 + dat$x1 + rnorm(m, sd = 0.8)[dat$area] + rnorm(nrow(dat), sd = 0.6) +
    rnorm(nrow(dat), sd = sqrt(dat$vardir))
  dat
}

bhf_data <- function() {
  pop <- cornsoybeanmeans
  names(pop)[names(pop) == "MeanCornPixPerSeg"] <- "CornPix"
  names(pop)[names(pop) == "MeanSoyBeansPixPerSeg"] <- "SoyBeansPix"
  units <- cornsoybean
  names(units)[names(units) == "County"] <- "CountyIndex"
  list(pop = pop, units = units)
}

# ---- Sampled design and vardir ------------------------------------------------

test_that("a missing vardir of a sampled area gives a clear error", {
  d <- mys
  d$vardir[sampled_rows[2]] <- NA
  expect_error(eblup_fh(y ~ x1, vardir = "vardir", data = d, print_result = FALSE),
               "finite and strictly positive")
  expect_error(eblup_sfh(y ~ x1, vardir = "vardir", data = d, W = mys_proxmat, print_result = FALSE),
               "finite and strictly positive")
  dat <- twofold_data()
  dat$vardir[3] <- NA
  expect_error(eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", data = dat,
                             print_result = FALSE),
               "finite and strictly positive")
})

test_that("linearly dependent covariates are rejected with the aliased column", {
  d <- mys
  d$k <- 2
  expect_error(eblup_fh(y ~ x1 + k, vardir = "vardir", data = d, print_result = FALSE),
               "linearly dependent")
  expect_error(eblup_sfh(y ~ x1 + k, vardir = "vardir", data = d, W = mys_proxmat, print_result = FALSE),
               "linearly dependent")
  p <- mys_panel
  p$k <- 1
  expect_error(eblup_stfh(y ~ x1 + k, data = p, vardir = ~vardir, domain = ~area, time = ~year,
                          W = mys_proxmat, mse_type = "none", print_result = FALSE),
               "linearly dependent")
})

test_that("too few sampled areas for the number of coefficients is an error", {
  d <- mys[sampled_rows[1:2], ]
  expect_error(eblup_fh(y ~ x1, vardir = "vardir", data = d, print_result = FALSE), "Too few")
  one_area <- data.frame(area = 1, x1 = 1:6, vardir = 0.5, y = c(2, 1, 4, 3, 6, 5))
  expect_error(eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", data = one_area,
                             print_result = FALSE),
               "At least 2 areas")
  b <- bhf_data()
  big <- names(which.max(table(b$units$CountyIndex)))
  one <- b$units[b$units$CountyIndex == big, ]
  expect_error(eblup_bhf(CornHec ~ CornPix + SoyBeansPix, unit_data = one, Xpop = b$pop,
                         domain_var = "CountyIndex", popsize_var = "PopnSegments",
                         mse_type = "none", print_result = FALSE),
               "At least 2 sampled domains")
})

# ---- Spatial weights --------------------------------------------------------

test_that("W with zero rows or rows not summing to 1 gives a warning", {
  W <- as.matrix(mys_proxmat)
  W0 <- W
  W0[1, ] <- 0
  expect_warning(eblup_sfh(y ~ x1, vardir = "vardir", data = mys, W = W0, print_result = FALSE),
                 "no neighbours")
  # rows summing to 0.5: still a valid SAR model for |rho| < 1, but not the usual one
  expect_warning(eblup_sfh(y ~ x1, vardir = "vardir", data = mys, W = W / 2, print_result = FALSE),
                 "not row-standardised")
  expect_warning(eblup_stfh(y ~ x1, data = mys_panel, vardir = ~vardir, domain = ~area, time = ~year,
                            W = W / 2, mse_type = "none", print_result = FALSE),
                 "not row-standardised")
  Wn <- W
  Wn[1, 2] <- -0.1
  expect_error(eblup_sfh(y ~ x1, vardir = "vardir", data = mys, W = Wn, print_result = FALSE),
               "non-negative")
  expect_no_warning(eblup_sfh(y ~ x1, vardir = "vardir", data = mys, W = W, print_result = FALSE))
})

# ---- Convergence -------------------------------------------------------------

test_that("non-convergence is a warning, not only a message", {
  expect_warning(
    fit <- eblup_fh(y ~ x1 + x2 + x3, vardir = "vardir", data = mys, maxiter = 1,
                    print_result = FALSE),
    "not converged"
  )
  expect_false(fit$convergence)
})

test_that("INLA convergence is read from the fit", {
  expect_true(.inla_convergence(list(ok = TRUE, mode = list(mode.status = 0))))
  expect_true(.inla_convergence(list()))
  expect_warning(ok <- .inla_convergence(list(ok = TRUE, mode = list(mode.status = 1))),
                 "mode status 1")
  expect_false(ok)
  expect_warning(ok <- .inla_convergence(list(ok = FALSE)), "posterior mode")
  expect_false(ok)
})

# ---- mse_type ----------------------------------------------------------------

test_that("mse_type = 'none' returns NA MSE in every model", {
  all_na <- function(fit) all(is.na(fit$df_eblup$mse)) && all(is.na(fit$df_eblup$rse))
  f1 <- eblup_fh(y ~ x1, vardir = "vardir", data = mys, mse_type = "none", print_result = FALSE)
  f2 <- eblup_sfh(y ~ x1, vardir = "vardir", data = mys, W = mys_proxmat, mse_type = "none",
                  print_result = FALSE)
  f3 <- eblup_stfh(y ~ x1, data = mys_panel, vardir = ~vardir, domain = ~area, time = ~year,
                   W = mys_proxmat, mse_type = "none", print_result = FALSE)
  f4 <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", data = twofold_data(),
                      mse_type = "none", print_result = FALSE)
  b <- bhf_data()
  f5 <- eblup_bhf(CornHec ~ CornPix + SoyBeansPix, unit_data = b$units, Xpop = b$pop,
                  domain_var = "CountyIndex", popsize_var = "PopnSegments",
                  mse_type = "none", print_result = FALSE)
  for (f in list(f1, f2, f3, f4, f5)) {
    expect_true(all_na(f))
    expect_identical(f$mse_type, "none")
  }
  # the point estimates do not depend on mse_type
  f1a <- eblup_fh(y ~ x1, vardir = "vardir", data = mys, print_result = FALSE)
  expect_identical(f1$df_eblup$eblup, f1a$df_eblup$eblup)
  expect_false(anyNA(f1a$df_eblup$mse[!is.na(mys$y)]))
})

test_that("the default mse_type computes an MSE in every model", {
  f3 <- eblup_stfh(y ~ x1, data = mys_panel, vardir = ~vardir, domain = ~area, time = ~year,
                   W = mys_proxmat, print_result = FALSE)
  expect_identical(f3$mse_type, "analytical")
  expect_false(anyNA(f3$df_eblup$mse))
  b <- bhf_data()
  f5 <- eblup_bhf(CornHec ~ CornPix + SoyBeansPix, unit_data = b$units, Xpop = b$pop,
                  domain_var = "CountyIndex", popsize_var = "PopnSegments", B = 5, seed = 1,
                  print_result = FALSE)
  expect_identical(f5$mse_type, "bootstrap")
  expect_true(all(f5$df_eblup$mse[f5$df_eblup$samp_size > 0] > 0))
})

test_that("the arguments of fastsae 0.1.0 still work, with a deprecation warning", {
  lifecycle::expect_deprecated(
    a <- eblup_sfh(y ~ x1, vardir = "vardir", data = mys, W = mys_proxmat,
                   mse_method = "analytical", print_result = FALSE)
  )
  expect_identical(a$mse_type, "analytical")

  args <- list(formula = y ~ x1, data = mys_panel, vardir = ~vardir, domain = ~area,
               time = ~year, W = mys_proxmat, print_result = FALSE)
  lifecycle::expect_deprecated(s0 <- do.call(eblup_stfh, c(args, compute_mse = FALSE)))
  expect_identical(s0$mse_type, "none")
  lifecycle::expect_deprecated(s1 <- do.call(eblup_stfh, c(args, compute_mse = TRUE)))
  expect_identical(s1$mse_type, "analytical")
  lifecycle::expect_deprecated(s2 <- do.call(eblup_stfh, c(args, compute_mse = TRUE, B = 3, seed = 1)))
  expect_identical(s2$mse_type, "bootstrap")

  b <- bhf_data()
  lifecycle::expect_deprecated(
    u <- eblup_bhf(CornHec ~ CornPix + SoyBeansPix, unit_data = b$units, Xpop = b$pop,
                   domain_var = "CountyIndex", popsize_var = "PopnSegments",
                   compute_mse = FALSE, print_result = FALSE)
  )
  expect_identical(u$mse_type, "none")
})
