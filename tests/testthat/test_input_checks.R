# Input checks, convergence warnings and mse_type shared by the EBLUP models

sampled_rows <- which(!is.na(mys$y))

twofold_data <- function(m = 15, seed = 1) {
  set.seed(seed)
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:4, 1)
    data.frame(
      area = d, subarea = paste0(d, "-", seq_len(nd)), x1 = rnorm(nd),
      vardir = runif(nd, 0.2, 0.8)
    )
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
  expect_error(
    eblup_fh(y ~ x1, vardir = "vardir", data = d, print_result = FALSE),
    "finite and strictly positive"
  )
  expect_error(
    eblup_sfh(y ~ x1, vardir = "vardir", data = d, W = mys_proxmat, print_result = FALSE),
    "finite and strictly positive"
  )
  dat <- twofold_data()
  dat$vardir[3] <- NA
  expect_error(
    eblup_twofold(y ~ x1,
      vardir = "vardir", domain = "area", data = dat,
      print_result = FALSE
    ),
    "finite and strictly positive"
  )
})

test_that("linearly dependent covariates are rejected with the aliased column", {
  d <- mys
  d$k <- 2
  expect_error(
    eblup_fh(y ~ x1 + k, vardir = "vardir", data = d, print_result = FALSE),
    "linearly dependent"
  )
  expect_error(
    eblup_sfh(y ~ x1 + k, vardir = "vardir", data = d, W = mys_proxmat, print_result = FALSE),
    "linearly dependent"
  )
  p <- mys_panel
  p$k <- 1
  expect_error(
    eblup_stfh(y ~ x1 + k,
      data = p, vardir = ~vardir, domain = ~area, time = ~year,
      W = mys_proxmat, mse_type = "none", print_result = FALSE
    ),
    "linearly dependent"
  )
})

test_that("too few sampled areas for the number of coefficients is an error", {
  d <- mys[sampled_rows[1:2], ]
  expect_error(eblup_fh(y ~ x1, vardir = "vardir", data = d, print_result = FALSE), "Too few")
  one_area <- data.frame(area = 1, x1 = 1:6, vardir = 0.5, y = c(2, 1, 4, 3, 6, 5))
  expect_error(
    eblup_twofold(y ~ x1,
      vardir = "vardir", domain = "area", data = one_area,
      print_result = FALSE
    ),
    "At least 2 areas"
  )
  b <- bhf_data()
  big <- names(which.max(table(b$units$CountyIndex)))
  one <- b$units[b$units$CountyIndex == big, ]
  expect_error(
    eblup_bhf(CornHec ~ CornPix + SoyBeansPix,
      unit_data = one, Xpop = b$pop,
      domain_var = "CountyIndex", popsize_var = "PopnSegments",
      mse_type = "none", print_result = FALSE
    ),
    "At least 2 sampled domains"
  )
})

# ---- Spatial weights --------------------------------------------------------

test_that("W with zero rows or rows not summing to 1 gives a warning", {
  W <- as.matrix(mys_proxmat)
  W0 <- W
  W0[1, ] <- 0
  expect_warning(
    eblup_sfh(y ~ x1, vardir = "vardir", data = mys, W = W0, print_result = FALSE),
    "no neighbours"
  )
  # rows summing to 0.5: still a valid SAR model for |rho| < 1, but not the usual one
  expect_warning(
    eblup_sfh(y ~ x1, vardir = "vardir", data = mys, W = W / 2, print_result = FALSE),
    "not row-standardised"
  )
  expect_warning(
    eblup_stfh(y ~ x1,
      data = mys_panel, vardir = ~vardir, domain = ~area, time = ~year,
      W = W / 2, mse_type = "none", print_result = FALSE
    ),
    "not row-standardised"
  )
  Wn <- W
  Wn[1, 2] <- -0.1
  expect_error(
    eblup_sfh(y ~ x1, vardir = "vardir", data = mys, W = Wn, print_result = FALSE),
    "non-negative"
  )
  expect_no_warning(eblup_sfh(y ~ x1, vardir = "vardir", data = mys, W = W, print_result = FALSE))
})

# ---- Convergence -------------------------------------------------------------

test_that("non-convergence is a warning, not only a message", {
  expect_warning(
    fit <- eblup_fh(y ~ x1 + x2 + x3,
      vardir = "vardir", data = mys, maxiter = 1,
      print_result = FALSE
    ),
    "not converged"
  )
  expect_false(fit$convergence)
})

test_that("INLA convergence is read from the fit", {
  expect_true(.inla_convergence(list(ok = TRUE, mode = list(mode.status = 0))))
  expect_true(.inla_convergence(list()))
  expect_warning(
    ok <- .inla_convergence(list(ok = TRUE, mode = list(mode.status = 1))),
    "mode status 1"
  )
  expect_false(ok)
  expect_warning(ok <- .inla_convergence(list(ok = FALSE)), "posterior mode")
  expect_false(ok)
})

# ---- mse_type ----------------------------------------------------------------

test_that("mse_type = 'none' returns NA MSE in every model", {
  all_na <- function(fit) all(is.na(fit$df_eblup$mse)) && all(is.na(fit$df_eblup$rse))
  f1 <- eblup_fh(y ~ x1, vardir = "vardir", data = mys, mse_type = "none", print_result = FALSE)
  f2 <- eblup_sfh(y ~ x1,
    vardir = "vardir", data = mys, W = mys_proxmat, mse_type = "none",
    print_result = FALSE
  )
  f3 <- eblup_stfh(y ~ x1,
    data = mys_panel, vardir = ~vardir, domain = ~area, time = ~year,
    W = mys_proxmat, mse_type = "none", print_result = FALSE
  )
  f4 <- eblup_twofold(y ~ x1,
    vardir = "vardir", domain = "area", data = twofold_data(),
    mse_type = "none", print_result = FALSE
  )
  b <- bhf_data()
  f5 <- eblup_bhf(CornHec ~ CornPix + SoyBeansPix,
    unit_data = b$units, Xpop = b$pop,
    domain_var = "CountyIndex", popsize_var = "PopnSegments",
    mse_type = "none", print_result = FALSE
  )
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
  f3 <- eblup_stfh(y ~ x1,
    data = mys_panel, vardir = ~vardir, domain = ~area, time = ~year,
    W = mys_proxmat, print_result = FALSE
  )
  expect_identical(f3$mse_type, "analytical")
  expect_false(anyNA(f3$df_eblup$mse))
  b <- bhf_data()
  f5 <- eblup_bhf(CornHec ~ CornPix + SoyBeansPix,
    unit_data = b$units, Xpop = b$pop,
    domain_var = "CountyIndex", popsize_var = "PopnSegments", B = 5, seed = 1,
    print_result = FALSE
  )
  expect_identical(f5$mse_type, "bootstrap")
  expect_true(all(f5$df_eblup$mse[f5$df_eblup$samp_size > 0] > 0))
})

test_that("the arguments of fastsae 0.1.0 still work, with a deprecation warning", {
  lifecycle::expect_deprecated(
    a <- eblup_sfh(y ~ x1,
      vardir = "vardir", data = mys, W = mys_proxmat,
      mse_method = "analytical", print_result = FALSE
    )
  )
  expect_identical(a$mse_type, "analytical")

  args <- list(
    formula = y ~ x1, data = mys_panel, vardir = ~vardir, domain = ~area,
    time = ~year, W = mys_proxmat, print_result = FALSE
  )
  lifecycle::expect_deprecated(s0 <- do.call(eblup_stfh, c(args, compute_mse = FALSE)))
  expect_identical(s0$mse_type, "none")
  lifecycle::expect_deprecated(s1 <- do.call(eblup_stfh, c(args, compute_mse = TRUE)))
  expect_identical(s1$mse_type, "analytical")
  lifecycle::expect_deprecated(s2 <- do.call(eblup_stfh, c(args, compute_mse = TRUE, B = 3, seed = 1)))
  expect_identical(s2$mse_type, "bootstrap")

  b <- bhf_data()
  lifecycle::expect_deprecated(
    u <- eblup_bhf(CornHec ~ CornPix + SoyBeansPix,
      unit_data = b$units, Xpop = b$pop,
      domain_var = "CountyIndex", popsize_var = "PopnSegments",
      compute_mse = FALSE, print_result = FALSE
    )
  )
  expect_identical(u$mse_type, "none")
})

# ---- Hierarchical Bayes inputs ------------------------------------------------
# INLA treats a missing covariate as 0 and the twofold model replaced a missing
# vardir by 1, so these inputs gave wrong estimates without any message.

test_that("missing covariates are an error in the HB models", {
  skip_if_not_installed("INLA")
  d <- mys
  d$x1[which(is.na(d$y))[1]] <- NA
  expect_error(
    hb_area(y ~ x1, data = d, vardir = "vardir", print_result = FALSE),
    "Auxiliary variables contain NA"
  )
  dat <- twofold_data()
  dat$x1[2] <- NA
  expect_error(
    hb_twofold(y ~ x1,
      vardir = "vardir", domain = "area", data = dat,
      print_result = FALSE
    ),
    "Auxiliary variables contain NA"
  )
  b <- bhf_data()
  b$pop$CornPix[1] <- NA
  expect_error(
    hb_unit(CornHec ~ CornPix + SoyBeansPix,
      unit_data = b$units, Xpop = b$pop,
      domain_var = "CountyIndex", popsize_var = "PopnSegments",
      print_result = FALSE
    ),
    "Xpop"
  )
})

test_that("a missing population mean is a clear error in eblup_bhf", {
  b <- bhf_data()
  b$pop$CornPix[1] <- NA
  expect_error(
    eblup_bhf(CornHec ~ CornPix + SoyBeansPix,
      unit_data = b$units, Xpop = b$pop,
      domain_var = "CountyIndex", popsize_var = "PopnSegments",
      mse_type = "none", print_result = FALSE
    ),
    "Xpop"
  )
})

test_that("an invalid vardir of a sampled area is an error in the HB models", {
  skip_if_not_installed("INLA")
  d <- mys
  d$vardir[sampled_rows[2]] <- NA
  expect_error(
    hb_area(y ~ x1, data = d, vardir = "vardir", print_result = FALSE),
    "finite and strictly positive"
  )
  dat <- twofold_data()
  for (v in list(NA, 0, -0.5)) {
    bad <- dat
    bad$vardir[2] <- v
    expect_error(
      hb_twofold(y ~ x1,
        vardir = "vardir", domain = "area", data = bad,
        print_result = FALSE
      ),
      "finite and strictly positive"
    )
  }
})

test_that("linearly dependent covariates are rejected in hb_area", {
  skip_if_not_installed("INLA")
  d <- mys
  d$k <- 2
  expect_error(
    hb_area(y ~ x1 + k, data = d, vardir = "vardir", print_result = FALSE),
    "linearly dependent"
  )
})

# ---- Domain identifiers -------------------------------------------------------

test_that("duplicated domain identifiers give a warning", {
  d <- mys
  d$dom <- rep(seq_len(nrow(d) / 2), each = 2)
  expect_warning(
    eblup_fh(y ~ x1, vardir = "vardir", domain = "dom", data = d, print_result = FALSE),
    "duplicated"
  )
  expect_warning(
    eblup_sfh(y ~ x1,
      vardir = "vardir", domain = "dom", data = d, W = mys_proxmat,
      print_result = FALSE
    ),
    "duplicated"
  )
  skip_if_not_installed("INLA")
  expect_warning(
    hb_area(y ~ x1, data = d, domain = "dom", vardir = "vardir", print_result = FALSE),
    "duplicated"
  )
})

# ---- Seeds --------------------------------------------------------------------

test_that("the seed argument does not change the random number stream of the session", {
  check_rng <- function(fit_call) {
    set.seed(11)
    before <- stats::runif(3)
    set.seed(11)
    fit_call()
    expect_identical(stats::runif(3), before)
  }
  check_rng(function() {
    eblup_sfh(y ~ x1,
      vardir = "vardir", data = mys, W = mys_proxmat,
      mse_type = "pbmse", B = 3, seed = 5, print_result = FALSE
    )
  })
  check_rng(function() {
    eblup_stfh(y ~ x1,
      data = mys_panel, vardir = ~vardir, domain = ~area,
      time = ~year, W = mys_proxmat, mse_type = "bootstrap", B = 2,
      seed = 5, print_result = FALSE
    )
  })
  tf <- twofold_data()
  check_rng(function() {
    eblup_twofold(y ~ x1,
      vardir = "vardir", domain = "area",
      data = tf, mse_type = "bootstrap", B = 3,
      seed = 5, print_result = FALSE
    )
  })
  b <- bhf_data()
  check_rng(function() {
    eblup_bhf(CornHec ~ CornPix + SoyBeansPix,
      unit_data = b$units,
      Xpop = b$pop, domain_var = "CountyIndex",
      popsize_var = "PopnSegments", B = 3, seed = 5,
      print_result = FALSE
    )
  })
  # a fixed seed still gives the same bootstrap MSE
  f <- function() {
    eblup_twofold(y ~ x1,
      vardir = "vardir", domain = "area", data = tf,
      mse_type = "bootstrap", B = 3, seed = 5, print_result = FALSE
    )
  }
  expect_identical(f()$df_eblup$mse, f()$df_eblup$mse)
})

test_that("eblup_twofold uses the bootstrap defaults of the other EBLUP functions", {
  fm <- formals(eblup_twofold)
  expect_identical(fm$B, formals(eblup_sfh)$B)
  expect_identical(fm$seed, formals(eblup_sfh)$seed)
})

# ---- INLA reproducibility --------------------------------------------------------
# hb_unit() with INLA's default "compact" mode returned different county means and
# SDs on identical calls (cornsoybean), up to 3 ha and 13% in the SDs; "classic"
# matches the exact posterior and differs by < 0.1% (means) and < 2% (SDs). Posterior draws used seed 0 (random), so aggregates changed too.

test_that("hb_unit() fits in INLA's classic mode and is reproducible", {
  skip_if_not_installed("INLA")
  b <- bhf_data()
  fit <- function() {
    hb_unit(CornHec ~ CornPix + SoyBeansPix,
      unit_data = b$units, Xpop = b$pop,
      domain_var = "CountyIndex", popsize_var = "PopnSegments",
      print_result = FALSE
    )
  }
  f1 <- fit()
  expect_identical(f1$fit$.args$inla.mode, "classic")
  num <- function(f) as.matrix(f$df_hb[c("hb", "sd")])
  for (i in 1:3) expect_equal(num(fit()), num(f1), tolerance = 0.03)
  # the mode can still be chosen
  f2 <- hb_unit(CornHec ~ CornPix + SoyBeansPix,
    unit_data = b$units, Xpop = b$pop,
    domain_var = "CountyIndex", popsize_var = "PopnSegments",
    engine = inla_engine(inla.mode = "compact"), print_result = FALSE
  )
  expect_identical(f2$fit$.args$inla.mode, "compact")
})

test_that("posterior draws of the INLA engine use a fixed seed by default", {
  expect_gt(inla_engine()$seed, 0)
  skip_if_not_installed("INLA")
  # Gaussian hb_twofold() aggregates by lincombs and stores no configurations
  fit <- hb_twofold(y ~ x1,
    vardir = "vardir", domain = "area", data = twofold_data(),
    engine = inla_engine(control.compute = list(config = TRUE)),
    print_result = FALSE
  )
  draw <- function() .inla_draws(fit$fit, 5, fit$engine_spec$seed)[[3]]$latent
  first <- draw()
  set.seed(11)
  before <- stats::runif(2)
  set.seed(11)
  # a seed given only to INLA::inla.posterior.sample() did not fix the draws
  for (i in 1:3) expect_identical(draw(), first)
  expect_identical(stats::runif(2), before)
})

# hb_twofold() in INLA's compact mode overstated sub-area SDs by up to 10% and
# moved posterior means by up to 0.1 against the exact posterior (numerical
# integration over the two variances); classic matched it.
test_that("hb_twofold() fits in INLA's classic mode and is reproducible", {
  skip_if_not_installed("INLA")
  dat <- twofold_data()
  fit <- function(...) {
    hb_twofold(y ~ x1,
      vardir = "vardir", domain = "area", data = dat,
      print_result = FALSE, ...
    )
  }
  f1 <- fit()
  expect_identical(f1$fit$.args$inla.mode, "classic")
  for (i in 1:2) expect_equal(fit()$df_hb[c("hb", "sd")], f1$df_hb[c("hb", "sd")], tolerance = 0.03)
  expect_identical(fit(engine = inla_engine(inla.mode = "compact"))$fit$.args$inla.mode, "compact")
  dat$cnt <- rpois(nrow(dat), 5)
  p <- hb_twofold(cnt ~ x1, domain = "area", data = dat, family = "poisson", print_result = FALSE)
  expect_identical(p$fit$.args$inla.mode, "compact")
})

test_that("missing or negative aggregation weights of hb_twofold() give a warning", {
  skip_if_not_installed("INLA")
  dat <- twofold_data()
  dat$pop <- 100
  dat$pop[1] <- NA
  dat$pop[2] <- -3
  expect_warning(
    hb_twofold(y ~ x1,
      vardir = "vardir", domain = "area", weight = "pop",
      data = dat, print_result = FALSE
    ),
    "2 subareas"
  )
})

test_that("compare_sae() does not call a smaller posterior variance more efficient", {
  skip_if_not_installed("INLA")
  dat <- twofold_data()
  e <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", data = dat, print_result = FALSE)
  h <- hb_twofold(y ~ x1, vardir = "vardir", domain = "area", data = dat, print_result = FALSE)
  out <- capture.output(print(compare_sae(e, h)), type = "message")
  out <- c(out, capture.output(print(compare_sae(e, h))))
  expect_false(any(grepl("more efficient", out)))
  expect_true(any(grepl("not the same quantity", out)))
  # two EBLUPs keep the efficiency wording
  e2 <- eblup_fh(y ~ x1, vardir = "vardir", data = dat, print_result = FALSE)
  expect_true(any(grepl("more efficient", compare_sae(e2, e)$metrics$Metric)))
})

# hb_unit() predicts the domain mean of a binomial or Poisson model at the mean
# covariates, inv_link(xbar' beta + u_d), which is biased for the mean of
# inv_link(x' beta + u_d) over the units: -8% (binomial) and -16% (Poisson) with
# 78-90% coverage in a simulation. The estimator is kept; the user is told.
test_that("hb_unit() warns that non-Gaussian domain means are evaluated at mean covariates", {
  skip_if_not_installed("INLA")
  set.seed(3)
  D <- 10
  units <- data.frame(dom = rep(seq_len(D), each = 8), x = rnorm(8 * D))
  units$y <- rbinom(nrow(units), 1, plogis(-0.3 + units$x))
  pop <- data.frame(dom = seq_len(D), x = tapply(units$x, units$dom, mean), N = 200)
  expect_warning(
    hb_unit(y ~ x,
      unit_data = units, Xpop = pop, domain_var = "dom",
      popsize_var = "N", family = "binomial", print_result = FALSE
    ),
    "mean covariates"
  )
  units$y <- rpois(nrow(units), exp(0.2 + 0.5 * units$x))
  expect_warning(
    hb_unit(y ~ x,
      unit_data = units, Xpop = pop, domain_var = "dom",
      popsize_var = "N", family = "poisson", print_result = FALSE
    ),
    "mean covariates"
  )
  b <- bhf_data()
  expect_no_warning(hb_unit(CornHec ~ CornPix + SoyBeansPix,
    unit_data = b$units, Xpop = b$pop,
    domain_var = "CountyIndex", popsize_var = "PopnSegments",
    print_result = FALSE
  ))
})

# hb_area() in INLA's compact mode overstated the posterior SDs of Gaussian
# models by up to 18% against Stan when the area variance is small relative to
# vardir; classic matched Stan at the same speed. For the other families both
# modes matched Stan and classic was much slower (beta, 2000 areas: 32 s
# against 1.3 s), so they keep INLA's default compact mode.
test_that("hb_area() fits Gaussian models in INLA's classic mode, other families in compact", {
  skip_if_not_installed("INLA")
  f <- hb_area(y ~ x1, data = mys, vardir = "vardir", print_result = FALSE)
  expect_identical(f$fit$.args$inla.mode, "classic")
  g <- hb_area(y ~ x1,
    data = mys, vardir = "vardir", engine = inla_engine(inla.mode = "compact"),
    print_result = FALSE
  )
  expect_identical(g$fit$.args$inla.mode, "compact")
  b <- hb_area(y_binomial ~ x1,
    data = sim_area, trials = "trials", family = "binomial",
    print_result = FALSE
  )
  expect_identical(b$fit$.args$inla.mode, "compact")
})

# A CAR model uses only the pattern W > 0: a dense matrix such as distance-decay
# weights (all 861 pairs linked) becomes a complete graph.
test_that("a dense W for a CAR model gives a warning", {
  skip_if_not_installed("INLA")
  xy <- cbind(sim_area$x_coord, sim_area$y_coord)
  W_dense <- 1 / as.matrix(stats::dist(xy))
  diag(W_dense) <- 0
  W_dense <- W_dense / rowSums(W_dense)
  expect_warning(
    hb_area(y ~ x1,
      data = mys, vardir = "vardir", spatial = "besag", W = W_dense,
      print_result = FALSE
    ),
    "pattern"
  )
  # mys_proxmat is a 4-nearest-neighbour graph (it was distance-decay weights)
  expect_no_warning(hb_area(y ~ x1,
    data = mys, vardir = "vardir", spatial = "besag", W = mys_proxmat,
    print_result = FALSE
  ))
  W_knn <- create_weights(xy, method = "knn", k = 4, style = "B")
  expect_no_warning(hb_area(y ~ x1,
    data = mys, vardir = "vardir", spatial = "besag", W = W_knn,
    print_result = FALSE
  ))
  # the SAR (slm) model uses the weights themselves
  expect_no_warning(.convert_spatial_weights(W_dense, n_domains = 42, spatial = "slm"))
})

# Knorr-Held (2000) interactions with a random walk in time are intrinsic in
# the direction "one level (rw1) or one linear trend (rw2) per domain". INLA's
# default constraints left those directions with a flat prior: a domain without
# direct estimates got a posterior SD of 42-155 (others about 0.4).
st_panel <- function() {
  D <- 25
  Tt <- 6
  set.seed(1)
  d <- expand.grid(area = 1:D, year = 1:Tt)
  d$x1 <- rnorm(nrow(d))
  d$vardir <- runif(nrow(d), 0.3, 1)
  d$y <- 5 + d$x1 + rnorm(D, sd = 0.5)[d$area] + cumsum(rnorm(Tt, sd = 0.3))[d$year] +
    rnorm(nrow(d), sd = 0.3) + rnorm(nrow(d), sd = sqrt(d$vardir))
  d$y[d$area == 7] <- NA
  list(d = d, W = sim_spatial_weights(D, type = "grid", seed = 1))
}

test_that("type IV and type II interactions are identified by Knorr-Held constraints", {
  skip_if_not_installed("INLA")
  p <- st_panel()
  for (cfg in list(c("type4", "rw1"), c("type4", "rw2"), c("type2", "rw2"))) {
    f <- suppressWarnings(hb_area(y ~ x1,
      data = p$d, vardir = "vardir", domain = "area", time = "year",
      spatial = "besag", temporal = cfg[2], st_interaction = cfg[1],
      W = p$W, print_result = FALSE
    ))
    sd7 <- mean(f$df_hb$sd[f$df_hb$domain == 7])
    expect_lt(sd7, 3 * stats::median(f$df_hb$sd[f$df_hb$domain != 7]))
  }
})

test_that("separable random walks and domain-specific rw2 are rejected", {
  p <- st_panel()
  for (tp in c("rw1", "rw2")) {
    expect_error(
      hb_area(y ~ x1,
        data = p$d, vardir = "vardir", domain = "area", time = "year",
        spatial = "besag", temporal = tp, st_interaction = "separable", W = p$W,
        print_result = FALSE
      ),
      "ar1"
    )
  }
  expect_error(
    hb_area(y ~ x1,
      data = p$d, vardir = "vardir", domain = "area", time = "year",
      temporal = "rw2", st_interaction = "domain-specific", print_result = FALSE
    ),
    "type2"
  )
})

test_that("the Knorr-Held constraints span the null space of the interaction", {
  skip_if_not_installed("INLA")
  W <- as.matrix(unclass(sim_spatial_weights(9, type = "grid", seed = 1)))
  W[9, ] <- W[, 9] <- 0 # an island
  for (cfg in list(c("type4", "rw1"), c("type4", "rw2"), c("type2", "rw2"))) {
    g <- .st_generic0(cfg[1], cfg[2],
      domain_id = rep(1:9, 4), time_id = rep(1:4, each = 9),
      n_domains = 9, n_times = 4, graph = (W > 0) * 1
    )
    A <- g$extraconstr$A
    C <- as.matrix(g$Cmatrix)
    # independent constraints, as many as the rank deficiency, and each null
    # direction of the structure matrix is removed by them
    expect_equal(qr(A)$rank, nrow(A))
    expect_equal(nrow(A), g$rankdef)
    expect_equal(ncol(C) - qr(C)$rank, g$rankdef)
    expect_equal(qr(rbind(C, A))$rank, ncol(C))
  }
})

# Against Stan: Leroux and SLM in compact mode moved Poisson posterior means by
# 0.12-0.15 posterior SD (classic: <= 0.03); Gaussian models with a Knorr-Held
# interaction overstated SDs by 6-13% in both modes, and only grid integration
# matched Stan (<= 0.04 SD); the negative binomial model was off in every mode.
test_that("hb_area() picks the INLA settings validated against Stan", {
  skip_if_not_installed("INLA")
  d <- sim_area
  for (sp in c("leroux", "slm")) {
    f <- suppressWarnings(hb_area(y_poisson ~ x1,
      data = d, exposure = "exposure", family = "poisson",
      spatial = sp, W = if (sp == "slm") mys_proxmat else mys_knn,
      print_result = FALSE
    ))
    expect_identical(f$fit$.args$inla.mode, "classic")
  }
  p <- st_panel()
  g <- suppressWarnings(hb_area(y ~ x1,
    data = p$d, vardir = "vardir", domain = "area", time = "year",
    spatial = "besag", temporal = "rw1", st_interaction = "type2", W = p$W,
    print_result = FALSE
  ))
  expect_identical(g$fit$.args$control.inla$int.strategy, "grid")
  expect_identical(g$fit$.args$inla.mode, "classic")
  # without an interaction the default integration is kept
  h <- suppressWarnings(hb_area(y ~ x1,
    data = p$d, vardir = "vardir", domain = "area", time = "year",
    spatial = "besag", temporal = "rw1", W = p$W, print_result = FALSE
  ))
  expect_false(identical(h$fit$.args$control.inla$int.strategy, "grid"))
})

test_that("the hyperparameter grid is used only for the three-hyperparameter models validated against Stan", {
  skip_if_not_installed("INLA")
  set.seed(42)
  d <- expand.grid(time = 1:2, domain = 1:6)
  d$x <- rnorm(12)
  d$y <- 5 + 0.5 * d$x + rnorm(12, 0, 0.2)
  d$vardir <- 0.1
  W <- sim_spatial_weights(6, type = "ring")
  fit <- function(...) {
    suppressWarnings(hb_area(y ~ x,
      data = d, domain = "domain", time = "time",
      print_result = FALSE, ...
    ))
  }
  # AR(1) adds a correlation to the temporal effect and to a type II/IV
  # interaction: a grid over 4 or 5 hyperparameters took 30 s and 13 min here
  cfgs <- list(
    list("none", "ar1", "type1", 4), list("none", "ar1", "type2", 5),
    list("bym2", "rw1", "type3", 4), list("besag", "ar1", "type4", 5)
  )
  for (cfg in cfgs) {
    sp <- cfg[[1]]
    f <- fit(
      vardir = "vardir", spatial = sp, temporal = cfg[[2]], st_interaction = cfg[[3]],
      W = if (sp != "none") W
    )
    expect_equal(nrow(f$fit$summary.hyperpar), cfg[[4]])
    expect_equal(.n_hyper_st(sp, cfg[[2]], cfg[[3]], vardir_known = TRUE), cfg[[4]])
    expect_false(identical(f$fit$.args$control.inla$int.strategy, "grid"))
  }
  # the validated case keeps the grid; an estimated residual variance is a fourth
  g <- fit(vardir = "vardir", spatial = "besag", temporal = "rw1", st_interaction = "type2", W = W)
  expect_equal(nrow(g$fit$summary.hyperpar), 3)
  expect_identical(g$fit$.args$control.inla$int.strategy, "grid")
  h <- fit(spatial = "besag", temporal = "rw1", st_interaction = "type2", W = W)
  expect_equal(nrow(h$fit$summary.hyperpar), .n_hyper_st("besag", "rw1", "type2", vardir_known = FALSE))
  expect_false(identical(h$fit$.args$control.inla$int.strategy, "grid"))
})

test_that("the negative binomial model of hb_area() warns about its accuracy", {
  skip_if_not_installed("INLA")
  expect_warning(
    hb_area(y_nbinomial ~ x1,
      data = sim_area, exposure = "exposure", family = "nbinomial",
      print_result = FALSE
    ),
    "weakly identified"
  )
})

test_that("a random effect variance estimated as 0 is reported with a message", {
  set.seed(1)
  d <- data.frame(x = rnorm(25), v = 4)
  # sampling variances far above the spread of y around the regression line
  d$y <- 1 + d$x + rnorm(25, sd = 0.5)
  W <- mys_proxmat[1:25, 1:25]
  W <- W / rowSums(W)
  expect_message(
    f <- eblup_fh(y ~ x, vardir = "v", data = d, print_result = FALSE),
    "random effect variance is estimated as 0"
  )
  expect_equal(f$random_effect_var, 0)
  expect_message(
    s <- eblup_sfh(y ~ x, vardir = "v", data = d, W = W, print_result = FALSE),
    "random effect variance is estimated as 0"
  )
  expect_equal(s$random_effect_var, 0)
  # a positive estimate gives no such message
  expect_no_message(f <- eblup_fh(y ~ x1 + x2 + x3, vardir = "vardir", data = mys, print_result = FALSE))
  expect_gt(f$random_effect_var, 0)
})
