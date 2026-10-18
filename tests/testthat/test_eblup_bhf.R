suppressMessages({
  library(testthat)
  library(fastsae)
  library(dplyr)
  library(sae)
})

skip_if_not_installed("sae")

# ============================================================================
# Helper
# ============================================================================
quiet <- function(expr) {
  result <- NULL
  invisible(capture.output(result <- suppressWarnings(suppressMessages(expr))))
  result
}

# ============================================================================
# Load data
# ============================================================================
data(cornsoybean, package = "sae")
data(cornsoybeanmeans, package = "sae")

# Attach formula variables to environment
CornHec <- cornsoybean$CornHec
CornPix <- cornsoybean$CornPix
SoyBeansPix <- cornsoybean$SoyBeansPix

# Prepare data for fastsae
df_cornsoybean <- cornsoybean
df_cornsoybean$CountyIndex <- df_cornsoybean$County
df_cornsoybean$County <- NULL

# Prepare Xmean from cornsoybeanmeans
Xmean <- cornsoybeanmeans[, c("CountyIndex", "MeanCornPixPerSeg", "MeanSoyBeansPixPerSeg")]
names(Xmean) <- c("CountyIndex", "CornPix", "SoyBeansPix")

# Prepare population sizes
Popn <- cornsoybeanmeans[, c("CountyIndex", "PopnSegments")]

# Merge Xmean with Popn
Xpop <- merge(Xmean, Popn, by = "CountyIndex")

# ============================================================================
# Fit using fastsae
# ============================================================================
fit_fast <- eblup_bhf(
  formula = CornHec ~ CornPix + SoyBeansPix,
  unit_data = df_cornsoybean,
  Xpop = Xpop,
  domain_var = "CountyIndex",
  popsize_var = "PopnSegments",
  method = "REML",
  print_result = FALSE
)

# ============================================================================
# Fit using sae package
# ============================================================================
# Prepare meanxpop and popnsize for sae
Xmean_sae <- data.frame(
  CountyIndex = cornsoybeanmeans$CountyIndex,
  MeanCornPixPerSeg = cornsoybeanmeans$MeanCornPixPerSeg,
  MeanSoyBeansPixPerSeg = cornsoybeanmeans$MeanSoyBeansPixPerSeg
)
Popn_sae <- data.frame(
  CountyIndex = cornsoybeanmeans$CountyIndex,
  PopnSegments = cornsoybeanmeans$PopnSegments
)

fit_sae <- quiet(sae::eblupBHF(
  CornHec ~ CornPix + SoyBeansPix,
  dom = cornsoybean$County,
  selectdom = unique(cornsoybean$County),
  meanxpop = Xmean_sae,
  popnsize = Popn_sae,
  method = "REML"
))

tol <- 1e-4

# ============================================================================
# Tests: EBLUP
# ============================================================================
test_that("eblup_bhf returns valid structure", {
  expect_s3_class(fit_fast, "fastsae_unit")
  expect_true("df_eblup" %in% names(fit_fast))
  expect_true("fit" %in% names(fit_fast))
  expect_true("random_effect_var" %in% names(fit_fast$fit))
  expect_true("sigma2_e" %in% names(fit_fast$fit))
})

test_that("EBLUP agrees with sae::eblupBHF", {
  # Order by domain
  fast_eblup <- fit_fast$df_eblup[order(fit_fast$df_eblup$domain), ]
  sae_eblup <- fit_sae$eblup[order(fit_sae$eblup$domain), ]

  expect_equal(
    fast_eblup$eblup,
    as.numeric(sae_eblup$eblup),
    tolerance = tol
  )
})

test_that("Variance components agree with sae::eblupBHF", {
  expect_equal(
    fit_fast$fit$random_effect_var,
    fit_sae$fit$refvar,
    tolerance = tol
  )
})

test_that("Fixed effects agree with sae::eblupBHF", {
  fast_beta <- as.numeric(fit_fast$fit$beta)
  sae_beta <- as.numeric(fit_sae$fit$fixed)

  expect_equal(
    fast_beta,
    sae_beta,
    tolerance = tol
  )
})

# ============================================================================
# Tests: Bootstrap MSE
# ============================================================================
test_that("pbmse_unit returns valid structure", {
  skip_on_cran()

  B_test <- 50
  mse_fast <- .pbmse_unit(
    formula = CornHec ~ CornPix + SoyBeansPix,
    unit_data = df_cornsoybean,
    Xpop = Xpop,
    domain_var = "CountyIndex",
    popsize_var = "PopnSegments",
    method = "REML",
    B = B_test,
    seed = 123
  )

  expect_type(mse_fast$mse, "double")
  expect_length(mse_fast$mse, length(unique(df_cornsoybean$CountyIndex)))
  expect_true(all(mse_fast$mse >= 0))
})

test_that("pbmse_unit produces valid MSE estimates", {
  skip_on_cran()

  B_test <- 50
  set.seed(123)

  mse_result <- .pbmse_unit(
    formula = CornHec ~ CornPix + SoyBeansPix,
    unit_data = df_cornsoybean,
    Xpop = Xpop,
    domain_var = "CountyIndex",
    popsize_var = "PopnSegments",
    method = "REML",
    B = B_test,
    seed = 123
  )

  expect_type(mse_result$mse, "double")
  expect_length(mse_result$mse, 12)
  expect_true(all(mse_result$mse > 0))

  # MSE should be in reasonable range (based on response scale)
  expect_true(all(mse_result$mse < 1000))
})

# ============================================================================
# Tests: Error handling
# ============================================================================
test_that("Invalid method throws error", {
  expect_error(
    eblup_bhf(
      formula = CornHec ~ CornPix + SoyBeansPix,
      unit_data = df_cornsoybean,
      Xpop = Xpop,
      domain_var = "CountyIndex",
      popsize_var = "PopnSegments",
      method = "INVALID"
    )
  )
})

test_that("eblup_bhf works with popnmean_xpop and handles all options", {
  # 1. popnmean_xpop without intercept (auto-prepended)
  mean_no_int <- as.matrix(Xmean[, c("CornPix", "SoyBeansPix")])
  fit_popmean <- eblup_bhf(
    formula = CornHec ~ CornPix + SoyBeansPix,
    unit_data = df_cornsoybean,
    Xpop = Xpop,
    popnmean_xpop = mean_no_int,
    domain_var = "CountyIndex",
    popsize_var = "PopnSegments",
    method = "ML",
    print_result = TRUE
  )
  expect_s3_class(fit_popmean, "fastsae")

  # 2. popnmean_xpop with wrong number of columns
  expect_error(
    eblup_bhf(
      formula = CornHec ~ CornPix + SoyBeansPix,
      unit_data = df_cornsoybean,
      Xpop = Xpop,
      popnmean_xpop = matrix(1:60, 12, 5),
      domain_var = "CountyIndex",
      popsize_var = "PopnSegments",
      print_result = FALSE
    ),
    "Number of columns in auxiliary population means"
  )

  # 3. Direct cpp call covering empty domain / warn_domains and error stops
  upred_df <- data.frame(upred = c(0.1, 0.2), row.names = c("A", "B"))
  Xs_test <- matrix(c(1, 1, 2, 3), 2, 2, dimnames = list(NULL, c("(Intercept)", "x")))
  meanx_test <- matrix(c(1, 1, 2, 3), 2, 2)
  res_cpp <- fastsae:::.eblup_bhf_cpp(
    selectdom = c("A", "C"), # C not in dom
    dom = c("A", "B"),
    Xs = Xs_test,
    meanxpop = meanx_test,
    ys = c(10, 20),
    popnsize = c(100, 100),
    betaest = matrix(c(1, 1), 2, 1),
    upred = upred_df
  )
  expect_equal(res_cpp$warn_domains, "C")

  # Error checks in .eblup_bhf_cpp
  expect_error(fastsae:::.eblup_bhf_cpp(c("A", "B"), c("A", "B"), Xs_test, matrix(1:6, 2, 3), c(1, 2), c(10, 10), matrix(1, 2, 1), upred_df), "Number of columns in meanxpop")
  expect_error(fastsae:::.eblup_bhf_cpp(c("A", "B"), c("A", "B"), Xs_test, matrix(1:2, 1, 2), c(1, 2), c(10, 10), matrix(1, 2, 1), upred_df), "Number of rows in meanxpop")
  expect_error(fastsae:::.eblup_bhf_cpp(c("A", "B"), c("A", "B"), Xs_test, meanx_test, c(1, 2), c(10), matrix(1, 2, 1), upred_df), "Length of popnsize")
})

test_that("eblup_bhf bootstrap MSE computation works", {
  fit_mse <- eblup_bhf(
    formula = CornHec ~ CornPix + SoyBeansPix,
    unit_data = df_cornsoybean,
    Xpop = Xpop,
    domain_var = "CountyIndex",
    popsize_var = "PopnSegments",
    method = "REML",
    mse_type = "bootstrap",
    B = 5,
    seed = 123,
    print_result = TRUE
  )

  expect_s3_class(fit_mse, "fastsae_unit")
  expect_true("mse" %in% names(fit_mse$df_eblup))
  expect_true("rse" %in% names(fit_mse$df_eblup))
  expect_true(all(!is.na(fit_mse$df_eblup$mse)))
  expect_true(all(fit_mse$df_eblup$mse > 0))
  expect_true(all(fit_mse$df_eblup$rse > 0))
})

test_that("eblup_bhf handles unsampled domains, full enumeration, and NA in unit data", {
  # 1. NA in unit data (leaves at least 1 obs for county 1)
  df_na <- df_cornsoybean
  df_na$CornPix[nrow(df_na)] <- NA

  # 2. Unsampled domain in Xpop
  Xpop_extra <- rbind(Xpop, data.frame(CountyIndex = 999, PopnSegments = 50, CornPix = 300, SoyBeansPix = 200))

  # 3. Full enumeration domain (nd == Ni)
  Xpop_extra$PopnSegments[1] <- sum(df_na$CountyIndex == Xpop_extra$CountyIndex[1], na.rm = TRUE)

  expect_warning(
    fit_warn <- eblup_bhf(
      formula = CornHec ~ CornPix + SoyBeansPix,
      unit_data = df_na,
      Xpop = Xpop_extra,
      domain_var = "CountyIndex",
      popsize_var = "PopnSegments",
      mse_type = "bootstrap",
      B = 2,
      seed = 42,
      print_result = FALSE
    ),
    "domain\\(s\\) have no sample units"
  )

  expect_s3_class(fit_warn, "fastsae_unit")
  expect_true(999 %in% fit_warn$df_eblup$domain)
  expect_equal(fit_warn$df_eblup$samp_size[fit_warn$df_eblup$domain == 999], 0)
  expect_true(fit_warn$df_eblup$mse[fit_warn$df_eblup$domain == 999] > 0)

  # 4. .pbmse_unit directly with NA in the unit data
  res_na <- fastsae:::.pbmse_unit(
    formula = CornHec ~ CornPix + SoyBeansPix,
    unit_data = df_na,
    Xpop = Xpop,
    domain_var = "CountyIndex",
    popsize_var = "PopnSegments",
    B = 2,
    seed = 42
  )
  expect_true(all(c("eblup", "mse", "B") %in% names(res_na)))
  expect_equal(res_na$B, 2L)
})


# ============================================================================
# Tests: compiled REML/ML fit and bootstrap
# ============================================================================
.lmer_ref <- function(y, x, dom, reml) {
  dd <- data.frame(y = y, x = x, dom = factor(dom))
  lme4::lmer(y ~ x + (1 | dom), dd,
    REML = reml,
    control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(rhoend = 1e-12))
  )
}

test_that("compiled BHF fit agrees with lme4 (REML, ML, and the s2u = 0 boundary)", {
  skip_if_not_installed("lme4")
  set.seed(11)
  D <- 40
  nd <- sample(2:12, D, TRUE)
  dom <- rep(seq_len(D), nd)
  n <- length(dom)
  x <- rnorm(n, 50, 10)
  for (su in c(1, 0)) {
    y <- 100 + 2 * x + rnorm(D, sd = su)[dom] + rnorm(n, sd = 2)
    for (reml in c(TRUE, FALSE)) {
      f <- fastsae:::.bhf_fit(cbind(1, x), y, dom - 1L, D, if (reml) "REML" else "ML")
      m <- suppressMessages(.lmer_ref(y, x, dom, reml))
      expect_equal(as.numeric(f$beta), as.numeric(lme4::fixef(m)), tolerance = 1e-6)
      expect_equal(f$sigma2_e, sigma(m)^2, tolerance = 1e-6)
      expect_equal(f$sigma2_u, lme4::VarCorr(m)$dom[1], tolerance = 1e-4)
      expect_equal(f$vcov, unname(as.matrix(vcov(m))), tolerance = 1e-5)
      if (su == 0) expect_equal(f$sigma2_u, 0)
      # the compiled maximiser is never worse than lme4's
      devf <- lme4::lmer(y ~ x + (1 | dom), data.frame(y, x, dom = factor(dom)), REML = reml, devFunOnly = TRUE)
      expect_lte(devf(sqrt(f$sigma2_u / f$sigma2_e)), devf(lme4::getME(m, "theta")) + 1e-8)
    }
  }
})

test_that("eblup_bhf returns lme4-equivalent coefficients and residuals", {
  skip_if_not_installed("lme4")
  m <- suppressMessages(.lmer_ref(df_cornsoybean$CornHec, df_cornsoybean$CornPix, df_cornsoybean$CountyIndex, TRUE))
  fit <- eblup_bhf(CornHec ~ CornPix,
    unit_data = df_cornsoybean, Xpop = Xpop,
    domain_var = "CountyIndex", popsize_var = "PopnSegments", print_result = FALSE
  )
  expect_equal(fit$estcoef$std.error, unname(sqrt(diag(as.matrix(vcov(m))))), tolerance = 1e-5)
  expect_equal(residuals(fit), unname(residuals(m)), tolerance = 1e-5)
  expect_true(fit$convergence)
})

test_that("bootstrap MSE matches a reference refit-based bootstrap and is thread invariant", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  args <- list(
    formula = CornHec ~ CornPix + SoyBeansPix, unit_data = df_cornsoybean, Xpop = Xpop,
    domain_var = "CountyIndex", popsize_var = "PopnSegments", B = 30, seed = 3
  )
  r1 <- do.call(fastsae:::.pbmse_unit, args)
  r2 <- do.call(fastsae:::.pbmse_unit, c(args, n_threads = 2))
  expect_identical(r1$mse, r2$mse)

  # Reference: same draws (same order), lme4 refit, true population mean
  formuladata <- model.frame(args$formula, df_cornsoybean)
  Xs <- model.matrix(args$formula, formuladata)
  dom <- df_cornsoybean$CountyIndex
  sel <- unique(c(dom, Xpop$CountyIndex))
  xp <- Xpop[match(sel, Xpop$CountyIndex), ]
  Xbar <- cbind(1, xp$CornPix, xp$SoyBeansPix)
  di <- match(dom, sel)
  m0 <- lme4::lmer(CornHec ~ CornPix + SoyBeansPix + (1 | di), data.frame(formuladata, di = factor(di)))
  b <- lme4::fixef(m0)
  s2u <- lme4::VarCorr(m0)$di[1]
  s2e <- sigma(m0)^2
  nd <- tabulate(di, length(sel))
  N <- xp$PopnSegments
  rd <- N - nd
  set.seed(3)
  mse <- 0
  for (k in seq_len(30)) {
    u <- rnorm(length(sel), sd = sqrt(s2u))
    e <- rnorm(length(di), sd = sqrt(s2e))
    er <- rnorm(length(sel), 0, sqrt(s2e / rd))
    yb <- drop(Xs %*% b) + u[di] + e
    mb <- lme4::refit(m0, newresp = yb)
    ub <- lme4::ranef(mb)$di[, 1]
    bb <- lme4::fixef(mb)
    fd <- nd / N
    est <- fd * tapply(yb, di, mean) + drop((Xbar - fd * rowsum(Xs, di) / nd) %*% bb) + (1 - fd) * ub
    truth <- drop(Xbar %*% b) + u + (nd * tapply(e, di, mean) + rd * er) / N
    mse <- mse + (est - truth)^2 / 30
  }
  expect_equal(unname(r1$mse), as.numeric(mse), tolerance = 1e-3)
})

test_that("compiled BHF fit rejects unidentifiable designs", {
  x <- rnorm(20)
  expect_error(fastsae:::.bhf_fit(cbind(1, x, 2 * x), rnorm(20), rep(0:3, 5), 4, "REML"), "linearly dependent")
  expect_error(fastsae:::.bhf_fit(cbind(1, x), rnorm(20), 0:19, 20, "REML"), "must exceed")
})
