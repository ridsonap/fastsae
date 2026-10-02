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
    compute_mse = TRUE,
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
      compute_mse = TRUE,
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

  # 4. .pbmse_unit directly with NA and fit_boot NULL
  # NOTE: with_mocked_bindings intercepts even qualified lme4::lmer calls,
  # so grab the real implementation first to avoid infinite recursion.
  real_lmer <- lme4::lmer
  lmer_env <- new.env(parent = emptyenv())
  lmer_env$call_count <- 0L
  testthat::with_mocked_bindings(
    lmer = function(formula, data, ...) {
      lmer_env$call_count <- lmer_env$call_count + 1L
      # First call is the initial fit; return real result
      if (lmer_env$call_count == 1L) return(real_lmer(formula, data = data, ...))
      # All subsequent (bootstrap) calls return NULL -> exercises `next` path
      NULL
    },
    .package = "lme4",
    {
      fit_lme4_na <- fastsae:::.pbmse_unit(
        formula = CornHec ~ CornPix + SoyBeansPix,
        unit_data = df_na,
        Xpop = Xpop,
        domain_var = "CountyIndex",
        popsize_var = "PopnSegments",
        B = 2,
        seed = 42
      )
      expect_type(fit_lme4_na, "list")
      expect_true(all(c("eblup", "mse", "B") %in% names(fit_lme4_na)))
    }
  )
})

