mys_noise <- function() {
  data("mys", package = "fastsae", envir = environment())
  mys <- as.data.frame(mys)
  set.seed(1)
  mys$x5 <- rnorm(nrow(mys))
  mys
}

test_that("step_sae() follows emdi::step() for the Fay-Herriot model", {
  d <- mys_noise()
  fit <- eblup_fh(y ~ x1 + x2 + x3 + x5, vardir = "vardir", data = d, print_result = FALSE)
  s <- step_sae(fit, criterion = "BIC", trace = FALSE)
  expect_s3_class(s, "fastsae_step")
  # emdi 2.2.3, fh(method = "ml") and step(criteria = "BIC") on the same data:
  # 150.45 -> 147.18 (- x1) -> 143.91 (- x5) -> 141.04 (- x3)
  expect_equal(s$path$Step, c("", "- x1", "- x5", "- x3"))
  expect_equal(s$path$BIC, c(150.4542, 147.1784, 143.9140, 141.0402), tolerance = 1e-4)
  expect_equal(deparse(s$final), "y ~ x2")
  # emdi::step() selected y ~ x2 with AIC and KIC as well
  for (cr in c("AIC", "KIC")) {
    expect_equal(deparse(step_sae(fit, criterion = cr, trace = FALSE)$final), "y ~ x2")
  }
})

test_that("step_sae() searches with ML and refits the selected model with the original method", {
  d <- mys_noise()
  fit <- eblup_fh(y ~ x1 + x2 + x5, vardir = "vardir", data = d, method = "REML", print_result = FALSE)
  s <- step_sae(fit, trace = FALSE)
  ml <- eblup_fh(y ~ x1 + x2 + x5, vardir = "vardir", data = d, method = "ML", print_result = FALSE)
  expect_equal(s$path$BIC[1], ml$goodness[["BIC"]])
  expect_equal(s$model$method, "REML")
  expect_equal(deparse(s$model$formula), deparse(s$final))
})

test_that("step_sae() keeps the spatial structure and supports forward selection", {
  data("mys", package = "fastsae")
  fit <- eblup_sfh(y ~ x1, vardir = "vardir", data = mys, W = mys_proxmat, print_result = FALSE)
  s <- step_sae(fit, scope = ~ x1 + x2 + x3, direction = "forward", criterion = "AIC", trace = FALSE)
  expect_equal(s$model$model, fit$model)
  expect_true(all(attr(stats::terms(s$final), "term.labels") %in% c("x1", "x2", "x3")))
  expect_true("x1" %in% attr(stats::terms(s$final), "term.labels"))
  # the start is the AIC of the spatial fit
  ml <- eblup_sfh(y ~ x1, vardir = "vardir", data = mys, W = mys_proxmat, method = "ML", print_result = FALSE)
  expect_equal(s$path$AIC[1], ml$goodness[["AIC"]])
})

test_that("step_sae() works for eblup_bhf() with KIC", {
  data("cornsoybean", package = "fastsae")
  data("cornsoybeanmeans", package = "fastsae")
  pop <- cornsoybeanmeans
  names(pop)[names(pop) == "MeanCornPixPerSeg"] <- "CornPix"
  names(pop)[names(pop) == "MeanSoyBeansPixPerSeg"] <- "SoyBeansPix"
  names(pop)[names(pop) == "CountyIndex"] <- "County"
  fit <- eblup_bhf(CornHec ~ CornPix + SoyBeansPix,
    unit_data = cornsoybean, Xpop = pop, domain_var = "County",
    popsize_var = "PopnSegments", print_result = FALSE
  )
  s <- step_sae(fit, criterion = "KIC", trace = FALSE)
  ml <- eblup_bhf(CornHec ~ CornPix + SoyBeansPix,
    unit_data = cornsoybean, Xpop = pop, domain_var = "County",
    popsize_var = "PopnSegments", method = "ML", print_result = FALSE
  )
  g <- ml$goodness
  expect_equal(s$path$KIC[1], -2 * g[["loglikelihood"]] + 3 * 5)
  expect_true("CornPix" %in% attr(stats::terms(s$final), "term.labels"))
})

test_that("step_sae() selects HB models by WAIC, DIC or LCPO", {
  skip_if_not_installed("INLA")
  d <- mys_noise()
  fit <- hb_area(y ~ x1 + x2 + x5, data = d, vardir = "vardir", print_result = FALSE)
  s <- step_sae(fit, criterion = "WAIC", trace = FALSE)
  expect_equal(s$path$WAIC[1], fit$goodness[["WAIC"]])
  # leave-one-out (LCPO) is the default with INLA: WAIC is unstable with one
  # observation per area (benchmarks/validation/step_sae.R)
  expect_equal(step_sae(fit, direction = "backward", trace = FALSE)$criterion, "LCPO")
  expect_s3_class(s$model, "fastsae_hb")
  s2 <- step_sae(fit, criterion = "LCPO", direction = "backward", trace = FALSE)
  cpo <- fit$fit$cpo$cpo
  expect_equal(s2$path$LCPO[1], -mean(log(cpo[is.finite(cpo) & cpo > 0])))
  expect_error(step_sae(fit, criterion = "AIC"), "LCPO")
  expect_error(step_sae(eblup_fh(y ~ x1, vardir = "vardir", data = d, print_result = FALSE), criterion = "WAIC"), "AIC")
  expect_message(print(s), "Stepwise covariate selection")
})

test_that("step_sae() validates its input", {
  expect_error(step_sae(lm(mpg ~ wt, mtcars)), "fastsae")
})

test_that("step_sae() works for eblup_stfh() and eblup_twofold()", {
  data("mys_panel", package = "fastsae")
  data("mys_proxmat", package = "fastsae")
  keep <- !mys_panel$area %in% unique(mys_panel$area)[c(21, 25)]
  panel <- as.data.frame(mys_panel[keep, ])
  W <- mys_proxmat[-c(21, 25), -c(21, 25)]
  fit <- eblup_stfh(y ~ x1 + x2,
    data = panel, vardir = ~vardir, domain = ~area, time = ~year, W = W,
    print_result = FALSE
  )
  s <- step_sae(fit, direction = "backward", criterion = "AIC", trace = FALSE)
  # no ML fit: compared through the likelihood at the REML estimates
  expect_equal(s$path$AIC[1], fit$goodness[["AIC"]])
  expect_identical(s$model$call$W, fit$call$W)

  set.seed(2)
  dat <- do.call(rbind, lapply(1:15, function(d) {
    nd <- sample(2:4, 1)
    data.frame(area = d, subarea = paste0(d, "_", seq_len(nd)), x1 = rnorm(nd), x2 = rnorm(nd))
  }))
  dat$vardir <- runif(nrow(dat), 0.5, 2)
  dat$y <- 10 + 2 * dat$x1 + rnorm(15, sd = 1)[dat$area] + rnorm(nrow(dat), sd = sqrt(dat$vardir))
  tf <- eblup_twofold(y ~ x1 + x2,
    vardir = "vardir", domain = "area", subarea = "subarea", data = dat, print_result = FALSE
  )
  st <- step_sae(tf, trace = FALSE)
  expect_true("x1" %in% attr(stats::terms(st$final), "term.labels"))
  expect_equal(st$model$method, tf$method)
})

test_that("step_sae() selects hb_unit() and hb_twofold() models by WAIC, DIC or LCPO", {
  skip_if_not_installed("INLA")
  data("cornsoybean", package = "fastsae")
  data("cornsoybeanmeans", package = "fastsae")
  pop <- cornsoybeanmeans
  names(pop)[names(pop) == "MeanCornPixPerSeg"] <- "CornPix"
  names(pop)[names(pop) == "MeanSoyBeansPixPerSeg"] <- "SoyBeansPix"
  names(pop)[names(pop) == "CountyIndex"] <- "County"
  units <- cornsoybean
  set.seed(3)
  units$noise <- rnorm(nrow(units))
  pop$noise <- 0
  fu <- hb_unit(CornHec ~ CornPix + noise,
    unit_data = units, Xpop = pop, domain_var = "County",
    popsize_var = "PopnSegments", print_result = FALSE
  )
  # hb_unit() and hb_twofold() report goodness as list(dic, waic, ...), not DIC, WAIC
  for (cr in c("WAIC", "DIC")) {
    s <- step_sae(fu, criterion = cr, direction = "backward", trace = FALSE)
    expect_equal(s$path[[cr]][1], fu$goodness[[tolower(cr)]])
    expect_s3_class(s$model, "fastsae_hb_unit")
  }
  s <- step_sae(fu, direction = "backward", trace = FALSE)
  expect_equal(s$criterion, "LCPO")
  # the leave-one-out criterion of the sampled units; the population rows have no CPO
  cpo <- fu$fit$cpo$cpo
  expect_equal(s$path$LCPO[1], -mean(log(cpo[is.finite(cpo) & cpo > 0])))
  expect_true("CornPix" %in% attr(stats::terms(s$final), "term.labels"))

  set.seed(2)
  dat <- do.call(rbind, lapply(1:15, function(d) {
    nd <- sample(2:4, 1)
    data.frame(area = d, subarea = paste0(d, "_", seq_len(nd)), x1 = rnorm(nd), x2 = rnorm(nd))
  }))
  dat$vardir <- runif(nrow(dat), 0.5, 2)
  dat$y <- 10 + 2 * dat$x1 + rnorm(15, sd = 1)[dat$area] + rnorm(nrow(dat), sd = sqrt(dat$vardir))
  ft <- hb_twofold(y ~ x1 + x2,
    vardir = "vardir", domain = "area", subarea = "subarea", data = dat, print_result = FALSE
  )
  for (cr in c("WAIC", "DIC", "LCPO")) {
    st <- step_sae(ft, criterion = cr, trace = FALSE)
    expect_s3_class(st$model, "fastsae_hb_twofold")
    expect_true("x1" %in% attr(stats::terms(st$final), "term.labels"))
    if (cr != "LCPO") expect_equal(st$path[[cr]][1], ft$goodness[[tolower(cr)]])
  }
})
