library(testthat)
library(fastsae)

# Two-fold data with non-sampled subareas in sampled areas and one area without
# any sampled subarea (area 12).
sim_twofold_area <- function(seed = 11, m = 12) {
  set.seed(seed)
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(2:5, 1)
    data.frame(area = d, subarea = paste0(d, "-", seq_len(nd)),
               x1 = rnorm(nd), vardir = runif(nd, 0.2, 0.9), pop = rpois(nd, 50) + 1)
  }))
  dat$y <- 1 + dat$x1 + rnorm(m, sd = 0.8)[dat$area] + rnorm(nrow(dat), sd = 0.5) +
    rnorm(nrow(dat), sd = sqrt(dat$vardir))
  dat$y[c(2, 9)] <- NA
  dat$y[dat$area == m] <- NA
  dat
}

test_that("eblup_twofold area estimates are the weighted means of the subarea EBLUPs", {
  dat <- sim_twofold_area()
  fit <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                       data = dat, weight = "pop", print_result = FALSE)

  expect_named(fit$df_area, c("domain", "eblup_area", "mse_area", "rse_area",
                              "n_subareas", "n_sampled"))
  expect_equal(fit$df_area$domain, sort(unique(dat$area)))
  w <- dat$pop / ave(dat$pop, dat$area, FUN = sum)
  ref <- as.numeric(rowsum(w * fit$df_eblup$eblup, dat$area))
  expect_equal(fit$df_area$eblup_area, ref, tolerance = 1e-10)
  expect_equal(fit$df_area$n_subareas, as.integer(table(dat$area)))
  expect_equal(fit$df_area$n_sampled, as.integer(tapply(!is.na(dat$y), dat$area, sum)))
  expect_true(all(fit$df_area$mse_area > 0))

  # equal weights by default
  fit_eq <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                          data = dat, print_result = FALSE)
  ref_eq <- as.numeric(tapply(fit_eq$df_eblup$eblup, dat$area, mean))
  expect_equal(fit_eq$df_area$eblup_area, ref_eq, tolerance = 1e-10)
})

test_that("eblup_twofold area MSE with a one-subarea weight equals that subarea's MSE", {
  dat <- sim_twofold_area()
  # weight 1 on one subarea per area: sampled, non-sampled (rows 2, 9) and
  # subareas of the area without a sample
  pick <- c(2, 9, which(dat$area == 12)[1], which(dat$area == 5)[2])
  first <- !duplicated(dat$area)
  sel <- first
  sel[dat$area %in% dat$area[pick]] <- FALSE
  sel[pick] <- TRUE
  dat$w1 <- as.numeric(sel)

  for (method in c("REML", "ML")) {
    fit <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                         data = dat, weight = "w1", method = method, print_result = FALSE)
    sub <- fit$df_eblup[sel, ]
    expect_equal(fit$df_area$eblup_area, sub$eblup, tolerance = 1e-10)
    expect_equal(fit$df_area$mse_area, sub$mse, tolerance = 1e-8, info = method)
  }
})

test_that("eblup_twofold analytical area MSE matches a dense pure-R g1 + g2 + 2 g3", {
  dat <- sim_twofold_area(seed = 3, m = 10)
  fit <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                       data = dat, weight = "pop", print_result = FALSE)
  s2v <- fit$random_effect_var[["sigma2_v"]]
  s2u <- fit$random_effect_var[["sigma2_u"]]

  X <- cbind(1, dat$x1)
  s <- !is.na(dat$y)
  Xs <- X[s, ]
  ar_s <- dat$area[s]
  w <- dat$pop / ave(dat$pop, dat$area, FUN = sum)
  areas <- sort(unique(dat$area))
  Zv <- outer(ar_s, areas, "==") * 1          # sampled rows x areas
  Vfun <- function(dl) dl[1] * Zv %*% t(Zv) + diag(dl[2] + dat$vardir[s])

  # REML information I_rt = tr(P V_r P V_t) / 2 at the estimates
  V <- Vfun(c(s2v, s2u))
  Vi <- solve(V)
  Q <- solve(t(Xs) %*% Vi %*% Xs)
  P <- Vi - Vi %*% Xs %*% Q %*% t(Xs) %*% Vi
  Vr <- list(Zv %*% t(Zv), diag(sum(s)))
  info <- matrix(0, 2, 2)
  for (r in 1:2) for (t in 1:2) info[r, t] <- 0.5 * sum(diag(P %*% Vr[[r]] %*% P %*% Vr[[t]]))
  info_inv <- solve(info)

  mse_ref <- numeric(length(areas))
  for (k in seq_along(areas)) {
    a <- dat$area == areas[k]
    wa <- ifelse(a, w, 0)
    kvec <- colSums(wa * X)
    # Cov(y_s, mu) and b(delta) = V^{-1} Cov(y_s, mu) for the area target
    s_fun <- function(dl) dl[1] * sum(wa) * (ar_s == areas[k]) + dl[2] * wa[s]
    b_fun <- function(dl) solve(Vfun(dl), s_fun(dl))
    sv <- s_fun(c(s2v, s2u))
    b <- b_fun(c(s2v, s2u))
    g1 <- s2v * sum(wa)^2 + s2u * sum(wa^2) - sum(sv * b)
    dvec <- kvec - colSums(b * Xs)
    g2 <- drop(t(dvec) %*% Q %*% dvec)
    h <- 1e-6
    db <- cbind((b_fun(c(s2v + h, s2u)) - b_fun(c(s2v - h, s2u))) / (2 * h),
                (b_fun(c(s2v, s2u + h)) - b_fun(c(s2v, s2u - h))) / (2 * h))
    g3 <- sum(diag(t(db) %*% V %*% db %*% info_inv))
    mse_ref[k] <- g1 + g2 + 2 * g3
  }
  expect_equal(fit$df_area$mse_area, mse_ref, tolerance = 1e-5)
})

test_that("eblup_twofold bootstrap area MSE agrees with the analytical area MSE", {
  set.seed(5)
  m <- 30
  dat <- do.call(rbind, lapply(seq_len(m), function(d) {
    nd <- sample(3:6, 1)
    data.frame(area = d, subarea = seq_len(nd), x1 = rnorm(nd), vardir = runif(nd, 0.3, 0.8))
  }))
  dat$y <- 2 + dat$x1 + rnorm(m)[dat$area] + rnorm(nrow(dat), sd = 0.5) +
    rnorm(nrow(dat), sd = sqrt(dat$vardir))

  fa <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                      data = dat, print_result = FALSE)
  fb <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                      data = dat, mse_type = "bootstrap", B = 400, seed = 1,
                      print_result = FALSE)
  expect_equal(fb$df_area$eblup_area, fa$df_area$eblup_area)
  ratio <- fb$df_area$mse_area / fa$df_area$mse_area
  expect_gt(median(ratio), 0.8)
  expect_lt(median(ratio), 1.25)
})

test_that("eblup_twofold compute_area, mse_type = 'none' and invalid weights", {
  dat <- sim_twofold_area()
  fit <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                       data = dat, compute_area = FALSE, print_result = FALSE)
  expect_null(fit$df_area)

  fit_none <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                            data = dat, mse_type = "none", print_result = FALSE)
  expect_true(all(is.na(fit_none$df_area$mse_area)))
  expect_false(anyNA(fit_none$df_area$eblup_area))

  dat$pop[1] <- -1
  expect_warning(
    eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                  data = dat, weight = "pop", print_result = FALSE),
    "missing or negative"
  )
  expect_error(
    eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                  data = dat, weight = dat$pop[-1], print_result = FALSE),
    "length does not match"
  )
})

test_that("eblup_twofold area estimates follow the self-benchmarked sub-area estimates", {
  dat <- sim_twofold_area()
  fit <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area", subarea = "subarea",
                       data = dat, weight = "pop", self_benchmark = TRUE,
                       benchmark_weight = "pop", print_result = FALSE)
  w <- dat$pop / ave(dat$pop, dat$area, FUN = sum)
  expect_equal(fit$df_area$eblup_area,
               as.numeric(rowsum(w * fit$df_eblup$eblup, dat$area)), tolerance = 1e-10)
  expect_true(all(fit$df_area$mse_area > 0))
})
