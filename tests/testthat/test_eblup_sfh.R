library(testthat)
library(fastsae)

skip_if_not_installed("sae")

quiet <- function(expr) {
  result <- NULL
  invisible(capture.output(result <- suppressWarnings(suppressMessages(expr))))
  result
}

# ------------------------------------------------------------------
# Shared objects
# ------------------------------------------------------------------

idx <- !is.na(mys$y)

mysnona <- mys[idx, ]
mys_proxmat_nona <- mys_proxmat[idx, idx]

fit_fast <- eblup_sfh(
  y ~ x1 + x2 + x3,
  vardir = "vardir",
  domain = "area",
  method = "REML",
  data = mysnona,
  W = mys_proxmat_nona,
  print_result = FALSE
)

fit_sae <- sae::mseSFH(
  mysnona$y ~ mysnona$x1 + mysnona$x2 + mysnona$x3,
  vardir = mysnona$vardir,
  proxmat = mys_proxmat_nona,
  method = "REML"
)


fit_fast_ml <- eblup_sfh(
  y ~ x1 + x2 + x3,
  vardir = "vardir",
  domain = "area",
  method = "ML",
  data = mysnona,
  W = mys_proxmat_nona,
  print_result = FALSE
)

fit_sae_ml <- sae::mseSFH(
  mysnona$y ~ mysnona$x1 + mysnona$x2 + mysnona$x3,
  vardir = mysnona$vardir,
  proxmat = mys_proxmat_nona,
  method = "ML"
)

tol <- 1e-6


# -----------------------------------------------------------------
# Structure
# ------------------------------------------------------------------

test_that("eblup_sfh returns valid structure", {
  expect_true(is.list(fit_fast))

  # Standardized fields
  expect_true("estcoef" %in% names(fit_fast))
  expect_true("random_effect_var" %in% names(fit_fast))
  expect_true("rho" %in% names(fit_fast))
  expect_true("estvarcomp" %in% names(fit_fast))
  expect_true("goodness" %in% names(fit_fast))
  expect_true("df_eblup" %in% names(fit_fast))
  expect_true("model" %in% names(fit_fast))
  expect_true("level" %in% names(fit_fast))
  expect_true("n_iter" %in% names(fit_fast))
  expect_true("convergence" %in% names(fit_fast))
  expect_true("method" %in% names(fit_fast))

  # Check estvarcomp structure
  expect_true("parameter" %in% names(fit_fast$estvarcomp))
  expect_true("estimate" %in% names(fit_fast$estvarcomp))
  expect_true("std.error" %in% names(fit_fast$estvarcomp))

  # Check values
  expect_equal(fit_fast$model, "SFH")
  expect_equal(fit_fast$level, "area")
  expect_true(fit_fast$convergence)
  expect_true(fit_fast$n_iter > 0)

  expect_length(
    fit_fast$df_eblup$eblup,
    nrow(mysnona)
  )

  expect_length(
    fit_fast$df_eblup$mse,
    nrow(mysnona)
  )
})

# ------------------------------------------------------------------
# EBLUP
# ------------------------------------------------------------------

test_that("EBLUP agrees with sae::mseSFH", {
  expect_equal(
    fit_fast$df_eblup$eblup,
    as.numeric(fit_sae$est$eblup),
    tolerance = tol
  )

  expect_equal(
    fit_fast_ml$df_eblup$eblup,
    as.numeric(fit_sae_ml$est$eblup),
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# MSE
# ------------------------------------------------------------------

# Independent R implementation of the SFH analytical MSE (sampled domains):
# g1 + g2 + 2 g3 - g4 (Singh et al. 2005 / sae::mseSFH), minus b_ML' grad g1
# for ML. Information: REML (P) for REML, ML (V^-1) for ML. A rho on the
# boundary is treated as fixed (only the sigma2 block is inverted).
.sfh_mse_R <- function(fit, X, y, psi, W) {
  m <- nrow(W); s2 <- fit$random_effect_var; rho <- fit$rho
  Sig <- solve(crossprod(diag(m) - rho * W))
  C <- 2 * rho * t(W) %*% W - W - t(W)
  der3 <- Sig %*% C %*% Sig
  G <- s2 * Sig; V <- G + diag(psi); Vi <- solve(V)
  Q <- solve(t(X) %*% Vi %*% X); P <- Vi - Vi %*% X %*% Q %*% t(X) %*% Vi
  Va <- list(Sig, -s2 * der3)
  M <- if (fit$method == "ML") Vi else P
  I <- outer(1:2, 1:2, Vectorize(function(a, b) 0.5 * sum(diag(M %*% Va[[a]] %*% M %*% Va[[b]]))))
  Ii <- if (abs(rho) >= 0.999) matrix(c(1 / I[1, 1], 0, 0, 0), 2) else solve(I)
  g1 <- diag(G - G %*% Vi %*% G)
  R <- X - G %*% Vi %*% X; g2 <- rowSums((R %*% Q) * R)
  L <- list(Vi %*% Va[[1]] - Vi %*% Va[[1]] %*% Vi %*% G, Vi %*% Va[[2]] - Vi %*% Va[[2]] %*% Vi %*% G)
  g3 <- 0; for (a in 1:2) for (b in 1:2) g3 <- g3 + Ii[a, b] * colSums(L[[a]] * (V %*% L[[b]]))
  D12 <- -der3; D22 <- 2 * s2 * der3 %*% C %*% Sig - 2 * s2 * Sig %*% t(W) %*% W %*% Sig
  Ps <- diag(psi)
  g4 <- 0.5 * diag(Ps %*% Vi %*% D12 %*% Vi %*% Ps) * (Ii[1, 2] + Ii[2, 1]) +
    0.5 * diag(Ps %*% Vi %*% D22 %*% Vi %*% Ps) * Ii[2, 2]
  mse <- g1 + g2 + 2 * g3 - g4
  if (fit$method == "ML") {
    h <- sapply(1:2, function(a) -sum(diag(Q %*% t(X) %*% Vi %*% Va[[a]] %*% Vi %*% X)))
    b <- Ii %*% h / 2
    dg <- sapply(1:2, function(a) diag(Va[[a]] - 2 * G %*% Vi %*% Va[[a]] + G %*% Vi %*% Va[[a]] %*% Vi %*% G))
    mse <- mse - as.vector(dg %*% b)
  }
  mse
}

test_that("MSE agrees with sae::mseSFH (REML) and with an independent computation (ML)", {
  # Interior estimates; on mysnona rho-hat is on the boundary (next test).
  skip_if_not_installed("MASS")
  pop <- simulate_safe_sfh(rho_true = 0.7)
  set.seed(1)
  d <- data.frame(y = pop$theta_true + rnorm(pop$m, sd = sqrt(pop$vardir)),
                  x1 = pop$X[, 2], x2 = pop$X[, 3], x3 = pop$X[, 4], vardir = pop$vardir)
  s <- sae::mseSFH(d$y ~ d$x1 + d$x2 + d$x3, vardir = d$vardir, proxmat = pop$W, method = "REML")
  f <- eblup_sfh(y ~ x1 + x2 + x3, vardir = "vardir", data = d, W = pop$W, method = "REML",
                 print_result = FALSE)
  expect_true(abs(f$rho) < 0.999)
  expect_equal(f$df_eblup$mse, s$mse, tolerance = tol)
  expect_equal(f$df_eblup$mse, .sfh_mse_R(f, pop$X, d$y, d$vardir, pop$W),
               tolerance = 1e-8, ignore_attr = TRUE)

  # ML uses the ML information (sae::mseSFH uses the REML one), so it is checked
  # against the independent computation instead
  f_ml <- eblup_sfh(y ~ x1 + x2 + x3, vardir = "vardir", data = d, W = pop$W, method = "ML",
                    print_result = FALSE)
  expect_true(abs(f_ml$rho) < 0.999)
  expect_equal(f_ml$df_eblup$mse, .sfh_mse_R(f_ml, pop$X, d$y, d$vardir, pop$W),
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("MSE with rho-hat on the boundary treats rho as fixed", {
  X <- model.matrix(~ x1 + x2 + x3, mysnona)
  for (fit in list(fit_fast, fit_fast_ml)) {
    expect_equal(fit$rho, -1)
    expect_true(is.na(fit$estvarcomp$std.error[2]))
    expect_equal(fit$df_eblup$mse, .sfh_mse_R(fit, X, mysnona$y, mysnona$vardir, mys_proxmat_nona),
                 tolerance = 1e-8, ignore_attr = TRUE)
  }
})

# ------------------------------------------------------------------
# Variance component
# ------------------------------------------------------------------

test_that("random effect variance agrees with sae::mseSFH", {
  expect_equal(
    fit_fast$random_effect_var,
    fit_sae$est$fit$refvar,
    tolerance = tol
  )

  expect_equal(
    fit_fast_ml$random_effect_var,
    fit_sae_ml$est$fit$refvar,
    tolerance = tol
  )
})

test_that("rho parameter is in reasonable range", {
  # Spatial autocorrelation should be in (-1, 1)
  expect_true(fit_fast$rho >= -1 && fit_fast$rho <= 1)
})

# ------------------------------------------------------------------
# Goodness of fit
# ------------------------------------------------------------------

test_that("goodness statistics agree with sae::mseSFH", {
  expect_equal(
    as.numeric(fit_fast$goodness),
    as.numeric(fit_sae$est$fit$goodness[-4]),
    tolerance = tol
  )

  expect_equal(
    as.numeric(fit_fast_ml$goodness),
    as.numeric(fit_sae_ml$est$fit$goodness[-4]),
    tolerance = tol
  )
})


# ------------------------------------------------------------------
# Regression coefficients
# ------------------------------------------------------------------

test_that("beta estimates agree with sae::mseSFH", {
  expect_identical(
    names(fit_fast$estcoef$beta),
    names(fit_sae$est$fit$estcoef$beta)
  )

  expect_equal(
    fit_fast$estcoef$beta,
    fit_sae$est$fit$estcoef$beta,
    tolerance = tol
  )

  expect_identical(
    names(fit_fast_ml$estcoef$beta),
    names(fit_sae_ml$est$fit$estcoef$beta)
  )

  expect_equal(
    fit_fast_ml$estcoef$beta,
    fit_sae_ml$est$fit$estcoef$beta,
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# Standard errors
# ------------------------------------------------------------------

test_that("beta standard errors agree with sae::mseSFH", {
  expect_equal(
    fit_fast$estcoef$stderr_beta,
    fit_sae$est$fit$estcoef$std.error,
    tolerance = tol
  )

  expect_equal(
    fit_fast_ml$estcoef$stderr_beta,
    fit_sae_ml$est$fit$estcoef$std.error,
    tolerance = tol
  )
})

# ------------------------------------------------------------------
# Error handling
# ------------------------------------------------------------------

test_that("non-square proximity matrix throws error", {
  W_bad <- mys_proxmat_nona[-1, ]

  expect_error(
    eblup_sfh(
      y ~ x1 + x2 + x3,
      vardir = "vardir",
      data = mysnona,
      W = W_bad,
      print_result = FALSE
    )
  )
})

test_that("negative sampling variance throws error", {
  dat_bad <- mysnona
  dat_bad$vardir[1] <- -1

  expect_error(
    eblup_sfh(
      y ~ x1 + x2 + x3,
      vardir = "vardir",
      data = dat_bad,
      W = mys_proxmat_nona,
      print_result = FALSE
    )
  )
})

test_that("eblup_sfh validation and edge cases", {
  # W is NULL
  expect_error(eblup_sfh(y ~ x1, vardir = "vardir", data = mysnona, print_result = FALSE), "must be provided")

  # domain = NULL defaults to 1:nrow(data)
  fit_nodom <- eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mysnona, W = mys_proxmat_nona, domain = NULL, print_result = FALSE)
  expect_equal(fit_nodom$df_eblup$domain, 1:nrow(mysnona))

  # NA in auxiliary variables
  bad_mys <- mysnona
  bad_mys$x1[1] <- NA
  expect_error(eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = bad_mys, W = mys_proxmat_nona, print_result = FALSE), "Auxiliary variables contain NA")

  # Length of vardir mismatch
  expect_error(eblup_sfh(y ~ x1 + x2, vardir = 1:5, data = mysnona, W = mys_proxmat_nona, print_result = FALSE), "length does not match")

  # Non-convergence with maxiter = 1
  quiet(
    fit_nc <- eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mysnona, W = mys_proxmat_nona, maxiter = 1, print_result = FALSE)
  )

  expect_false(fit_nc$convergence)

  # print_result = TRUE
  expect_output(eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mysnona, W = mys_proxmat_nona, print_result = TRUE))

  # Bootstrap MSE with unsampled domains
  fit_pb_unsampled <- eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mys, W = mys_proxmat, mse_method = "pbmse", B = 5, print_result = FALSE)
  expect_s3_class(fit_pb_unsampled, "fastsae")

  fit_npb_unsampled <- eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mys, W = mys_proxmat, mse_method = "npbmse", B = 5, print_result = FALSE)
  expect_s3_class(fit_npb_unsampled, "fastsae")

  # Direct nrow(mf) != length(vardir) check in eblup_sfh
  testthat::with_mocked_bindings(
    .get_variable = function(data, variable) c(1, 2),
    .package = "fastsae",
    {
      expect_error(
        eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mysnona, W = mys_proxmat_nona, print_result = FALSE),
        "Length of 'vardir' must equal number of observations"
      )
    }
  )
})

# ------------------------------------------------------------------
# Unsampled domains: SAR is defined on ALL domains of W, so the sampled
# domains have covariance s2 * Sigma_full[s,s] (not the SAR built from W[s,s])
# and kriging uses the same model.
# ------------------------------------------------------------------
test_that("SFH with unsampled domains uses the full-domain SAR marginal covariance", {
  fit <- eblup_sfh(y ~ x1 + x2 + x3, vardir = "vardir", data = mys, W = mys_proxmat,
                   print_result = FALSE)
  s <- !is.na(mys$y)
  D <- nrow(mys_proxmat)
  rho <- fit$rho
  s2 <- fit$random_effect_var
  Sig <- solve(t(diag(D) - rho * mys_proxmat) %*% (diag(D) - rho * mys_proxmat))
  X <- model.matrix(~ x1 + x2 + x3, mys)
  V <- s2 * Sig[s, s] + diag(mys$vardir[s])
  Vi <- solve(V)
  beta <- solve(t(X[s, ]) %*% Vi %*% X[s, ], t(X[s, ]) %*% Vi %*% mys$y[s])
  r <- mys$y[s] - X[s, ] %*% beta
  expect_equal(fit$estcoef$beta, as.vector(beta), tolerance = 1e-6)
  expect_equal(fit$df_eblup$eblup[s], as.vector(X[s, ] %*% beta + s2 * Sig[s, s] %*% Vi %*% r),
               tolerance = 1e-6)
  expect_equal(fit$df_eblup$eblup[!s], as.vector(X[!s, ] %*% beta + s2 * Sig[!s, s] %*% Vi %*% r),
               tolerance = 1e-6)
})

test_that("SFH bootstrap validity does not depend on the scale of y", {
  skip_if_not_installed("MASS")
  sim <- simulate_safe_sfh(rho_true = 0.7)
  set.seed(11)
  k <- 1e4 # y ~ 1e4-1e5 => eblup well above the old absolute 1e6 cut-off
  d <- data.frame(
    y = (sim$theta_true + rnorm(sim$m, sd = sqrt(sim$vardir))) * k,
    x1 = sim$X[, "x1"], x2 = sim$X[, "x2"], x3 = sim$X[, "x3"],
    vardir = sim$vardir * k^2
  )
  fit_s <- suppressWarnings(eblup_sfh(y ~ x1 + x2 + x3, vardir = "vardir", data = d, W = sim$W,
                                      mse_method = "pbmse", B = 10, seed = 1,
                                      print_result = FALSE))
  expect_true(all(is.finite(fit_s$df_eblup$mse_pb)))
  expect_gt(max(fit_s$df_eblup$mse_pb), 1e6)
})

# ------------------------------------------------------------------
# Constrained (RE)ML at the boundary
# ------------------------------------------------------------------
.sfh_grid_data <- function(seed) {
  W <- make_grid_W(6, 6)
  m <- 36
  set.seed(seed)
  x1 <- rnorm(m)
  vd <- runif(m, 0.5, 1.5)
  u <- as.numeric(MASS::mvrnorm(1, rep(0, m), 0.05 * solve(crossprod(diag(m) + 0.5 * W))))
  list(data = data.frame(y = 1 + x1 + u + rnorm(m, sd = sqrt(vd)), x1, vardir = vd), W = W)
}

.sfh_reml_ll <- function(s2, rho, d, W) {
  m <- nrow(W)
  X <- cbind(1, d$x1)
  Sig <- solve(crossprod(diag(m) - rho * W))
  V <- s2 * Sig + diag(d$vardir)
  Vi <- solve(V)
  XtViX <- t(X) %*% Vi %*% X
  P <- Vi - Vi %*% X %*% solve(XtViX, t(X) %*% Vi)
  -0.5 * (as.numeric(determinant(V)$modulus) + as.numeric(determinant(XtViX)$modulus) +
            drop(t(d$y) %*% P %*% d$y))
}

test_that("SFH returns the constrained REML maximiser when sigma2 is at 0", {
  skip_if_not_installed("MASS")
  g <- .sfh_grid_data(27)
  fit <- eblup_sfh(y ~ x1, vardir = "vardir", data = g$data, W = g$W, print_result = FALSE)
  expect_true(fit$convergence)
  expect_equal(fit$random_effect_var, 0)
  # KKT: sigma2 = 0 with no feasible improvement in sigma2; rho is then unidentified
  ll0 <- .sfh_reml_ll(0, fit$rho, g$data, g$W)
  for (h in c(1e-3, 1e-2)) expect_lt(.sfh_reml_ll(h, fit$rho, g$data, g$W), ll0 + 1e-8)
})

test_that("rho = -0.999 is not mapped to -1 when I + W is singular (bipartite grid)", {
  skip_if_not_installed("MASS")
  g <- .sfh_grid_data(30)
  fit <- eblup_sfh(y ~ x1, vardir = "vardir", data = g$data, W = g$W, print_result = FALSE)
  expect_equal(fit$rho, -0.999)
  expect_true(all(is.finite(fit$df_eblup$mse)))
  expect_lt(max(fit$df_eblup$mse), max(g$data$vardir) * 2)
})

test_that("SFH MSE of kriged (unsampled) domains matches an independent second-order computation", {
  skip_if_not_installed("MASS")
  pop <- simulate_safe_sfh(rho_true = 0.6)
  set.seed(3)
  y <- pop$theta_true + rnorm(pop$m, sd = sqrt(pop$vardir))
  ns <- c(5, 18, 33, 47, 60)
  y[ns] <- NA
  d <- data.frame(y, x1 = pop$X[, 2], x2 = pop$X[, 3], x3 = pop$X[, 4], vardir = pop$vardir)
  W <- pop$W; Dn <- nrow(W); s <- setdiff(seq_len(Dn), ns)
  X <- pop$X; Xs <- X[s, ]; Xr <- X[ns, ]; psi <- pop$vardir[s]
  for (meth in c("REML", "ML")) {
    f <- eblup_sfh(y ~ x1 + x2 + x3, vardir = "vardir", data = d, W = W, method = meth,
                   print_result = FALSE)
    th <- c(f$random_effect_var, f$rho)
    Gf <- function(t) t[1] * solve(crossprod(diag(Dn) - t[2] * W))
    G <- Gf(th); Vi <- solve(G[s, s] + diag(psi)); Q <- solve(t(Xs) %*% Vi %*% Xs)
    P <- Vi - Vi %*% Xs %*% Q %*% t(Xs) %*% Vi
    A <- G[ns, s]; Rr <- Xr - A %*% Vi %*% Xs
    g1 <- function(t) {
      Gt <- Gf(t)
      diag(Gt[ns, ns] - Gt[ns, s] %*% solve(Gt[s, s] + diag(psi), Gt[s, ns]))
    }
    # first/second derivatives by finite differences (independent of the C++ algebra)
    e <- function(a, h) replace(c(0, 0), a, h)
    dG <- lapply(1:2, function(a) (Gf(th + e(a, 1e-5)) - Gf(th - e(a, 1e-5))) / 2e-5)
    Va <- lapply(dG, function(M) M[s, s]); Aa <- lapply(dG, function(M) M[ns, s])
    Mi <- if (meth == "ML") Vi else P # ML information for ML, REML information for REML
    Fm <- outer(1:2, 1:2, Vectorize(function(a, b) 0.5 * sum(diag(Mi %*% Va[[a]] %*% Mi %*% Va[[b]]))))
    Fi <- solve(Fm)
    Bm <- A %*% Vi + Rr %*% Q %*% t(Xs) %*% Vi
    L <- lapply(1:2, function(a) Aa[[a]] - Bm %*% Va[[a]])
    g3 <- 0
    for (a in 1:2) for (b in 1:2) g3 <- g3 + Fi[a, b] * rowSums((L[[a]] %*% P) * L[[b]])
    hh <- 1e-4
    H <- function(j, k) (g1(th + e(j, hh) + e(k, hh)) - g1(th + e(j, hh) - e(k, hh)) -
                           g1(th - e(j, hh) + e(k, hh)) + g1(th - e(j, hh) - e(k, hh))) / (4 * hh^2)
    half_tr <- 0
    for (j in 1:2) for (k in 1:2) half_tr <- half_tr + 0.5 * Fi[j, k] * H(j, k)
    mse <- g1(th) + rowSums((Rr %*% Q) * Rr) + g3 - half_tr
    if (meth == "ML") {
      hvec <- sapply(1:2, function(a) -sum(diag(Q %*% t(Xs) %*% Vi %*% Va[[a]] %*% Vi %*% Xs)))
      grad <- sapply(1:2, function(a) (g1(th + e(a, hh)) - g1(th - e(a, hh))) / (2 * hh))
      mse <- mse - as.vector(grad %*% (Fi %*% hvec / 2))
    }
    expect_equal(f$df_eblup$mse[ns], mse, tolerance = 1e-6)
  }
})

test_that("SFH with a dense W (dense factorisation path) matches the independent computation", {
  # inverse-distance weights: every pair is a neighbour, so the C++ engine
  # uses dense instead of sparse Cholesky factorisations
  set.seed(29)
  D <- 50
  xy <- cbind(runif(D), runif(D))
  W <- 1 / as.matrix(dist(xy)); diag(W) <- 0; W <- W / rowSums(W)
  X <- cbind(1, rnorm(D), rnorm(D))
  psi <- runif(D, 0.3, 1)
  v <- solve(diag(D) - 0.5 * W, rnorm(D, sd = 0.8))
  d <- data.frame(y = drop(X %*% c(1, 1, -1)) + v + rnorm(D, sd = sqrt(psi)),
                  x1 = X[, 2], x2 = X[, 3], vardir = psi)
  for (meth in c("REML", "ML")) {
    f <- eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = d, W = W, method = meth, print_result = FALSE)
    expect_true(f$convergence)
    expect_true(abs(f$rho) < 0.999)
    expect_equal(f$df_eblup$mse, .sfh_mse_R(f, X, d$y, psi, W), tolerance = 1e-8, ignore_attr = TRUE)
    # GLS beta, EBLUP and log-likelihood from dense algebra
    G <- f$random_effect_var * solve(crossprod(diag(D) - f$rho * W))
    V <- G + diag(psi); Vi <- solve(V)
    b <- solve(t(X) %*% Vi %*% X, t(X) %*% Vi %*% d$y)
    expect_equal(f$estcoef$beta, as.numeric(b), tolerance = 1e-10)
    expect_equal(f$df_eblup$eblup, as.numeric(X %*% b + G %*% Vi %*% (d$y - X %*% b)), tolerance = 1e-10)
    r <- d$y - X %*% b
    ll <- -0.5 * (D * log(2 * pi) + as.numeric(determinant(V)$modulus) + sum(r * (Vi %*% r)))
    expect_equal(unname(f$goodness["loglikelihood"]), ll, tolerance = 1e-10)
  }
})
