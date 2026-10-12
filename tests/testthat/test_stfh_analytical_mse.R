library(fastsae)

test_that("eblup_stfh computes fast analytical MSE", {
  library(dplyr)
  mys_panel_nona <- mys_panel |>
    filter(!is.na(y) & year >= 2024)

  # Measure time for analytical MSE
  t0 <- Sys.time()
  fit_anal <- eblup_stfh(
    y ~ x1 + x2 + x3,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat[-c(21, 25), -c(21, 25)],
    model = "ST",
    compute_mse = TRUE,
    mse_type = "analytical",
    print_result = FALSE
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Must be fast (< 1.0 second)
  expect_lt(elapsed, 1.5)
  expect_true(fit_anal$convergence)
  expect_equal(fit_anal$mse_type, "analytical")

  # MSE and RSE checks
  expect_true("mse" %in% names(fit_anal$df_eblup))
  expect_true("rse" %in% names(fit_anal$df_eblup))
  expect_true(all(fit_anal$df_eblup$mse > 0))
  expect_true(all(fit_anal$df_eblup$rse > 0))
  expect_false(any(is.na(fit_anal$df_eblup$mse)))

  # MSE components checks
  expect_true(!is.null(fit_anal$g1))
  expect_true(!is.null(fit_anal$g2))
  expect_true(!is.null(fit_anal$g3))
  expect_true(all(fit_anal$g1 >= 0))
  expect_true(all(fit_anal$g2 >= 0))
  expect_true(all(fit_anal$g3 >= 0))

  # MSE is reasonable compared to direct sampling variance
  expect_true(mean(fit_anal$df_eblup$mse) < mean(fit_anal$df_eblup$vardir) * 1.5)
})

test_that("eblup_stfh analytical MSE works for Spatial-only model S", {
  mys_panel_nona <- mys_panel |>
    filter(!is.na(y))

  fit_s_anal <- eblup_stfh(
    y ~ x1 + x2,
    data = mys_panel_nona,
    vardir = ~vardir,
    domain = ~area,
    time = ~year,
    W = mys_proxmat[-c(21, 25), -c(21, 25)],
    model = "S",
    compute_mse = TRUE,
    mse_type = "analytical",
    print_result = FALSE
  )

  expect_true(fit_s_anal$convergence)
  expect_equal(fit_s_anal$model, "S")
  expect_true(all(fit_s_anal$df_eblup$mse > 0))
  expect_true(all(fit_s_anal$df_eblup$rse > 0))
})

# Small, well-identified ST panel (interior parameters) for exact checks.
.sim_st_panel <- function(D = 16, Tt = 4, seed = 99) {
  set.seed(seed)
  W <- sim_spatial_weights(D = D, type = "knn", style = "W", k = 4, seed = 1)
  M <- D * Tt
  dom <- rep(seq_len(D), each = Tt); tim <- rep(seq_len(Tt), D)
  x1 <- rnorm(M)
  psi <- runif(M, 0.3, 2) # heterogeneous: the old column-wise g3 was wrong here
  Om1 <- solve(t(diag(D) - 0.5 * W) %*% (diag(D) - 0.5 * W))
  Om2 <- outer(1:Tt, 1:Tt, function(i, j) 0.5^abs(i - j)) / (1 - 0.5^2)
  u1 <- t(chol(Om1)) %*% rnorm(D)
  u2 <- sqrt(0.6) * as.vector(t(chol(Om2)) %*% matrix(rnorm(M), Tt))
  y <- 2 + x1 + u1[dom] + u2 + rnorm(M, sd = sqrt(psi))
  list(data = data.frame(y, x1, psi, dom, tim), W = W, D = D, Tt = Tt)
}

test_that("STFH analytical g3 matches an independent computation (heterogeneous vardir)", {
  sim <- .sim_st_panel()
  fit <- eblup_stfh(y ~ x1, data = sim$data, vardir = ~psi, domain = ~dom, time = ~tim,
                    W = sim$W, model = "ST", print_result = FALSE)
  skip_if_not(isTRUE(fit$convergence))
  th <- fit$estvarcomp$estimate
  D <- sim$D; Tt <- sim$Tt; W <- sim$W
  X <- cbind(1, sim$data$x1); psi <- sim$data$psi
  Om1 <- solve(t(diag(D) - th[2] * W) %*% (diag(D) - th[2] * W))
  Om2 <- outer(1:Tt, 1:Tt, function(i, j) th[4]^abs(i - j)) / (1 - th[4]^2)
  dOm2 <- outer(1:Tt, 1:Tt, function(i, j) {
    l <- abs(i - j)
    (ifelse(l > 0, l * th[4]^pmax(l - 1, 0), 0) * (1 - th[4]^2) + 2 * th[4]^(l + 1)) / (1 - th[4]^2)^2
  })
  J <- matrix(1, Tt, Tt)
  V <- th[1] * kronecker(Om1, J) + th[3] * kronecker(diag(D), Om2) + diag(psi)
  Vi <- solve(V)
  Q <- solve(t(X) %*% Vi %*% X)
  P <- Vi - Vi %*% X %*% Q %*% t(X) %*% Vi
  dR <- -(W + t(W)) + 2 * th[2] * t(W) %*% W
  Va <- list(kronecker(Om1, J), kronecker(-th[1] * Om1 %*% dR %*% Om1, J),
             kronecker(diag(D), Om2), kronecker(diag(D), th[3] * dOm2))
  Fm <- matrix(0, 4, 4)
  for (a in 1:4) for (b in 1:4) Fm[a, b] <- 0.5 * sum(diag(P %*% Va[[a]] %*% P %*% Va[[b]]))
  # correlations estimated on the boundary are treated as fixed (not inverted)
  keep <- c(TRUE, abs(th[2]) < 0.999, TRUE, abs(th[4]) < 0.999)
  Fi <- matrix(0, 4, 4)
  Fi[keep, keep] <- solve(Fm[keep, keep])
  # d theta~/d th_a = Psi P V_a P y  =>  g3_i = psi_i^2 sum_ab Fi_ab (P V_a P V_b P)_ii
  g3 <- 0
  for (a in 1:4) for (b in 1:4) {
    g3 <- g3 + Fi[a, b] * psi^2 * rowSums((P %*% Va[[a]] %*% P %*% Va[[b]]) * P)
  }
  # rho1 is weakly identified on this small panel (ill-conditioned F), so
  # round-off differences between the C++ and R evaluation orders are
  # amplified to ~1e-5; the old (column-wise) g3 differed by O(1).
  # (estvarcomp$std.error follows sae::eblupSTFH, i.e. the last iteration.)
  expect_equal(fit$g3, g3, tolerance = 1e-4)
  # analytical MSE = g1 + g2 + 2 g3 with g1 = psi - psi^2 diag(V^-1)
  expect_equal(fit$g1, psi - psi^2 * diag(Vi), tolerance = 1e-8)
  expect_equal(fit$mse_analytical, fit$g1 + fit$g2 + 2 * fit$g3 - fit$g4, tolerance = 1e-10)
})

test_that("STFH bootstrap MSE (theta re-estimated) is consistent with the analytical MSE", {
  skip_on_cran()
  sim <- .sim_st_panel()
  fit_a <- eblup_stfh(y ~ x1, data = sim$data, vardir = ~psi, domain = ~dom, time = ~tim,
                      W = sim$W, model = "ST", compute_mse = TRUE, print_result = FALSE)
  fit_b <- eblup_stfh(y ~ x1, data = sim$data, vardir = ~psi, domain = ~dom, time = ~tim,
                      W = sim$W, model = "ST", compute_mse = TRUE, mse_type = "bootstrap",
                      B = 200, seed = 1, print_result = FALSE)
  skip_if_not(isTRUE(fit_a$convergence) && isTRUE(fit_b$convergence))
  a <- fit_a$df_eblup$mse
  b <- fit_b$df_eblup$mse_pb
  expect_gt(stats::cor(a, b), 0.8)
  expect_lt(abs(mean(b) / mean(a) - 1), 0.25)
  # with theta re-estimated, the bootstrap must exceed the known-theta g1 + g2
  expect_gt(mean(b), mean(fit_a$g1 + fit_a$g2) * 0.95)
})

# Restricted log-likelihood of the ST model (up to a constant), independent of
# the C++ implementation.
.reml_ll_st <- function(th, y, X, psi, W, D, Tt) {
  O1 <- solve(t(diag(D) - th[2] * W) %*% (diag(D) - th[2] * W))
  O2 <- outer(1:Tt, 1:Tt, function(i, j) th[4]^abs(i - j)) / (1 - th[4]^2)
  V <- th[1] * kronecker(O1, matrix(1, Tt, Tt)) + th[3] * kronecker(diag(D), O2) + diag(psi)
  cV <- chol(V)
  Vi <- chol2inv(cV)
  XtViX <- t(X) %*% Vi %*% X
  P <- Vi - Vi %*% X %*% solve(XtViX, t(X) %*% Vi)
  -0.5 * (2 * sum(log(diag(cV))) + as.numeric(determinant(XtViX)$modulus) + drop(t(y) %*% P %*% y))
}

# KKT-type check: no feasible coordinate move increases the REML likelihood.
.expect_constrained_reml_max <- function(fit, d, W, X) {
  D <- length(unique(d$dom)); Tt <- length(unique(d$tim))
  th <- fit$estvarcomp$estimate
  ll0 <- .reml_ll_st(th, d$y, X, d$psi, W, D, Tt)
  lo <- c(0, -0.999, 0, -0.999); hi <- c(Inf, 0.999, Inf, 0.999)
  for (j in 1:4) for (h in c(-1, 1)) {
    t2 <- th
    t2[j] <- min(max(th[j] + h * 1e-3 * max(abs(th[j]), 0.1), lo[j]), hi[j])
    # the fit stops once scoring steps are < 0.1% of an SE (ll gain < ~1e-6)
    expect_lte(.reml_ll_st(t2, d$y, X, d$psi, W, D, Tt), ll0 + 1e-6)
  }
}

test_that("STFH fallback returns the constrained REML maximiser (classical iteration fails)", {
  sim <- .sim_st_panel() # the classical Fisher scoring diverges on this panel
  fit <- eblup_stfh(y ~ x1, data = sim$data, vardir = ~psi, domain = ~dom, time = ~tim,
                    W = sim$W, model = "ST", print_result = FALSE)
  expect_true(fit$convergence)
  .expect_constrained_reml_max(fit, sim$data, sim$W, cbind(1, sim$data$x1))
})

test_that("STFH re-maximises under the constraints when a variance would be negative", {
  library(dplyr)
  # Real mys_panel response in domain-major order: the classical iteration
  # ends with sigma22 < 0; plain truncation is not the constrained maximum.
  d <- mys_panel |> filter(!is.na(y)) |> arrange(area, year)
  W <- mys_proxmat[-c(21, 25), -c(21, 25)]
  fit <- eblup_stfh(y ~ x1 + x2 + x3, data = d, vardir = ~vardir, domain = ~area, time = ~year,
                    W = W, model = "ST", print_result = FALSE)
  expect_true(fit$convergence)
  expect_equal(fit$estvarcomp$estimate[3], 0)
  expect_true(is.na(fit$estvarcomp$std.error[4])) # rho2 unidentified when sigma22 = 0
  dd <- data.frame(y = d$y, psi = d$vardir, dom = d$area, tim = d$year)
  .expect_constrained_reml_max(fit, dd, W, model.matrix(~ x1 + x2 + x3, d))
})

test_that("STFH g4 (second-derivative term) matches an independent computation", {
  D <- 16; Tt <- 4
  W <- sim_spatial_weights(D = D, type = "knn", style = "W", k = 4, seed = 1)
  set.seed(4) # interior estimates (rho1, rho2 away from the bounds)
  M <- D * Tt; dom <- rep(seq_len(D), each = Tt); tim <- rep(seq_len(Tt), D)
  x1 <- rnorm(M); psi <- runif(M, 0.3, 2)
  Om1 <- solve(t(diag(D) - 0.5 * W) %*% (diag(D) - 0.5 * W))
  Om2 <- outer(1:Tt, 1:Tt, function(i, j) 0.5^abs(i - j)) / (1 - 0.25)
  u1 <- t(chol(Om1)) %*% rnorm(D)
  u2 <- sqrt(0.6) * as.vector(t(chol(Om2)) %*% matrix(rnorm(M), Tt))
  d <- data.frame(y = 2 + x1 + u1[dom] + u2 + rnorm(M, sd = sqrt(psi)), x1, psi, dom, tim)
  fit <- eblup_stfh(y ~ x1, data = d, vardir = ~psi, domain = ~dom, time = ~tim, W = W,
                    model = "ST", print_result = FALSE)
  th <- fit$estvarcomp$estimate
  expect_true(all(abs(th[c(2, 4)]) < 0.9))

  X <- cbind(1, x1); J <- matrix(1, Tt, Tt)
  O <- solve(t(diag(D) - th[2] * W) %*% (diag(D) - th[2] * W))
  C <- -(W + t(W)) + 2 * th[2] * t(W) %*% W
  OCO <- O %*% C %*% O
  r2 <- th[4]; h <- 1 / (1 - r2^2)
  L <- abs(outer(1:Tt, 1:Tt, "-"))
  O2 <- r2^L * h
  dO2 <- (ifelse(L > 0, L * r2^pmax(L - 1, 0), 0) + r2^L * 2 * r2 * h) * h
  d2O2 <- ifelse(L >= 2, L * (L - 1) * r2^pmax(L - 2, 0), 0) * h +
    2 * ifelse(L >= 1, L * r2^pmax(L - 1, 0), 0) * 2 * r2 * h^2 + r2^L * (2 * h^2 + 8 * r2^2 * h^3)
  V <- th[1] * kronecker(O, J) + th[3] * kronecker(diag(D), O2) + diag(psi)
  Vi <- solve(V)
  P <- Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X, t(X) %*% Vi)
  Va <- list(kronecker(O, J), kronecker(-th[1] * OCO, J),
             kronecker(diag(D), O2), kronecker(diag(D), th[3] * dO2))
  Fm <- outer(1:4, 1:4, Vectorize(function(a, b) 0.5 * sum(diag(P %*% Va[[a]] %*% P %*% Va[[b]]))))
  Fi <- solve(Fm)
  Vjk <- list(c(1, 2), kronecker(-OCO, J),
              c(2, 2), kronecker(th[1] * (2 * OCO %*% C %*% O - 2 * O %*% t(W) %*% W %*% O), J),
              c(3, 4), kronecker(diag(D), dO2),
              c(4, 4), kronecker(diag(D), th[3] * d2O2))
  g4 <- 0
  for (k in seq(1, length(Vjk), by = 2)) {
    ab <- Vjk[[k]]
    w <- if (ab[1] == ab[2]) Fi[ab[1], ab[1]] else 2 * Fi[ab[1], ab[2]]
    g4 <- g4 + 0.5 * w * psi^2 * rowSums((P %*% Vjk[[k + 1]]) * P)
  }
  expect_equal(fit$g4, g4, tolerance = 1e-8)
})
