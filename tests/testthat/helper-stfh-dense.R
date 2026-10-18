# Dense reference implementation of the spatio-temporal Fay-Herriot model on a
# D x T grid with unobserved cells (NA y). Used only by the tests: every
# quantity is formed from the full covariance matrices by direct algebra.

.stfh_dense_parts <- function(th, W, D, T, model = "ST") {
  I_D <- diag(D)
  Om1 <- solve(t(I_D - th[2] * W) %*% (I_D - th[2] * W))
  if (model == "ST") {
    lag <- abs(outer(seq_len(T), seq_len(T), "-"))
    Om2 <- th[4]^lag / (1 - th[4]^2)
  } else {
    Om2 <- diag(T)
  }
  Z <- kronecker(I_D, matrix(1, T, 1))
  # covariance of the signal u = Z u1 + u2 on the full grid
  G <- th[1] * Z %*% Om1 %*% t(Z) + th[3] * kronecker(I_D, Om2)
  list(G = G, Om1 = Om1, Om2 = Om2, Z = Z)
}

# Fit-free quantities at given th: BLUP of every cell and g1, g2
.stfh_dense_predict <- function(th, X, y, vardir, W, D, T, model = "ST") {
  s <- !is.na(y)
  pr <- .stfh_dense_parts(th, W, D, T, model)
  V <- pr$G[s, s] + diag(vardir[s], sum(s))
  Vi <- solve(V)
  Xs <- X[s, , drop = FALSE]
  Q <- solve(t(Xs) %*% Vi %*% Xs)
  beta <- Q %*% t(Xs) %*% Vi %*% y[s]
  A <- pr$G[, s] # Cov(theta_i, y_s), all cells
  pred <- drop(X %*% beta + A %*% Vi %*% (y[s] - Xs %*% beta))
  g1 <- diag(pr$G) - rowSums((A %*% Vi) * A)
  R <- X - A %*% Vi %*% Xs
  g2 <- rowSums((R %*% Q) * R)
  list(pred = pred, beta = drop(beta), g1 = g1, g2 = g2, Vi = Vi, Q = Q, V = V, s = s)
}

# Restricted log-likelihood of the observed cells
.stfh_dense_reml <- function(th, X, y, vardir, W, D, T, model = "ST") {
  s <- !is.na(y)
  pr <- .stfh_dense_parts(th, W, D, T, model)
  V <- pr$G[s, s] + diag(vardir[s], sum(s))
  Xs <- X[s, , drop = FALSE]
  Vi <- solve(V)
  XtViX <- t(Xs) %*% Vi %*% Xs
  P <- Vi - Vi %*% Xs %*% solve(XtViX) %*% t(Xs) %*% Vi
  -0.5 * (determinant(V)$modulus + determinant(XtViX)$modulus + drop(t(y[s]) %*% P %*% y[s]))
}

# Analytical MSE pieces at th by dense algebra (numerical derivatives of the
# covariance in th): g3 and g4 for every cell, with the REML information
# F_ab = 1/2 tr(P V_a P V_b) on the parameters in `free`.
.stfh_dense_mse <- function(th, X, y, vardir, W, D, T, model = "ST", free = seq_along(th), h = 1e-5) {
  s <- !is.na(y)
  np <- length(th)
  Gfun <- function(t) .stfh_dense_parts(t, W, D, T, model)$G
  G0 <- Gfun(th)
  dG <- lapply(seq_len(np), function(a) {
    e <- replace(numeric(np), a, h)
    (Gfun(th + e) - Gfun(th - e)) / (2 * h)
  })
  d2G <- function(a, b) {
    ea <- replace(numeric(np), a, h)
    eb <- replace(numeric(np), b, h)
    (Gfun(th + ea + eb) - Gfun(th + ea - eb) - Gfun(th - ea + eb) + Gfun(th - ea - eb)) / (4 * h^2)
  }
  V <- G0[s, s] + diag(vardir[s], sum(s))
  Vi <- solve(V)
  Xs <- X[s, , drop = FALSE]
  Q <- solve(t(Xs) %*% Vi %*% Xs)
  P <- Vi - Vi %*% Xs %*% Q %*% t(Xs) %*% Vi
  Va <- lapply(dG, function(M) M[s, s])
  Fm <- matrix(0, np, np)
  for (a in seq_len(np)) for (b in seq_len(np)) Fm[a, b] <- 0.5 * sum(diag(P %*% Va[[a]] %*% P %*% Va[[b]]))
  Finv <- matrix(0, np, np)
  Finv[free, free] <- solve(Fm[free, free, drop = FALSE])
  psi <- vardir
  n_all <- nrow(X)
  g3 <- g4 <- numeric(n_all)
  d2Gs <- list()
  for (a in seq_len(np)) for (b in seq_len(np)) if (Finv[a, b] != 0) d2Gs[[paste(a, b)]] <- d2G(a, b)
  # sampled cells: Prasad-Rao form in P (as the C++ code)
  g3s <- g4s <- numeric(sum(s))
  for (a in seq_len(np)) {
    for (b in seq_len(np)) {
      if (Finv[a, b] == 0) next
      g3s <- g3s + Finv[a, b] * diag(P %*% Va[[a]] %*% P %*% Va[[b]] %*% P)
      g4s <- g4s + 0.5 * Finv[a, b] * diag(P %*% d2Gs[[paste(a, b)]][s, s] %*% P)
    }
  }
  g3[s] <- psi[s]^2 * g3s
  g4[s] <- psi[s]^2 * g4s
  # unsampled cells: g3 from the exact derivative of the BLUP, g4 from d2 g1
  g1fun <- function(t) {
    G <- Gfun(t)
    Vt <- G[s, s] + diag(vardir[s], sum(s))
    diag(G)[!s] - rowSums((G[!s, s, drop = FALSE] %*% solve(Vt)) * G[!s, s, drop = FALSE])
  }
  idx_u <- which(!s)
  if (length(idx_u)) {
    A <- G0[!s, s, drop = FALSE]
    R <- X[!s, , drop = FALSE] - A %*% Vi %*% Xs
    Bm <- A %*% Vi + R %*% Q %*% t(Xs) %*% Vi
    L <- lapply(seq_len(np), function(a) dG[[a]][!s, s, drop = FALSE] - Bm %*% Va[[a]])
    g3u <- numeric(length(idx_u))
    for (a in seq_len(np)) for (b in seq_len(np)) g3u <- g3u + Finv[a, b] * rowSums((L[[a]] %*% P) * L[[b]])
    g4u <- numeric(length(idx_u))
    for (a in seq_len(np)) {
      for (b in seq_len(np)) {
        if (Finv[a, b] == 0) next
        ea <- replace(numeric(np), a, h)
        eb <- replace(numeric(np), b, h)
        d2 <- (g1fun(th + ea + eb) - g1fun(th + ea - eb) - g1fun(th - ea + eb) + g1fun(th - ea - eb)) / (4 * h^2)
        g4u <- g4u + 0.5 * Finv[a, b] * d2
      }
    }
    g3[idx_u] <- g3u
    g4[idx_u] <- g4u
  }
  list(g3 = g3, g4 = g4, Finv = Finv)
}
