# Dense reference for the goodness-of-fit measures of the Gaussian area-level
# models: full (ML-form) and restricted log-likelihood of y ~ N(X beta, V) at
# the GLS estimate of beta.
.loglik_dense <- function(y, X, V) {
  n <- length(y)
  p <- ncol(X)
  Vi <- solve(V)
  XtViX <- t(X) %*% Vi %*% X
  beta <- solve(XtViX, t(X) %*% Vi %*% y)
  r <- y - X %*% beta
  ld <- function(A) as.numeric(determinant(A)$modulus)
  quad <- drop(t(r) %*% Vi %*% r)
  c(
    ml = -0.5 * (n * log(2 * pi) + ld(V) + quad),
    reml = -0.5 * ((n - p) * log(2 * pi) + ld(V) + ld(XtViX) - ld(t(X) %*% X) + quad)
  )
}
