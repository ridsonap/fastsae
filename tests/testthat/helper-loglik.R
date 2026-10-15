# Dense reference for the goodness-of-fit measures of the Gaussian area-level
# models: full (ML-form) log-likelihood of y ~ N(X beta, V) at the GLS
# estimate of beta.
.loglik_dense <- function(y, X, V) {
  Vi <- solve(V)
  beta <- solve(t(X) %*% Vi %*% X, t(X) %*% Vi %*% y)
  r <- y - X %*% beta
  -0.5 * (length(y) * log(2 * pi) + as.numeric(determinant(V)$modulus) +
            drop(t(r) %*% Vi %*% r))
}
