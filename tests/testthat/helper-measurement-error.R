# Area data for hb_area(var_x = ): the covariate w = x + e is observed with error
# variance vx, z is observed without error, and the first n_unsampled areas have no
# direct estimate. Used by test_hb_area_me.R and test_hb_area_engine.R.
sim_me_df <- function(D = 50, beta = 0.8, seed = 1, n_unsampled = 0) {
  set.seed(seed)
  x <- rnorm(D, 5, 2)
  vx <- runif(D, 0.2, 1.5)
  vardir <- runif(D, 0.5, 1.5)
  theta <- 1 + beta * x + rnorm(D, 0, 0.7)
  d <- data.frame(
    area = paste0("a", seq_len(D)),
    w = x + rnorm(D, 0, sqrt(vx)),
    vx = vx,
    z = rnorm(D),
    vardir = vardir,
    theta = theta
  )
  d$y <- theta + rnorm(D, 0, sqrt(vardir))
  if (n_unsampled > 0) d$y[seq_len(n_unsampled)] <- NA
  d
}
