# Gamma hb_area() with and without vardir: y ~ Gamma(shape 20, mean mu), vardir the true
# sampling variance mu^2 / 20, default priors. 15/30/60 areas (1/8 without a direct
# estimate), area SD 0.3/0.8/1.5 on the log scale, 100 replicates; ~20 min on 6 cores.
#
# Run from the package root: Rscript benchmarks/validation/gamma_vardir.R
#
# Result (2026-10-19): the 95% intervals covered 83.4-88.1% of the true means without
# vardir and 91.2-94.6% with it (median relative error 0.138-0.188 and 0.137-0.162).
# Without vardir the gamma precision is estimated from one observation per area and
# trades off against the area effect; hb_area() warns.
devtools::load_all(quiet = TRUE)
q <- function(e) suppressWarnings(suppressMessages(e))
one <- function(m, su, r) {
  set.seed(7e6 + 1e3 * m + 100 * su + r)
  x <- rnorm(m)
  theta <- exp(-0.5 + 0.5 * x + rnorm(m, sd = su))
  y <- rgamma(m, shape = 20, rate = 20 / theta)
  miss <- seq_len(max(2, round(m / 8)))
  y[miss] <- NA
  d <- data.frame(y = y, x = x, v = theta^2 / 20)
  out <- NULL
  for (vd in c(FALSE, TRUE)) {
    f <- tryCatch(q(hb_area(y ~ x,
      data = d, family = "gamma", vardir = if (vd) "v" else NULL,
      print_result = FALSE, engine = inla_engine(num.threads = "1:1")
    )), error = function(e) NULL)
    if (is.null(f)) next
    h <- f$df_hb
    out <- rbind(out, data.frame(
      m = m, sigma_u = su, vardir = vd,
      cover = theta >= h$ci_lower & theta <= h$ci_upper, rel_err = abs(h$hb - theta) / theta
    ))
  }
  out
}
grid <- expand.grid(r = 1:100, su = c(0.3, 0.8, 1.5), m = c(15, 30, 60))
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  tryCatch(one(g$m, g$su, g$r), error = function(e) NULL)
}, mc.cores = 6))
a <- aggregate(cbind(cover, rel_err) ~ m + sigma_u + vardir, res, function(v) c(mean(v), median(v)))
a <- data.frame(a[1:3], cover = round(a$cover[, 1], 3), med_rel_err = round(a$rel_err[, 2], 3))
print(a[order(a$vardir, a$sigma_u, a$m), ], row.names = FALSE)
