# Prior on the gamma precision (shape) without vardir: INLA's default loggamma(1, 0.01)
# on log precision vs half-normal(10) on the precision, the prior of jax_engine(). Raw INLA
# with the model of hb_area(): iid area effect with a half-normal(1) SD, gamma likelihood
# with a common precision. 30 areas (4 without a direct estimate), true shape 5/20/100,
# area SD 0.3/0.8 on the log scale, 100 replicates; results are saved per replicate in
# out_dir.
#
# Run from the package root:
#   Rscript benchmarks/validation/gamma_shape_prior.R [out_dir]
#
# Result (2026-10-22), coverage of 95% intervals (area SD 0.3 / 0.8):
#   true shape  sampled: HN(10)  sampled: loggamma  unsampled: HN(10)  unsampled: loggamma
#   5           0.930 / 0.925    0.593 / 0.522      0.972 / 0.938      0.995 / 0.955
#   20          0.959 / 0.986    0.862 / 0.826      0.907 / 0.932      0.973 / 0.955
#   100         0.950 / 0.999    0.979 / 0.993      0.847 / 0.925      0.935 / 0.948
# The estimated shape follows its prior: median 8-13 with HN(10), 76-84 with loggamma, for
# every true shape. HN(10) costs wider intervals and larger errors for nearly noiseless
# responses (median relative error 0.14-0.17 against 0.07 for shape 100).
suppressMessages(library(INLA))
q <- function(e) suppressWarnings(suppressMessages(e))
out_dir <- if (length(commandArgs(TRUE)) > 0) commandArgs(TRUE)[1] else file.path(tempdir(), "gamma_shape_prior")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
hn_sd1 <- list(prior = "expression: s = 1; sigma = exp(-theta/2); logdens = log(2) - log(s) - 0.5*log(2*pi) - 0.5*(sigma/s)^2 + log(sigma/2); return(logdens)")
shape_priors <- list(
  loggamma = list(prior = "loggamma", param = c(1, 0.01)),
  hn10 = list(prior = "expression: size = exp(theta); logdens = log(2) - log(10) - 0.5*log(2*pi) - 0.5*(size/10)^2 + theta; return(logdens)")
)
one <- function(shape, su, r, m = 30) {
  set.seed(4e6 + 1e4 * shape + 100 * su + r)
  x <- rnorm(m)
  theta <- exp(1 + 0.5 * x + rnorm(m, sd = su))
  y <- rgamma(m, shape = shape, rate = shape / theta)
  miss <- 1:4
  y[miss] <- NA
  d <- data.frame(y = y, x = x, id = 1:m)
  out <- NULL
  for (p in names(shape_priors)) {
    fit <- q(inla(y ~ x + f(id, model = "iid", hyper = list(prec = hn_sd1)),
      data = d, family = "gamma",
      control.family = list(hyper = list(theta = shape_priors[[p]])),
      control.predictor = list(compute = TRUE, link = 1), num.threads = "1:1"
    ))
    h <- fit$summary.fitted.values
    out <- rbind(out, data.frame(
      true_shape = shape, sigma_u = su, prior = p,
      unsamp = seq_len(m) %in% miss, cover = theta >= h[["0.025quant"]] & theta <= h[["0.975quant"]],
      rel_err = abs(h$mean - theta) / theta, shape_hat = fit$summary.hyperpar[1, "0.5quant"]
    ))
  }
  out
}
grid <- expand.grid(r = 1:100, su = c(0.3, 0.8), shape = c(5, 20, 100))
invisible(parallel::mclapply(seq_len(nrow(grid)), function(i) {
  f <- file.path(out_dir, sprintf("rep_%04d.rds", i))
  if (file.exists(f)) {
    return(NULL)
  }
  res <- tryCatch(one(grid$shape[i], grid$su[i], grid$r[i]), error = function(e) NULL)
  if (!is.null(res)) saveRDS(res, f)
  NULL
}, mc.cores = 6))
res <- do.call(rbind, lapply(list.files(out_dir, full.names = TRUE), readRDS))
cat("replicates done:", length(list.files(out_dir)), "of", nrow(grid), "\n")
a <- aggregate(cbind(cover, rel_err) ~ true_shape + sigma_u + prior + unsamp, res, function(v) c(mean(v), median(v)))
a <- data.frame(a[1:4], cover = round(a$cover[, 1], 3), med_rel_err = round(a$rel_err[, 2], 3))
print(a[order(a$unsamp, a$true_shape, a$sigma_u, a$prior), ], row.names = FALSE)
sh <- aggregate(shape_hat ~ true_shape + sigma_u + prior, unique(res[, c("true_shape", "sigma_u", "prior", "shape_hat")]), median)
cat("\nmedian estimated shape:\n")
print(sh, digits = 3, row.names = FALSE)
