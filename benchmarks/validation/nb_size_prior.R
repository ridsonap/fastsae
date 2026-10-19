# Prior on the negative binomial size: INLA's pc.mgamma(7) vs half-normal(10), with one
# estimated size (no vardir) and with vardir (two stages: the size of each area from the
# posterior mean of a first fit, .nb_size_from_vardir). Raw INLA with the model of
# hb_area(): iid area effect with a half-normal(1) SD, NB with exposure 10-100.
# 30 areas (4 without a direct estimate), true size 3/10/100, area SD 0.3/0.8 on the log
# scale, 100 replicates; ~10 min on 6 cores. Results are saved per replicate in out_dir.
#
# Run from the package root:
#   Rscript benchmarks/validation/nb_size_prior.R [out_dir]
#
# Result (2026-10-19), coverage of 95% intervals for sampled areas (area SD 0.3 / 0.8):
#   true size  no vardir: HN(10)  no vardir: PC(7)  two stages: HN(10)  two stages: PC(7)
#   3          0.928 / 0.896      0.816 / 0.756     0.906 / 0.903       0.871 / 0.893
#   10         0.949 / 0.965      0.913 / 0.915     0.938 / 0.931       0.927 / 0.930
#   100        0.939 / 0.990      0.965 / 0.980     0.954 / 0.950       0.958 / 0.951
# The estimated size follows its prior: median 8-12 with HN(10), 110-180 with PC(7), for
# every true size. With two stages HN(10) was as good or better (median relative error
# 0.31-0.33 against 0.42-0.49 for size 3); without vardir it was safer when the
# overdispersion was real, at the cost of wider intervals and larger errors (0.15-0.20
# against 0.12-0.15) for nearly Poisson counts. Unsampled areas: 0.89-0.99.
suppressMessages(library(INLA))
q <- function(e) suppressWarnings(suppressMessages(e))
out_dir <- if (length(commandArgs(TRUE)) > 0) commandArgs(TRUE)[1] else file.path(tempdir(), "nb_size_prior")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
hn_sd1 <- list(prior = "expression: s = 1; sigma = exp(-theta/2); logdens = log(2) - log(s) - 0.5*log(2*pi) - 0.5*(sigma/s)^2 + log(sigma/2); return(logdens)")
size_priors <- list(
  pc_mgamma7 = list(prior = "pc.mgamma", param = 7),
  hn10 = list(prior = "expression: size = exp(theta); logdens = log(2) - log(10) - 0.5*log(2*pi) - 0.5*(size/10)^2 + theta; return(logdens)")
)
fit_nb <- function(d, size_prior, scale = NULL) {
  cf <- if (is.null(scale)) {
    list(hyper = list(theta = size_prior))
  } else {
    list(variant = 2, hyper = list(theta = list(initial = 0, fixed = TRUE)))
  }
  args <- list(y ~ x + f(id, model = "iid", hyper = list(prec = hn_sd1)),
    data = d, family = "nbinomial", E = d$E, control.family = cf,
    control.predictor = list(compute = TRUE, link = 1), num.threads = "1:1"
  )
  if (!is.null(scale)) args$scale <- scale
  r <- q(do.call(inla, args))
  list(fv = r$summary.fitted.values, size = if (is.null(scale)) r$summary.hyperpar[1, "0.5quant"] else NA)
}
nb_size <- function(mu, V) {
  s <- rep(1e8, length(mu))
  ok <- V > mu & mu > 0
  s[ok] <- mu[ok]^2 / (V[ok] - mu[ok])
  pmin(s, 1e8)
}
one <- function(size, su, r, m = 30) {
  set.seed(3e6 + 1e4 * size + 100 * su + r)
  x <- rnorm(m)
  E <- runif(m, 10, 100)
  theta <- exp(-0.5 + 0.5 * x + rnorm(m, sd = su))
  mu <- E * theta
  y <- rnbinom(m, mu = mu, size = size)
  V <- mu + mu^2 / size
  miss <- 1:4
  y[miss] <- NA
  d <- data.frame(y = y, x = x, E = E, id = 1:m)
  out <- NULL
  for (p in names(size_priors)) {
    s1 <- fit_nb(d, size_priors[[p]])
    s2 <- fit_nb(d, size_priors[[p]], scale = nb_size(s1$fv$mean * E, V))
    for (meth in c("no_vardir", "two_stage")) {
      h <- if (meth == "no_vardir") s1$fv else s2$fv
      out <- rbind(out, data.frame(
        true_size = size, sigma_u = su, prior = p, method = meth,
        unsamp = seq_len(m) %in% miss, cover = theta >= h[["0.025quant"]] & theta <= h[["0.975quant"]],
        rel_err = abs(h$mean - theta) / theta, size_hat = s1$size
      ))
    }
  }
  out
}
grid <- expand.grid(r = 1:100, su = c(0.3, 0.8), size = c(3, 10, 100))
invisible(parallel::mclapply(seq_len(nrow(grid)), function(i) {
  f <- file.path(out_dir, sprintf("rep_%04d.rds", i))
  if (file.exists(f)) {
    return(NULL)
  }
  res <- tryCatch(one(grid$size[i], grid$su[i], grid$r[i]), error = function(e) NULL)
  if (!is.null(res)) saveRDS(res, f)
  NULL
}, mc.cores = 6))
res <- do.call(rbind, lapply(list.files(out_dir, full.names = TRUE), readRDS))
cat("replicates done:", length(list.files(out_dir)), "of", nrow(grid), "\n")
a <- aggregate(cbind(cover, rel_err) ~ true_size + sigma_u + prior + method + unsamp, res, function(v) c(mean(v), median(v)))
a <- data.frame(a[1:5], cover = round(a$cover[, 1], 3), med_rel_err = round(a$rel_err[, 2], 3))
cat("\nsampled areas:\n")
print(a[!a$unsamp, -5][order(a$true_size[!a$unsamp], a$sigma_u[!a$unsamp], a$method[!a$unsamp]), ], row.names = FALSE)
cat("\nunsampled areas:\n")
print(a[a$unsamp, -5][order(a$true_size[a$unsamp], a$sigma_u[a$unsamp], a$method[a$unsamp]), c("true_size", "sigma_u", "prior", "method", "cover")], row.names = FALSE)
sz <- aggregate(size_hat ~ true_size + sigma_u + prior, unique(res[res$method == "no_vardir", c("true_size", "sigma_u", "prior", "size_hat")]), median)
cat("\nmedian estimated size (one size):\n")
print(sz, digits = 3, row.names = FALSE)
