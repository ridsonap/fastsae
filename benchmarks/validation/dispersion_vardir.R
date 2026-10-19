# Gamma and negative binomial hb_area() with and without vardir. Without vardir the
# dispersion (gamma shape, NB size) is estimated from one observation per area and is
# hardly identified against the area effect; with vardir it is fixed per area (gamma:
# y^2 / vardir; NB: y^2 / (vardir - y)). Truth: gamma shape 20, NB size 10, vardir the
# true sampling variance. 15/30/60 areas (1/8 without a direct estimate), area SD
# 0.3/0.8/1.5 on the log scale, exposure 10-100 for NB, default priors, 100 replicates.
#
# Run from the package root: Rscript benchmarks/validation/dispersion_vardir.R (~1 hour)
#
# Result (2026-10-19), coverage of 95% intervals (median relative error):
#                  no vardir                       vardir
#   gamma          0.824-0.890 (0.139-0.184)       0.916-0.942 (0.135-0.163)
#   neg. binomial  0.940-0.961 (0.170-0.390)       0.908-0.953 (0.161-0.265)
# Gamma without vardir: intervals too short for every prior on the area SD (hb_area()
# warns). Negative binomial (true size 10) without vardir: covered thanks to the
# half-normal(10) size prior, whose median estimate stays near 8-12 whatever the true
# size (benchmarks/validation/nb_size_prior.R), with errors up to 50% larger than with
# vardir.
devtools::load_all(quiet = TRUE)
q <- function(e) suppressWarnings(suppressMessages(e))
one <- function(fam, m, su, r) {
  set.seed(7e6 + 1e5 * (fam == "nbinomial") + 1e3 * m + 100 * su + r)
  x <- rnorm(m)
  E <- runif(m, 10, 100)
  theta <- exp(-0.5 + 0.5 * x + rnorm(m, sd = su))
  if (fam == "gamma") {
    y <- rgamma(m, shape = 20, rate = 20 / theta)
    v <- theta^2 / 20
  } else {
    mu <- E * theta
    y <- rnbinom(m, mu = mu, size = 10)
    v <- mu + mu^2 / 10
  }
  miss <- seq_len(max(2, round(m / 8)))
  y[miss] <- NA
  d <- data.frame(y = y, x = x, E = E, v = v)
  out <- NULL
  for (vd in c(FALSE, TRUE)) {
    args <- list(y ~ x,
      data = d, family = fam, print_result = FALSE,
      engine = inla_engine(num.threads = "1:1")
    )
    if (vd) args$vardir <- "v"
    if (fam == "nbinomial") args$exposure <- "E"
    f <- tryCatch(q(do.call(hb_area, args)), error = function(e) NULL)
    if (is.null(f)) next
    h <- f$df_hb
    out <- rbind(out, data.frame(
      family = fam, m = m, sigma_u = su, vardir = vd, unsamp = seq_len(m) %in% miss,
      cover = theta >= h$ci_lower & theta <= h$ci_upper, rel_err = abs(h$hb - theta) / theta
    ))
  }
  out
}
grid <- expand.grid(
  r = 1:100, su = c(0.3, 0.8, 1.5), m = c(15, 30, 60), fam = c("gamma", "nbinomial"),
  stringsAsFactors = FALSE
)
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  tryCatch(one(g$fam, g$m, g$su, g$r), error = function(e) NULL)
}, mc.cores = 6))
a <- aggregate(cbind(cover, rel_err) ~ family + m + sigma_u + vardir, res, function(v) c(mean(v), median(v)))
a <- data.frame(a[1:4], cover = round(a$cover[, 1], 3), med_rel_err = round(a$rel_err[, 2], 3))
print(a[order(a$family, a$vardir, a$sigma_u, a$m), ], row.names = FALSE)
s <- aggregate(cover ~ family + vardir + unsamp, res, mean)
print(s, digits = 3, row.names = FALSE)
