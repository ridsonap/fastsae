# Default prior on the area SD of binomial and beta hb_area() (.default_prior_prec):
# frequentist coverage of the 95% intervals of p_d for INLA's PC(1, 0.01), half-normal(1)
# (the default) and half-normal(2.5) (as tipsae), on the logit scale.
# m = 15/30/60 areas (1/8 without a direct estimate), sigma_u = 0.3/0.8/1/1.5,
# logit p = -1 + 0.5 x + u, n_d ~ U{20, ..., 200}, 100 replicates; 6 cores, ~45 min.
#
# Run from the package root: Rscript benchmarks/validation/default_sd_prior_glm.R
#
# Result (2026-10-19; 24 cells = 2 families x 3 m x 4 sigma_u):
#   prior          coverage      RMSE of p_d   median area SD estimate
#   PC(1, 0.01)    0.918-0.944   equal         1-18% too small (1.23 for 1.5 with 15 areas)
#   HN(1)          0.927-0.950   equal         -3% to +15% (+15%: sigma_u = 0.3, 15 areas)
#   HN(2.5)        0.927-0.952   equal         -1% to +18% (1.69 for 1.5 with 15 areas)
# HN(2.5) covers at least as often as HN(1) in every cell (0.2 points more on average,
# at most 0.5); HN(1) estimates the area SD with less bias.
devtools::load_all(quiet = TRUE)
q <- function(e) suppressWarnings(suppressMessages(e))
priors <- list(
  "PC(1, 0.01)" = list(prior = "pc.prec", param = c(1, 0.01)),
  "HN(1)" = .half_normal_prec(1),
  "HN(2.5)" = .half_normal_prec(2.5)
)
R <- 100
one_rep <- function(fam, m, su, r) {
  set.seed(1e5 * (fam == "beta") + 1e3 * m + 100 * su + r)
  x <- rnorm(m)
  n <- sample(20:200, m, replace = TRUE)
  p <- plogis(-1 + 0.5 * x + rnorm(m, sd = su))
  miss <- seq_len(max(2, round(m / 8)))
  if (fam == "binomial") {
    y <- rbinom(m, n, p)
  } else {
    y <- rbeta(m, p * (n - 1), (1 - p) * (n - 1))
    y <- pmin(pmax(y, 1e-3), 1 - 1e-3)
  }
  y[miss] <- NA
  d <- data.frame(y = y, x = x, n = n)
  out <- NULL
  for (nm in names(priors)) {
    f <- tryCatch(q(hb_area(y ~ x,
      data = d, family = fam, trials = "n", print_result = FALSE,
      engine = inla_engine(prior_prec = priors[[nm]], num.threads = "1:1")
    )), error = function(e) NULL)
    if (is.null(f)) next
    h <- f$df_hb
    out <- rbind(out, data.frame(
      family = fam, m = m, sigma_u = su, prior = nm, unsamp = seq_len(m) %in% miss,
      cover = p >= h$ci_lower & p <= h$ci_upper, w = h$ci_upper - h$ci_lower, err2 = (h$hb - p)^2,
      s2u = f$random_effect_var
    ))
  }
  out
}
grid <- expand.grid(
  r = seq_len(R), su = c(0.3, 0.8, 1, 1.5), m = c(15, 30, 60), fam = c("binomial", "beta"),
  stringsAsFactors = FALSE
)
t0 <- Sys.time()
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  one_rep(g$fam, g$m, g$su, g$r)
}, mc.cores = 6))
a <- aggregate(cbind(cover, w, err2) ~ family + m + sigma_u + prior, res, mean)
a$rmse <- sqrt(a$err2)
a$err2 <- NULL
a$sigma_u_hat <- aggregate(sqrt(s2u) ~ family + m + sigma_u + prior, unique(res[, c("family", "m", "sigma_u", "prior", "s2u")]), median)[, 5]
print(a[order(a$family, a$sigma_u, a$m, a$prior), ], digits = 3, row.names = FALSE)
cat("elapsed:", format(Sys.time() - t0), "\n")
