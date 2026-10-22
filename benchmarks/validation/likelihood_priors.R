# Priors of the likelihood parameters that INLA and jax_engine() set differently without
# vardir: the residual SD of a Gaussian area-level model (INLA's default loggamma(1, 5e-5)
# on the precision vs half-normal(sd y), the prior of jax_engine()) and the beta precision
# (INLA's default loggamma(1, 0.1) vs half-normal(10)). With one observation per area both
# trade off against the area effect. Raw INLA with the model of hb_area(): iid area effect
# with the default half-normal SD prior. 30 areas (4 without a direct estimate),
# 60 replicates; results are saved per replicate in out_dir.
#
# Run from the package root:
#   Rscript benchmarks/validation/likelihood_priors.R [out_dir]
#
# Result (2026-10-21), coverage of 95% intervals for sampled areas:
#   Gaussian (area SD x residual SD 0.5/1 x 0.5/1, y in its unit / times 10):
#     half-normal(sd y) 0.938-0.968, INLA's default 0.001-0.038; INLA's default put the
#     residual precision at ~1.5e4, so the means collapsed onto the direct estimates.
#   Beta (area SD 0.3 / 0.8):
#     true phi   loggamma(1, 0.1)   half-normal(10)
#     10         0.960 / 0.937      0.965 / 0.938
#     50         0.931 / 0.903      0.913 / 0.869
#     200        0.875 / 0.925      0.841 / 0.871
#   Unsampled areas: Gaussian 0.92-1.00 with half-normal; beta 0.73-0.99 with either prior.
suppressMessages(library(INLA))
q <- function(e) suppressWarnings(suppressMessages(e))
out_dir <- if (length(commandArgs(TRUE)) > 0) commandArgs(TRUE)[1] else file.path(tempdir(), "likelihood_priors")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
hn_prec <- function(s) {
  list(prior = paste0("expression: s = ", s, "; sigma = exp(-theta/2); logdens = log(2) - log(s) - 0.5*log(2*pi) - 0.5*(sigma/s)^2 + log(sigma/2); return(logdens)"))
}
hn_size10 <- list(prior = "expression: size = exp(theta); logdens = log(2) - log(10) - 0.5*log(2*pi) - 0.5*(size/10)^2 + theta; return(logdens)")
fit_one <- function(d, family, s_u, lik_prior) {
  # the Gaussian precision is INLA's hyperparameter "prec", the beta precision "theta"
  cf <- if (is.null(lik_prior)) list() else list(hyper = setNames(list(lik_prior), if (family == "gaussian") "prec" else "theta"))
  r <- q(inla(y ~ x + f(id, model = "iid", hyper = list(prec = hn_prec(s_u))),
    data = d, family = family, control.family = cf,
    control.predictor = list(compute = TRUE, link = 1), num.threads = "1:1"
  ))
  list(fv = r$summary.fitted.values, disp = r$summary.hyperpar[1, "0.5quant"])
}
one <- function(g, r, m = 30) {
  set.seed(5e6 + 1000 * g$cell + r)
  x <- rnorm(m)
  miss <- 1:4
  if (g$family == "gaussian") {
    theta <- g$unit * (1 + 0.5 * x + rnorm(m, sd = g$su))
    y <- theta + rnorm(m, sd = g$unit * g$se)
    y[miss] <- NA
    s_y <- sd(y, na.rm = TRUE)
    priors <- list(inla_default = NULL, half_normal = hn_prec(s_y))
    s_u <- s_y
  } else {
    theta <- plogis(-1 + 0.5 * x + rnorm(m, sd = g$su))
    y <- rbeta(m, theta * g$phi, (1 - theta) * g$phi)
    y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
    y[miss] <- NA
    priors <- list(inla_default = NULL, half_normal = hn_size10)
    s_u <- 1
  }
  d <- data.frame(y = y, x = x, id = 1:m)
  out <- NULL
  for (p in names(priors)) {
    f <- fit_one(d, g$family, s_u, priors[[p]])
    h <- f$fv
    out <- rbind(out, data.frame(g, prior = p,
      unsamp = seq_len(m) %in% miss, cover = theta >= h[["0.025quant"]] & theta <= h[["0.975quant"]],
      rel_err = abs(h$mean - theta) / abs(theta), disp_hat = f$disp
    ))
  }
  out
}
cells <- rbind(
  expand.grid(family = "gaussian", su = c(0.5, 1), se = c(0.5, 1), unit = c(1, 10), phi = NA, stringsAsFactors = FALSE),
  expand.grid(family = "beta", su = c(0.3, 0.8), se = NA, unit = 1, phi = c(10, 50, 200), stringsAsFactors = FALSE)
)
cells$cell <- seq_len(nrow(cells))
grid <- expand.grid(r = 1:60, cell = cells$cell)
invisible(parallel::mclapply(seq_len(nrow(grid)), function(i) {
  f <- file.path(out_dir, sprintf("rep_%04d.rds", i))
  if (file.exists(f)) {
    return(NULL)
  }
  res <- tryCatch(one(cells[grid$cell[i], ], grid$r[i]), error = function(e) NULL)
  if (!is.null(res)) saveRDS(res, f)
  NULL
}, mc.cores = 6))
res <- do.call(rbind, lapply(list.files(out_dir, full.names = TRUE), readRDS))
cat("replicates done:", length(list.files(out_dir)), "of", nrow(grid), "\n")
res$se[is.na(res$se)] <- 0
res$phi[is.na(res$phi)] <- 0
a <- aggregate(cbind(cover, rel_err, disp_hat) ~ family + su + se + unit + phi + prior + unsamp, res, function(v) c(mean(v), median(v)))
a <- data.frame(a[1:7], cover = round(a$cover[, 1], 3), med_rel_err = round(a$rel_err[, 2], 3), disp_hat = signif(a$disp_hat[, 2], 3))
print(a[order(a$family, a$unsamp, a$su, a$se, a$unit, a$phi, a$prior), ], row.names = FALSE)
