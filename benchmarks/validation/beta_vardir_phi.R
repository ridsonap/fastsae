# Beta precision of each area from vardir. The beta variance is mu (1 - mu) / (phi + 1), so
# a known sampling variance V_d fixes phi_d = mu_d (1 - mu_d) / V_d - 1 at the true mean.
# Compared: the plug-in phi_d of the direct estimate y_d (hb_area() before), two stages
# with the posterior mean of a first fit (hb_area() now, both engines), phi_d at the model
# mean inside the model (tipsae::fit_sae(type_disp = "var"), with its own priors:
# half-normal(2.5) area SD, N(0, 2.5^2) coefficients) and the true phi_d (oracle).
# y_d ~ Beta with phi_d = n_d - 1, n_d 10-80, V_d = p_d (1 - p_d) / n_d;
# 30 areas (4 without a direct estimate), logit p = b0 + 0.5 x + u, 60 replicates.
#
# Run from the package root (INLA part on 6 cores, then tipsae one fit at a time):
#   Rscript benchmarks/validation/beta_vardir_phi.R [out_dir]
#
# Result (2026-10-22), coverage of 95% intervals for sampled areas (median relative error):
#   b0    area SD  plug-in y_d     two stages      tipsae          oracle
#   -2.5  0.3      0.844 (0.257)   0.938 (0.188)   0.938 (0.175)   0.953 (0.179)
#   -2.5  0.8      0.812 (0.332)   0.896 (0.269)   0.931 (0.246)   0.918 (0.263)
#   -1.2  0.3      0.922 (0.122)   0.948 (0.116)   0.951 (0.116)   0.950 (0.119)
#   -1.2  0.8      0.931 (0.163)   0.933 (0.161)   0.954 (0.161)   0.944 (0.163)
# Unsampled areas: plug-in 0.83-0.95, two stages 0.97-0.98, tipsae 0.93-0.97, oracle
# 0.95-0.97. tipsae: 11 of 240 fits with R-hat > 1.05 (max 8.7): the precision as a
# function of the model mean gives a second, low-precision mode of the likelihood.
suppressMessages(library(INLA))
q <- function(e) suppressWarnings(suppressMessages(e))
out_dir <- if (length(commandArgs(TRUE)) > 0) commandArgs(TRUE)[1] else file.path(tempdir(), "beta_vardir_phi")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
hn1 <- list(prior = "expression: s = 1; sigma = exp(-theta/2); logdens = log(2) - log(s) - 0.5*log(2*pi) - 0.5*(sigma/s)^2 + log(sigma/2); return(logdens)")
sim_data <- function(b0, su, r, m = 30) {
  set.seed(6e6 + 1e4 * abs(b0) + 100 * su + r)
  x <- rnorm(m)
  p <- plogis(b0 + 0.5 * x + rnorm(m, sd = su))
  n <- sample(10:80, m, replace = TRUE)
  y <- rbeta(m, p * (n - 1), (1 - p) * (n - 1))
  y <- pmin(pmax(y, 1e-4), 1 - 1e-4)
  V <- p * (1 - p) / n
  miss <- 1:4
  y[miss] <- NA
  V[miss] <- NA
  list(d = data.frame(y = y, x = x, id = 1:m, V = V), p = p, n = n, miss = miss)
}
fit_inla <- function(d, phi) {
  phi[is.na(phi)] <- 1
  r <- q(inla(y ~ x + f(id, model = "iid", hyper = list(prec = hn1)),
    data = d, family = "beta", scale = pmax(phi, 1),
    control.family = list(hyper = list(theta = list(initial = 0, fixed = TRUE))),
    control.predictor = list(compute = TRUE, link = 1), num.threads = "1:1"
  ))
  r$summary.fitted.values
}
summ <- function(h, s, method, cell) {
  data.frame(cell, method = method, unsamp = seq_along(s$p) %in% s$miss,
    cover = s$p >= h$lo & s$p <= h$hi, rel_err = abs(h$mean - s$p) / s$p, width = h$hi - h$lo
  )
}
cells <- expand.grid(b0 = c(-1.2, -2.5), su = c(0.3, 0.8))
grid <- expand.grid(r = 1:60, cell = seq_len(nrow(cells)))
invisible(parallel::mclapply(seq_len(nrow(grid)), function(i) {
  f <- file.path(out_dir, sprintf("inla_%04d.rds", i))
  if (file.exists(f)) {
    return(NULL)
  }
  cl <- cells[grid$cell[i], ]
  s <- sim_data(cl$b0, cl$su, grid$r[i])
  y <- s$d$y
  h1 <- fit_inla(s$d, y * (1 - y) / s$d$V - 1)
  mu <- h1$mean
  h2 <- fit_inla(s$d, mu * (1 - mu) / s$d$V - 1)
  h3 <- fit_inla(s$d, s$n - 1)
  tidy <- function(h) data.frame(mean = h$mean, lo = h[["0.025quant"]], hi = h[["0.975quant"]])
  res <- rbind(
    summ(tidy(h1), s, "plug_in_y", cl), summ(tidy(h2), s, "two_stage", cl),
    summ(tidy(h3), s, "oracle", cl)
  )
  saveRDS(res, f)
  NULL
}, mc.cores = 6))
for (i in seq_len(nrow(grid))) {
  f <- file.path(out_dir, sprintf("tipsae_%04d.rds", i))
  if (file.exists(f)) next
  cl <- cells[grid$cell[i], ]
  s <- sim_data(cl$b0, cl$su, grid$r[i])
  fit <- tryCatch(q(tipsae::fit_sae(y ~ x,
    data = s$d, domains = "id", disp_direct = "V", type_disp = "var",
    chains = 4, iter = 2000, seed = i, refresh = 0
  )), error = function(e) NULL)
  if (is.null(fit)) next
  e <- tipsae::extract(summary(fit))
  e <- rbind(e$in_sample, e$out_of_sample)
  e <- e[match(s$d$id, e$Domains), ]
  h <- data.frame(mean = e[["HB est."]], lo = e[["2.5%"]], hi = e[["97.5%"]])
  rhat <- max(rstan::summary(fit$stanfit)$summary[, "Rhat"], na.rm = TRUE)
  saveRDS(cbind(summ(h, s, "tipsae_model_mean", cl), rhat = rhat), f)
}
res <- do.call(rbind, lapply(list.files(out_dir, full.names = TRUE), function(f) readRDS(f)[, c("b0", "su", "method", "unsamp", "cover", "rel_err", "width")]))
cat("files:", length(list.files(out_dir)), "of", 2 * nrow(grid), "\n")
a <- aggregate(cbind(cover, rel_err, width) ~ b0 + su + method + unsamp, res, function(v) c(mean(v), median(v)))
a <- data.frame(a[1:4], cover = round(a$cover[, 1], 3), med_rel_err = round(a$rel_err[, 2], 3), mean_width = round(a$width[, 1], 4))
print(a[order(a$unsamp, a$b0, a$su, a$method), ], row.names = FALSE)
rh <- sapply(list.files(out_dir, "^tipsae", full.names = TRUE), function(f) readRDS(f)$rhat[1])
cat("\ntipsae: max R-hat", round(max(rh), 3), "; fits with R-hat > 1.05:", sum(rh > 1.05), "of", length(rh), "\n")
