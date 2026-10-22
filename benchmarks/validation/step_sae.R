# Does step_sae() find the true covariates? x1, x2 enter the model; x3, x4 are noise.
#   fh:        Fay-Herriot (eblup_fh), 40 areas, BIC and AIC.
#   binomial:  hb_area(family = "binomial"), 40 areas, n_d 20-100, WAIC and LCPO.
#   spatial:   area effect with a spatial (SAR) component, noise covariate x4 spatially
#              smooth; selection with eblup_sfh() (the spatial model) against eblup_fh().
# 50 replicates (30 for binomial); share of replicates selecting exactly {x1, x2}, and the
# share keeping each noise covariate.
#
# Run from the package root: Rscript benchmarks/validation/step_sae.R
#
# Result (2026-10-22), share of replicates (exact = exactly {x1, x2}):
#   setting                    exact   x1     x2     x3     x4
#   FH, BIC                    0.72    1.00   0.88   0.10   0.12
#   FH, AIC                    0.60    1.00   1.00   0.16   0.30
#   binomial HB, WAIC          0.27    0.87   0.80   0.47   0.33
#   binomial HB, LCPO          0.70    1.00   1.00   0.10   0.23
#   spatial truth, SFH, BIC    0.74    0.98   0.90   0.10   0.08
#   spatial truth, FH, BIC     0.68    0.96   0.90   0.06   0.20
# WAIC is unstable with one observation per area (hence LCPO as the HB default with INLA);
# the non-spatial model keeps the spatially smooth noise covariate x4 more often.
suppressMessages(devtools::load_all(quiet = TRUE))
q <- function(e) suppressWarnings(suppressMessages(e))
W <- sim_spatial_weights(D = 40, type = "knn", style = "W", k = 4, seed = 1)
sel_terms <- function(s) attr(stats::terms(s$final), "term.labels")
summarise <- function(sets, label) {
  data.frame(
    setting = label, exact = mean(vapply(sets, function(t) setequal(t, c("x1", "x2")), logical(1))),
    x1 = mean(vapply(sets, function(t) "x1" %in% t, logical(1))),
    x2 = mean(vapply(sets, function(t) "x2" %in% t, logical(1))),
    x3 = mean(vapply(sets, function(t) "x3" %in% t, logical(1))),
    x4 = mean(vapply(sets, function(t) "x4" %in% t, logical(1)))
  )
}
sim_area <- function(r, spatial = FALSE) {
  set.seed(100 + r)
  m <- 40
  d <- data.frame(area = 1:m, x1 = rnorm(m), x2 = rnorm(m), x3 = rnorm(m))
  if (spatial) {
    s <- as.numeric(solve(diag(m) - 0.8 * W, rnorm(m)))
    d$x4 <- as.numeric(solve(diag(m) - 0.8 * W, rnorm(m)))
    u <- 0.6 * s / sd(s) + rnorm(m, sd = 0.3)
  } else {
    d$x4 <- rnorm(m)
    u <- rnorm(m, sd = 0.5)
  }
  d$theta <- 5 + 0.5 * d$x1 - 0.4 * d$x2 + u
  d$vardir <- runif(m, 0.1, 0.4)
  d$y <- d$theta + rnorm(m, sd = sqrt(d$vardir))
  d
}
res <- list()
fh_sets <- list(BIC = list(), AIC = list())
for (r in 1:50) {
  d <- sim_area(r)
  fit <- q(eblup_fh(y ~ x1 + x2 + x3 + x4, vardir = "vardir", data = d, print_result = FALSE))
  for (cr in names(fh_sets)) fh_sets[[cr]][[r]] <- sel_terms(q(step_sae(fit, criterion = cr, trace = FALSE)))
}
res$fh_bic <- summarise(fh_sets$BIC, "FH, BIC")
res$fh_aic <- summarise(fh_sets$AIC, "FH, AIC")

bin_sets <- list(WAIC = list(), LCPO = list())
for (r in 1:30) {
  set.seed(200 + r)
  m <- 40
  d <- data.frame(area = 1:m, x1 = rnorm(m), x2 = rnorm(m), x3 = rnorm(m), x4 = rnorm(m))
  d$n <- sample(20:100, m, replace = TRUE)
  p <- plogis(-1 + 0.5 * d$x1 - 0.4 * d$x2 + rnorm(m, sd = 0.4))
  d$y <- rbinom(m, d$n, p)
  fit <- q(hb_area(y ~ x1 + x2 + x3 + x4, data = d, family = "binomial", trials = "n", print_result = FALSE))
  for (cr in names(bin_sets)) bin_sets[[cr]][[r]] <- sel_terms(q(step_sae(fit, criterion = cr, trace = FALSE)))
}
res$bin_waic <- summarise(bin_sets$WAIC, "binomial HB, WAIC")
res$bin_lcpo <- summarise(bin_sets$LCPO, "binomial HB, LCPO")

sp_sets <- list(SFH = list(), FH = list())
for (r in 1:50) {
  d <- sim_area(r, spatial = TRUE)
  f_sfh <- q(eblup_sfh(y ~ x1 + x2 + x3 + x4, vardir = "vardir", data = d, W = W, print_result = FALSE))
  f_fh <- q(eblup_fh(y ~ x1 + x2 + x3 + x4, vardir = "vardir", data = d, print_result = FALSE))
  sp_sets$SFH[[r]] <- sel_terms(q(step_sae(f_sfh, criterion = "BIC", trace = FALSE)))
  sp_sets$FH[[r]] <- sel_terms(q(step_sae(f_fh, criterion = "BIC", trace = FALSE)))
}
res$sp_sfh <- summarise(sp_sets$SFH, "spatial truth, SFH, BIC")
res$sp_fh <- summarise(sp_sets$FH, "spatial truth, FH, BIC")
print(do.call(rbind, res), row.names = FALSE, digits = 3)
