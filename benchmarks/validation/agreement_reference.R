#!/usr/bin/env Rscript

# ==============================================================================
# Agreement of the fastsae EBLUP functions with sae and emdi on the same data.
#
# For each model and fitting method, the largest relative difference
# max |a - b| / max |b| between fastsae (a) and the reference package (b) of the
# regression coefficients, the variance parameters, the EBLUPs and the analytical
# MSEs. Writes inst/extdata/agreement_reference.rds (read by the R Journal article).
#
# Usage (from the package root): Rscript benchmarks/validation/agreement_reference.R
# ==============================================================================

suppressPackageStartupMessages({
  library(fastsae)
  library(sae)
  library(emdi)
})

rel <- function(a, b) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  if (length(a) != length(b) || anyNA(a) || anyNA(b)) return(NA_real_)
  max(abs(a - b)) / max(abs(b))
}
row <- function(model, data, method, ref, beta, var, eblup, mse, note = "") {
  data.frame(Model = model, Data = data, Method = method, Reference = ref,
             Coefficients = beta, Variances = var, EBLUP = eblup, MSE = mse, Note = note)
}

# Simulated area-level data as in benchmarks/bench_frequentist_child.R
p_cov <- 5L
T_periods <- 5L
sim_fh <- function(D, W = NULL, rho = 0, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(D * p_cov), D, p_cov, dimnames = list(NULL, paste0("x", 1:p_cov)))
  beta <- c(1, seq(0.8, -0.8, length.out = p_cov))
  u <- rnorm(D, sd = 0.7)
  v <- if (is.null(W)) u else solve(diag(D) - rho * W, u)
  vardir <- runif(D, 0.3, 1.2)
  data.frame(area = seq_len(D), y = drop(cbind(1, X) %*% beta) + v + rnorm(D, sd = sqrt(vardir)),
             X, vardir = vardir)
}
sim_stfh <- function(D, W, rho1 = 0.5, rho2 = 0.6, seed = 1) {
  set.seed(seed)
  M <- D * T_periods
  X <- matrix(rnorm(M * p_cov), M, p_cov, dimnames = list(NULL, paste0("x", 1:p_cov)))
  beta <- c(1, seq(0.8, -0.8, length.out = p_cov))
  v1 <- solve(diag(D) - rho1 * W, rnorm(D, sd = 0.6))
  v2 <- as.vector(sapply(seq_len(D), function(d) {
    e <- numeric(T_periods)
    e[1] <- rnorm(1, sd = 0.4 / sqrt(1 - rho2^2))
    for (t in 2:T_periods) e[t] <- rho2 * e[t - 1] + rnorm(1, sd = 0.4)
    e
  }))
  vardir <- runif(M, 0.3, 1.2)
  data.frame(
    area = rep(seq_len(D), each = T_periods), year = rep(seq_len(T_periods), D),
    y = drop(cbind(1, X) %*% beta) + rep(v1, each = T_periods) + v2 + rnorm(M, sd = sqrt(vardir)),
    X, vardir = vardir
  )
}
fml5 <- y ~ x1 + x2 + x3 + x4 + x5
rows <- list()

# ---- Fay-Herriot: sampled areas of mys ---------------------------------------
mys_s <- as.data.frame(na.omit(mys))
for (m in c("REML", "ML")) {
  f <- eblup_fh(y ~ x1 + x2 + x3, vardir = "vardir", domain = "area", data = mys_s,
                method = m, print_result = FALSE)
  s <- sae::mseFH(y ~ x1 + x2 + x3, vardir, method = m, data = mys_s)
  rows[[length(rows) + 1]] <- row("Fay-Herriot", "mys", m, "sae",
    rel(f$estcoef$beta, s$est$fit$estcoef$beta), rel(f$random_effect_var, s$est$fit$refvar),
    rel(f$df_eblup$eblup, s$est$eblup), rel(f$df_eblup$mse, s$mse))
  e <- emdi::fh(fixed = y ~ x1 + x2 + x3, vardir = "vardir", combined_data = mys_s,
                domains = "area", method = tolower(m), MSE = TRUE)
  o <- match(f$df_eblup$domain, e$ind$Domain)
  rows[[length(rows) + 1]] <- row("Fay-Herriot", "mys", m, "emdi",
    rel(f$estcoef$beta, e$model$coefficients$coefficients), rel(f$random_effect_var, e$model$variance),
    rel(f$df_eblup$eblup, e$ind$FH[o]), rel(f$df_eblup$mse, e$MSE$FH[o]))
}

# ---- Spatial Fay-Herriot: simulated, D = 100, rho = 0.5 ----------------------
D <- 100
W <- as.matrix(fastsae:::sim_spatial_weights(D = D, type = "knn", style = "W", k = 4, seed = D))
d_sfh <- sim_fh(D, W = W, rho = 0.5, seed = D)
for (m in c("REML", "ML")) {
  f <- eblup_sfh(fml5, vardir = "vardir", domain = "area", data = d_sfh, W = W, method = m,
                 print_result = FALSE)
  s <- sae::mseSFH(fml5, vardir, W, method = m, data = d_sfh)
  rows[[length(rows) + 1]] <- row("Spatial Fay-Herriot", "simulated, D = 100", m, "sae",
    rel(f$estcoef$beta, s$est$fit$estcoef$beta),
    rel(c(f$random_effect_var, f$rho), c(s$est$fit$refvar, s$est$fit$spatialcorr)),
    rel(f$df_eblup$eblup, s$est$eblup), rel(f$df_eblup$mse, s$mse),
    if (m == "ML") "ML information in g3 and the bias term (deliberate)" else "")
  e <- emdi::fh(fixed = fml5, vardir = "vardir", combined_data = d_sfh, domains = "area",
                method = tolower(m), MSE = TRUE, correlation = "spatial", corMatrix = W)
  o <- match(f$df_eblup$domain, e$ind$Domain)
  rows[[length(rows) + 1]] <- row("Spatial Fay-Herriot", "simulated, D = 100", m, "emdi",
    rel(f$estcoef$beta, e$model$coefficients$coefficients),
    rel(c(f$random_effect_var, f$rho), c(e$model$variance$variance, e$model$variance$correlation)),
    rel(f$df_eblup$eblup, e$ind$FH[o]), rel(f$df_eblup$mse, e$MSE$FH[o]),
    if (m == "ML") "ML information in g3 and the bias term (deliberate)" else "")
}

# ---- Spatio-temporal Fay-Herriot: simulated, D = 100, T = 5 (REML only) -------
d_st <- sim_stfh(D, W, seed = D)
f <- eblup_stfh(fml5, data = d_st, vardir = ~vardir, domain = ~area, time = ~year, W = W,
                print_result = FALSE)
s <- sae::eblupSTFH(fml5, D = D, T = T_periods, vardir = vardir, proxmat = W, data = d_st)
stopifnot(isTRUE(s$fit$convergence))
rows[[length(rows) + 1]] <- row("Spatio-temporal Fay-Herriot", "simulated, D = 100, T = 5", "REML", "sae",
  rel(f$estcoef$beta, s$fit$estcoef$beta), rel(f$estvarcomp$estimate, s$fit$estvarcomp$estimate),
  rel(f$df_eblup$eblup, s$eblup), NA_real_, "sae has no analytical MSE for this model")

# ---- Battese-Harter-Fuller: cornsoybean (REML) --------------------------------
data(cornsoybean, package = "sae")
data(cornsoybeanmeans, package = "sae")
unit <- cornsoybean
names(unit)[names(unit) == "County"] <- "CountyIndex"
Xpop <- data.frame(CountyIndex = cornsoybeanmeans$CountyIndex,
                   CornPix = cornsoybeanmeans$MeanCornPixPerSeg,
                   SoyBeansPix = cornsoybeanmeans$MeanSoyBeansPixPerSeg,
                   PopnSegments = cornsoybeanmeans$PopnSegments)
f <- eblup_bhf(CornHec ~ CornPix + SoyBeansPix, unit_data = unit, Xpop = Xpop,
               domain_var = "CountyIndex", popsize_var = "PopnSegments", method = "REML",
               mse_type = "none", print_result = FALSE)
s <- sae::eblupBHF(CornHec ~ CornPix + SoyBeansPix, dom = County,
                   selectdom = unique(cornsoybean$County),
                   meanxpop = cornsoybeanmeans[, c("CountyIndex", "MeanCornPixPerSeg", "MeanSoyBeansPixPerSeg")],
                   popnsize = cornsoybeanmeans[, c("CountyIndex", "PopnSegments")],
                   method = "REML", data = cornsoybean)
fe <- f$df_eblup[order(f$df_eblup$domain), ]
se <- s$eblup[order(s$eblup$domain), ]
rows[[length(rows) + 1]] <- row("Battese-Harter-Fuller", "cornsoybean", "REML", "sae",
  rel(f$fit$beta, s$fit$fixed), rel(c(f$fit$random_effect_var, f$fit$sigma2_e), c(s$fit$refvar, s$fit$errorvar)),
  rel(fe$eblup, se$eblup), NA_real_, "MSE by bootstrap in both packages")

agreement <- do.call(rbind, rows)
print(agreement, digits = 3)
saveRDS(list(agreement = agreement, date = Sys.Date(),
             versions = vapply(c("fastsae", "sae", "emdi"), function(p) as.character(packageVersion(p)), "")),
        "inst/extdata/agreement_reference.rds")
