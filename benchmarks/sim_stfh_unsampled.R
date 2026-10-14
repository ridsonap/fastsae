#!/usr/bin/env Rscript

# ==============================================================================
# Model-based Monte Carlo: eblup_stfh() with unsampled domains and missing
# periods.
#
# A fixed design (X, W, vardir and the pattern of missing cells) comes from
# sim_series_data(D = 50, T = 5, n_unsampled = 4, prop_intermittent = 0.10).
# In each replicate the area parameters are drawn from the STFH model
#   theta = X beta + u1 (SAR(rho1) on all D domains) + u2 (AR(1)(rho2) on all T),
# the direct estimates y = theta + e only in the sampled cells, and the model
# is refitted. Per cell type (sampled, gap in a sampled domain, unsampled domain)
# the script reports
#   * the empirical MSE of the EBLUP and of the synthetic estimator x'beta^;
#   * the relative bias of the analytical (and bootstrap) MSE estimator,
#     RB = mean(mse^) / MSE_emp - 1, averaged over the cells of that type;
#   * the coverage of theta^ -/+ 1.96 sqrt(mse^).
#
# Usage: Rscript benchmarks/sim_stfh_unsampled.R [--R=200] [--B=100] [--R_boot=50] [--cores=6]
# Output: benchmarks/output/sim_stfh_unsampled.rds and a printed summary.
# ==============================================================================

suppressPackageStartupMessages(library(fastsae))
args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else as.numeric(sub(paste0("^--", name, "="), "", hit[1]))
}
R <- arg_value("R", 200)
B <- arg_value("B", 100)
R_boot <- arg_value("R_boot", 50)       # replicates that also get the bootstrap MSE
cores <- arg_value("cores", max(1, parallel::detectCores() - 2))

# ---- fixed design --------------------------------------------------------------
D <- 50; T <- 5
sim <- sim_series_data(D = D, T = T, n_unsampled = 4, prop_intermittent = 0.10, seed = 2026)
d <- sim$data
d <- d[order(d$area, d$year), ]
W <- unclass(sim$W)
W <- W / rowSums(W)
attributes(W)[c("coords", "type", "style")] <- NULL
X <- cbind(1, d$x1, d$x2)
beta <- c(1, 0.8, -0.5)
pars <- c(sigma21 = 0.5, rho1 = 0.5, sigma22 = 0.2, rho2 = 0.5)
sampled <- !is.na(d$y_gaussian)
dom_never <- tapply(!sampled, d$area, all)
type <- ifelse(sampled, "sampled",
               ifelse(d$area %in% names(dom_never)[dom_never], "unsampled domain", "gap"))
vardir <- d$vardir

L1 <- t(chol(solve(crossprod(diag(D) - pars[["rho1"]] * W))))   # Omega1 = L1 L1'
draw_theta <- function() {
  u1 <- sqrt(pars[["sigma21"]]) * drop(L1 %*% rnorm(D))
  u2 <- unlist(lapply(seq_len(D), function(i) {
    x <- numeric(T)
    x[1] <- rnorm(1, 0, sqrt(pars[["sigma22"]] / (1 - pars[["rho2"]]^2)))
    for (t in 2:T) x[t] <- pars[["rho2"]] * x[t - 1] + rnorm(1, 0, sqrt(pars[["sigma22"]]))
    x
  }))
  drop(X %*% beta) + rep(u1, each = T) + u2
}

one_rep <- function(r) {
  set.seed(1000 + r)
  theta <- draw_theta()
  dd <- d
  dd$y <- ifelse(sampled, theta + rnorm(nrow(d), 0, sqrt(vardir)), NA_real_)
  dd$vd <- ifelse(sampled, vardir, NA_real_)
  fit <- function(...) suppressMessages(eblup_stfh(y ~ x1 + x2, data = dd, vardir = ~vd, domain = ~area,
                                                  time = ~year, W = W, compute_mse = TRUE,
                                                  print_result = FALSE, ...))
  m <- tryCatch(fit(), error = function(e) NULL)
  if (is.null(m) || !isTRUE(m$convergence)) return(NULL)
  out <- data.frame(rep = r, cell = seq_len(nrow(d)), type = type, theta = theta,
                    eblup = m$df_eblup$eblup, mse_an = m$df_eblup$mse,
                    synth = drop(X %*% m$estcoef$beta), mse_pb = NA_real_)
  if (r <= R_boot) {
    b <- tryCatch(fit(mse_type = "bootstrap", B = B, seed = r), error = function(e) NULL)
    if (!is.null(b)) out$mse_pb <- b$df_eblup$mse
  }
  out
}

t0 <- Sys.time()
res <- parallel::mclapply(seq_len(R), one_rep, mc.cores = cores)
res <- do.call(rbind, Filter(Negate(is.null), res))
cat(sprintf("%d of %d replicates converged; %.1f minutes\n",
            length(unique(res$rep)), R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- summaries per cell, then per cell type ----------------------------------
per_cell <- do.call(rbind, lapply(split(res, res$cell), function(x) {
  emp <- mean((x$eblup - x$theta)^2)
  xb <- x[!is.na(x$mse_pb), ]
  data.frame(
    cell = x$cell[1], type = x$type[1],
    mse_emp = emp,
    mse_synth = mean((x$synth - x$theta)^2),
    rb_an = mean(x$mse_an) / emp - 1,
    rb_pb = if (nrow(xb)) mean(xb$mse_pb) / mean((xb$eblup - xb$theta)^2) - 1 else NA_real_,
    cov_an = mean(abs(x$eblup - x$theta) <= 1.96 * sqrt(x$mse_an))
  )
}))
summ <- do.call(rbind, lapply(split(per_cell, per_cell$type), function(x) data.frame(
  type = x$type[1], cells = nrow(x),
  mse_eblup = mean(x$mse_emp), mse_synthetic = mean(x$mse_synth),
  rb_analytical_pct = 100 * mean(x$rb_an), rb_bootstrap_pct = 100 * mean(x$rb_pb),
  coverage_pct = 100 * mean(x$cov_an)
)))
summ <- summ[match(c("sampled", "gap", "unsampled domain"), summ$type), ]
rownames(summ) <- NULL
print(summ, digits = 3)

dir.create(file.path("benchmarks", "output"), showWarnings = FALSE)
saveRDS(list(summary = summ, per_cell = per_cell,
             settings = list(D = D, T = T, R = R, B = B, R_boot = R_boot, pars = pars, beta = beta,
                             date = format(Sys.Date()))),
        file.path("benchmarks", "output", "sim_stfh_unsampled.rds"))
