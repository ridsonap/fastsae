# Repeated-sampling check of hb_area(var_x = ): coverage of the 95% intervals,
# RMSE and the coefficient of the covariate measured with error, for the naive
# Fay-Herriot model, the structural and functional measurement error models, and
# the Ybarra-Lohr EBLUP (emdi::fh(method = "me"), point estimates only).
# Run from the package root: Rscript benchmarks/validation/me_area_coverage.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
vdir <- "benchmarks/validation"
n_rep <- as.integer(Sys.getenv("ME_REPS", "100"))
D <- 60
n_out <- 10
beta_true <- 0.8

scenarios <- list(
  moderate = list(x = function(D) rnorm(D, 5, 2), C = function(D) runif(D, 0.2, 1.5)),
  heavy = list(x = function(D) rnorm(D, 5, 1), C = function(D) runif(D, 0.5, 2)),
  skewed = list(
    x = function(D) {
      z <- rlnorm(D, 0, 0.8)
      5 + 2 * (z - mean(z)) / stats::sd(z)
    },
    C = function(D) runif(D, 0.2, 1.5)
  ),
  # some values without error (e.g. from a census), others very noisy
  mixed = list(x = function(D) rnorm(D, 5, 2), C = function(D) sample(c(0, 0.1, 4), D, replace = TRUE))
)

one_rep <- function(scen, r) {
  INLA::inla.setOption(num.threads = "1:1")
  set.seed(1000 * match(scen, names(scenarios)) + r)
  x <- scenarios[[scen]]$x(D)
  C <- scenarios[[scen]]$C(D)
  d <- data.frame(area = seq_len(D), w = x + rnorm(D, 0, sqrt(C)), vx = C, vardir = runif(D, 0.5, 1.5))
  theta <- 1 + beta_true * x + rnorm(D, 0, 0.7)
  d$y <- theta + rnorm(D, 0, sqrt(d$vardir))
  d$y[sort(sample(D, n_out))] <- NA
  s <- !is.na(d$y)
  fit <- function(...) {
    suppressWarnings(suppressMessages(hb_area(y ~ w,
      data = d, domain = "area", vardir = "vardir", print_result = FALSE, ...
    )))
  }
  fits <- list(
    naive = fit(),
    structural = fit(var_x = c(w = "vx")),
    functional = fit(var_x = c(w = "vx"), me_model = "functional")
  )
  rows <- lapply(names(fits), function(m) {
    f <- fits[[m]]
    data.frame(
      scenario = scen, rep = r, model = m, area = d$area, sampled = s, C = C, theta = theta,
      est = f$df_hb$hb, lower = f$df_hb$ci_lower, upper = f$df_hb$ci_upper,
      beta = f$estcoef["w", "beta"], beta_lower = f$estcoef["w", "ci_lower"],
      beta_upper = f$estcoef["w", "ci_upper"]
    )
  })
  Ci <- array(0, c(2, 2, sum(s)))
  Ci[2, 2, ] <- d$vx[s]
  yl <- tryCatch(suppressMessages(suppressWarnings(emdi::fh(
    fixed = y ~ w, vardir = "vardir", combined_data = d[s, ], domains = "area",
    method = "me", Ci = Ci, MSE = FALSE
  ))), error = function(e) NULL)
  if (!is.null(yl)) {
    i <- match(d$area[s], as.numeric(as.character(yl$ind$Domain)))
    rows[[length(rows) + 1]] <- data.frame(
      scenario = scen, rep = r, model = "ybarra_lohr", area = d$area[s], sampled = TRUE, C = C[s],
      theta = theta[s], est = yl$ind$FH[i], lower = NA_real_, upper = NA_real_,
      beta = yl$model$coefficients[2, "coefficients"], beta_lower = NA_real_, beta_upper = NA_real_
    )
  }
  do.call(rbind, rows)
}

# Coverage and RMSE by model and scenario (sampled and unsampled areas), by size of
# the error variance, and the mean and coverage of the coefficient
summarise_coverage <- function(res) {
  res$covered <- res$theta >= res$lower & res$theta <= res$upper
  res$c_group <- cut(res$C, c(-Inf, 0, 0.75, Inf), labels = c("C = 0", "C <= 0.75", "C > 0.75"))
  summ <- function(df) {
    data.frame(coverage = mean(df$covered), rmse = sqrt(mean((df$est - df$theta)^2)), n = nrow(df))
  }
  by_group <- function(keys) {
    out <- do.call(rbind, lapply(split(res, res[keys], drop = TRUE), function(df) cbind(df[1, keys], summ(df))))
    out[do.call(order, out[keys]), ]
  }
  coef <- unique(res[c("scenario", "rep", "model", "beta", "beta_lower", "beta_upper")])
  beta_tab <- do.call(rbind, lapply(split(coef, coef[c("scenario", "model")], drop = TRUE), function(df) {
    data.frame(
      scenario = df$scenario[1], model = df$model[1], beta_mean = mean(df$beta),
      beta_coverage = mean(beta_true >= df$beta_lower & beta_true <= df$beta_upper)
    )
  }))
  list(
    by_model = by_group(c("scenario", "sampled", "model")),
    by_error_variance = by_group(c("scenario", "c_group", "model")),
    coefficient = beta_tab[order(beta_tab$scenario, beta_tab$model), ]
  )
}

cores <- max(1L, min(4L, parallel::detectCores() - 1L))
res <- do.call(rbind, lapply(names(scenarios), function(scen) {
  do.call(rbind, parallel::mclapply(seq_len(n_rep), function(r) one_rep(scen, r), mc.cores = cores))
}))
tabs <- summarise_coverage(res)
saveRDS(tabs, file.path(vdir, "results_me_coverage.rds"))
options(width = 160)
for (tab in tabs) print(tab, digits = 3, row.names = FALSE)
