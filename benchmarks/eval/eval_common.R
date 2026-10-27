# ==============================================================================
# Shared pieces of the evaluation study (vignettes/articles/fastsae-evaluation.Rmd):
# scenario definitions, data-generating mechanisms and one fit per method.
# Sourced by benchmarks/eval/sim_sfh_accuracy.R, sim_beta_accuracy.R and
# bench_eval_resources.R. Uses the *installed* fastsae.
# ==============================================================================

# ---- maps ----------------------------------------------------------------------

# Regular grid of n cells (queen contiguity), as in bench_hb_area_engines.R
make_grid_sf <- function(n) {
  factors <- which(n %% 1:n == 0)
  nx <- factors[which.min(abs(factors - sqrt(n)))]
  ny <- n / nx
  if (nx > ny) {
    tmp <- nx; nx <- ny; ny <- tmp
  }
  poly <- sf::st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE)))
  grid <- sf::st_make_grid(poly, n = c(nx, ny))
  sf::st_sf(domain = as.character(seq_len(n)), geometry = grid[seq_len(n)])
}

# Irregular map: Voronoi tessellation of n random points in the unit square
make_voronoi_sf <- function(n, seed = 1) {
  set.seed(seed)
  pts <- sf::st_multipoint(cbind(stats::runif(n), stats::runif(n)))
  box <- sf::st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE)))
  vor <- sf::st_collection_extract(sf::st_voronoi(pts, sf::st_sfc(box)), "POLYGON")
  vor <- sf::st_intersection(vor, box)
  sf::st_sf(domain = as.character(seq_len(n)), geometry = vor)
}

# Binary adjacency matrix with domain names
adjacency_from_sf <- function(map) {
  W <- spdep::nb2mat(spdep::poly2nb(map), style = "B", zero.policy = TRUE)
  dimnames(W) <- list(map$domain, map$domain)
  W
}

# The map of a scenario: "grid" or "voronoi" (polygons, needed by tipsae) or "mys"
# (mys_proxmat, no polygons; SFH only)
scenario_map <- function(map, D) {
  if (map == "mys") {
    W <- get("mys_proxmat", envir = asNamespace("fastsae"))
    W <- as.matrix(W)
    dimnames(W) <- list(as.character(seq_len(nrow(W))), as.character(seq_len(nrow(W))))
    return(list(sf = NULL, W = (W > 0) * 1, W_std = as.matrix(W) / rowSums(W)))
  }
  m <- if (map == "grid") make_grid_sf(D) else make_voronoi_sf(D)
  W <- adjacency_from_sf(m)
  list(sf = m, W = W, W_std = W / rowSums(W))
}

# ---- SFH (SAR) data --------------------------------------------------------------
# y_d = x_d'beta + u_d + e_d, u = (I - rho W)^{-1} eps, eps ~ N(0, sigma2 I),
# e_d ~ N(0, psi_d). The design (x, psi) is fixed within a scenario; u and e are
# redrawn in every replicate. Truth: theta_d = x_d'beta + u_d.
sfh_design <- function(sc, seed) {
  mp <- scenario_map(sc$map, sc$D)
  D <- nrow(mp$W)
  set.seed(seed)
  x1 <- stats::rnorm(D)
  x2 <- stats::runif(D)
  # psi_d / sigma2 has mean `ratio` and ranges over [0.5, 1.5] x ratio
  psi <- sc$ratio * sc$sigma2 * stats::runif(D, 0.5, 1.5)
  list(D = D, W_std = mp$W_std, x1 = x1, x2 = x2, psi = psi, beta = c(1, 0.8, -0.5))
}

sfh_replicate <- function(des, sc, seed) {
  set.seed(seed)
  eps <- stats::rnorm(des$D, sd = sqrt(sc$sigma2))
  u <- solve(diag(des$D) - sc$rho * des$W_std, eps)
  theta <- des$beta[1] + des$beta[2] * des$x1 + des$beta[3] * des$x2 + u
  y <- theta + stats::rnorm(des$D, sd = sqrt(des$psi))
  data.frame(domain = seq_len(des$D), y = y, x1 = des$x1, x2 = des$x2, vardir = des$psi, theta = theta)
}

# One SFH fit. Returns estimates (est, mse) and parameters (beta, sigma2, rho)
fit_sfh <- function(method, dat, W) {
  fml <- y ~ x1 + x2
  if (method == "fastsae") {
    f <- fastsae::eblup_sfh(fml, vardir = "vardir", data = dat, W = W, method = "REML",
                            print_result = FALSE)
    return(list(est = f$df_eblup$eblup, mse = f$df_eblup$mse,
                beta = f$estcoef$beta, sigma2 = f$random_effect_var, rho = f$rho,
                converged = isTRUE(f$convergence)))
  }
  if (method == "sae") {
    f <- sae::mseSFH(fml, vardir, W, method = "REML", data = dat)
    return(list(est = as.numeric(f$est$eblup), mse = as.numeric(f$mse),
                beta = f$est$fit$estcoef$beta, sigma2 = f$est$fit$refvar,
                rho = f$est$fit$spatialcorr, converged = isTRUE(f$est$fit$convergence)))
  }
  if (method == "emdi") {
    f <- suppressMessages(emdi::fh(fixed = fml, vardir = "vardir", combined_data = dat,
                                   domains = "domain", method = "reml", correlation = "spatial",
                                   corMatrix = W, MSE = TRUE, mse_type = "analytical"))
    return(list(est = f$ind$FH, mse = f$MSE$FH,
                beta = f$model$coefficients$coefficients, sigma2 = f$model$variance[["variance"]],
                rho = f$model$variance[["correlation"]], converged = TRUE))
  }
  stop("unknown method ", method)
}

# ---- spatial Beta data ------------------------------------------------------------
# logit(mu_d) = x_d'beta + v_d + s_d, v ~ N(0, sigma_v^2), s = sigma_s * scaled ICAR
# (Besag, York & Mollie 1991; scaled to a geometric mean marginal variance of 1 as in
# Riebler et al. 2016 and in INLA/tipsae). Direct estimate y_d ~ Beta(mu_d kappa_d,
# (1 - mu_d) kappa_d); its sampling variance is estimated from y_d as
# y_d (1 - y_d) / (kappa_d + 1), as a survey would. Truth: mu_d.
beta_design <- function(sc, seed) {
  mp <- scenario_map(sc$map, sc$D)
  D <- nrow(mp$W)
  set.seed(seed)
  x1 <- stats::rnorm(D)
  x2 <- stats::runif(D, -1, 1)
  kappa <- sc$kappa * exp(stats::runif(D, -0.5, 0.5))
  Q <- diag(rowSums(mp$W)) - mp$W
  eg <- eigen(Q, symmetric = TRUE)
  pos <- eg$values > 1e-8
  gm <- exp(mean(log(diag(eg$vectors[, pos] %*% (t(eg$vectors[, pos]) / eg$values[pos])))))
  list(D = D, sf = mp$sf, W = mp$W, x1 = x1, x2 = x2, kappa = kappa,
       U = eg$vectors[, pos], lam = eg$values[pos], gm = gm, beta = c(-1, 0.5, -0.4))
}

beta_replicate <- function(des, sc, seed) {
  set.seed(seed)
  s <- as.vector(des$U %*% (stats::rnorm(length(des$lam)) / sqrt(des$lam))) / sqrt(des$gm)
  s <- sc$sigma_s * s
  v <- stats::rnorm(des$D, sd = sc$sigma_v)
  mu <- stats::plogis(des$beta[1] + des$beta[2] * des$x1 + des$beta[3] * des$x2 + v + s)
  y <- stats::rbeta(des$D, mu * des$kappa, (1 - mu) * des$kappa)
  y <- pmin(pmax(y, 1e-3), 1 - 1e-3)
  data.frame(domain = as.character(seq_len(des$D)), y = y, x1 = des$x1, x2 = des$x2,
             vardir = y * (1 - y) / (des$kappa + 1), mu = mu)
}

# One spatial Beta fit: fastsae (bym: iid + scaled ICAR, the structure of tipsae's
# spatial_error = TRUE) with INLA or NumPyro, or tipsae. Returns area estimates
# (mean, sd, 95% interval), coefficients and convergence diagnostics.
fit_beta <- function(method, dat, des, seed = 1) {
  if (method == "tipsae") {
    f <- tipsae::fit_sae(
      formula_fixed = y ~ x1 + x2, data = dat, domains = "domain",
      disp_direct = "vardir", type_disp = "var", likelihood = "beta",
      spatial_error = TRUE, spatial_df = des$sf, domains_spatial_df = "domain",
      temporal_error = FALSE, refresh = 0, seed = seed
    )
    sm <- rstan::summary(f$stanfit)$summary
    th <- sm[grep("^theta\\[", rownames(sm)), , drop = FALSE]
    est <- summary(f)$model_estimates
    ord <- match(dat$domain, as.character(est$Domains))
    bt <- sm[grep("^beta0\\[|^beta\\[", rownames(sm)), "mean"]
    sp <- rstan::get_sampler_params(f$stanfit, inc_warmup = FALSE)
    return(list(
      est = est$mean[ord], sd = est$sd[ord],
      lower = th[ord, "2.5%"], upper = th[ord, "97.5%"],
      beta = unname(bt), sigma_v = sm["sigma_v[1]", "50%"], sigma_s = sm["sigma_s[1]", "50%"],
      rhat_max = max(sm[, "Rhat"], na.rm = TRUE), ess_min = min(sm[, "n_eff"], na.rm = TRUE),
      divergences = sum(vapply(sp, function(x) sum(x[, "divergent__"]), 0))
    ))
  }
  eng <- if (method == "inla") fastsae::inla_engine() else
    fastsae::jax_engine(device = if (method == "metal") "metal" else "cpu", progress = FALSE, seed = seed)
  f <- fastsae::hb_area(y ~ x1 + x2, vardir = "vardir", data = dat, domain = "domain",
                        family = "beta", spatial = "bym", W = des$W, engine = eng,
                        print_result = FALSE)
  h <- f$df_hb[match(dat$domain, as.character(f$df_hb$domain)), ]
  # posterior medians of the two SDs (the median is invariant to the precision -> SD map)
  if (method == "inla") {
    hp <- f$hyperpar
    med <- function(pat) 1 / sqrt(hp[grep(pat, rownames(hp)), "0.5quant"])
    sig <- c(sigma_v = med("iid component"), sigma_s = med("spatial component"))
  } else {
    tr <- f$trace
    med <- function(p) if (!is.null(tr) && p %in% dimnames(tr)$parameter) stats::median(sqrt(tr[, , p])) else NA_real_
    sig <- c(sigma_v = med("sigma2_iid"), sigma_s = med("sigma2_spatial"))
  }
  out <- list(est = h$hb, sd = h$sd, lower = h$ci_lower, upper = h$ci_upper,
              beta = unname(f$estcoef$beta), sigma_v = sig[["sigma_v"]], sigma_s = sig[["sigma_s"]],
              rhat_max = NA_real_, ess_min = NA_real_, divergences = NA_real_)
  if (method != "inla") {
    out$rhat_max <- f$rhat_max %||% NA_real_
    out$ess_min <- f$ess_min %||% NA_real_
    out$divergences <- f$divergences %||% NA_real_
  }
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a
