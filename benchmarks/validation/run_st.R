# Validation of hb_area() space-time interactions (Knorr-Held types I-IV) against Stan.
# Run from the package root: Rscript benchmarks/validation/run_st.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
vdir <- "benchmarks/validation"
mod <- cmdstanr::cmdstan_model(file.path(vdir, "st_interaction.stan"), quiet = TRUE)
dir.create(file.path(vdir, "out"), showWarnings = FALSE)

D <- 25
Tt <- 6
W <- as.matrix(unclass(sim_spatial_weights(D, type = "grid", seed = 1)))
A <- (W > 0 | t(W) > 0) * 1
diag(A) <- 0
Rs <- as.matrix(INLA::inla.scale.model(diag(rowSums(A)) - A, constr = list(A = matrix(1, 1, D), e = 0)))
rw_struct <- function(k) {
  R <- crossprod(diff(diag(Tt), differences = k))
  Ac <- if (k == 1) matrix(1, 1, Tt) else rbind(rep(1, Tt), seq_len(Tt) - (Tt + 1) / 2)
  as.matrix(INLA::inla.scale.model(R, constr = list(A = Ac, e = rep(0, nrow(Ac)))))
}

simulate <- function(fam, seed) {
  set.seed(seed)
  d <- expand.grid(area = 1:D, year = 1:Tt)
  d$x1 <- rnorm(nrow(d))
  s_eff <- rnorm(D, sd = 0.4)[d$area]
  t_eff <- cumsum(rnorm(Tt, sd = 0.2))[d$year]
  i_eff <- rnorm(nrow(d), sd = 0.2)
  if (fam == "gaussian") {
    d$vardir <- runif(nrow(d), 0.3, 1)
    d$y <- 5 + d$x1 + s_eff + t_eff + i_eff + rnorm(nrow(d), sd = sqrt(d$vardir))
  } else {
    d$E <- runif(nrow(d), 5, 40)
    d$y <- rpois(nrow(d), d$E * exp(-1 + 0.4 * d$x1 + s_eff + t_eff + i_eff))
  }
  d$y[sample(nrow(d), 15)] <- NA
  d$y[d$area == 7] <- NA
  d
}

basis <- function(Q, tol = 1e-8) {
  e <- eigen((Q + t(Q)) / 2, symmetric = TRUE)
  pos <- e$values > tol * max(e$values)
  e$vectors[, pos, drop = FALSE] %*% diag(1 / sqrt(e$values[pos]), sum(pos))
}

interaction_spec <- function(type, temporal, d) {
  N <- nrow(d)
  if (type == "type1") {
    return(list(Nint = N, iint = seq_len(N), Q = diag(N), Aint = matrix(0, 0, N)))
  }
  if (type == "type3") {
    Q <- kronecker(diag(Tt), Rs)
    Aint <- t(sapply(seq_len(Tt), function(t) {
      a <- numeric(D * Tt)
      a[(t - 1) * D + 1:D] <- 1
      a
    }))
    return(list(Nint = D * Tt, iint = (d$year - 1) * D + d$area, Q = Q, Aint = Aint))
  }
  if (type == "type2" && temporal == "rw1") {
    Q <- kronecker(diag(D), rw_struct(1))
    Aint <- t(sapply(seq_len(D), function(dd) {
      a <- numeric(D * Tt)
      a[(dd - 1) * Tt + 1:Tt] <- 1
      a
    }))
    return(list(Nint = D * Tt, iint = (d$area - 1) * Tt + d$year, Q = Q, Aint = Aint))
  }
  g <- .st_generic0(type, temporal, domain_id = d$area, time_id = d$year, n_domains = D, n_times = Tt, graph = A)
  list(Nint = D * Tt, iint = g$index, Q = as.matrix(g$Cmatrix), Aint = g$extraconstr$A)
}

pc1 <- list(prior = "pc.prec", param = c(1, 0.01))
compare <- function(type, temporal, fam, seed) {
  d <- simulate(fam, seed)
  s <- !is.na(d$y)
  fits <- list()
  times <- c()
  for (mode in c("default", "compact", "classic")) {
    args <- list(
      formula = y ~ x1, data = d, domain = "area", time = "year", family = fam, spatial = "besag",
      temporal = temporal, st_interaction = type, W = W,
      # the PC priors of st_interaction.stan, not the half-normal default of fastsae >= 1.0.0
      engine = do.call(inla_engine, c(
        list(prior_prec = pc1, prior_prec_time = pc1),
        if (mode != "default") list(inla.mode = mode)
      )),
      print_result = FALSE
    )
    if (fam == "gaussian") args$vardir <- "vardir" else args$exposure <- "E"
    t0 <- Sys.time()
    fits[[mode]] <- suppressWarnings(do.call(hb_area, args))
    times[mode] <- as.numeric(Sys.time() - t0, units = "secs")
  }
  it <- interaction_spec(type, temporal, d)
  Bs <- basis(Rs)
  Bt <- basis(rw_struct(if (temporal == "rw2") 2 else 1))
  Bi <- basis(it$Q)
  stopifnot(ncol(Bi) == it$Nint - nrow(it$Aint)) # constraints span the null space
  sd_ <- list(
    N = nrow(d), p = 2, X = cbind(1, d$x1), D = D, Tt = Tt, dom = d$area, tim = d$year,
    Ns = sum(s), idx = which(s), fam = if (fam == "gaussian") 1L else 2L,
    y = if (fam == "gaussian") d$y[s] else rep(0, sum(s)),
    yi = if (fam == "gaussian") rep(0L, sum(s)) else as.integer(d$y[s]),
    psi = if (fam == "gaussian") d$vardir[s] else rep(1, sum(s)),
    logE = if (fam == "gaussian") rep(0, sum(s)) else log(d$E[s]),
    rs = ncol(Bs), Bs = Bs, rt = ncol(Bt), Bt = Bt, trend = as.integer(temporal == "rw2"),
    tc = seq_len(Tt) - (Tt + 1) / 2, Nint = it$Nint, iint = it$iint, ri = ncol(Bi), Bi = Bi,
    lam = -log(0.01)
  )
  fit <- mod$sample(
    data = sd_, chains = 4, parallel_chains = 4, iter_warmup = 2000, iter_sampling = 4000,
    refresh = 0, show_messages = FALSE, seed = 1, adapt_delta = 0.99,
    output_dir = file.path(vdir, "out")
  )
  dr <- fit$draws("theta", format = "matrix")
  sm <- colMeans(dr)
  ss <- apply(dr, 2, sd)
  div <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
  ess <- min(posterior::summarise_draws(dr, "ess_bulk")$ess_bulk)
  res <- do.call(rbind, lapply(names(fits), function(mode) {
    f <- fits[[mode]]
    z <- (f$df_hb$hb - sm) / ss
    r <- f$df_hb$sd / ss
    data.frame(
      config = paste(type, temporal, fam), seed = seed, mode = mode, max_dmean_sd = max(abs(z)),
      mean_dmean_sd = mean(abs(z)), sd_ratio_min = min(r), sd_ratio_max = max(r),
      sd_ratio_dom7 = mean(r[d$area == 7]), time_s = times[[mode]], divergent = div, ess_min = ess
    )
  }))
  print(format(res, digits = 3), row.names = FALSE)
  res
}

cfgs <- list(
  c("type1", "rw1"), c("type2", "rw1"), c("type3", "rw1"), c("type4", "rw1"),
  c("type2", "rw2"), c("type4", "rw2")
)
res <- do.call(rbind, c(
  lapply(cfgs, function(cf) compare(cf[1], cf[2], "gaussian", 1)),
  lapply(cfgs[1:4], function(cf) compare(cf[1], cf[2], "poisson", 1))
))
saveRDS(res, file.path(vdir, "results_st.rds"))
