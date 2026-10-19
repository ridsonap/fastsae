# Default prior on the area SD of non-Gaussian hb_area() (.default_prior_prec), beyond
# default_sd_prior_glm.R: PC(1, 0.01) vs half-normal(1) (the default) vs half-normal(2.5)
# on the link scale, for Poisson, negative binomial (size 10) and gamma (shape 20, no
# vardir) responses with exposure 10-100, binomial and beta responses from small samples
# (n_d 5-20) and a binomial BYM2 model (4-nearest-neighbour graph, phi = 0.6).
# m = 15/30/60 areas (30/60 for BYM2; 1/8 without a direct estimate), area SD 0.3/0.8
# (and 1.5 for the counts and gamma), 100 replicates. ~3.5 hours on 6 cores.
#
# Run from the package root: Rscript benchmarks/validation/default_sd_prior_glm_more.R
#
# Result (2026-10-19), mean coverage of the 95% intervals over the cells (worst cell):
#   scenario             HN(1)          HN(2.5)        PC(1, 0.01)
#   poisson              0.947 (0.944)  0.948 (0.944)  0.943 (0.934)
#   binomial BYM2        0.947 (0.941)  0.947 (0.942)  0.942 (0.934)
#   binomial, n_d 5-20   0.950 (0.932)  0.953 (0.938)  0.919 (0.864)
#   beta, n_d 5-20       0.938 (0.925)  0.940 (0.929)  0.908 (0.873)
#   negative binomial    0.906 (0.874)  0.902 (0.870)  0.907 (0.880)
#   gamma, no vardir     0.855 (0.823)  0.849 (0.815)  0.877 (0.857)
# No fit failed. The median relative error was equal for the three priors (gamma and
# negative binomial: 0.159 and 0.245 with HN(1), 0.167 and 0.264 with PC). With small
# samples the PC prior estimated an area SD of 0.8 as 0.45 (binomial, 15 areas). The
# negative binomial and gamma intervals are too short with every prior: with one
# observation per area the dispersion and the area effect are only weakly identified.
devtools::load_all(quiet = TRUE)
q <- function(e) suppressWarnings(suppressMessages(e))
priors <- list(
  "PC(1, 0.01)" = list(prior = "pc.prec", param = c(1, 0.01)),
  "HN(1)" = .half_normal_prec(1),
  "HN(2.5)" = .half_normal_prec(2.5)
)
R <- 100

knn_graph <- function(m, seed) {
  set.seed(seed)
  xy <- cbind(runif(m), runif(m))
  dd <- as.matrix(stats::dist(xy))
  diag(dd) <- Inf
  A <- matrix(0, m, m)
  for (i in seq_len(m)) A[i, order(dd[i, ])[1:4]] <- 1
  pmax(A, t(A))
}
graphs <- list("30" = knn_graph(30, 1), "60" = knn_graph(60, 2))
spec <- lapply(graphs, .icar_spectral)

one_rep <- function(scn, m, su, r) {
  set.seed(1e6 * match(scn, scenarios) + 1e3 * m + 100 * su + r)
  x <- rnorm(m)
  miss <- seq_len(max(2, round(m / 8)))
  args <- list(formula = y ~ x, print_result = FALSE)
  if (scn == "bym2") {
    sp <- spec[[as.character(m)]]
    pos <- sp$eig_values > 1e-6
    s_icar <- drop(sp$eig_vectors[, pos] %*% (rnorm(sum(pos)) / sqrt(sp$eig_values[pos] * sp$scale_factor)))
    u <- su * (sqrt(0.4) * rnorm(m) + sqrt(0.6) * s_icar)
    args$spatial <- "bym2"
    args$W <- graphs[[as.character(m)]]
  } else {
    u <- rnorm(m, sd = su)
  }
  if (scn %in% c("poisson", "nbinomial", "gamma")) {
    theta <- exp(-0.5 + 0.5 * x + u)
    E <- runif(m, 10, 100)
    y <- switch(scn,
      poisson = rpois(m, E * theta),
      nbinomial = rnbinom(m, mu = E * theta, size = 10),
      gamma = rgamma(m, shape = 20, rate = 20 / theta)
    )
    if (scn != "gamma") args$exposure <- "E"
    fam <- scn
  } else {
    theta <- plogis(-1 + 0.5 * x + u)
    n <- if (scn %in% c("binomial_small", "beta_small")) sample(5:20, m, replace = TRUE) else sample(20:200, m, replace = TRUE)
    fam <- if (startsWith(scn, "beta")) "beta" else "binomial"
    y <- if (fam == "binomial") rbinom(m, n, theta) else pmin(pmax(rbeta(m, theta * (n - 1), (1 - theta) * (n - 1)), 1e-3), 1 - 1e-3)
    args$trials <- "n"
    E <- n
  }
  y[miss] <- NA
  d <- data.frame(y = y, x = x, E = E, n = E)
  args$data <- d
  args$family <- fam
  out <- NULL
  for (nm in names(priors)) {
    f <- tryCatch(q(do.call(hb_area, c(args, list(
      engine = inla_engine(prior_prec = priors[[nm]], num.threads = "1:1")
    )))), error = function(e) NULL)
    if (is.null(f)) {
      out <- rbind(out, data.frame(
        scenario = scn, m = m, sigma_u = su, prior = nm, ok = FALSE,
        cover = NA, w = NA, rel_err2 = NA, s2u = NA
      ))
      next
    }
    h <- f$df_hb
    out <- rbind(out, data.frame(
      scenario = scn, m = m, sigma_u = su, prior = nm, ok = TRUE,
      cover = theta >= h$ci_lower & theta <= h$ci_upper, w = (h$ci_upper - h$ci_lower) / theta,
      rel_err2 = ((h$hb - theta) / theta)^2, s2u = f$random_effect_var
    ))
  }
  out
}
scenarios <- c("poisson", "nbinomial", "gamma", "binomial_small", "beta_small", "bym2")
grid <- rbind(
  expand.grid(r = seq_len(R), su = c(0.3, 0.8, 1.5), m = c(15, 30, 60), scn = c("poisson", "nbinomial", "gamma"), stringsAsFactors = FALSE),
  expand.grid(r = seq_len(R), su = c(0.3, 0.8), m = c(15, 30, 60), scn = c("binomial_small", "beta_small"), stringsAsFactors = FALSE),
  expand.grid(r = seq_len(R), su = c(0.3, 0.8), m = c(30, 60), scn = "bym2", stringsAsFactors = FALSE)
)
t0 <- Sys.time()
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  tryCatch(one_rep(g$scn, g$m, g$su, g$r), error = function(e) NULL)
}, mc.cores = 6))
ok <- res[res$ok, ]
cat("failed fits by prior:\n")
print(tapply(!res$ok, res$prior, sum))
a <- aggregate(cbind(cover, w, rel_err2) ~ scenario + m + sigma_u + prior, ok, mean)
a$rel_rmse <- sqrt(a$rel_err2)
a$rel_err2 <- NULL
s <- aggregate(s2u ~ scenario + m + sigma_u + prior, unique(ok[, c("scenario", "m", "sigma_u", "prior", "s2u")]), function(v) median(sqrt(v)))
a <- merge(a, setNames(s, c("scenario", "m", "sigma_u", "prior", "sigma_u_hat")))
print(a[order(a$scenario, a$sigma_u, a$m, a$prior), ], digits = 3, row.names = FALSE)
cat("elapsed:", format(Sys.time() - t0), "\n")
