#' @title Simulate Spatio-Temporal Multi-Distribution Small Area Data
#' @description Generates synthetic spatio-temporal panel datasets for small area estimation.
#'   Combines domain-level spatial autocorrelation (SAR) with domain-specific first-order
#'   autoregressive AR(1) temporal dynamics across \eqn{T} time periods, producing multi-distribution
#'   responses (Gaussian Fay-Herriot, Poisson, Binomial, Beta, Negative Binomial, Gamma).
#'
#' @param D Integer. Number of domains (areas). Default is 42.
#' @param T Integer. Number of time periods (\eqn{T \ge 2}). Default is 5.
#' @param time_start Integer. Starting year / time label. Default is 2022.
#' @param W Optional spatial proximity or adjacency matrix of dimension \eqn{D \times D}.
#'   If \code{NULL}, a spatial matrix is automatically simulated using \code{\link{sim_spatial_weights}}.
#' @param spatial_type Character. Topology for simulated \code{W} if \code{W = NULL}:
#'   \code{"knn"}, \code{"grid"}, or \code{"ring"}. Default is \code{"knn"}.
#' @param rho_s Numeric. Domain spatial autoregressive autocorrelation parameter (\eqn{|\rho_s| < 1},
#'   corresponding to \code{rho1} in \code{\link{eblup_stfh}}). Default is 0.5.
#' @param rho_t Numeric. Temporal AR(1) autocorrelation parameter (\eqn{|\rho_t| < 1},
#'   corresponding to \code{rho2} in \code{\link{eblup_stfh}}). Default is 0.6.
#' @param sigma_s Numeric. Standard deviation of domain spatial random effects (\eqn{u_1}). Default is 0.4.
#' @param sigma_t Numeric. Standard deviation of stationary spatio-temporal effects (\eqn{u_2}). Default is 0.3.
#' @param trend Character. Deterministric temporal trend type: \code{"linear"}, \code{"random_walk"},
#'   \code{"ar1"}, or \code{"none"}. Default is \code{"linear"}.
#' @param trend_slope Numeric. Slope / scale of the temporal trend. Default is 0.05.
#' @param beta Numeric vector of length 3 giving fixed regression coefficients \eqn{(\beta_0, \beta_1, \beta_2)}.
#'   Default is \code{c(1.0, 0.8, -0.5)}.
#' @param n_unsampled Integer. Number of domains that are completely unsampled across all time periods. Default is 4.
#' @param prop_intermittent Numeric. Proportion of domain-by-time observations to mark as intermittently missing
#'   (excluding persistent unsampled domains). Default is 0.05.
#' @param sort_order Character. Row sorting order: \code{"domain-major"} (all time periods for domain 1,
#'   then domain 2, required by \code{\link{eblup_stfh}}) or \code{"time-major"}. Default is \code{"domain-major"}.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return An object of class \code{c("fastsae_sim_series", "list")} containing:
#'   \describe{
#'     \item{data}{A data frame of dimension \eqn{(D \cdot T) \times 15} with columns:
#'       \code{area}, \code{year}, \code{x1}, \code{x2}, \code{y_gaussian}, \code{vardir},
#'       \code{y_poisson}, \code{exposure}, \code{y_binomial}, \code{trials}, \code{y_beta},
#'       \code{y_nbinomial}, \code{y_gamma}, \code{x_coord}, \code{y_coord}.}
#'     \item{W}{Binary adjacency matrix \eqn{D \times D}.}
#'     \item{W_std}{Row-standardized spatial weights matrix \eqn{D \times D}.}
#'     \item{coords}{Centroid coordinates of the \eqn{D} domains.}
#'     \item{u_spatial}{True domain spatial random effect vector of length \eqn{D}.}
#'     \item{u_temporal}{True spatio-temporal random effect matrix of dimension \eqn{D \times T}.}
#'     \item{parameters}{List of true simulation parameters (\code{rho_s}, \code{rho_t}, \code{sigma_s}, \code{sigma_t}, \code{beta}, \code{trend}).}
#'   }
#'
#' @examples
#' # 1. Simulate a 30-domain panel over 4 years
#' sim_panel <- sim_series_data(D = 30, T = 4, time_start = 2021, seed = 123)
#' print(sim_panel)
#' head(sim_panel$data)
#'
#' # 2. Directly estimate Spatio-Temporal Fay-Herriot model (eblup_stfh)
#' \dontrun{
#' fit_stfh <- eblup_stfh(
#'   y_gaussian ~ x1 + x2,
#'   data = sim_panel$data,
#'   domain = ~area,
#'   time = ~year,
#'   vardir = ~vardir,
#'   W = sim_panel$W_std
#' )
#' summary(fit_stfh)
#' }
#'
#' @export
sim_series_data <- function(D = 42,
                            T = 5,
                            time_start = 2022,
                            W = NULL,
                            spatial_type = c("knn", "grid", "ring"),
                            rho_s = 0.5,
                            rho_t = 0.6,
                            sigma_s = 0.4,
                            sigma_t = 0.3,
                            trend = c("linear", "random_walk", "ar1", "none"),
                            trend_slope = 0.05,
                            beta = c(1.0, 0.8, -0.5),
                            n_unsampled = 4,
                            prop_intermittent = 0.05,
                            sort_order = c("domain-major", "time-major"),
                            seed = NULL) {
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  T <- as.integer(T)
  if (is.na(T) || T < 2L) stop("'T' must be an integer >= 2.", call. = FALSE)

  time_start <- as.integer(time_start)
  time_vec <- time_start + seq_len(T) - 1L

  spatial_type <- match.arg(spatial_type)
  trend <- match.arg(trend)
  sort_order <- match.arg(sort_order)

  # 1. Process spatial matrix
  if (is.null(W)) {
    D <- as.integer(D)
    if (is.na(D) || D < 2L) stop("'D' must be an integer >= 2.", call. = FALSE)
    W <- sim_spatial_weights(D = D, type = spatial_type, style = "B")
    coords <- attr(W, "coords")
  } else {
    W <- as.matrix(W)
    D <- nrow(W)
    if (ncol(W) != D) stop("'W' must be a square matrix.", call. = FALSE)
    coords <- attr(W, "coords")
    if (is.null(coords)) {
      d_dist <- stats::as.dist(1 - (W / max(W + 1e-6)))
      cmd_res <- tryCatch(stats::cmdscale(d_dist, k = 2), error = function(e) NULL)
      if (!is.null(cmd_res) && nrow(cmd_res) == D && ncol(cmd_res) == 2L) {
        coords <- cmd_res
      } else {
        coords <- matrix(stats::runif(2L * D), ncol = 2L)
      }
      colnames(coords) <- c("x", "y")
    }
  }

  # Ensure domain names / area identifiers
  area_ids <- seq_len(D)
  domain_names <- rownames(W)
  if (is.null(domain_names)) {
    domain_names <- as.character(area_ids)
    rownames(W) <- colnames(W) <- domain_names
  }
  rownames(coords) <- domain_names

  # Row-standardized matrix for SAR process
  rs <- rowSums(W)
  rs[rs == 0] <- 1
  W_std <- W / rs

  # 2. Random Effects Simulation
  # a. Domain spatial effect u1 (SAR)
  rho_s <- max(-0.95, min(0.95, as.numeric(rho_s)))
  sigma_s <- max(0.01, as.numeric(sigma_s))
  I_D <- diag(D)
  A_sar <- I_D - rho_s * W_std
  u1_raw <- solve(A_sar, stats::rnorm(D))
  u1_sd <- stats::sd(u1_raw)
  u1 <- if (u1_sd > 0) sigma_s * (u1_raw - mean(u1_raw)) / u1_sd else u1_raw

  # b. Domain-specific temporal AR(1) effects u2 (D x T matrix)
  rho_t <- max(-0.95, min(0.95, as.numeric(rho_t)))
  sigma_t <- max(0.01, as.numeric(sigma_t))
  innov_sd <- sigma_t * sqrt(max(0.01, 1 - rho_t^2))

  u2 <- matrix(0, nrow = D, ncol = T)
  for (d in seq_len(D)) {
    u2[d, 1] <- stats::rnorm(1, mean = 0, sd = sigma_t)
    for (t in 2:T) {
      u2[d, t] <- rho_t * u2[d, t - 1] + stats::rnorm(1, mean = 0, sd = innov_sd)
    }
  }

  # c. Deterministic temporal trend
  trend_vec <- switch(
    trend,
    "linear" = trend_slope * (seq_len(T) - 1),
    "random_walk" = cumsum(c(0, stats::rnorm(T - 1, mean = 0, sd = trend_slope))),
    "ar1" = {
      rw <- numeric(T)
      for (t in 2:T) rw[t] <- 0.7 * rw[t - 1] + stats::rnorm(1, mean = 0, sd = trend_slope)
      rw
    },
    "none" = rep(0, T)
  )

  # 3. Covariates and Response Generation over Domain x Time
  # Generate baseline covariates per domain with dynamic temporal variations
  x1_base <- stats::rnorm(D, mean = 2, sd = 0.8)
  x2_base <- stats::runif(D, min = 1, max = 4)

  # Construct observation indices according to sort_order
  if (sort_order == "domain-major") {
    # All T periods for domain 1, then all T for domain 2, etc.
    grid_idx <- expand.grid(time_idx = seq_len(T), d_idx = seq_len(D))
  } else {
    # All D domains for time 1, then all D for time 2, etc.
    grid_idx <- expand.grid(d_idx = seq_len(D), time_idx = seq_len(T))
  }

  N <- nrow(grid_idx)
  d_vec <- grid_idx$d_idx
  t_vec <- grid_idx$time_idx

  area_col <- area_ids[d_vec]
  year_col <- time_vec[t_vec]

  # Dynamic covariates
  x1 <- round(x1_base[d_vec] + stats::rnorm(N, mean = 0, sd = 0.25), 4)
  x2 <- round(x2_base[d_vec] + stats::runif(N, min = -0.4, max = 0.4), 4)

  # Linear predictor eta_{d, t}
  eta_fixed <- beta[1] + beta[2] * x1 + beta[3] * x2 + trend_vec[t_vec]
  u1_rep <- u1[d_vec]
  u2_rep <- u2[cbind(d_vec, t_vec)]
  eta <- eta_fixed + u1_rep + u2_rep

  # 4. Generate Multi-Distribution Responses
  # a. Gaussian (Fay-Herriot)
  vardir <- round(stats::runif(N, min = 0.04, max = 0.16), 5)
  e_gauss <- stats::rnorm(N, mean = 0, sd = sqrt(vardir))
  y_gaussian <- round(eta + e_gauss, 4)

  # b. Poisson (counts with dynamic exposure)
  exposure <- round(stats::runif(N, min = 90, max = 420) * (1 + 0.02 * (t_vec - 1)))
  lambda_rate <- exp(pmin(pmax(eta / 2 - 0.5, -4), 4))
  y_poisson <- stats::rpois(N, lambda = exposure * lambda_rate)

  # c. Binomial (success counts with trials)
  trials <- round(stats::runif(N, min = 60, max = 280))
  prob_bin <- stats::plogis(eta - 1)
  y_binomial <- stats::rbinom(N, size = trials, prob = prob_bin)

  # d. Negative Binomial (overdispersed counts)
  theta_nb <- 3.0
  mu_nb <- exp(pmin(pmax(eta / 2, -3), 4)) * 15
  rate_gamma <- stats::rgamma(N, shape = theta_nb, scale = mu_nb / theta_nb)
  y_nbinomial <- stats::rpois(N, lambda = rate_gamma)

  # e. Beta (continuous proportion in (0, 1))
  mu_beta <- stats::plogis(eta - 1)
  mu_beta <- pmin(pmax(mu_beta, 0.02), 0.98)
  kappa_prec <- 25
  y_beta <- round(stats::rbeta(N, shape1 = mu_beta * kappa_prec, shape2 = (1 - mu_beta) * kappa_prec), 5)
  y_beta <- pmin(pmax(y_beta, 0.001), 0.999)

  # f. Gamma (skewed positive continuous)
  shape_gamma <- 5.0
  mu_gamma <- exp(pmin(pmax(eta / 2, -3), 3))
  y_gamma <- round(stats::rgamma(N, shape = shape_gamma, scale = mu_gamma / shape_gamma), 4)

  # 5. Handle Unsampled Domains & Intermittent Missingness
  # Persistent unsampled domains across all years
  n_unsampled <- max(0L, min(as.integer(n_unsampled), D - 1L))
  if (n_unsampled > 0L) {
    persistent_domains <- sample(area_ids, size = n_unsampled)
    pers_idx <- which(area_col %in% persistent_domains)
    y_gaussian[pers_idx] <- NA_real_
    y_poisson[pers_idx] <- NA_integer_
    y_binomial[pers_idx] <- NA_integer_
    y_nbinomial[pers_idx] <- NA_integer_
    y_beta[pers_idx] <- NA_real_
    y_gamma[pers_idx] <- NA_real_
  } else {
    persistent_domains <- integer(0)
  }

  # Intermittent missing observations in sampled domains
  prop_intermittent <- max(0, min(0.3, as.numeric(prop_intermittent)))
  sampled_indices <- which(!(area_col %in% persistent_domains))
  if (prop_intermittent > 0 && length(sampled_indices) > 0) {
    n_inter <- round(length(sampled_indices) * prop_intermittent)
    if (n_inter > 0) {
      inter_idx <- sample(sampled_indices, size = n_inter)
      y_gaussian[inter_idx] <- NA_real_
      y_poisson[inter_idx] <- NA_integer_
      y_binomial[inter_idx] <- NA_integer_
      y_nbinomial[inter_idx] <- NA_integer_
      y_beta[inter_idx] <- NA_real_
      y_gamma[inter_idx] <- NA_real_
    }
  }

  # Coordinates matching domain
  coords_x <- round(coords[d_vec, 1], 4)
  coords_y <- round(coords[d_vec, 2], 4)

  df <- data.frame(
    area = area_col,
    year = year_col,
    x1 = x1,
    x2 = x2,
    y_gaussian = y_gaussian,
    vardir = vardir,
    y_poisson = y_poisson,
    exposure = exposure,
    y_binomial = y_binomial,
    trials = trials,
    y_beta = y_beta,
    y_nbinomial = y_nbinomial,
    y_gamma = y_gamma,
    x_coord = coords_x,
    y_coord = coords_y,
    stringsAsFactors = FALSE
  )

  out <- list(
    data = df,
    W = W,
    W_std = W_std,
    coords = coords,
    u_spatial = u1,
    u_temporal = u2,
    parameters = list(
      D = D,
      T = T,
      rho_s = rho_s,
      rho_t = rho_t,
      sigma_s = sigma_s,
      sigma_t = sigma_t,
      beta = beta,
      trend = trend,
      trend_slope = trend_slope,
      sort_order = sort_order
    )
  )

  class(out) <- c("fastsae_sim_series", "list")
  return(out)
}

#' @export
print.fastsae_sim_series <- function(x, ...) {
  params <- x$parameters
  D <- params$D
  T <- params$T
  N <- nrow(x$data)
  n_na <- sum(is.na(x$data$y_gaussian))
  n_obs <- N - n_na

  cat("=================================================================\n")
  cat("  fastsae Simulated Spatio-Temporal Multi-Distribution Panel\n")
  cat("=================================================================\n")
  cat(sprintf("Domains (D): %d | Time periods (T): %d | Total observations: %d\n", D, T, N))
  cat(sprintf("Observations: %d observed, %d missing/unsampled (%.1f%%)\n",
              n_obs, n_na, 100 * n_na / N))
  cat(sprintf("Spatio-temporal parameters:\n"))
  cat(sprintf("  - Spatial correlation (rho_s / rho1):  %.2f (sigma_s = %.2f)\n", params$rho_s, params$sigma_s))
  cat(sprintf("  - Temporal AR(1)     (rho_t / rho2):  %.2f (sigma_t = %.2f)\n", params$rho_t, params$sigma_t))
  cat(sprintf("  - Temporal trend:                     %s (slope = %.3f)\n", params$trend, params$trend_slope))
  cat(sprintf("  - Row sort order:                     %s\n", params$sort_order))
  cat("Responses generated:\n")
  cat("  - Gaussian (FH):     y_gaussian (vardir)\n")
  cat("  - Poisson:           y_poisson (exposure)\n")
  cat("  - Binomial:          y_binomial (trials)\n")
  cat("  - Beta:              y_beta\n")
  cat("  - Negative Binomial: y_nbinomial\n")
  cat("  - Gamma:             y_gamma\n")
  cat(sprintf("Spatial proximity:     %d x %d (W binary, W_std row-standardized)\n", D, D))
  cat("=================================================================\n")
  invisible(x)
}
