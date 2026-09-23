#' @title Simulate Multi-Distribution Small Area Data with Spatial Autocorrelation
#' @description Generates synthetic area-level datasets for small area estimation,
#'   including continuous (Gaussian, Gamma), count (Poisson, Negative Binomial),
#'   and proportion/rate responses (Binomial, Beta), driven by shared covariates
#'   and structured spatial random effects (BYM2 / SAR).
#'
#' @param D Integer. Number of domains (areas). Default is 42.
#' @param W Optional spatial proximity or adjacency matrix of dimension \eqn{D \times D}.
#'   If \code{NULL}, a spatial matrix is automatically simulated using \code{\link{sim_spatial_weights}}.
#' @param spatial_type Character. Topology for simulated \code{W} if \code{W = NULL}:
#'   \code{"knn"}, \code{"grid"}, or \code{"ring"}. Default is \code{"knn"}.
#' @param rho Numeric. Spatial autocorrelation autoregressive parameter (\eqn{|\rho| < 1}). Default is 0.5.
#' @param phi Numeric. Proportion of spatial marginal variance relative to total area variance (BYM2 parameter, \eqn{0 \le \phi \le 1}). Default is 0.6.
#' @param sigma_u Numeric. Overall standard deviation of domain random effects. Default is 0.5.
#' @param beta Numeric vector of length 3 giving true regression coefficients \eqn{(\beta_0, \beta_1, \beta_2)}.
#'   Default is \code{c(1.0, 0.8, -0.5)}.
#' @param n_unsampled Integer. Number of domains to mark as unsampled (where response variables are set to \code{NA}). Default is 6.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return An object of class \code{c("fastsae_sim_data", "list")} containing:
#'   \describe{
#'     \item{data}{A data frame with columns: \code{domain}, \code{x1}, \code{x2},
#'       \code{y_gaussian}, \code{vardir}, \code{y_poisson}, \code{exposure},
#'       \code{y_binomial}, \code{trials}, \code{y_beta}, \code{y_nbinomial},
#'       \code{y_gamma}, \code{x_coord}, \code{y_coord}.}
#'     \item{W}{Binary adjacency matrix \eqn{D \times D}.}
#'     \item{W_std}{Row-standardized proximity matrix \eqn{D \times D}.}
#'     \item{coords}{Coordinates of domain centroids.}
#'     \item{u_spatial}{True structured spatial random effect vector.}
#'     \item{u_iid}{True unstructured random effect vector.}
#'     \item{u_total}{True total random effect vector.}
#'     \item{phi}{True spatial variance fraction.}
#'     \item{rho}{True spatial autoregressive parameter.}
#'   }
#'
#' @examples
#' # 1. Simulate dataset with default KNN spatial structure
#' sim <- sim_area_data(D = 40, spatial_type = "knn", seed = 123)
#' print(sim)
#' head(sim$data)
#'
#' # 2. Using an existing spatial matrix
#' my_W <- sim_spatial_weights(D = 30, type = "grid")
#' sim2 <- sim_area_data(W = my_W, seed = 456)
#'
#' @export
sim_area_data <- function(D = 42,
                          W = NULL,
                          spatial_type = c("knn", "grid", "ring"),
                          rho = 0.5,
                          phi = 0.6,
                          sigma_u = 0.5,
                          beta = c(1.0, 0.8, -0.5),
                          n_unsampled = 6,
                          seed = NULL) {
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  spatial_type <- match.arg(spatial_type)

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
      # Fallback 2D coordinates via MDS or random
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

  # Ensure domain names
  domain_names <- rownames(W)
  if (is.null(domain_names)) {
    domain_names <- paste0("Area_", sprintf("%02d", seq_len(D)))
    rownames(W) <- colnames(W) <- domain_names
  }
  rownames(coords) <- domain_names

  # Row-standardized matrix for SAR process
  rs <- rowSums(W)
  rs[rs == 0] <- 1
  W_std <- W / rs

  # 2. Covariates and Linear Predictor
  x1 <- stats::rnorm(D, mean = 2, sd = 1)
  x2 <- stats::runif(D, min = 0, max = 5)
  eta_fix <- beta[1] + beta[2] * x1 + beta[3] * x2

  # 3. Spatial and IID Random Effects (SAR + BYM2 combination)
  rho <- max(-0.95, min(0.95, as.numeric(rho)))
  phi <- max(0, min(1, as.numeric(phi)))
  sigma_u <- max(0.01, as.numeric(sigma_u))

  I_D <- diag(D)
  A_sar <- I_D - rho * W_std
  v_raw <- solve(A_sar, stats::rnorm(D))
  v_sd <- stats::sd(v_raw)
  v <- if (v_sd > 0) (v_raw - mean(v_raw)) / v_sd else v_raw

  u_raw <- stats::rnorm(D)
  u_sd <- stats::sd(u_raw)
  u_iid <- if (u_sd > 0) (u_raw - mean(u_raw)) / u_sd else u_raw

  u_total <- sigma_u * (sqrt(phi) * v + sqrt(1 - phi) * u_iid)
  eta <- eta_fix + u_total

  # 4. Generate Responses
  # a. Gaussian (Fay-Herriot)
  vardir <- stats::runif(D, min = 0.04, max = 0.16)
  e_gauss <- stats::rnorm(D, mean = 0, sd = sqrt(vardir))
  y_gaussian <- eta + e_gauss

  # b. Poisson (count response with exposure offset)
  exposure <- round(stats::runif(D, min = 80, max = 400))
  lambda_rate <- exp(pmin(pmax(eta / 2 - 0.5, -4), 4))
  y_poisson <- stats::rpois(D, lambda = exposure * lambda_rate)

  # c. Binomial (success count with trials)
  trials <- round(stats::runif(D, min = 50, max = 250))
  prob_bin <- stats::plogis(eta - 1)
  y_binomial <- stats::rbinom(D, size = trials, prob = prob_bin)

  # d. Negative Binomial (overdispersed count)
  theta_nb <- 3.0
  mu_nb <- exp(pmin(pmax(eta / 2, -3), 4)) * 15
  rate_gamma <- stats::rgamma(D, shape = theta_nb, scale = mu_nb / theta_nb)
  y_nbinomial <- stats::rpois(D, lambda = rate_gamma)

  # e. Beta (continuous proportion in (0, 1))
  mu_beta <- stats::plogis(eta - 1)
  mu_beta <- pmin(pmax(mu_beta, 0.02), 0.98)
  kappa_prec <- 25
  y_beta <- stats::rbeta(D, shape1 = mu_beta * kappa_prec, shape2 = (1 - mu_beta) * kappa_prec)
  y_beta <- pmin(pmax(y_beta, 0.001), 0.999)

  # f. Gamma (skewed positive continuous)
  shape_gamma <- 5.0
  mu_gamma <- exp(pmin(pmax(eta / 2, -3), 3))
  y_gamma <- stats::rgamma(D, shape = shape_gamma, scale = mu_gamma / shape_gamma)

  # 5. Handle Unsampled Areas
  n_unsampled <- max(0L, min(as.integer(n_unsampled), D - 1L))
  if (n_unsampled > 0L) {
    unsampled_idx <- sample(seq_len(D), size = n_unsampled)
    y_gaussian[unsampled_idx] <- NA_real_
    y_poisson[unsampled_idx] <- NA_integer_
    y_binomial[unsampled_idx] <- NA_integer_
    y_nbinomial[unsampled_idx] <- NA_integer_
    y_beta[unsampled_idx] <- NA_real_
    y_gamma[unsampled_idx] <- NA_real_
  }

  df <- data.frame(
    domain = domain_names,
    x1 = round(x1, 4),
    x2 = round(x2, 4),
    y_gaussian = round(y_gaussian, 4),
    vardir = round(vardir, 5),
    y_poisson = y_poisson,
    exposure = exposure,
    y_binomial = y_binomial,
    trials = trials,
    y_beta = round(y_beta, 5),
    y_nbinomial = y_nbinomial,
    y_gamma = round(y_gamma, 4),
    x_coord = round(coords[, 1], 4),
    y_coord = round(coords[, 2], 4),
    stringsAsFactors = FALSE
  )

  out <- list(
    data = df,
    W = W,
    W_std = W_std,
    coords = coords,
    u_spatial = v,
    u_iid = u_iid,
    u_total = u_total,
    phi = phi,
    rho = rho
  )

  class(out) <- c("fastsae_sim_data", "list")
  return(out)
}

#' @export
print.fastsae_sim_data <- function(x, ...) {
  D <- nrow(x$data)
  n_na <- sum(is.na(x$data$y_gaussian))
  n_sampled <- D - n_na

  cat("======================================================\n")
  cat("  fastsae Simulated Multi-Distribution Area Data\n")
  cat("======================================================\n")
  cat(sprintf("Domains (D): %d (Sampled: %d, Unsampled: %d)\n", D, n_sampled, n_na))
  cat(sprintf("Spatial parameters: rho = %.2f, phi = %.2f\n", x$rho, x$phi))
  cat("Responses generated:\n")
  cat("  - Gaussian (FH):     y_gaussian (vardir)\n")
  cat("  - Poisson:           y_poisson (exposure)\n")
  cat("  - Binomial:          y_binomial (trials)\n")
  cat("  - Beta:              y_beta\n")
  cat("  - Negative Binomial: y_nbinomial\n")
  cat("  - Gamma:             y_gamma\n")
  cat(sprintf("Spatial matrices:      %d x %d (W binary, W_std row-standardized)\n", D, D))
  cat("======================================================\n")
  invisible(x)
}
