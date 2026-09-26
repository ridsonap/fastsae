#' Empirical Best Prediction for Area-Level Small Area Estimation
#'
#' @description
#' Estimates small area parameters using area-level models under various distributions
#' (Gaussian, Binomial, Poisson, Negative Binomial, Beta, Gamma) and spatial random
#' effect structures (non-spatial, BYM2, BYM, Besag/ICAR, Leroux) using
#' Integrated Nested Laplace Approximations (INLA) or frequentist Laplace Approximation.
#'
#' @param formula An object of class \code{formula} specifying the fixed-effects model
#'   (e.g., \code{y ~ x1 + x2}).
#' @param data A \code{data.frame} containing the area-level data (one row per area).
#' @param domain Vector, column name, or one-sided formula referencing the area/domain
#'   identifier in \code{data}. If \code{NULL}, domains are numbered consecutively.
#' @param family Character string specifying the response distribution / likelihood.
#'   Options:
#'   \itemize{
#'     \item \code{"gaussian"}: Continuous response (Fay-Herriot model). If \code{vardir}
#'       is provided, sampling variances are treated as known. If \code{vardir = NULL},
#'       residual variance is estimated.
#'     \item \code{"binomial"}: Binary / proportion response. Requires \code{trials}.
#'     \item \code{"poisson"}: Count response / disease rate. Supports \code{exposure}.
#'     \item \code{"nbinomial"}: Negative Binomial for overdispersed counts.
#'     \item \code{"beta"}: Continuous proportions strictly in (0, 1). If \code{vardir}
#'       is provided, area-specific precision is calculated via Janicki (2020) formula
#'       \eqn{\phi_i = y_i(1 - y_i)/V_i - 1} and injected via INLA's \code{scale} argument.
#'       If \code{trials} is provided, precision is scaled as \eqn{n_i - 1}.
#'     \item \code{"gamma"}: Skewed positive continuous response.
#'   }
#' @param spatial Character string specifying the spatial random effect structure:
#'   \itemize{
#'     \item \code{"none"}: Non-spatial IID area random intercept (default).
#'     \item \code{"bym2"}: Scaled Besag-York-Mollié 2 model (Riebler et al., 2016),
#'       decomposing area variance into spatial (\eqn{\phi}) and unstructured components.
#'     \item \code{"bym"}: Classic Besag-York-Mollié (1991) model.
#'     \item \code{"besag"}: Intrinsic Conditional Autoregressive (ICAR) model.
#'     \item \code{"generic1"}: Leroux spatial autoregressive model.
#'   }
#' @param W Proximity or spatial adjacency matrix. Can be a square \code{matrix},
#'   \code{Matrix}, or \code{spdep} \code{nb} or \code{listw} object. Dimensions must
#'   match the total number of domains in \code{data}. Required when \code{spatial != "none"}.
#' @param vardir Vector, column name, or formula specifying the sampling variances of the
#'   direct estimator (for \code{family = "gaussian"} or \code{family = "beta"}).
#' @param trials Vector, column name, or formula specifying sample sizes / total trials
#'   per area (for \code{family = "binomial"}).
#' @param exposure Vector, column name, or formula specifying expected counts or population
#'   exposure offsets (for \code{family = "poisson"} or \code{"nbinomial"}).
#' @param method Estimation method: \code{"inla"} (Bayesian INLA, default) or
#'   \code{"laplace"} (Frequentist GLMM with Laplace approximation via \code{lme4}).
#' @param strategy INLA approximation strategy: \code{"simplified.laplace"} (fast default),
#'   \code{"laplace"} (full Laplace approximation), or \code{"gaussian"}.
#' @param link Optional character string for link function. If \code{NULL}, default canonical
#'   link is used (identity for Gaussian, logit for Binomial/Beta, log for Poisson/NegBinom/Gamma).
#' @param scale_model Logical. If \code{TRUE} (default), scales the spatial graph so the
#'   marginal variance of the structured effect is approximately 1 (recommended for BYM2/Besag).
#' @param prior_prec List specifying the prior for random effect precision. Default is a
#'   Penalized Complexity (PC) prior: \code{list(prior = "pc.prec", param = c(1, 0.01))}.
#' @param prior_phi List specifying the PC-prior for the spatial mixing parameter \eqn{\phi}
#'   in the BYM2 model. Default is \code{list(prior = "pc", param = c(0.5, 0.5))}.
#' @param prior_rho Optional list specifying prior for spatial autocorrelation parameter (generic1/slm).
#' @param time Vector, column name, or one-sided formula referencing the time period
#'   identifier in \code{data} (e.g. \code{time = "year"} or \code{~year}). Required
#'   when \code{temporal != "none"}.
#' @param temporal Character string specifying the temporal random effect structure:
#'   \itemize{
#'     \item \code{"none"}: No temporal random effect (default cross-sectional model).
#'     \item \code{"rw1"}: First-order random walk across time periods.
#'     \item \code{"rw2"}: Second-order random walk across time periods.
#'     \item \code{"ar1"}: First-order autoregressive process across time periods.
#'     \item \code{"iid"}: Unstructured independent time effects.
#'   }
#' @param st_interaction Character string specifying the spatio-temporal structure:
#'   \itemize{
#'     \item \code{"none"}: Additive main spatial and temporal effects (default).
#'     \item \code{"domain-specific"}: Domain-specific temporal random walk / AR(1) dynamics with
#'       shared variance parameter (directly corresponds to the \pkg{tipsae} spatio-temporal model).
#'     \item \code{"separable"}: Grouped dynamic spatial field evolving over time via AR(1) / RW
#'       (corresponds to classical Spatio-Temporal Fay-Herriot models, Marhuenda et al. 2013).
#'     \item \code{"type1"}: Knorr-Held (2000) Type I interaction (unstructured space \eqn{\times} unstructured time).
#'     \item \code{"type2"}: Knorr-Held Type II interaction (unstructured space \eqn{\times} structured time).
#'     \item \code{"type3"}: Knorr-Held Type III interaction (structured space \eqn{\times} unstructured time).
#'     \item \code{"type4"}: Knorr-Held Type IV interaction (structured space \eqn{\times} structured time).
#'   }
#' @param prior_prec_time Optional list specifying the prior for temporal precision. Default is
#'   a PC prior: \code{list(prior = "pc.prec", param = c(1, 0.01))}.
#' @param prior_rho_time Optional list specifying prior for temporal autocorrelation parameter in AR(1).
#' @param print_result Logical. If \code{TRUE} (default), prints a summary of results.
#' @param ... Additional arguments passed to \code{INLA::inla()}.
#'
#' @returns An object of class \code{c("fastsae_ebp_area", "fastsae")} containing:
#' \itemize{
#'   \item \code{df_ebp}: Data frame with domain estimates, including \code{domain},
#'     \code{time} (if specified), observed \code{y}, predicted \code{ebp}, linear predictor \code{linear_pred},
#'     posterior standard error \code{sd}, \code{mse}, relative error \code{rse} (\%),
#'     95\% credible interval (\code{ci_lower}, \code{ci_upper}), and \code{random_effect}.
#'   \item \code{estcoef}: Data frame of estimated regression coefficients.
#'   \item \code{hyperpar}: Data frame of hyperparameter posterior estimates.
#'   \item \code{random_effect_var}: Estimated area random effect variance.
#'   \item \code{random_effect_var_time}: Estimated temporal random effect variance (if temporal).
#'   \item \code{phi}: Estimated spatial variance proportion (for BYM2) or spatial autocorrelation.
#'   \item \code{rho_time}: Estimated temporal autocorrelation parameter (for AR1).
#'   \item \code{goodness}: Model fit metrics (DIC, WAIC, Marginal Log-Likelihood).
#'   \item \code{family}: Response family used.
#'   \item \code{spatial}: Spatial model type used.
#'   \item \code{temporal}: Temporal model type used.
#'   \item \code{st_interaction}: Spatio-temporal interaction type used.
#'   \item \code{fit}: Raw fitted model object.
#'   \item \code{call}: Matched function call.
#' }
#'
#' @references
#' \enumerate{
#'   \item Rao, J. N. K., and Molina, I. (2015). \emph{Small Area Estimation}. John Wiley & Sons.
#'   \item Riebler, A., Sørbye, S. H., Simpson, D., and Rue, H. (2016). An intuitive Bayesian spatial model
#'     for disease mapping that accounts for scaling. \emph{Statistical Methods in Medical Research}, 25(4), 1145-1165.
#'   \item Marhuenda, Y., Molina, I., and Morales, D. (2013). Small area estimation with spatio-temporal Fay-Herriot models.
#'     \emph{Computational Statistics & Data Analysis}, 58, 308-325.
#'   \item De Nicolò, S., and Gardini, A. (2024). The R Package tipsae: Tools for Mapping Proportions and Indicators on the Unit Interval.
#'     \emph{Journal of Statistical Software}, 108(1), 1-36.
#'   \item Rue, H., Martino, S., and Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian
#'     models by using integrated nested Laplace approximations. \emph{Journal of the Royal Statistical Society: Series B}, 71(2), 319-392.
#' }
#'
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   library(fastsae)
#'
#'   # 1. Non-spatial Gaussian Fay-Herriot with INLA
#'   m_norm <- ebp_area(
#'     y ~ x1 + x2 + x3,
#'     data = mys,
#'     vardir = "vardir",
#'     family = "gaussian"
#'   )
#'
#'   # 2. Spatial BYM2 Gaussian Fay-Herriot with INLA
#'   m_bym2 <- ebp_area(
#'     y ~ x1 + x2 + x3,
#'     data = mys,
#'     vardir = "vardir",
#'     family = "gaussian",
#'     spatial = "bym2",
#'     W = mys_proxmat
#'   )
#' }
#' }
ebp_area <- function(
  formula,
  data,
  domain = NULL,
  time = NULL,
  family = c("gaussian", "binomial", "poisson", "nbinomial", "beta", "gamma"),
  spatial = c("none", "bym2", "bym", "besag", "generic1", "slm"),
  temporal = c("none", "rw1", "rw2", "ar1", "iid"),
  st_interaction = c("none", "domain-specific", "separable", "type1", "type2", "type3", "type4"),
  W = NULL,
  vardir = NULL,
  trials = NULL,
  exposure = NULL,
  method = c("inla", "laplace"),
  strategy = c("simplified.laplace", "laplace", "gaussian"),
  link = NULL,
  scale_model = TRUE,
  prior_prec = list(prior = "pc.prec", param = c(1, 0.01)),
  prior_phi = list(prior = "pc", param = c(0.5, 0.5)),
  prior_rho = NULL,
  prior_prec_time = NULL,
  prior_rho_time = NULL,
  print_result = TRUE,
  ...
) {
  call_matched <- match.call()
  family <- match.arg(tolower(family), choices = c("gaussian", "binomial", "poisson", "nbinomial", "beta", "gamma"))
  spatial <- match.arg(tolower(spatial), choices = c("none", "bym2", "bym", "besag", "generic1", "slm"))
  temporal <- match.arg(tolower(temporal), choices = c("none", "rw1", "rw2", "ar1", "iid"))
  st_interaction <- match.arg(tolower(st_interaction), choices = c("none", "domain-specific", "separable", "type1", "type2", "type3", "type4"))
  method <- match.arg(tolower(method), choices = c("inla", "laplace"))
  strategy <- match.arg(tolower(strategy), choices = c("simplified.laplace", "laplace", "gaussian"))

  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame or tibble.")
  }

  n_obs <- nrow(data)

  # 1. Extract domain & time identifiers
  if (is.null(domain)) {
    if (is.null(time)) {
      domain_vec <- seq_len(n_obs)
    } else {
      cli::cli_abort("When {.arg time} is specified, {.arg domain} must also be specified to identify domains across time periods.")
    }
  } else {
    domain_vec <- .get_variable(data, domain)
  }

  unique_domains <- unique(domain_vec)
  n_unique_domains <- length(unique_domains)

  time_vec <- if (!is.null(time)) .get_variable(data, time) else NULL
  if (temporal != "none" && is.null(time_vec)) {
    cli::cli_abort("When {.code temporal != 'none'}, {.arg time} (time/period column) must be specified.")
  }

  unique_times <- if (!is.null(time_vec)) sort(unique(time_vec)) else NULL
  n_unique_times <- if (!is.null(unique_times)) length(unique_times) else 0

  if (temporal != "none" && n_unique_times < 2) {
    cli::cli_abort("Temporal modeling requires at least 2 distinct time periods in {.arg time}.")
  }

  # 2. Model frame & Response validation
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  y <- stats::model.response(mf)

  # Validate optional inputs
  vardir_vec <- if (!is.null(vardir)) .get_variable(data, vardir) else NULL
  trials_vec <- if (!is.null(trials)) .get_variable(data, trials) else NULL
  exposure_vec <- if (!is.null(exposure)) .get_variable(data, exposure) else NULL

  response_var <- all.vars(formula)[1]

  if (family == "binomial") {
    if (is.null(trials_vec)) {
      cli::cli_abort("For {.code family = 'binomial'}, {.arg trials} (sample sizes / total trials per area) must be specified.")
    }
    # Check if response contains proportions (values in [0, 1] with non-integers)
    y_non_na <- y[!is.na(y)]
    if (length(y_non_na) > 0 && all(y_non_na >= 0 & y_non_na <= 1) && any(y_non_na %% 1 != 0)) {
      cli::cli_alert_info("Response {.arg y} appears to be proportions; converting to integer counts: {.code round(y * trials)}.")
      y_counts <- as.integer(round(y * trials_vec))
      data[[response_var]] <- y_counts
      y <- y_counts
    } else {
      y_counts <- ifelse(is.na(y), NA_integer_, as.integer(round(y)))
      data[[response_var]] <- y_counts
      y <- y_counts
    }
  }

  if (family == "gamma") {
    y_non_na <- y[!is.na(y)]
    if (any(y_non_na <= 0, na.rm = TRUE)) {
      cli::cli_abort("For {.code family = 'gamma'}, response {.arg y} must be strictly positive (> 0).")
    }
  }

  # 3. Spatial weights handling
  W_obj <- NULL
  if (spatial != "none") {
    W_obj <- .convert_spatial_weights(W, n_domains = n_unique_domains, spatial = spatial, domain_names = unique_domains)
  }

  # --------------------------------------------------------------------------
  # Route estimation: Frequentist Laplace (lme4) vs Bayesian INLA
  # --------------------------------------------------------------------------
  if (method == "laplace" && spatial == "none" && temporal == "none" && family %in% c("binomial", "poisson")) {
    out <- .fit_glmm_laplace(
      formula = formula,
      data = data,
      domain = domain_vec,
      family = family,
      trials = trials_vec,
      exposure = exposure_vec,
      call = call_matched
    )
  } else {
    if (method == "laplace") {
      if (family == "gaussian") {
        cli::cli_alert_info("For Gaussian area-level models with known sampling variance, utilizing INLA's full Laplace approximation strategy.")
      } else if (spatial != "none" || temporal != "none") {
        cli::cli_alert_info("Spatial and temporal models require INLA; utilizing INLA's full Laplace approximation strategy.")
      }
      strategy <- "laplace"
    }

    .check_inla_installed()

    out <- .fit_inla_area(
      formula = formula,
      data = data,
      y = y,
      domain = domain_vec,
      time = time_vec,
      family = family,
      spatial = spatial,
      temporal = temporal,
      st_interaction = st_interaction,
      W_obj = W_obj,
      vardir = vardir_vec,
      trials = trials_vec,
      exposure = exposure_vec,
      strategy = strategy,
      link = link,
      scale_model = scale_model,
      prior_prec = prior_prec,
      prior_phi = prior_phi,
      prior_rho = prior_rho,
      prior_prec_time = prior_prec_time,
      prior_rho_time = prior_rho_time,
      call = call_matched,
      ...
    )
  }

  out$W <- W

  if (print_result) {
    print(out)
  }

  return(out)
}

#' Internal worker to fit area-level model via INLA
#' @noRd
.fit_inla_area <- function(
  formula,
  data,
  y,
  domain,
  time = NULL,
  family,
  spatial = "none",
  temporal = "none",
  st_interaction = "none",
  W_obj = NULL,
  vardir = NULL,
  trials = NULL,
  exposure = NULL,
  strategy = "simplified.laplace",
  link = NULL,
  scale_model = TRUE,
  prior_prec = list(prior = "pc.prec", param = c(1, 0.01)),
  prior_phi = list(prior = "pc", param = c(0.5, 0.5)),
  prior_rho = NULL,
  prior_prec_time = NULL,
  prior_rho_time = NULL,
  call = NULL,
  ...
) {
  n_obs <- nrow(data)

  unique_domains <- unique(domain)
  n_domains <- length(unique_domains)
  domain_id <- as.integer(factor(domain, levels = unique_domains))

  inla_data <- data
  inla_data$..domain_id.. <- domain_id
  inla_data$..obs_id.. <- seq_len(n_obs)

  unique_times <- NULL
  n_times <- 1
  time_id <- NULL
  if (!is.null(time)) {
    unique_times <- sort(unique(time))
    n_times <- length(unique_times)
    time_id <- as.integer(factor(time, levels = unique_times))
    inla_data$..time_id.. <- time_id
  }

  # Map family name to INLA family
  inla_family <- switch(family,
    "gaussian"  = "gaussian",
    "binomial"  = "binomial",
    "poisson"   = "poisson",
    "nbinomial" = "nbinomial",
    "beta"      = "beta",
    "gamma"     = "gamma",
    "gaussian"
  )

  # Construct INLA formula with random effect
  response_var <- all.vars(formula)[1]
  terms_fixed <- attr(stats::terms(formula), "term.labels")
  fixed_str <- if (length(terms_fixed) > 0) paste(terms_fixed, collapse = " + ") else "1"

  X_mat <- NULL
  rho_range <- NULL

  # Setup temporal prior hyperparameters
  if (is.null(prior_prec_time)) {
    prior_prec_time <- list(prior = "pc.prec", param = c(1, 0.01))
  }
  hyper_time <- list(prec = prior_prec_time)
  if (temporal == "ar1" && !is.null(prior_rho_time)) {
    hyper_time$rho <- prior_rho_time
  }

  rand_terms <- c()

  # 1. Spatial random effect specification
  if (spatial == "none") {
    if (temporal == "none" || st_interaction %in% c("none", "domain-specific", "type1", "type2")) {
      rand_terms <- c(rand_terms, "f(..domain_id.., model = 'iid', hyper = list(prec = prior_prec))")
    }
  } else if (spatial == "bym2") {
    if (st_interaction != "separable") {
      rand_terms <- c(rand_terms, paste0(
        "f(..domain_id.., model = 'bym2', graph = W_obj$graph, scale.model = ",
        scale_model,
        ", hyper = list(prec = prior_prec, phi = prior_phi))"
      ))
    }
  } else if (spatial == "bym") {
    if (st_interaction != "separable") {
      rand_terms <- c(rand_terms, paste0(
        "f(..domain_id.., model = 'bym', graph = W_obj$graph, scale.model = ",
        scale_model,
        ", hyper = list(prec.unstruct = prior_prec, prec.spatial = prior_prec))"
      ))
    }
  } else if (spatial == "besag") {
    if (st_interaction != "separable") {
      rand_terms <- c(rand_terms, paste0(
        "f(..domain_id.., model = 'besag', graph = W_obj$graph, scale.model = ",
        scale_model,
        ", hyper = list(prec = prior_prec))"
      ))
    }
  } else if (spatial == "generic1") {
    adj <- W_obj$graph
    R_mat <- diag(rowSums(adj)) - adj
    C_mat <- diag(nrow(adj)) - R_mat
    hyper_generic1 <- list(prec = prior_prec)
    if (!is.null(prior_rho)) {
      hyper_generic1$beta <- prior_rho
    } else if (!is.null(prior_phi) && !identical(prior_phi$prior, "pc")) {
      hyper_generic1$beta <- prior_phi
    }
    if (st_interaction != "separable") {
      rand_terms <- c(rand_terms, "f(..domain_id.., model = 'generic1', Cmatrix = C_mat, hyper = hyper_generic1)")
    }
  } else if (spatial == "slm") {
    W_mat <- W_obj$W_mat
    rs <- rowSums(W_mat)
    rs_inv <- ifelse(rs > 0, 1 / rs, 0)
    W_std <- W_mat * rs_inv
    W_sparse <- Matrix::Matrix(W_std, sparse = TRUE)
    e <- eigen(W_std, only.values = TRUE)$values
    re_e <- Re(e[abs(Im(e)) < 1e-5])
    rho_min <- 1 / min(re_e)
    rho_max <- 1 / max(re_e)
    if (is.infinite(rho_min) || rho_min < -1) rho_min <- -0.999
    if (is.infinite(rho_max) || rho_max > 1) rho_max <- 0.999
    rho_range <- c(rho_min, rho_max)

    terms_no_y <- stats::delete.response(stats::terms(formula))
    mf_X <- stats::model.frame(terms_no_y, data = data, na.action = stats::na.pass)
    X_mat <- stats::model.matrix(terms_no_y, data = mf_X)
    p_cov <- ncol(X_mat)
    Q_beta <- Matrix::Diagonal(p_cov, 1e-4)

    args_slm <- list(
      rho.min = rho_min,
      rho.max = rho_max,
      W = W_sparse,
      X = X_mat,
      Q.beta = Q_beta
    )
    hyper_slm <- list(prec = prior_prec)
    if (!is.null(prior_rho)) {
      hyper_slm$rho <- prior_rho
    } else if (!is.null(prior_phi) && !identical(prior_phi$prior, "pc")) {
      hyper_slm$rho <- prior_phi
    } else {
      hyper_slm$rho <- list(prior = "logitbeta", param = c(1, 1))
    }

    rand_terms <- c(rand_terms, "f(..domain_id.., model = 'slm', args.slm = args_slm, hyper = hyper_slm)")
  }

  # 2. Temporal & Space-Time random effect specification
  if (temporal != "none") {
    scale_model_time <- if (temporal %in% c("rw1", "rw2")) paste0(", scale.model = ", scale_model) else ""

    if (st_interaction == "none") {
      # Additive temporal main effect
      rand_terms <- c(rand_terms, paste0(
        "f(..time_id.., model = '", temporal, "'", scale_model_time, ", hyper = hyper_time)"
      ))
    } else if (st_interaction == "domain-specific") {
      # Domain-specific temporal dynamics (directly corresponds to tipsae)
      rand_terms <- c(rand_terms, paste0(
        "f(..time_id.., model = '", temporal, "', group = ..domain_id.., control.group = list(model = 'iid')",
        scale_model_time, ", hyper = hyper_time)"
      ))
    } else if (st_interaction == "separable") {
      # Dynamic spatial field evolving over time
      if (spatial == "bym2") {
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id.., model = 'bym2', graph = W_obj$graph, group = ..time_id.., control.group = list(model = '",
          temporal, "'), scale.model = ", scale_model, ", hyper = list(prec = prior_prec, phi = prior_phi))"
        ))
      } else if (spatial == "besag") {
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id.., model = 'besag', graph = W_obj$graph, group = ..time_id.., control.group = list(model = '",
          temporal, "'), scale.model = ", scale_model, ", hyper = list(prec = prior_prec))"
        ))
      } else if (spatial == "bym") {
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id.., model = 'bym', graph = W_obj$graph, group = ..time_id.., control.group = list(model = '",
          temporal, "'), scale.model = ", scale_model, ", hyper = list(prec.unstruct = prior_prec, prec.spatial = prior_prec))"
        ))
      } else if (spatial == "generic1") {
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id.., model = 'generic1', Cmatrix = C_mat, group = ..time_id.., control.group = list(model = '",
          temporal, "'), hyper = hyper_generic1)"
        ))
      } else {
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id.., model = 'iid', group = ..time_id.., control.group = list(model = '",
          temporal, "'), hyper = list(prec = prior_prec))"
        ))
      }
    } else if (st_interaction == "type1") {
      rand_terms <- c(rand_terms, paste0(
        "f(..time_id.., model = '", temporal, "'", scale_model_time, ", hyper = hyper_time)"
      ))
      rand_terms <- c(rand_terms, "f(..obs_id.., model = 'iid', hyper = list(prec = prior_prec))")
    } else if (st_interaction == "type2") {
      rand_terms <- c(rand_terms, paste0(
        "f(..time_id.., model = '", temporal, "'", scale_model_time, ", hyper = hyper_time)"
      ))
      rand_terms <- c(rand_terms, paste0(
        "f(..time_id.., model = '", temporal, "', group = ..domain_id.., control.group = list(model = 'iid')",
        scale_model_time, ", hyper = hyper_time)"
      ))
    } else if (st_interaction == "type3") {
      rand_terms <- c(rand_terms, paste0(
        "f(..time_id.., model = '", temporal, "'", scale_model_time, ", hyper = hyper_time)"
      ))
      if (spatial != "none") {
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id.., model = '", spatial, "', graph = W_obj$graph, group = ..time_id.., control.group = list(model = 'iid'), scale.model = ",
          scale_model, ", hyper = list(prec = prior_prec))"
        ))
      }
    } else if (st_interaction == "type4") {
      rand_terms <- c(rand_terms, paste0(
        "f(..time_id.., model = '", temporal, "'", scale_model_time, ", hyper = hyper_time)"
      ))
      if (spatial != "none") {
        rand_terms <- c(rand_terms, paste0(
          "f(..domain_id.., model = '", spatial, "', graph = W_obj$graph, group = ..time_id.., control.group = list(model = '",
          temporal, "'), scale.model = ", scale_model, ", hyper = list(prec = prior_prec))"
        ))
      }
    }
  }

  if (spatial == "slm") {
    inla_formula_str <- paste(response_var, "~ -1 +", paste(rand_terms, collapse = " + "))
  } else {
    inla_formula_str <- paste(response_var, "~", fixed_str, "+", paste(rand_terms, collapse = " + "))
  }

  inla_formula <- stats::as.formula(inla_formula_str)

  # INLA controls
  control_compute <- list(
    dic = TRUE,
    waic = TRUE,
    cpo = TRUE,
    mlik = TRUE,
    config = TRUE
  )

  # Predictor control: compute fitted values on original scale
  control_predictor <- list(compute = TRUE, link = 1)

  # INLA strategy control
  control_inla <- list(strategy = strategy)

  # Family specific settings
  inla_args <- list(
    formula = inla_formula,
    family = inla_family,
    data = inla_data,
    control.compute = control_compute,
    control.predictor = control_predictor,
    control.inla = control_inla,
    ...
  )

  # Gaussian with known sampling variances (Fay-Herriot)
  if (family == "gaussian" && !is.null(vardir)) {
    if (any(vardir[!is.na(y)] <= 0, na.rm = TRUE)) {
      cli::cli_abort("{.arg vardir} must be strictly positive for sampled domains.")
    }
    inla_args$scale <- 1 / vardir
    inla_args$control.family <- list(
      hyper = list(prec = list(initial = 0, fixed = TRUE))
    )
  }

  # Beta with area-specific sampling variances or sample sizes (Janicki 2020)
  if (family == "beta") {
    if (!is.null(vardir)) {
      if (any(vardir[!is.na(y)] <= 0, na.rm = TRUE)) {
        cli::cli_abort("{.arg vardir} must be strictly positive for sampled domains.")
      }
      phi_dir <- (y * (1 - y) / vardir) - 1
      phi_dir[phi_dir < 1] <- 1
      inla_args$scale <- phi_dir
      inla_args$control.family <- list(
        hyper = list(theta = list(initial = 0, fixed = TRUE))
      )
    } else if (!is.null(trials)) {
      phi_trials <- trials - 1
      phi_trials[phi_trials < 1] <- 1
      inla_args$scale <- phi_trials
      inla_args$control.family <- list(
        hyper = list(theta = list(initial = 0, fixed = TRUE))
      )
    }
  }

  # Gamma with area-specific sampling variances or CV (Gamma SAE)
  if (family == "gamma" && !is.null(vardir)) {
    if (any(vardir[!is.na(y)] <= 0, na.rm = TRUE)) {
      cli::cli_abort("{.arg vardir} must be strictly positive for sampled domains.")
    }
    s_gamma <- rep(1, n_obs)
    idx_valid <- which(!is.na(y) & !is.na(vardir) & vardir > 0)
    s_gamma[idx_valid] <- (y[idx_valid]^2) / vardir[idx_valid]
    s_gamma[s_gamma <= 0 | is.na(s_gamma)] <- 1
    inla_args$scale <- s_gamma
    inla_args$control.family <- list(
      hyper = list(theta = list(initial = 0, fixed = TRUE))
    )
  }

  # Binomial with trials
  if (family == "binomial" && !is.null(trials)) {
    inla_args$Ntrials <- trials
  }

  # Poisson / Negative Binomial with exposure/offset
  if (family %in% c("poisson", "nbinomial") && !is.null(exposure)) {
    inla_args$E <- exposure
  }

  # Run INLA fitting
  fit <- tryCatch(
    do.call(INLA::inla, inla_args),
    error = function(e) {
      cli::cli_abort(c(
        "Fitting model with INLA failed.",
        "x" = e$message
      ))
    }
  )

  # Extract formatted results
  res <- .extract_inla_results(
    fit = fit,
    data = data,
    y = y,
    domain = domain,
    time = time,
    family = family,
    spatial = spatial,
    temporal = temporal,
    st_interaction = st_interaction,
    vardir = vardir,
    trials = trials,
    exposure = exposure,
    X_mat = X_mat,
    rho_range = rho_range,
    domain_id = domain_id,
    time_id = time_id,
    unique_domains = unique_domains,
    unique_times = unique_times,
    call = call
  )

  return(res)
}

#' Internal worker to fit non-spatial model via frequentist GLMM (Laplace approximation)
#' @noRd
.fit_glmm_laplace <- function(
  formula,
  data,
  domain,
  family,
  trials = NULL,
  exposure = NULL,
  call = NULL
) {
  n_domains <- nrow(data)
  lmer_data <- data
  lmer_data$..domain_id.. <- factor(domain)

  response_var <- all.vars(formula)[1]
  terms_fixed <- attr(stats::terms(formula), "term.labels")
  fixed_str <- if (length(terms_fixed) > 0) paste(terms_fixed, collapse = " + ") else "1"
  formula_lmer <- stats::as.formula(paste(response_var, "~", fixed_str, "+ (1 | ..domain_id..)"))

  y <- stats::model.response(stats::model.frame(formula, data, na.action = stats::na.pass))

  fit <- NULL
  if (family == "gaussian") {
    fit <- lme4::lmer(formula_lmer, data = lmer_data, REML = FALSE)
  } else if (family == "binomial") {
    if (!is.null(trials)) {
      lmer_data$..successes.. <- y
      lmer_data$..failures.. <- trials - y
      formula_bin <- stats::as.formula(paste("cbind(..successes.., ..failures..) ~", fixed_str, "+ (1 | ..domain_id..)"))
      fit <- lme4::glmer(formula_bin, data = lmer_data, family = stats::binomial(link = "logit"), nAGQ = 1)
    } else {
      fit <- lme4::glmer(formula_lmer, data = lmer_data, family = stats::binomial(link = "logit"), nAGQ = 1)
    }
  } else if (family == "poisson") {
    if (!is.null(exposure)) {
      lmer_data$..offset.. <- log(exposure)
      formula_pois <- stats::as.formula(paste(response_var, "~", fixed_str, "+ offset(..offset..) + (1 | ..domain_id..)"))
      fit <- lme4::glmer(formula_pois, data = lmer_data, family = stats::poisson(link = "log"), nAGQ = 1)
    } else {
      fit <- lme4::glmer(formula_lmer, data = lmer_data, family = stats::poisson(link = "log"), nAGQ = 1)
    }
  }

  # Extract coefficients
  sum_fit <- summary(fit)
  coef_mat <- as.data.frame(sum_fit$coefficients)
  estcoef <- data.frame(
    beta = coef_mat[, 1],
    std.error = coef_mat[, 2],
    zvalue = coef_mat[, 3],
    pvalue = 2 * stats::pnorm(abs(coef_mat[, 3]), lower.tail = FALSE),
    ci_lower = coef_mat[, 1] - 1.96 * coef_mat[, 2],
    ci_upper = coef_mat[, 1] + 1.96 * coef_mat[, 2],
    row.names = rownames(coef_mat)
  )

  # Random effect variance
  vc <- lme4::VarCorr(fit)
  sigma2_u <- as.numeric(vc[[1]])

  # Predictions
  preds <- stats::predict(fit, newdata = lmer_data, type = "response", allow.new.levels = TRUE)
  linpred <- stats::predict(fit, newdata = lmer_data, type = "link", allow.new.levels = TRUE)
  ranefs <- lme4::ranef(fit)$..domain_id..[, 1]

  df_ebp <- data.frame(
    domain = domain,
    y = y,
    ebp = preds,
    linear_pred = linpred,
    sd = NA_real_,
    mse = NA_real_,
    rse = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    random_effect = if (length(ranefs) == n_domains) ranefs else rep(NA_real_, n_domains),
    stringsAsFactors = FALSE
  )

  goodness <- c(
    AIC = stats::AIC(fit),
    BIC = stats::BIC(fit),
    LogLik = as.numeric(stats::logLik(fit))
  )

  out <- list(
    df_ebp = df_ebp,
    ebp = df_ebp,
    df_eblup = df_ebp,
    estcoef = estcoef,
    hyperpar = data.frame(Parameter = "sigma2_u", Estimate = sigma2_u),
    random_effect_var = sigma2_u,
    phi = NULL,
    goodness = goodness,
    family = family,
    spatial = "none",
    level = "area",
    model = paste0("EBP-", toupper(family), " (Laplace GLMM)"),
    convergence = TRUE,
    fit = fit,
    call = call
  )

  class(out) <- c("fastsae_ebp_area", "fastsae")
  return(out)
}
