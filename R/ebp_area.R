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
#'     \item \code{"beta"}: Continuous proportions strictly in (0, 1).
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
#'   direct estimator (for \code{family = "gaussian"}).
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
#' @param print_result Logical. If \code{TRUE} (default), prints a summary of results.
#' @param ... Additional arguments passed to \code{INLA::inla()}.
#'
#' @returns An object of class \code{c("fastsae_ebp_area", "fastsae")} containing:
#' \itemize{
#'   \item \code{df_ebp}: Data frame with domain estimates, including \code{domain},
#'     observed \code{y}, predicted \code{ebp}, linear predictor \code{linear_pred},
#'     posterior standard error \code{sd}, \code{mse}, relative error \code{rse} (\%),
#'     95\% credible interval (\code{ci_lower}, \code{ci_upper}), and \code{random_effect}.
#'   \item \code{estcoef}: Data frame of estimated regression coefficients.
#'   \item \code{hyperpar}: Data frame of hyperparameter posterior estimates.
#'   \item \code{random_effect_var}: Estimated area random effect variance.
#'   \item \code{phi}: Estimated spatial variance proportion (for BYM2).
#'   \item \code{goodness}: Model fit metrics (DIC, WAIC, Marginal Log-Likelihood).
#'   \item \code{family}: Response family used.
#'   \item \code{spatial}: Spatial model type used.
#'   \item \code{fit}: Raw fitted model object.
#'   \item \code{call}: Matched function call.
#' }
#'
#' @references
#' \enumerate{
#'   \item Rao, J. N. K., and Molina, I. (2015). \emph{Small Area Estimation}. John Wiley & Sons.
#'   \item Riebler, A., Sørbye, S. H., Simpson, D., and Rue, H. (2016). An intuitive Bayesian spatial model
#'     for disease mapping that accounts for scaling. \emph{Statistical Methods in Medical Research}, 25(4), 1145-1165.
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
  family = c("gaussian", "binomial", "poisson", "nbinomial", "beta", "gamma"),
  spatial = c("none", "bym2", "bym", "besag", "generic1", "slm"),
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
  print_result = TRUE,
  ...
) {
  call_matched <- match.call()
  family <- match.arg(tolower(family), choices = c("gaussian", "binomial", "poisson", "nbinomial", "beta", "gamma"))
  spatial <- match.arg(tolower(spatial), choices = c("none", "bym2", "bym", "besag", "generic1", "slm"))
  method <- match.arg(tolower(method), choices = c("inla", "laplace"))
  strategy <- match.arg(tolower(strategy), choices = c("simplified.laplace", "laplace", "gaussian"))

  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame or tibble.")
  }

  n_domains <- nrow(data)

  # 1. Extract domain identifiers
  if (is.null(domain)) {
    domain_vec <- seq_len(n_domains)
  } else {
    domain_vec <- .get_variable(data, domain)
  }

  # 2. Model frame & Response validation
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  y <- stats::model.response(mf)

  # Validate optional inputs
  vardir_vec <- if (!is.null(vardir)) .get_variable(data, vardir) else NULL
  trials_vec <- if (!is.null(trials)) .get_variable(data, trials) else NULL
  exposure_vec <- if (!is.null(exposure)) .get_variable(data, exposure) else NULL

  if (family == "binomial" && is.null(trials_vec)) {
    cli::cli_abort("For {.code family = 'binomial'}, {.arg trials} (sample sizes / total trials per area) must be specified.")
  }

  # 3. Spatial weights handling
  W_obj <- NULL
  if (spatial != "none") {
    W_obj <- .convert_spatial_weights(W, n_domains = n_domains, spatial = spatial)
  }

  # --------------------------------------------------------------------------
  # Route estimation: Frequentist Laplace (lme4) vs Bayesian INLA
  # --------------------------------------------------------------------------
  if (method == "laplace" && spatial == "none" && family %in% c("binomial", "poisson")) {
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
      } else if (spatial != "none") {
        cli::cli_alert_info("Spatial models require INLA; utilizing INLA's full Laplace approximation strategy.")
      }
      strategy <- "laplace"
    }

    .check_inla_installed()

    out <- .fit_inla_area(
      formula = formula,
      data = data,
      y = y,
      domain = domain_vec,
      family = family,
      spatial = spatial,
      W_obj = W_obj,
      vardir = vardir_vec,
      trials = trials_vec,
      exposure = exposure_vec,
      strategy = strategy,
      link = link,
      scale_model = scale_model,
      prior_prec = prior_prec,
      prior_phi = prior_phi,
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
  family,
  spatial,
  W_obj,
  vardir,
  trials,
  exposure,
  strategy,
  link,
  scale_model,
  prior_prec,
  prior_phi,
  call,
  ...
) {
  n_domains <- nrow(data)

  # Create an internal domain index 1:n_domains for INLA
  inla_data <- data
  inla_data$..domain_id.. <- seq_len(n_domains)

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

  # Spatial random effect specification
  if (spatial == "none") {
    rand_term <- "f(..domain_id.., model = 'iid', hyper = list(prec = prior_prec))"
  } else if (spatial == "bym2") {
    rand_term <- paste0(
      "f(..domain_id.., model = 'bym2', graph = W_obj$graph, scale.model = ",
      scale_model,
      ", hyper = list(prec = prior_prec, phi = prior_phi))"
    )
  } else if (spatial == "bym") {
    rand_term <- paste0(
      "f(..domain_id.., model = 'bym', graph = W_obj$graph, scale.model = ",
      scale_model,
      ", hyper = list(prec.unstruct = prior_prec, prec.spatial = prior_prec))"
    )
  } else if (spatial == "besag") {
    rand_term <- paste0(
      "f(..domain_id.., model = 'besag', graph = W_obj$graph, scale.model = ",
      scale_model,
      ", hyper = list(prec = prior_prec))"
    )
  } else if (spatial == "generic1") {
    # Leroux model: Q = tau * (rho * (diag(degree) - adj) + (1 - rho) * I)
    adj <- W_obj$graph
    R_mat <- diag(rowSums(adj)) - adj
    rand_term <- "f(..domain_id.., model = 'generic1', Cmatrix = R_mat, hyper = list(prec = prior_prec))"
  } else {
    rand_term <- "f(..domain_id.., model = 'iid')"
  }

  inla_formula_str <- paste(response_var, "~", fixed_str, "+", rand_term)
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
    # Fix Gaussian observation precision at scale
    inla_args$scale <- 1 / vardir
    inla_args$control.family <- list(
      hyper = list(prec = list(initial = 0, fixed = TRUE))
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
    family = family,
    spatial = spatial,
    vardir = vardir,
    trials = trials,
    exposure = exposure,
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
