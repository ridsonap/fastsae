#' Two-Fold Hierarchical Bayes for Sub-Area Level Small Area Estimation
#'
#' @description
#' Estimates small area parameters under a two-fold sub-area level model
#' (Torabi & Rao, 2014 Bayesian hierarchical framework) using Integrated Nested
#' Laplace Approximations (\pkg{INLA}) or, with \code{engine = jax_engine()}, MCMC
#' (NumPyro's No-U-Turn sampler on JAX). The model accounts for nested random effects
#' at both the primary area level (\eqn{v_i}) and the nested sub-area level (\eqn{u_{ij}}),
#' supporting Gaussian, Binomial (logistic), and Poisson responses, with optional
#' spatial correlation structures (BYM2, Besag) across areas and aggregate area-level
#' estimation.
#'
#' @param formula An object of class \code{formula} specifying the fixed-effects model
#'   (e.g., \code{y ~ x1 + x2}).
#' @param vardir Character string or vector specifying direct sampling variances
#'   (\eqn{\psi_{ij}}) for continuous (Gaussian) responses. Required when \code{family = "gaussian"}.
#' @param domain Character string, column name, or vector referencing the primary area identifier
#'   in \code{data}. Must be supplied.
#' @param subarea Optional character string, column name, or vector referencing the sub-area identifier
#'   in \code{data}. If \code{NULL}, sub-areas are numbered consecutively per domain.
#' @param data A \code{data.frame} or tibble containing one row per sub-area.
#' @param family Character string specifying the response likelihood. Options:
#'   \itemize{
#'     \item \code{"gaussian"}: Continuous response with known sampling variance (default).
#'     \item \code{"binomial"}: Binary / count proportion with sample sizes specified by \code{trials}.
#'     \item \code{"poisson"}: Count response with expected exposures specified by \code{exposure}.
#'   }
#' @param spatial Character string specifying the spatial structure across primary areas:
#'   \itemize{
#'     \item \code{"none"}: Non-spatial independent and identically distributed (IID) area effects (default).
#'     \item \code{"bym2"}: Scaled Besag-York-Molli\enc{é}{e} 2 spatial model (Riebler et al., 2016).
#'     \item \code{"besag"}: Intrinsic Conditional Autoregressive (ICAR) spatial model.
#'   }
#' @param W Proximity or spatial adjacency matrix across unique primary areas. Can be a square \code{matrix},
#'   \code{Matrix}, or \code{spdep} \code{nb} or \code{listw} object. Dimensions must
#'   match the total number of unique primary areas in \code{domain}. Required when \code{spatial != "none"}.
#' @param weight Optional character string or numeric vector specifying sub-area population weights
#'   (\eqn{w_{ij}}) used to aggregate sub-area predictions to primary area estimates
#'   \eqn{\hat{\theta}_i = \sum_j w_{ij} \hat{\theta}_{ij}}. If \code{NULL}, equal weights per area are used.
#' @param trials Optional character string or vector specifying total number of trials per sub-area.
#'   Required when \code{family = "binomial"}.
#' @param exposure Optional character string or vector specifying expected baseline exposures per sub-area.
#'   Required when \code{family = "poisson"}.
#' @param engine Estimation engine: \code{"inla"} (default) or an \code{\link{inla_engine}()}
#'   object, or \code{"jax"} or a \code{\link{jax_engine}()} object. \code{inla_engine()} holds the approximation strategy, the priors (\code{prior_prec} for the area
#'   effect \eqn{\tau_v} and \code{prior_prec_subarea} for the sub-area effect \eqn{\tau_u},
#'   both by default \code{list(prior = "pc.prec", param = c(1, 0.01))}; \code{prior_phi} for
#'   BYM2), \code{scale_model}, the number of posterior \code{draws} for the area aggregates
#'   (default 200) and their \code{seed}, and further \code{INLA::inla()} arguments. The
#'   Gaussian model is fitted in INLA's \code{"classic"} mode, the binomial and Poisson models
#'   in INLA's default \code{"compact"} mode, as in \code{\link{hb_area}};
#'   \code{inla_engine(inla.mode = )} selects a mode. \code{jax_engine()} samples the same
#'   model by MCMC (see the Engines section).
#' @param compute_area Logical. If \code{TRUE} (default), computes aggregate primary area-level
#'   estimates and posterior standard errors from joint posterior draws (sampled from the INLA
#'   fit, or the MCMC draws with \code{jax_engine()}).
#' @param print_result Logical. If \code{TRUE} (default), prints a summary of results.
#' @param ... Deprecated. \code{strategy}, \code{scale_model}, \code{prior_prec_area} (now
#'   \code{prior_prec}), \code{prior_prec_subarea}, \code{prior_phi}, \code{n_samples} (now
#'   \code{draws}) and arguments of \code{INLA::inla()} are now given to
#'   \code{\link{inla_engine}()}; passed here they are forwarded to it with a warning.
#'
#' @returns An object of class \code{c("fastsae_hb_twofold", "fastsae_hb", "fastsae")} containing:
#' \itemize{
#'   \item \code{df_hb}: Data frame with sub-area level estimates:
#'     \itemize{
#'       \item \code{domain}: Primary area identifier.
#'       \item \code{subarea}: Sub-area identifier.
#'       \item \code{y}: Direct estimate / observed response (NA for unsampled).
#'       \item \code{hb}: Posterior mean estimate of sub-area mean.
#'       \item \code{linear_pred}: Posterior mean of linear predictor.
#'       \item \code{vardir}: Direct sampling variance (for Gaussian).
#'       \item \code{sd}: Posterior standard deviation (standard error).
#'       \item \code{mse}: Posterior Mean Squared Error (\eqn{\text{sd}^2}).
#'       \item \code{rse}: Relative Standard Error (\%).
#'       \item \code{ci_lower}: 2.5\% quantile of posterior credible interval.
#'       \item \code{ci_upper}: 97.5\% quantile of posterior credible interval.
#'       \item \code{random_effect_area}: Posterior mean of area random effect (\eqn{\hat{v}_i}).
#'       \item \code{random_effect_subarea}: Posterior mean of sub-area random effect (\eqn{\hat{u}_{ij}}).
#'     }
#'   \item \code{df_area}: Data frame with aggregate primary area-level estimates:
#'     \itemize{
#'       \item \code{domain}: Primary area identifier.
#'       \item \code{hb_area}: Aggregate area mean estimate.
#'       \item \code{sd_area}: Posterior standard error of area estimate.
#'       \item \code{mse_area}: Posterior Mean Squared Error of area estimate.
#'       \item \code{rse_area}: Relative Standard Error of area estimate (\%).
#'       \item \code{ci_lower_area}: 2.5\% quantile of area credible interval.
#'       \item \code{ci_upper_area}: 97.5\% quantile of area credible interval.
#'       \item \code{n_subareas}: Number of sub-areas in the area.
#'     }
#'   \item \code{df_subarea}: Alias pointing to \code{df_hb}.
#'   \item \code{estcoef}: Data frame of estimated regression coefficients (\eqn{\beta}, standard error, z-value, p-value, credible intervals).
#'   \item \code{hyperpar}: Data frame of hyperparameter posterior estimates.
#'   \item \code{random_effect_var}: Named vector with the posterior means of the area
#'     (\code{sigma2_v}) and sub-area (\code{sigma2_u}) variances. When a variance is close to
#'     0 its posterior is skewed and INLA's summary can be biased (on simulated data
#'     \code{sigma2_u} was 0.047 against an exact posterior mean of 0.072); the sub-area and
#'     area estimates are the reliable output.
#'   \item \code{phi}: Estimated spatial variance proportion (for BYM2).
#'   \item \code{goodness}: Model fit metrics (DIC, WAIC, Marginal Log-Likelihood).
#'   \item \code{family}: Response family used.
#'   \item \code{spatial}: Spatial model type used.
#'   \item \code{model}: Label of model (\code{"HB-TWOFOLD"}).
#'   \item \code{method}: Description of method.
#'   \item \code{convergence}: Logical indicating convergence.
#'   \item \code{level}: \code{"subarea"}.
#'   \item \code{fit}: Raw fitted INLA model object, or the list returned by the Python backend
#'     for the JAX engine.
#'   \item \code{engine}, \code{engine_spec}: Engine name (\code{"inla"} or \code{"jax"}) and its
#'     settings.
#'   \item \code{draws}, \code{rhat_max}, \code{device}, \code{elapsed_seconds}: With the JAX
#'     engine, the posterior draws of the sub-area parameters (rows of \code{df_hb} by draws, kept
#'     when \code{jax_engine(keep_draws = TRUE)}), the largest split R-hat (\code{convergence} is
#'     \code{rhat_max < 1.05}), the JAX device and the sampling time.
#'   \item \code{call}: Matched function call.
#' }
#'
#' @section Engines:
#' Both engines fit the same model and return the same structure. With
#' \code{jax_engine()} the area effect (iid, Besag or BYM2) and the iid sub-area effect have
#' half-normal priors on their standard deviations (on the standardised response scale for
#' \code{family = "gaussian"}) instead of the PC priors of \code{inla_engine()}, so the
#' estimates differ slightly. On simulated data (20 areas, 56 sub-areas; Gaussian, binomial,
#' Poisson and Gaussian with BYM2) the area means of the two engines correlated above 0.999 and
#' differed by a median of 0.02-0.07 posterior SD, with area SDs 0.85-1.14 times those of INLA.
#' The posterior means of the variances were larger with MCMC (area variance 0.93 against 0.68
#' for the Gaussian model): their posteriors are skewed, and INLA's summary is biased downwards
#' for small variances. \code{jax_engine(robust = TRUE)} gives Student-t area and sub-area
#' effects and \code{prior_beta = "horseshoe"} a regularized horseshoe prior on the coefficients.
#'
#' @references
#' \enumerate{
#'   \item Torabi, M., & Rao, J. N. K. (2014). On small area estimation under a sub-area
#'     level model. \emph{Journal of Multivariate Analysis}, 127, 36-55.
#'   \item Rao, J. N. K., & Molina, I. (2015). \emph{Small Area Estimation} (2nd ed.). John Wiley & Sons.
#'   \item Riebler, A., S\enc{ø}{o}rbye, S. H., Simpson, D., & Rue, H. (2016). An intuitive Bayesian spatial model
#'     for disease mapping that accounts for scaling. \emph{Statistical Methods in Medical Research}, 25(4), 1145-1165.
#' }
#'
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   library(fastsae)
#'   set.seed(42)
#'   m <- 10
#'   dat <- do.call(rbind, lapply(1:m, function(d) {
#'     nd <- sample(2:4, 1)
#'     data.frame(
#'       area = d,
#'       subarea = paste0(d, "-", seq_len(nd)),
#'       x1 = rnorm(nd),
#'       vardir = runif(nd, 0.2, 0.8)
#'     )
#'   }))
#'   dat$y <- 1 + dat$x1 + rnorm(m)[dat$area] + rnorm(nrow(dat), sd = 0.5) +
#'     rnorm(nrow(dat), sd = sqrt(dat$vardir))
#'
#'   fit <- hb_twofold(y ~ x1, vardir = "vardir", domain = "area",
#'                     subarea = "subarea", data = dat)
#'   summary(fit)
#' }
#' }
hb_twofold <- function(
  formula,
  vardir = NULL,
  domain,
  subarea = NULL,
  data,
  family = c("gaussian", "binomial", "poisson"),
  spatial = c("none", "bym2", "besag"),
  W = NULL,
  weight = NULL,
  trials = NULL,
  exposure = NULL,
  engine = "inla",
  compute_area = TRUE,
  print_result = TRUE,
  ...
) {
  call_matched <- match.call()
  pc_prec <- list(prior = "pc.prec", param = c(1, 0.01))
  engine <- .resolve_inla_engine(
    engine, missing(engine), list(...), "hb_twofold",
    defaults = list(
      prior_prec = pc_prec,
      prior_prec_subarea = pc_prec,
      prior_phi = list(prior = "pc", param = c(0.5, 0.5)),
      draws = 200L
    ),
    used = c("scale_model", "prior_prec", "prior_prec_subarea", "prior_phi", "draws")
  )
  strategy <- engine$strategy
  scale_model <- engine$scale_model
  # referenced by name in the INLA formula
  prior_prec_area <- engine$prior_prec
  prior_prec_subarea <- engine$prior_prec_subarea
  prior_phi <- engine$prior_phi

  family <- match.arg(tolower(family), choices = c("gaussian", "binomial", "poisson"))
  spatial <- match.arg(tolower(spatial), choices = c("none", "bym2", "besag"))

  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame or tibble.")
  }

  n_subareas <- nrow(data)
  if (n_subareas == 0) {
    cli::cli_abort("{.arg data} contains 0 rows.")
  }

  # 1. Extract Domain & Subarea Identifiers
  if (missing(domain) || is.null(domain)) {
    cli::cli_abort("Primary area identifier {.arg domain} must be specified.")
  }
  domain_vec <- .get_variable(data, domain)

  if (is.null(subarea)) {
    subarea_vec <- seq_len(n_subareas)
  } else {
    subarea_vec <- .get_variable(data, subarea)
  }

  unique_domains <- unique(domain_vec)
  n_domains <- length(unique_domains)

  # Integer indices for INLA
  area_id <- match(domain_vec, unique_domains)
  subarea_id <- seq_len(n_subareas)

  # 2. Extract Response and Design Matrix
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  term_labels <- attr(stats::terms(formula), "term.labels")

  # Optional variables
  vardir_vec <- if (!is.null(vardir)) .get_variable(data, vardir) else NULL
  weight_vec <- if (!is.null(weight)) .get_variable(data, weight) else NULL
  trials_vec <- if (!is.null(trials)) .get_variable(data, trials) else NULL
  exposure_vec <- if (!is.null(exposure)) .get_variable(data, exposure) else NULL

  # INLA treats a missing covariate as 0, so every subarea needs complete covariates
  .check_covariates_na(stats::model.matrix(attr(mf, "terms"), mf), "subarea")

  # 3. Likelihood Validation
  if (family == "gaussian") {
    if (is.null(vardir_vec)) {
      cli::cli_abort("Direct sampling variance {.arg vardir} must be provided for {.code family = 'gaussian'}.")
    }
    .check_vardir(y, vardir_vec, "subarea")
    # Scale for Gaussian: 1 / vardir
    scale_vec <- 1 / as.numeric(vardir_vec)
    # Unsampled subareas (y = NA) do not enter the likelihood
    scale_vec[is.na(scale_vec) | is.infinite(scale_vec) | scale_vec <= 0] <- 1
  } else if (family == "binomial") {
    if (is.null(trials_vec)) {
      cli::cli_abort("Sample size / trial count {.arg trials} must be provided for {.code family = 'binomial'}.")
    }
    scale_vec <- NULL
  } else if (family == "poisson") {
    scale_vec <- NULL
  }

  # 4. Spatial Matrix Conversion
  adj_graph <- NULL
  if (spatial != "none") {
    adj_res <- .convert_spatial_weights(W, n_domains = n_domains, spatial = spatial, domain_names = unique_domains)
    adj_graph <- adj_res$graph
  }

  if (engine$name == "jax") {
    result <- .fit_jax_twofold(
      formula = formula, mf = mf, y = y, domain_vec = domain_vec, subarea_vec = subarea_vec,
      unique_domains = unique_domains, family = family, spatial = spatial, adj_graph = adj_graph,
      vardir = vardir_vec, trials = trials_vec, exposure = exposure_vec, weight = weight_vec,
      compute_area = compute_area, engine = engine, call = call_matched
    )
    if (print_result) print(result)
    return(invisible(result))
  }
  .check_inla_installed()

  # 5. Build Data for INLA
  inla_data <- as.data.frame(mf[, term_labels, drop = FALSE])
  inla_data$y <- y
  inla_data$..area_id.. <- area_id
  inla_data$..subarea_id.. <- subarea_id

  # 6. Construct INLA Formula
  fixed_str <- if (length(term_labels) > 0) paste(term_labels, collapse = " + ") else "1"

  # Area random effect
  rand_area_str <- switch(spatial,
    "none" = "f(..area_id.., model = 'iid', hyper = list(prec = prior_prec_area))",
    "bym2" = paste0(
      "f(..area_id.., model = 'bym2', graph = adj_graph, scale.model = ",
      scale_model, ", hyper = list(prec = prior_prec_area, phi = prior_phi))"
    ),
    "besag" = paste0(
      "f(..area_id.., model = 'besag', graph = adj_graph, scale.model = ",
      scale_model, ", hyper = list(prec = prior_prec_area))"
    )
  )

  # Sub-area random effect (IID)
  rand_subarea_str <- "f(..subarea_id.., model = 'iid', hyper = list(prec = prior_prec_subarea))"

  inla_formula_str <- paste("y ~", fixed_str, "+", rand_area_str, "+", rand_subarea_str)
  inla_form <- stats::as.formula(inla_formula_str)

  # Thread control for numerical determinism
  old_inla_threads <- tryCatch(INLA::inla.getOption("num.threads"), error = function(e) NULL)
  if (!is.null(old_inla_threads)) {
    INLA::inla.setOption("num.threads", "1:1")
    on.exit(INLA::inla.setOption("num.threads", old_inla_threads), add = TRUE)
  }

  # 7. INLA Arguments
  inla_args <- list(
    formula = inla_form,
    data = inla_data,
    family = family,
    control.predictor = list(compute = TRUE, link = 1),
    control.inla = list(strategy = strategy),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, mlik = TRUE, config = compute_area),
    num.threads = 1
  )

  if (family == "gaussian") {
    inla_args$scale <- scale_vec
    inla_args$control.family <- list(hyper = list(prec = list(initial = 0, fixed = TRUE)))
  } else if (family == "binomial") {
    inla_args$Ntrials <- as.integer(round(trials_vec))
  } else if (family == "poisson") {
    if (!is.null(exposure_vec)) {
      inla_args$E <- as.numeric(exposure_vec)
    }
  }

  # inla_engine(...) settings override; control.* lists are merged
  # INLA's default "compact" mode overstated sub-area SDs by up to 10% and moved
  # posterior means by up to 0.1 against the exact posterior (numerical
  # integration over the two variances, simulated data) and varied slightly
  # between identical calls; "classic" matched it. inla_engine(inla.mode = ) overrides.
  # Gaussian only, as in hb_area(): classic is much slower for non-Gaussian data.
  if (family == "gaussian" && "inla.mode" %in% names(formals(INLA::inla))) inla_args$inla.mode <- "classic"
  inla_args <- utils::modifyList(inla_args, engine$inla_args)

  # 8. Fit INLA Model
  fit <- tryCatch(
    do.call(INLA::inla, inla_args),
    error = function(e) {
      cli::cli_abort(c(
        "Fitting two-fold model with INLA failed.",
        "x" = e$message
      ))
    }
  )

  # 9. Extract Sub-area Predictions
  pred_summary <- fit$summary.fitted.values[seq_len(n_subareas), , drop = FALSE]
  linear_pred_summary <- fit$summary.linear.predictor[seq_len(n_subareas), , drop = FALSE]

  hb_pred <- pred_summary$mean
  hb_sd <- pred_summary$sd
  hb_lower <- pred_summary[["0.025quant"]]
  hb_upper <- pred_summary[["0.975quant"]]
  eta_mean <- linear_pred_summary$mean

  # Extract random effects
  area_re_sum <- fit$summary.random$..area_id..
  area_re_mean <- if (!is.null(area_re_sum)) area_re_sum$mean[area_id] else rep(0, n_subareas)

  subarea_re_sum <- fit$summary.random$..subarea_id..
  subarea_re_mean <- if (!is.null(subarea_re_sum)) subarea_re_sum$mean[subarea_id] else rep(0, n_subareas)

  # Build subarea table
  mse_val <- hb_sd^2
  rse_val <- ifelse(abs(hb_pred) < .Machine$double.eps, NA_real_, (hb_sd / abs(hb_pred)) * 100)

  df_hb <- data.frame(
    domain = domain_vec,
    subarea = subarea_vec,
    y = y,
    hb = hb_pred,
    linear_pred = eta_mean,
    vardir = if (!is.null(vardir_vec)) vardir_vec else rep(NA_real_, n_subareas),
    sd = hb_sd,
    mse = mse_val,
    rse = rse_val,
    ci_lower = hb_lower,
    ci_upper = hb_upper,
    random_effect_area = area_re_mean,
    random_effect_subarea = subarea_re_mean,
    stringsAsFactors = FALSE
  )

  # 10. Area-Level Aggregation
  df_area <- NULL
  if (compute_area) {
    norm_weight <- .twofold_weights(weight_vec, domain_vec, unique_domains, subarea_vec)
    # Joint posterior draws of the sub-area parameters on the response scale
    mat_samples <- tryCatch({
      ps <- .inla_draws(fit, engine$draws, engine$seed)
      pred_names <- paste0("Predictor:", seq_len(n_subareas))
      latent_rows <- rownames(ps[[1]]$latent)
      pred_indices <- match(pred_names, latent_rows)

      if (!any(is.na(pred_indices))) {
        mat_samples <- sapply(ps, function(s) s$latent[pred_indices, 1])

        # Link function transformation if binomial or poisson
        if (family == "binomial") {
          mat_samples <- 1 / (1 + exp(-mat_samples))
        } else if (family == "poisson") {
          mat_samples <- exp(mat_samples)
        }
        mat_samples
      } else {
        NULL
      }
    }, error = function(e) NULL)
    df_area <- .twofold_area_table(mat_samples, norm_weight, domain_vec, unique_domains, hb_pred, hb_sd)
  }

  # 11. Fixed Effects Coefficients
  fixed_sum <- fit$summary.fixed
  estcoef <- data.frame(
    beta = fixed_sum$mean,
    std.error = fixed_sum$sd,
    zvalue = fixed_sum$mean / fixed_sum$sd,
    pvalue = 2 * stats::pnorm(abs(fixed_sum$mean / fixed_sum$sd), lower.tail = FALSE),
    ci_lower = fixed_sum[["0.025quant"]],
    ci_upper = fixed_sum[["0.975quant"]],
    row.names = rownames(fixed_sum)
  )

  # 12. Hyperparameters
  hyp_sum <- fit$summary.hyperpar
  s2v <- NA_real_  # area variance
  s2u <- NA_real_  # subarea variance
  phi_val <- NA_real_

  if (!is.null(hyp_sum) && nrow(hyp_sum) > 0) {
    hyp_names <- rownames(hyp_sum)

    # Area precision
    area_idx <- grep("..area_id..", hyp_names)
    if (length(area_idx) > 0) {
      prec_v <- hyp_sum$mean[area_idx[1]]
      s2v <- .post_mean_variance(fit, hyp_names[area_idx[1]], prec_v)
    }

    # Sub-area precision
    subarea_idx <- grep("..subarea_id..", hyp_names)
    if (length(subarea_idx) > 0) {
      prec_u <- hyp_sum$mean[subarea_idx[1]]
      s2u <- .post_mean_variance(fit, hyp_names[subarea_idx[1]], prec_u)
    }

    # Spatial mixing parameter (BYM2)
    phi_idx <- grep("Phi", hyp_names, ignore.case = TRUE)
    if (length(phi_idx) > 0) {
      phi_val <- hyp_sum$mean[phi_idx[1]]
    }
  }

  random_effect_var <- c(sigma2_v = s2v, sigma2_u = s2u)

  # Goodness of fit
  goodness <- list(
    dic = fit$dic$dic,
    p_eff_dic = fit$dic$p.eff,
    waic = fit$waic$waic,
    p_eff_waic = fit$waic$p.eff,
    mlik = fit$mlik[1, 1]
  )

  model_label <- paste0(
    "HB-TWOFOLD-", toupper(family),
    if (spatial != "none") paste0(" (", toupper(spatial), ")") else " (Non-spatial)"
  )

  result <- list(
    df_hb = df_hb,
    df_subarea = df_hb,
    df_area = df_area,
    estcoef = estcoef,
    hyperpar = hyp_sum,
    random_effect_var = random_effect_var,
    phi = phi_val,
    goodness = goodness,
    family = family,
    spatial = spatial,
    model = model_label,
    method = paste0("INLA (", strategy, ")"),
    convergence = .inla_convergence(fit),
    inla_mode_status = fit$mode$mode.status,
    level = "subarea",
    formula = formula,
    fit = fit,
    engine = "inla",
    engine_spec = engine,
    call = call_matched
  )

  class(result) <- c("fastsae_hb_twofold", "fastsae_hb", "fastsae")

  if (print_result) {
    print(result)
  }

  return(invisible(result))
}

#' Sub-area weights normalised within each area
#'
#' Missing or negative weights get weight 0 (with a warning); an area whose
#' weights sum to 0, and every area when `weight_vec` is NULL, gets equal weights.
#' @noRd
.twofold_weights <- function(weight_vec, domain_vec, unique_domains, subarea_vec) {
  if (!is.null(weight_vec)) {
    bad_w <- which(is.na(weight_vec) | weight_vec < 0)
    if (length(bad_w) > 0) {
      cli::cli_warn(c(
        "{.arg weight} is missing or negative for {length(bad_w)} subarea{?s}: {.val {utils::head(subarea_vec[bad_w], 10)}}.",
        "i" = "They get weight 0 in the area aggregates of {.code df_area}."
      ))
    }
  }
  norm_weight <- numeric(length(domain_vec))
  for (d in unique_domains) {
    idx <- which(domain_vec == d)
    if (is.null(weight_vec)) {
      norm_weight[idx] <- 1 / length(idx)
    } else {
      w_sub <- weight_vec[idx]
      w_sub[is.na(w_sub) | w_sub < 0] <- 0
      sum_w <- sum(w_sub)
      norm_weight[idx] <- if (sum_w > 0) w_sub / sum_w else 1 / length(idx)
    }
  }
  norm_weight
}

#' Area aggregates of a two-fold fit
#'
#' `draws` is the n x S matrix of joint posterior draws of the sub-area
#' parameters on the response scale; the area SD and interval come from the
#' weighted sums of the draws, which keep the covariance the sub-areas of an
#' area share through v_i. Without draws (`NULL`) they are approximated from
#' independent sub-area errors.
#' @noRd
.twofold_area_table <- function(draws, norm_weight, domain_vec, unique_domains, hb_pred, hb_sd) {
  n_domains <- length(unique_domains)
  area_hb <- as.numeric(tapply(norm_weight * hb_pred, domain_vec, sum, na.rm = TRUE)[as.character(unique_domains)])
  n_subs <- as.integer(table(factor(domain_vec, levels = unique_domains)))
  area_sd <- area_lower <- area_upper <- rep(NA_real_, n_domains)

  for (i in seq_along(unique_domains)) {
    idx <- which(domain_vec == unique_domains[i])
    w_i <- norm_weight[idx]
    if (!is.null(draws)) {
      area_draws <- colSums(w_i * draws[idx, , drop = FALSE])
      area_sd[i] <- stats::sd(area_draws)
      area_lower[i] <- stats::quantile(area_draws, probs = 0.025)
      area_upper[i] <- stats::quantile(area_draws, probs = 0.975)
    } else {
      area_sd[i] <- sqrt(sum((w_i * hb_sd[idx])^2))
      area_lower[i] <- area_hb[i] - 1.96 * area_sd[i]
      area_upper[i] <- area_hb[i] + 1.96 * area_sd[i]
    }
  }

  data.frame(
    domain = unique_domains,
    hb_area = area_hb,
    sd_area = area_sd,
    mse_area = area_sd^2,
    rse_area = ifelse(abs(area_hb) < .Machine$double.eps, NA_real_, (area_sd / abs(area_hb)) * 100),
    ci_lower_area = area_lower,
    ci_upper_area = area_upper,
    n_subareas = n_subs,
    stringsAsFactors = FALSE
  )
}
