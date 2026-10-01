#' Hierarchical Bayes for Unit-Level Small Area Estimation
#'
#' @description
#' Estimates small area parameters using unit-level Hierarchical Bayesian models
#' (Battese-Harter-Fuller model and generalized linear mixed models)
#' with Integrated Nested Laplace Approximations (\pkg{INLA}). Supports Gaussian,
#' Binomial (logistic regression for binary unit responses), and Poisson likelihoods,
#' with optional spatial random effect structures (BYM2, Besag) and finite
#' population adjustments.
#'
#' @param formula An object of class \code{formula} specifying the unit-level
#'   fixed-effects model (e.g., \code{y ~ x1 + x2}).
#' @param unit_data A \code{data.frame} containing the unit-level survey sample data.
#' @param Xpop A \code{data.frame} containing auxiliary population information for all
#'   domains (either domain population means or unit-level population records). Must contain
#'   the domain identifier specified by \code{domain_var}.
#' @param domain_var Character string specifying the column name for the domain identifier
#'   in both \code{unit_data} and \code{Xpop}.
#' @param popsize_var Optional character string specifying the column name for domain population
#'   sizes (\eqn{N_d}) in \code{Xpop}. If provided, finite population adjustment
#'   \eqn{\hat{\bar{Y}}_d = f_d \bar{y}_{d,s} + (1 - f_d) \hat{\mu}_{d,r}} is computed,
#'   where \eqn{f_d = n_d / N_d}. If \code{NULL}, the superpopulation expectation
#'   \eqn{\hat{\bar{Y}}_d = \bar{X}_d^\top \hat{\beta} + \hat{u}_d} is reported.
#' @param family Character string specifying the response likelihood. Options:
#'   \itemize{
#'     \item \code{"gaussian"}: Continuous response (Battese-Harter-Fuller model, default).
#'     \item \code{"binomial"}: Binary / Bernoulli response (0 or 1, e.g. poverty indicator).
#'     \item \code{"poisson"}: Count response (non-negative integer).
#'   }
#' @param spatial Character string specifying the spatial random effect structure across domains:
#'   \itemize{
#'     \item \code{"none"}: Non-spatial independent and identically distributed (IID) domain effects (default).
#'     \item \code{"bym2"}: Scaled Besag-York-Molli\enc{é}{e} 2 spatial model (Riebler et al., 2016).
#'     \item \code{"besag"}: Intrinsic Conditional Autoregressive (ICAR) spatial model.
#'   }
#' @param W Proximity or spatial adjacency matrix. Can be a square \code{matrix},
#'   \code{Matrix}, or \code{spdep} \code{nb} or \code{listw} object. Dimensions must
#'   match the total number of unique domains in \code{Xpop}. Required when \code{spatial != "none"}.
#' @param popnmean_xpop Optional matrix or data frame of auxiliary population means per domain.
#'   If \code{NULL} (default), means are extracted directly from \code{Xpop}.
#' @param strategy INLA approximation strategy: \code{"laplace"} (default, highest accuracy) or
#'   \code{"simplified.laplace"} (faster).
#' @param scale_model Logical. If \code{TRUE} (default), scales the spatial graph so the
#'   marginal variance of the structured effect is approximately 1 (recommended for BYM2/Besag).
#' @param prior_prec List specifying the prior for domain random effect precision. Default is
#'   \code{list(prior = "loggamma", param = c(0.01, 0.01))}.
#' @param prior_phi List specifying the PC-prior for the spatial mixing parameter \eqn{\phi}
#'   in the BYM2 model. Default is \code{list(prior = "pc", param = c(0.5, 0.5))}.
#' @param print_result Logical. If \code{TRUE} (default), prints a summary of results.
#' @param ... Additional arguments passed to \code{INLA::inla()}.
#'
#' @returns An object of class \code{c("fastsae_hb_unit", "fastsae_hb", "fastsae")} containing:
#' \itemize{
#'   \item \code{df_hb}: Data frame with domain estimates, including:
#'     \itemize{
#'       \item \code{domain}: Domain identifier.
#'       \item \code{hb}: Estimated domain mean or proportion.
#'       \item \code{linear_pred}: Posterior mean of domain linear predictor.
#'       \item \code{sd}: Posterior standard deviation (standard error).
#'       \item \code{mse}: Posterior Mean Squared Error (\eqn{\text{sd}^2}).
#'       \item \code{rse}: Relative Standard Error (\%).
#'       \item \code{ci_lower}: 2.5\% quantile of posterior credible interval.
#'       \item \code{ci_upper}: 97.5\% quantile of posterior credible interval.
#'       \item \code{samp_size}: Sample size (\eqn{n_d}) in the domain (0 for unsampled).
#'       \item \code{pop_size}: Population size (\eqn{N_d}), if \code{popsize_var} is provided.
#'       \item \code{estimated_total}: Estimated domain population total (if \code{popsize_var} is provided).
#'       \item \code{sample_mean}: Direct sample mean (\eqn{\bar{y}_{d,s}}).
#'       \item \code{random_effect}: Posterior mean of domain random effect (\eqn{\hat{u}_d}).
#'     }
#'   \item \code{estcoef}: Data frame of estimated regression coefficients (posterior mean, sd, z-value, p-value, credible intervals).
#'   \item \code{hyperpar}: Data frame of hyperparameter posterior estimates.
#'   \item \code{random_effect_var}: Estimated domain random effect variance (\eqn{\sigma_u^2}).
#'   \item \code{residual_var}: Estimated residual variance (\eqn{\sigma_e^2}, for Gaussian).
#'   \item \code{phi}: Estimated spatial variance proportion (for BYM2).
#'   \item \code{goodness}: Model fit metrics (DIC, WAIC, Marginal Log-Likelihood).
#'   \item \code{family}: Response family used.
#'   \item \code{spatial}: Spatial model type used.
#'   \item \code{unsampled_domains}: Character vector of unsampled domain identifiers.
#'   \item \code{fit}: Raw fitted INLA model object.
#'   \item \code{call}: Matched function call.
#' }
#'
#' @references
#' \enumerate{
#'   \item Battese, G. E., Harter, R. M., & Fuller, W. A. (1988). An error-components model
#'     for prediction of county crop areas using survey and satellite data.
#'     \emph{Journal of the American Statistical Association}, 83(401), 28-36.
#'   \item Rao, J. N. K., & Molina, I. (2015). \emph{Small Area Estimation} (2nd ed.). John Wiley & Sons.
#'   \item Riebler, A., S\enc{ø}{o}rbye, S. H., Simpson, D., & Rue, H. (2016). An intuitive Bayesian spatial model
#'     for disease mapping that accounts for scaling. \emph{Statistical Methods in Medical Research}, 25(4), 1145-1165.
#'   \item Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian
#'     models by using integrated nested Laplace approximations.
#'     \emph{Journal of the Royal Statistical Society: Series B}, 71(2), 319-392.
#' }
#'
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   library(fastsae)
#'   data(cornsoybean)
#'   data(cornsoybeanmeans)
#'
#'   # Prepare auxiliary population means
#'   df_pop <- cornsoybeanmeans
#'   names(df_pop)[names(df_pop) == "MeanCornPixPerSeg"] <- "CornPix"
#'   names(df_pop)[names(df_pop) == "MeanSoyBeansPixPerSeg"] <- "SoyBeansPix"
#'
#'   # Unit survey data
#'   df_sample <- cornsoybean
#'   names(df_sample)[names(df_sample) == "County"] <- "CountyIndex"
#'
#'   # Fit unit-level Hierarchical Bayes model
#'   fit_hb_u <- hb_unit(
#'     formula = CornHec ~ CornPix + SoyBeansPix,
#'     unit_data = df_sample,
#'     Xpop = df_pop,
#'     domain_var = "CountyIndex",
#'     popsize_var = "PopnSegments"
#'   )
#'   print(fit_hb_u)
#' }
#' }
hb_unit <- function(
  formula,
  unit_data,
  Xpop,
  domain_var,
  popsize_var = NULL,
  family = c("gaussian", "binomial", "poisson"),
  spatial = c("none", "bym2", "besag"),
  W = NULL,
  popnmean_xpop = NULL,
  strategy = c("laplace", "simplified.laplace"),
  scale_model = TRUE,
  prior_prec = list(prior = "loggamma", param = c(0.01, 0.01)),
  prior_phi = list(prior = "pc", param = c(0.5, 0.5)),
  print_result = TRUE,
  ...
) {
  # 1. Dependency checks
  .check_inla_installed()

  family <- match.arg(family)
  spatial <- match.arg(spatial)
  strategy <- match.arg(strategy)

  # 2. Input validation
  if (!inherits(formula, "formula")) {
    cli::cli_abort("{.arg formula} must be a valid formula object.")
  }
  if (!is.data.frame(unit_data)) {
    cli::cli_abort("{.arg unit_data} must be a data frame.")
  }
  if (!is.data.frame(Xpop)) {
    cli::cli_abort("{.arg Xpop} must be a data frame.")
  }
  if (!is.character(domain_var) || length(domain_var) != 1) {
    cli::cli_abort("{.arg domain_var} must be a single character string.")
  }
  if (!domain_var %in% names(unit_data)) {
    cli::cli_abort("Domain identifier {.val {domain_var}} not found in {.arg unit_data}.")
  }
  if (!domain_var %in% names(Xpop)) {
    cli::cli_abort("Domain identifier {.val {domain_var}} not found in {.arg Xpop}.")
  }
  if (!is.null(popsize_var) && !popsize_var %in% names(Xpop)) {
    cli::cli_abort("Population size variable {.val {popsize_var}} not found in {.arg Xpop}.")
  }

  # 3. Process unit sample data
  formuladata_s <- stats::model.frame(formula, na.action = stats::na.omit, data = unit_data)
  dom_s <- unit_data[[domain_var]]
  if (!is.null(attr(formuladata_s, "na.action"))) {
    dom_s <- dom_s[-attr(formuladata_s, "na.action")]
  }

  y_s <- stats::model.response(formuladata_s)
  n_sample <- nrow(formuladata_s)
  term_labels <- attr(stats::terms(formula), "term.labels")

  if (length(term_labels) == 0) {
    cli::cli_abort("Formula must contain at least one predictor variable.")
  }

  # 4. Process population data (domain-level aggregation)
  all_domains <- unique(Xpop[[domain_var]])
  n_domains <- length(all_domains)

  # Check if any sample domain is missing in Xpop
  missing_in_pop <- setdiff(unique(dom_s), all_domains)
  if (length(missing_in_pop) > 0) {
    cli::cli_abort(
      "The following domain(s) in {.arg unit_data} are missing from {.arg Xpop}: {paste(missing_in_pop, collapse = ', ')}"
    )
  }

  # Unsampled domains
  unsampled_domains <- setdiff(all_domains, unique(dom_s))

  # If Xpop is unit-level (more rows than unique domains), aggregate to domain means
  if (nrow(Xpop) > n_domains) {
    agg_formula <- stats::as.formula(paste("cbind(", paste(term_labels, collapse = ", "), ") ~", domain_var))
    Xpop_domain <- stats::aggregate(agg_formula, data = Xpop, FUN = mean, na.rm = TRUE)
    if (!is.null(popsize_var)) {
      pop_sizes <- as.integer(table(factor(Xpop[[domain_var]], levels = all_domains)))
      Xpop_domain[[popsize_var]] <- pop_sizes
    }
  } else {
    Xpop_domain <- Xpop
  }

  # Sort Xpop_domain to match all_domains
  Xpop_domain <- Xpop_domain[match(all_domains, Xpop_domain[[domain_var]]), , drop = FALSE]

  # Matrix of population means
  if (is.null(popnmean_xpop)) {
    formula_noy <- stats::reformulate(term_labels)
    meanxpop <- stats::model.matrix(formula_noy, data = Xpop_domain)
  } else {
    meanxpop <- as.matrix(popnmean_xpop)
    if (ncol(meanxpop) == length(term_labels)) {
      meanxpop <- cbind(`(Intercept)` = 1, meanxpop)
    }
  }

  # 5. Spatial matrix preparation
  adj_graph <- NULL
  if (spatial != "none") {
    adj_res <- .convert_spatial_weights(W, n_domains = n_domains, spatial = spatial, domain_names = all_domains)
    adj_graph <- adj_res$graph
  }

  # 6. Build combined INLA data frame
  dom_levels <- as.character(all_domains)
  unit_dom_id <- match(as.character(dom_s), dom_levels)
  pred_dom_id <- 1:n_domains

  # Sample covariates data frame
  df_s <- as.data.frame(formuladata_s[, term_labels, drop = FALSE])
  df_s$y <- y_s
  df_s$..dom_id.. <- unit_dom_id

  # Prediction rows (1 per domain, response = NA)
  df_p <- as.data.frame(Xpop_domain[, term_labels, drop = FALSE])
  df_p$y <- rep(NA_real_, n_domains)
  df_p$..dom_id.. <- pred_dom_id

  combined_df <- rbind(df_s, df_p)

  # 7. Construct INLA formula
  fixed_str <- paste(term_labels, collapse = " + ")
  rand_str <- switch(spatial,
    "none" = "f(..dom_id.., model = 'iid', hyper = list(prec = prior_prec))",
    "bym2" = paste0(
      "f(..dom_id.., model = 'bym2', graph = adj_graph, scale.model = ",
      scale_model, ", hyper = list(prec = prior_prec, phi = prior_phi))"
    ),
    "besag" = paste0(
      "f(..dom_id.., model = 'besag', graph = adj_graph, scale.model = ",
      scale_model, ", hyper = list(prec = prior_prec))"
    )
  )

  inla_form <- stats::as.formula(paste("y ~", fixed_str, "+", rand_str))

  # Thread control for numerical determinism
  old_inla_threads <- tryCatch(INLA::inla.getOption("num.threads"), error = function(e) NULL)
  if (!is.null(old_inla_threads)) {
    INLA::inla.setOption("num.threads", "1:1")
    on.exit(INLA::inla.setOption("num.threads", old_inla_threads), add = TRUE)
  }

  inla_args <- list(
    formula = inla_form,
    data = combined_df,
    family = family,
    control.predictor = list(compute = TRUE, link = 1),
    control.inla = list(strategy = strategy),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, mlik = TRUE),
    num.threads = 1
  )
  extra_args <- list(...)
  inla_args[names(extra_args)] <- extra_args

  # 8. Fit INLA model
  fit <- tryCatch(
    do.call(INLA::inla, inla_args),
    error = function(e) {
      cli::cli_abort(c(
        "Fitting model with INLA failed.",
        "x" = e$message
      ))
    }
  )

  # 9. Extract domain predictions
  pred_indices <- (n_sample + 1):nrow(combined_df)
  pred_summary <- fit$summary.fitted.values[pred_indices, , drop = FALSE]
  linear_pred_summary <- fit$summary.linear.predictor[pred_indices, , drop = FALSE]

  hb_pred <- pred_summary$mean
  hb_sd <- pred_summary$sd
  hb_lower <- pred_summary[["0.025quant"]]
  hb_upper <- pred_summary[["0.975quant"]]
  eta_mean <- linear_pred_summary$mean

  # Sample sizes and sample means per domain
  dom_factor <- factor(dom_s, levels = all_domains)
  samp_size <- as.integer(table(dom_factor))
  samp_mean_list <- tapply(y_s, dom_factor, mean, na.rm = TRUE)
  sample_mean_y <- as.numeric(samp_mean_list)

  # Finite population adjustment if popsize_var is present
  if (!is.null(popsize_var)) {
    pop_size <- as.numeric(Xpop_domain[[popsize_var]])
    f_d <- samp_size / pop_size
    f_d[is.na(f_d) | f_d < 0] <- 0
    f_d[f_d > 1] <- 1

    # Adjusted domain estimator: f_d * y_bar_s + (1 - f_d) * pred_mean
    # For domains with samp_size == 0, f_d == 0 so hb_adj = hb_pred
    hb_adj <- ifelse(samp_size > 0, f_d * sample_mean_y + (1 - f_d) * hb_pred, hb_pred)

    # Variance adjustment
    sd_adj <- ifelse(samp_size > 0, (1 - f_d) * hb_sd, hb_sd)
    # Ensure minimum standard deviation
    sd_adj[sd_adj < .Machine$double.eps] <- hb_sd[sd_adj < .Machine$double.eps]

    ci_lower_adj <- hb_adj - 1.96 * sd_adj
    ci_upper_adj <- hb_adj + 1.96 * sd_adj

    if (family == "binomial") {
      ci_lower_adj <- pmax(0, ci_lower_adj)
      ci_upper_adj <- pmin(1, ci_upper_adj)
    } else if (family == "poisson") {
      ci_lower_adj <- pmax(0, ci_lower_adj)
    }

    final_hb <- hb_adj
    final_sd <- sd_adj
    final_lower <- ci_lower_adj
    final_upper <- ci_upper_adj
  } else {
    pop_size <- rep(NA_real_, n_domains)
    final_hb <- hb_pred
    final_sd <- hb_sd
    final_lower <- hb_lower
    final_upper <- hb_upper
  }

  # 10. Extract random effects
  u_summary <- fit$summary.random$..dom_id..
  u_mean <- if (!is.null(u_summary)) u_summary$mean[1:n_domains] else rep(0, n_domains)

  # 11. Extract fixed coefficients
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
  sigma2_u <- NA_real_
  sigma2_e <- NA_real_
  phi_val <- NA_real_

  if (!is.null(hyp_sum) && nrow(hyp_sum) > 0) {
    hyp_names <- rownames(hyp_sum)

    # Random effect variance
    re_idx <- grep("..dom_id..", hyp_names)
    if (length(re_idx) > 0) {
      prec_u <- hyp_sum$mean[re_idx[1]]
      sigma2_u <- if (prec_u > 0) 1 / prec_u else NA_real_
    }

    # Observation / residual variance (Gaussian)
    if (family == "gaussian") {
      obs_idx <- grep("Gaussian", hyp_names)
      if (length(obs_idx) > 0) {
        prec_e <- hyp_sum$mean[obs_idx[1]]
        sigma2_e <- if (prec_e > 0) 1 / prec_e else NA_real_
      }
    }

    # Spatial mixing parameter (BYM2)
    phi_idx <- grep("Phi", hyp_names, ignore.case = TRUE)
    if (length(phi_idx) > 0) {
      phi_val <- hyp_sum$mean[phi_idx[1]]
    }
  }

  # Direct sampling variance for Gaussian unit level: sigma2_e / n_d
  vardir_val <- if (!is.na(sigma2_e)) {
    ifelse(samp_size > 0, sigma2_e / samp_size, NA_real_)
  } else if (family == "binomial") {
    ifelse(samp_size > 0 & !is.na(sample_mean_y), (sample_mean_y * (1 - sample_mean_y)) / samp_size, NA_real_)
  } else if (family == "poisson") {
    ifelse(samp_size > 0 & !is.na(sample_mean_y), sample_mean_y / samp_size, NA_real_)
  } else {
    rep(NA_real_, n_domains)
  }

  # 13. Build df_hb results table
  mse_val <- final_sd^2
  rse_val <- ifelse(abs(final_hb) < .Machine$double.eps, NA_real_, (final_sd / abs(final_hb)) * 100)

  df_hb <- data.frame(
    domain = all_domains,
    y = sample_mean_y,
    hb = final_hb,
    linear_pred = eta_mean,
    vardir = vardir_val,
    sd = final_sd,
    mse = mse_val,
    rse = rse_val,
    ci_lower = final_lower,
    ci_upper = final_upper,
    samp_size = samp_size,
    sample_mean = sample_mean_y,
    random_effect = u_mean,
    stringsAsFactors = FALSE
  )

  if (!is.null(popsize_var)) {
    df_hb$pop_size <- pop_size
    df_hb$estimated_total <- final_hb * pop_size
  }

  # Model goodness of fit
  goodness <- list(
    dic = fit$dic$dic,
    p_eff_dic = fit$dic$p.eff,
    waic = fit$waic$waic,
    p_eff_waic = fit$waic$p.eff,
    mlik = fit$mlik[1, 1]
  )

  model_label <- paste0(
    "HB-UNIT-", toupper(family),
    if (spatial != "none") paste0(" (", toupper(spatial), ")") else " (Non-spatial)"
  )

  result <- list(
    df_hb = df_hb,
    estcoef = estcoef,
    hyperpar = hyp_sum,
    random_effect_var = sigma2_u,
    residual_var = sigma2_e,
    phi = phi_val,
    goodness = goodness,
    family = family,
    spatial = spatial,
    unsampled_domains = unsampled_domains,
    model = model_label,
    method = paste0("INLA (", strategy, ")"),
    convergence = TRUE,
    level = "unit",
    formula = formula,
    fit = fit,
    call = match.call()
  )

  class(result) <- c("fastsae_hb_unit", "fastsae_hb", "fastsae")

  if (print_result) {
    print(result)
    if (length(unsampled_domains) > 0) {
      cli::cli_alert_info("{length(unsampled_domains)} domain(s) are unsampled (synthetic prediction applied).")
    }
  }

  return(invisible(result))
}
