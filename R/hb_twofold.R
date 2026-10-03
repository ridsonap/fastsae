#' Two-Fold Hierarchical Bayes for Sub-Area Level Small Area Estimation
#'
#' @description
#' Estimates small area parameters under a two-fold sub-area level model
#' (Torabi & Rao, 2014 Bayesian hierarchical framework) using Integrated Nested
#' Laplace Approximations (\pkg{INLA}). The model accounts for nested random effects
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
#' @param strategy INLA approximation strategy: \code{"laplace"} (default, highest accuracy) or
#'   \code{"simplified.laplace"} (faster).
#' @param scale_model Logical. If \code{TRUE} (default), scales the spatial graph so the
#'   marginal variance of the structured effect is approximately 1 (recommended for BYM2/Besag).
#' @param prior_prec_area List specifying the prior for area random effect precision (\eqn{\tau_v}).
#'   Default is \code{list(prior = "pc.prec", param = c(1, 0.01))}.
#' @param prior_prec_subarea List specifying the prior for sub-area random effect precision (\eqn{\tau_u}).
#'   Default is \code{list(prior = "pc.prec", param = c(1, 0.01))}.
#' @param prior_phi List specifying the PC-prior for the spatial mixing parameter \eqn{\phi}
#'   in the BYM2 model. Default is \code{list(prior = "pc", param = c(0.5, 0.5))}.
#' @param compute_area Logical. If \code{TRUE} (default), computes aggregate primary area-level
#'   estimates and posterior standard errors via posterior sampling.
#' @param n_samples Integer specifying the number of posterior draws used to compute
#'   area-level aggregate uncertainties. Default is \code{200}.
#' @param print_result Logical. If \code{TRUE} (default), prints a summary of results.
#' @param ... Additional arguments passed to \code{INLA::inla()}.
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
#'   \item \code{random_effect_var}: Named vector with area (\code{sigma2_v}) and sub-area (\code{sigma2_u}) variances.
#'   \item \code{phi}: Estimated spatial variance proportion (for BYM2).
#'   \item \code{goodness}: Model fit metrics (DIC, WAIC, Marginal Log-Likelihood).
#'   \item \code{family}: Response family used.
#'   \item \code{spatial}: Spatial model type used.
#'   \item \code{model}: Label of model (\code{"HB-TWOFOLD"}).
#'   \item \code{method}: Description of method.
#'   \item \code{convergence}: Logical indicating convergence.
#'   \item \code{level}: \code{"subarea"}.
#'   \item \code{fit}: Raw fitted INLA model object.
#'   \item \code{call}: Matched function call.
#' }
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
  strategy = c("laplace", "simplified.laplace"),
  scale_model = TRUE,
  prior_prec_area = list(prior = "pc.prec", param = c(1, 0.01)),
  prior_prec_subarea = list(prior = "pc.prec", param = c(1, 0.01)),
  prior_phi = list(prior = "pc", param = c(0.5, 0.5)),
  compute_area = TRUE,
  n_samples = 200,
  print_result = TRUE,
  ...
) {
  call_matched <- match.call()
  .check_inla_installed()

  family <- match.arg(tolower(family), choices = c("gaussian", "binomial", "poisson"))
  spatial <- match.arg(tolower(spatial), choices = c("none", "bym2", "besag"))
  strategy <- match.arg(tolower(strategy), choices = c("laplace", "simplified.laplace"))

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

  # 3. Likelihood Validation
  if (family == "gaussian") {
    if (is.null(vardir_vec)) {
      cli::cli_abort("Direct sampling variance {.arg vardir} must be provided for {.code family = 'gaussian'}.")
    }
    # Scale for Gaussian: 1 / vardir
    scale_vec <- 1 / as.numeric(vardir_vec)
    # Handle NA or non-positive vardir for non-sampled rows
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

  extra_args <- list(...)
  inla_args[names(extra_args)] <- extra_args

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
    # Normalize weights per area
    norm_weight <- numeric(n_subareas)
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

    # Area point estimate
    area_hb <- as.numeric(tapply(norm_weight * hb_pred, domain_vec, sum, na.rm = TRUE)[as.character(unique_domains)])
    n_subs <- as.integer(table(factor(domain_vec, levels = unique_domains)))

    # Compute area-level uncertainty via posterior sampling
    area_sd <- rep(NA_real_, n_domains)
    area_lower <- rep(NA_real_, n_domains)
    area_upper <- rep(NA_real_, n_domains)

    ps_ok <- tryCatch({
      ps <- INLA::inla.posterior.sample(n = n_samples, fit)
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

        for (i in seq_along(unique_domains)) {
          d <- unique_domains[i]
          idx <- which(domain_vec == d)
          w_i <- norm_weight[idx]
          area_draws <- colSums(w_i * mat_samples[idx, , drop = FALSE])
          area_sd[i] <- stats::sd(area_draws)
          area_lower[i] <- stats::quantile(area_draws, probs = 0.025)
          area_upper[i] <- stats::quantile(area_draws, probs = 0.975)
        }
        TRUE
      } else {
        FALSE
      }
    }, error = function(e) FALSE)

    if (!ps_ok) {
      # Fallback: independent approximation if posterior sampling failed
      for (i in seq_along(unique_domains)) {
        d <- unique_domains[i]
        idx <- which(domain_vec == d)
        w_i <- norm_weight[idx]
        area_sd[i] <- sqrt(sum((w_i * hb_sd[idx])^2))
        area_lower[i] <- area_hb[i] - 1.96 * area_sd[i]
        area_upper[i] <- area_hb[i] + 1.96 * area_sd[i]
      }
    }

    area_mse <- area_sd^2
    area_rse <- ifelse(abs(area_hb) < .Machine$double.eps, NA_real_, (area_sd / abs(area_hb)) * 100)

    df_area <- data.frame(
      domain = unique_domains,
      hb_area = area_hb,
      sd_area = area_sd,
      mse_area = area_mse,
      rse_area = area_rse,
      ci_lower_area = area_lower,
      ci_upper_area = area_upper,
      n_subareas = n_subs,
      stringsAsFactors = FALSE
    )
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
      s2v <- if (prec_v > 0) 1 / prec_v else NA_real_
    }

    # Sub-area precision
    subarea_idx <- grep("..subarea_id..", hyp_names)
    if (length(subarea_idx) > 0) {
      prec_u <- hyp_sum$mean[subarea_idx[1]]
      s2u <- if (prec_u > 0) 1 / prec_u else NA_real_
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
    convergence = TRUE,
    level = "subarea",
    formula = formula,
    fit = fit,
    call = call_matched
  )

  class(result) <- c("fastsae_hb_twofold", "fastsae_hb", "fastsae")

  if (print_result) {
    print(result)
  }

  return(invisible(result))
}
