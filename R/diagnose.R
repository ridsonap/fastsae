#' @title Diagnostic Evaluation of Small Area Estimation Models
#' @description Conducts a comprehensive diagnostic evaluation of small area estimation
#'   models (HB or EBLUP). Evaluates estimation precision, efficiency gains over direct
#'   estimators, external calibration and bias diagnostics (Brown et al., 2001), goodness-of-fit
#'   tests, and residual spatial autocorrelation (Moran's I).
#'
#' @param object A fitted model object of class \code{"fastsae"} (e.g., from \code{\link{hb_area}},
#'   \code{\link{eblup_fh}}, \code{\link{eblup_sfh}}, or \code{\link{eblup_stfh}}).
#' @param W Optional spatial proximity or adjacency matrix of dimension \eqn{D \times D}.
#'   If \code{NULL} and the fitted model contains a spatial matrix (e.g. \code{eblup_sfh} or spatial \code{hb_area}),
#'   it is automatically extracted from \code{object}.
#' @param truth Optional numeric vector containing true parameter values \eqn{\theta_d} (useful for simulation studies).
#' @param rse_threshold Numeric. Threshold for reliable Relative Standard Error (default is 25\%).
#' @param alpha_level Numeric. Significance level for hypothesis tests (default is 0.05).
#'
#' @return An object of class \code{c("fastsae_diagnose", "list")} containing:
#'   \describe{
#'     \item{model_info}{Basic model characteristics (model name, sample size, unsampled count).}
#'     \item{precision}{Summary of direct vs SAE RSE, proportion of areas with RSE < threshold, and MSE reduction ratios.}
#'     \item{brown_test}{Brown et al. (2001) diagnostics: \code{bias_regression} (descriptive OLS
#'       fit of the direct estimates on the model estimates, no test; see Details) and
#'       \code{goodness_of_fit} (the \eqn{W} statistic with its upper-tail \eqn{\chi^2} p-value).}
#'     \item{spatial_test}{Moran's I test for spatial autocorrelation in model residuals (if \code{W} is available).}
#'     \item{bayesian_metrics}{For Bayesian fits: DIC, WAIC, marginal log-likelihood and, when they are
#'       available per domain, CPO summaries and a Kolmogorov-Smirnov test of the PIT values.}
#'     \item{simulation_metrics}{Mean relative bias, mean absolute relative error, relative RMSE across
#'       domains, and coverage rate (if \code{truth} is supplied).}
#'     \item{df_diag}{Data frame combining direct estimates, SAE predictions, MSE, RSE, and residuals for diagnostics.}
#'     \item{status}{Overall assessment summary string.}
#'   }
#'
#' @details
#' \strong{Goodness of fit.} \eqn{W = \sum_d (y_d - \hat\theta_d)^2 / (\psi_d + \mathrm{MSE}_d)} over the
#' sampled domains is compared with \eqn{\chi^2_m} (Brown et al., 2001); only large values indicate
#' that the model estimates differ from the direct estimates more than their uncertainty allows, so
#' the p-value is upper-tailed (as in \code{emdi}). The statistic ignores the positive correlation
#' between \eqn{y_d} and \eqn{\hat\theta_d}, so it is conservative: under a correct model \eqn{W}
#' is typically well below its degrees of freedom.
#'
#' \strong{Bias regression.} The OLS fit \eqn{y_d = \alpha + \beta\hat\theta_d} is reported
#' without a formal test of \eqn{\alpha = 0, \beta = 1}: for shrinkage estimators such as the EBLUP
#' the slope exceeds 1 even when the model is correct (about \eqn{1/\gamma_d} for the part of the
#' variation not explained by the covariates), so that test would reject correct models. Inspect
#' the regression (and \code{autoplot(x, type = "calibration")}) for systematic deviations instead.
#'
#' \strong{Moran's I} on the residuals \eqn{y_d - \hat\theta_d} of the sampled domains uses the
#' randomisation variance of Cliff and Ord (1981), two-sided.
#'
#' @references
#' \enumerate{
#'   \item Brown, G., Chambers, R., Heady, P., & Heasman, D. (2001). Evaluation of small area estimation methods:
#'     An application to the British Labour Force Survey. \emph{ONS Internal Report}.
#'   \item Rao, J. N. K., & Molina, I. (2015). \emph{Small Area Estimation} (2nd ed.). John Wiley & Sons.
#'   \item Cliff, A. D., & Ord, J. K. (1981). \emph{Spatial Processes: Models & Applications}. Pion London.
#' }
#'
#' @examples
#' library(fastsae)
#' data(mys)
#' data(mys_proxmat)
#'
#' # 1. Fit Fay-Herriot model
#' fit_fh <- eblup_fh(y ~ x1 + x2, data = mys, vardir = ~vardir, domain = ~area)
#' diag_fh <- diagnose(fit_fh)
#' print(diag_fh)
#'
#' # 2. Fit Spatial HB model (requires INLA)
#' if (requireNamespace("INLA", quietly = TRUE)) {
#'   fit_hb <- hb_area(y ~ x1 + x2, data = mys, vardir = "vardir")
#'   diag_hb <- diagnose(fit_hb)
#'   print(diag_hb)
#' }
#'
#' @export
diagnose <- function(object,
                     W = NULL,
                     truth = NULL,
                     rse_threshold = 25,
                     alpha_level = 0.05) {
  if (!inherits(object, "fastsae")) {
    cli::cli_abort("{.arg object} must be a fitted model of class {.cls fastsae}.")
  }

  # 1. Extract data frame of estimates
  df <- .sae_table(object)
  if (is.null(df)) {
    cli::cli_abort("Could not find estimation results table ({.code df_hb} or {.code df_eblup}) in {.arg object}.")
  }

  est_col <- .sae_est_col(df)
  if (is.na(est_col)) {
    cli::cli_abort("Prediction column ({.val hb} or {.val eblup}) not found in results table.")
  }

  sae_pred <- df[[est_col]]
  direct_y <- if ("y" %in% names(df)) df$y else rep(NA_real_, nrow(df))
  vardir <- if ("vardir" %in% names(df)) df$vardir else rep(NA_real_, nrow(df))
  mse <- if ("mse" %in% names(df)) df$mse else rep(NA_real_, nrow(df))
  rse <- if ("rse" %in% names(df)) df$rse else rep(NA_real_, nrow(df))

  N_total <- nrow(df)
  sampled_mask <- !is.na(direct_y)
  N_sampled <- sum(sampled_mask)
  N_unsampled <- N_total - N_sampled

  # Direct RSE calculation where direct_y is available
  direct_rse <- rep(NA_real_, N_total)
  valid_direct_rse <- sampled_mask & !is.na(vardir) & vardir > 0 & !is.na(direct_y) & direct_y != 0
  direct_rse[valid_direct_rse] <- (sqrt(vardir[valid_direct_rse]) / abs(direct_y[valid_direct_rse])) * 100

  # Efficiency Gain: vardir / mse
  eff_ratio <- rep(NA_real_, N_total)
  valid_eff <- sampled_mask & !is.na(vardir) & !is.na(mse) & mse > 0
  eff_ratio[valid_eff] <- vardir[valid_eff] / mse[valid_eff]

  # Residuals y - theta^ and their standardised version. Var(y - theta^) is
  # not vardir: y - theta^ = e - (theta^ - theta), and the prediction error is
  # positively correlated with e, so Var(y - theta^) ~ vardir - MSE (for the
  # EBLUP vardir (1 - gamma)). Dividing by sqrt(vardir) gave an sd of ~0.76 under
  # a correct FH model (1.01 with vardir - MSE).
  residuals <- direct_y - sae_pred
  std_residuals <- rep(NA_real_, N_total)
  res_var <- vardir - mse
  valid_res <- sampled_mask & !is.na(res_var) & res_var > 0
  std_residuals[valid_res] <- residuals[valid_res] / sqrt(res_var[valid_res])

  # ----------------------------------------------------------------------------
  # 2. Precision and Efficiency Gain Summary
  # ----------------------------------------------------------------------------
  rse_clean <- rse[!is.na(rse)]
  prop_reliable <- if (length(rse_clean) > 0) mean(rse_clean < rse_threshold) * 100 else NA_real_

  eff_clean <- eff_ratio[!is.na(eff_ratio)]
  gain_pct <- if (length(eff_clean) > 0) mean(eff_clean > 1.0) * 100 else NA_real_

  precision_summary <- list(
    rse_threshold = rse_threshold,
    prop_reliable = prop_reliable,
    mean_direct_rse = if (any(!is.na(direct_rse))) mean(direct_rse, na.rm = TRUE) else NA_real_,
    median_direct_rse = if (any(!is.na(direct_rse))) stats::median(direct_rse, na.rm = TRUE) else NA_real_,
    mean_sae_rse = if (length(rse_clean) > 0) mean(rse_clean) else NA_real_,
    median_sae_rse = if (length(rse_clean) > 0) stats::median(rse_clean) else NA_real_,
    eff_ratio_summary = if (length(eff_clean) > 0) {
      c(Min = min(eff_clean),
        Q1 = stats::quantile(eff_clean, 0.25),
        Median = stats::median(eff_clean),
        Mean = mean(eff_clean),
        Q3 = stats::quantile(eff_clean, 0.75),
        Max = max(eff_clean))
    } else NULL,
    prop_gain = gain_pct
  )

  # ----------------------------------------------------------------------------
  # 3. Brown et al. (2001) Calibration & Goodness-of-Fit Tests
  # ----------------------------------------------------------------------------
  brown_res <- list(bias_regression = NULL, goodness_of_fit = NULL)

  calib_mask <- sampled_mask & !is.na(direct_y) & !is.na(sae_pred)
  if (sum(calib_mask) >= 5) {
    y_sub <- direct_y[calib_mask]
    pred_sub <- sae_pred[calib_mask]
    n_sub <- length(y_sub)

    # a. Bias regression y_dir = alpha + beta * sae_pred (descriptive only).
    # No test of alpha = 0, beta = 1: for shrinkage estimators the slope exceeds 1
    # under a correct model (y = theta + e and theta^ = gamma y + (1 - gamma) x'b give
    # slope 1 / gamma for the part of y not explained by x), so the F test rejects
    # correct models (93% at the 5% level for a correct FH model with D = 100).
    fit_ols <- stats::lm(y_sub ~ pred_sub)
    coefs <- unname(stats::coef(fit_ols))
    se <- unname(sqrt(diag(stats::vcov(fit_ols))))
    brown_res$bias_regression <- list(
      alpha = coefs[1],
      beta = coefs[2],
      se_alpha = se[1],
      se_beta = se[2],
      r_squared = summary(fit_ols)$r.squared,
      n = n_sub
    )

    # b. Chi-Square goodness-of-fit statistic W
    vd_sub <- vardir[calib_mask]
    mse_sub <- mse[calib_mask]
    valid_w <- !is.na(vd_sub) & !is.na(mse_sub) & (vd_sub + mse_sub) > 0
    if (sum(valid_w) >= 3) {
      denom <- vd_sub[valid_w] + mse_sub[valid_w]
      w_stat <- sum(((y_sub[valid_w] - pred_sub[valid_w])^2) / denom)
      df_w <- sum(valid_w)
      # Upper-tail test (as emdi::compare.fh): only a large W indicates a lack of
      # fit. W is conservative (y and theta^ are positively correlated), so a small
      # W is expected under a correct model and is no evidence against it.
      p_val_w <- stats::pchisq(w_stat, df = df_w, lower.tail = FALSE)
      brown_res$goodness_of_fit <- list(
        w_stat = w_stat,
        df = df_w,
        crit = stats::qchisq(1 - alpha_level, df = df_w),
        p_value = p_val_w,
        is_good_fit = (p_val_w >= alpha_level)
      )
    }
  }

  # ----------------------------------------------------------------------------
  # 4. Spatial Autocorrelation Test on Residuals (Moran's I)
  # ----------------------------------------------------------------------------
  spatial_res <- NULL
  W_mat <- if (!is.null(W)) {
    as.matrix(W)
  } else if (!is.null(object$W)) {
    as.matrix(object$W)
  } else {
    NULL
  }

  if (!is.null(W_mat)) {
    # Match sampled domains
    if (nrow(W_mat) == N_total) {
      res_clean <- residuals[sampled_mask]
      W_sub <- W_mat[sampled_mask, sampled_mask, drop = FALSE]
      n_spat <- length(res_clean)

      if (n_spat >= 5 && all(!is.na(res_clean))) {
        # Standard Moran's I calculation
        z <- res_clean - mean(res_clean)
        s0 <- sum(W_sub)
        if (s0 > 0) {
          num <- sum(z * (W_sub %*% z))
          denom <- sum(z^2)
          moran_I <- (n_spat / s0) * (num / denom)

          # Expected value and variance under randomisation (Cliff & Ord, 1981)
          E_I <- -1 / (n_spat - 1)
          s1 <- sum(W_sub^2) + sum(W_sub * t(W_sub))
          s2 <- sum((rowSums(W_sub) + colSums(W_sub))^2)
          kurt <- (n_spat * sum(z^4)) / (sum(z^2)^2)

          var_I <- (n_spat * ((n_spat^2 - 3 * n_spat + 3) * s1 - n_spat * s2 + 3 * s0^2) -
                      kurt * ((n_spat^2 - n_spat) * s1 - 2 * n_spat * s2 + 6 * s0^2)) /
            ((n_spat - 1) * (n_spat - 2) * (n_spat - 3) * s0^2) - E_I^2

          var_I <- max(1e-8, var_I)
          z_stat <- (moran_I - E_I) / sqrt(var_I)
          p_val_moran <- 2 * stats::pnorm(-abs(z_stat))

          spatial_res <- list(
            moran_I = moran_I,
            expected_I = E_I,
            sd_I = sqrt(var_I),
            z_stat = z_stat,
            p_value = p_val_moran,
            no_residual_autocorrelation = (p_val_moran >= alpha_level)
          )
        }
      }
    }
  }

  # ----------------------------------------------------------------------------
  # 5. Simulation Validation Metrics (if ground truth is supplied)
  # ----------------------------------------------------------------------------
  sim_metrics <- NULL
  if (!is.null(truth)) {
    truth_vec <- as.numeric(truth)
    if (length(truth_vec) == N_total) {
      # One realisation per domain: per-domain relative errors, and the relative
      # RMSE across domains (a per-domain RRMSE needs repeated samples).
      err <- sae_pred - truth_vec
      rb <- (err / truth_vec) * 100
      ok <- !is.na(err) & !is.na(truth_vec)
      rrmse <- if (any(ok)) 100 * sqrt(mean(err[ok]^2)) / mean(abs(truth_vec[ok])) else NA_real_

      # Coverage rate if confidence/credible intervals are present
      has_ci <- all(c("ci_lower", "ci_upper") %in% names(df))
      cr <- if (has_ci) {
        mean(truth_vec >= df$ci_lower & truth_vec <= df$ci_upper, na.rm = TRUE) * 100
      } else NA_real_

      sim_metrics <- list(
        mean_rb = mean(rb, na.rm = TRUE),
        mean_abs_rb = mean(abs(rb), na.rm = TRUE),
        rrmse = rrmse,
        coverage_rate = cr
      )
    }
  }

  # ----------------------------------------------------------------------------
  # 5b. Bayesian Diagnostics (WAIC, DIC, CPO, PIT) if available
  # ----------------------------------------------------------------------------
  bayesian_metrics <- NULL
  cpo_obj <- if (!is.null(object$fit) && !is.null(object$fit$cpo)) object$fit$cpo else NULL
  # CPO/PIT are per INLA data row; they match the domains of df only when the rows
  # coincide (hb_area, hb_twofold), not for hb_unit (unit rows + prediction rows)
  if (!is.null(cpo_obj) && length(cpo_obj$cpo) != N_total) cpo_obj <- NULL
  cpo_vec <- rep(NA_real_, N_total)
  pit_vec <- rep(NA_real_, N_total)

  gd <- object$goodness
  # hb_area: named vector (WAIC, pWAIC, DIC, pD, Marginal_LogLik);
  # hb_unit / hb_twofold: list(dic, p_eff_dic, waic, p_eff_waic, mlik)
  gd_get <- function(...) {
    for (nm in c(...)) {
      if (!is.null(gd) && nm %in% names(gd)) return(as.numeric(unname(gd[[nm]]))[1])
    }
    NA_real_
  }
  is_bayesian <- inherits(object, c("fastsae_hb_area", "fastsae_hb")) ||
    any(c("WAIC", "DIC", "waic", "dic") %in% names(gd))
  if (is_bayesian) {
    waic_val <- gd_get("WAIC", "waic")
    p_waic <- gd_get("pWAIC", "p_eff_waic")
    dic_val <- gd_get("DIC", "dic")
    p_d <- gd_get("pD", "p_eff_dic")
    log_mlik <- gd_get("Marginal_LogLik", "mlik")

    cpo_summary <- NULL
    pit_test <- NULL

    if (!is.null(cpo_obj)) {
      if (!is.null(cpo_obj$cpo)) {
        cpo_raw <- cpo_obj$cpo[seq_len(min(length(cpo_obj$cpo), N_total))]
        cpo_vec[seq_along(cpo_raw)] <- cpo_raw
      }
      if (!is.null(cpo_obj$pit)) {
        pit_raw <- cpo_obj$pit[seq_len(min(length(cpo_obj$pit), N_total))]
        pit_vec[seq_along(pit_raw)] <- pit_raw
      }

      failures <- if (!is.null(cpo_obj$failure)) sum(cpo_obj$failure > 0, na.rm = TRUE) else 0

      cpo_clean <- cpo_vec[!is.na(cpo_vec)]
      cpo_summary <- list(
        mean_cpo = if (length(cpo_clean) > 0) mean(cpo_clean) else NA_real_,
        min_cpo = if (length(cpo_clean) > 0) min(cpo_clean) else NA_real_,
        failures = failures
      )

      pit_clean <- pit_vec[!is.na(pit_vec) & pit_vec > 0 & pit_vec < 1]
      if (length(pit_clean) >= 5) {
        ks_res <- suppressWarnings(stats::ks.test(pit_clean, "punif", 0, 1))
        pit_test <- list(
          d_stat = unname(ks_res$statistic),
          p_value = ks_res$p.value,
          is_calibrated = (ks_res$p.value >= alpha_level)
        )
      }
    }

    bayesian_metrics <- list(
      waic = waic_val,
      pWAIC = p_waic,
      dic = dic_val,
      pD = p_d,
      log_mlik = log_mlik,
      cpo_summary = cpo_summary,
      pit_test = pit_test
    )
  }

  # ----------------------------------------------------------------------------
  # 6. Overall Evaluation Status
  # ----------------------------------------------------------------------------
  is_fit <- is.null(brown_res$goodness_of_fit) || brown_res$goodness_of_fit$is_good_fit
  has_high_precision <- is.na(prop_reliable) || prop_reliable >= 80

  status <- if (is_fit && has_high_precision) {
    "PASS: Model estimates are consistent with the direct estimates (Brown W test) and reliable."
  } else if (!is_fit) {
    "CAUTION: Model estimates differ from the direct estimates more than their uncertainty allows (Brown W test)."
  } else {
    "WARNING: High proportions of domains with RSE >= threshold."
  }

  df_diag <- data.frame(
    domain = if ("domain" %in% names(df)) df$domain else if ("area" %in% names(df)) df$area else seq_len(N_total),
    direct_y = direct_y,
    sae_pred = sae_pred,
    vardir = vardir,
    mse = mse,
    direct_rse = direct_rse,
    sae_rse = rse,
    eff_ratio = eff_ratio,
    residuals = residuals,
    std_residuals = std_residuals,
    is_sampled = sampled_mask,
    stringsAsFactors = FALSE
  )
  if (!is.null(cpo_obj)) {
    df_diag$cpo <- cpo_vec
    df_diag$pit <- pit_vec
  }

  out <- list(
    model_info = list(
      model_type = if (!is.null(object$model)) object$model else class(object)[1],
      N_total = N_total,
      N_sampled = N_sampled,
      N_unsampled = N_unsampled
    ),
    precision = precision_summary,
    brown_test = brown_res,
    spatial_test = spatial_res,
    bayesian_metrics = bayesian_metrics,
    simulation_metrics = sim_metrics,
    df_diag = df_diag,
    status = status
  )

  class(out) <- c("fastsae_diagnose", "list")
  return(out)
}

#' @export
print.fastsae_diagnose <- function(x, ...) {
  info <- x$model_info
  prec <- x$precision
  brown <- x$brown_test
  spat <- x$spatial_test
  sim <- x$simulation_metrics

  cli::cli_rule(left = cli::style_bold("fastsae Small Area Estimation Diagnostic Report"))
  cli::cli_text("{.strong Model:} {.val {info$model_type}} | {.strong Domains:} {info$N_total} (Sampled: {info$N_sampled}, Unsampled: {info$N_unsampled})")

  # 1. Precision & Efficiency Gain
  cli::cli_h2("1. Precision & Efficiency Gain")
  cli::cli_ul()
  if (!is.na(prec$prop_reliable)) {
    if (prec$prop_reliable >= 80) {
      cli::cli_alert_success("Domains with RSE < {prec$rse_threshold}%: {.strong {round(prec$prop_reliable, 1)}%} (High precision)")
    } else {
      cli::cli_alert_warning("Domains with RSE < {prec$rse_threshold}%: {.strong {round(prec$prop_reliable, 1)}%} (Caution: low precision)")
    }
  }
  if (!is.na(prec$mean_direct_rse) && !is.na(prec$mean_sae_rse)) {
    diff_rse <- prec$mean_direct_rse - prec$mean_sae_rse
    if (diff_rse > 0) {
      cli::cli_alert_success("Average RSE reduction: Direct {.val {round(prec$mean_direct_rse, 2)}%} -> SAE {.val {round(prec$mean_sae_rse, 2)}%} (Gain: {.strong {round(diff_rse, 2)}%})")
    } else {
      cli::cli_alert_info("Average Direct RSE: {.val {round(prec$mean_direct_rse, 2)}%} | Average SAE RSE: {.val {round(prec$mean_sae_rse, 2)}%}")
    }
  }
  if (!is.null(prec$eff_ratio_summary)) {
    s <- prec$eff_ratio_summary
    if (!is.na(prec$prop_gain) && prec$prop_gain >= 80) {
      cli::cli_alert_success("Variance reduction in {.strong {round(prec$prop_gain, 1)}%} of areas (MSE ratio median: {.val {round(s['Median'], 2)}}, max: {.val {round(s['Max'], 2)}})")
    } else if (!is.na(prec$prop_gain)) {
      cli::cli_alert_info("Variance reduction in {.strong {round(prec$prop_gain, 1)}%} of areas (MSE ratio median: {.val {round(s['Median'], 2)}})")
    }
  }
  cli::cli_end()

  # 2. Calibration & Goodness-of-Fit Tests (Brown et al., 2001)
  cli::cli_h2("2. Brown et al. (2001) Diagnostics")
  cli::cli_ul()
  if (!is.null(brown$goodness_of_fit)) {
    gof <- brown$goodness_of_fit
    if (gof$is_good_fit) {
      cli::cli_alert_success("Goodness-of-Fit Statistic W: W = {.val {round(gof$w_stat, 2)}} (df = {gof$df}), p-value = {.val {round(gof$p_value, 4)}} [Consistent with direct estimates]")
    } else {
      cli::cli_alert_warning("Goodness-of-Fit Statistic W: W = {.val {round(gof$w_stat, 2)}} (df = {gof$df}), p-value = {.val {round(gof$p_value, 4)}} [Differs from direct estimates]")
    }
  }
  if (!is.null(brown$bias_regression)) {
    br <- brown$bias_regression
    cli::cli_alert_info("Bias regression (direct ~ model, descriptive): {.field alpha} = {round(br$alpha, 4)} (SE {round(br$se_alpha, 4)}), {.field beta} = {round(br$beta, 4)} (SE {round(br$se_beta, 4)}), R^2 = {round(br$r_squared, 3)}")
    cli::cli_text("  A slope above 1 is expected for shrinkage estimators; see ?diagnose.")
  } else {
    cli::cli_alert_info("Bias regression: insufficient sampled observations.")
  }
  cli::cli_end()

  # 3. Residual Spatial Autocorrelation
  if (!is.null(spat)) {
    cli::cli_h2("3. Residual Spatial Autocorrelation")
    cli::cli_ul()
    if (spat$no_residual_autocorrelation) {
      cli::cli_alert_success("Moran's I on Residuals: I = {.val {round(spat$moran_I, 4)}} (Expected: {.val {round(spat$expected_I, 4)}}, p-value = {.val {round(spat$p_value, 4)}}) [No Residual Spatial Autocorrelation]")
    } else {
      cli::cli_alert_warning("Moran's I on Residuals: I = {.val {round(spat$moran_I, 4)}} (Expected: {.val {round(spat$expected_I, 4)}}, p-value = {.val {round(spat$p_value, 4)}}) [Residual Spatial Pattern Remains]")
    }
    cli::cli_end()
  }

  # 4. Simulation Ground Truth (if provided)
  if (!is.null(sim)) {
    cli::cli_h2("4. Simulation Ground Truth Evaluation")
    cli::cli_ul()
    if (abs(sim$mean_rb) < 5) {
      cli::cli_alert_success("Mean Relative Bias (MRB): {.val {round(sim$mean_rb, 2)}%} (Low bias < 5%)")
    } else {
      cli::cli_alert_info("Mean Relative Bias (MRB): {.val {round(sim$mean_rb, 2)}%}")
    }
    cli::cli_alert_info("Mean absolute relative error: {.val {round(sim$mean_abs_rb, 2)}%} | Relative RMSE across domains: {.val {round(sim$rrmse, 2)}%}")
    if (!is.na(sim$coverage_rate)) {
      if (sim$coverage_rate >= 90) {
        cli::cli_alert_success("95% CI Coverage Rate: {.val {round(sim$coverage_rate, 1)}%} (Nominal target met)")
      } else {
        cli::cli_alert_warning("95% CI Coverage Rate: {.val {round(sim$coverage_rate, 1)}%} (Below nominal 90%)")
      }
    }
    cli::cli_end()
  }

  # 5. Bayesian Model Diagnostics (WAIC, DIC, CPO/PIT)
  if (!is.null(x$bayesian_metrics)) {
    bm <- x$bayesian_metrics
    cli::cli_h2("5. Bayesian Information Criteria & Predictive Diagnostics")
    cli::cli_ul()
    if (!is.na(bm$waic) || !is.na(bm$dic)) {
      cli::cli_alert_info("WAIC: {.val {round(bm$waic, 2)}} (p_eff: {round(bm$pWAIC, 2)}) | DIC: {.val {round(bm$dic, 2)}} (p_eff: {round(bm$pD, 2)})")
    }
    if (!is.na(bm$log_mlik)) {
      cli::cli_alert_info("Marginal Log-Likelihood: {.val {round(bm$log_mlik, 2)}}")
    }
    if (!is.null(bm$pit_test)) {
      if (bm$pit_test$is_calibrated) {
        cli::cli_alert_success("PIT Calibration Test vs Uniform(0,1): D = {.val {round(bm$pit_test$d_stat, 3)}}, p-value = {.val {round(bm$pit_test$p_value, 4)}} [Well-calibrated predictive distribution]")
      } else {
        cli::cli_alert_warning("PIT Calibration Test vs Uniform(0,1): D = {.val {round(bm$pit_test$d_stat, 3)}}, p-value = {.val {round(bm$pit_test$p_value, 4)}} [Potential dispersion/skewness deviation]")
      }
    }
    if (!is.null(bm$cpo_summary)) {
      if (bm$cpo_summary$failures == 0) {
        cli::cli_alert_success("Leave-One-Out CPO: No numerical approximation issues (min CPO = {.val {round(bm$cpo_summary$min_cpo, 4)}})")
      } else {
        cli::cli_alert_warning("Leave-One-Out CPO: {bm$cpo_summary$failures} observation(s) flagged with numerical caution.")
      }
    }
    cli::cli_end()
  }

  cli::cli_rule()
  if (grepl("^PASS", x$status)) {
    cli::cli_alert_success("{.strong Final Assessment:} {x$status}")
  } else if (grepl("^CAUTION", x$status)) {
    cli::cli_alert_warning("{.strong Final Assessment:} {x$status}")
  } else {
    cli::cli_alert_danger("{.strong Final Assessment:} {x$status}")
  }
  cli::cli_rule()

  invisible(x)
}
