#' Pre-Modeling Spatial Diagnostic and Autocorrelation Tests
#'
#' @description
#' Evaluates whether spatial autocorrelation is present in direct domain estimates
#' and regression residuals, conducts Lagrange Multiplier (LM) tests (Anselin, 1988),
#' and performs comparative model fit tests between the standard Fay-Herriot (\code{\link{eblup_fh}})
#' and Spatial Fay-Herriot (\code{\link{eblup_sfh}}) models. Provides an explicit,
#' automated recommendation on whether a spatial model is statistically necessary.
#'
#' @param formula Regression formula describing the model (e.g. \code{y ~ x1 + x2}).
#' @param vardir Optional vector, column name, or formula referencing sampling variances
#'   of the direct estimator. Required for full model comparison (Standard FH vs Spatial FH).
#' @param data A \code{data.frame} or \code{sf} polygon object containing the variables in \code{formula}.
#' @param W Proximity or spatial adjacency matrix. If \code{NULL} and \code{data} is an \code{sf}
#'   polygon object, it is automatically constructed via \code{\link{create_weights}}.
#' @param alpha Significance level for hypothesis tests (default is 0.05).
#' @param alternative Alternative hypothesis for Moran's I (\code{"greater"}, \code{"two.sided"}, \code{"less"}).
#'
#' @return An object of class \code{c("fastsae_spatial_test", "list")} containing:
#'   \describe{
#'     \item{response_test}{Global Moran's I test results on direct response variable \eqn{y}.}
#'     \item{residual_test}{Moran's I, Anselin LM-Error, and Anselin LM-Lag test results on regression residuals.}
#'     \item{model_comparison}{Log-Likelihood, AIC, BIC, and Likelihood Ratio Test (LRT) comparing standard vs spatial FH (if \code{vardir} provided).}
#'     \item{recommendation}{Character string with diagnosis and logical flag \code{spatial_recommended}.}
#'   }
#'
#' @references
#' \enumerate{
#'   \item Anselin, L. (1988). \emph{Spatial Econometrics: Methods and Models}. Kluwer Academic Publishers.
#'   \item Cliff, A. D., & Ord, J. K. (1981). \emph{Spatial Processes: Models & Applications}. Pion London.
#'   \item Petrucci, A., & Salvati, N. (2006). Small area estimation for spatial data in R. \emph{Statistical Methods & Applications}, 14(3), 391-407.
#' }
#'
#' @examples
#' library(fastsae)
#' data(mys)
#' data(mys_proxmat)
#'
#' # Run spatial diagnostics on simulated survey data
#' st <- spatial_test(
#'   formula = y ~ x1 + x2,
#'   vardir = "vardir",
#'   data = mys,
#'   W = mys_proxmat
#' )
#' print(st)
#'
#' @export
spatial_test <- function(
  formula,
  vardir = NULL,
  data,
  W = NULL,
  alpha = 0.05,
  alternative = c("greater", "two.sided", "less")
) {
  call <- match.call()
  alternative <- match.arg(alternative)

  # 1. Handle W matrix
  if (is.null(W)) {
    if (inherits(data, "sf")) {
      cli::cli_alert_info("Argument {.arg W} not provided; automatically constructing Queen spatial weights from {.cls sf} polygon geometry.")
      W <- create_weights(data, method = "queen", style = "W")
    } else {
      cli::cli_abort("Spatial weight matrix {.arg W} must be provided when {.arg data} is not an {.cls sf} polygon object.")
    }
  }

  W_mat <- as.matrix(W)
  n_total <- nrow(data)

  if (nrow(W_mat) != n_total || ncol(W_mat) != n_total) {
    cli::cli_abort(c(
      "Dimensions of {.arg W} ({nrow(W_mat)} x {ncol(W_mat)}) must match the number of rows in {.arg data} ({n_total})."
    ))
  }

  # 2. Extract model frame and variables
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  y <- stats::model.response(mf, "numeric")
  X <- stats::model.matrix(attr(mf, "terms"), mf)

  sampled_mask <- !is.na(y)
  n <- sum(sampled_mask)

  if (n < 5) {
    cli::cli_abort("At least 5 sampled observations are required to conduct spatial tests.")
  }

  y_sub <- y[sampled_mask]
  X_sub <- X[sampled_mask, , drop = FALSE]
  W_sub <- W_mat[sampled_mask, sampled_mask, drop = FALSE]

  if (anyNA(X_sub)) {
    cli::cli_abort("Auxiliary variables in {.arg formula} contain NA values for sampled domains.")
  }

  # Precompute spatial weight matrix constants once (O(N^2) instead of repeated O(N^3))
  s0 <- sum(W_sub)
  # tr(W^2 + W'W) == sum(W^2) + sum(W * t(W)) == s1 in Cliff & Ord (1981)
  s1 <- sum(W_sub^2) + sum(W_sub * t(W_sub))
  s2 <- sum((rowSums(W_sub) + colSums(W_sub))^2)
  E_I <- -1 / (n - 1)

  # 3. Helper: Global Moran's I computation (fast matrix-vector implementation)
  .calc_moran <- function(vec, W_mat, s0, s1, s2, E_I, alt = "greater") {
    n_pts <- length(vec)
    z <- vec - mean(vec)
    denom <- sum(z^2)
    if (s0 <= 0 || denom <= 0) {
      return(list(I = NA_real_, expected = NA_real_, sd = NA_real_, z = NA_real_, p_value = NA_real_))
    }

    # Matrix-vector multiplication O(N^2), zero large matrix allocation
    num <- sum(z * (W_mat %*% z))
    moran_I <- (n_pts / s0) * (num / denom)

    kurt <- (n_pts * sum(z^4)) / (denom^2)

    var_I <- (n_pts * ((n_pts^2 - 3 * n_pts + 3) * s1 - n_pts * s2 + 3 * s0^2) -
                kurt * ((n_pts^2 - n_pts) * s1 - 2 * n_pts * s2 + 6 * s0^2)) /
      ((n_pts - 1) * (n_pts - 2) * (n_pts - 3) * s0^2) - E_I^2

    var_I <- max(1e-10, var_I)
    z_stat <- (moran_I - E_I) / sqrt(var_I)

    p_val <- switch(
      alt,
      "greater" = stats::pnorm(z_stat, lower.tail = FALSE),
      "less" = stats::pnorm(z_stat, lower.tail = TRUE),
      "two.sided" = 2 * stats::pnorm(-abs(z_stat))
    )

    list(
      I = moran_I,
      expected = E_I,
      sd = sqrt(var_I),
      z = z_stat,
      p_value = p_val
    )
  }

  # Step A: Direct Response Test (ESDA on y)
  moran_y <- .calc_moran(y_sub, W_sub, s0, s1, s2, E_I, alt = alternative)

  # Step B: Model Residual Test (OLS on auxiliary variables)
  ols_fit <- stats::lm.fit(X_sub, y_sub)
  res_ols <- stats::residuals(ols_fit)
  moran_res <- .calc_moran(res_ols, W_sub, s0, s1, s2, E_I, alt = alternative)

  # Anselin (1988) Lagrange Multiplier (LM) Tests
  # LM-Error Test:
  # LM_err = ( (e' W e) / (e'e / n) )^2 / tr(W^2 + W'W)
  tr_T <- s1
  s2_e <- sum(res_ols^2) / n
  eWe <- as.numeric(crossprod(res_ols, W_sub %*% res_ols))

  lm_error_stat <- if (tr_T > 0 && s2_e > 0) {
    ((eWe / s2_e)^2) / tr_T
  } else {
    0.0
  }
  lm_error_pval <- stats::pchisq(lm_error_stat, df = 1, lower.tail = FALSE)

  # LM-Lag Test:
  # LM_lag = ( (e' W y) / (e'e / n) )^2 / D_lag
  # Note: Uses QR residual projection to guarantee numerical stability and robustness to rank-deficiency
  eWy <- as.numeric(crossprod(res_ols, W_sub %*% y_sub))
  WXbeta <- W_sub %*% ols_fit$fitted.values
  M_WXbeta <- qr.resid(ols_fit$qr, WXbeta)
  D_lag <- tr_T + sum(M_WXbeta^2) / s2_e

  lm_lag_stat <- if (is.finite(D_lag) && D_lag > 0 && s2_e > 0) {
    ((eWy / s2_e)^2) / D_lag
  } else {
    0.0
  }
  lm_lag_pval <- stats::pchisq(lm_lag_stat, df = 1, lower.tail = FALSE)

  # Step C: Model Comparison (Standard FH vs Spatial FH) if vardir supplied
  comp_res <- NULL
  if (!is.null(vardir)) {
    vardir_vec <- if (inherits(vardir, "gvf_smooth")) {
      vardir$smooth_vardir
    } else {
      .get_variable(data, vardir)
    }

    tryCatch({
      fit_fh <- eblup_fh(
        formula = formula,
        vardir = vardir_vec,
        data = data,
        method = "REML",
        print_result = FALSE
      )

      fit_sfh <- eblup_sfh(
        formula = formula,
        vardir = vardir_vec,
        data = data,
        W = W_mat,
        method = "REML",
        mse_method = "analytical",
        print_result = FALSE
      )

      ll_fh <- as.numeric(fit_fh$goodness["loglikelihood"])
      ll_sfh <- as.numeric(fit_sfh$goodness["loglikelihood"])

      # LRT Statistic: 2 * (ll_sfh - ll_fh)
      lrt_stat <- max(0.0, 2 * (ll_sfh - ll_fh))
      lrt_pval <- stats::pchisq(lrt_stat, df = 1, lower.tail = FALSE)

      aic_fh <- as.numeric(fit_fh$goodness["AIC"])
      aic_sfh <- as.numeric(fit_sfh$goodness["AIC"])
      delta_aic <- aic_sfh - aic_fh

      bic_fh <- as.numeric(fit_fh$goodness["BIC"])
      bic_sfh <- as.numeric(fit_sfh$goodness["BIC"])
      delta_bic <- bic_sfh - bic_fh

      rho_est <- as.numeric(fit_sfh$rho)
      sigma2_u_fh <- as.numeric(fit_fh$random_effect_var)
      sigma2_u_sfh <- as.numeric(fit_sfh$random_effect_var)

      comp_res <- list(
        ll_fh = ll_fh,
        ll_sfh = ll_sfh,
        lrt_stat = lrt_stat,
        lrt_pval = lrt_pval,
        aic_fh = aic_fh,
        aic_sfh = aic_sfh,
        delta_aic = delta_aic,
        bic_fh = bic_fh,
        bic_sfh = bic_sfh,
        delta_bic = delta_bic,
        rho = rho_est,
        sigma2_u_fh = sigma2_u_fh,
        sigma2_u_sfh = sigma2_u_sfh
      )
    }, error = function(e) {
      cli::cli_alert_warning("Could not complete model comparison (FH vs SFH): {e$message}")
    })
  }

  # Step D: Formulate Recommendation
  is_moran_y_sig <- !is.na(moran_y$p_value) && moran_y$p_value < alpha
  is_moran_res_sig <- !is.na(moran_res$p_value) && moran_res$p_value < alpha
  is_lm_err_sig <- !is.na(lm_error_pval) && lm_error_pval < alpha
  is_lrt_sig <- !is.null(comp_res) && !is.na(comp_res$lrt_pval) && comp_res$lrt_pval < alpha

  # Spatial model is recommended if residual spatial autocorrelation persists
  # OR if the comparative LRT test shows significant improvement
  spatial_rec <- (is_moran_res_sig || is_lm_err_sig || is_lrt_sig)

  reason_text <- if (spatial_rec) {
    if (is_lrt_sig && is_lm_err_sig) {
      "Residual spatial autocorrelation is significant and the Spatial Fay-Herriot model yields a statistically significant improvement over standard FH."
    } else if (is_lm_err_sig || is_moran_res_sig) {
      "Model residuals exhibit significant spatial autocorrelation; auxiliary variables do not fully capture spatial dependence."
    } else {
      "Likelihood Ratio Test and information criteria favor the Spatial Fay-Herriot model."
    }
  } else {
    if (is_moran_y_sig) {
      "Although the direct response exhibits spatial clustering, auxiliary variables adequately capture the spatial pattern (residuals are spatially random). A standard model is sufficient."
    } else {
      "Neither the response variable nor regression residuals exhibit significant spatial autocorrelation. A standard model is sufficient."
    }
  }

  out <- list(
    response_test = moran_y,
    residual_test = list(
      moran = moran_res,
      lm_error = list(statistic = lm_error_stat, p_value = lm_error_pval),
      lm_lag = list(statistic = lm_lag_stat, p_value = lm_lag_pval)
    ),
    model_comparison = comp_res,
    spatial_recommended = spatial_rec,
    reason = reason_text,
    alpha = alpha,
    n = n,
    W = W_mat,
    formula = formula,
    call = call
  )

  class(out) <- c("fastsae_spatial_test", "list")
  return(out)
}

#' @export
print.fastsae_spatial_test <- function(x, ...) {
  cli::cli_h1("Pre-Modeling Spatial Diagnostic Tests (fastsae)")

  # 1. Direct response test
  cli::cli_h2("1. Direct Response Variable (ESDA)")
  y_t <- x$response_test
  if (!is.na(y_t$I)) {
    sig_str_y <- if (y_t$p_value < 0.001) "***" else if (y_t$p_value < 0.01) "**" else if (y_t$p_value < 0.05) "*" else "ns"
    cli::cli_bullets(c(
      "*" = "Global Moran's I: {.val {round(y_t$I, 4)}} (Expected: {.val {round(y_t$expected, 4)}}, z-stat: {.val {round(y_t$z, 3)}})",
      "*" = "p-value: {.val {format.pval(y_t$p_value, digits = 4)}} [{sig_str_y}]",
      "i" = if (y_t$p_value < x$alpha) {
        "Direct estimates exhibit significant spatial clustering."
      } else {
        "Direct estimates show no significant spatial clustering."
      }
    ))
  }

  # 2. Residual tests
  cli::cli_h2("2. Regression Residuals (Controlling for Auxiliary Variables)")
  res_t <- x$residual_test
  sig_str_res <- if (res_t$moran$p_value < 0.001) "***" else if (res_t$moran$p_value < 0.01) "**" else if (res_t$moran$p_value < 0.05) "*" else "ns"
  sig_str_lm <- if (res_t$lm_error$p_value < 0.001) "***" else if (res_t$lm_error$p_value < 0.01) "**" else if (res_t$lm_error$p_value < 0.05) "*" else "ns"

  cli::cli_bullets(c(
    "*" = "Residual Moran's I: {.val {round(res_t$moran$I, 4)}} (z-stat: {.val {round(res_t$moran$z, 3)}}, p-value: {.val {format.pval(res_t$moran$p_value, digits = 4)}}) [{sig_str_res}]",
    "*" = "Anselin LM-Error Test: {.val {round(res_t$lm_error$statistic, 3)}} (df = 1, p-value: {.val {format.pval(res_t$lm_error$p_value, digits = 4)}}) [{sig_str_lm}]",
    "*" = "Anselin LM-Lag Test:   {.val {round(res_t$lm_lag$statistic, 3)}} (df = 1, p-value: {.val {format.pval(res_t$lm_lag$p_value, digits = 4)}})"
  ))

  # 3. Model comparison if available
  if (!is.null(x$model_comparison)) {
    mc <- x$model_comparison
    cli::cli_h2("3. Model Comparison (Standard FH vs Spatial FH)")
    sig_str_lrt <- if (mc$lrt_pval < 0.001) "***" else if (mc$lrt_pval < 0.01) "**" else if (mc$lrt_pval < 0.05) "*" else "ns"

    cli::cli_bullets(c(
      "*" = "Standard FH:       sigma2_u = {.val {round(mc$sigma2_u_fh, 4)}}, REML LogLik = {.val {round(mc$ll_fh, 2)}}, AIC = {.val {round(mc$aic_fh, 2)}}",
      "*" = "Spatial FH (SAR):  sigma2_u = {.val {round(mc$sigma2_u_sfh, 4)}}, rho = {.val {round(mc$rho, 4)}}, REML LogLik = {.val {round(mc$ll_sfh, 2)}}, AIC = {.val {round(mc$aic_sfh, 2)}}",
      "*" = "Likelihood Ratio:  Chisq = {.val {round(mc$lrt_stat, 3)}} (df = 1, p-value = {.val {format.pval(mc$lrt_pval, digits = 4)}}) [{sig_str_lrt}]",
      "*" = "Information Delta: Delta AIC = {.val {round(mc$delta_aic, 2)}}, Delta BIC = {.val {round(mc$delta_bic, 2)}}"
    ))
  }

  # 4. Final Recommendation
  cli::cli_h2("4. Final Recommendation")
  if (x$spatial_recommended) {
    cli::cli_alert_success("SPATIAL MODEL RECOMMENDED: Significant spatial spillover detected.")
    cli::cli_bullets(c(
      "v" = "{x$reason}",
      "i" = "Recommended functions: {.code eblup_sfh()} (Frequentist) or {.code hb_area(..., spatial = 'bym2')} (Bayesian)."
    ))
  } else {
    cli::cli_alert_info("STANDARD MODEL SUFFICIENT: No significant spatial spillover.")
    cli::cli_bullets(c(
      "v" = "{x$reason}",
      "i" = "Recommended function: {.code eblup_fh()} (avoids unnecessary g3 MSE inflation)."
    ))
  }

  invisible(x)
}
