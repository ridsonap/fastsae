#' Generalized Variance Function (GVF) Smoothing for Sampling Variances
#'
#' @description Fits Generalized Variance Function (GVF) models to smooth unstable
#' direct sampling variance estimates (\eqn{D_i}) in Small Area Estimation.
#' Supports log-linear regression (Wolter, 2007), Gamma generalized linear models
#' (GLM with log link), and CV power models.
#'
#' @param vardir vector, column name, or one-sided formula referencing the direct sampling variances in \code{data}.
#' @param y optional vector, column name, or one-sided formula referencing the direct survey estimates.
#' @param n optional vector, column name, or one-sided formula referencing domain sample sizes.
#' @param method character string specifying the GVF smoothing model:
#'   \describe{
#'     \item{\code{"log_linear"}}{Log-linear model: \eqn{\log(D_i) = \alpha_0 + \alpha_1 \log(y_i) + \alpha_2 \log(n_i) + \varepsilon_i}.
#'       Includes second-order bias correction upon exponentiation.}
#'     \item{\code{"gamma_glm"}}{Gamma GLM with log link: \eqn{\mathrm{E}[D_i] = \mu_i, \log(\mu_i) = \alpha_0 + \dots}.
#'       Naturally suited for strictly positive variances without back-transformation bias.}
#'     \item{\code{"cv_power"}}{Power model on direct coefficients of variation: \eqn{\log(CV_i^2) = \alpha_0 + \dots}.}
#'   }
#'   Default is \code{"log_linear"}.
#' @param data optional data frame containing the variables.
#' @param formula optional user-supplied formula overriding the default model specification.
#' @param bias_correction logical. For \code{"log_linear"}, whether to apply log-normal bias
#'   correction \eqn{\tilde{D}_i = \exp(\hat{\mu}_i + \hat{\sigma}^2 / 2)}. Default is \code{TRUE}.
#'
#' @returns An object of class \code{"gvf_smooth"} and \code{"fastsae_gvf"}, containing:
#'   \describe{
#'     \item{\code{smooth_vardir}}{numeric vector of smoothed sampling variances \eqn{\tilde{D}_i}.}
#'     \item{\code{vardir}}{original raw sampling variances.}
#'     \item{\code{fitted_model}}{the fitted \code{lm} or \code{glm} model object.}
#'     \item{\code{method}}{method used.}
#'     \item{\code{r_squared}}{coefficient of determination \eqn{R^2} (or pseudo-\eqn{R^2} for GLM).}
#'     \item{\code{df}}{data frame with raw \code{vardir}, \code{smooth_vardir}, and covariates.}
#'     \item{\code{call}}{the matched function call.}
#'   }
#'
#' @references
#' \enumerate{
#'   \item Wolter, K. M. (2007). \emph{Introduction to Variance Estimation} (2nd ed.). Springer.
#'   \item Rao, J. N. K., & Molina, I. (2015). \emph{Small Area Estimation} (2nd ed.). John Wiley & Sons.
#'   \item Rivest, L. P., & Belmonte, E. (2000). A conditional mean squared error of small area estimators. \emph{Survey Methodology}, 26(1), 67-78.
#' }
#'
#' @export
#' @examples
#' library(fastsae)
#'
#' # Smooth direct variances using direct estimates y and sample size n
#' gvf_res <- gvf_smooth(
#'   vardir = "vardir",
#'   y = "y",
#'   n = "n",
#'   data = mys,
#'   method = "log_linear"
#' )
#'
#' # Inspect summary
#' summary(gvf_res)
#'
#' # Use smoothed variances directly in Fay-Herriot model
#' fit_fh <- eblup_fh(
#'   y ~ x1 + x2 + x3,
#'   data = mys,
#'   vardir = gvf_res$smooth_vardir
#' )
#'
#' @md
gvf_smooth <- function(
  vardir,
  y = NULL,
  n = NULL,
  method = c("log_linear", "gamma_glm", "cv_power"),
  data = NULL,
  formula = NULL,
  bias_correction = TRUE
) {
  method <- match.arg(method, choices = c("log_linear", "gamma_glm", "cv_power"))

  # Extract variables
  if (!is.null(data)) {
    vardir_vec <- .get_variable(data, vardir)
    y_vec <- if (!is.null(y)) .get_variable(data, y) else NULL
    n_vec <- if (!is.null(n)) .get_variable(data, n) else NULL
  } else {
    vardir_vec <- as.numeric(vardir)
    y_vec <- if (!is.null(y)) as.numeric(y) else NULL
    n_vec <- if (!is.null(n)) as.numeric(n) else NULL
  }

  m <- length(vardir_vec)
  if (!is.null(y_vec) && length(y_vec) != m) {
    cli::cli_abort("Length of 'y' ({length(y_vec)}) must match length of 'vardir' ({m}).")
  }
  if (!is.null(n_vec) && length(n_vec) != m) {
    cli::cli_abort("Length of 'n' ({length(n_vec)}) must match length of 'vardir' ({m}).")
  }

  # Identify valid rows for model fitting (strictly positive vardir, positive y if present)
  valid_fit <- !is.na(vardir_vec) & (vardir_vec > 0)
  if (!is.null(y_vec)) valid_fit <- valid_fit & !is.na(y_vec) & (y_vec > 0)
  if (!is.null(n_vec)) valid_fit <- valid_fit & !is.na(n_vec) & (n_vec > 0)

  if (sum(valid_fit) < 3) {
    cli::cli_abort("At least 3 valid observations with strictly positive 'vardir' are required to fit a GVF model.")
  }

  # Build data frame for fitting
  fit_df <- data.frame(vardir = vardir_vec)

  if (!is.null(y_vec)) {
    fit_df$y <- y_vec
    fit_df$log_y <- ifelse(is.na(y_vec) | y_vec <= 0, NA_real_, log(y_vec))
  }

  if (!is.null(n_vec)) {
    if (length(n_vec) != m) {
      cli::cli_abort("Length of 'n' ({length(n_vec)}) must match length of 'vardir' ({m}).")
    }
    fit_df$n <- n_vec
    fit_df$log_n <- ifelse(is.na(n_vec) | n_vec <= 0, NA_real_, log(n_vec))
    fit_df$inv_n <- ifelse(is.na(n_vec) | n_vec <= 0, NA_real_, 1 / n_vec)
  }

  # Check if at least one covariate or formula is given
  if (is.null(formula) && is.null(y_vec) && is.null(n_vec)) {
    cli::cli_abort("At least one explanatory variable ('y' or 'n') or a custom 'formula' must be provided for GVF smoothing.")
  }

  # Construct model formula if not user-supplied
  if (is.null(formula)) {
    pred_terms <- character(0)
    if (!is.null(y_vec)) pred_terms <- c(pred_terms, "log_y")
    if (!is.null(n_vec)) pred_terms <- c(pred_terms, "log_n")

    rhs <- paste(pred_terms, collapse = " + ")

    if (method == "log_linear") {
      fit_df$log_vardir <- ifelse(fit_df$vardir > 0, log(fit_df$vardir), NA_real_)
      model_formula <- stats::as.formula(paste("log_vardir ~", rhs))
    } else if (method == "gamma_glm") {
      model_formula <- stats::as.formula(paste("vardir ~", rhs))
    } else if (method == "cv_power") {
      if (is.null(y_vec)) {
        cli::cli_abort("Method 'cv_power' requires direct estimate 'y' to compute direct CV.")
      }
      fit_df$log_cv2 <- ifelse(fit_df$vardir > 0 & !is.na(fit_df$y) & fit_df$y > 0,
                               log(fit_df$vardir / (fit_df$y^2)), NA_real_)
      model_formula <- stats::as.formula(paste("log_cv2 ~", rhs))
    }
  } else {
    model_formula <- stats::as.formula(formula)
  }

  # Fit GVF model on valid observations
  train_df <- fit_df[valid_fit, , drop = FALSE]

  if (method %in% c("log_linear", "cv_power")) {
    fit_mod <- stats::lm(model_formula, data = train_df)
    r_sq <- summary(fit_mod)$r.squared
    s2 <- summary(fit_mod)$sigma^2
    bc <- if (isTRUE(bias_correction)) exp(s2 / 2) else 1.0

    pred_log <- stats::predict(fit_mod, newdata = fit_df)
    if (method == "log_linear") {
      smooth_d <- exp(pred_log) * bc
    } else {
      # cv_power
      smooth_cv2 <- exp(pred_log) * bc
      smooth_d <- smooth_cv2 * (fit_df$y^2)
    }
  } else if (method == "gamma_glm") {
    fit_mod <- stats::glm(model_formula, data = train_df, family = stats::Gamma(link = "log"))
    smooth_d <- stats::predict(fit_mod, newdata = fit_df, type = "response")
    # Pseudo R-squared: 1 - deviance / null.deviance
    r_sq <- if (!is.null(fit_mod$null.deviance) && fit_mod$null.deviance > 0) {
      1 - (fit_mod$deviance / fit_mod$null.deviance)
    } else {
      NA_real_
    }
  }

  smooth_d[is.na(smooth_d) & vardir_vec == 0] <- 0

  # Clean return data frame
  fit_df$smooth_vardir <- as.numeric(smooth_d)

  out <- list(
    smooth_vardir = as.numeric(smooth_d),
    vardir = vardir_vec,
    fitted_model = fit_mod,
    formula = model_formula,
    method = method,
    r_squared = r_sq,
    df = fit_df,
    call = match.call()
  )

  class(out) <- c("gvf_smooth", "fastsae_gvf")
  return(out)
}

#' @export
print.gvf_smooth <- function(x, ...) {
  cat("-- Generalized Variance Function (GVF) Smoothing --\n")
  cat("Method: \"", x$method, "\"\n", sep = "")
  cat("Model formula: ", deparse(x$formula), "\n", sep = "")
  cat("R-squared / Pseudo R-squared: ", round(x$r_squared, 4), "\n", sep = "")
  cat("Number of domains: ", length(x$smooth_vardir), "\n", sep = "")
  cat("Raw vardir range: [", round(min(x$vardir), 5), ", ", round(max(x$vardir), 5), "]\n", sep = "")
  cat("Smoothed vardir range: [", round(min(x$smooth_vardir), 5), ", ", round(max(x$smooth_vardir), 5), "]\n", sep = "")
  invisible(x)
}

#' @export
summary.gvf_smooth <- function(object, ...) {
  print(object)
  cat("\n-- Model Coefficients --\n")
  print(summary(object$fitted_model)$coefficients)
  invisible(object)
}

#' @export
plot.gvf_smooth <- function(x, ...) {
  graphics::plot(
    x$vardir,
    x$smooth_vardir,
    xlab = "Direct Sampling Variance (vardir)",
    ylab = "Smoothed Variance (smooth_vardir)",
    main = paste("GVF Smoothing:", x$method),
    pch = 19,
    col = "steelblue",
    ...
  )
  graphics::abline(0, 1, col = "red", lty = 2)
  invisible(x)
}
