#' Empirical Best Linear Unbiased Prediction based on a Fay-Herriot Model.
#'
#' @description This function gives the Empirical Best Linear Unbiased Prediction (EBLUP) or Empirical Best (EB) predictor under normality based on a Fay-Herriot model.
#'
#' @references
#' \enumerate{
#'  \item Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley & Sons.
#' }
#'
#' @param formula an object of class formula that contains a description of the model to be fitted.
#' @param data a data frame or a data frame extension (e.g. a tibble).
#' @param vardir vector or column names from data that contain variance sampling from the direct estimator.
#' @param domain vector, column name or one-sided formula referencing a domain names column
#'   in \code{data}. If NULL, the domains are numbered consecutively.
#' @param method Character string specifying the estimation method: \code{"REML"} (default) or \code{"ML"}.
#' @param transform Character string specifying data transformation for the response:
#'   \code{"none"} (default) or \code{"log"} (Log-Fay-Herriot model with Slud & Maiti 2006
#'   second-order bias-corrected back-transformation).
#' @param maxiter maximum number of iterations allowed in the Fisher-scoring algorithm.
#' @param precision convergence tolerance limit for the Fisher-scoring algorithm.
#' @param print_result print coefficient or not, default value is TRUE.
#' @param self_benchmark Logical. If \code{TRUE}, fits a self-benchmarking model via the
#'   augmented model approach (Wang, Fuller, and Qu, 2008; Rao and Molina, 2015, Sec 10.5.1),
#'   ensuring that the weighted aggregate of small area estimates automatically matches the direct
#'   survey aggregate or specified benchmark targets without post-hoc adjustment.
#' @param benchmark_weight Optional vector, column name, or one-sided formula referencing the domain
#'   benchmarking weights in \code{data} (e.g. population sizes or shares). Defaults to equal weights.
#' @param benchmark_target Optional numeric scalar or named vector specifying fixed aggregate
#'   benchmark target(s). If \code{NULL} (default), calibrates to the weighted aggregate of the direct estimates.
#' @param benchmark_group Optional vector, column name, or one-sided formula referencing the
#'   grouping / stratum column in \code{data} (e.g. province) for group-specific self-benchmarking.
#'
#' @returns The function returns a list with the following objects:
#' \code{estcoef} a data frame with the estimated model coefficients,
#' \code{random_effect_var} estimated random effect variance,
#' \code{goodness} vector containing several goodness-of-fit measures,
#' \code{df_eblup} a data frame that contains y, eblup, random_effect, vardir, mse, and rse. \cr
#'    * \code{y} variable response \cr
#'    * \code{eblup} estimated results for each area \cr
#'    * \code{random_effect} random effect for each area \cr
#'    * \code{vardir} variance sampling from the direct estimator for each area \cr
#'    * \code{mse} Mean Square Error \cr
#'    * \code{rse} Relative Standart Error (%) \cr
#' \code{self_benchmark} logical indicating whether self-benchmarking was active, \cr
#' \code{benchmark_summary} summary table of benchmark calibration (if self-benchmarking). \cr
#'
#' @details
#' The model has a form that is response ~ auxiliary variables.
#' where numeric type response variables can contain NA.
#' When the response variable contains NA it will be estimated with synthetic estimator.
#'
#' @references
#' \enumerate{
#'  \item Rao, J. N. K., & Molina, I. (2015). Small area estimation. John Wiley & Sons.
#'  \item Slud, E. V., & Maiti, T. (2006). Small-area estimation via Fay-Herriot models with logit and log transformations. \emph{Communications in Statistics - Theory and Methods}, 35(10), 1957-1971.
#'  \item Wang, J., Fuller, W. A., and Qu, Y. (2008). Small area estimation under a restriction.
#'    \emph{Survey Methodology}, 34(1), 29-36.
#' }
#'
#' @export
#' @examples
#' library(fastsae)
#'
#' # Standard Fay-Herriot model
#' m1 <- eblup_fh(
#'   y ~ x1 + x2 + x3,
#'   data = mys,
#'   vardir = "vardir"
#' )
#'
#' # Log-Fay-Herriot with Slud-Maiti bias correction
#' m_log <- eblup_fh(
#'   y ~ x1 + x2 + x3,
#'   data = mys,
#'   vardir = "vardir",
#'   transform = "log"
#' )
#'
#' # Self-benchmarking Fay-Herriot model
#' m1_bench <- eblup_fh(
#'   y ~ x1 + x2 + x3,
#'   data = mys,
#'   vardir = "vardir",
#'   self_benchmark = TRUE,
#'   benchmark_weight = "n"
#' )
#'
#' @md
eblup_fh <- function(
  formula,
  vardir,
  domain = NULL,
  data,
  method = c("REML", "ML"),
  transform = c("none", "log"),
  maxiter = 100,
  precision = 1e-4,
  print_result = TRUE,
  self_benchmark = FALSE,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL
) {
  method <- match.arg(method, choices = c("REML", "ML"))
  transform <- match.arg(transform, choices = c("none", "log"))

  if (is.null(domain)) {
    domain <- 1:nrow(data)
  } else {
    domain <- .get_variable(data, domain)
  }

  # model frame & validasi
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)

  if (inherits(vardir, "gvf_smooth") || inherits(vardir, "fastsae_gvf")) {
    vardir <- vardir$smooth_vardir
  } else {
    vardir <- .get_variable(data, vardir)
  }

  if (nrow(mf) != length(vardir)) {
    cli::cli_abort("Length of 'vardir' must equal number of observations in data ({nrow(mf)} vs {length(vardir)}).")
  }

  y <- stats::model.response(mf, "numeric")
  y_orig <- y
  vardir_orig <- vardir

  if (transform == "log") {
    if (any(y[!is.na(y)] <= 0)) {
      cli::cli_abort("Log transformation requires strictly positive response values (y > 0).")
    }
    y_log <- ifelse(is.na(y), NA_real_, log(y))
    vardir_log <- ifelse(is.na(y), vardir, vardir / (y^2))
    y <- y_log
    vardir <- vardir_log
  }

  # Self-benchmarking setup (Wang, Fuller, and Qu, 2008; Rao and Molina, 2015, Sec. 10.5.1)
  # ponytail: for log-FH benchmark target is on original scale (user-facing), but z-weight uses vardir_log;
  # Slud-Maiti back-transform is applied after fitting, so benchmark holds approximately on original scale.
  sb_info <- NULL
  if (isTRUE(self_benchmark)) {
    sb_info <- .setup_self_benchmark(
      data = data,
      formula = formula,
      y = if (transform == "log") y_orig else y,
      vardir = vardir,
      weight = benchmark_weight,
      target = benchmark_target,
      group = benchmark_group
    )
    data <- sb_info$data_aug
    formula <- sb_info$formula_aug
    mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
    y <- if (transform == "log") ifelse(is.na(stats::model.response(mf, "numeric")), NA_real_, log(stats::model.response(mf, "numeric"))) else stats::model.response(mf, "numeric")
  }

  X <- stats::model.matrix(attr(mf, "terms"), mf)

  # Check auxiliary variables for NA values
  if (anyNA(X)) {
    cli::cli_abort("Auxiliary variables contain NA values.")
  }

  res <- .eblup_core(
    Xall = X,
    yall = y,
    vardirall = vardir,
    method = method,
    maxiter = maxiter,
    precision = precision
  )

  # attach metadata
  row.names(res$estcoef) <- colnames(X)
  res$formula <- formula
  res$model <- if (transform == "log") "Log-FH" else "FH"
  res$transform <- transform

  if (transform == "log") {
    # Slud-Maiti (2006) second-order bias-corrected back-transformation
    sigma2_u <- res$random_effect_var
    mu_star <- res$df_eblup$eblup
    mse_star <- res$df_eblup$mse
    gamma_star <- ifelse(is.na(vardir) | vardir <= 0, 0.0, sigma2_u / (sigma2_u + vardir))

    eblup_orig <- exp(mu_star + 0.5 * sigma2_u * (1.0 - gamma_star))
    mse_orig <- (eblup_orig^2) * mse_star
    rse_orig <- ifelse(eblup_orig < .Machine$double.eps, NA_real_, sqrt(mse_orig) / eblup_orig * 100)

    res$df_eblup$eblup_log <- mu_star
    res$df_eblup$mse_log <- mse_star
    res$df_eblup$eblup <- eblup_orig
    res$df_eblup$mse <- mse_orig
    res$df_eblup$rse <- rse_orig
    res$df_eblup$y <- y_orig
    res$df_eblup$vardir <- vardir_orig
  }

  # Add domain identifier as first column in df_eblup
  res$df_eblup$domain <- domain
  cols_order <- c("domain", "y", "eblup", "vardir", "random_effect", "mse", "rse", "eblup_log", "mse_log")
  res$df_eblup <- res$df_eblup[, intersect(cols_order, names(res$df_eblup))]

  if (!is.null(sb_info)) {
    res$self_benchmark <- TRUE
    res$benchmark_weight <- benchmark_weight
    res$benchmark_target <- benchmark_target
    res$benchmark_group <- benchmark_group
    res$benchmark_summary <- .compute_benchmark_summary(
      unique_groups = sb_info$unique_groups,
      group_vec = sb_info$group_vec,
      w_norm = sb_info$w_norm,
      y = sb_info$original_y,
      target_vec = sb_info$target_vec,
      estimates = res$df_eblup$eblup
    )
    res$df_eblup$y <- sb_info$original_y
  }

  res$call <- match.call()
  res$data <- data
  class(res) <- "fastsae"

  if (!res$convergence) {
    cli::cli_alert_danger("After {res$n_iter} iterations, there is no convergence.")
    return(res)
  }

  if (print_result) {
    print(res)
  }
  return(res)
}
