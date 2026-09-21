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
#' @param method Fitting method can be chosen between 'ML' and 'REML'.
#' @param maxiter maximum number of iterations allowed in the Fisher-scoring algorithm.
#' @param precision convergence tolerance limit for the Fisher-scoring algorithm.
#' @param print_result print coefficient or not, default value is TRUE.
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
#'
#' @details
#' The model has a form that is response ~ auxiliary variables.
#' where numeric type response variables can contain NA.
#' When the response variable contains NA it will be estimated with synthetic estimator.
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
#' @md
eblup_fh <- function(
  formula,
  vardir,
  domain = NULL,
  data,
  method = c("REML", "ML"),
  maxiter = 100,
  precision = 1e-4,
  print_result = TRUE
) {
  method <- match.arg(method, choices = c("REML", "ML"))
  if (is.null(domain)) {
    domain <- 1:nrow(data)
  } else {
    domain <- .get_variable(data, domain)
  }

  # model frame & validasi
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  vardir <- .get_variable(data, vardir)

  if (nrow(mf) != length(vardir)) {
    cli::cli_abort("Length of 'vardir' must equal number of observations in data ({nrow(mf)} vs {length(vardir)}).")
  }

  y <- stats::model.response(mf, "numeric")
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
  res$model <- "FH"

  # Add domain identifier as first column in df_eblup
  res$df_eblup$domain <- domain
  res$df_eblup <- res$df_eblup[, c("domain", "y", "eblup", "vardir", "random_effect","mse", "rse")]

  res$call <- match.call()
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
