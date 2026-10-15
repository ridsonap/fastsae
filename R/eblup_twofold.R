#' Empirical Best Linear Unbiased Prediction under a Two-fold Fay-Herriot Model.
#'
#' @description This function gives the Empirical Best Linear Unbiased Prediction (EBLUP)
#'   under a two-fold Fay-Herriot model (Torabi & Rao, 2014): subareas nested within areas
#'   with an area random effect \eqn{v_i} and a subarea random effect \eqn{u_{ij}}.
#'
#' @references
#' \enumerate{
#'  \item Torabi, M., & Rao, J. N. K. (2014). On small area estimation under a sub-area
#'    level model. \emph{Journal of Multivariate Analysis}, 127, 36--55.
#'  \item Rao, J. N. K., & Molina, I. (2015). \emph{Small Area Estimation}. John Wiley & Sons.
#' }
#'
#' @param formula an object of class formula that contains a description of the model to be fitted.
#' @param data a data frame or a data frame extension (e.g. a tibble), one row per subarea.
#' @param vardir vector or column names from data that contain variance sampling from the direct estimator.
#' @param domain vector, column name or one-sided formula referencing the area identifier column
#'   in \code{data}. Must be supplied.
#' @param subarea vector, column name or one-sided formula referencing the subarea identifier column
#'   in \code{data}. If NULL, rows are numbered consecutively.
#' @param method Fitting method can be chosen between 'REML' and 'ML'.
#' @param mse_type MSE estimator: 'analytical' (default; Prasad-Rao g1+g2+2*g3, fast; ML adds
#'   the bias correction), 'bootstrap' (parametric bootstrap, slower but robust) or 'none'
#'   (no MSE; \code{mse} and \code{rse} are \code{NA}).
#' @param B number of parametric bootstrap replicates (only used if mse_type = 'bootstrap').
#' @param seed integer seed for bootstrap resampling (NULL = no seed set).
#' @param maxiter maximum number of iterations allowed in the Fisher-scoring algorithm.
#' @param precision convergence tolerance limit for the Fisher-scoring algorithm.
#' @param self_benchmark Logical. If \code{TRUE}, fits a self-benchmarking model via the
#'   augmented model approach (Wang, Fuller, and Qu, 2008; Rao and Molina, 2015, Sec 10.5.1),
#'   ensuring that the weighted aggregate of small area estimates automatically matches the direct
#'   survey aggregate or specified benchmark targets without post-hoc adjustment.
#' @param benchmark_weight Optional vector, column name, or one-sided formula referencing the domain
#'   benchmarking weights in \code{data} (e.g. population sizes or shares). Defaults to equal weights.
#' @param benchmark_target Optional numeric scalar or named vector specifying fixed aggregate
#'   benchmark target(s). If \code{NULL} (default), calibrates to the weighted aggregate of the direct estimates.
#' @param benchmark_group Optional vector, column name, or one-sided formula referencing the
#'   grouping / stratum column in \code{data} (e.g. area identifier for subarea-to-area benchmarking)
#'   for group-specific self-benchmarking.
#' @param print_result print coefficient or not, default value is TRUE.
#'
#' @returns The function returns a list with the following objects:
#' \code{estcoef} a data frame with the estimated model coefficients (\code{beta}), their
#'    standard errors (\code{std.error}, repeated as \code{stderr_beta} for compatibility
#'    with fastsae 0.1.0), z-statistics (\code{zvalue}) and normal p-values (\code{pvalue}),
#' \code{random_effect_var} named vector with the estimated area (\code{sigma2_v}) and
#' subarea (\code{sigma2_u}) variances (paper notation),
#' \code{goodness} named vector: \code{loglikelihood} (the full Gaussian log-likelihood at the
#'    estimates, for ML and REML fits alike), \code{AIC} and \code{BIC} computed from it,
#'    following the convention of \pkg{sae} for the Fay-Herriot models; BIC uses the number
#'    of sampled subareas. The restricted (REML) likelihood is not reported: it depends on the
#'    covariates through X, so it cannot compare models with different covariates.,
#' \code{df_eblup} a data frame that contains domain, subarea, y, eblup, vardir,
#' random_effect_area, random_effect_subarea, mse, and rse. \cr
#'
#' @details
#' The model has a form that is response ~ auxiliary variables, with one row per subarea.
#' When the response contains NA the subarea is treated as non-sampled and predicted
#' with the synthetic estimator. MSE is estimated either analytically (Prasad-Rao
#' g1+g2+2*g3 for REML, minus the ML bias term b'grad(g1) for ML, with an exact matrix
#' derivation for g3) or by parametric bootstrap.
#'
#' @export
#' @examples
#' library(fastsae)
#'
#' # Two-fold Fay-Herriot model (simulated illustration)
#' set.seed(1)
#' m <- 20
#' dat <- do.call(rbind, lapply(1:m, function(d) {
#'   nd <- sample(2:5, 1)
#'   data.frame(
#'     area = d, subarea = paste0(d, "-", seq_len(nd)),
#'     x1 = rnorm(nd), y = NA_real_, vardir = runif(nd, 0.2, 1)
#'   )
#' }))
#' dat$y <- 1 + dat$x1 + rnorm(m)[dat$area] + rnorm(nrow(dat), sd = 0.5) +
#'   rnorm(nrow(dat), sd = sqrt(dat$vardir))
#' m1 <- eblup_twofold(y ~ x1, vardir = "vardir", domain = "area",
#'                     subarea = "subarea", data = dat)
#'
#' @md
eblup_twofold <- function(
  formula,
  vardir,
  domain = NULL,
  subarea = NULL,
  data,
  method = c("REML", "ML"),
  mse_type = c("analytical", "bootstrap", "none"),
  B = 200,
  seed = NULL,
  maxiter = 100,
  precision = 1e-4,
  self_benchmark = FALSE,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL,
  print_result = TRUE
) {
  method <- match.arg(method, choices = c("REML", "ML"))
  mse_type <- match.arg(mse_type, choices = c("analytical", "bootstrap", "none"))
  if (is.null(domain)) {
    cli::cli_abort("`domain` (area identifier) must be supplied for the two-fold model.")
  }
  domain <- .get_variable(data, domain)
  if (is.null(subarea)) {
    subarea <- seq_len(nrow(data))
  } else {
    subarea <- .get_variable(data, subarea)
  }

  # model frame & validasi
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  vardir <- .get_variable(data, vardir)

  if (nrow(mf) != length(vardir)) {
    cli::cli_abort("Length of 'vardir' must equal number of observations in data ({nrow(mf)} vs {length(vardir)}).")
  }

  y <- stats::model.response(mf, "numeric")

  sb_info <- NULL
  if (isTRUE(self_benchmark)) {
    sb_info <- .setup_self_benchmark(
      data = data,
      formula = formula,
      y = y,
      vardir = vardir,
      weight = benchmark_weight,
      target = benchmark_target,
      group = benchmark_group
    )
    data <- sb_info$data_aug
    formula <- sb_info$formula_aug
    mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
    y <- stats::model.response(mf, "numeric")
  }

  X <- stats::model.matrix(attr(mf, "terms"), mf)

  # Check auxiliary variables for NA values
  if (anyNA(X)) {
    cli::cli_abort("Auxiliary variables contain NA values.")
  }
  .check_area_inputs(X, y, vardir, unit = "subarea")

  area_idx <- as.integer(factor(domain)) - 1L
  if (length(unique(area_idx[!is.na(y)])) < 2) {
    cli::cli_abort("At least 2 areas with sampled subareas are needed to estimate the area variance.")
  }

  if (!is.null(seed)) set.seed(seed)

  res <- .eblup_twofold_core(
    Xall = X,
    yall = y,
    vardirall = vardir,
    area = area_idx,
    method = method,
    mse_type = if (mse_type == "none") "analytical" else mse_type,
    B = as.integer(B),
    maxiter = maxiter,
    precision = precision
  )

  # attach metadata
  row.names(res$estcoef) <- colnames(X)
  res$formula <- formula
  res$model <- "TWOFOLD"
  res$method <- method
  res$mse_type <- mse_type

  # Add domain/subarea identifiers as first columns in df_eblup
  res$df_eblup$domain <- domain
  res$df_eblup$subarea <- subarea
  res$df_eblup <- res$df_eblup[, c("domain", "subarea", "y", "eblup", "vardir",
                                   "random_effect_area", "random_effect_subarea",
                                   "mse", "rse")]

  # A negative second-order MSE estimate (possible with the ML bias correction) is not a
  # usable MSE: report NA rather than a NaN RSE
  if (mse_type == "none") res$df_eblup <- .drop_mse(res$df_eblup)
  neg_mse <- which(res$df_eblup$mse < 0)
  if (length(neg_mse) > 0) {
    cli::cli_warn(c(
      "The MSE estimate is negative for {length(neg_mse)} subarea{?s}; set to {.code NA}.",
      "i" = "Consider {.code method = \"REML\"} or {.code mse = \"bootstrap\"}."
    ))
    res$df_eblup$mse[neg_mse] <- NA_real_
    res$df_eblup$rse[neg_mse] <- NA_real_
  }

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
    .warn_no_convergence(res$n_iter)
    return(res)
  }

  if (print_result) {
    print(res)
  }
  return(res)
}
