#' Empirical Best Linear Unbiased Prediction based on a Spatial Fay-Herriot Model.
#'
#' @description This function gives the Spatial Empirical Best Linear Unbiased Prediction (EBLUP) or Empirical Best (EB) predictor under normality based on a Fay-Herriot model.
#'
#' @references
#' \enumerate{
#'  \item Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley & Sons.
#' }
#'
#' @param formula an object of class formula that contains a description of the model to be fitted. The variables included in the formula must be contained in the data.
#' @param data a data frame or a data frame extension (e.g. a tibble).
#' @param vardir vector or column names from data that contain variance sampling from the direct estimator for each area.
#' @param domain vector, column name or one-sided formula referencing a domain names column
#'   in \code{data}. If NULL, the domains are numbered consecutively.
#' @param method Fitting method can be chosen between 'ML' and 'REML'.
#' @param W A square matrix with dimension equal to the TOTAL number of domains in `data`
#' (including any unsampled domains where the response is `NA`). It should contain the
#' row-standardized spatial weights (proximities) between ALL domains, with values ranging
#' from 0 to 1. Rows and columns must be ordered consistently with the domain identifiers in
#' `data`. Do NOT pre-subset `W` to sampled domains only -- unsampled domains' spatial MSE/
#' prediction relies on their relation (in `W`) to sampled neighbors.
#' @param mse_method a character string determining the estimation method of the MSE.
#'   Methods that can be chosen: "analytical", "pbmse", or "npbmse".
#' When there are unsampled domains, "pbmse"/"npbmse" bootstrap MSE is only defined for the
#' sampled domains; unsampled domains automatically get a full-spatial synthetic (kriging)
#' prediction and analytical MSE merged into the result regardless of `mse_method`.
#' @param B Number of bootstrap replications when mse_method = "pbmse" or "npbmse".
#' @param n_threads Number of threads used in parallel computation (default 1).
#'   Values less than or equal to 0 use all available cores.
#' @param seed Integer seed for bootstrap resampling.
#'   A value of -1 leaves the current R RNG state unchanged.
#' @param maxiter maximum number of iterations allowed in the Fisher-scoring algorithm. Default is 100 iterations.
#' @param precision convergence tolerance limit for the Fisher-scoring algorithm. Default value is 0.0001.
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
#' @param print_result print coefficient or not, default value is TRUE.
#'
#' @returns The function returns a list with the following objects:
#' \code{estcoef} a data frame with the estimated model coefficients in the first column (beta),
#'    their asymptotic standard errors in the second column (std.error),
#'    the t-statistics in the third column (tvalue) and the p-values of the significance of each coefficient
#'    in last column (pvalue) \cr
#' \code{random_effect_var} estimated random effect variance \cr
#' \code{rho} estimated spatial autocorrelation parameter \cr
#' \code{estvarcomp} a data frame with parameter, estimate, and std.error \cr
#' \code{goodness} vector containing several goodness-of-fit measures: loglikelihood, AIC, and BIC \cr
#' \code{df_eblup} a data frame that contains the following columns: \cr
#'    * \code{y} variable response \cr
#'    * \code{eblup} estimated results for each area \cr
#'    * \code{random_effect} random effect for each area \cr
#'    * \code{vardir} variance sampling from the direct estimator for each area \cr
#'    * \code{mse} Mean Square Error \cr
#'    * \code{rse} Relative Standart Error (%) \cr
#'    * \code{mse_pb} or \code{mse_npb} Parametric / Non Parametric bootstrap MSE\cr
#'    * \code{mse_pbbc} or \code{mse_npbbc} Bias Corrected Parametric / Non Parametric bootstrap MSE \cr
#'
#' @details
#' The model has a form that is response ~ auxiliary variables.
#' where numeric type response variables can contain NA.
#' When the response variable contains NA, that domain is treated as unsampled: it is
#' predicted with a full-spatial synthetic (kriging) estimator that borrows strength from
#' sampled neighbors via `W`. The SAR process is defined on all domains of `W`, so the
#' sampled domains have covariance \eqn{\sigma^2 \Sigma[s,s]} with
#' \eqn{\Sigma = [(I - \rho W')(I - \rho W)]^{-1}}.
#'
#' Analytical MSE. Sampled domains: \eqn{g_1 + g_2 + 2 g_3 - g_4} (Singh et al., 2005;
#' as \code{sae::mseSFH}), minus \eqn{b_{ML}'\nabla g_1} for ML. Unsampled domains:
#' \eqn{g_1 + g_2 + g_3 - \frac{1}{2}\mathrm{tr}(I^{-1}\nabla^2 g_1)} (minus
#' \eqn{b_{ML}'\nabla g_1} for ML), with \eqn{g_3} from the exact derivative of the
#' kriging predictor; for sampled domains this second-order form agrees, up to o(1/m), with
#' \eqn{g_1 + g_2 + 2 g_3 - g_4}. The inverse information used for \eqn{g_3} and
#' \eqn{b_{ML}} is the one of the fitting method: REML information for REML and ML
#' information for ML. (\code{sae::mseSFH} uses the REML information also for ML; the two
#' agree asymptotically, but with few domains, or \eqn{\hat\sigma^2 \approx 0} with
#' \eqn{\hat\rho \approx 1}, the REML choice can make the ML bias term explode. Hence, for
#' \code{method = "ML"}, the MSE differs from \code{sae::mseSFH}; REML results agree.)
#'
#' When \eqn{\hat\rho} lies on the boundary (\eqn{|\hat\rho| \ge 0.999}, or reported as
#' \eqn{-1}), the Prasad-Rao expansion in \eqn{\rho} is not valid and can explode, so
#' \eqn{\rho} is treated as fixed: only the \eqn{\sigma^2} block of the information is
#' inverted, which removes the \eqn{\rho} terms from \eqn{g_3}, \eqn{g_4} and
#' \eqn{b_{ML}} (the standard error of \eqn{\rho} is then \code{NA}). In this case the
#' MSE differs from \code{sae::mseSFH}, which keeps those terms.
#'
#' If the classical (sae-compatible) Fisher scoring fails, or converges to a negative
#' \eqn{\sigma^2} (which \code{sae::eblupSFH} truncates to 0 without re-estimating
#' \eqn{\rho}), the (RE)ML likelihood is maximised under the constraints
#' \eqn{\sigma^2 \ge 0}, \eqn{|\rho| \le 0.999}. A lower-bound estimate
#' \eqn{\rho = -0.999} is reported as \eqn{-1} (as \code{sae}) only when \eqn{I + W} is
#' non-singular (not for bipartite graphs such as a rook grid).
#'
#' Computation. No \eqn{D \times D} inverse is formed: with
#' \eqn{A = (I - \rho W')(I - \rho W)} and \eqn{M = A + \sigma^2 E'\Psi^{-1}E}
#' (\eqn{E} selects the sampled domains), \eqn{V^{-1} = \Psi^{-1} - \sigma^2\Psi^{-1}EM^{-1}E'\Psi^{-1}},
#' \eqn{g_1 = \sigma^2\,\mathrm{diag}(M^{-1})} and the remaining terms follow from Cholesky
#' solves with \eqn{A} and \eqn{M}. These are sparse for a sparse \code{W} (contiguity or
#' k-nearest-neighbour weights), which makes models with a few thousand domains feasible;
#' a dense \code{W} (e.g. inverse-distance weights) uses dense factorisations.
#'
#' @export
#' @examples
#' library(fastsae)
#'
#' # Spatial Fay-Herriot model
#' m1 <- eblup_sfh(
#'   y ~ x1 + x2 + x3,
#'   data = mys,
#'   vardir = ~vardir,
#'   W = mys_proxmat
#' )
#'
#' # Spatial Fay-Herriot model with Parametric Bootstrap MSE
#' m2 <- eblup_sfh(
#'   y ~ x1 + x2 + x3,
#'   data = mys,
#'   vardir = ~vardir,
#'   mse_method = "pbmse",
#'   B = 50,
#'   W = mys_proxmat
#' )
#'
#' @md
eblup_sfh <- function(
  formula,
  vardir,
  domain = NULL,
  data,
  method = c("REML", "ML"),
  mse_method = c("analytical", "pbmse", "npbmse"),
  W = NULL,
  B = 100,
  n_threads = 1,
  seed = -1,
  maxiter = 100,
  precision = 1e-4,
  self_benchmark = FALSE,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL,
  print_result = TRUE
) {
  method <- match.arg(toupper(method), choices = c("REML", "ML"))
  mse_method <- match.arg(tolower(mse_method), choices = c("analytical", "pbmse", "npbmse"))
  if (is.null(domain)) {
    domain <- 1:nrow(data)
  } else {
    domain <- .get_variable(data, domain)
  }

  # model frame & validasi
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  vardir <- .get_variable(data, vardir)


  if (nrow(mf) != length(vardir)) {
    cli::cli_abort("Length of 'vardir' must equal number of observations in data.")
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

  # Validasi vardir untuk area tersampel saja (y tidak NA)
  y_valid <- !is.na(y)
  if (any(vardir[y_valid] <= 0, na.rm = TRUE)) {
    cli::cli_abort("vardir must be strictly positive for sampled areas.")
  }

  if (is.null(W)) {
    if (inherits(data, "sf")) {
      cli::cli_alert_info("Argument {.arg W} not provided; automatically constructing Queen spatial weights from {.cls sf} polygon geometry.")
      W <- create_weights(data, method = "queen", style = "W")
    } else {
      cli::cli_abort("Argument `W` (spatial weight matrix) must be provided.")
    }
  }

  n_total <- nrow(X)
  if (!is.matrix(W) || nrow(W) != n_total || ncol(W) != n_total) {
    cli::cli_abort(c(
      "`W` must be a square matrix with dimension equal to the total number of domains in `data` ({n_total}), including any unsampled ones.",
      "i" = "Got a {nrow(W)}x{ncol(W)} matrix. Do not pre-subset `W` to sampled domains only."
    ))
  }


  idx_s <- !is.na(y)
  has_unsampled <- any(!idx_s)

  if (mse_method %in% c("pbmse", "npbmse")) {
    if (has_unsampled) {
      cli::cli_alert_info(
        "{sum(!idx_s)} unsampled domain(s) detected. Bootstrap MSE ({mse_method}) is computed for the {sum(idx_s)} sampled domain(s); unsampled domain(s) get a full-spatial synthetic (kriging) prediction and analytical MSE instead."
      )
    }

    # Data penuh (termasuk area tak-tersampel, y = NA) dan W penuh: proses SAR
    # didefinisikan pada SEMUA domain, sehingga fit, replikasi bootstrap, dan
    # prediksi kriging area tak-tersampel memakai model yang sama. Kolom
    # bootstrap (mse_pb/mse_pbbc atau mse_npb/mse_npbbc) bernilai NA untuk
    # area tak-tersampel; mse/rse mereka berasal dari kriging analitis.
    boot_fun <- if (mse_method == "pbmse") .seblup_pbmse else .seblup_npbmse
    res <- boot_fun(
      Xall = X, yall = y, vardirall = vardir, method = method, W = W,
      maxiter = maxiter, precision = precision, B = B,
      n_threads = n_threads, seed = seed
    )
  } else {
    res <- .seblup_core(
      Xall = X,
      yall = y,
      vardirall = vardir,
      method = method,
      Wall = W,
      maxiter = maxiter,
      precision = precision
    )
  }

  # attach metadata
  row.names(res$estcoef) <- colnames(X)
  res$formula <- formula
  res$model <- "SFH"

  # Tambahkan domain identifier ke df_eblup sebagai kolom pertama
  res$df_eblup$domain <- domain
  res$df_eblup <- res$df_eblup[, c("domain", setdiff(names(res$df_eblup), "domain"))]

  if (inherits(data, "sf") && requireNamespace("sf", quietly = TRUE)) {
    res$df_eblup <- sf::st_set_geometry(res$df_eblup, sf::st_geometry(data))
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
  res$W <- W
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
