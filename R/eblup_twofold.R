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
#' @param vardir sampling variances of the sub-area direct estimates: a numeric vector, the name
#'   of a column of \code{data} or a one-sided formula.
#' @param domain vector, column name or one-sided formula referencing the area identifier column
#'   in \code{data}. Must be supplied.
#' @param subarea vector, column name or one-sided formula referencing the subarea identifier column
#'   in \code{data}. If NULL, rows are numbered consecutively.
#' @param method Fitting method can be chosen between 'REML' and 'ML'.
#' @param mse_type MSE estimator: 'analytical' (default; Prasad-Rao g1+g2+2*g3, fast; ML adds
#'   the bias correction), 'bootstrap' (parametric bootstrap, slower but robust) or 'none'
#'   (no MSE; \code{mse} and \code{rse} are \code{NA}).
#' @param B number of parametric bootstrap replicates (only used if mse_type = 'bootstrap').
#' @param seed integer seed for the bootstrap; \code{-1} (default) or \code{NULL} sets no seed.
#'   A fixed seed does not change the random number stream of the session.
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
#' @param weight Optional vector, column name or one-sided formula referencing the sub-area
#'   weights \eqn{w_{ij}} (e.g. population sizes) used for the area means
#'   \eqn{\theta_i = \sum_j w_{ij}\theta_{ij}}. They are normalised within each area; if
#'   \code{NULL}, every sub-area of an area gets the same weight. Missing or negative weights
#'   get weight 0 (with a warning), as in \code{\link{hb_twofold}}.
#' @param compute_area Logical. If \code{TRUE} (default), adds the EBLUPs of the area means and
#'   their MSE in \code{df_area} (see Details).
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
#' \code{df_area} (if \code{compute_area = TRUE}) a data frame with one row per area: domain,
#' \code{eblup_area}, \code{mse_area}, \code{rse_area}, \code{n_subareas} (sub-areas in
#' \code{data}) and \code{n_sampled} (those with a direct estimate). \cr
#'
#' @details
#' The model has a form that is response ~ auxiliary variables, with one row per subarea.
#' When the response contains NA the subarea is treated as non-sampled. It is predicted by
#' \eqn{x_{ij}'\hat\beta + \hat v_i}, borrowing the area effect from the sampled sub-areas of
#' its area; only the sub-areas of an area without any sampled sub-area get the synthetic
#' estimator \eqn{x_{ij}'\hat\beta}. MSE is estimated either analytically (Prasad-Rao
#' g1+g2+2*g3 for REML, minus the ML bias term b'grad(g1) for ML, with an exact matrix
#' derivation for g3) or by parametric bootstrap.
#'
#' Area means. The target of area \eqn{i} is \eqn{\theta_i = \sum_j w_{ij}\theta_{ij}} over
#' the sub-areas of the area in \code{data} (sampled or not), with the normalised
#' \code{weight}. Its EBLUP is \eqn{\sum_j w_{ij}\tilde\theta_{ij}}. The analytical MSE is the
#' same second-order estimator as for the sub-areas, applied to this linear target (Torabi
#' and Rao, 2014, derive it for population-weighted area means): it includes the covariance
#' of the sub-area prediction errors through the shared area effect and \eqn{\hat\beta}, so
#' it is not \eqn{\sum_j w_{ij}^2\,\mathrm{mse}_{ij}}. An area without sampled sub-areas gets
#' the synthetic estimator \eqn{\bar x_i'\hat\beta} with MSE
#' \eqn{\bar x_i' Q \bar x_i + \hat\sigma_v^2 + \hat\sigma_u^2\sum_j w_{ij}^2}. With
#' \code{mse_type = "bootstrap"} the area MSE comes from the same bootstrap replicates.
#' Sub-areas of the population that are not rows of \code{data} are not part of the area mean.
#'
#' @seealso \code{\link{hb_twofold}} for the hierarchical Bayes counterpart, and
#'   \code{vignette("twofold-subarea", package = "fastsae")}.
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
#' m1 <- eblup_twofold(y ~ x1,
#'   vardir = "vardir", domain = "area",
#'   subarea = "subarea", data = dat
#' )
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
  B = 100,
  seed = -1,
  maxiter = 100,
  precision = 1e-4,
  self_benchmark = FALSE,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL,
  weight = NULL,
  compute_area = TRUE,
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

  .check_covariates_na(X, "subarea")
  .check_area_inputs(X, y, vardir, unit = "subarea")

  area_idx <- as.integer(factor(domain)) - 1L
  unique_areas <- sort(unique(domain))

  norm_weight <- NULL
  if (isTRUE(compute_area)) {
    weight_vec <- if (!is.null(weight)) .get_variable(data, weight) else NULL
    norm_weight <- .twofold_weights(weight_vec, domain, unique_areas, subarea)
  }
  if (length(unique(area_idx[!is.na(y)])) < 2) {
    cli::cli_abort("At least 2 areas with sampled subareas are needed to estimate the area variance.")
  }

  res <- .with_seed(seed, .eblup_twofold_core(
    Xall = X,
    yall = y,
    vardirall = vardir,
    area = area_idx,
    method = method,
    mse_type = if (mse_type == "none") "analytical" else mse_type,
    B = as.integer(B),
    maxiter = maxiter,
    precision = precision,
    wall_ = norm_weight
  ))
  area_eblup <- res$area_eblup
  area_mse <- res$area_mse
  res$area_eblup <- res$area_mse <- NULL

  # attach metadata
  row.names(res$estcoef) <- colnames(X)
  res$formula <- formula
  res$model <- "TWOFOLD"
  res$method <- method
  res$mse_type <- mse_type

  # Add domain/subarea identifiers as first columns in df_eblup
  res$df_eblup$domain <- domain
  res$df_eblup$subarea <- subarea
  res$df_eblup <- res$df_eblup[, c(
    "domain", "subarea", "y", "eblup", "vardir",
    "random_effect_area", "random_effect_subarea",
    "mse", "rse"
  )]

  # Area means (see Details); a negative area MSE is NA, as for the sub-areas below
  if (isTRUE(compute_area)) {
    sampled <- !is.na(y)
    if (mse_type == "none") area_mse[] <- NA_real_
    area_mse[which(area_mse < 0)] <- NA_real_
    res$df_area <- data.frame(
      domain = unique_areas,
      eblup_area = area_eblup,
      mse_area = area_mse,
      rse_area = ifelse(abs(area_eblup) < .Machine$double.eps, NA_real_, sqrt(area_mse) / abs(area_eblup) * 100),
      n_subareas = as.integer(tabulate(area_idx + 1L, length(unique_areas))),
      n_sampled = as.integer(tabulate(area_idx[sampled] + 1L, length(unique_areas))),
      stringsAsFactors = FALSE
    )
  }

  # A negative second-order MSE estimate (possible with the ML bias correction) is not a
  # usable MSE: report NA rather than a NaN RSE
  if (mse_type == "none") res$df_eblup <- .drop_mse(res$df_eblup)
  neg_mse <- which(res$df_eblup$mse < 0)
  if (length(neg_mse) > 0) {
    cli::cli_warn(c(
      "The MSE estimate is negative for {length(neg_mse)} subarea{?s}; set to {.code NA}.",
      "i" = "Consider {.code method = \"REML\"} or {.code mse_type = \"bootstrap\"}."
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
