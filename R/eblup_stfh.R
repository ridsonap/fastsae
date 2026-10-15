#' Empirical Best Linear Unbiased Prediction based on a Spatio-Temporal Fay-Herriot Model.
#'
#' @description This function gives the Spatio-Temporal Empirical Best Linear
#' Unbiased Prediction (EBLUP) under normality based on a spatio-temporal
#' Fay-Herriot model. It reimplements the same Fisher-scoring algorithm as
#' \code{eblupSTFH()} (package \pkg{sae}, Marhuenda, Molina & Morales 2013),
#' but the estimation loop runs in compiled C++/Armadillo, making it substantially faster and more memory-efficient than the original
#' R implementation -- especially for a large number of domains/time periods.
#'
#' @references
#' \enumerate{
#'  \item Marhuenda, Y., Molina, I., & Morales, D. (2013). Small area estimation
#'    with spatio-temporal Fay-Herriot models. Computational Statistics & Data
#'    Analysis, 58, 308-325.
#'  \item Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley & Sons.
#'  \item Saei, A., & Chambers, R. (2005). Out of sample estimation for small areas using
#'    area level data. Southampton Statistical Sciences Research Institute, Methodology
#'    Working Paper M05/11.
#'  \item Benedetti, R., Piersimoni, F., Pratesi, M., Salvati, N., & Suesse, T. (2024).
#'    Handling out-of-sample areas to estimate the unemployment rate at local labour market
#'    areas in Italy. International Statistical Review. \doi{10.1111/insr.12596}
#' }
#'
#' @param formula an object of class formula describing the model to fit
#'   (response ~ auxiliary variables). Variables must be present in \code{data}.
#' @param data a data frame (or extension) with one row per domain and time period
#'   (\code{domain * time} rows), including the cells without a direct estimate: for
#'   those, the response (and \code{vardir}) is \code{NA}, while the auxiliary
#'   variables must be available (see Details). Rows may be in any order: they are
#'   sorted internally into the domain-major order required by the model (domains in
#'   order of first appearance, periods increasing), and the results are returned in
#'   the original row order. Note that \code{sae::eblupSTFH()} instead assumes a
#'   complete panel already in domain-major order.
#' @param vardir vector, column name or one-sided formula referencing a column
#'   in \code{data}, with the sampling variances of the direct estimator. Must be
#'   positive in the sampled cells; it is ignored (may be \code{NA}) where the
#'   response is \code{NA}.
#' @param domain vector, column name or one-sided formula referencing a domain names column
#'   in \code{data}. If NULL, the domains are numbered consecutively.
#' @param time vector, column name, or one-sided formula referencing a time names column
#'   in \code{data}.
#' @param W a square proximity/spatial weights matrix of dimension
#'   \code{domain x domain} (row-standardized, values typically in \eqn{[0,1]}), with rows
#'   and columns in the order in which the domains first appear in \code{data}. It must
#'   cover all domains, including the ones that are never sampled: do not drop them from
#'   \code{W}, since their predictions borrow strength from their neighbours through it.
#' @param model character, either \code{"ST"} (spatio-temporal, default) or
#'   \code{"S"} (spatial only, no AR(1) temporal component).
#' @param maxiter maximum number of Fisher-scoring iterations. Default 100.
#' @param precision convergence tolerance for the Fisher-scoring algorithm.
#'   Default \code{1e-4}.
#' @param sigma21_start,rho1_start,sigma22_start,rho2_start starting values for
#'   the variance/autocorrelation components. Defaults mirror \code{eblupSTFH()}:
#'   \code{0.5 * median(vardir)} for the variances and \code{0.5} for the
#'   autocorrelations. \code{rho2_start} is ignored when \code{model = "S"}.
#'   The REML fit first uses the same Fisher-scoring iteration as \code{sae::eblupSTFH()};
#'   if that fails (typically a variance iterate turns negative and the covariance matrix
#'   becomes indefinite) it falls back to constrained REML (variances \eqn{\ge 0},
#'   \eqn{|\rho| \le 0.999}) with step halving. When a variance is estimated at 0, its
#'   autocorrelation is not identified and its standard error is reported as \code{NA}.
#' @param compute_mse logical, if \code{TRUE} computes MSE. Default \code{FALSE}.
#' @param mse_type character, method for MSE estimation: \code{"analytical"} (fast
#'   Prasad-Rao type estimator g1 + g2 + 2 g3 - g4, default when \code{compute_mse = TRUE}; g3 uses
#'   the exact derivative of the EBLUP, \eqn{\Psi P V_a P y}, and the REML information at the
#'   final estimates; correlations estimated on the boundary \eqn{|\rho| \ge 0.999} are
#'   treated as fixed, so they do not enter \eqn{g_3}/\eqn{g_4}) or \code{"bootstrap"}
#'   (parametric bootstrap in which the variance components are re-estimated by REML in every
#'   replicate, as in Marhuenda et al., 2013).
#' @param B number of bootstrap replicates for MSE computation. Only used when
#'   \code{compute_mse = TRUE} and \code{mse_type = "bootstrap"}. Default \code{100}.
#' @param n_threads number of threads for parallel bootstrap MSE computation.
#'   Use \code{0} for all available cores. Default \code{1}.
#' @param seed random seed for bootstrap. Use \code{-1} for no seed. Default \code{-1}.
#' @param self_benchmark Logical. If \code{TRUE}, fits a self-benchmarking model via the
#'   augmented model approach (Wang, Fuller, and Qu, 2008; Rao and Molina, 2015, Sec 10.5.1),
#'   ensuring that the weighted aggregate of small area estimates automatically matches the direct
#'   survey aggregate or specified benchmark targets without post-hoc adjustment.
#' @param benchmark_weight Optional vector, column name, or one-sided formula referencing the domain
#'   benchmarking weights in \code{data} (e.g. population sizes or shares). Defaults to equal weights.
#' @param benchmark_target Optional numeric scalar or named vector specifying fixed aggregate
#'   benchmark target(s). If \code{NULL} (default), calibrates to the weighted aggregate of the direct estimates.
#' @param benchmark_group Optional vector, column name, or one-sided formula referencing the
#'   grouping / stratum column in \code{data} (e.g. time period / year) for group-specific self-benchmarking.
#' @param print_result print the estimated coefficients or not. Default \code{TRUE}.
#'
#' @returns A list with the same structure as \code{seblup_area()}, with additional
#'   \code{estvarcomp} for spatio-temporal variance/autocorrelation components:
#'   \describe{
#'     \item{\code{estcoef}}{data frame with \code{beta}, \code{std.error}, \code{tvalue},
#'       \code{pvalue}, as in \code{sae::eblupSTFH()}, plus \code{stderr_beta} and \code{zvalue}
#'       (copies, for the column names of the other fastsae models). The p-values use the
#'       normal distribution, so \code{tvalue} is a z-statistic.}
#'     \item{\code{estvarcomp}}{data frame with estimate, std.error for sigma21, rho1, sigma22, rho2.}
#'     \item{\code{goodness}}{named vector: \code{loglikelihood}, the full Gaussian
#'       log-likelihood at the REML estimates (as \code{loglike} in \code{sae::eblupSTFH()}),
#'       \code{AIC} and \code{BIC} computed from it, and \code{loglikelihood_reml}, the
#'       restricted log-likelihood that REML maximises. The full log-likelihood is used for AIC
#'       and BIC because the restricted one depends on the covariates.}
#'     \item{\code{df_eblup}}{data frame with domain, time, y, eblup, vardir, mse, rse,
#'       random_effect_u1, random_effect_u2 and \code{sampled} (\code{FALSE} for the cells
#'       predicted without a direct estimate).}
#'       When \code{compute_mse = FALSE}, mse and rse are \code{NA}.
#'     \item{\code{model}}{model type ("ST" or "S").}
#'     \item{\code{convergence}}{logical, whether the algorithm converged.}
#'     \item{\code{n_iter}}{number of iterations.}
#'     \item{\code{B}}{number of bootstrap replicates (only when compute_mse = TRUE and mse_type = "bootstrap").}
#'   }
#'
#' @details
#' \strong{Estimation method.} The variance components are estimated by REML only, as in
#' Marhuenda, Molina and Morales (2013) and \code{sae::eblupSTFH()}, which derive the
#' Fisher-scoring algorithm, the analytical MSE and the bootstrap for REML; there is no
#' \code{method} argument. The standard errors of the variance components are those of
#' \code{sae::eblupSTFH()}: the inverse REML information of the last Fisher-scoring
#' iteration.
#'
#' \strong{Unsampled domains and missing periods.} Cells whose response is \code{NA}
#' (a domain without a sample in some period, or in every period) are not dropped. The
#' model is defined on the full \eqn{D \times T} grid: the SAR effect \eqn{u_1} lives on
#' all \eqn{D} domains of \code{W} and the AR(1) effect \eqn{u_2} on all \eqn{T} periods.
#' The variance components are estimated by REML from the sampled cells, i.e. from the
#' marginal likelihood with the missing cells integrated out (missing at random), and every
#' cell gets the best linear unbiased predictor of \eqn{x_{dt}'\beta + u_{1d} + u_{2dt}}
#' given the sampled data (Saei and Chambers, 2005; Benedetti et al., 2024 for the
#' spatial effect):
#' \deqn{\tilde\theta_{dt} = x_{dt}'\tilde\beta + a_{dt}' V^{-1} (y - X\tilde\beta),}
#' with \eqn{a_{dt}} the covariance between the cell and the sampled direct estimates.
#' A missing period of a sampled domain borrows from the other periods of the domain (AR(1))
#' and from the neighbours (SAR); a domain that is never sampled is predicted from its
#' neighbours only (spatial kriging). With \code{model = "S"} the period effect of a missing
#' cell is independent of the data, so only the spatial part is borrowed.
#'
#' \strong{MSE.} For sampled cells the analytical MSE is \eqn{g_1 + g_2 + 2 g_3 - g_4}
#' (Marhuenda et al., 2013). For the cells predicted without a direct estimate it is
#' \eqn{g_1 + g_2 + g_3}, where \eqn{g_1} and \eqn{g_2} are the MSE of the predictor with known
#' variance components and \eqn{g_3} comes from the exact derivative of the predictor in those
#' components (\code{g4} is 0 for these cells). The further second-order correction of the
#' plug-in \eqn{g_1} is not applied: it is numerically unstable when \eqn{\rho_1} approaches
#' \eqn{\pm 1} or a variance component is weakly identified, while in simulations
#' \eqn{g_1 + g_2 + g_3} has a relative bias of a few percent and close to nominal coverage
#' for these cells. The parametric bootstrap generates the
#' random effects on the full grid and the sampling errors on the sampled cells, re-estimates
#' the model on the sampled cells and predicts every cell, so it covers both kinds of cells.
#' The auxiliary variables must be known for every cell, since they enter the prediction.
#' Self-benchmarking calibrates the weighted aggregate of the sampled cells of each group.
#'
#' @export
#' @examples
#' library(fastsae)
#'
#' # mys_panel: 42 domains x 11 years; domains 21 and 25 are never sampled (y = NA).
#' # They stay in the data and in W, and are predicted from their neighbours.
#' m1 <- eblup_stfh(
#'   y ~ x1 + x2 + x3,
#'   data = mys_panel,
#'   vardir = ~vardir,
#'   domain = ~area,
#'   time = ~year,
#'   W = mys_proxmat,
#'   model = "ST"
#' )
#'
#' # Analytical MSE for every cell, sampled or not
#' m2 <- eblup_stfh(
#'   y ~ x1 + x2 + x3,
#'   data = mys_panel,
#'   vardir = ~vardir,
#'   domain = ~area,
#'   time = ~year,
#'   W = mys_proxmat,
#'   model = "ST",
#'   compute_mse = TRUE,
#'   mse_type = "analytical"
#' )
#' head(m2$df_eblup[!m2$df_eblup$sampled, ])
#'
#' @md
eblup_stfh <- function(
  formula,
  vardir,
  data,
  domain,
  time,
  W,
  model = c("ST", "S"),
  maxiter = 100,
  precision = 1e-4,
  sigma21_start = NULL,
  rho1_start = 0.5,
  sigma22_start = NULL,
  rho2_start = 0.5,
  compute_mse = FALSE,
  mse_type = c("analytical", "bootstrap"),
  B = 100,
  n_threads = 1,
  seed = -1,
  self_benchmark = FALSE,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL,
  print_result = TRUE
) {
  model <- match.arg(model, choices = c("ST", "S"))
  if (missing(mse_type)) {
    if (!missing(B)) {
      mse_type <- "bootstrap"
    } else {
      mse_type <- "analytical"
    }
  } else {
    mse_type <- match.arg(mse_type, choices = c("analytical", "bootstrap"))
  }

  # ---- model frame & matrix design ----
  # Rows with NA in the response are cells without a direct estimate: they are
  # kept and predicted (see Details); only the covariates must be complete.

  # ---- domain-major order ----
  # The C++ core requires all periods of domain 1, then domain 2, ... with the
  # same (increasing) time order within every domain. Sort internally (domains
  # in order of first appearance, which must also be the row order of W) and
  # return the results in the caller's row order. Previously unsorted input,
  # e.g. time-major panels such as `mys_panel`, was silently mis-specified.
  data_orig <- data
  .resolve_rows <- function(v) {
    if (is.null(v)) return(NULL)
    if ((is.character(v) && length(v) == 1) || inherits(v, "formula")) {
      return(.get_variable(data, v))
    }
    v
  }
  dom0 <- .get_variable(data, domain)
  tim0 <- .get_variable(data, time)
  vardir <- .get_variable(data, vardir)
  n_rows <- nrow(data)
  ord <- if (length(dom0) == n_rows && length(tim0) == n_rows) {
    order(match(as.character(dom0), unique(as.character(dom0))), tim0)
  } else {
    seq_len(n_rows) # malformed input: leave it to the validation below
  }
  data <- data[ord, , drop = FALSE]
  domain <- if (length(dom0) == n_rows) dom0[ord] else dom0
  time <- if (length(tim0) == n_rows) tim0[ord] else tim0
  if (length(vardir) == n_rows) vardir <- vardir[ord]
  if (!is.null(benchmark_weight)) {
    bw <- .resolve_rows(benchmark_weight)
    if (length(bw) == length(ord)) benchmark_weight <- bw[ord]
  }
  if (!is.null(benchmark_group)) {
    bg <- .resolve_rows(benchmark_group)
    if (length(bg) == length(ord)) benchmark_group <- bg[ord]
  }

  # All cells (NA response = unsampled cell)
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
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
  if (anyNA(X)) {
    cli::cli_abort(c(
      "The auxiliary variables have missing values in {sum(!stats::complete.cases(X))} row{?s}.",
      "i" = "They are needed in every domain x time cell, including the unsampled ones (response {.code NA}), to predict them."
    ))
  }
  sampled <- !is.na(y)
  if (!any(sampled)) cli::cli_abort("The response is {.code NA} in every cell: there is nothing to fit.")
  if (nrow(mf) != length(vardir)) {
    cli::cli_abort("Length of 'vardir' must equal number of observations in data ({nrow(mf)} vs {length(vardir)}).")
  }

  if (anyNA(vardir[sampled]) || any(vardir[sampled] <= 0)) {
    cli::cli_abort("vardir must be available and strictly positive in every sampled cell.")
  }

  # domain & time
  n_domain <- length(unique(domain))
  n_time <- length(unique(time))

  # After sorting, every domain must have the same sequence of periods
  # (balanced panel, one row per domain x time).
  if (length(domain) == n_domain * n_time) {
    dom_blocks <- rep(unique(domain), each = n_time)
    time_blocks <- rep(time[seq_len(n_time)], n_domain)
    if (!identical(as.character(domain), as.character(dom_blocks)) ||
        !identical(as.character(time), as.character(time_blocks))) {
      cli::cli_abort(c(
        "{.arg data} must be a balanced panel: one row per domain and time period, with the same periods in every domain.",
        "i" = "Check for duplicated or missing domain x time combinations."
      ))
    }
  }

  # check dimension
  M <- n_domain * n_time
  if (nrow(X) != M || length(y) != M || length(vardir) != M) {
    cli::cli_abort(
      c("Dimensions mismatch!",
        "i" = "Formula produces {nrow(X)} rows, but domain * time = {M}.",
        "i" = "vardir has {length(vardir)} values."
      )
    )
  }

  # check W
  if (anyNA(W)) cli::cli_abort("Argument W contains NA values.")
  if (nrow(W) != n_domain || ncol(W) != n_domain) {
    cli::cli_abort(
      c("Argument W must be a square matrix of size Domain = {n_domain}.",
        "i" = "Got a {nrow(W)}x{ncol(W)} matrix."
      )
    )
  }

  # check sigma and rho
  if (!is.null(sigma21_start) && sigma21_start < 0) {
    cli::cli_abort("sigma21_start must be >= 0.")
  }
  if (!is.null(sigma22_start) && sigma22_start < 0) {
    cli::cli_abort("sigma22_start must be >= 0.")
  }
  if (rho1_start <= -1 || rho1_start >= 1) {
    cli::cli_abort("rho1_start must be in the interval (-1, 1).")
  }
  if (model == "ST" && (rho2_start <= -1 || rho2_start >= 1)) {
    cli::cli_abort("rho2_start must be in the interval (-1, 1).")
  }

  if (!all(sampled)) {
    dom_unsampled <- tapply(!sampled, factor(domain, levels = unique(domain)), all)
    n_dom0 <- sum(dom_unsampled)
    n_gap <- sum(!sampled) - n_dom0 * n_time
    if (print_result) {
      cli::cli_alert_info(paste0(
        "Predicting {sum(!sampled)} cell{?s} without a direct estimate: ",
        "{n_dom0} unsampled domain{?s} (spatial borrowing) and {n_gap} missing period{?s} of sampled domains."
      ))
    }
  }

  # ---- call C++ core ----
  res <- .eblup_stfh_core(
    Xall = X,
    yall = y,
    vardirall = vardir,
    proxmat = W,
    D = as.integer(n_domain),
    Tt = as.integer(n_time),
    model = model,
    maxiter = as.integer(maxiter),
    precision = precision,
    sigma21_start = if (is.null(sigma21_start)) -1 else sigma21_start,
    rho1_start = rho1_start,
    sigma22_start = if (is.null(sigma22_start)) -1 else sigma22_start,
    rho2_start = rho2_start
  )

  # ---- attach additional info (matching fastsae standardized structure) ----
  # Row names for estcoef
  if (!is.null(res$estcoef)) {
    row.names(res$estcoef) <- colnames(X)
    res$estcoef$stderr_beta <- res$estcoef$std.error
    res$estcoef$zvalue <- res$estcoef$tvalue
  }

  # Attach formula and metadata
  res$formula <- formula
  res$level <- "area"
  res$method <- "REML"

  # Add class
  res$call <- match.call()
  res$data <- data
  res$W <- W
  class(res) <- "fastsae"

  # Map results back to the caller's row order (see the sorting above)
  .restore_order <- function(res) {
    inv <- order(ord)
    if (!is.null(res$df_eblup) && nrow(res$df_eblup) == length(inv)) {
      res$df_eblup <- res$df_eblup[inv, , drop = FALSE]
      rownames(res$df_eblup) <- NULL
    }
    for (nm in c("mse_analytical", "rse_analytical", "g1", "g2", "g3", "g4")) {
      if (!is.null(res[[nm]]) && length(res[[nm]]) == length(inv)) res[[nm]] <- res[[nm]][inv]
    }
    if (!is.null(res$data) && nrow(res$data) == length(inv)) res$data <- res$data[inv, , drop = FALSE]
    res
  }

  # Convergence check
  if (!isTRUE(res$convergence)) {
    cli::cli_alert_danger(
      "After {res$n_iter} iteration(s), there is no convergence."
    )
    return(.restore_order(res))
  }

  # ---- compute MSE if requested ----
  if (compute_mse) {
    if (mse_type == "analytical") {
      res$df_eblup$mse <- as.numeric(res$mse_analytical)
      res$df_eblup$rse <- as.numeric(res$rse_analytical)
      res$df_eblup$mse_pb <- NA_real_
      res$mse_type <- "analytical"
      res$B <- NA_integer_
    } else {
      if (print_result) {
        cli::cli_alert_info("Computing parametric bootstrap MSE with B = {B} replicates...")
      }

      # Call C++ PBMSE function
      pbmse_res <- .pbmse_stfh(
        Xall = X,
        yall = y,
        vardirall = vardir,
        proxmat = W,
        D = as.integer(n_domain),
        Tt = as.integer(n_time),
        model = model,
        maxiter = as.integer(maxiter),
        precision = precision,
        B = as.integer(B),
        n_threads = as.integer(n_threads),
        seed = as.integer(seed)
      )

      # Update df_eblup with MSE
      res$df_eblup$mse_pb <- as.numeric(pbmse_res$mse_pb)
      res$df_eblup$mse <- as.numeric(pbmse_res$mse_pb)
      res$df_eblup$rse <- ifelse(abs(res$df_eblup$eblup) < .Machine$double.eps, NA_real_, sqrt(as.numeric(pbmse_res$mse_pb)) / abs(res$df_eblup$eblup) * 100)
      res$mse_type <- "bootstrap"
      res$B <- pbmse_res$B
    }
  } else {
    # Add NA columns for mse and rse
    res$df_eblup$mse_pb <- NA_real_
    res$df_eblup$mse <- NA_real_
    res$df_eblup$rse <- NA_real_
    res$B <- NA_integer_
  }

  # Standardize df_eblup columns
  res$df_eblup$domain <- domain
  res$df_eblup$time <- time
  res$df_eblup$y <- y
  res$df_eblup$vardir <- ifelse(sampled, vardir, NA_real_)
  cols_order <- c(
    "domain", "time", "y", "eblup", "vardir", "mse", "rse", "mse_pb",
    "random_effect_u1", "random_effect_u2", "sampled"
  )
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

  res <- .restore_order(res)

  # Print results
  if (print_result) {
    print(res)
  }

  return(res)
}
