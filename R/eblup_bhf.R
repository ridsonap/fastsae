#' Empirical Best Linear Unbiased Prediction (EBLUP) for the Battese-Harter-Fuller Model
#'
#' This function estimates small area means using the unit-level model
#' proposed by Battese, Harter, and Fuller (1988), which combines survey data
#' (sample units) and auxiliary population information.
#'
#' @param formula An object of class `formula` describing the model.
#' @param unit_data A `data.frame` containing the unit-level survey data.
#' @param Xpop A `data.frame` with one row per domain: the domain identifier
#'   (`domain_var`), the population size (`popsize_var`) and the population means of the
#'   covariates, in columns named as the covariates in `formula`. Domains in `Xpop` without
#'   sampled units are predicted synthetically.
#' @param domain_var A character string giving the column name for domain identifier.
#' @param popsize_var A character string for population size variable.
#' @param method Fitting method: "REML" (default) or "ML".
#' @param popnmean_xpop Optional matrix or data frame of the population means of the
#'   covariates (one row per domain, in the order of `Xpop`, columns as in the design matrix,
#'   with or without the intercept); if `NULL`, they are taken from `Xpop`.
#' @param B Number of bootstrap replicates for MSE (if mse_type = "bootstrap").
#' @param mse_type MSE estimator: \code{"bootstrap"} (default; parametric bootstrap of
#'   Gonzalez-Manteiga et al. (2008), as in \code{sae::pbmseBHF()}) or \code{"none"} (no MSE;
#'   \code{mse} and \code{rse} are \code{NA}).
#' @param compute_mse Deprecated (fastsae 1.0.0); use \code{mse_type}.
#'   \code{compute_mse = TRUE} is \code{mse_type = "bootstrap"} and \code{FALSE} is
#'   \code{"none"}.
#' @param n_threads Number of OpenMP threads for the bootstrap MSE. Results depend on
#'   \code{seed} only, not on the number of threads.
#' @param seed Random seed for reproducibility; \code{-1} (default) sets no seed. A fixed seed
#'   does not change the random number stream of the session.
#' @param print_result Print results (default TRUE).
#'
#' @details
#' The nested error model \eqn{y_{dj} = x_{dj}'\beta + u_d + e_{dj}} is fitted by
#' REML or ML in compiled code. With \eqn{\gamma = \sigma_u^2/\sigma_e^2},
#' \eqn{\beta} and \eqn{\sigma_e^2} are profiled out, so the fit is a
#' one-dimensional maximisation over \eqn{\gamma \ge 0} (including the boundary
#' \eqn{\sigma_u^2 = 0}); each evaluation uses only within-domain cross-products
#' and domain means. The estimates are those of \code{lme4::lmer} and
#' \code{sae::eblupBHF} up to the optimiser tolerance of those packages. The
#' standard errors in \code{estcoef} are from \eqn{(X'\hat V^{-1}X)^{-1}}.
#'
#' The bootstrap MSE refits the model in each replicate with the same code and
#' runs replicates in parallel (\code{n_threads}); random draws are generated
#' sequentially from R's RNG, so a given \code{seed} gives the same result for
#' any number of threads.
#'
#' @returns An object of class \code{c("fastsae", "fastsae_unit")}, a list with
#'   \code{df_eblup} (one row per domain: \code{domain}, \code{eblup}, \code{samp_size},
#'   \code{mse}, \code{rse}), \code{estcoef} (coefficients, standard errors, z-statistics and
#'   p-values), \code{random_effect_var}, \code{fit} (the fitting method, the variance
#'   components \code{random_effect_var} and \code{sigma2_e}, \code{beta}, the predicted random
#'   effects and the unit-level conditional residuals), \code{convergence}, \code{mse_type},
#'   \code{formula}, \code{model}, \code{data} and \code{call}.
#'   \code{goodness} is a named vector: \code{loglikelihood} (the full Gaussian
#'   log-likelihood of the sampled units at the estimates, for ML and REML fits alike, as in
#'   \code{\link{eblup_fh}()}), and \code{AIC} and \code{BIC} computed from it with
#'   \eqn{p + 2} parameters and the number of sampled units.
#'
#' @references
#' Battese, G. E., Harter, R. M., and Fuller, W. A. (1988). An error-components
#' model for prediction of county crop areas using survey and satellite data.
#' *Journal of the American Statistical Association*, 83(401), 28-36.
#'
#' Gonzalez-Manteiga, W., Lombardia, M. J., Molina, I., Morales, D., and Santamaria, L.
#' (2008). Bootstrap mean squared error of a small-area EBLUP. *Journal of Statistical
#' Computation and Simulation*, 78(5), 443-462.
#'
#' @seealso \code{\link{hb_unit}} for the hierarchical Bayes counterpart.
#'
#' @examples
#'
#' # the population means of the covariates must carry the names used in the formula
#' df_meanpop <- cornsoybeanmeans
#' names(df_meanpop)[match(c("MeanCornPixPerSeg", "MeanSoyBeansPixPerSeg"), names(df_meanpop))] <-
#'   c("CornPix", "SoyBeansPix")
#' df_cornsoybean <- cornsoybean
#' names(df_cornsoybean)[names(df_cornsoybean) == "County"] <- "CountyIndex"
#'
#' res <- eblup_bhf(
#'   formula = CornHec ~ CornPix + SoyBeansPix,
#'   Xpop = df_meanpop,
#'   unit_data = df_cornsoybean,
#'   domain_var = "CountyIndex",
#'   popsize_var = "PopnSegments"
#' )
#'
#' @md
#' @export
eblup_bhf <- function(
  formula,
  unit_data,
  Xpop,
  domain_var,
  popsize_var,
  method = c("REML", "ML"),
  popnmean_xpop = NULL,
  B = 100,
  mse_type = c("bootstrap", "none"),
  n_threads = 1,
  seed = -1,
  print_result = TRUE,
  compute_mse = lifecycle::deprecated()
) {
  method <- match.arg(method, choices = c("REML", "ML"))
  if (lifecycle::is_present(compute_mse)) {
    lifecycle::deprecate_warn("1.0.0", "eblup_bhf(compute_mse)", "eblup_bhf(mse_type)")
    mse_type <- if (isTRUE(compute_mse)) "bootstrap" else "none"
  }
  mse_type <- match.arg(mse_type, choices = c("bootstrap", "none"))
  compute_mse <- mse_type == "bootstrap"

  # --- data preparation ---
  formuladata <- stats::model.frame(formula, na.action = stats::na.omit, unit_data)
  dom <- .get_variable(unit_data, domain_var)
  if (!is.null(attr(formuladata, "na.action"))) {
    dom <- dom[-attr(formuladata, "na.action")]
  }
  selectdom <- unique(dom)

  # --- build design matrix and response ---
  ys <- stats::model.response(formuladata)
  Xs <- stats::model.matrix(formula, formuladata)

  # --- prepare population data sorted to match selectdom ---
  pop_dom <- .get_variable(Xpop, domain_var)
  selectdom <- unique(c(selectdom, pop_dom))
  dom_int <- match(dom, selectdom) - 1L

  # --- fit the nested error model (REML/ML, profiled in s2u/s2e; C++) ---
  fit <- .bhf_fit(Xs, ys, dom_int, length(selectdom), method)
  betaest <- matrix(fit$beta, ncol = 1, dimnames = list(colnames(Xs), NULL))
  upred <- .bhf_upred(fit, selectdom)

  Xpop <- Xpop[order(match(pop_dom, selectdom)), , drop = FALSE]
  popnsize <- .get_variable(Xpop, popsize_var)

  if (is.null(popnmean_xpop)) {
    formula_noy <- stats::reformulate(attr(stats::terms(formula), "term.labels"))
    meanxpop <- .pop_means(formula_noy, Xpop)
  } else {
    meanxpop <- as.matrix(popnmean_xpop)
    if (ncol(meanxpop) == ncol(Xs) - 1 && "(Intercept)" %in% colnames(Xs)) {
      meanxpop <- cbind(`(Intercept)` = 1, meanxpop)
    }
  }

  if (ncol(meanxpop) != ncol(Xs)) {
    cli::cli_abort("Number of columns in auxiliary population means ({ncol(meanxpop)}) must match design matrix ({ncol(Xs)}).")
  }
  .check_xpop_na(meanxpop, .get_variable(Xpop, domain_var))

  sigma2_u <- fit$sigma2_u
  sigma2_e <- fit$sigma2_e

  # --- compute EBLUP ---
  result <- .eblup_bhf_cpp(
    selectdom = as.character(selectdom),
    dom = as.character(dom),
    Xs = Xs,
    meanxpop = meanxpop,
    ys = ys,
    popnsize = popnsize,
    betaest = betaest,
    upred = upred
  )

  eblup_df <- data.frame(
    domain = selectdom,
    eblup = result$eblup,
    samp_size = result$samp_size
  )

  if (length(result$warn_domains) > 0) {
    cli::cli_warn(
      "{length(result$warn_domains)} domain(s) have no sample units: {paste(result$warn_domains, collapse = ', ')}"
    )
  }

  # --- bootstrap MSE ---
  if (compute_mse) {
    if (print_result) {
      cli::cli_alert_info("Computing bootstrap MSE with B = {B} replicates...")
    }

    mse_result <- .with_seed(seed, .pbmse_unit(
      formula = formula,
      unit_data = unit_data,
      Xpop = Xpop,
      domain_var = domain_var,
      popsize_var = popsize_var,
      meanxpop = meanxpop,
      method = method,
      B = B,
      n_threads = n_threads,
      seed = seed
    ))

    mse_df <- data.frame(
      domain = selectdom,
      mse = mse_result$mse
    )

    eblup_df <- merge(eblup_df, mse_df, by = "domain", sort = FALSE)
  }

  df_eblup <- eblup_df
  if (compute_mse && "mse" %in% names(df_eblup)) {
    df_eblup$rse <- ifelse(abs(df_eblup$eblup) < .Machine$double.eps, NA_real_,
      sqrt(df_eblup$mse) / abs(df_eblup$eblup) * 100
    )
  } else {
    df_eblup$mse <- NA_real_
    df_eblup$rse <- NA_real_
  }

  # Build standardized estcoef (Var(beta) = (X' V^-1 X)^-1 at the estimates)
  se_beta <- sqrt(diag(fit$vcov))
  tval <- fit$beta / se_beta
  estcoef <- data.frame(
    beta = fit$beta,
    std.error = se_beta,
    stderr_beta = se_beta,
    tvalue = tval,
    zvalue = tval,
    pvalue = 2 * stats::pnorm(abs(tval), lower.tail = FALSE),
    row.names = colnames(Xs)
  )

  # --- assemble output ---
  out <- list(
    df_eblup = df_eblup,
    eblup = df_eblup, # backward compatibility
    estcoef = estcoef,
    random_effect_var = sigma2_u,
    fit = list(
      method = method,
      random_effect_var = sigma2_u,
      sigma2_e = sigma2_e,
      beta = betaest,
      random_effect = upred,
      residuals = as.numeric(ys - Xs %*% betaest - fit$u[dom_int + 1L])
    ),
    formula = formula,
    model = "BHF",
    level = "unit",
    mse_type = mse_type,
    goodness = .bhf_goodness(ys, Xs, betaest, dom_int, sigma2_u, sigma2_e),
    convergence = fit$convergence,
    data = unit_data,
    call = match.call()
  )

  class(out) <- c("fastsae", "fastsae_unit")

  if (!isTRUE(fit$convergence)) .warn_no_convergence(NULL)

  if (print_result) {
    print(out)
  }

  return(out)
}

#' Full Gaussian log-likelihood of the nested error model at the estimates, with AIC
#' and BIC (p + 2 parameters, n sampled units). Per domain V_d = s2e I + s2u 11', so
#' log|V_d| = (n_d - 1) log s2e + log(s2e + n_d s2u) and
#' r' V_d^-1 r = (sum r^2 - s2u (sum r)^2 / (s2e + n_d s2u)) / s2e.
#' @noRd
.bhf_goodness <- function(ys, Xs, beta, dom_int, s2u, s2e) {
  r <- as.numeric(ys - Xs %*% beta)
  n_d <- tabulate(dom_int + 1L)
  sum_r <- vapply(split(r, dom_int), sum, numeric(1))
  sum_r2 <- vapply(split(r^2, dom_int), sum, numeric(1))
  n_d <- n_d[n_d > 0]
  logdet <- sum((n_d - 1) * log(s2e) + log(s2e + n_d * s2u))
  quad <- sum((sum_r2 - s2u * sum_r^2 / (s2e + n_d * s2u)) / s2e)
  n <- length(r)
  loglik <- -0.5 * (n * log(2 * pi) + logdet + quad)
  k <- ncol(Xs) + 2
  c(loglikelihood = loglik, AIC = -2 * loglik + 2 * k, BIC = -2 * loglik + k * log(n))
}

#' Fit the nested error model by REML/ML (C++, profiled likelihood)
#' @noRd
.bhf_fit <- function(Xs, ys, dom_int, n_dom, method) {
  n <- nrow(Xs)
  .check_design(Xs, "sampled units")
  if (length(unique(dom_int)) >= n || n <= ncol(Xs)) {
    cli::cli_abort("The number of sampled units must exceed both the number of sampled domains and the number of coefficients.")
  }
  if (length(unique(dom_int)) < 2) {
    cli::cli_abort("At least 2 sampled domains are needed to estimate the domain variance.")
  }
  .bhf_fit_cpp(Xs, as.numeric(ys), as.integer(dom_int), as.integer(n_dom), method == "REML")
}

#' Predicted random effects of the sampled domains as the data frame used by
#' .eblup_bhf_cpp
#' @noRd
.bhf_upred <- function(fit, selectdom) {
  sampled <- !is.na(fit$u)
  data.frame(
    `(Intercept)` = fit$u[sampled], row.names = as.character(selectdom[sampled]),
    check.names = FALSE
  )
}

#' @noRd
.pbmse_unit <- function(
  formula,
  unit_data,
  Xpop,
  domain_var,
  popsize_var,
  method = c("REML", "ML"),
  B = 100,
  n_threads = 1,
  seed = -1,
  meanxpop = NULL
) {
  method <- match.arg(method, choices = c("REML", "ML"))

  # --- data preparation ---
  formuladata <- stats::model.frame(formula, na.action = stats::na.omit, unit_data)
  dom <- .get_variable(unit_data, domain_var)
  if (!is.null(attr(formuladata, "na.action"))) {
    dom <- dom[-attr(formuladata, "na.action")]
  }
  selectdom <- unique(dom)

  # --- build design matrix and response ---
  ys <- stats::model.response(formuladata)
  Xs <- stats::model.matrix(formula, formuladata)

  # --- prepare population data sorted to match selectdom ---
  pop_dom <- .get_variable(Xpop, domain_var)
  selectdom <- unique(c(selectdom, pop_dom))
  Xpop <- Xpop[order(match(pop_dom, selectdom)), , drop = FALSE]
  popnsize <- .get_variable(Xpop, popsize_var)
  dom_int <- match(dom, selectdom) - 1L

  # Use the same population means as the point estimate (incl. a user-supplied
  # popnmean_xpop), so that the bootstrap targets the same estimator.
  if (is.null(meanxpop)) {
    formula_noy <- stats::reformulate(attr(stats::terms(formula), "term.labels"))
    meanxpop <- .pop_means(formula_noy, Xpop)
    if (ncol(meanxpop) == ncol(Xs) - 1 && "(Intercept)" %in% colnames(Xs)) {
      meanxpop <- cbind(`(Intercept)` = 1, meanxpop)
    }
  }
  meanxpop <- as.matrix(meanxpop)
  storage.mode(meanxpop) <- "double"

  # --- initial fit and EBLUP ---
  fit <- .bhf_fit(Xs, ys, dom_int, length(selectdom), method)
  eblup_init <- .eblup_bhf_cpp(
    selectdom = as.character(selectdom),
    dom = as.character(dom),
    Xs = Xs,
    meanxpop = meanxpop,
    ys = ys,
    popnsize = popnsize,
    betaest = matrix(fit$beta, ncol = 1),
    upred = .bhf_upred(fit, selectdom)
  )$eblup
  names(eblup_init) <- as.character(selectdom)

  # --- bootstrap (C++): y* = X beta + u* + e*, refit, EBLUP vs. true mean ---
  # Draws come from R's RNG in a fixed order, so results depend on the seed
  # only, not on n_threads.
  if (seed >= 0) set.seed(seed)
  boot <- .pbmse_bhf_cpp(
    Xs, as.integer(dom_int), length(selectdom), meanxpop, as.numeric(popnsize),
    fit$beta, fit$sigma2_u, fit$sigma2_e, as.integer(B), method == "REML",
    as.integer(n_threads)
  )
  mse <- as.numeric(boot$mse)
  names(mse) <- as.character(selectdom)

  return(list(
    eblup = eblup_init,
    mse = mse,
    method = method,
    B = boot$B
  ))
}
