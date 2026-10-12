#' Empirical Best Linear Unbiased Prediction (EBLUP) for the Battese-Harter-Fuller Model
#'
#' This function estimates small area means or totals using the unit-level model
#' proposed by Battese, Harter, and Fuller (1988), which combines survey data
#' (sample units) and auxiliary population information.
#'
#' @param formula An object of class `formula` describing the model.
#' @param unit_data A `data.frame` containing the unit-level survey data.
#' @param Xpop A `data.frame` containing auxiliary variables and domain info.
#' @param domain_var A character string giving the column name for domain identifier.
#' @param popsize_var A character string for population size variable.
#' @param method Fitting method: "REML" (default) or "ML".
#' @param popnmean_xpop Population mean of auxiliary variables per domain.
#' @param B Number of bootstrap replicates for MSE (if compute_mse = TRUE).
#' @param compute_mse If TRUE, compute bootstrap MSE.
#' @param n_threads Currently unused: the bootstrap MSE runs sequentially (each replicate
#'   is a \code{lme4} refit). Kept for interface compatibility.
#' @param seed Random seed for reproducibility.
#' @param print_result Print results (default TRUE).
#'
#' @returns List containing EBLUP estimates, fit, and optionally MSE.
#'
#' @references
#' Battese, G. E., Harter, R. M., and Fuller, W. A. (1988). An error-components
#' model for prediction of county crop areas using survey and satellite data.
#' *Journal of the American Statistical Association*, 83(401), 28-36.
#'
#' @examples
#'
#' library(dplyr)
#' df_meanpop <- cornsoybeanmeans |>
#'   rename(CornPix = MeanCornPixPerSeg, SoyBeansPix = MeanSoyBeansPixPerSeg)
#' df_cornsoybean <- cornsoybean |>
#'   rename(CountyIndex = County)
#'
#' res <- eblup_bhf(
#'   formula = CornHec ~ CornPix + SoyBeansPix,
#'   Xpop = df_meanpop,
#'   unit_data = df_cornsoybean,
#'   domain_var = "CountyIndex",
#'   popsize_var = "PopnSegments"
#' )
#'
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
  compute_mse = FALSE,
  n_threads = 1,
  seed = -1,
  print_result = TRUE
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

  # --- fit linear mixed model with method (REML vs ML) ---
  term_labels <- attr(stats::terms(formula), "term.labels")
  lmer_data <- as.data.frame(formuladata)
  lmer_data$dom <- dom
  lmer_data$..y_var.. <- ys
  fixed_part <- if (length(term_labels) > 0) paste(term_labels, collapse = " + ") else "1"
  formula_lmer <- stats::as.formula(paste("..y_var.. ~", fixed_part, "+ (1 | dom)"))

  fit <- lme4::lmer(formula_lmer, data = lmer_data, REML = (method == "REML"))
  betaest <- matrix(lme4::fixef(fit), ncol = 1)
  upred <- lme4::ranef(fit)$dom

  # --- prepare population data sorted to match selectdom ---
  pop_dom <- .get_variable(Xpop, domain_var)
  selectdom <- unique(c(selectdom, pop_dom))
  Xpop <- Xpop[order(match(pop_dom, selectdom)), , drop = FALSE]
  popnsize <- .get_variable(Xpop, popsize_var)

  if (is.null(popnmean_xpop)) {
    formula_noy <- stats::reformulate(attr(stats::terms(formula), "term.labels"))
    meanxpop <- stats::model.matrix(formula_noy, Xpop)
  } else {
    meanxpop <- as.matrix(popnmean_xpop)
    if (ncol(meanxpop) == ncol(Xs) - 1 && "(Intercept)" %in% colnames(Xs)) {
      meanxpop <- cbind(`(Intercept)` = 1, meanxpop)
    }
  }

  if (ncol(meanxpop) != ncol(Xs)) {
    cli::cli_abort("Number of columns in auxiliary population means ({ncol(meanxpop)}) must match design matrix ({ncol(Xs)}).")
  }

  # --- extract variance components ---
  sigma2_u <- lme4::VarCorr(fit)$dom[1, 1]
  sigma2_e <- attr(lme4::VarCorr(fit), "sc")^2

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

    mse_result <- .pbmse_unit(
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
    )

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

  # Build standardized estcoef
  lmer_sum <- summary(fit)
  coef_mat <- as.data.frame(lmer_sum$coefficients)
  estcoef <- data.frame(
    beta = coef_mat[, 1],
    std.error = coef_mat[, 2],
    stderr_beta = coef_mat[, 2],
    tvalue = coef_mat[, 3],
    zvalue = coef_mat[, 3],
    pvalue = 2 * stats::pnorm(abs(coef_mat[, 3]), lower.tail = FALSE),
    row.names = rownames(coef_mat)
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
      lme = fit
    ),
    formula = formula,
    model = "BHF",
    level = "unit",
    convergence = TRUE,
    data = unit_data,
    call = match.call()
  )

  class(out) <- c("fastsae", "fastsae_unit")

  if (print_result) {
    print(out)
  }

  return(out)
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

  # Use the same population means as the point estimate (incl. a user-supplied
  # popnmean_xpop), so that the bootstrap targets the same estimator.
  if (is.null(meanxpop)) {
    formula_noy <- stats::reformulate(attr(stats::terms(formula), "term.labels"))
    meanxpop <- stats::model.matrix(formula_noy, Xpop)
    if (ncol(meanxpop) == ncol(Xs) - 1 && "(Intercept)" %in% colnames(Xs)) {
      meanxpop <- cbind(`(Intercept)` = 1, meanxpop)
    }
  }

  # --- initial fit ---
  term_labels <- attr(stats::terms(formula), "term.labels")
  lmer_data <- as.data.frame(formuladata)
  lmer_data$dom <- dom
  lmer_data$..y_boot.. <- ys
  fixed_part <- if (length(term_labels) > 0) paste(term_labels, collapse = " + ") else "1"
  formula_lmer <- stats::as.formula(paste("..y_boot.. ~", fixed_part, "+ (1 | dom)"))

  fit <- lme4::lmer(formula_lmer, data = lmer_data, REML = (method == "REML"))
  betaest <- matrix(lme4::fixef(fit), ncol = 1)
  upred <- lme4::ranef(fit)$dom

  sigma2_u <- lme4::VarCorr(fit)$dom[1, 1]
  sigma2_e <- attr(lme4::VarCorr(fit), "sc")^2

  # --- initial EBLUP ---
  init_result <- .eblup_bhf_cpp(
    selectdom = as.character(selectdom),
    dom = as.character(dom),
    Xs = Xs,
    meanxpop = meanxpop,
    ys = ys,
    popnsize = popnsize,
    betaest = betaest,
    upred = upred
  )

  eblup_init <- init_result$eblup
  names(eblup_init) <- as.character(selectdom)

  # --- bootstrap loop ---
  if (seed >= 0) set.seed(seed)

  mse <- numeric(length(selectdom))
  valid_B <- 0L

  # Domain bookkeeping, computed once (the loop below is vectorised over domains)
  dom_factor <- factor(dom, levels = as.character(selectdom))
  Xbeta_pop <- as.numeric(Xs %*% betaest)
  n_d <- as.numeric(table(dom_factor))                 # sample size per domain
  N_d <- as.numeric(popnsize)
  r_d <- pmax(N_d - n_d, 0)                            # non-sampled units
  mu_d <- as.numeric(meanxpop %*% betaest)             # Xbar_d' beta
  # population mean error: ebar_d = (n_d ebar_{d,s} + r_d ebar_{d,r}) / N_d,
  # ebar_{d,r} ~ N(0, s2e / r_d); fully sampled domains (r_d = 0): ebar_{d,s};
  # unsampled domains (n_d = 0): ebar_d ~ N(0, s2e / N_d)
  w_s <- ifelse(n_d > 0, ifelse(r_d > 0, n_d / N_d, 1), 0)
  w_r <- ifelse(n_d > 0, ifelse(r_d > 0, r_d / N_d, 0), 1)
  sd_r <- ifelse(n_d > 0, ifelse(r_d > 0, sqrt(sigma2_e / r_d), 0), sqrt(sigma2_e / N_d))
  dom_int <- as.integer(dom_factor)

  for (b in seq_len(B)) {
    # Generate bootstrap sample
    u_boot <- stats::rnorm(length(selectdom), sd = sqrt(sigma2_u))

    e_boot <- stats::rnorm(length(ys), sd = sqrt(sigma2_e))

    # Fast vectorized y_boot: X * beta + u_boot[dom] + e_boot
    y_boot <- Xbeta_pop + u_boot[dom_int] + e_boot   # integer index (no name matching)

    # Re-fit the bootstrap response with the model structure of the initial
    # fit (lme4::refit): same estimator as lmer() on the new response, without
    # re-parsing the formula and rebuilding the model matrices every time
    fit_boot <- tryCatch(
      lme4::refit(fit, newresp = y_boot),
      error = function(e) NULL
    )
    if (is.null(fit_boot)) next

    beta_boot <- matrix(lme4::fixef(fit_boot), ncol = 1)
    upred_boot <- lme4::ranef(fit_boot)$dom

    # Compute bootstrap EBLUP
    eblup_boot <- .eblup_bhf_cpp(
      selectdom = as.character(selectdom),
      dom = as.character(dom),
      Xs = Xs,
      meanxpop = meanxpop,
      ys = y_boot,
      popnsize = popnsize,
      betaest = beta_boot,
      upred = upred_boot
    )$eblup

    # True population mean in the bootstrap population (vectorised over domains)
    esd <- numeric(length(n_d))
    has_s <- n_d > 0
    esd[has_s] <- rowsum(e_boot, dom_int)[, 1] / n_d[has_s]   # groups sorted by domain code
    e_r <- stats::rnorm(length(n_d), 0, sd_r)
    truemean_boot <- mu_d + as.numeric(u_boot) + w_s * esd + w_r * e_r

    # Accumulate MSE: (eblup_boot - truemean_boot)^2
    mse <- mse + (eblup_boot - truemean_boot)^2
    valid_B <- valid_B + 1L
  }

  if (valid_B > 0L) {
    mse <- mse / valid_B
  } else {
    mse <- rep(NA_real_, length(selectdom))
  }
  names(mse) <- as.character(selectdom)

  return(list(
    eblup = eblup_init,
    mse = mse,
    method = method,
    B = valid_B
  ))
}
