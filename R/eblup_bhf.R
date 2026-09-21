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
#' @param n_threads Number of threads for parallel computation.
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
  lmer_data <- data.frame(
    y = ys,
    formuladata[, term_labels, drop = FALSE],
    dom = dom
  )
  fixed_part <- paste(term_labels, collapse = " + ")
  formula_lmer <- stats::as.formula(paste("y ~", fixed_part, "+ (1 | dom)"))

  fit <- lme4::lmer(formula_lmer, data = lmer_data, REML = (method == "REML"))
  betaest <- matrix(lme4::fixef(fit), ncol = 1)
  upred <- lme4::ranef(fit)$dom

  # --- prepare population data sorted to match selectdom ---
  pop_dom <- .get_variable(Xpop, domain_var)
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
  seed = -1
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
  Xpop <- Xpop[order(match(pop_dom, selectdom)), , drop = FALSE]
  popnsize <- .get_variable(Xpop, popsize_var)

  formula_noy <- stats::reformulate(attr(stats::terms(formula), "term.labels"))
  meanxpop <- stats::model.matrix(formula_noy, Xpop)
  if (ncol(meanxpop) == ncol(Xs) - 1 && "(Intercept)" %in% colnames(Xs)) {
    meanxpop <- cbind(`(Intercept)` = 1, meanxpop)
  }

  # --- initial fit ---
  term_labels <- attr(stats::terms(formula), "term.labels")
  lmer_data <- data.frame(
    y = ys,
    formuladata[, term_labels, drop = FALSE],
    dom = dom
  )
  fixed_part <- paste(term_labels, collapse = " + ")
  formula_lmer <- stats::as.formula(paste("y ~", fixed_part, "+ (1 | dom)"))

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

  # Group indices by domain
  dom_factor <- factor(dom, levels = as.character(selectdom))
  dom_idx <- split(seq_along(dom), dom_factor)

  for (b in seq_len(B)) {
    # Generate bootstrap sample
    u_boot <- stats::rnorm(length(selectdom), sd = sqrt(sigma2_u))
    names(u_boot) <- as.character(selectdom)

    e_boot <- stats::rnorm(length(ys), sd = sqrt(sigma2_e))

    # y_boot = X * beta + u_boot[dom] + e_boot
    y_boot <- numeric(length(ys))
    for (i in seq_along(ys)) {
      d <- as.character(dom[i])
      y_boot[i] <- as.numeric(Xs[i, ] %*% betaest) + u_boot[d] + e_boot[i]
    }

    # Fit bootstrap model
    lmer_data_boot <- data.frame(
      y = y_boot,
      formuladata[, term_labels, drop = FALSE],
      dom = dom
    )
    fit_boot <- tryCatch(
      lme4::lmer(formula_lmer, data = lmer_data_boot, REML = (method == "REML")),
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

    # Compute true population mean in bootstrap population
    truemean_boot <- numeric(length(selectdom))
    for (i in seq_along(selectdom)) {
      d <- as.character(selectdom[i])
      idx_d <- dom_idx[[d]]
      nd <- length(idx_d)
      Ni <- popnsize[i]
      rd <- max(Ni - nd, 0)
      mud <- sum(meanxpop[i, ] * betaest)
      if (nd > 0 && rd > 0) {
        esdmean <- mean(e_boot[idx_d])
        erdmean <- stats::rnorm(1, 0, sqrt(sigma2_e / rd))
        edmean <- esdmean * (nd / Ni) + erdmean * (rd / Ni)
        truemean_boot[i] <- mud + u_boot[d] + edmean
      } else if (nd > 0 && rd == 0) {
        esdmean <- mean(e_boot[idx_d])
        truemean_boot[i] <- mud + u_boot[d] + esdmean
      } else {
        truemean_boot[i] <- mud + u_boot[d] + stats::rnorm(1, 0, sqrt(sigma2_e / Ni))
      }
    }

    # Accumulate MSE: (eblup_boot - truemean_boot)^2
    mse <- mse + (eblup_boot - truemean_boot)^2
  }

  mse <- mse / B
  names(mse) <- as.character(selectdom)

  return(list(
    eblup = eblup_init,
    mse = mse,
    method = method,
    B = B
  ))
}
