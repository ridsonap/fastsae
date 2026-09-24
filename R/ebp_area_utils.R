# ============================================================================
# Area-level EBLUP/EBP using INLA (Integrated Nested Laplace Approximation)
# ============================================================================

#' Check if INLA package is available and provide installation guidance
#'
#' @noRd
.check_inla_available <- function() {
  if (!requireNamespace("INLA", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg INLA} is required but not installed.",
      "i" = "INLA is not available on CRAN and must be installed from a special repository.",
      " " = "",
      "To install INLA, run:",
      "  install.packages('INLA', repos = 'https://inla.r-inla-download.org/R/testing')",
      " " = "",
      "For more information, visit: https://www.r-inla.org/download/install"
    ))
  }
  invisible(TRUE)
}

#' Build INLA formula for area-level models
#'
#' @param formula Model formula
#' @param domain Domain identifier
#' @param spatial Spatial model type: "none", "bym", "bym2"
#' @param W Spatial weight matrix (required if spatial != "none")
#' @param family Distribution family
#'
#' @noRd
.build_inla_formula <- function(formula, domain, spatial, W = NULL, family = "gaussian") {
  # Extract response and fixed effects from formula
  formula_parts <- as.list(formula)
  response <- deparse(formula_parts[[2]])
  fixed_effects <- deparse(formula_parts[[3]])

  # Build random effects based on spatial model
  random_effect <- switch(spatial,
    "none" = sprintf("f(%s, model = 'iid')", domain),
    "bym" = {
      if (is.null(W)) {
        cli::cli_abort("Spatial weight matrix {.arg W} is required when spatial = 'bym' or 'bym2'")
      }
      sprintf("f(%s, model = 'bym', graph = W)", domain)
    },
    "bym2" = {
      if (is.null(W)) {
        cli::cli_abort("Spatial weight matrix {.arg W} is required when spatial = 'bym' or 'bym2'")
      }
      sprintf("f(%s, model = 'bym2', graph = W, scale.model = TRUE)", domain)
    }
  )

  # Construct the full formula
  full_formula_str <- sprintf("%s ~ %s + %s", response, fixed_effects, random_effect)
  formula(as.formula(full_formula_str))
}

#' Convert INLA result to fastsae object
#'
#' @param result INLA result object
#' @param data Original data frame
#' @param domain Domain column
#' @param y Response variable
#' @param family Distribution family
#' @param spatial Spatial model type
#' @param formula Original formula
#'
#' @noRd
.convert_inla_to_fastsae <- function(result, data, domain, y, family, spatial, formula) {

  # Extract fixed effects
  fixed_summary <- result$summary.fixed
  beta <- fixed_summary$mean
  se <- fixed_summary$sd
  zvalue <- beta / se
  pvalue <- 2 * pnorm(abs(zvalue), lower.tail = FALSE)

  estcoef <- data.frame(
    beta = beta,
    std.error = se,
    zvalue = zvalue,
    pvalue = pvalue,
    stringsAsFactors = FALSE
  )
  rownames(estcoef) <- rownames(fixed_summary)

  # Extract variance components from hyperparameters
  hyper_summary <- result$summary.hyperpar

  # Extract random effect variance
  if (spatial == "none") {
    # For iid model: Precision for iid
    random_var <- 1 / hyper_summary$mean[1]
  } else if (spatial == "bym") {
    # For BYM: Sum of iid and spatial variances
    # hyperparam order: Precision for iid, Precision for spatial
    random_var <- sum(1 / hyper_summary$mean[1:2])
  } else if (spatial == "bym2") {
    # For BYM2: Total variance
    # hyperparam: Theta (on log scale)
    random_var <- exp(hyper_summary$mean[1])
  } else if (spatial == "besagproper") {
    # For Besag proper: precision parameter
    random_var <- 1 / hyper_summary$mean[1]
  } else {
    random_var <- NA
  }

  # Extract mixing/spatial parameter for spatial models
  rho <- NULL
  if (spatial == "bym") {
    # rho = precision_iid / (precision_iid + precision_spatial)
    prec_iid <- hyper_summary$mean[1]
    prec_spatial <- hyper_summary$mean[2]
    rho <- prec_iid / (prec_iid + prec_spatial)
  } else if (spatial == "bym2") {
    # phi = mixing parameter (0 = iid, 1 = spatial)
    rho <- INLA::inla.emarginal(function(x) x, result$marginals.hyperpar[[1]])
  }

  # Extract EBLUP (fitted values)
  fitted_summary <- result$summary.fitted.values
  eblup <- fitted_summary$mean
  mse <- fitted_summary$sd^2
  rse <- sqrt(mse) / abs(eblup) * 100

  # Build df_eblup
  domain_vals <- data[[domain]]
  df_eblup <- data.frame(
    domain = domain_vals,
    y = y,
    eblup = eblup,
    mse = mse,
    rse = rse,
    stringsAsFactors = FALSE
  )

  # Compute goodness of fit (DIC, WAIC)
  goodness <- c(
    dic = result$dic$dic,
    waic = result$waic$waic,
    CPO = -mean(log(result$cpo$cpo), na.rm = TRUE)
  )

  # Build output object
  out <- list(
    estcoef = estcoef,
    estvarcomp = NULL,
    random_effect_var = random_var,
    goodness = goodness,
    df_eblup = df_eblup,
    convergence = result$converged,
    n_iter = NA,
    model = "INLA",
    method = "INLA",
    family = family,
    spatial = spatial,
    formula = formula,
    rho = rho,
    call = match.call()
  )

  class(out) <- "fastsae"
  return(out)
}

#' Prepare data for INLA
#'
#' @param data Data frame
#' @param formula Model formula
#' @param domain Domain column name
#'
#' @noRd
.prepare_inla_data <- function(data, formula, domain) {
  # Create model frame
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)

  # Get response
  y <- stats::model.response(mf, "numeric")

  # Get design matrix
  X <- stats::model.matrix(attr(mf, "terms"), mf)

  # Domain indicator
  domain_col <- data[[domain]]

  # Build data list for INLA
  data_list <- list(
    y = y,
    intercept = rep(1, length(y))
  )

  # Add domain index
  domain_idx <- as.numeric(factor(domain_col))
  data_list[[domain]] <- domain_idx

  # Add covariates (excluding intercept)
  covar_names <- setdiff(colnames(X), "(Intercept)")
  for (nm in covar_names) {
    data_list[[nm]] <- X[, nm]
  }

  attr(data_list, "domain") <- domain
  attr(data_list, "y") <- y

  return(data_list)
}

#' Get graph from spatial weight matrix for INLA
#'
#' @param W Row-standardized spatial weight matrix
#'
#' @noRd
.get_inla_graph <- function(W) {
  # Convert to binary adjacency matrix (0/1)
  # INLA accepts either a matrix or uses the graph= argument
  W_binary <- (W != 0) * 1
  diag(W_binary) <- 0  # Remove self-loops

  return(W_binary)
}
