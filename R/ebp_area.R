#' Area-level EBP using INLA
#'
#' @description
#' Area-level Empirical Best Predictor (EBP) using Integrated Nested Laplace
#' Approximation (INLA) for Bayesian inference. Supports multiple distributions
#' and spatial random effects.
#'
#' @references
#' \enumerate{
#'   \item Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
#'     inference for latent Gaussian models using integrated nested Laplace
#'     approximations (with discussion). Journal of the Royal Statistical
#'     Society B, 71(2), 319-392.
#'   \item Lindgren, F., & Rue, H. (2015). Bayesian spatial modelling with
#'     R-INLA. Journal of Statistical Software, 63(19), 1-25.
#' }
#'
#' @param formula an object of class formula that contains a description of
#'   the model to be fitted.
#' @param data a data frame or a data frame extension (e.g. a tibble).
#' @param family distribution family for the response. Options include:
#'   \itemize{
#'     \item \code{"gaussian"} - Normal distribution for continuous data
#'     \item \code{"binomial"} - Binomial distribution for proportion/binary data
#'     \item \code{"poisson"} - Poisson distribution for count data
#'     \item \code{"zip"} - Zero-inflated Poisson for overdispersed counts
#'     \item \code{"nbinomial"} - Negative binomial for overdispersed counts
#'   }
#' @param domain vector, column name or one-sided formula referencing a domain
#'   names column in \code{data}. If NULL, the domains are numbered consecutively.
#' @param spatial spatial random effects model. Options include:
#'   \itemize{
#'     \item \code{"none"} - No spatial random effects (iid)
#'     \item \code{"bym"} - Besag-York-Mollie model (iid + spatial)
#'     \item \code{"bym2"} - Proper BYM2 model with normalized spatial effects
#'     \item \code{"besagproper"} - Besag proper model (proper CAR)
#'   }
#' @param W A square matrix with dimension equal to the number of domains.
#'   Required when \code{spatial} is not "none". Should be a row-standardized
#'   proximity matrix.
#' @param N.trials number of trials for binomial family. Can be a vector or
#'   column name. Defaults to 1 for binary data.
#' @param E expected values for Poisson rate model (offset term).
#' @param control.compute list of options controlling computation. Set to
#'   \code{list(cpo = TRUE, dic = TRUE, waic = TRUE)} for model comparison.
#' @param control.predictor compute fitted values and linear predictors.
#' @param ... additional arguments passed to \code{\link[INLA]{inla}}.
#'
#' @returns The function returns a list with class \code{"fastsae"} containing:
#'   \itemize{
#'     \item \code{df_eblup} - Data frame with domain, y, eblup, mse, rse
#'     \item \code{estcoef} - Fixed effects with estimates, SE, z-values, p-values
#'     \item \code{random_effect_var} - Estimated variance of random effects
#'     \item \code{goodness} - Model fit statistics (DIC, WAIC, CPO)
#'     \item \code{convergence} - Whether INLA converged
#'     \item \code{family} - Distribution family used
#'     \item \code{spatial} - Spatial model used
#'   }
#'
#' @examples
#' \dontrun{
#' # Note: These examples require INLA to be installed.
#' # Install from: install.packages('INLA', repos =
#' #   'https://inla.r-inla-download.org/R/testing')
#'
#' # Gaussian model (non-spatial)
#' m1 <- ebp_area(
#'   y ~ x1 + x2 + x3,
#'   data = mys,
#'   family = "gaussian",
#'   spatial = "none"
#' )
#'
#' # Spatial Gaussian with BYM2
#' m2 <- ebp_area(
#'   y ~ x1 + x2 + x3,
#'   data = mys,
#'   family = "gaussian",
#'   spatial = "bym2",
#'   W = mys_proxmat
#' )
#'
#' # Binomial model
#' m3 <- ebp_area(
#'   y ~ x1 + x2,
#'   data = proportion_data,
#'   family = "binomial",
#'   N.trials = n,
#'   spatial = "none"
#' )
#'
#' # Poisson model with offset
#' m4 <- ebp_area(
#'   y ~ x1 + x2,
#'   data = count_data,
#'   family = "poisson",
#'   E = expected_counts,
#'   spatial = "bym"
#' )
#' }
#'
#' @export
ebp_area <- function(
  formula,
  data,
  family = c("gaussian", "binomial", "poisson", "zip", "nbinomial"),
  domain = NULL,
  spatial = c("none", "bym", "bym2", "besagproper"),
  W = NULL,
  N.trials = NULL,
  E = NULL,
  control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
  control.predictor = list(link = TRUE),
  print_result = TRUE,
  ...
) {
  # Check INLA availability
  .check_inla_available()

  # Match arguments
  family <- match.arg(family)
  spatial <- match.arg(spatial)

  # Validate domain
  if (is.null(domain)) {
    domain <- "area"
    data$area <- seq_len(nrow(data))
  } else if (is.character(domain) && length(domain) == 1) {
    domain <- domain
  } else {
    domain <- .get_variable(data, domain)
    if (length(domain) != nrow(data)) {
      cli::cli_abort("Length of 'domain' must equal nrow(data)")
    }
    domain_name <- if (is.character(substitute(domain))) {
      as.character(substitute(domain))
    } else {
      "area"
    }
    data$area <- domain
    domain <- "area"
  }

  # Create model frame
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  y <- stats::model.response(mf, "numeric")

  # Get design matrix
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  fixed_effects <- colnames(X)

  # Check auxiliary variables
  if (anyNA(X)) {
    cli::cli_abort("Auxiliary variables contain NA values.")
  }

  # Validate spatial requirements
  spatial_needs_W <- spatial != "none"
  if (spatial_needs_W && is.null(W)) {
    cli::cli_abort("Spatial weight matrix {.arg W} is required when spatial != 'none'.")
  }

  if (!is.null(W)) {
    n_domains <- length(unique(data[[domain]]))
    if (!is.matrix(W) || nrow(W) != n_domains || ncol(W) != n_domains) {
      cli::cli_abort(c(
        "`W` must be a square matrix with dimension equal to the number of domains ({n_domains}).",
        "i" = "Got a {nrow(W)}x{ncol(W)} matrix."
      ))
    }
  }

  # Prepare data for INLA
  domain_idx <- as.numeric(factor(data[[domain]]))
  n_domains <- length(unique(domain_idx))

  # Build data list with all variables from the model frame
  # INLA will use the formula to extract variables
  data_list <- as.list(mf)

  # Add domain index (use "idx" as the variable name for INLA::f)
  data_list[["idx"]] <- domain_idx

  # Add N.trials for binomial
  if (family == "binomial") {
    if (!is.null(N.trials)) {
      if (is.character(N.trials) && length(N.trials) == 1) {
        data_list$Ntrials <- data[[N.trials]]
      } else {
        data_list$Ntrials <- N.trials
      }
    } else {
      data_list$Ntrials <- rep(1, length(y))
    }
  }

  # Add offset for Poisson
  if (family == "poisson" && !is.null(E)) {
    if (is.character(E) && length(E) == 1) {
      data_list$E <- data[[E]]
    } else {
      data_list$E <- E
    }
  }

  # Build INLA formula
  # Fixed effects part - use original formula directly
  # INLA handles intercept automatically

  # Create temporary graph file for spatial models
  graph_file <- NULL
  if (spatial != "none" && !is.null(W)) {
    # Convert to binary adjacency matrix
    W_binary <- (W != 0) * 1
    diag(W_binary) <- 0

    # Save as graph file for INLA
    graph_file <- tempfile(fileext = ".graph")
    write.graph <- function(Wmat, file) {
      n <- nrow(Wmat)
      lines <- character(n)
      for (i in 1:n) {
        idx <- which(Wmat[i, ] > 0) - 1  # 0-indexed
        lines[i] <- paste(i - 1, length(idx), paste(idx, collapse = " "))
      }
      writeLines(c(as.character(n), lines), file)
    }
    write.graph(W_binary, graph_file)
  }

  # Random effects part - build formula string for all models
  random_part <- switch(spatial,
    "none" = sprintf("f(idx, model = 'iid')"),
    "bym" = sprintf("f(idx, model = 'bym', graph = '%s')", graph_file),
    "bym2" = sprintf("f(idx, model = 'bym2', graph = '%s', scale.model = TRUE)", graph_file),
    "besagproper" = sprintf("f(idx, model = 'besagproper', graph = '%s')", graph_file)
  )

  # Construct formula - use original formula and add random effects
  random_formula <- as.formula(paste0("~", random_part))
  inla_formula <- update.formula(formula, random_formula)

  # Map family names for INLA
  inla_family <- switch(family,
    "gaussian"   = "gaussian",
    "binomial"   = "binomial",
    "poisson"    = "poisson",
    "zip"        = "zeroinflatedpoisson0",
    "nbinomial"  = "nbinomial",
    family
  )

  # Build INLA call arguments
  inla_args <- list(
    formula = inla_formula,
    data = data_list,
    family = inla_family,
    control.compute = control.compute,
    control.predictor = control.predictor
  )

  # Add Ntrials for binomial
  if (family == "binomial") {
    inla_args$Ntrials <- data_list$Ntrials
  }

  # Call INLA
  result <- do.call(INLA::inla, inla_args)

  # Convert to fastsae object
  out <- .convert_inla_to_fastsae(
    result = result,
    data = data,
    domain = domain,
    y = y,
    family = family,
    spatial = spatial,
    formula = formula
  )

  # Add call metadata
  out$call <- match.call()

  # Print result if requested
  if (print_result) {
    print(out)
  }

  return(out)
}
