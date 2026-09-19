# ============================================================================
# S3 Methods for fastsae Objects
# ============================================================================

# Helper operator for default value
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Print a fastsae object
#'
#' @param x An object of class \code{fastsae}.
#' @param ... Additional arguments passed to print methods.
#'
#' @return The original object invisibly.
#' @export
print.fastsae <- function(x, ...) {
  cli::cli_h1("Fast Small Area Estimation (fastsae)")

  if (!is.null(x$call)) {
    cli::cli_par()
    cli::cli_text("{.strong Call}:")
    cli::cli_text(deparse(x$call))
    cli::cli_end()
  }

  if (!is.null(x$convergence)) {
    if (isTRUE(x$convergence)) {
      n_it <- x$n_iter %||% "-"
      cli::cli_alert_success("Convergence: Yes (in {n_it} iterations)")
    } else {
      cli::cli_alert_danger("Convergence: No")
    }
  }

  model_type <- switch(
    x$model %||% "FH",
    "FH"  = "Fay-Herriot (Area-level)",
    "SFH" = "Spatial Fay-Herriot (Area-level SAR)",
    "ST"  = "Spatio-Temporal Fay-Herriot (ST-FH)",
    "S"   = "Spatial Fay-Herriot (S-FH)",
    "BHF" = "Battese-Harter-Fuller (Unit-level)",
    x$model %||% "Fay-Herriot"
  )

  cli::cli_text("{.strong Model}: {model_type}")

  if (!is.null(x$method)) {
    cli::cli_text("{.strong Method}: {x$method}")
  }


  if (!is.null(x$random_effect_var)) {
    cat("Random effect variance (sigma2_u):", round(x$random_effect_var, 6), "\n")
  }
  if (!is.null(x$rho)) {
    cat("Spatial autocorrelation (rho):", round(x$rho, 4), "\n")
  }

  if (!is.null(x$estcoef)) {
    cat("\nFixed Effects Coefficients:\n")
    cols_to_print <- intersect(c("beta", "std.error", "zvalue", "pvalue"), names(x$estcoef))
    if (length(cols_to_print) == 0) cols_to_print <- names(x$estcoef)
    stats::printCoefmat(as.matrix(x$estcoef[, cols_to_print, drop = FALSE]),
                        signif.stars = TRUE, ...)
  }

  if (!is.null(x$df_eblup)) {
    cat("\nEBLUP Estimates (First 6 domains):\n")
    print(utils::head(x$df_eblup, 6), ...)
    if (nrow(x$df_eblup) > 6) {
      cat("... and", nrow(x$df_eblup) - 6, "more rows.\n")
    }
  }

  cat("\n")
  invisible(x)
}

#' Summarize a fastsae object
#'
#' @param object An object of class \code{fastsae}.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{summary.fastsae}.
#' @export
summary.fastsae <- function(object, ...) {
  ans <- list(
    call = object$call,
    model = object$model,
    method = object$method,
    convergence = object$convergence,
    n_iter = object$n_iter,
    coefficients = object$estcoef,
    estvarcomp = object$estvarcomp,
    random_effect_var = object$random_effect_var,
    rho = object$rho,
    goodness = object$goodness,
    df_eblup = object$df_eblup,
    level = object$level
  )
  class(ans) <- "summary.fastsae"
  ans
}

#' Print summary of a fastsae object
#'
#' @param x An object of class \code{summary.fastsae}.
#' @param ... Additional arguments.
#'
#' @return The original object invisibly.
#' @export
print.summary.fastsae <- function(x, ...) {
  cli::cli_h1("Summary of fastsae Fit")

  if (!is.null(x$call)) {
    cli::cli_par()
    cli::cli_text("Call :")
    cli::cli_text(deparse(x$call))
    cli::cli_end()
  }

  model_type <- switch(
    x$model %||% "FH",
    "FH"  = "Fay-Herriot (Area-level)",
    "SFH" = "Spatial Fay-Herriot (Area-level SAR)",
    "ST"  = "Spatio-Temporal Fay-Herriot (ST-FH)",
    "S"   = "Spatial Fay-Herriot (S-FH)",
    "BHF" = "Battese-Harter-Fuller (Unit-level)",
    x$model %||% "Fay-Herriot"
  )

  if (!is.null(x$convergence)) {
    if (isTRUE(x$convergence)) {
      cli::cli_alert_success("Convergence: Yes (in {x$n_iter} iterations)")
    } else {
      cli::cli_alert_danger("Convergence: No")
    }
  }

  cli::cli_text("{.strong Model}: {model_type}")

  if (!is.null(x$method)) {
    cli::cli_text("{.strong Method}: {x$method}")
  }

  if (!is.null(x$estvarcomp)) {
    cat("\nVariance & Correlation Components:\n")
    print(x$estvarcomp, row.names = FALSE)
  } else if (!is.null(x$random_effect_var)) {
    cat("\nVariance Components:\n")
    cat("sigma2_u:", round(x$random_effect_var, 6), "\n")
  }

  if (!is.null(x$coefficients)) {
    cat("\nCoefficients:\n")
    cols_to_print <- intersect(c("beta", "std.error", "zvalue", "pvalue"), names(x$coefficients))
    if (length(cols_to_print) == 0) cols_to_print <- names(x$coefficients)
    stats::printCoefmat(as.matrix(x$coefficients[, cols_to_print, drop = FALSE]),
                        signif.stars = TRUE, ...)
  }

  if (!is.null(x$goodness)) {
    cat("\nGoodness of Fit:\n")
    print(x$goodness)
  }

  if (!is.null(x$df_eblup)) {
    cat("\nEBLUP Summary Statistics:\n")
    summary_cols <- intersect(c("eblup", "mse", "rse"), names(x$df_eblup))
    if (length(summary_cols) > 0) {
      print(summary(x$df_eblup[, summary_cols, drop = FALSE]))
    }
  }

  cat("\n")
  invisible(x)
}

#' Extract coefficients from a fastsae object
#'
#' @param object An object of class \code{fastsae}.
#' @param ... Additional arguments.
#'
#' @return Named vector of estimated coefficients.
#' @export
coef.fastsae <- function(object, ...) {
  if (!is.null(object$estcoef) && "beta" %in% names(object$estcoef)) {
    return(stats::setNames(object$estcoef$beta, rownames(object$estcoef)))
  }
  if (!is.null(object$fit$beta)) {
    return(as.vector(object$fit$beta))
  }
  NULL
}

#' Extract fitted values (EBLUP) from a fastsae object
#'
#' @param object An object of class \code{fastsae}.
#' @param ... Additional arguments.
#'
#' @return Vector of fitted EBLUP estimates.
#' @export
fitted.fastsae <- function(object, ...) {
  if (!is.null(object$df_eblup) && "eblup" %in% names(object$df_eblup)) {
    return(object$df_eblup$eblup)
  }
  NULL
}

#' Extract residuals from a fastsae object
#'
#' @param object An object of class \code{fastsae}.
#' @param ... Additional arguments.
#'
#' @return Vector of residuals (y - eblup) for sampled areas.
#' @export
residuals.fastsae <- function(object, ...) {
  if (!is.null(object$df_eblup) && all(c("y", "eblup") %in% names(object$df_eblup))) {
    return(object$df_eblup$y - object$df_eblup$eblup)
  }
  if (!is.null(object$fit$lme)) {
    return(stats::residuals(object$fit$lme))
  }
  NULL
}
