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

  model_type <- switch(x$model %||% "FH",
    "FH"  = "Fay-Herriot (Area-level)",
    "SFH" = "Spatial Fay-Herriot (Area-level SAR)",
    "ST"  = "Spatio-Temporal Fay-Herriot (ST-FH)",
    "S"   = "Spatial Fay-Herriot (S-FH)",
    "BHF" = "Battese-Harter-Fuller (Unit-level)",
    "TFH" = "Two-fold Fay-Herriot (Area/Subarea-level)",
    "TWOFOLD" = "Two-fold Fay-Herriot (Area/Subarea-level)",
    x$model %||% "Fay-Herriot"
  )

  cli::cli_text("{.strong Model}: {model_type}")

  if (!is.null(x$temporal) && x$temporal != "none") {
    cli::cli_text("{.strong Temporal}: {toupper(x$temporal)}")
    if (!is.null(x$st_interaction) && x$st_interaction != "none") {
      cli::cli_text("{.strong Spatio-Temporal Structure}: {toupper(x$st_interaction)}")
    }
  }

  if (!is.null(x$method)) {
    cli::cli_text("{.strong Method}: {x$method}")
  }

  if (!is.null(x$random_effect_var)) {
    rev <- x$random_effect_var
    if (length(rev) == 2L && !is.null(names(rev))) {
      cat(paste0("Area effect variance (", names(rev)[1], "):"), round(rev[1], 6), "\n")
      cat(paste0("Subarea effect variance (", names(rev)[2], "):"), round(rev[2], 6), "\n")
    } else {
      cat("Random effect variance (sigma2_u):", round(rev, 6), "\n")
    }
  }
  if (!is.null(x$random_effect_var_time)) {
    cat("Temporal effect variance (sigma2_t):", round(x$random_effect_var_time, 6), "\n")
  }
  if (!is.null(x$rho)) {
    cat("Spatial autocorrelation (rho):", round(x$rho, 4), "\n")
  }
  if (!is.null(x$rho_time)) {
    cat("Temporal autocorrelation (rho_t):", round(x$rho_time, 4), "\n")
  }
  if (!is.null(x$phi)) {
    cat("Spatial mixing fraction (phi):", round(x$phi, 4), "\n")
  }

  if (!is.null(x$estcoef)) {
    cat("\nFixed Effects Coefficients:\n")
    cols_to_print <- intersect(c("beta", "std.error", "zvalue", "pvalue", "ci_lower", "ci_upper"), names(x$estcoef))
    if (length(cols_to_print) == 0) cols_to_print <- names(x$estcoef)
    stats::printCoefmat(as.matrix(x$estcoef[, cols_to_print, drop = FALSE]),
      signif.stars = TRUE, ...
    )
  }

  est_df <- x$df_hb %||% x$df_ebp %||% x$df_eblup
  if (!is.null(est_df)) {
    label_est <- if (!is.null(x$df_hb)) "HB Estimates" else if (!is.null(x$df_ebp)) "EBP Estimates" else "EBLUP Estimates"
    cat("\n", label_est, " (First 6 domains):\n", sep = "")
    print(utils::head(est_df, 6), ...)
    if (nrow(est_df) > 6) {
      cat("... and", nrow(est_df) - 6, "more rows.\n")
    }
  }

  if (!is.null(x$df_area)) {
    cat("\nArea-Level Aggregates (First 6 areas):\n")
    print(utils::head(x$df_area, 6), ...)
    if (nrow(x$df_area) > 6) {
      cat("... and", nrow(x$df_area) - 6, "more areas.\n")
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
    family = object$family,
    spatial = object$spatial,
    temporal = object$temporal,
    st_interaction = object$st_interaction,
    convergence = object$convergence,
    n_iter = object$n_iter,
    coefficients = object$estcoef,
    estvarcomp = object$estvarcomp,
    hyperpar = object$hyperpar,
    random_effect_var = object$random_effect_var,
    random_effect_var_time = object$random_effect_var_time,
    rho = object$rho,
    rho_time = object$rho_time,
    phi = object$phi,
    goodness = object$goodness,
    df_eblup = object$df_eblup,
    df_hb = object$df_hb,
    df_area = object$df_area,
    df_ebp = object$df_ebp,
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

  model_type <- switch(x$model %||% "FH",
    "FH"  = "Fay-Herriot (Area-level)",
    "SFH" = "Spatial Fay-Herriot (Area-level SAR)",
    "ST"  = "Spatio-Temporal Fay-Herriot (ST-FH)",
    "S"   = "Spatial Fay-Herriot (S-FH)",
    "BHF" = "Battese-Harter-Fuller (Unit-level)",
    "TFH" = "Two-fold Fay-Herriot (Area/Subarea-level)",
    "TWOFOLD" = "Two-fold Fay-Herriot (Area/Subarea-level)",
    x$model %||% "Fay-Herriot"
  )

  if (!is.null(x$convergence)) {
    if (isTRUE(x$convergence)) {
      n_it <- x$n_iter %||% "-"
      cli::cli_alert_success("Convergence: Yes (in {n_it} iterations)")
    } else {
      cli::cli_alert_danger("Convergence: No")
    }
  }

  cli::cli_text("{.strong Model}: {model_type}")

  if (!is.null(x$temporal) && x$temporal != "none") {
    cli::cli_text("{.strong Temporal}: {toupper(x$temporal)}")
    if (!is.null(x$st_interaction) && x$st_interaction != "none") {
      cli::cli_text("{.strong Spatio-Temporal Structure}: {toupper(x$st_interaction)}")
    }
  }

  if (!is.null(x$method)) {
    cli::cli_text("{.strong Method}: {x$method}")
  }

  if (!is.null(x$estvarcomp)) {
    cat("\nVariance & Correlation Components:\n")
    print(x$estvarcomp, row.names = FALSE)
  } else if (!is.null(x$random_effect_var)) {
    cat("\nVariance Components:\n")
    rev <- x$random_effect_var
    if (length(rev) == 2L && !is.null(names(rev))) {
      cat(paste0(names(rev)[1], " (area):"), round(rev[1], 6), "\n")
      cat(paste0(names(rev)[2], " (subarea):"), round(rev[2], 6), "\n")
    } else {
      cat("sigma2_u:", round(rev, 6), "\n")
    }
  }
  if (!is.null(x$random_effect_var_time)) {
    cat("sigma2_t (temporal):", round(x$random_effect_var_time, 6), "\n")
  }
  if (!is.null(x$rho)) {
    cat("rho (spatial):", round(x$rho, 4), "\n")
  }
  if (!is.null(x$rho_time)) {
    cat("rho_t (temporal):", round(x$rho_time, 4), "\n")
  }
  if (!is.null(x$phi)) {
    cat("phi (spatial fraction):", round(x$phi, 4), "\n")
  }

  if (!is.null(x$coefficients)) {
    cat("\nCoefficients:\n")
    cols_to_print <- intersect(c("beta", "std.error", "zvalue", "pvalue", "ci_lower", "ci_upper"), names(x$coefficients))
    if (length(cols_to_print) == 0) cols_to_print <- names(x$coefficients)
    stats::printCoefmat(as.matrix(x$coefficients[, cols_to_print, drop = FALSE]),
      signif.stars = TRUE, ...
    )
  }

  if (!is.null(x$hyperpar) && nrow(x$hyperpar) > 0) {
    cat("\nHyperparameters:\n")
    print(x$hyperpar)
  }

  if (!is.null(x$goodness)) {
    cat("\nGoodness of Fit:\n")
    print(x$goodness)
  }

  est_df <- x$df_hb %||% x$df_ebp %||% x$df_eblup
  if (!is.null(est_df)) {
    label_sum <- if (!is.null(x$df_hb)) "HB Summary Statistics" else if (!is.null(x$df_ebp)) "EBP Summary Statistics" else "EBLUP Summary Statistics"
    cat("\n", label_sum, ":\n", sep = "")
    summary_cols <- intersect(c("hb", "ebp", "linear_pred", "sd", "mse", "rse"), names(est_df))
    if (length(summary_cols) == 0) summary_cols <- intersect(c("eblup", "mse", "rse"), names(est_df))
    if (length(summary_cols) > 0) {
      print(summary(est_df[, summary_cols, drop = FALSE]))
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

#' Extract fitted values (EBLUP or HB) from a fastsae object
#'
#' @param object An object of class \code{fastsae}.
#' @param ... Additional arguments.
#'
#' @return Vector of fitted estimates.
#' @export
fitted.fastsae <- function(object, ...) {
  if (!is.null(object$df_hb) && "hb" %in% names(object$df_hb)) {
    return(object$df_hb$hb)
  }
  if (!is.null(object$df_ebp) && "ebp" %in% names(object$df_ebp)) {
    return(object$df_ebp$ebp)
  }
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
#' @return Vector of residuals for sampled areas.
#' @export
residuals.fastsae <- function(object, ...) {
  if (!is.null(object$df_hb) && all(c("y", "hb") %in% names(object$df_hb))) {
    return(object$df_hb$y - object$df_hb$hb)
  }
  if (!is.null(object$df_ebp) && all(c("y", "ebp") %in% names(object$df_ebp))) {
    return(object$df_ebp$y - object$df_ebp$ebp)
  }
  if (!is.null(object$df_eblup) && all(c("y", "eblup") %in% names(object$df_eblup))) {
    return(object$df_eblup$y - object$df_eblup$eblup)
  }
  if (!is.null(object$fit$lme)) {
    return(stats::residuals(object$fit$lme))
  }
  NULL
}

#' Plot method for fastsae objects
#'
#' @description
#' Provides standard R \code{plot()} dispatch for \code{fastsae} models.
#' Automatically routes to \code{\link{map_sae}} if spatial geometry is supplied,
#' to two-model comparison if a second \code{fastsae} model is provided, or
#' to \code{\link[ggplot2]{autoplot}} for diagnostic and estimate plots.
#'
#' @param x An object of class \code{fastsae}.
#' @param y Optional second object (e.g. another \code{fastsae} model for comparison).
#' @param ... Additional arguments passed to \code{\link{map_sae}} or \code{\link[ggplot2]{autoplot}}.
#'
#' @return A \code{ggplot} object.
#' @export
plot.fastsae <- function(x, y = NULL, ...) {
  dots <- list(...)
  if (!is.null(y) && inherits(y, "fastsae")) {
    if (!is.null(dots$sf_geom) || (!is.null(x$data) && inherits(x$data, "sf"))) {
      return(map_sae(x, model2 = y, ...))
    }
    return(ggplot2::autoplot(list("Model 1" = x, "Model 2" = y), ...))
  }
  if (!is.null(dots$sf_geom) || (!is.null(x$data) && inherits(x$data, "sf"))) {
    return(map_sae(x, ...))
  }
  ggplot2::autoplot(x, ...)
}

