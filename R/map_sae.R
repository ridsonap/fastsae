# ============================================================================
# Unified Spatial Choropleth Mapping Interface for fastsae (map_sae)
# ============================================================================

#' Unified Spatial Choropleth Mapping for Small Area Estimates
#'
#' @description
#' Provides a comprehensive, single-gateway ("satu pintu") spatial choropleth mapping
#' interface for Small Area Estimation models fitted with \pkg{fastsae} (e.g. \code{eblup_fh},
#' \code{eblup_sfh}, \code{eblup_stfh}, \code{eblup_bhf}, \code{ebp_area}, and \code{benchmark}).
#'
#' Intelligently handles:
#' \itemize{
#'   \item \strong{Single Model Visualizations}: Point estimates (\code{"estimate"}), Relative
#'     Standard Error (\code{"rse"}), and 3-color official statistics reliability classification
#'     (\code{"reliability"}: green for reliable < 20\%, amber for caution 20-30\%, red for unreliable \eqn{\ge 30\%}).
#'   \item \strong{Direct vs Model Comparisons}: Side-by-side facet maps contrasting erratic direct survey
#'     estimates against smoothed small area model estimates (\code{"comparison"}).
#'   \item \strong{Two-Model Comparisons}: Side-by-side maps or spatial difference maps (\code{"difference"})
#'     contrasting two competing models (e.g. Classical FH vs Spatial FH, or FastSAE vs Stan).
#'   \item \strong{Benchmarked Calibrations}: Visualizing pre- vs post-benchmarked calibrations.
#'   \item \strong{Smart Key Matching}: Automatically detects the common domain identifier between
#'     the model and \code{sf_geom}, with optional manual override via \code{key}.
#' }
#'
#' @param object A fitted \code{fastsae} model object, a \code{fastsae_benchmark} object,
#'   a named list of \code{fastsae} objects, or a data frame containing domain estimates.
#' @param model2 Optional second \code{fastsae} model object for direct model-to-model comparison.
#' @param sf_geom Spatial polygon geometry of class \code{sf}. Optional if the model was
#'   fitted on data that already inherits from \code{sf}.
#' @param key Optional character string or named vector specifying the column(s) used to match
#'   model domain identifiers with the spatial geometry in \code{sf_geom}.
#'   \itemize{
#'     \item If \code{NULL} (default), automatically detects the matching key by evaluating
#'       common column names or the highest set intersection overlap between domain IDs and \code{sf_geom} columns.
#'     \item If a single string (e.g. \code{key = "kd_kab"}), specifies the column in \code{sf_geom}
#'       to match with model domain IDs.
#'     \item If a named string (e.g. \code{key = c("domain" = "kd_kab")}), explicitly pairs the model's
#'       domain column with the spatial column in \code{sf_geom}.
#'   }
#' @param type Character string specifying the visualization type:
#'   \itemize{
#'     \item \code{"estimate"}: Choropleth map of small area point estimates (EBP or EBLUP).
#'     \item \code{"rse"}: Choropleth map of Relative Standard Errors / CV (\%).
#'     \item \code{"reliability"}: Official statistics 3-color traffic light map (< 20\% Reliable,
#'       20-30\% Use with Caution, \eqn{\ge 30\%} Unreliable).
#'     \item \code{"comparison"}: Side-by-side facet choropleth map (Direct vs Model, or Model 1 vs Model 2).
#'     \item \code{"difference"}: Spatial difference map (\eqn{\hat{\theta}^{(2)} - \hat{\theta}^{(1)}}
#'       or \eqn{\hat{\theta}^{\text{Model}} - \hat{\theta}^{\text{Direct}}}) with divergent color palette.
#'   }
#' @param indicator Optional alias for \code{type}. If supplied, overrides \code{type}.
#' @param palette Character string specifying a color palette. Default is automatically selected
#'   based on \code{type}: \code{"viridis"} for estimates, \code{"magma"} for RSE, official traffic-light
#'   colors for reliability, and \code{"RdBu"} / \code{"PuOr"} for differences.
#' @param thresholds Numeric vector of length 2 defining the RSE (\%) thresholds for reliability
#'   classification. Default is \code{c(20, 30)} according to official statistical guidelines (BPS / Eurostat).
#' @param facet_scales Character string passed to \code{ggplot2::facet_wrap(scales = ...)}. Default is \code{"fixed"}.
#' @param title Optional title for the plot. If \code{NULL}, a descriptive title is generated automatically.
#' @param subtitle Optional subtitle for the plot.
#' @param ... Additional arguments passed to \code{ggplot2} layers.
#'
#' @return A \code{ggplot} object containing the thematic choropleth map.
#'
#' @export
#' @examples
#' library(fastsae)
#' data(mys)
#'
#' # Create a synthetic spatial grid for demonstration
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   grid_sf <- sf::st_make_grid(
#'     sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
#'     cellsize = c(1, 1),
#'     what = "polygons"
#'   )[seq_len(nrow(mys))]
#'   mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)
#'
#'   # 1. Fit Fay-Herriot model
#'   fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
#'
#'   # 2. Single model point estimate map (auto-detected key)
#'   p1 <- map_sae(fit_fh, sf_geom = mys_sf)
#'
#'   # 3. Explicit key matching
#'   p2 <- map_sae(fit_fh, sf_geom = mys_sf, key = "area")
#'
#'   # 4. Reliability classification map (BPS standard <20%, 20-30%, >=30%)
#'   p3 <- map_sae(fit_fh, sf_geom = mys_sf, type = "reliability")
#'
#'   # 5. Side-by-side comparison map: Direct Survey vs Model SAE
#'   p4 <- map_sae(fit_fh, sf_geom = mys_sf, type = "comparison")
#' }
map_sae <- function(object, ...) {
  UseMethod("map_sae")
}

#' @rdname map_sae
#' @export
map_sae.fastsae <- function(
  object,
  model2 = NULL,
  sf_geom = NULL,
  key = NULL,
  type = c("estimate", "rse", "reliability", "comparison", "difference"),
  indicator = NULL,
  palette = NULL,
  thresholds = c(20, 30),
  facet_scales = "fixed",
  title = NULL,
  subtitle = NULL,
  ...
) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg sf} is required for spatial choropleth mapping.")
  }

  type_chosen <- if (!is.null(indicator)) indicator else match.arg(type)

  # Check if model2 is provided -> delegate to comparison mapping
  if (!is.null(model2)) {
    return(.map_two_models(
      model1 = object,
      model2 = model2,
      sf_geom = sf_geom,
      key = key,
      type = type_chosen,
      palette = palette,
      facet_scales = facet_scales,
      title = title,
      subtitle = subtitle,
      ...
    ))
  }

  # Extract single model estimation data
  df_model <- .extract_model_data_for_map(object)

  # Check if sf_geom is embedded in object$data
  sf_obj <- sf_geom
  if (is.null(sf_obj)) {
    if (!is.null(object$data) && inherits(object$data, "sf")) {
      sf_obj <- object$data
    } else if (inherits(object$df_ebp, "sf")) {
      sf_obj <- object$df_ebp
    } else if (inherits(object$df_eblup, "sf")) {
      sf_obj <- object$df_eblup
    }
  }

  if (is.null(sf_obj)) {
    cli::cli_abort(c(
      "Spatial polygon geometry {.arg sf_geom} must be provided when model data is not an {.cls sf} object.",
      "i" = "Pass your spatial polygons via {.code map_sae(fit, sf_geom = your_sf_polygons)}."
    ))
  }

  # Perform smart key matching
  keys <- .match_spatial_keys(model_df = df_model, sf_geom = sf_obj, key = key, domain_col = "domain")
  matched_sf <- .merge_model_with_sf(df_model, sf_obj, keys)

  # Plot according to type
  .generate_sae_map(
    matched_sf = matched_sf,
    type = type_chosen,
    model_name = attr(df_model, "model_name") %||% "SAE Model",
    palette = palette,
    thresholds = thresholds,
    title = title,
    subtitle = subtitle,
    ...
  )
}

#' @rdname map_sae
#' @export
map_sae.fastsae_benchmark <- function(
  object,
  sf_geom = NULL,
  key = NULL,
  type = c("comparison", "estimate", "difference"),
  indicator = NULL,
  palette = NULL,
  title = NULL,
  subtitle = NULL,
  ...
) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg sf} is required for spatial choropleth mapping.")
  }

  type_chosen <- if (!is.null(indicator)) indicator else match.arg(type)

  # Prepare data from benchmarked object
  df_bm <- as.data.frame(object)
  df_map <- data.frame(
    domain = df_bm$domain,
    direct = NA_real_,
    estimate = df_bm$benchmarked,
    original = df_bm$original,
    adjustment = df_bm$adjustment,
    rse = NA_real_,
    stringsAsFactors = FALSE
  )
  attr(df_map, "model_name") <- "Benchmarked"

  if (is.null(sf_geom)) {
    cli::cli_abort("Please provide spatial polygons via {.arg sf_geom} for benchmarked objects.")
  }

  keys <- .match_spatial_keys(model_df = df_map, sf_geom = sf_geom, key = key, domain_col = "domain")
  matched_sf <- .merge_model_with_sf(df_map, sf_geom, keys)

  if (type_chosen == "comparison") {
    # Faceted: Original vs Benchmarked
    long_data <- rbind(
      data.frame(matched_sf[, c("domain", "geometry")], Source = "Original Model", Value = matched_sf$original),
      data.frame(matched_sf[, c("domain", "geometry")], Source = "Benchmarked (Calibrated)", Value = matched_sf$estimate)
    )
    long_sf <- sf::st_as_sf(long_data)

    p <- ggplot2::ggplot(long_sf) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$Value), color = "gray30", linewidth = 0.2) +
      ggplot2::facet_wrap(~ Source) +
      ggplot2::scale_fill_viridis_c(option = palette %||% "viridis", name = "Value") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = title %||% "Benchmark Calibration Comparison",
        subtitle = subtitle %||% "Original Model Estimate vs Aggregate-Calibrated Estimate"
      ) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
    return(p)

  } else if (type_chosen == "difference") {
    matched_sf$diff_val <- matched_sf$adjustment
    p <- ggplot2::ggplot(matched_sf) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$diff_val), color = "gray30", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(
        low = "#2b83ba", mid = "#ffffbf", high = "#d7191c", midpoint = 0,
        name = "Adjustment"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = title %||% "Benchmarking Calibration Adjustment",
        subtitle = subtitle %||% "Difference: (Benchmarked - Original)"
      ) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
    return(p)

  } else {
    # Single estimate map of benchmarked values
    matched_sf$plot_val <- matched_sf$estimate
    p <- ggplot2::ggplot(matched_sf) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$plot_val), color = "gray30", linewidth = 0.2) +
      ggplot2::scale_fill_viridis_c(option = palette %||% "viridis", name = "Benchmarked") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = title %||% "Benchmarked Small Area Estimates",
        subtitle = subtitle %||% "Aggregate-Calibrated Values"
      ) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
    return(p)
  }
}

#' @rdname map_sae
#' @export
map_sae.list <- function(
  object,
  sf_geom = NULL,
  key = NULL,
  type = c("comparison", "difference"),
  indicator = NULL,
  palette = NULL,
  title = NULL,
  subtitle = NULL,
  ...
) {
  if (length(object) < 2) {
    if (length(object) == 1) return(map_sae(object[[1]], sf_geom = sf_geom, key = key, ...))
    cli::cli_abort("Empty list provided to {.fn map_sae}.")
  }

  type_chosen <- if (!is.null(indicator)) indicator else match.arg(type)

  model1 <- object[[1]]
  model2 <- object[[2]]
  names_list <- names(object) %||% c("Model 1", "Model 2")

  .map_two_models(
    model1 = model1,
    model2 = model2,
    model1_name = names_list[1],
    model2_name = names_list[2],
    sf_geom = sf_geom,
    key = key,
    type = type_chosen,
    palette = palette,
    title = title,
    subtitle = subtitle,
    ...
  )
}

#' @rdname map_sae
#' @export
map_sae.default <- function(
  object,
  sf_geom = NULL,
  key = NULL,
  type = c("estimate", "rse", "reliability"),
  indicator = NULL,
  palette = NULL,
  title = NULL,
  subtitle = NULL,
  ...
) {
  if (!is.data.frame(object)) {
    cli::cli_abort("{.arg object} must be a {.cls fastsae} object, a benchmark object, or a data frame.")
  }

  type_chosen <- if (!is.null(indicator)) indicator else match.arg(type)

  # Check required columns
  df <- object
  domain_candidates <- c("domain", "area", "id", "code", "kd_kab")
  found_dom <- intersect(domain_candidates, names(df))
  if (length(found_dom) == 0) {
    df$domain <- seq_len(nrow(df))
  } else {
    df$domain <- df[[found_dom[1]]]
  }

  est_candidates <- c("ebp", "eblup", "estimate", "est", "y_hat", "pred")
  found_est <- intersect(est_candidates, names(df))
  if (length(found_est) == 0) {
    cli::cli_abort("Data frame must contain an estimate column (e.g. {.code ebp}, {.code eblup}, or {.code estimate}).")
  }
  df$estimate <- df[[found_est[1]]]

  if (!"rse" %in% names(df)) {
    if ("mse" %in% names(df)) {
      df$rse <- (sqrt(df$mse) / abs(df$estimate)) * 100
    } else {
      df$rse <- rep(NA_real_, nrow(df))
    }
  }

  attr(df, "model_name") <- "Estimates"

  sf_obj <- sf_geom %||% (if (inherits(df, "sf")) df else NULL)
  if (is.null(sf_obj)) {
    cli::cli_abort("Please provide spatial polygons via {.arg sf_geom}.")
  }

  keys <- .match_spatial_keys(model_df = df, sf_geom = sf_obj, key = key, domain_col = "domain")
  matched_sf <- .merge_model_with_sf(df, sf_obj, keys)

  .generate_sae_map(
    matched_sf = matched_sf,
    type = type_chosen,
    model_name = "Estimates",
    palette = palette,
    title = title,
    subtitle = subtitle,
    ...
  )
}

# ------------------------------------------------------------------------------
# Internal Helpers for Map Processing
# ------------------------------------------------------------------------------

#' Internal worker to map two competing models
#' @noRd
.map_two_models <- function(
  model1,
  model2,
  model1_name = "Model 1",
  model2_name = "Model 2",
  sf_geom = NULL,
  key = NULL,
  type = "comparison",
  palette = NULL,
  facet_scales = "fixed",
  title = NULL,
  subtitle = NULL,
  ...
) {
  df1 <- .extract_model_data_for_map(model1)
  df2 <- .extract_model_data_for_map(model2)

  # Check sf_geom
  sf_obj <- sf_geom
  if (is.null(sf_obj)) {
    if (!is.null(model1$data) && inherits(model1$data, "sf")) sf_obj <- model1$data
    else if (!is.null(model2$data) && inherits(model2$data, "sf")) sf_obj <- model2$data
  }
  if (is.null(sf_obj)) {
    cli::cli_abort("Please provide spatial polygons via {.arg sf_geom} for two-model comparison.")
  }

  # Match keys
  keys1 <- .match_spatial_keys(model_df = df1, sf_geom = sf_obj, key = key, domain_col = "domain")
  matched_sf <- .merge_model_with_sf(df1, sf_obj, keys1)

  # Match model2
  rownames(df2) <- as.character(df2$domain)
  matched_doms <- as.character(matched_sf$domain)
  matched_sf$est2 <- df2[matched_doms, "estimate"]
  matched_sf$est1 <- matched_sf$estimate

  if (type == "difference") {
    matched_sf$diff_val <- matched_sf$est2 - matched_sf$est1
    p <- ggplot2::ggplot(matched_sf) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$diff_val), color = "gray30", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(
        low = "#2b83ba", mid = "#ffffbf", high = "#d7191c", midpoint = 0,
        name = "Difference"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = title %||% paste("Spatial Difference:", model2_name, "minus", model1_name),
        subtitle = subtitle %||% "Relative shift in domain estimates between models"
      ) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
    return(p)

  } else {
    # Default comparison: Faceted side-by-side
    long_data <- rbind(
      data.frame(matched_sf[, c("domain", "geometry")], Model = model1_name, Value = matched_sf$est1),
      data.frame(matched_sf[, c("domain", "geometry")], Model = model2_name, Value = matched_sf$est2)
    )
    long_sf <- sf::st_as_sf(long_data)

    p <- ggplot2::ggplot(long_sf) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$Value), color = "gray30", linewidth = 0.2) +
      ggplot2::facet_wrap(~ Model, scales = facet_scales) +
      ggplot2::scale_fill_viridis_c(option = palette %||% "viridis", name = "Estimate") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = title %||% paste("Model Comparison:", model1_name, "vs", model2_name),
        subtitle = subtitle %||% "Side-by-side geographic distribution of estimates"
      ) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
    return(p)
  }
}

#' Internal extractor of standardized prediction data frame
#' @noRd
.extract_model_data_for_map <- function(object) {
  df_est <- object$df_ebp %||% object$df_eblup
  if (is.null(df_est)) {
    cli::cli_abort("Model object does not contain estimation data ({.code df_ebp} or {.code df_eblup}).")
  }

  dom_col <- intersect(c("domain", "area", "id"), names(df_est))[1]
  domain_vals <- df_est[[dom_col]] %||% seq_len(nrow(df_est))

  y_raw <- if ("y" %in% names(df_est)) df_est$y else NA_real_
  y_hat <- df_est$ebp %||% df_est$eblup %||% df_est$est
  mse_val <- df_est$mse %||% (if (!is.null(df_est$sd)) df_est$sd^2 else NA_real_)

  rse_val <- df_est$rse
  if (is.null(rse_val) && !all(is.na(mse_val))) {
    rse_val <- (sqrt(mse_val) / abs(y_hat)) * 100
  }

  res <- data.frame(
    domain = as.character(domain_vals),
    direct = as.numeric(y_raw),
    estimate = as.numeric(y_hat),
    mse = as.numeric(mse_val),
    rse = as.numeric(rse_val),
    stringsAsFactors = FALSE
  )

  # Reliability flag (BPS official standards: <20%, 20-30%, >=30%)
  res$reliability <- cut(
    res$rse,
    breaks = c(-Inf, 20, 30, Inf),
    labels = c("Reliable (< 20%)", "Use with Caution (20-30%)", "Unreliable (\u2265 30%)"),
    right = FALSE
  )

  model_label <- switch(object$model %||% "SAE",
    "FH"  = "Fay-Herriot EBLUP",
    "SFH" = "Spatial FH EBLUP",
    "ST"  = "Spatio-Temporal FH EBLUP",
    "BHF" = "Battese-Harter-Fuller",
    "EBP" = paste0("EBP (", toupper(object$family %||% "Beta"), ")"),
    "SAE Model"
  )
  if (!is.null(object$family)) {
    model_label <- paste0("EBP-", toupper(object$family))
  }
  attr(res, "model_name") <- model_label

  return(res)
}

#' Smart Key Matching between model data and spatial geometry
#' @noRd
.match_spatial_keys <- function(model_df, sf_geom, key = NULL, domain_col = "domain") {
  sf_cols <- setdiff(names(sf_geom), attr(sf_geom, "sf_column"))
  model_doms <- as.character(model_df[[domain_col]])

  matched_sf_col <- NULL
  matched_model_col <- domain_col

  if (!is.null(key)) {
    if (is.character(key) && length(key) == 1 && is.null(names(key))) {
      # Single string specifies column in sf_geom
      if (!key %in% names(sf_geom)) {
        cli::cli_abort("Specified key column {.val {key}} was not found in {.arg sf_geom}.")
      }
      matched_sf_col <- key
    } else if (is.character(key) && !is.null(names(key))) {
      # Named key c(model_col = sf_col)
      m_col <- names(key)[1]
      s_col <- unname(key)[1]
      if (!m_col %in% names(model_df)) {
        cli::cli_abort("Model domain column {.val {m_col}} was not found in model results.")
      }
      if (!s_col %in% names(sf_geom)) {
        cli::cli_abort("Spatial key column {.val {s_col}} was not found in {.arg sf_geom}.")
      }
      matched_model_col <- m_col
      matched_sf_col <- s_col
    }
  }

  # Automatic key detection if not provided
  if (is.null(matched_sf_col)) {
    # 1. Exact name match
    if (domain_col %in% sf_cols) {
      matched_sf_col <- domain_col
    } else {
      # 2. Check candidate standard names
      cand_names <- c("domain", "area", "id", "code", "kd_kab", "kode", "kabupaten", "provinsi", "district", "region")
      found_cand <- intersect(cand_names, tolower(sf_cols))
      if (length(found_cand) > 0) {
        idx <- which(tolower(sf_cols) == found_cand[1])
        matched_sf_col <- sf_cols[idx[1]]
      }

      # 3. Best overlap matching (set intersection)
      best_overlap <- 0
      best_col <- NULL
      for (col in sf_cols) {
        col_vals <- as.character(sf_geom[[col]])
        common <- length(intersect(model_doms, col_vals))
        if (common > best_overlap) {
          best_overlap <- common
          best_col <- col
        }
      }

      # If best overlap matches at least 40% of domains, adopt it
      if (best_overlap >= (0.4 * length(model_doms))) {
        matched_sf_col <- best_col
      }
    }
  }

  if (is.null(matched_sf_col)) {
    cli::cli_abort(c(
      "Could not automatically match model domain IDs with any column in {.arg sf_geom}.",
      "i" = "Please specify the matching key explicitly via {.code key = 'column_name_in_sf'}."
    ))
  }

  list(model_col = matched_model_col, sf_col = matched_sf_col)
}

#' Merge model results with sf geometry
#' @noRd
.merge_model_with_sf <- function(model_df, sf_geom, keys) {
  m_col <- keys$model_col
  s_col <- keys$sf_col

  sf_work <- sf_geom
  sf_work$..merge_key.. <- as.character(sf_work[[s_col]])
  model_df$..merge_key.. <- as.character(model_df[[m_col]])

  # Merge
  merged <- merge(sf_work, model_df, by = "..merge_key..", all.x = TRUE)
  merged$..merge_key.. <- NULL
  return(merged)
}

#' Plotting engine for map_sae
#' @noRd
.generate_sae_map <- function(
  matched_sf,
  type,
  model_name = "SAE Model",
  palette = NULL,
  thresholds = c(20, 30),
  title = NULL,
  subtitle = NULL,
  ...
) {
  if (type == "estimate") {
    matched_sf$plot_val <- matched_sf$estimate
    pal <- palette %||% "viridis"

    p <- ggplot2::ggplot(matched_sf) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$plot_val), color = "gray30", linewidth = 0.2) +
      ggplot2::scale_fill_viridis_c(option = pal, name = "Estimate", na.value = "gray90") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = title %||% paste("Small Area Estimates (", model_name, ")", sep = ""),
        subtitle = subtitle %||% "Model-based predictions across domains"
      ) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

  } else if (type == "rse") {
    matched_sf$plot_val <- matched_sf$rse
    pal <- palette %||% "magma"

    p <- ggplot2::ggplot(matched_sf) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$plot_val), color = "gray30", linewidth = 0.2) +
      ggplot2::scale_fill_viridis_c(option = pal, direction = -1, name = "RSE (%)", na.value = "gray90") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = title %||% paste("Relative Standard Error (%) \u2014", model_name),
        subtitle = subtitle %||% "Sampling error and precision distribution across domains"
      ) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

  } else if (type == "reliability") {
    # 3-color traffic light scale (BPS / Eurostat guidelines)
    colors_traffic <- c(
      "Reliable (< 20%)"          = "#198754",  # Green
      "Use with Caution (20-30%)" = "#ffc107",  # Amber/Yellow
      "Unreliable (\u2265 30%)"   = "#dc3545"   # Red
    )

    p <- ggplot2::ggplot(matched_sf) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$reliability), color = "gray30", linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = colors_traffic, na.value = "gray90", name = "Status") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = title %||% paste("Statistical Reliability Classification \u2014", model_name),
        subtitle = subtitle %||% "Evaluation against official statistical publication thresholds"
      ) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

  } else if (type == "comparison") {
    # Faceted: Direct Survey vs Model SAE
    if (all(is.na(matched_sf$direct))) {
      cli::cli_warn("Direct estimates are not available in model object; displaying model estimate map instead.")
      return(.generate_sae_map(matched_sf, type = "estimate", model_name = model_name, palette = palette, title = title, subtitle = subtitle))
    }

    long_data <- rbind(
      data.frame(matched_sf[, c("domain", "geometry")], Source = "Direct Survey", Value = matched_sf$direct),
      data.frame(matched_sf[, c("domain", "geometry")], Source = paste("SAE (", model_name, ")", sep = ""), Value = matched_sf$estimate)
    )
    long_sf <- sf::st_as_sf(long_data)

    p <- ggplot2::ggplot(long_sf) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$Value), color = "gray30", linewidth = 0.2) +
      ggplot2::facet_wrap(~ Source) +
      ggplot2::scale_fill_viridis_c(option = palette %||% "viridis", name = "Estimate", na.value = "gray90") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = title %||% "Direct Survey vs Small Area Estimation",
        subtitle = subtitle %||% "Side-by-side geographical comparison of domain values"
      ) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

  } else if (type == "difference") {
    # Difference: Model - Direct
    if (all(is.na(matched_sf$direct))) {
      cli::cli_abort("Direct estimates are required for difference mapping against a single model.")
    }
    matched_sf$diff_val <- matched_sf$estimate - matched_sf$direct

    p <- ggplot2::ggplot(matched_sf) +
      ggplot2::geom_sf(ggplot2::aes(fill = .data$diff_val), color = "gray30", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(
        low = "#2b83ba", mid = "#ffffbf", high = "#d7191c", midpoint = 0,
        name = "Difference", na.value = "gray90"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::labs(
        title = title %||% "Estimation Difference: (Model - Direct)",
        subtitle = subtitle %||% "Spatial smoothing adjustments produced by small area model"
      ) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

  } else {
    cli::cli_abort("Unknown map type: {.val {type}}.")
  }

  return(p)
}
