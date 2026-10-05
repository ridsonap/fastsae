#' Create Spatial Proximity and Weights Matrix
#'
#' @description
#' Constructs spatial proximity, contiguity, and weight matrices (\eqn{\mathbf{W}})
#' from spatial data (\code{sf} polygon or point objects, coordinate matrices, or neighbor lists).
#' Supports Queen and Rook contiguity, \eqn{k}-nearest neighbors (\eqn{k}-NN),
#' distance bands, and inverse distance weighting (IDW), with automatic fallback
#' handling for isolated islands and disconnected domains.
#'
#' @param data An \code{sf} polygon or point object, a numeric matrix or data frame of
#'   2D coordinates (with columns \code{x}/\code{y} or longitude/latitude), or an \code{spdep}
#'   \code{nb} or \code{listw} object.
#' @param method Character string specifying the spatial weight method:
#'   \itemize{
#'     \item \code{"queen"}: Queen contiguity (sharing at least one vertex or boundary edge).
#'     \item \code{"rook"}: Rook contiguity (sharing a common linear boundary edge).
#'     \item \code{"knn"}: \eqn{k}-nearest neighbors based on centroid/point distances.
#'     \item \code{"distance"}: Distance band where domains within \code{d_max} are neighbors.
#'     \item \code{"idw"}: Inverse distance weighting \eqn{w_{ij} = 1 / d_{ij}^\alpha}.
#'   }
#' @param style Character string specifying the normalization style:
#'   \itemize{
#'     \item \code{"W"}: Row-standardized (row sums equal 1, default for SAR/SFH models).
#'     \item \code{"B"}: Binary adjacency (0/1 values, standard for INLA BYM2/Besag models).
#'     \item \code{"C"}: Globally standardized.
#'     \item \code{"U"}: Equal weighting (\eqn{1 / m}).
#'     \item \code{"minmax"}: Min-max spectral normalization (divided by largest eigenvalue).
#'   }
#' @param k Integer. Number of nearest neighbors when \code{method = "knn"} (default is 4).
#' @param d_max Numeric. Distance threshold when \code{method = "distance"}. If \code{NULL},
#'   it is automatically set to the minimum distance ensuring no disconnected domains.
#' @param alpha Numeric. Distance decay power for \code{method = "idw"} (default is 1).
#' @param island_fallback Logical. If \code{TRUE} (default), disconnected domains or islands
#'   with zero contiguous neighbors are automatically connected to their nearest neighbor(s)
#'   based on centroid distance.
#' @param k_fallback Integer. Number of nearest neighbors to connect for isolated islands
#'   when \code{island_fallback = TRUE} (default is 1).
#' @param domain Optional vector, column name, or one-sided formula referencing domain IDs.
#'   If provided, used as row and column names for the resulting matrix.
#'
#' @return A square numeric matrix of class \code{c("fastsae_weights", "matrix")} of dimension
#'   \eqn{m \times m}. Attributes include \code{method}, \code{style}, \code{n_domains},
#'   \code{n_islands}, and \code{call}.
#'
#' @details
#' When working with real administrative boundaries (such as Indonesian provinces, regencies,
#' or municipalities), island regencies (e.g. Kepulauan Seribu, Natuna, Sabang, Mentawai,
#' Sangihe) lack contiguous land borders. If \code{island_fallback = TRUE}, \code{create_weights}
#' detects these zero-neighbor domains and connects each to its closest \code{k_fallback}
#' neighbor(s) based on centroid distances, ensuring a connected spatial graph without
#' computational singularities.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   # Assuming 'data_kab' is an sf polygon object:
#'   W_queen <- create_weights(data_kab, method = "queen", style = "W")
#'   W_knn   <- create_weights(data_kab, method = "knn", k = 5, style = "W")
#'   W_bin   <- create_weights(data_kab, method = "queen", style = "B")
#' }
#' }
#'
#' # Create spatial weights from 2D coordinates
#' coords <- matrix(runif(40), ncol = 2)
#' W_coord <- create_weights(coords, method = "knn", k = 3, style = "W")
#' dim(W_coord)
#'
#' @export
create_weights <- function(
  data,
  method = c("queen", "rook", "knn", "distance", "idw"),
  style = c("W", "B", "C", "U", "minmax"),
  k = 4,
  d_max = NULL,
  alpha = 1,
  island_fallback = TRUE,
  k_fallback = 1,
  domain = NULL
) {
  call <- match.call()
  method <- match.arg(method)
  style <- match.arg(style)

  # 1. Determine domain names if provided
  dom_names <- NULL
  if (!is.null(domain)) {
    if (is.character(domain) && length(domain) == 1 && is.data.frame(data) && domain %in% names(data)) {
      dom_names <- as.character(data[[domain]])
    } else if (inherits(domain, "formula")) {
      dom_names <- as.character(stats::model.frame(domain, data = as.data.frame(data))[[1]])
    } else if (length(domain) == nrow(data)) {
      dom_names <- as.character(domain)
    }
  }

  n <- 0
  adj_list <- NULL
  dist_mat <- NULL

  # 2. Extract spatial adjacency or coordinates based on data type
  if (inherits(data, "sf")) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg sf} is required to process {.cls sf} objects.")
    }

    n <- nrow(data)
    geom_types <- unique(as.character(sf::st_geometry_type(data)))
    is_polygon <- any(geom_types %in% c("POLYGON", "MULTIPOLYGON"))

    if (method %in% c("queen", "rook")) {
      if (!is_polygon) {
        cli::cli_abort("Method {.val {method}} requires polygon geometry, but got {.val {geom_types}}.")
      }

      if (method == "queen") {
        adj_raw <- sf::st_touches(data, sparse = TRUE)
        adj_list <- lapply(seq_len(n), function(i) setdiff(adj_raw[[i]], i))
      } else {
        if (requireNamespace("spdep", quietly = TRUE)) {
          nb_rook <- spdep::poly2nb(data, queen = FALSE)
          adj_list <- lapply(seq_len(n), function(i) {
            nbrs <- nb_rook[[i]]
            if (length(nbrs) == 1 && nbrs[1] == 0) integer(0) else setdiff(as.integer(nbrs), i)
          })
        } else {
          rel_mat <- sf::st_relate(data, data, pattern = "F***1****", sparse = TRUE)
          adj_list <- lapply(seq_len(n), function(i) setdiff(rel_mat[[i]], i))
        }
      }

      suppressWarnings({
        centr <- sf::st_centroid(sf::st_geometry(data))
      })
      dist_mat <- as.matrix(sf::st_distance(centr, centr))
      storage.mode(dist_mat) <- "double"
      # ponytail: strip units if present (sf >=1.0 returns units matrix)
      if (inherits(dist_mat[1,1], "units") || inherits(dist_mat, "units")) {
        dist_mat <- units::drop_units(dist_mat)
        storage.mode(dist_mat) <- "double"
      }
      # fallback if storage.mode stripped failed
      if (is.character(dist_mat)) storage.mode(dist_mat) <- "double"

    } else {
      suppressWarnings({
        centr <- if (is_polygon) sf::st_centroid(sf::st_geometry(data)) else sf::st_geometry(data)
      })
      dist_mat <- as.matrix(sf::st_distance(centr, centr))
      storage.mode(dist_mat) <- "double"
      if (inherits(dist_mat[1,1], "units") || inherits(dist_mat, "units")) {
        dist_mat <- units::drop_units(dist_mat)
        storage.mode(dist_mat) <- "double"
      }
    }

  } else if (inherits(data, "listw")) {
    if (requireNamespace("spdep", quietly = TRUE)) {
      W_mat <- spdep::listw2mat(data)
      n <- nrow(W_mat)
      adj_list <- lapply(seq_len(n), function(i) which(W_mat[i, ] > 0))
    } else {
      cli::cli_abort("Package {.pkg spdep} is required to process {.cls listw} objects.")
    }

  } else if (inherits(data, "nb")) {
    n <- length(data)
    adj_list <- lapply(seq_len(n), function(i) {
      nbrs <- data[[i]]
      if (length(nbrs) == 1 && nbrs[1] == 0) integer(0) else setdiff(as.integer(nbrs), i)
    })

  } else if (is.matrix(data) || is.data.frame(data)) {
    df_coords <- as.data.frame(data)
    coord_cols <- intersect(tolower(names(df_coords)), c("x", "y", "lon", "lat", "longitude", "latitude"))
    if (length(coord_cols) >= 2) {
      coords_mat <- as.matrix(df_coords[, coord_cols[1:2]])
    } else if (ncol(df_coords) >= 2) {
      coords_mat <- as.matrix(df_coords[, 1:2])
    } else {
      cli::cli_abort("Coordinate data must have at least 2 numeric columns.")
    }

    n <- nrow(coords_mat)
    if (method %in% c("queen", "rook")) {
      cli::cli_alert_warning("Contiguity methods (Queen/Rook) require polygon geometry. Automatically falling back to k-NN (k = {k}).")
      method <- "knn"
    }

    dist_mat <- as.matrix(stats::dist(coords_mat))
    storage.mode(dist_mat) <- "double"

  } else {
    cli::cli_abort("Unsupported input type for {.arg data}. Expected an {.cls sf} object, coordinate matrix, or {.cls nb}/{.cls listw} object.")
  }

  if (is.null(dom_names)) {
    if (!is.null(rownames(data)) && !all(rownames(data) == as.character(seq_len(n)))) {
      dom_names <- rownames(data)
    } else {
      dom_names <- as.character(seq_len(n))
    }
  }

  # 3. Handle Distance-Based Methods (ponytail: C++ for n>500, else stay in R — no new dep)
  if (method == "knn") {
    k_eff <- min(max(1L, as.integer(k)), n - 1L)
    if (n > 500) {
      A_dir <- .knn_adj_cpp(dist_mat, k_eff)
      adj_list <- lapply(seq_len(n), function(i) which(A_dir[i, ] > 0.5))
    } else {
      adj_list <- vector("list", n)
      for (i in seq_len(n)) {
        dists_i <- dist_mat[i, ]
        dists_i[i] <- Inf
        nn_idx <- order(dists_i)[seq_len(k_eff)]
        adj_list[[i]] <- nn_idx
      }
    }

  } else if (method == "distance") {
    if (is.null(d_max)) {
      min_dists <- apply(dist_mat + diag(Inf, n), 1, min)
      d_max <- max(min_dists) * 1.001
      cli::cli_alert_info("Distance threshold {.code d_max} set automatically to {.val {round(d_max, 2)}} to guarantee connectivity.")
    }
    if (n > 500) {
      A_band <- .dist_adj_cpp(dist_mat, d_max)
      adj_list <- lapply(seq_len(n), function(i) which(A_band[i, ] > 0.5))
    } else {
      adj_list <- vector("list", n)
      for (i in seq_len(n)) {
        dists_i <- dist_mat[i, ]
        nbrs <- which(dists_i <= d_max & dists_i > 0)
        adj_list[[i]] <- setdiff(nbrs, i)
      }
    }

  } else if (method == "idw") {
    if (n > 500) {
      W_raw <- .idw_mat_cpp(dist_mat, alpha)
    } else {
      W_raw <- 1 / (dist_mat ^ alpha)
      diag(W_raw) <- 0
      W_raw[!is.finite(W_raw)] <- 0
    }
  }

  # 4. Handle Isolated Islands (Zero Neighbors) if Contiguity / Distance
  n_islands <- 0
  island_idx <- integer(0)

  if (method != "idw") {
    degrees <- lengths(adj_list)
    island_idx <- which(degrees == 0)
    n_islands <- length(island_idx)

    if (n_islands > 0) {
      if (isTRUE(island_fallback) && !is.null(dist_mat)) {
        cli::cli_alert_info(
          "{n_islands} island/isolated domain(s) detected with zero neighbors ({paste(dom_names[island_idx], collapse = ', ')}). Automatically connecting to nearest neighbor(s) via centroid distance."
        )
        for (idx in island_idx) {
          dists_i <- dist_mat[idx, ]
          dists_i[idx] <- Inf
          nearest_k <- order(dists_i)[seq_len(min(k_fallback, n - 1L))]
          adj_list[[idx]] <- nearest_k
          for (k_neighbor in nearest_k) {
            adj_list[[k_neighbor]] <- sort(unique(c(adj_list[[k_neighbor]], idx)))
          }
        }
      } else {
        cli::cli_alert_warning(
          "{n_islands} domain(s) have zero neighbors in the spatial matrix (isolated nodes). Consider setting {.code island_fallback = TRUE} or using {.code method = 'knn'}."
        )
      }
    }

    W_raw <- matrix(0, nrow = n, ncol = n)
    for (i in seq_len(n)) {
      nbrs <- adj_list[[i]]
      if (length(nbrs) > 0) {
        W_raw[i, nbrs] <- 1
      }
    }
  }

  diag(W_raw) <- 0

  # 5. Apply Normalization Style
  W_out <- matrix(0, nrow = n, ncol = n)

  if (style == "W") {
    rs <- rowSums(W_raw)
    nonzero_rows <- rs > 0
    W_out[nonzero_rows, ] <- W_raw[nonzero_rows, ] / rs[nonzero_rows]

  } else if (style == "B") {
    W_out <- (W_raw > 0) * 1

  } else if (style == "C") {
    total_w <- sum(W_raw)
    if (total_w > 0) {
      W_out <- (n / total_w) * W_raw
    }

  } else if (style == "U") {
    W_out <- W_raw / n

  } else if (style == "minmax") {
    # Guard dense eigen O(n^3): defer for large n (ponytail: use power iteration / RSpectra when n>500)
    if (n > 500) {
      cli::cli_alert_warning("Style {.val minmax} with {.val n}={n} would require dense eigen O(n^3); falling back to row-standardized {.val W} and capping spectral radius via row-sum bound.")
      rs <- rowSums(abs(W_raw))
      max_rs <- max(rs)
      if (max_rs > 0) W_out <- W_raw / max_rs else W_out <- W_raw
    } else {
      ev <- eigen(W_raw, only.values = TRUE)$values
      max_ev <- max(abs(Re(ev)))
      if (max_ev > 0) {
        W_out <- W_raw / max_ev
      } else {
        W_out <- W_raw
      }
    }
  }

  rownames(W_out) <- dom_names
  colnames(W_out) <- dom_names

  attr(W_out, "method") <- method
  attr(W_out, "style") <- style
  attr(W_out, "k") <- if (method == "knn") k else NULL
  attr(W_out, "d_max") <- if (method == "distance") d_max else NULL
  attr(W_out, "alpha") <- if (method == "idw") alpha else NULL
  attr(W_out, "n_domains") <- n
  attr(W_out, "n_islands") <- n_islands
  attr(W_out, "island_domains") <- if (n_islands > 0) dom_names[island_idx] else character(0)
  attr(W_out, "call") <- call

  class(W_out) <- c("fastsae_weights", "matrix")
  return(W_out)
}

#' @export
print.fastsae_weights <- function(x, ...) {
  n <- attr(x, "n_domains") %||% nrow(x)
  meth <- attr(x, "method") %||% "unknown"
  sty <- attr(x, "style") %||% "unknown"
  n_isl <- attr(x, "n_islands") %||% 0

  cli::cli_h1("Spatial Weights Matrix (fastsae)")
  cli::cli_bullets(c(
    "*" = "Dimensions: {.val {n}} x {.val {n}} domains",
    "*" = "Construction Method: {.field {toupper(meth)}}",
    "*" = "Standardization Style: {.val {sty}} ({.emph {switch(sty, 'W'='Row-standardized', 'B'='Binary', 'C'='Globally standardized', 'U'='Equal weights', 'minmax'='Spectral min-max', sty)}})",
    "*" = "Sparsity: {.val {round((1 - sum(x != 0) / (n * n)) * 100, 1)}\\%} zeros"
  ))

  if (n_isl > 0) {
    isl_names <- attr(x, "island_domains")
    cli::cli_alert_info("Island domains connected: {.val {length(isl_names)}} ({paste(utils::head(isl_names, 5), collapse = ', ')}{if (length(isl_names) > 5) '...' else ''})")
  }

  cat("\nMatrix preview (top-left 5x5):\n")
  p_rows <- seq_len(min(5, n))
  print(round(x[p_rows, p_rows, drop = FALSE], 4))
  if (n > 5) {
    cat(sprintf("... and %d more rows/columns.\n", n - 5))
  }
  invisible(x)
}

#' @export
summary.fastsae_weights <- function(object, ...) {
  n <- nrow(object)
  degree <- rowSums(object > 0)

  cli::cli_h2("Spatial Neighbor Distribution Summary")
  cat("Neighbors per domain (non-zero entries):\n")
  print(summary(degree))

  isolated <- sum(degree == 0)
  if (isolated > 0) {
    cli::cli_alert_warning("{isolated} domain(s) have zero neighbors.")
  } else {
    cli::cli_alert_success("All domains are connected (minimum {min(degree)} neighbor(s)).")
  }
  invisible(summary(degree))
}

#' @export
plot.fastsae_weights <- function(x, ...) {
  degree <- rowSums(x > 0)
  df_deg <- data.frame(Neighbors = degree)

  p <- ggplot2::ggplot(df_deg, ggplot2::aes(x = .data$Neighbors)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "#2C3E50", color = "white", alpha = 0.85) +
    ggplot2::geom_vline(xintercept = mean(degree), color = "#E74C3C", linetype = "dashed", linewidth = 1) +
    ggplot2::labs(
      title = "Distribution of Spatial Neighbors per Domain",
      subtitle = sprintf("Mean: %.1f neighbors (Method: %s, Style: %s)", mean(degree), toupper(attr(x, "method")), attr(x, "style")),
      x = "Number of Neighbors",
      y = "Domain Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))

  print(p)
  invisible(p)
}
