#' @title Simulate Spatial Proximity and Adjacency Matrices
#' @description Generates synthetic spatial proximity or adjacency matrices for
#'   small area estimation models. Supports k-nearest neighbors (KNN) on random
#'   2D coordinates, regular grid contiguity (lattice), and 1D circular rings.
#'
#' @param D Integer. Number of small areas (domains). Default is 40.
#' @param type Character. Topology type: \code{"knn"} (k-nearest neighbors on 2D coordinates),
#'   \code{"grid"} (Rook contiguity on a regular 2D lattice), or \code{"ring"} (1D circular chain).
#'   Default is \code{"knn"}.
#' @param style Character. Weight style: \code{"B"} for symmetric binary adjacency (0/1) with zero
#'   diagonal (standard for INLA graph models such as BYM2, BYM, Besag), or \code{"W"} for
#'   row-standardized weights (sum of each row equals 1, standard for SAR/SLM models).
#'   Default is \code{"B"}.
#' @param k Integer. Number of nearest neighbors per domain when \code{type = "knn"}.
#'   Default is 4. Must be between 1 and \code{D - 1}.
#' @param coords Optional numeric matrix of dimension \code{c(D, 2)} containing 2D spatial coordinates.
#'   If provided and \code{type = "knn"}, distances are computed from these coordinates.
#' @param seed Optional integer. Random seed for reproducible coordinate generation.
#'
#' @return A square numeric matrix of dimension \eqn{D \times D} with row and column names
#'   set to domain identifiers (\code{"1"}, \code{"2"}, ..., \code{"D"}).
#'   Attributes include \code{"coords"} (2D coordinates) and \code{"type"}.
#'
#' @details
#' When \code{type = "knn"}, coordinates are sampled from \eqn{\text{Uniform}(0, 1)^2} if not provided.
#' Each domain is connected to its \eqn{k} nearest neighbors. To ensure a symmetric adjacency graph
#' required by GMRF spatial priors (INLA), edges are made mutual: domain \eqn{i} and \eqn{j} are
#' connected if \eqn{j} is among the \eqn{k} nearest neighbors of \eqn{i} or vice versa.
#'
#' When \code{type = "grid"}, domains are placed on a regular \eqn{r \times c} lattice where
#' \eqn{r = \lfloor\sqrt{D}\rfloor} and \eqn{c = \lceil D / r \rceil}. If \eqn{r \times c > D},
#' the first \eqn{D} lattice points are retained. Rook contiguity (sharing a common edge) is used.
#'
#' When \code{type = "ring"}, domain \eqn{i} is connected to \eqn{i-1} and \eqn{i+1} in a closed ring.
#'
#' @examples
#' # 1. Binary adjacency matrix for 20 domains using KNN
#' W_bin <- sim_spatial_weights(D = 20, type = "knn", k = 3, seed = 123)
#' dim(W_bin)
#' table(W_bin)
#'
#' # 2. Row-standardized matrix on a regular grid
#' W_std <- sim_spatial_weights(D = 25, type = "grid", style = "W")
#' rowSums(W_std)
#'
#' @export
sim_spatial_weights <- function(D = 40,
                                type = c("knn", "grid", "ring"),
                                style = c("B", "W"),
                                k = 4,
                                coords = NULL,
                                seed = NULL) {
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  D <- as.integer(D)
  if (is.na(D) || D < 2L) {
    stop("'D' must be an integer >= 2.", call. = FALSE)
  }

  type <- match.arg(type)
  style <- match.arg(style)

  domain_ids <- as.character(seq_len(D))

  if (type == "knn") {
    if (is.null(coords)) {
      coords <- matrix(stats::runif(2L * D), ncol = 2L)
    } else {
      coords <- as.matrix(coords)
      if (nrow(coords) != D || ncol(coords) < 2L) {
        stop(sprintf("'coords' must be a matrix or data.frame with %d rows and >= 2 columns.", D),
             call. = FALSE)
      }
      coords <- coords[, 1:2, drop = FALSE]
    }
    colnames(coords) <- c("x", "y")
    rownames(coords) <- domain_ids

    k_eff <- min(max(1L, as.integer(k)), D - 1L)
    dmat <- as.matrix(stats::dist(coords))

    A <- matrix(0, nrow = D, ncol = D)
    for (i in seq_len(D)) {
      # Find k nearest neighbors (excluding self at distance 0)
      order_idx <- order(dmat[i, ])
      nn <- order_idx[order_idx != i][seq_len(k_eff)]
      A[i, nn] <- 1
    }
    # Mutual connectivity to ensure symmetry
    A <- pmax(A, t(A))
    diag(A) <- 0

  } else if (type == "grid") {
    r <- max(1L, as.integer(floor(sqrt(D))))
    c <- as.integer(ceiling(D / r))

    grid_df <- expand.grid(x = seq_len(c), y = seq_len(r))
    coords <- as.matrix(grid_df[seq_len(D), c("x", "y")])
    rownames(coords) <- domain_ids

    # Rook contiguity: Manhattan distance == 1
    dmat_m <- as.matrix(stats::dist(coords, method = "manhattan"))
    A <- matrix(0, nrow = D, ncol = D)
    A[dmat_m == 1] <- 1
    diag(A) <- 0

    # Ensure no isolated domains if D is not a multiple of r
    rs <- rowSums(A)
    isolated <- which(rs == 0)
    if (length(isolated) > 0) {
      dmat_e <- as.matrix(stats::dist(coords, method = "euclidean"))
      for (iso in isolated) {
        order_idx <- order(dmat_e[iso, ])
        nearest <- order_idx[order_idx != iso][1]
        A[iso, nearest] <- 1
        A[nearest, iso] <- 1
      }
    }

  } else if (type == "ring") {
    angles <- seq(0, 2 * pi, length.out = D + 1L)[seq_len(D)]
    coords <- cbind(x = cos(angles), y = sin(angles))
    rownames(coords) <- domain_ids

    A <- matrix(0, nrow = D, ncol = D)
    for (i in seq_len(D)) {
      left <- if (i == 1L) D else i - 1L
      right <- if (i == D) 1L else i + 1L
      A[i, left] <- 1
      A[i, right] <- 1
    }
    diag(A) <- 0
  }

  if (style == "W") {
    rs <- rowSums(A)
    rs[rs == 0] <- 1
    W <- A / rs
  } else {
    W <- A
  }

  dimnames(W) <- list(domain_ids, domain_ids)
  attr(W, "coords") <- coords
  attr(W, "type") <- type
  attr(W, "style") <- style

  return(W)
}
