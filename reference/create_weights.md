# Create Spatial Proximity and Weights Matrix

Constructs spatial proximity, contiguity, and weight matrices
(\\\mathbf{W}\\) from spatial data (`sf` polygon or point objects,
coordinate matrices, or neighbor lists). Supports Queen and Rook
contiguity, \\k\\-nearest neighbors (\\k\\-NN), distance bands, and
inverse distance weighting (IDW), with automatic fallback handling for
isolated islands and disconnected domains.

## Usage

``` r
create_weights(
  data,
  method = c("queen", "rook", "knn", "distance", "idw"),
  style = c("W", "B", "C", "U", "minmax"),
  k = 4,
  d_max = NULL,
  alpha = 1,
  island_fallback = TRUE,
  k_fallback = 1,
  domain = NULL
)
```

## Arguments

- data:

  An `sf` polygon or point object, a numeric matrix or data frame of 2D
  coordinates (with columns `x`/`y` or longitude/latitude), or an
  `spdep` `nb` or `listw` object.

- method:

  Character string specifying the spatial weight method:

  - `"queen"`: Queen contiguity (sharing at least one vertex or boundary
    edge).

  - `"rook"`: Rook contiguity (sharing a common linear boundary edge).

  - `"knn"`: \\k\\-nearest neighbors based on centroid/point distances.

  - `"distance"`: Distance band where domains within `d_max` are
    neighbors.

  - `"idw"`: Inverse distance weighting \\w\_{ij} = 1 /
    d\_{ij}^\alpha\\.

- style:

  Character string specifying the normalization style:

  - `"W"`: Row-standardized (row sums equal 1, default for SAR/SFH
    models).

  - `"B"`: Binary adjacency (0/1 values, standard for INLA BYM2/Besag
    models).

  - `"C"`: Globally standardized.

  - `"U"`: Equal weighting (\\1 / m\\).

  - `"minmax"`: Min-max spectral normalization (divided by largest
    eigenvalue).

- k:

  Integer. Number of nearest neighbors when `method = "knn"` (default is
  4).

- d_max:

  Numeric. Distance threshold when `method = "distance"`. If `NULL`, it
  is automatically set to the minimum distance ensuring no disconnected
  domains.

- alpha:

  Numeric. Distance decay power for `method = "idw"` (default is 1).

- island_fallback:

  Logical. If `TRUE` (default), disconnected domains or islands with
  zero contiguous neighbors are automatically connected to their nearest
  neighbor(s) based on centroid distance.

- k_fallback:

  Integer. Number of nearest neighbors to connect for isolated islands
  when `island_fallback = TRUE` (default is 1).

- domain:

  Optional vector, column name, or one-sided formula referencing domain
  IDs. If provided, used as row and column names for the resulting
  matrix.

## Value

A square numeric matrix of class `c("fastsae_weights", "matrix")` of
dimension \\m \times m\\. Attributes include `method`, `style`,
`n_domains`, `n_islands`, and `call`.

## Details

When working with real administrative boundaries (such as Indonesian
provinces, regencies, or municipalities), island regencies (e.g.
Kepulauan Seribu, Natuna, Sabang, Mentawai, Sangihe) lack contiguous
land borders. If `island_fallback = TRUE`, `create_weights` detects
these zero-neighbor domains and connects each to its closest
`k_fallback` neighbor(s) based on centroid distances, ensuring a
connected spatial graph without computational singularities.

## Examples

``` r
if (requireNamespace("sf", quietly = TRUE) && requireNamespace("spdep", quietly = TRUE)) {
  # North Carolina counties, an sf polygon data set shipped with 'sf'
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  W_queen <- create_weights(nc, method = "queen", style = "W")
  W_knn   <- create_weights(nc, method = "knn", k = 5, style = "W")
  W_bin   <- create_weights(nc, method = "queen", style = "B")
  W_queen
}
#> 
#> ── Spatial Weights Matrix (fastsae) ────────────────────────────────────────────
#> • Dimensions: 100 x 100 domains
#> • Construction Method: QUEEN
#> • Standardization Style: "W" (Row-standardized)
#> • Sparsity: 95.1% zeros
#> 
#> Matrix preview (top-left 5x5):
#>        1      2      3 4 5
#> 1 0.0000 0.3333 0.0000 0 0
#> 2 0.3333 0.0000 0.3333 0 0
#> 3 0.0000 0.2000 0.0000 0 0
#> 4 0.0000 0.0000 0.0000 0 0
#> 5 0.0000 0.0000 0.0000 0 0
#> ... and 95 more rows/columns.

# Create spatial weights from 2D coordinates
coords <- matrix(runif(40), ncol = 2)
W_coord <- create_weights(coords, method = "knn", k = 3, style = "W")
dim(W_coord)
#> [1] 20 20
```
