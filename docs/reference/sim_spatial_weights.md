# Simulate Spatial Proximity and Adjacency Matrices

Generates synthetic spatial proximity or adjacency matrices for small
area estimation models. Supports k-nearest neighbors (KNN) on random 2D
coordinates, regular grid contiguity (lattice), and 1D circular rings.

## Usage

``` r
sim_spatial_weights(
  D = 40,
  type = c("knn", "grid", "ring"),
  style = c("B", "W"),
  k = 4,
  coords = NULL,
  seed = NULL
)
```

## Arguments

- D:

  Integer. Number of small areas (domains). Default is 40.

- type:

  Character. Topology type: `"knn"` (k-nearest neighbors on 2D
  coordinates), `"grid"` (Rook contiguity on a regular 2D lattice), or
  `"ring"` (1D circular chain). Default is `"knn"`.

- style:

  Character. Weight style: `"B"` for symmetric binary adjacency (0/1)
  with zero diagonal (standard for INLA graph models such as BYM2, BYM,
  Besag), or `"W"` for row-standardized weights (sum of each row equals
  1, standard for SAR/SLM models). Default is `"B"`.

- k:

  Integer. Number of nearest neighbors per domain when `type = "knn"`.
  Default is 4. Must be between 1 and `D - 1`.

- coords:

  Optional numeric matrix of dimension `c(D, 2)` containing 2D spatial
  coordinates. If provided and `type = "knn"`, distances are computed
  from these coordinates.

- seed:

  Optional integer. Random seed for reproducible coordinate generation.

## Value

A square numeric matrix of dimension \\D \times D\\ with row and column
names set to domain identifiers (`"1"`, `"2"`, ..., `"D"`). Attributes
include `"coords"` (2D coordinates) and `"type"`.

## Details

When `type = "knn"`, coordinates are sampled from \\\text{Uniform}(0,
1)^2\\ if not provided. Each domain is connected to its \\k\\ nearest
neighbors. To ensure a symmetric adjacency graph required by GMRF
spatial priors (INLA), edges are made mutual: domain \\i\\ and \\j\\ are
connected if \\j\\ is among the \\k\\ nearest neighbors of \\i\\ or vice
versa.

When `type = "grid"`, domains are placed on a regular \\r \times c\\
lattice where \\r = \lfloor\sqrt{D}\rfloor\\ and \\c = \lceil D / r
\rceil\\. If \\r \times c \> D\\, the first \\D\\ lattice points are
retained. Rook contiguity (sharing a common edge) is used.

When `type = "ring"`, domain \\i\\ is connected to \\i-1\\ and \\i+1\\
in a closed ring.

## Examples

``` r
# 1. Binary adjacency matrix for 20 domains using KNN
W_bin <- sim_spatial_weights(D = 20, type = "knn", k = 3, seed = 123)
dim(W_bin)
#> [1] 20 20
table(W_bin)
#> W_bin
#>   0   1 
#> 322  78 

# 2. Row-standardized matrix on a regular grid
W_std <- sim_spatial_weights(D = 25, type = "grid", style = "W")
rowSums(W_std)
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 
#>  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1 
```
