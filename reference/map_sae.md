# Unified Spatial Choropleth Mapping for Small Area Estimates

Provides a comprehensive, single-gateway ("satu pintu") spatial
choropleth mapping interface for Small Area Estimation models fitted
with fastsae (e.g. `eblup_fh`, `eblup_sfh`, `eblup_stfh`, `eblup_bhf`,
`hb_area`, and `benchmark_sae`).

Intelligently handles:

- **Single Model Visualizations**: Point estimates (`"estimate"`),
  Relative Standard Error (`"rse"`), and 3-color official statistics
  reliability classification (`"reliability"`: green for reliable \<
  20%, amber for caution 20-30%, red for unreliable \\\ge 30\\\\).

- **Direct vs Model Comparisons**: Side-by-side facet maps contrasting
  erratic direct survey estimates against smoothed small area model
  estimates (`"comparison"`).

- **Two-Model Comparisons**: Side-by-side maps or spatial difference
  maps (`"difference"`) contrasting two competing models (e.g. Classical
  FH vs Spatial FH, or FastSAE vs Stan).

- **Benchmarked Calibrations**: Visualizing pre- vs post-benchmarked
  calibrations.

- **Smart Key Matching**: Automatically detects the common domain
  identifier between the model and `sf_geom`, with optional manual
  override via `key`.

## Usage

``` r
map_sae(object, ...)

# S3 method for class 'fastsae'
map_sae(
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
)

# S3 method for class 'fastsae_benchmark'
map_sae(
  object,
  sf_geom = NULL,
  key = NULL,
  type = c("comparison", "estimate", "difference"),
  indicator = NULL,
  palette = NULL,
  title = NULL,
  subtitle = NULL,
  ...
)

# S3 method for class 'list'
map_sae(
  object,
  sf_geom = NULL,
  key = NULL,
  type = c("comparison", "difference"),
  indicator = NULL,
  palette = NULL,
  title = NULL,
  subtitle = NULL,
  ...
)

# Default S3 method
map_sae(
  object,
  sf_geom = NULL,
  key = NULL,
  type = c("estimate", "rse", "reliability"),
  indicator = NULL,
  palette = NULL,
  title = NULL,
  subtitle = NULL,
  ...
)
```

## Arguments

- object:

  A fitted `fastsae` model object, a `fastsae_benchmark` object, a named
  list of `fastsae` objects, or a data frame containing domain
  estimates.

- ...:

  Additional arguments passed to `ggplot2` layers.

- model2:

  Optional second `fastsae` model object for direct model-to-model
  comparison.

- sf_geom:

  Spatial polygon geometry of class `sf`. Optional if the model was
  fitted on data that already inherits from `sf`.

- key:

  Optional character string or named vector specifying the column(s)
  used to match model domain identifiers with the spatial geometry in
  `sf_geom`.

  - If `NULL` (default), automatically detects the matching key by
    evaluating common column names or the highest set intersection
    overlap between domain IDs and `sf_geom` columns.

  - If a single string (e.g. `key = "kd_kab"`), specifies the column in
    `sf_geom` to match with model domain IDs.

  - If a named string (e.g. `key = c("domain" = "kd_kab")`), explicitly
    pairs the model's domain column with the spatial column in
    `sf_geom`.

- type:

  Character string specifying the visualization type:

  - `"estimate"`: Choropleth map of small area point estimates (EBP or
    EBLUP).

  - `"rse"`: Choropleth map of Relative Standard Errors / CV (%).

  - `"reliability"`: Official statistics 3-color traffic light map (\<
    20% Reliable, 20-30% Use with Caution, \\\ge 30\\\\ Unreliable).

  - `"comparison"`: Side-by-side facet choropleth map (Direct vs Model,
    or Model 1 vs Model 2).

  - `"difference"`: Spatial difference map (\\\hat{\theta}^{(2)} -
    \hat{\theta}^{(1)}\\ or \\\hat{\theta}^{\text{Model}} -
    \hat{\theta}^{\text{Direct}}\\) with divergent color palette.

- indicator:

  Optional alias for `type`. If supplied, overrides `type`.

- palette:

  Character string specifying a color palette. Default is automatically
  selected based on `type`: `"viridis"` for estimates, `"magma"` for
  RSE, official traffic-light colors for reliability, and `"RdBu"` /
  `"PuOr"` for differences.

- thresholds:

  Numeric vector of length 2 defining the RSE (%) thresholds for
  reliability classification. Default is `c(20, 30)` according to
  official statistical guidelines (BPS / Eurostat).

- facet_scales:

  Character string passed to `ggplot2::facet_wrap(scales = ...)`.
  Default is `"fixed"`.

- title:

  Optional title for the plot. If `NULL`, a descriptive title is
  generated automatically.

- subtitle:

  Optional subtitle for the plot.

## Value

A `ggplot` object containing the thematic choropleth map.

## Examples

``` r
library(fastsae)
data(mys)

# Create a synthetic spatial grid for demonstration
if (requireNamespace("sf", quietly = TRUE)) {
  grid_sf <- sf::st_make_grid(
    sf::st_polygon(list(matrix(c(0,0, 6,0, 6,7, 0,7, 0,0), ncol = 2, byrow = TRUE))),
    cellsize = c(1, 1),
    what = "polygons"
  )[seq_len(nrow(mys))]
  mys_sf <- sf::st_sf(area = mys$area, geometry = grid_sf)

  # 1. Fit Fay-Herriot model
  fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)

  # 2. Single model point estimate map (auto-detected key)
  p1 <- map_sae(fit_fh, sf_geom = mys_sf)

  # 3. Explicit key matching
  p2 <- map_sae(fit_fh, sf_geom = mys_sf, key = "area")

  # 4. Reliability classification map (BPS standard <20%, 20-30%, >=30%)
  p3 <- map_sae(fit_fh, sf_geom = mys_sf, type = "reliability")

  # 5. Side-by-side comparison map: Direct Survey vs Model SAE
  p4 <- map_sae(fit_fh, sf_geom = mys_sf, type = "comparison")
}
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_fh(formula = y ~ x1 + x2, vardir = "vardir", data = mys)
#> 
#> ✔ Convergence: Yes (in 7 iterations)
#> Model: Fay-Herriot (Area-level)
#> Method: eblup
#> Random effect variance (sigma2_u): 2.569182 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.0689476  0.7631917  4.0212018 0.0001
#> x1          -0.0041073  0.0090768 -0.4525074 0.6509
#> x2           0.0861173  0.0309576  2.7817847 0.0054
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    eblup    vardir random_effect       mse       rse
#> 1      1 8.359527 7.594807 0.6618838    2.96835204 0.5602338  9.855256
#> 2      2 7.599650 6.777056 0.8374691    2.52354833 0.6766501 12.137828
#> 3      3 5.514137 5.247512 0.8822257    0.77645304 0.7048881 15.999508
#> 4      4 3.869326 4.247072 0.6581716   -1.47453645 0.5556741 17.551750
#> 5      5 6.305063 6.353632 1.2788021   -0.09757891 0.9308379 15.185005
#> 6      6 3.926807 4.090117 0.3878004   -1.08192989 0.3497100 14.458337
#> ... and 36 more rows.
#> 
```
