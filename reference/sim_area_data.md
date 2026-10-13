# Simulate Multi-Distribution Small Area Data with Spatial Autocorrelation

Generates synthetic area-level datasets for small area estimation,
including continuous (Gaussian, Gamma), count (Poisson, Negative
Binomial), and proportion/rate responses (Binomial, Beta), driven by
shared covariates and structured spatial random effects (BYM2 / SAR).

## Usage

``` r
sim_area_data(
  D = 42,
  W = NULL,
  spatial_type = c("knn", "grid", "ring"),
  rho = 0.5,
  phi = 0.6,
  sigma_u = 0.5,
  beta = c(1, 0.8, -0.5),
  n_unsampled = 6,
  seed = NULL
)
```

## Arguments

- D:

  Integer. Number of domains (areas). Default is 42.

- W:

  Optional spatial proximity or adjacency matrix of dimension \\D \times
  D\\. If `NULL`, a spatial matrix is automatically simulated using
  [`sim_spatial_weights`](https://ridsonap.github.io/fastsae/reference/sim_spatial_weights.md).

- spatial_type:

  Character. Topology for simulated `W` if `W = NULL`: `"knn"`,
  `"grid"`, or `"ring"`. Default is `"knn"`.

- rho:

  Numeric. Spatial autocorrelation autoregressive parameter (\\\|\rho\|
  \< 1\\). Default is 0.5.

- phi:

  Numeric. Proportion of spatial marginal variance relative to total
  area variance (BYM2 parameter, \\0 \le \phi \le 1\\). Default is 0.6.

- sigma_u:

  Numeric. Overall standard deviation of domain random effects. Default
  is 0.5.

- beta:

  Numeric vector of length 3 giving true regression coefficients
  \\(\beta_0, \beta_1, \beta_2)\\. Default is `c(1.0, 0.8, -0.5)`.

- n_unsampled:

  Integer. Number of domains to mark as unsampled (where response
  variables are set to `NA`). Default is 6.

- seed:

  Optional integer. Random seed for reproducibility.

## Value

An object of class `c("fastsae_sim_data", "list")` containing:

- data:

  A data frame with columns: `domain`, `x1`, `x2`, `y_gaussian`,
  `vardir`, `y_poisson`, `exposure`, `y_binomial`, `trials`, `y_beta`,
  `y_nbinomial`, `y_gamma`, `x_coord`, `y_coord`.

- W:

  Binary adjacency matrix \\D \times D\\.

- W_std:

  Row-standardized proximity matrix \\D \times D\\.

- coords:

  Coordinates of domain centroids.

- u_spatial:

  True structured spatial random effect vector.

- u_iid:

  True unstructured random effect vector.

- u_total:

  True total random effect vector.

- phi:

  True spatial variance fraction.

- rho:

  True spatial autoregressive parameter.

## Examples

``` r
# 1. Simulate dataset with default KNN spatial structure
sim <- sim_area_data(D = 40, spatial_type = "knn", seed = 123)
print(sim)
#> ======================================================
#>   fastsae Simulated Multi-Distribution Area Data
#> ======================================================
#> Domains (D): 40 (Sampled: 34, Unsampled: 6)
#> Spatial parameters: rho = 0.50, phi = 0.60
#> Responses generated:
#>   - Gaussian (FH):     y_gaussian (vardir)
#>   - Poisson:           y_poisson (exposure)
#>   - Binomial:          y_binomial (trials)
#>   - Beta:              y_beta
#>   - Negative Binomial: y_nbinomial
#>   - Gamma:             y_gamma
#> Spatial matrices:      40 x 40 (W binary, W_std row-standardized)
#> ======================================================
head(sim$data)
#>   domain     x1     x2 y_gaussian  vardir y_poisson exposure y_binomial trials
#> 1      1 1.3053 2.5115     1.4564 0.05726       137      149         92    205
#> 2      2 1.7921 1.7695     1.7618 0.06314       385      340         72    125
#> 3      3 0.7346 3.2499    -0.3993 0.14761       105      178         17     58
#> 4      4 4.1690 1.8736     2.8760 0.07697       771      300        108    123
#> 5      5 3.2080 1.7772     1.7168 0.08360       611      378         72    105
#> 6      6 0.8769 2.6684     0.0746 0.13407        79      117         64    220
#>    y_beta y_nbinomial y_gamma truth_gaussian x_coord y_coord
#> 1 0.48992          32  2.0911         0.9302  0.2876  0.1428
#> 2 0.47039          41  2.0245         1.4321  0.7883  0.4145
#> 3 0.18907           5  0.5666        -0.2974  0.4090  0.4137
#> 4 0.93489           9  3.6981         2.7253  0.8830  0.3688
#> 5 0.72223          44  1.2174         1.8366  0.9405  0.1524
#> 6 0.29830           8  0.5810         0.2490  0.0456  0.1388

# 2. Using an existing spatial matrix
my_W <- sim_spatial_weights(D = 30, type = "grid")
sim2 <- sim_area_data(W = my_W, seed = 456)
```
