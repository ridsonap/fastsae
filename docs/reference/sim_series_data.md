# Simulate Spatio-Temporal Multi-Distribution Small Area Data

Generates synthetic spatio-temporal panel datasets for small area
estimation. Combines domain-level spatial autocorrelation (SAR) with
domain-specific first-order autoregressive AR(1) temporal dynamics
across \\T\\ time periods, producing multi-distribution responses
(Gaussian Fay-Herriot, Poisson, Binomial, Beta, Negative Binomial,
Gamma).

## Usage

``` r
sim_series_data(
  D = 42,
  T = 5,
  time_start = 2022,
  W = NULL,
  spatial_type = c("knn", "grid", "ring"),
  rho_s = 0.5,
  rho_t = 0.6,
  sigma_s = 0.4,
  sigma_t = 0.3,
  trend = c("linear", "random_walk", "ar1", "none"),
  trend_slope = 0.05,
  beta = c(1, 0.8, -0.5),
  n_unsampled = 4,
  prop_intermittent = 0.05,
  sort_order = c("domain-major", "time-major"),
  seed = NULL
)
```

## Arguments

- D:

  Integer. Number of domains (areas). Default is 42.

- T:

  Integer. Number of time periods (\\T \ge 2\\). Default is 5.

- time_start:

  Integer. Starting year / time label. Default is 2022.

- W:

  Optional spatial proximity or adjacency matrix of dimension \\D \times
  D\\. If `NULL`, a spatial matrix is automatically simulated using
  [`sim_spatial_weights`](https://ridsonap.github.io/fastsae/reference/sim_spatial_weights.md).

- spatial_type:

  Character. Topology for simulated `W` if `W = NULL`: `"knn"`,
  `"grid"`, or `"ring"`. Default is `"knn"`.

- rho_s:

  Numeric. Domain spatial autoregressive autocorrelation parameter
  (\\\|\rho_s\| \< 1\\, corresponding to `rho1` in
  [`eblup_stfh`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)).
  Default is 0.5.

- rho_t:

  Numeric. Temporal AR(1) autocorrelation parameter (\\\|\rho_t\| \<
  1\\, corresponding to `rho2` in
  [`eblup_stfh`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)).
  Default is 0.6.

- sigma_s:

  Numeric. Standard deviation of domain spatial random effects
  (\\u_1\\). Default is 0.4.

- sigma_t:

  Numeric. Standard deviation of stationary spatio-temporal effects
  (\\u_2\\). Default is 0.3.

- trend:

  Character. Deterministric temporal trend type: `"linear"`,
  `"random_walk"`, `"ar1"`, or `"none"`. Default is `"linear"`.

- trend_slope:

  Numeric. Slope / scale of the temporal trend. Default is 0.05.

- beta:

  Numeric vector of length 3 giving fixed regression coefficients
  \\(\beta_0, \beta_1, \beta_2)\\. Default is `c(1.0, 0.8, -0.5)`.

- n_unsampled:

  Integer. Number of domains that are completely unsampled across all
  time periods. Default is 4.

- prop_intermittent:

  Numeric. Proportion of domain-by-time observations to mark as
  intermittently missing (excluding persistent unsampled domains).
  Default is 0.05.

- sort_order:

  Character. Row sorting order: `"domain-major"` (all time periods for
  domain 1, then domain 2, required by
  [`eblup_stfh`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md))
  or `"time-major"`. Default is `"domain-major"`.

- seed:

  Optional integer. Random seed for reproducibility.

## Value

An object of class `c("fastsae_sim_series", "list")` containing:

- data:

  A data frame of dimension \\(D \cdot T) \times 15\\ with columns:
  `area`, `year`, `x1`, `x2`, `y_gaussian`, `vardir`, `y_poisson`,
  `exposure`, `y_binomial`, `trials`, `y_beta`, `y_nbinomial`,
  `y_gamma`, `x_coord`, `y_coord`.

- W:

  Binary adjacency matrix \\D \times D\\.

- W_std:

  Row-standardized spatial weights matrix \\D \times D\\.

- coords:

  Centroid coordinates of the \\D\\ domains.

- u_spatial:

  True domain spatial random effect vector of length \\D\\.

- u_temporal:

  True spatio-temporal random effect matrix of dimension \\D \times T\\.

- parameters:

  List of true simulation parameters (`rho_s`, `rho_t`, `sigma_s`,
  `sigma_t`, `beta`, `trend`).

## Examples

``` r
# 1. Simulate a 30-domain panel over 4 years
sim_panel <- sim_series_data(D = 30, T = 4, time_start = 2021, seed = 123)
print(sim_panel)
#> =================================================================
#>   fastsae Simulated Spatio-Temporal Multi-Distribution Panel
#> =================================================================
#> Domains (D): 30 | Time periods (T): 4 | Total observations: 120
#> Observations: 99 observed, 21 missing/unsampled (17.5%)
#> Spatio-temporal parameters:
#>   - Spatial correlation (rho_s / rho1):  0.50 (sigma_s = 0.40)
#>   - Temporal AR(1)     (rho_t / rho2):  0.60 (sigma_t = 0.30)
#>   - Temporal trend:                     linear (slope = 0.050)
#>   - Row sort order:                     domain-major
#> Responses generated:
#>   - Gaussian (FH):     y_gaussian (vardir)
#>   - Poisson:           y_poisson (exposure)
#>   - Binomial:          y_binomial (trials)
#>   - Beta:              y_beta
#>   - Negative Binomial: y_nbinomial
#>   - Gamma:             y_gamma
#> Spatial proximity:     30 x 30 (W binary, W_std row-standardized)
#> =================================================================
head(sim_panel$data)
#>   area year     x1     x2 y_gaussian  vardir y_poisson exposure y_binomial
#> 1    1 2021 1.1643 2.2742     1.4259 0.04593       258      225         98
#> 2    1 2022 0.9732 2.7502     0.9198 0.12341       179      227         47
#> 3    1 2023 0.9700 2.6740     1.1845 0.08359       261      316         77
#> 4    1 2024 1.3705 2.3545     1.0411 0.14610       439      416         45
#> 5    2 2021 2.7566 2.7597     1.7132 0.13304       160      161         60
#> 6    2 2022 3.4994 2.9991     2.4992 0.05670       482      276        188
#>   trials  y_beta y_nbinomial y_gamma x_coord y_coord
#> 1    223 0.71535           9  2.0850  0.2876  0.9630
#> 2    125 0.36109          24  1.4626  0.2876  0.9630
#> 3    201 0.49363          34  1.8408  0.2876  0.9630
#> 4     87 0.48010          13  0.9031  0.2876  0.9630
#> 5    116 0.41287          52  0.7981  0.7883  0.9023
#> 6    241 0.74524          67  3.2009  0.7883  0.9023

# 2. Directly estimate Spatio-Temporal Fay-Herriot model (eblup_stfh)
if (FALSE) { # \dontrun{
fit_stfh <- eblup_stfh(
  y_gaussian ~ x1 + x2,
  data = sim_panel$data,
  domain = ~area,
  time = ~year,
  vardir = ~vardir,
  W = sim_panel$W_std
)
summary(fit_stfh)
} # }
```
