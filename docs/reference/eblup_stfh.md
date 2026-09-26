# Empirical Best Linear Unbiased Prediction based on a Spatio-Temporal Fay-Herriot Model.

This function gives the Spatio-Temporal Empirical Best Linear Unbiased
Prediction (EBLUP) under normality based on a spatio-temporal
Fay-Herriot model. It reimplements the same Fisher-scoring algorithm as
`eblupSTFH()` (package sae, Marhuenda, Molina & Morales 2013), but the
estimation loop runs in compiled C++/Armadillo, making it substantially
faster and more memory-efficient than the original R implementation –
especially for a large number of domains/time periods.

## Usage

``` r
eblup_stfh(
  formula,
  vardir,
  data,
  domain,
  time,
  W,
  model = c("ST", "S"),
  maxiter = 100,
  precision = 1e-04,
  sigma21_start = NULL,
  rho1_start = 0.5,
  sigma22_start = NULL,
  rho2_start = 0.5,
  compute_mse = FALSE,
  B = 100,
  n_threads = 1,
  seed = -1,
  print_result = TRUE
)
```

## Arguments

- formula:

  an object of class formula describing the model to fit (response ~
  auxiliary variables). Variables must be present in `data`.

- vardir:

  vector, column name or one-sided formula referencing a column in
  `data`, with the sampling variances of the direct estimator.

- data:

  a data frame (or extension) with `domain * time` rows, sorted so that
  all `time` periods of domain 1 come first, then all periods of domain
  2, and so on (i.e. domain-major order) – exactly as required by
  `eblupSTFH()`.

- domain:

  vector, column name or one-sided formula referencing a domain names
  column in `data`. If NULL, the domains are numbered consecutively.

- time:

  vector, column name, or one-sided formula referencing a time names
  column in `data`.

- W:

  a square proximity/spatial weights matrix of dimension
  `domain x domain` (row-standardized, values typically in \\\[0,1\]\\).

- model:

  character, either `"ST"` (spatio-temporal, default) or `"S"` (spatial
  only, no AR(1) temporal component).

- maxiter:

  maximum number of Fisher-scoring iterations. Default 100.

- precision:

  convergence tolerance for the Fisher-scoring algorithm. Default
  `1e-4`.

- sigma21_start, rho1_start, sigma22_start, rho2_start:

  starting values for the variance/autocorrelation components. Defaults
  mirror `eblupSTFH()`: `0.5 * median(vardir)` for the variances and
  `0.5` for the autocorrelations. `rho2_start` is ignored when
  `model = "S"`.

- compute_mse:

  logical, if `TRUE` computes parametric bootstrap MSE using `B`
  bootstrap replicates. Default `FALSE`.

- B:

  number of bootstrap replicates for MSE computation. Only used when
  `compute_mse = TRUE`. Default `100`.

- n_threads:

  number of threads for parallel bootstrap MSE computation. Use `0` for
  all available cores. Default `1`.

- seed:

  random seed for bootstrap. Use `-1` for no seed. Default `-1`.

- print_result:

  print the estimated coefficients or not. Default `TRUE`.

## Value

A list with the same structure as `seblup_area()`, with additional
`estvarcomp` for spatio-temporal variance/autocorrelation components:

- `estcoef`:

  data frame with beta, std.error, tvalue, pvalue.

- `estvarcomp`:

  data frame with estimate, std.error for sigma21, rho1, sigma22, rho2.

- `goodness`:

  vector with loglike, AIC, BIC.

- `df_eblup`:

  data frame with eblup, mse, rse, random_effect_u1, random_effect_u2.

- `model`:

  model type ("ST" or "S").

- `convergence`:

  logical, whether the algorithm converged.

- `n_iter`:

  number of iterations.

- `B`:

  number of bootstrap replicates (only when compute_mse = TRUE).

## Details

This function requires a complete panel (no NA values in the response
variable). If your data contains NA values (unsampled areas/time
periods), please filter them out before calling this function. Future
versions may support automatic handling of unsampled areas.

## References

1.  Marhuenda, Y., Molina, I., & Morales, D. (2013). Small area
    estimation with spatio-temporal Fay-Herriot models. Computational
    Statistics & Data Analysis, 58, 308-325.

2.  Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley &
    Sons.

## Examples

``` r
library(fastsae)
library(dplyr)

mys_panel_nona <- mys_panel |>
  filter(!is.na(y) & year >= 2024)

# Basic EBLUP without MSE
m1 <- eblup_stfh(
  y ~ x1 + x2 + x3,
  data = mys_panel_nona,
  vardir = ~vardir,
  domain = ~area,
  time = ~year,
  W = mys_proxmat[-c(21, 25), -c(21, 25)],
  model = "ST"
)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_stfh(formula = y ~ x1 + x2 + x3, vardir = ~vardir, data = mys_panel_nona,
#> domain = ~area, time = ~year, W = mys_proxmat[-c(21, 25), -c(21, 25)], model =
#> "ST")
#> 
#> ✔ Convergence: Yes (in 12 iterations)
#> Model: Spatio-Temporal Fay-Herriot (ST-FH)
#> Method: REML
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.7498891  0.4538487  8.2624211 0.0000
#> x1          -0.0053032  0.0046392 -1.1431391 0.2530
#> x2           0.0606926  0.0258702  2.3460387 0.0190
#> x3           0.0173309  0.0257334  0.6734780 0.5006
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain time        y    eblup    vardir mse rse mse_pb random_effect_u1
#> 1      1 2024 8.408797 7.967392 0.6684322  NA  NA     NA                0
#> 2      2 2024 7.746937 7.295214 0.8457547  NA  NA     NA                0
#> 3      3 2024 5.902558 6.026662 0.8909541  NA  NA     NA                0
#> 4      4 2024 4.326688 4.591993 0.6646833  NA  NA     NA                0
#> 5      5 2024 6.186035 5.939631 1.2914540  NA  NA     NA                0
#> 6      6 2024 4.455711 4.535236 0.3916371  NA  NA     NA                0
#>   random_effect_u2
#> 1        3.0829442
#> 2        2.7888151
#> 3        1.2876152
#> 4       -1.1934477
#> 5       -0.5181050
#> 6       -0.7532067
#> ... and 114 more rows.
#> 

# EBLUP with parametric bootstrap MSE
m2 <- eblup_stfh(
  y ~ x1 + x2 + x3,
  data = mys_panel_nona,
  vardir = ~vardir,
  domain = ~area,
  time = ~year,
  W = mys_proxmat[-c(21, 25), -c(21, 25)],
  model = "ST",
  compute_mse = TRUE,
  B = 100,
  seed = 42
)
#> ℹ Computing parametric bootstrap MSE with B = 100 replicates...
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_stfh(formula = y ~ x1 + x2 + x3, vardir = ~vardir, data = mys_panel_nona,
#> domain = ~area, time = ~year, W = mys_proxmat[-c(21, 25), -c(21, 25)], model =
#> "ST", compute_mse = TRUE, B = 100, seed = 42)
#> 
#> ✔ Convergence: Yes (in 12 iterations)
#> Model: Spatio-Temporal Fay-Herriot (ST-FH)
#> Method: REML
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.7498891  0.4538487  8.2624211 0.0000
#> x1          -0.0053032  0.0046392 -1.1431391 0.2530
#> x2           0.0606926  0.0258702  2.3460387 0.0190
#> x3           0.0173309  0.0257334  0.6734780 0.5006
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain time        y    eblup    vardir       mse       rse    mse_pb
#> 1      1 2024 8.408797 7.967392 0.6684322 0.5375116  9.201904 0.5375116
#> 2      2 2024 7.746937 7.295214 0.8457547 0.4766552  9.463769 0.4766552
#> 3      3 2024 5.902558 6.026662 0.8909541 0.7094756 13.976290 0.7094756
#> 4      4 2024 4.326688 4.591993 0.6646833 0.4692162 14.917131 0.4692162
#> 5      5 2024 6.186035 5.939631 1.2914540 0.9631173 16.522667 0.9631173
#> 6      6 2024 4.455711 4.535236 0.3916371 0.2681783 11.418575 0.2681783
#>   random_effect_u1 random_effect_u2
#> 1                0        3.0829442
#> 2                0        2.7888151
#> 3                0        1.2876152
#> 4                0       -1.1934477
#> 5                0       -0.5181050
#> 6                0       -0.7532067
#> ... and 114 more rows.
#> 
```
