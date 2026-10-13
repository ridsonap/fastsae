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
  mse_type = c("analytical", "bootstrap"),
  B = 100,
  n_threads = 1,
  seed = -1,
  self_benchmark = FALSE,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL,
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

  a data frame (or extension) with `domain * time` rows (balanced
  panel). Rows may be in any order: they are sorted internally into the
  domain-major order required by the model (domains in order of first
  appearance, periods increasing), and the results are returned in the
  original row order. Note that
  [`sae::eblupSTFH()`](https://rdrr.io/pkg/sae/man/eblupSTFH.html)
  instead assumes the data are already in domain-major order.

- domain:

  vector, column name or one-sided formula referencing a domain names
  column in `data`. If NULL, the domains are numbered consecutively.

- time:

  vector, column name, or one-sided formula referencing a time names
  column in `data`.

- W:

  a square proximity/spatial weights matrix of dimension
  `domain x domain` (row-standardized, values typically in \\\[0,1\]\\),
  with rows and columns in the order in which the domains first appear
  in `data`.

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
  `model = "S"`. The REML fit first uses the same Fisher-scoring
  iteration as
  [`sae::eblupSTFH()`](https://rdrr.io/pkg/sae/man/eblupSTFH.html); if
  that fails (typically a variance iterate turns negative and the
  covariance matrix becomes indefinite) it falls back to constrained
  REML (variances \\\ge 0\\, \\\|\rho\| \le 0.999\\) with step halving.
  When a variance is estimated at 0, its autocorrelation is not
  identified and its standard error is reported as `NA`.

- compute_mse:

  logical, if `TRUE` computes MSE. Default `FALSE`.

- mse_type:

  character, method for MSE estimation: `"analytical"` (fast Prasad-Rao
  type estimator g1 + g2 + 2 g3 - g4, default when `compute_mse = TRUE`;
  g3 uses the exact derivative of the EBLUP, \\\Psi P V_a P y\\, and the
  REML information at the final estimates; correlations estimated on the
  boundary \\\|\rho\| \ge 0.999\\ are treated as fixed, so they do not
  enter \\g_3\\/\\g_4\\) or `"bootstrap"` (parametric bootstrap in which
  the variance components are re-estimated by REML in every replicate,
  as in Marhuenda et al., 2013).

- B:

  number of bootstrap replicates for MSE computation. Only used when
  `compute_mse = TRUE` and `mse_type = "bootstrap"`. Default `100`.

- n_threads:

  number of threads for parallel bootstrap MSE computation. Use `0` for
  all available cores. Default `1`.

- seed:

  random seed for bootstrap. Use `-1` for no seed. Default `-1`.

- self_benchmark:

  Logical. If `TRUE`, fits a self-benchmarking model via the augmented
  model approach (Wang, Fuller, and Qu, 2008; Rao and Molina, 2015, Sec
  10.5.1), ensuring that the weighted aggregate of small area estimates
  automatically matches the direct survey aggregate or specified
  benchmark targets without post-hoc adjustment.

- benchmark_weight:

  Optional vector, column name, or one-sided formula referencing the
  domain benchmarking weights in `data` (e.g. population sizes or
  shares). Defaults to equal weights.

- benchmark_target:

  Optional numeric scalar or named vector specifying fixed aggregate
  benchmark target(s). If `NULL` (default), calibrates to the weighted
  aggregate of the direct estimates.

- benchmark_group:

  Optional vector, column name, or one-sided formula referencing the
  grouping / stratum column in `data` (e.g. time period / year) for
  group-specific self-benchmarking.

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

  number of bootstrap replicates (only when compute_mse = TRUE and
  mse_type = "bootstrap").

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
#> ✔ Convergence: Yes (in 7 iterations)
#> Model: Spatio-Temporal Fay-Herriot (ST-FH)
#> Method: REML
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.7761285  0.5612917  6.7275685 0.0000
#> x1          -0.0021031  0.0069943 -0.3006916 0.7636
#> x2           0.0394660  0.0362503  1.0887079 0.2763
#> x3           0.0344259  0.0349429  0.9852041 0.3245
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain time        y    eblup    vardir mse rse mse_pb random_effect_u1
#> 1      1 2024 8.408797 8.218580 0.6684322  NA  NA     NA        3.1764860
#> 2      2 2024 7.746937 7.535062 0.8457547  NA  NA     NA        2.9153948
#> 3      3 2024 5.902558 5.848560 0.8909541  NA  NA     NA        1.1749562
#> 4      4 2024 4.326688 4.223733 0.6646833  NA  NA     NA       -1.4343861
#> 5      5 2024 6.186035 6.394829 1.2914540  NA  NA     NA       -0.1137771
#> 6      6 2024 4.455711 4.343114 0.3916371  NA  NA     NA       -0.7942036
#>   random_effect_u2
#> 1                0
#> 2                0
#> 3                0
#> 4                0
#> 5                0
#> 6                0
#> ... and 114 more rows.
#> 

# Fast Analytical Prasad-Rao MSE (< 0.05 seconds)
m2 <- eblup_stfh(
  y ~ x1 + x2 + x3,
  data = mys_panel_nona,
  vardir = ~vardir,
  domain = ~area,
  time = ~year,
  W = mys_proxmat[-c(21, 25), -c(21, 25)],
  model = "ST",
  compute_mse = TRUE,
  mse_type = "analytical"
)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_stfh(formula = y ~ x1 + x2 + x3, vardir = ~vardir, data = mys_panel_nona,
#> domain = ~area, time = ~year, W = mys_proxmat[-c(21, 25), -c(21, 25)], model =
#> "ST", compute_mse = TRUE, mse_type = "analytical")
#> 
#> ✔ Convergence: Yes (in 7 iterations)
#> Model: Spatio-Temporal Fay-Herriot (ST-FH)
#> Method: REML
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.7761285  0.5612917  6.7275685 0.0000
#> x1          -0.0021031  0.0069943 -0.3006916 0.7636
#> x2           0.0394660  0.0362503  1.0887079 0.2763
#> x3           0.0344259  0.0349429  0.9852041 0.3245
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain time        y    eblup    vardir       mse       rse mse_pb
#> 1      1 2024 8.408797 8.218580 0.6684322 0.2460853  6.035955     NA
#> 2      2 2024 7.746937 7.535062 0.8457547 0.2890231  7.134761     NA
#> 3      3 2024 5.902558 5.848560 0.8909541 0.2994582  9.356625     NA
#> 4      4 2024 4.326688 4.223733 0.6646833 0.2461954 11.747447     NA
#> 5      5 2024 6.186035 6.394829 1.2914540 0.3990277  9.878080     NA
#> 6      6 2024 4.455711 4.343114 0.3916371 0.1919143 10.086777     NA
#>   random_effect_u1 random_effect_u2
#> 1        3.1764860                0
#> 2        2.9153948                0
#> 3        1.1749562                0
#> 4       -1.4343861                0
#> 5       -0.1137771                0
#> 6       -0.7942036                0
#> ... and 114 more rows.
#> 
```
