# Multi-Distribution Spatio-Temporal Panel Dataset

A synthetic panel dataset containing 42 domains observed over 5 time
periods (2022 to 2026, total 210 observations). Features multiple
response variables representing various probability distributions
(Gaussian, Poisson, Binomial, Beta, Negative Binomial, and Gamma) driven
by domain-level spatial autocorrelation and first-order autoregressive
AR(1) temporal dynamics. It shares the 42-domain spatial structure of
[`mys_proxmat`](https://ridsonap.github.io/fastsae/reference/mys_proxmat.md)
and is formatted in domain-major order for direct use in
[`eblup_stfh`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md).

## Usage

``` r
data(sim_panel)
```

## Format

A data frame with 210 rows and 15 variables in domain-major order:

- area:

  Integer domain identifier (1 to 42).

- year:

  Year of observation (2022 to 2026).

- x1:

  Dynamic explanatory covariate with spatial & temporal variation.

- x2:

  Dynamic uniform explanatory covariate.

- y_gaussian:

  Gaussian direct estimator response with unsampled domains as NA.

- vardir:

  Direct sampling variance for Gaussian Fay-Herriot model.

- y_poisson:

  Count response (Poisson) with unsampled domains as NA.

- exposure:

  Expected population count / exposure offset for Poisson model.

- y_binomial:

  Number of successes (Binomial) with unsampled domains as NA.

- trials:

  Sample size / number of trials for Binomial model.

- y_beta:

  Continuous proportion response in (0, 1) for Beta regression.

- y_nbinomial:

  Overdispersed count response for Negative Binomial model.

- y_gamma:

  Skewed positive continuous response for Gamma regression.

- x_coord:

  Centroid x coordinate.

- y_coord:

  Centroid y coordinate.

## Source

Simulated using
[`sim_series_data`](https://ridsonap.github.io/fastsae/reference/sim_series_data.md)
based on the spatial proximity structure of
[`mys_proxmat`](https://ridsonap.github.io/fastsae/reference/mys_proxmat.md).

## Examples

``` r
data(sim_panel)
head(sim_panel)
#>   area year     x1     x2 y_gaussian  vardir y_poisson exposure y_binomial
#> 1    1 2022 1.0967 2.6650     0.1591 0.04717       268      332         84
#> 2    1 2023 1.3387 2.7798         NA 0.15793        NA      363         NA
#> 3    1 2024 1.3375 2.3328     1.2971 0.06945       164      178         99
#> 4    1 2025 1.8651 2.1846     1.8610 0.04719       671      406         76
#> 5    1 2026 1.7945 2.6695     1.9166 0.10708       496      350         52
#> 6    2 2022 2.3614 2.2745         NA 0.11363        NA      258         NA
#>   trials  y_beta y_nbinomial y_gamma x_coord y_coord
#> 1    246 0.26165           2  1.8890  0.0354  0.2881
#> 2    126      NA          NA      NA  0.0354  0.2881
#> 3    239 0.30555           9  0.6434  0.0354  0.2881
#> 4    115 0.80090          73  3.1109  0.0354  0.2881
#> 5     73 0.47064          59  2.3906  0.0354  0.2881
#> 6    108      NA          NA      NA -0.2697  0.0128
```
