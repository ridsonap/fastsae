# Multi-Distribution Synthetic Area-Level Dataset

A synthetic dataset containing 42 domains with multiple response
variables representing various probability distributions (Gaussian,
Poisson, Binomial, Beta, Negative Binomial, and Gamma). Designed for
testing and benchmarking generalized linear and spatial small area
estimation models (EBP / INLA). It shares the 42-domain spatial
structure of
[`mys_proxmat`](https://ridsonap.github.io/fastsae/reference/mys_proxmat.md).

## Usage

``` r
data(sim_area)
```

## Format

A data frame with 42 rows and 14 variables:

- area:

  Integer domain identifier (1 to 42).

- x1:

  Continuous explanatory covariate.

- x2:

  Uniform explanatory covariate.

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
[`sim_area_data`](https://ridsonap.github.io/fastsae/reference/sim_area_data.md)
based on the spatial proximity structure of
[`mys_proxmat`](https://ridsonap.github.io/fastsae/reference/mys_proxmat.md).

## Examples

``` r
data(sim_area)
head(sim_area)
#>   area      x1     x2 y_gaussian  vardir y_poisson exposure y_binomial trials
#> 1    1  2.5206 0.5713     2.5100 0.13373       350      141         96    113
#> 2    2  0.9203 4.8878    -1.8859 0.11747        26       97         11    196
#> 3    3  2.1392 3.1023     1.3503 0.11125       427      330        137    212
#> 4    4  1.9153 4.6371    -0.9947 0.10081        83      216         35    233
#> 5    5  1.3334 1.0112     1.8221 0.04713       223      166         89    129
#> 6    6 -0.5161 1.0776         NA 0.09209        NA      174         NA     61
#>    y_beta y_nbinomial y_gamma x_coord y_coord
#> 1 0.85038          37  2.5916  0.0354  0.2881
#> 2 0.12609           1  0.5153 -0.2697  0.0128
#> 3 0.54495          37  1.1038  0.2234  0.0135
#> 4 0.17570           6  1.1678 -0.0070  0.0953
#> 5 0.61406          16  1.5058  0.1958  0.2200
#> 6      NA          NA      NA -0.2878 -0.0202
```
