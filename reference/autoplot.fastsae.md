# Autoplot Method for fastsae Objects

Creates diagnostic and comparison plots for Small Area Estimation (SAE)
models fitted with fastsae. This extends the generic
[`autoplot`](https://ggplot2.tidyverse.org/reference/autoplot.html) from
ggplot2.

## Usage

``` r
# S3 method for class 'fastsae'
autoplot(
  object,
  type = c("comparison", "ribbon", "mse", "estimates", "scatter", "rse", "map"),
  ...
)

# S3 method for class 'list'
autoplot(object, type = c("comparison", "mse", "scatter", "map"), ...)
```

## Arguments

- object:

  An object of class `fastsae`, or a (named) list of `fastsae` objects.

- type:

  Type of plot to create.

  - For a single `fastsae` object: `"comparison"`, `"ribbon"`, `"mse"`,
    `"estimates"`, `"scatter"`, `"rse"`, or `"map"`.

  - For a list of `fastsae` objects: `"comparison"`, `"mse"`,
    `"scatter"`, or `"map"`.

- ...:

  Additional arguments passed to internal plotting helpers (e.g.
  `title`) or to ggplot2 layers.

## Value

A `ggplot` object.

## Examples

``` r
library(fastsae)

# Single model plot
fit_fh <- eblup_fh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir")
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_fh(formula = y ~ x1 + x2 + x3, vardir = "vardir", data = mys)
#> 
#> ✔ Convergence: Yes (in 6 iterations)
#> Model: Fay-Herriot (Area-level)
#> Method: eblup
#> Random effect variance (sigma2_u): 2.608103 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.1077510  0.7697687  4.0372527 0.0001
#> x1          -0.0019323  0.0098886 -0.1954019 0.8451
#> x2           0.0555184  0.0614129  0.9040187 0.3660
#> x3           0.0335344  0.0580013  0.5781663 0.5632
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    eblup    vardir random_effect       mse       rse
#> 1      1 8.359527 7.612738 0.6618838    2.94266461 0.5616180  9.844182
#> 2      2 7.599650 6.782316 0.8374691    2.54539588 0.6784873 12.144869
#> 3      3 5.514137 5.187060 0.8822257    0.96692722 0.7190378 16.347620
#> 4      4 3.869326 4.201545 0.6581716   -1.31646412 0.5619314 17.841554
#> 5      5 6.305063 6.323679 1.2788021   -0.03796677 0.9370273 15.307573
#> 6      6 3.926807 4.048590 0.3878004   -0.81904048 0.3548310 14.713194
#> ... and 36 more rows.
#> 
autoplot(fit_fh, type = "estimates")

autoplot(fit_fh, type = "rse")
#> Warning: Removed 10 rows containing missing values or values outside the scale range
#> (`geom_point()`).


# Compare two models
fit_sfh <- eblup_sfh(y ~ x1 + x2 + x3, data = mys, vardir = "vardir", W = mys_proxmat)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_sfh(formula = y ~ x1 + x2 + x3, vardir = "vardir", data = mys, W =
#> mys_proxmat)
#> 
#> ✔ Convergence: Yes (in 9 iterations)
#> Model: Spatial Fay-Herriot (Area-level SAR)
#> Method: REML
#> Random effect variance (sigma2_u): 1.396229 
#> Spatial autocorrelation (rho): -1 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue pvalue
#> (Intercept)  3.02111730  0.61861801  4.88365555 0.0000
#> x1          -0.00072863  0.00805609 -0.09044507 0.9279
#> x2           0.04916379  0.05175135  0.95000019 0.3421
#> x3           0.04007654  0.04917325  0.81500691 0.4151
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    vardir    eblup random_effect       mse       rse
#> 1      1 8.359527 0.6618838 7.229608   2.557838297 0.5180135  9.955333
#> 2      2 7.599650 0.8374691 6.369903   2.167997274 0.6092215 12.253347
#> 3      3 5.514137 0.8822257 4.996497   0.882416743 0.6504225 16.141065
#> 4      4 3.869326 0.6581716 4.343987  -1.095198025 0.5213006 16.620942
#> 5      5 6.305063 1.2788021 6.344684   0.002549394 0.7994159 14.092123
#> 6      6 3.926807 0.3878004 4.094482  -0.683202000 0.3422232 14.287484
#> ... and 36 more rows.
#> 
autoplot(list("Fay-Herriot" = fit_fh, "Spatial FH" = fit_sfh), type = "comparison")

```
