# Empirical Best Linear Unbiased Prediction (EBLUP) for the Battese-Harter-Fuller Model

This function estimates small area means or totals using the unit-level
model proposed by Battese, Harter, and Fuller (1988), which combines
survey data (sample units) and auxiliary population information.

## Usage

``` r
eblup_bhf(
  formula,
  unit_data,
  Xpop,
  domain_var,
  popsize_var,
  method = c("REML", "ML"),
  popnmean_xpop = NULL,
  B = 100,
  mse = FALSE,
  n_threads = 1,
  seed = -1,
  print_result = TRUE
)

pbmse_unit(
  formula,
  unit_data,
  Xpop,
  domain_var,
  popsize_var,
  method = c("REML", "ML"),
  B = 100,
  n_threads = 1,
  seed = -1
)
```

## Arguments

- formula:

  An object of class \`formula\` describing the model.

- unit_data:

  A \`data.frame\` containing the unit-level survey data.

- Xpop:

  A \`data.frame\` containing auxiliary variables and domain info.

- domain_var:

  A character string giving the column name for domain identifier.

- popsize_var:

  A character string for population size variable.

- method:

  Fitting method: "REML" (default) or "ML".

- popnmean_xpop:

  Population mean of auxiliary variables per domain.

- B:

  Number of bootstrap replicates for MSE (if mse = TRUE).

- mse:

  If TRUE, compute bootstrap MSE.

- n_threads:

  Number of threads for parallel computation.

- seed:

  Random seed for reproducibility.

- print_result:

  Print results (default TRUE).

## Value

List containing EBLUP estimates, fit, and optionally MSE.

## References

Battese, G. E., Harter, R. M., and Fuller, W. A. (1988). An
error-components model for prediction of county crop areas using survey
and satellite data. \*Journal of the American Statistical Association\*,
83(401), 28-36.

## Examples

``` r

library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
df_meanpop <- cornsoybeanmeans |>
  rename(CornPix = MeanCornPixPerSeg, SoyBeansPix = MeanSoyBeansPixPerSeg)
df_cornsoybean <- cornsoybean |>
  rename(CountyIndex = County)

res <- eblup_bhf(
  formula = CornHec ~ CornPix + SoyBeansPix,
  Xpop = df_meanpop,
  unit_data = df_cornsoybean,
  domain_var = "CountyIndex",
  popsize_var = "PopnSegments"
)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_bhf(formula = CornHec ~ CornPix + SoyBeansPix, unit_data =
#> df_cornsoybean, Xpop = df_meanpop, domain_var = "CountyIndex", popsize_var =
#> "PopnSegments")
#> 
#> ✔ Convergence: Yes (in - iterations)
#> Model: Battese-Harter-Fuller (Unit-level)
#> Random effect variance (sigma2_u): 63.3149 
#> 
#> Fixed Effects Coefficients:
#>                  beta std.error    zvalue pvalue
#> (Intercept) 17.963979 30.974505  0.579960 0.5619
#> CornPix      0.366335  0.064959  5.639511 0.0000
#> SoyBeansPix -0.030364  0.067576 -0.449327 0.6532
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain    eblup samp_size mse rse
#> 1      1 122.5825         1  NA  NA
#> 2      2 123.5274         1  NA  NA
#> 3      3 113.0343         1  NA  NA
#> 4      4 114.9901         2  NA  NA
#> 5      5 137.2660         3  NA  NA
#> 6      6 108.9807         3  NA  NA
#> ... and 6 more rows.
#> 
```
