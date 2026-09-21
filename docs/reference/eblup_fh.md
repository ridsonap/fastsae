# EBLUPs based on a Fay-Herriot Model.

This function gives the Empirical Best Linear Unbiased Prediction
(EBLUP) or Empirical Best (EB) predictor under normality based on a
Fay-Herriot model.

## Usage

``` r
eblup_fh(
  formula,
  vardir,
  data,
  method = c("REML", "ML"),
  maxiter = 100,
  precision = 1e-04,
  print_result = TRUE
)
```

## Arguments

- formula:

  an object of class formula that contains a description of the model to
  be fitted.

- vardir:

  vector or column names from data that contain variance sampling from
  the direct estimator.

- data:

  a data frame or a data frame extension (e.g. a tibble).

- method:

  Fitting method can be chosen between 'ML' and 'REML'.

- maxiter:

  maximum number of iterations allowed in the Fisher-scoring algorithm.

- precision:

  convergence tolerance limit for the Fisher-scoring algorithm.

- print_result:

  print coefficient or not, default value is TRUE.

## Value

The function returns a list with the following objects: `estcoef` a data
frame with the estimated model coefficients, `random_effect_var`
estimated random effect variance, `goodness` vector containing several
goodness-of-fit measures, `df_eblup` a data frame that contains y,
eblup, vardir, mse, and rse.

## Details

The model has a form that is response ~ auxiliary variables. where
numeric type response variables can contain NA. When the response
variable contains NA it will be estimated with synthetic estimator.

## References

1.  Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley &
    Sons.

## Examples

``` r
library(fastsae)

# Standard Fay-Herriot model
m1 <- eblup_fh(
  y ~ x1 + x2 + x3,
  data = mys,
  vardir = "vardir"
)
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
#>   domain        y    eblup    vardir       mse       rse
#> 1      1 8.359527 7.612738 0.6618838 0.5616180  9.844182
#> 2      2 7.599650 6.782316 0.8374691 0.6784873 12.144869
#> 3      3 5.514137 5.187060 0.8822257 0.7190378 16.347620
#> 4      4 3.869326 4.201545 0.6581716 0.5619314 17.841554
#> 5      5 6.305063 6.323679 1.2788021 0.9370273 15.307573
#> 6      6 3.926807 4.048590 0.3878004 0.3548310 14.713194
#> ... and 36 more rows.
#> 
```
