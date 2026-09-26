# Empirical Best Linear Unbiased Prediction based on a Spatial Fay-Herriot Model.

This function gives the Spatial Empirical Best Linear Unbiased
Prediction (EBLUP) or Empirical Best (EB) predictor under normality
based on a Fay-Herriot model.

## Usage

``` r
eblup_sfh(
  formula,
  vardir,
  domain = NULL,
  data,
  method = c("REML", "ML"),
  mse_method = c("analytical", "pbmse", "npbmse"),
  W = NULL,
  B = 100,
  n_threads = 1,
  seed = -1,
  maxiter = 100,
  precision = 1e-04,
  print_result = TRUE
)
```

## Arguments

- formula:

  an object of class formula that contains a description of the model to
  be fitted. The variables included in the formula must be contained in
  the data.

- vardir:

  vector or column names from data that contain variance sampling from
  the direct estimator for each area.

- domain:

  vector, column name or one-sided formula referencing a domain names
  column in `data`. If NULL, the domains are numbered consecutively.

- data:

  a data frame or a data frame extension (e.g. a tibble).

- method:

  Fitting method can be chosen between 'ML' and 'REML'.

- mse_method:

  a character string determining the estimation method of the MSE.
  Methods that can be chosen: "analytical", "pbmse", or "npbmse". When
  there are unsampled domains, "pbmse"/"npbmse" bootstrap MSE is only
  defined for the sampled domains; unsampled domains automatically get a
  full-spatial synthetic (kriging) prediction and analytical MSE merged
  into the result regardless of `mse_method`.

- W:

  A square matrix with dimension equal to the TOTAL number of domains in
  `data` (including any unsampled domains where the response is `NA`).
  It should contain the row-standardized spatial weights (proximities)
  between ALL domains, with values ranging from 0 to 1. Rows and columns
  must be ordered consistently with the domain identifiers in `data`. Do
  NOT pre-subset `W` to sampled domains only – unsampled domains'
  spatial MSE/ prediction relies on their relation (in `W`) to sampled
  neighbors.

- B:

  Number of bootstrap replications when mse_method = "pbmse" or
  "npbmse".

- n_threads:

  Number of threads used in parallel computation (default 1). Values
  less than or equal to 0 use all available cores.

- seed:

  Integer seed for bootstrap resampling. A value of -1 leaves the
  current R RNG state unchanged.

- maxiter:

  maximum number of iterations allowed in the Fisher-scoring algorithm.
  Default is 100 iterations.

- precision:

  convergence tolerance limit for the Fisher-scoring algorithm. Default
  value is 0.0001.

- print_result:

  print coefficient or not, default value is TRUE.

## Value

The function returns a list with the following objects: `estcoef` a data
frame with the estimated model coefficients in the first column (beta),
their asymptotic standard errors in the second column (std.error), the
t-statistics in the third column (tvalue) and the p-values of the
significance of each coefficient in last column (pvalue)\
`random_effect_var` estimated random effect variance\
`rho` estimated spatial autocorrelation parameter\
`estvarcomp` a data frame with parameter, estimate, and std.error\
`goodness` vector containing several goodness-of-fit measures:
loglikelihood, AIC, and BIC\
`df_eblup` a data frame that contains the following columns:\

- `y` variable response\

- `eblup` estimated results for each area\

- `random_effect` random effect for each area\

- `vardir` variance sampling from the direct estimator for each area\

- `mse` Mean Square Error\

- `rse` Relative Standart Error (%)\

- `mse_pb` or `mse_npb` Parametric / Non Parametric bootstrap MSE\

- `mse_pbbc` or `mse_npbbc` Bias Corrected Parametric / Non Parametric
  bootstrap MSE\

## Details

The model has a form that is response ~ auxiliary variables. where
numeric type response variables can contain NA. When the response
variable contains NA, that domain is treated as unsampled: it is
predicted with a full-spatial synthetic (kriging) estimator that borrows
strength from sampled neighbors via `W`, with an analytical MSE
approximation.

## References

1.  Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley &
    Sons.

## Examples

``` r
library(fastsae)

# Spatial Fay-Herriot model
m1 <- eblup_sfh(
  y ~ x1 + x2 + x3,
  data = mys,
  vardir = ~vardir,
  W = mys_proxmat
)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_sfh(formula = y ~ x1 + x2 + x3, vardir = ~vardir, data = mys, W =
#> mys_proxmat)
#> 
#> ✔ Convergence: Yes (in 8 iterations)
#> Model: Spatial Fay-Herriot (Area-level SAR)
#> Method: REML
#> Random effect variance (sigma2_u): 1.493692 
#> Spatial autocorrelation (rho): -1 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue pvalue
#> (Intercept)  3.03659649  0.62669962  4.84537788 0.0000
#> x1          -0.00085935  0.00824230 -0.10426151 0.9170
#> x2           0.04983809  0.05273618  0.94504545 0.3446
#> x3           0.03929798  0.05006522  0.78493574 0.4325
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    vardir    eblup random_effect       mse       rse
#> 1      1 8.359527 0.6618838 7.282308   2.605985416 0.4797992  9.511756
#> 2      2 7.599650 0.8374691 6.423856   2.212974942 0.6156692 12.214562
#> 3      3 5.514137 0.8822257 5.023421   0.892393688 0.5519476 14.789360
#> 4      4 3.869326 0.6581716 4.324716  -1.127251956 0.4811841 16.039766
#> 5      5 6.305063 1.2788021 6.344991  -0.002629647 0.7315241 13.479795
#> 6      6 3.926807 0.3878004 4.087497  -0.704885769 0.3403445 14.272560
#> ... and 36 more rows.
#> 

# Spatial Fay-Herriot model with Parametric Bootstrap MSE
m2 <- eblup_sfh(
  y ~ x1 + x2 + x3,
  data = mys,
  vardir = ~vardir,
  mse_method = "pbmse",
  B = 50,
  W = mys_proxmat
)
#> ℹ 10 unsampled domain(s) detected. Bootstrap MSE (pbmse) is computed for the 32 sampled domain(s); unsampled domain(s) get a full-spatial synthetic (kriging) prediction and analytical MSE instead.
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_sfh(formula = y ~ x1 + x2 + x3, vardir = ~vardir, data = mys, mse_method
#> = "pbmse", W = mys_proxmat, B = 50)
#> 
#> ✔ Convergence: Yes (in 8 iterations)
#> Model: Spatial Fay-Herriot (Area-level SAR)
#> Method: REML
#> Random effect variance (sigma2_u): 1.493692 
#> Spatial autocorrelation (rho): -1 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue pvalue
#> (Intercept)  3.03659649  0.62669962  4.84537788 0.0000
#> x1          -0.00085935  0.00824230 -0.10426151 0.9170
#> x2           0.04983809  0.05273618  0.94504545 0.3446
#> x3           0.03929798  0.05006522  0.78493574 0.4325
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    vardir    eblup random_effect       mse       rse
#> 1      1 8.359527 0.6618838 7.282308   2.605985416 0.4797992  9.511756
#> 2      2 7.599650 0.8374691 6.423856   2.212974942 0.6156692 12.214562
#> 3      3 5.514137 0.8822257 5.023421   0.892393688 0.5519476 14.789360
#> 4      4 3.869326 0.6581716 4.324716  -1.127251956 0.4811841 16.039766
#> 5      5 6.305063 1.2788021 6.344991  -0.002629647 0.7315241 13.479795
#> 6      6 3.926807 0.3878004 4.087497  -0.704885769 0.3403445 14.272560
#>      mse_pb  mse_pbbc
#> 1 0.5249267 0.5478423
#> 2 0.6548961 0.6589172
#> 3 0.6786196 0.6981591
#> 4 0.4839694 0.5513615
#> 5 0.6318688 0.8754840
#> 6 0.2709669 0.3564423
#> ... and 36 more rows.
#> 
```
