# Empirical Best Linear Unbiased Prediction based on a Fay-Herriot Model.

This function gives the Empirical Best Linear Unbiased Prediction
(EBLUP) or Empirical Best (EB) predictor under normality based on a
Fay-Herriot model.

## Usage

``` r
eblup_fh(
  formula,
  vardir,
  domain = NULL,
  data,
  method = c("REML", "ML"),
  transform = c("none", "log"),
  mse_type = c("analytical", "none"),
  maxiter = 100,
  precision = 1e-04,
  print_result = TRUE,
  self_benchmark = FALSE,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL
)
```

## Arguments

- formula:

  an object of class formula that contains a description of the model to
  be fitted.

- vardir:

  vector or column names from data that contain variance sampling from
  the direct estimator.

- domain:

  vector, column name or one-sided formula referencing a domain names
  column in `data`. If NULL, the domains are numbered consecutively.

- data:

  a data frame or a data frame extension (e.g. a tibble).

- method:

  Fitting method can be chosen between 'ML' and 'REML'.

- transform:

  Character string specifying data transformation for the response:
  `"none"` (default) or `"log"` (Log-Fay-Herriot model). With `"log"`
  the estimate is back-transformed as \\\exp(\tilde\mu + \hat g_1/2)\\
  with \\g_1 = \sigma^2_u(1-\gamma)\\ (\\\sigma^2_u/2\\ for unsampled
  domains), i.e. the conditional-mean correction under the log-normal
  model (cf. Slud & Maiti, 2006); it does not include the additional
  O(1/m) terms for the estimation of \\\beta\\ and \\\sigma^2_u\\. The
  MSE on the original scale is the delta-method approximation
  \\\hat\theta^2 \cdot mse\_{log}\\.

- mse_type:

  MSE estimator: `"analytical"` (default; Prasad-Rao, with the
  Datta-Lahiri bias correction for ML, as in
  [`sae::mseFH()`](https://rdrr.io/pkg/sae/man/mseFH.html)) or `"none"`
  (no MSE; `mse` and `rse` are `NA`).

- maxiter:

  maximum number of iterations allowed in the Fisher-scoring algorithm.

- precision:

  convergence tolerance limit for the Fisher-scoring algorithm.

- print_result:

  print coefficient or not, default value is TRUE.

- self_benchmark:

  Logical. If `TRUE`, fits a self-benchmarking model via the augmented
  model approach (Wang, Fuller, and Qu, 2008; Rao and Molina, 2015, Sec
  10.5.1), ensuring that the weighted aggregate of small area estimates
  automatically matches the direct survey aggregate or specified
  benchmark targets without post-hoc adjustment. With
  `transform = "log"` that constraint would hold on the log scale only,
  so instead the back-transformed estimates are ratio-benchmarked on the
  original scale (exact constraint; MSE scaled by the squared group
  ratio; the unbenchmarked estimate is kept in `eblup_unbenchmarked`).

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
  grouping / stratum column in `data` (e.g. province) for group-specific
  self-benchmarking.

## Value

The function returns a list with the following objects: `estcoef` a data
frame with the estimated model coefficients (`beta`), their standard
errors (`std.error`, repeated as `stderr_beta` for compatibility with
fastsae 0.1.0), z-statistics (`zvalue`) and normal p-values (`pvalue`),
`random_effect_var` estimated random effect variance, `goodness` named
vector: `loglikelihood` (the full Gaussian log-likelihood at the
estimates, for ML and REML fits alike), `AIC` and `BIC` computed from
it, as sae does and as in Marhuenda, Morales and Pardo (2014). The
restricted (REML) likelihood is not reported: it depends on the
covariates through X, so it cannot compare models with different
covariates., `df_eblup` a data frame that contains y, eblup,
random_effect, vardir, mse, and rse.  

- `y` variable response  

- `eblup` estimated results for each area  

- `random_effect` random effect for each area  

- `vardir` variance sampling from the direct estimator for each area  

- `mse` Mean Square Error  

- `rse` Relative Standart Error (%)  
  `self_benchmark` logical indicating whether self-benchmarking was
  active,  
  `benchmark_summary` summary table of benchmark calibration (if
  self-benchmarking).  

## Details

The model has a form that is response ~ auxiliary variables. where
numeric type response variables can contain NA. When the response
variable contains NA it will be estimated with synthetic estimator.

## References

1.  Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley &
    Sons.

&nbsp;

1.  Marhuenda, Y., Morales, D., & Pardo, M. C. (2014). Information
    criteria for Fay-Herriot model selection. *Computational Statistics
    & Data Analysis*, 70, 268-280.

2.  Rao, J. N. K., & Molina, I. (2015). Small area estimation. John
    Wiley & Sons.

3.  Slud, E. V., & Maiti, T. (2006). Small-area estimation via
    Fay-Herriot models with logit and log transformations.
    *Communications in Statistics - Theory and Methods*, 35(10),
    1957-1971.

4.  Wang, J., Fuller, W. A., and Qu, Y. (2008). Small area estimation
    under a restriction. *Survey Methodology*, 34(1), 29-36.

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
#> Method: REML
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

# Log-Fay-Herriot with Slud-Maiti bias correction
m_log <- eblup_fh(
  y ~ x1 + x2 + x3,
  data = mys,
  vardir = "vardir",
  transform = "log"
)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_fh(formula = y ~ x1 + x2 + x3, vardir = "vardir", data = mys, transform =
#> "log")
#> 
#> ✔ Convergence: Yes (in 4 iterations)
#> Model: Fay-Herriot on the log scale (Area-level)
#> Method: REML
#> Random effect variance (sigma2_u): 0.049822 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue pvalue
#> (Intercept)  1.49622687  0.13824490 10.82301643 0.0000
#> x1          -0.00050299  0.00148560 -0.33857582 0.7349
#> x2           0.00489446  0.00955320  0.51233758 0.6084
#> x3           0.00666673  0.00902326  0.73883885 0.4600
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    eblup    vardir random_effect       mse       rse
#> 1      1 8.359527 7.774951 0.6618838   0.402270166 0.5187454  9.263592
#> 2      2 7.599650 6.938675 0.8374691   0.331933127 0.5970620 11.136102
#> 3      3 5.514137 5.337390 0.8822257   0.071683488 0.6181146 14.730097
#> 4      4 3.869326 4.715172 0.6581716  -0.210829864 0.6116069 16.585887
#> 5      5 6.305063 6.373069 1.2788021  -0.001478566 0.9059793 14.935181
#> 6      6 3.926807 4.362456 0.3878004  -0.191868095 0.3701073 13.945458
#>   eblup_log     mse_log
#> 1  2.046928 0.008581414
#> 2  1.931495 0.012401277
#> 3  1.665568 0.021697575
#> 4  1.539108 0.027509165
#> 5  1.842308 0.022305964
#> 6  1.464679 0.019447580
#> ... and 36 more rows.
#> 

# Self-benchmarking Fay-Herriot model
m1_bench <- eblup_fh(
  y ~ x1 + x2 + x3,
  data = mys,
  vardir = "vardir",
  self_benchmark = TRUE,
  benchmark_weight = "n"
)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_fh(formula = y ~ x1 + x2 + x3, vardir = "vardir", data = mys,
#> self_benchmark = TRUE, benchmark_weight = "n")
#> 
#> ✔ Convergence: Yes (in 6 iterations)
#> Model: Fay-Herriot (Area-level)
#> Method: REML
#> Self-Benchmarking: Active (Wang-Fuller-Qu augmented model)
#> Random effect variance (sigma2_u): 2.687629 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.1312345  0.7794439  4.0172672 0.0001
#> x1          -0.0015269  0.0100761 -0.1515406 0.8795
#> x2           0.0460540  0.0669946  0.6874290 0.4918
#> x3           0.0544386  0.0803391  0.6776108 0.4980
#> bench_aug   -4.8008951 12.4500804 -0.3856116 0.6998
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    eblup    vardir random_effect       mse       rse
#> 1      1 8.359527 7.634955 0.6618838    2.94217773 0.5640221  9.836522
#> 2      2 7.599650 6.805927 0.8374691    2.54723889 0.6821558 12.135412
#> 3      3 5.514137 5.178251 0.8822257    1.02324748 0.7246515 16.439230
#> 4      4 3.869326 4.208629 0.6581716   -1.38553580 0.5654438 17.867100
#> 5      5 6.305063 6.270749 1.2788021    0.07211539 0.9628948 15.648401
#> 6      6 3.926807 4.042689 0.3878004   -0.80311506 0.3556387 14.751431
#> ... and 36 more rows.
#> 
```
