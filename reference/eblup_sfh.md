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
  self_benchmark = FALSE,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL,
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
  grouping / stratum column in `data` (e.g. province) for group-specific
  self-benchmarking.

- print_result:

  print coefficient or not, default value is TRUE.

## Value

The function returns a list with the following objects: `estcoef` a data
frame with the estimated model coefficients (`beta`), their asymptotic
standard errors (`std.error`, repeated as `stderr_beta` for
compatibility with fastsae 0.1.0), the z-statistics (`zvalue`) and the
normal p-values (`pvalue`)  
`random_effect_var` estimated random effect variance  
`rho` estimated spatial autocorrelation parameter  
`estvarcomp` a data frame with parameter, estimate, and std.error  
`goodness` named vector: `loglikelihood` (the full Gaussian
log-likelihood at the estimates), `AIC` and `BIC` computed from it and,
for `method = "REML"`, `loglikelihood_reml` (the restricted
log-likelihood that REML maximises). AIC and BIC use the full
log-likelihood for ML and REML fits alike, as sae does and as in
Marhuenda, Morales and Pardo (2014): the restricted likelihood depends
on the covariates, so only the full one compares models with different
covariates.  
`df_eblup` a data frame that contains the following columns:  

- `y` variable response  

- `eblup` estimated results for each area  

- `random_effect` random effect for each area  

- `vardir` variance sampling from the direct estimator for each area  

- `mse` Mean Square Error  

- `rse` Relative Standart Error (%)  

- `mse_pb` or `mse_npb` Parametric / Non Parametric bootstrap MSE  

- `mse_pbbc` or `mse_npbbc` Bias Corrected Parametric / Non Parametric
  bootstrap MSE  

## Details

The model has a form that is response ~ auxiliary variables. where
numeric type response variables can contain NA. When the response
variable contains NA, that domain is treated as unsampled: it is
predicted with a full-spatial synthetic (kriging) estimator that borrows
strength from sampled neighbors via `W`. The SAR process is defined on
all domains of `W`, so the sampled domains have covariance \\\sigma^2
\Sigma\[s,s\]\\ with \\\Sigma = \[(I - \rho W')(I - \rho W)\]^{-1}\\.

Analytical MSE. Sampled domains: \\g_1 + g_2 + 2 g_3 - g_4\\ (Singh et
al., 2005; as [`sae::mseSFH`](https://rdrr.io/pkg/sae/man/mseSFH.html)),
minus \\b\_{ML}'\nabla g_1\\ for ML. Unsampled domains: \\g_1 + g_2 +
g_3 - \frac{1}{2}\mathrm{tr}(I^{-1}\nabla^2 g_1)\\ (minus
\\b\_{ML}'\nabla g_1\\ for ML), with \\g_3\\ from the exact derivative
of the kriging predictor; for sampled domains this second-order form
agrees, up to o(1/m), with \\g_1 + g_2 + 2 g_3 - g_4\\. The inverse
information used for \\g_3\\ and \\b\_{ML}\\ is the one of the fitting
method: REML information for REML and ML information for ML.
([`sae::mseSFH`](https://rdrr.io/pkg/sae/man/mseSFH.html) uses the REML
information also for ML; the two agree asymptotically, but with few
domains, or \\\hat\sigma^2 \approx 0\\ with \\\hat\rho \approx 1\\, the
REML choice can make the ML bias term explode. Hence, for
`method = "ML"`, the MSE differs from
[`sae::mseSFH`](https://rdrr.io/pkg/sae/man/mseSFH.html); REML results
agree.)

When \\\hat\rho\\ lies on the boundary (\\\|\hat\rho\| \ge 0.999\\, or
reported as \\-1\\), the Prasad-Rao expansion in \\\rho\\ is not valid
and can explode, so \\\rho\\ is treated as fixed: only the \\\sigma^2\\
block of the information is inverted, which removes the \\\rho\\ terms
from \\g_3\\, \\g_4\\ and \\b\_{ML}\\ (the standard error of \\\rho\\ is
then `NA`). In this case the MSE differs from
[`sae::mseSFH`](https://rdrr.io/pkg/sae/man/mseSFH.html), which keeps
those terms.

If the classical (sae-compatible) Fisher scoring fails, or converges to
a negative \\\sigma^2\\ (which
[`sae::eblupSFH`](https://rdrr.io/pkg/sae/man/eblupSFH.html) truncates
to 0 without re-estimating \\\rho\\), the (RE)ML likelihood is maximised
under the constraints \\\sigma^2 \ge 0\\, \\\|\rho\| \le 0.999\\. A
lower-bound estimate \\\rho = -0.999\\ is reported as \\-1\\ (as `sae`)
only when \\I + W\\ is non-singular (not for bipartite graphs such as a
rook grid).

Computation. No \\D \times D\\ inverse is formed: with \\A = (I - \rho
W')(I - \rho W)\\ and \\M = A + \sigma^2 E'\Psi^{-1}E\\ (\\E\\ selects
the sampled domains), \\V^{-1} = \Psi^{-1} -
\sigma^2\Psi^{-1}EM^{-1}E'\Psi^{-1}\\, \\g_1 =
\sigma^2\\\mathrm{diag}(M^{-1})\\ and the remaining terms follow from
Cholesky solves with \\A\\ and \\M\\. These are sparse for a sparse `W`
(contiguity or k-nearest-neighbour weights), which makes models with a
few thousand domains feasible; a dense `W` (e.g. inverse-distance
weights) uses dense factorisations.

## References

1.  Marhuenda, Y., Morales, D., & Pardo, M. C. (2014). Information
    criteria for Fay-Herriot model selection. *Computational Statistics
    & Data Analysis*, 70, 268-280.

2.  Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley &
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
#>      mse_pb  mse_pbbc
#> 1 0.3731986 0.5732823
#> 2 0.6002193 0.6787076
#> 3 0.4744233 0.7493375
#> 4 0.4001204 0.5693914
#> 5 0.7692070 0.8829821
#> 6 0.3054537 0.3742460
#> ... and 36 more rows.
#> 
```
