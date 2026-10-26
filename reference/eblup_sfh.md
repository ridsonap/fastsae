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
  mse_type = c("analytical", "pbmse", "npbmse", "none"),
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
  print_result = TRUE,
  mse_method = lifecycle::deprecated()
)
```

## Arguments

- formula:

  an object of class formula that contains a description of the model to
  be fitted. The variables included in the formula must be contained in
  the data.

- vardir:

  sampling variances of the direct estimates: a numeric vector, the name
  of a column of `data` or a one-sided formula (e.g. `~vardir`).

- domain:

  vector, column name or one-sided formula referencing a domain names
  column in `data`. If NULL, the domains are numbered consecutively.

- data:

  a data frame or a data frame extension (e.g. a tibble).

- method:

  Fitting method can be chosen between 'ML' and 'REML'.

- mse_type:

  MSE estimator: "analytical" (default), "pbmse" (parametric bootstrap),
  "npbmse" (nonparametric bootstrap) or "none" (no MSE; `mse` and `rse`
  are `NA`). When there are unsampled domains, "pbmse"/"npbmse"
  bootstrap MSE is only defined for the sampled domains; unsampled
  domains automatically get a full-spatial synthetic (kriging)
  prediction and analytical MSE merged into the result.

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

  Number of bootstrap replications when mse_type = "pbmse" or "npbmse".

- n_threads:

  Number of threads used in parallel computation (default 1). Values
  less than or equal to 0 use all available cores.

- seed:

  Integer seed for bootstrap resampling; `-1` (default) sets no seed. A
  fixed seed does not change the random number stream of the session.

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

- mse_method:

  Deprecated (fastsae 1.0.0); use `mse_type`.

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
log-likelihood at the estimates, for ML and REML fits alike), `AIC` and
`BIC` computed from it, as sae does and as in Marhuenda, Morales and
Pardo (2014). The restricted (REML) likelihood is not reported: it
depends on the covariates through X, so it cannot compare models with
different covariates.  
`df_eblup` a data frame with one row per domain:  

- `domain` domain identifier  

- `y` direct estimate (`NA` for unsampled domains)  

- `eblup` estimated results for each area  

- `random_effect` random effect for each area  

- `vardir` sampling variance of the direct estimate  

- `mse` Mean Square Error (the bootstrap MSE for `"pbmse"`/`"npbmse"`)  

- `rse` Relative Standard Error (%)  

- `mse_pb` or `mse_npb` Parametric / Non Parametric bootstrap MSE  

- `mse_pbbc` or `mse_npbbc` Bias Corrected Parametric / Non Parametric
  bootstrap MSE  
  `convergence`, `n_iter` whether the fit converged, and its
  iterations,  
  `model`, `method`, `mse_type`, `formula`, `call`, `data`, `W` the
  settings, the matched call, the data and the weight matrix used.  

## Details

The spatial Fay-Herriot model (Petrucci and Salvati, 2006; Pratesi and
Salvati, 2008) is \\y = X\beta + u + e\\ with a simultaneous
autoregressive area effect \\u = \rho W u + \varepsilon\\, \\\varepsilon
\sim N(0, \sigma^2 I)\\, and sampling errors \\e \sim N(0,
\mathrm{diag}(\psi_d))\\ with known \\\psi_d\\ (`vardir`). When the
response is `NA`, that domain is treated as unsampled: it is predicted
with a full-spatial synthetic (kriging) estimator that borrows strength
from sampled neighbors via `W`. The SAR process is defined on all
domains of `W`, so the sampled domains have covariance \\\sigma^2
\Sigma\[s,s\]\\ with \\\Sigma = \[(I - \rho W')(I - \rho W)\]^{-1}\\.

An estimate \\\hat\sigma^2 = 0\\ makes the EBLUPs synthetic estimates,
and the MSE estimates can then understate the error; the function says
so with a message (see Li and Lahiri, 2010).

Analytical MSE. Sampled domains: \\g_1 + g_2 + 2 g_3 - g_4\\ (Singh et
al., 2005; as [`sae::mseSFH`](https://rdrr.io/pkg/sae/man/mseSFH.html);
the bootstrap MSEs are those of Molina et al., 2009), minus
\\b\_{ML}'\nabla g_1\\ for ML. Unsampled domains: \\g_1 + g_2 + g_3 -
\frac{1}{2}\mathrm{tr}(I^{-1}\nabla^2 g_1)\\ (minus \\b\_{ML}'\nabla
g_1\\ for ML), with \\g_3\\ from the exact derivative of the kriging
predictor; for sampled domains this second-order form agrees, up to
o(1/m), with \\g_1 + g_2 + 2 g_3 - g_4\\. The inverse information used
for \\g_3\\ and \\b\_{ML}\\ is the one of the fitting method: REML
information for REML and ML information for ML.
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

1.  Li, H., and Lahiri, P. (2010). An adjusted maximum likelihood method
    for solving small area estimation problems. *Journal of Multivariate
    Analysis*, 101(4), 882-892.

2.  Marhuenda, Y., Morales, D., & Pardo, M. C. (2014). Information
    criteria for Fay-Herriot model selection. *Computational Statistics
    & Data Analysis*, 70, 268-280.

3.  Molina, I., Salvati, N., & Pratesi, M. (2009). Bootstrap for
    estimating the MSE of the spatial EBLUP. *Computational Statistics*,
    24(3), 441-458.

4.  Petrucci, A., & Salvati, N. (2006). Small area estimation for
    spatial correlation in watershed erosion assessment. *Journal of
    Agricultural, Biological, and Environmental Statistics*, 11(2),
    169-182.

5.  Pratesi, M., & Salvati, N. (2008). Small area estimation: the EBLUP
    estimator based on spatially correlated random area effects.
    *Statistical Methods & Applications*, 17(1), 113-141.

6.  Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley &
    Sons.

7.  Singh, B. B., Shukla, G. K., & Kundu, D. (2005). Spatio-temporal
    models in small area estimation. *Survey Methodology*, 31(2),
    183-195.

## See also

[`eblup_fh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_fh.md),
[`eblup_stfh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_stfh.md)
for repeated surveys,
[`spatial_test`](https://ridsonap.is-a.dev/fastsae/reference/spatial_test.md)
to check whether a spatial model is needed,
[`create_weights`](https://ridsonap.is-a.dev/fastsae/reference/create_weights.md)
to build `W`.

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
#> ✔ Convergence: Yes (in 18 iterations)
#> Model: Spatial Fay-Herriot (Area-level SAR)
#> Method: REML
#> Random effect variance (sigma2_u): 0.989226 
#> Spatial autocorrelation (rho): -1 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.2570714  0.6265689  5.1982651 0.0000
#> x1          -0.0020340  0.0080925 -0.2513457 0.8015
#> x2           0.0568516  0.0525488  1.0818827 0.2793
#> x3           0.0252179  0.0497283  0.5071123 0.6121
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    vardir    eblup random_effect       mse       rse
#> 1      1 8.359527 0.6618838 7.535193     2.8128506 0.4999462  9.383551
#> 2      2 7.599650 0.8374691 6.564385     2.2380505 0.5831410 11.633028
#> 3      3 5.514137 0.8822257 5.109460     0.7620572 0.6095214 15.279864
#> 4      4 3.869326 0.6581716 4.273391    -1.2736332 0.5018535 16.577378
#> 5      5 6.305063 1.2788021 6.791422     0.4746535 0.7644430 12.873951
#> 6      6 3.926807 0.3878004 4.016701    -0.9459807 0.3398919 14.514460
#> ... and 36 more rows.
#> 

# Spatial Fay-Herriot model with Parametric Bootstrap MSE
m2 <- eblup_sfh(
  y ~ x1 + x2 + x3,
  data = mys,
  vardir = ~vardir,
  mse_type = "pbmse",
  B = 50,
  W = mys_proxmat
)
#> ℹ 10 unsampled domain(s) detected. Bootstrap MSE (pbmse) is computed for the 32 sampled domain(s); unsampled domain(s) get a full-spatial synthetic (kriging) prediction and analytical MSE instead.
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_sfh(formula = y ~ x1 + x2 + x3, vardir = ~vardir, data = mys, mse_type =
#> "pbmse", W = mys_proxmat, B = 50)
#> 
#> ✔ Convergence: Yes (in 18 iterations)
#> Model: Spatial Fay-Herriot (Area-level SAR)
#> Method: REML
#> Random effect variance (sigma2_u): 0.989226 
#> Spatial autocorrelation (rho): -1 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.2570714  0.6265689  5.1982651 0.0000
#> x1          -0.0020340  0.0080925 -0.2513457 0.8015
#> x2           0.0568516  0.0525488  1.0818827 0.2793
#> x3           0.0252179  0.0497283  0.5071123 0.6121
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    vardir    eblup random_effect       mse       rse
#> 1      1 8.359527 0.6618838 7.535193     2.8128506 0.4999462  9.383551
#> 2      2 7.599650 0.8374691 6.564385     2.2380505 0.5831410 11.633028
#> 3      3 5.514137 0.8822257 5.109460     0.7620572 0.6095214 15.279864
#> 4      4 3.869326 0.6581716 4.273391    -1.2736332 0.5018535 16.577378
#> 5      5 6.305063 1.2788021 6.791422     0.4746535 0.7644430 12.873951
#> 6      6 3.926807 0.3878004 4.016701    -0.9459807 0.3398919 14.514460
#>      mse_pb  mse_pbbc
#> 1 0.5737887 0.5475450
#> 2 0.5894710 0.6517877
#> 3 0.9277117 0.6470065
#> 4 0.4324647 0.5414496
#> 5 0.8710300 0.8592546
#> 6 0.3221574 0.3562749
#> ... and 36 more rows.
#> 
```
