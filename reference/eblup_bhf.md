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
  mse_type = c("bootstrap", "none"),
  n_threads = 1,
  seed = -1,
  print_result = TRUE,
  compute_mse = lifecycle::deprecated()
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

  Number of bootstrap replicates for MSE (if mse_type = "bootstrap").

- mse_type:

  MSE estimator: `"bootstrap"` (default; parametric bootstrap of
  Gonzalez-Manteiga et al. (2008), as in
  [`sae::pbmseBHF()`](https://rdrr.io/pkg/sae/man/pbmseBHF.html)) or
  `"none"` (no MSE; `mse` and `rse` are `NA`).

- n_threads:

  Number of OpenMP threads for the bootstrap MSE. Results depend on
  `seed` only, not on the number of threads.

- seed:

  Random seed for reproducibility; `-1` (default) sets no seed. A fixed
  seed does not change the random number stream of the session.

- print_result:

  Print results (default TRUE).

- compute_mse:

  Deprecated (fastsae 1.0.0); use `mse_type`. `compute_mse = TRUE` is
  `mse_type = "bootstrap"` and `FALSE` is `"none"`.

## Value

List containing EBLUP estimates, fit, and optionally MSE. The `fit`
element holds the variance components, `beta`, the predicted random
effects and the unit-level conditional residuals.

## Details

The nested error model \\y\_{dj} = x\_{dj}'\beta + u_d + e\_{dj}\\ is
fitted by REML or ML in compiled code. With \\\gamma =
\sigma_u^2/\sigma_e^2\\, \\\beta\\ and \\\sigma_e^2\\ are profiled out,
so the fit is a one-dimensional maximisation over \\\gamma \ge 0\\
(including the boundary \\\sigma_u^2 = 0\\); each evaluation uses only
within-domain cross-products and domain means. The estimates are those
of [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) and
[`sae::eblupBHF`](https://rdrr.io/pkg/sae/man/eblupBHF.html) up to the
optimiser tolerance of those packages. The standard errors in `estcoef`
are from \\(X'\hat V^{-1}X)^{-1}\\.

The bootstrap MSE refits the model in each replicate with the same code
and runs replicates in parallel (`n_threads`); random draws are
generated sequentially from R's RNG, so a given `seed` gives the same
result for any number of threads.

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
#> ℹ Computing bootstrap MSE with B = 100 replicates...
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_bhf(formula = CornHec ~ CornPix + SoyBeansPix, unit_data =
#> df_cornsoybean, Xpop = df_meanpop, domain_var = "CountyIndex", popsize_var =
#> "PopnSegments")
#> 
#> ✔ Convergence: Yes (in - iterations)
#> Model: Battese-Harter-Fuller (Unit-level)
#> Random effect variance (sigma2_u): 63.31489 
#> 
#> Fixed Effects Coefficients:
#>                  beta std.error    zvalue pvalue
#> (Intercept) 17.963979 30.974505  0.579960 0.5619
#> CornPix      0.366335  0.064959  5.639511 0.0000
#> SoyBeansPix -0.030364  0.067576 -0.449327 0.6532
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain    eblup samp_size      mse      rse
#> 1      1 122.5825         1 64.29390 6.541184
#> 2      2 123.5274         1 61.19039 6.332544
#> 3      3 113.0343         1 87.53427 8.277114
#> 4      4 114.9901         2 64.10314 6.962726
#> 5      5 137.2660         3 77.81487 6.426408
#> 6      6 108.9807         3 59.14858 7.057040
#> ... and 6 more rows.
#> 
```
