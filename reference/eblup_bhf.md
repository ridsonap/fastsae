# Empirical Best Linear Unbiased Prediction (EBLUP) for the Battese-Harter-Fuller Model

This function estimates small area means using the unit-level model
proposed by Battese, Harter, and Fuller (1988), which combines survey
data (sample units) and auxiliary population information.

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

  An object of class `formula` describing the model.

- unit_data:

  A `data.frame` containing the unit-level survey data.

- Xpop:

  A `data.frame` with one row per domain: the domain identifier
  (`domain_var`), the population size (`popsize_var`) and the population
  means of the covariates, in columns named as the covariates in
  `formula`. Domains in `Xpop` without sampled units are predicted
  synthetically.

- domain_var:

  A character string giving the column name for domain identifier.

- popsize_var:

  A character string for population size variable.

- method:

  Fitting method: "REML" (default) or "ML".

- popnmean_xpop:

  Optional matrix or data frame of the population means of the
  covariates (one row per domain, in the order of `Xpop`, columns as in
  the design matrix, with or without the intercept); if `NULL`, they are
  taken from `Xpop`.

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

An object of class `c("fastsae", "fastsae_unit")`, a list with
`df_eblup` (one row per domain: `domain`, `eblup`, `samp_size`, `mse`,
`rse`), `estcoef` (coefficients, standard errors, z-statistics and
p-values), `random_effect_var`, `fit` (the fitting method, the variance
components `random_effect_var` and `sigma2_e`, `beta`, the predicted
random effects and the unit-level conditional residuals), `convergence`,
`mse_type`, `formula`, `model`, `data` and `call`. `goodness` is a named
vector: `loglikelihood` (the full Gaussian log-likelihood of the sampled
units at the estimates, for ML and REML fits alike, as in
[`eblup_fh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_fh.md)),
and `AIC` and `BIC` computed from it with \\p + 2\\ parameters and the
number of sampled units.

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
and satellite data. *Journal of the American Statistical Association*,
83(401), 28-36.

Gonzalez-Manteiga, W., Lombardia, M. J., Molina, I., Morales, D., and
Santamaria, L. (2008). Bootstrap mean squared error of a small-area
EBLUP. *Journal of Statistical Computation and Simulation*, 78(5),
443-462.

## See also

[`hb_unit`](https://ridsonap.is-a.dev/fastsae/reference/hb_unit.md) for
the hierarchical Bayes counterpart.

## Examples

``` r

# the population means of the covariates must carry the names used in the formula
df_meanpop <- cornsoybeanmeans
names(df_meanpop)[match(c("MeanCornPixPerSeg", "MeanSoyBeansPixPerSeg"), names(df_meanpop))] <-
  c("CornPix", "SoyBeansPix")
df_cornsoybean <- cornsoybean
names(df_cornsoybean)[names(df_cornsoybean) == "County"] <- "CountyIndex"

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
