# Empirical Best Linear Unbiased Prediction based on a Spatio-Temporal Fay-Herriot Model.

This function gives the Spatio-Temporal Empirical Best Linear Unbiased
Prediction (EBLUP) under normality based on a spatio-temporal
Fay-Herriot model. It reimplements the same Fisher-scoring algorithm as
`eblupSTFH()` (package sae, Marhuenda, Molina & Morales 2013), but the
estimation loop runs in compiled C++/Armadillo, making it substantially
faster and more memory-efficient than the original R implementation –
especially for a large number of domains/time periods.

## Usage

``` r
eblup_stfh(
  formula,
  vardir,
  data,
  domain,
  time,
  W,
  model = c("ST", "S"),
  maxiter = 100,
  precision = 1e-04,
  sigma21_start = NULL,
  rho1_start = 0.5,
  sigma22_start = NULL,
  rho2_start = 0.5,
  mse_type = c("analytical", "bootstrap", "none"),
  B = 100,
  n_threads = 1,
  seed = -1,
  self_benchmark = FALSE,
  benchmark_weight = NULL,
  benchmark_target = NULL,
  benchmark_group = NULL,
  print_result = TRUE,
  compute_mse = lifecycle::deprecated()
)
```

## Arguments

- formula:

  an object of class formula describing the model to fit (response ~
  auxiliary variables). Variables must be present in `data`.

- vardir:

  vector, column name or one-sided formula referencing a column in
  `data`, with the sampling variances of the direct estimator. Must be
  positive in the sampled cells; it is ignored (may be `NA`) where the
  response is `NA`.

- data:

  a data frame (or extension) with one row per domain and time period
  (`domain * time` rows), including the cells without a direct estimate:
  for those, the response (and `vardir`) is `NA`, while the auxiliary
  variables must be available (see Details). Rows may be in any order:
  they are sorted internally into the domain-major order required by the
  model (domains in order of first appearance, periods increasing), and
  the results are returned in the original row order. Note that
  [`sae::eblupSTFH()`](https://rdrr.io/pkg/sae/man/eblupSTFH.html)
  instead assumes a complete panel already in domain-major order.

- domain:

  vector, column name or one-sided formula referencing a domain names
  column in `data`. If NULL, the domains are numbered consecutively.

- time:

  vector, column name, or one-sided formula referencing a time names
  column in `data`.

- W:

  a square proximity/spatial weights matrix of dimension
  `domain x domain` (row-standardized, values typically in \\\[0,1\]\\),
  with rows and columns in the order in which the domains first appear
  in `data`. It must cover all domains, including the ones that are
  never sampled: do not drop them from `W`, since their predictions
  borrow strength from their neighbours through it.

- model:

  character, either `"ST"` (spatio-temporal, default) or `"S"` (spatial
  only, no AR(1) temporal component).

- maxiter:

  maximum number of Fisher-scoring iterations. Default 100.

- precision:

  convergence tolerance for the Fisher-scoring algorithm. Default
  `1e-4`.

- sigma21_start, rho1_start, sigma22_start, rho2_start:

  starting values for the variance/autocorrelation components. Defaults
  mirror `eblupSTFH()`: `0.5 * median(vardir)` for the variances and
  `0.5` for the autocorrelations. `rho2_start` is ignored when
  `model = "S"`. The REML fit first uses the same Fisher-scoring
  iteration as
  [`sae::eblupSTFH()`](https://rdrr.io/pkg/sae/man/eblupSTFH.html); if
  that fails (typically a variance iterate turns negative and the
  covariance matrix becomes indefinite) it falls back to constrained
  REML (variances \\\ge 0\\, \\\|\rho\| \le 0.999\\) with step halving.
  When a variance is estimated at 0, its autocorrelation is not
  identified and its standard error is reported as `NA`.

- mse_type:

  character, method for MSE estimation: `"analytical"` (default; fast
  Prasad-Rao type estimator g1 + g2 + 2 g3 - g4; g3 uses the exact
  derivative of the EBLUP, \\\Psi P V_a P y\\, and the REML information
  at the final estimates; correlations estimated on the boundary
  \\\|\rho\| \ge 0.999\\ are treated as fixed, so they do not enter
  \\g_3\\/\\g_4\\) or `"bootstrap"` (parametric bootstrap in which the
  variance components are re-estimated by REML in every replicate, as in
  Marhuenda et al., 2013) or `"none"` (no MSE; `mse` and `rse` are
  `NA`).

- B:

  number of bootstrap replicates for MSE computation. Only used when
  `mse_type = "bootstrap"`. Default `100`.

- n_threads:

  number of threads for parallel bootstrap MSE computation. Use `0` for
  all available cores. Default `1`.

- seed:

  random seed for bootstrap. Use `-1` for no seed. Default `-1`. A fixed
  seed does not change the random number stream of the session.

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
  grouping / stratum column in `data` (e.g. time period / year) for
  group-specific self-benchmarking.

- print_result:

  print the estimated coefficients or not. Default `TRUE`.

- compute_mse:

  Deprecated (fastsae 1.0.0); use `mse_type`. `compute_mse = FALSE` is
  `mse_type = "none"`; `compute_mse = TRUE` is `mse_type = "bootstrap"`
  when `B` is given and `"analytical"` otherwise.

## Value

A list with the same structure as the result of
[`eblup_sfh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md),
with `estvarcomp` for the spatio-temporal variance and autocorrelation
components:

- `estcoef`:

  data frame with `beta`, `std.error`, `tvalue`, `pvalue`, as in
  [`sae::eblupSTFH()`](https://rdrr.io/pkg/sae/man/eblupSTFH.html), plus
  `stderr_beta` and `zvalue` (copies, for the column names of the other
  fastsae models). The p-values use the normal distribution, so `tvalue`
  is a z-statistic.

- `estvarcomp`:

  data frame with estimate, std.error for sigma21, rho1, sigma22, rho2.

- `goodness`:

  named vector: `loglikelihood`, the full Gaussian log-likelihood at the
  REML estimates (as `loglike` in
  [`sae::eblupSTFH()`](https://rdrr.io/pkg/sae/man/eblupSTFH.html)), and
  `AIC` and `BIC` computed from it. The restricted (REML) likelihood is
  not reported: it depends on the covariates through X, so it cannot
  compare models with different covariates.

- `df_eblup`:

  data frame with domain, time, y, eblup, vardir, mse, rse,
  random_effect_u1, random_effect_u2 and `sampled` (`FALSE` for the
  cells predicted without a direct estimate).

- `model`:

  model type ("ST" or "S").

- `convergence`:

  logical, whether the algorithm converged.

- `n_iter`:

  number of iterations.

- `B`:

  number of bootstrap replicates (only when mse_type = "bootstrap").

## Details

**Estimation method.** The variance components are estimated by REML
only, as in Marhuenda, Molina and Morales (2013) and
[`sae::eblupSTFH()`](https://rdrr.io/pkg/sae/man/eblupSTFH.html), which
derive the Fisher-scoring algorithm, the analytical MSE and the
bootstrap for REML; there is no `method` argument. The standard errors
of the variance components are those of
[`sae::eblupSTFH()`](https://rdrr.io/pkg/sae/man/eblupSTFH.html): the
inverse REML information of the last Fisher-scoring iteration.

**Unsampled domains and missing periods.** Cells whose response is `NA`
(a domain without a sample in some period, or in every period) are not
dropped. The model is defined on the full \\D \times T\\ grid: the SAR
effect \\u_1\\ lives on all \\D\\ domains of `W` and the AR(1) effect
\\u_2\\ on all \\T\\ periods. The variance components are estimated by
REML from the sampled cells, i.e. from the marginal likelihood with the
missing cells integrated out (missing at random), and every cell gets
the best linear unbiased predictor of \\x\_{dt}'\beta + u\_{1d} +
u\_{2dt}\\ given the sampled data (Saei and Chambers, 2005; Benedetti et
al., 2024 for the spatial effect): \$\$\tilde\theta\_{dt} =
x\_{dt}'\tilde\beta + a\_{dt}' V^{-1} (y - X\tilde\beta),\$\$ with
\\a\_{dt}\\ the covariance between the cell and the sampled direct
estimates. A missing period of a sampled domain borrows from the other
periods of the domain (AR(1)) and from the neighbours (SAR); a domain
that is never sampled is predicted from its neighbours only (spatial
kriging). With `model = "S"` the period effect of a missing cell is
independent of the data, so only the spatial part is borrowed.

**MSE.** For sampled cells the analytical MSE is \\g_1 + g_2 + 2 g_3 -
g_4\\ (Marhuenda et al., 2013). For the cells predicted without a direct
estimate it is \\g_1 + g_2 + g_3\\, where \\g_1\\ and \\g_2\\ are the
MSE of the predictor with known variance components and \\g_3\\ comes
from the exact derivative of the predictor in those components (`g4` is
0 for these cells). The further second-order correction of the plug-in
\\g_1\\ is not applied: it is numerically unstable when \\\rho_1\\
approaches \\\pm 1\\ or a variance component is weakly identified, while
in simulations \\g_1 + g_2 + g_3\\ has a relative bias of a few percent
and close to nominal coverage for these cells. The parametric bootstrap
generates the random effects on the full grid and the sampling errors on
the sampled cells, re-estimates the model on the sampled cells and
predicts every cell, so it covers both kinds of cells. The auxiliary
variables must be known for every cell, since they enter the prediction.
Self-benchmarking calibrates the weighted aggregate of the sampled cells
of each group.

## References

1.  Marhuenda, Y., Molina, I., & Morales, D. (2013). Small area
    estimation with spatio-temporal Fay-Herriot models. Computational
    Statistics & Data Analysis, 58, 308-325.

2.  Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley &
    Sons.

3.  Saei, A., & Chambers, R. (2005). Out of sample estimation for small
    areas using area level data. Southampton Statistical Sciences
    Research Institute, Methodology Working Paper M05/11.

4.  Benedetti, R., Piersimoni, F., Pratesi, M., Salvati, N., &
    Suesse, T. (2024). Handling out-of-sample areas to estimate the
    unemployment rate at local labour market areas in Italy.
    International Statistical Review.
    [doi:10.1111/insr.12596](https://doi.org/10.1111/insr.12596)

## See also

[`eblup_sfh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md),
[`hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
(with `time`) for the hierarchical Bayes counterpart.

## Examples

``` r
library(fastsae)

# mys_panel: 42 domains x 11 years; domains 21 and 25 are never sampled (y = NA).
# They stay in the data and in W, and are predicted from their neighbours.
m1 <- eblup_stfh(
  y ~ x1 + x2 + x3,
  data = mys_panel,
  vardir = ~vardir,
  domain = ~area,
  time = ~year,
  W = mys_proxmat,
  model = "ST"
)
#> ℹ Predicting 22 cells without a direct estimate: 2 unsampled domains (spatial borrowing) and 0 missing periods of sampled domains.
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_stfh(formula = y ~ x1 + x2 + x3, vardir = ~vardir, data = mys_panel,
#> domain = ~area, time = ~year, W = mys_proxmat, model = "ST")
#> 
#> ✔ Convergence: Yes (in 11 iterations)
#> Model: Spatio-Temporal Fay-Herriot (SAR + AR(1) time effects)
#> Method: REML
#> Cells without a direct estimate: 22 of 462 (predicted from the model)
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  4.4442780  0.4445132  9.9980779 0.0000
#> x1          -0.0008817  0.0060240 -0.1463646 0.8836
#> x2           0.0254871  0.0277726  0.9177081 0.3588
#> x3           0.0163716  0.0256365  0.6386035 0.5231
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain time        y    eblup    vardir        mse      rse mse_pb
#> 1      1 2016 8.260507 8.377887 0.6618838 0.07413741 3.250006     NA
#> 2      2 2016 7.483862 7.669905 0.8374691 0.08422617 3.783847     NA
#> 3      3 2016 5.553658 5.561797 0.8822257 0.08902309 5.364586     NA
#> 4      4 2016 3.426513 4.008761 0.6581716 0.07340455 6.758515     NA
#> 5      5 2016 6.457152 6.361538 1.2788021 0.12313756 5.516113     NA
#> 6      6 2016 3.608120 4.050605 0.3878004 0.05601510 5.842956     NA
#>   random_effect_u1 random_effect_u2 sampled
#> 1        3.1977947                0    TRUE
#> 2        2.7040777                0    TRUE
#> 3        0.6227519                0    TRUE
#> 4       -1.5624851                0    TRUE
#> 5        0.4308669                0    TRUE
#> 6       -1.2267558                0    TRUE
#> ... and 456 more rows.
#> 

# Analytical MSE for every cell, sampled or not
m2 <- eblup_stfh(
  y ~ x1 + x2 + x3,
  data = mys_panel,
  vardir = ~vardir,
  domain = ~area,
  time = ~year,
  W = mys_proxmat,
  model = "ST",
  mse_type = "analytical"
)
#> ℹ Predicting 22 cells without a direct estimate: 2 unsampled domains (spatial borrowing) and 0 missing periods of sampled domains.
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_stfh(formula = y ~ x1 + x2 + x3, vardir = ~vardir, data = mys_panel,
#> domain = ~area, time = ~year, W = mys_proxmat, model = "ST", mse_type =
#> "analytical")
#> 
#> ✔ Convergence: Yes (in 11 iterations)
#> Model: Spatio-Temporal Fay-Herriot (SAR + AR(1) time effects)
#> Method: REML
#> Cells without a direct estimate: 22 of 462 (predicted from the model)
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  4.4442780  0.4445132  9.9980779 0.0000
#> x1          -0.0008817  0.0060240 -0.1463646 0.8836
#> x2           0.0254871  0.0277726  0.9177081 0.3588
#> x3           0.0163716  0.0256365  0.6386035 0.5231
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain time        y    eblup    vardir        mse      rse mse_pb
#> 1      1 2016 8.260507 8.377887 0.6618838 0.07413741 3.250006     NA
#> 2      2 2016 7.483862 7.669905 0.8374691 0.08422617 3.783847     NA
#> 3      3 2016 5.553658 5.561797 0.8822257 0.08902309 5.364586     NA
#> 4      4 2016 3.426513 4.008761 0.6581716 0.07340455 6.758515     NA
#> 5      5 2016 6.457152 6.361538 1.2788021 0.12313756 5.516113     NA
#> 6      6 2016 3.608120 4.050605 0.3878004 0.05601510 5.842956     NA
#>   random_effect_u1 random_effect_u2 sampled
#> 1        3.1977947                0    TRUE
#> 2        2.7040777                0    TRUE
#> 3        0.6227519                0    TRUE
#> 4       -1.5624851                0    TRUE
#> 5        0.4308669                0    TRUE
#> 6       -1.2267558                0    TRUE
#> ... and 456 more rows.
#> 
head(m2$df_eblup[!m2$df_eblup$sampled, ])
#>     domain time  y    eblup vardir      mse      rse mse_pb random_effect_u1
#> 21      21 2016 NA 6.356729     NA 2.389320 24.31664     NA         1.624999
#> 25      25 2016 NA 3.225462     NA 2.405876 48.08889     NA        -1.771180
#> 63      21 2017 NA 6.367907     NA 2.385248 24.25326     NA         1.624999
#> 67      25 2017 NA 3.213502     NA 2.416602 48.37534     NA        -1.771180
#> 105     21 2018 NA 6.361026     NA 2.389905 24.30319     NA         1.624999
#> 109     25 2018 NA 3.174470     NA 2.438947 49.19602     NA        -1.771180
#>     random_effect_u2 sampled
#> 21                 0   FALSE
#> 25                 0   FALSE
#> 63                 0   FALSE
#> 67                 0   FALSE
#> 105                0   FALSE
#> 109                0   FALSE
```
