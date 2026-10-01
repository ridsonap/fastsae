# Hierarchical Bayes for Area-Level Small Area Estimation

Estimates small area parameters using area-level Hierarchical Bayesian
models under various distributions (Gaussian, Binomial, Poisson,
Negative Binomial, Beta, Gamma) and spatial random effect structures
(non-spatial, BYM2, BYM, Besag/ICAR, Leroux) using Integrated Nested
Laplace Approximations (INLA) or frequentist Laplace Approximation.

## Usage

``` r
hb_area(
  formula,
  data,
  domain = NULL,
  time = NULL,
  family = c("gaussian", "binomial", "poisson", "nbinomial", "beta", "gamma"),
  spatial = c("none", "bym2", "bym", "besag", "generic1", "slm"),
  temporal = c("none", "rw1", "rw2", "ar1", "iid"),
  st_interaction = c("none", "domain-specific", "separable", "type1", "type2", "type3",
    "type4"),
  W = NULL,
  vardir = NULL,
  trials = NULL,
  exposure = NULL,
  method = c("inla", "laplace"),
  strategy = c("laplace", "simplified.laplace", "gaussian"),
  link = NULL,
  scale_model = TRUE,
  prior_prec = list(prior = "pc.prec", param = c(1, 0.01)),
  prior_phi = list(prior = "pc", param = c(0.5, 0.5)),
  prior_rho = NULL,
  prior_prec_time = NULL,
  prior_rho_time = NULL,
  print_result = TRUE,
  ...
)
```

## Arguments

- formula:

  An object of class `formula` specifying the fixed-effects model (e.g.,
  `y ~ x1 + x2`).

- data:

  A `data.frame` containing the area-level data (one row per area).

- domain:

  Vector, column name, or one-sided formula referencing the area/domain
  identifier in `data`. If `NULL`, domains are numbered consecutively.

- time:

  Vector, column name, or one-sided formula referencing the time period
  identifier in `data` (e.g. `time = "year"` or `~year`). Required when
  `temporal != "none"`.

- family:

  Character string specifying the response distribution / likelihood.
  Options:

  - `"gaussian"`: Continuous response (Fay-Herriot model). If `vardir`
    is provided, sampling variances are treated as known. If
    `vardir = NULL`, residual variance is estimated.

  - `"binomial"`: Binary / proportion response. Requires `trials`.

  - `"poisson"`: Count response / disease rate. Supports `exposure`.

  - `"nbinomial"`: Negative Binomial for overdispersed counts.

  - `"beta"`: Continuous proportions strictly in (0, 1). If `vardir` is
    provided, area-specific precision is calculated via Janicki (2020)
    formula \\\phi_i = y_i(1 - y_i)/V_i - 1\\ and injected via INLA's
    `scale` argument. If `trials` is provided, precision is scaled as
    \\n_i - 1\\.

  - `"gamma"`: Skewed positive continuous response.

- spatial:

  Character string specifying the spatial random effect structure:

  - `"none"`: Non-spatial IID area random intercept (default).

  - `"bym2"`: Scaled Besag-York-Mollié 2 model (Riebler et al., 2016),
    decomposing area variance into spatial (\\\phi\\) and unstructured
    components.

  - `"bym"`: Classic Besag-York-Mollié (1991) model.

  - `"besag"`: Intrinsic Conditional Autoregressive (ICAR) model.

  - `"generic1"`: Leroux spatial autoregressive model.

- temporal:

  Character string specifying the temporal random effect structure:

  - `"none"`: No temporal random effect (default cross-sectional model).

  - `"rw1"`: First-order random walk across time periods.

  - `"rw2"`: Second-order random walk across time periods.

  - `"ar1"`: First-order autoregressive process across time periods.

  - `"iid"`: Unstructured independent time effects.

- st_interaction:

  Character string specifying the spatio-temporal structure:

  - `"none"`: Additive main spatial and temporal effects (default).

  - `"domain-specific"`: Domain-specific temporal random walk / AR(1)
    dynamics with shared variance parameter (directly corresponds to the
    tipsae spatio-temporal model).

  - `"separable"`: Grouped dynamic spatial field evolving over time via
    AR(1) / RW (corresponds to classical Spatio-Temporal Fay-Herriot
    models, Marhuenda et al. 2013).

  - `"type1"`: Knorr-Held (2000) Type I interaction (unstructured space
    \\\times\\ unstructured time).

  - `"type2"`: Knorr-Held Type II interaction (unstructured space
    \\\times\\ structured time).

  - `"type3"`: Knorr-Held Type III interaction (structured space
    \\\times\\ unstructured time).

  - `"type4"`: Knorr-Held Type IV interaction (structured space
    \\\times\\ structured time).

- W:

  Proximity or spatial adjacency matrix. Can be a square `matrix`,
  `Matrix`, or `spdep` `nb` or `listw` object. Dimensions must match the
  total number of domains in `data`. Required when `spatial != "none"`.

- vardir:

  Vector, column name, or formula specifying the sampling variances of
  the direct estimator (for `family = "gaussian"` or `family = "beta"`).

- trials:

  Vector, column name, or formula specifying sample sizes / total trials
  per area (for `family = "binomial"`).

- exposure:

  Vector, column name, or formula specifying expected counts or
  population exposure offsets (for `family = "poisson"` or
  `"nbinomial"`).

- method:

  Estimation method: `"inla"` (Bayesian INLA, default) or `"laplace"`
  (Frequentist GLMM with Laplace approximation via `lme4`).

- strategy:

  INLA approximation strategy: `"laplace"` (default, full Laplace
  approximation for highest accuracy), `"simplified.laplace"` (faster
  Taylor approximation), or `"gaussian"`.

- link:

  Optional character string for link function. If `NULL`, default
  canonical link is used (identity for Gaussian, logit for
  Binomial/Beta, log for Poisson/NegBinom/Gamma).

- scale_model:

  Logical. If `TRUE` (default), scales the spatial graph so the marginal
  variance of the structured effect is approximately 1 (recommended for
  BYM2/Besag).

- prior_prec:

  List specifying the prior for random effect precision. Default is a
  Penalized Complexity (PC) prior:
  `list(prior = "pc.prec", param = c(1, 0.01))`.

- prior_phi:

  List specifying the PC-prior for the spatial mixing parameter \\\phi\\
  in the BYM2 model. Default is
  `list(prior = "pc", param = c(0.5, 0.5))`.

- prior_rho:

  Optional list specifying prior for spatial autocorrelation parameter
  (generic1/slm).

- prior_prec_time:

  Optional list specifying the prior for temporal precision. Default is
  a PC prior: `list(prior = "pc.prec", param = c(1, 0.01))`.

- prior_rho_time:

  Optional list specifying prior for temporal autocorrelation parameter
  in AR(1).

- print_result:

  Logical. If `TRUE` (default), prints a summary of results.

- ...:

  Additional arguments passed to
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html).

## Value

An object of class `c("fastsae_hb_area", "fastsae")` containing:

- `df_hb`: Data frame with domain estimates, including `domain`, `time`
  (if specified), observed `y`, predicted `hb`, linear predictor
  `linear_pred`, posterior standard error `sd`, `mse`, relative error
  `rse` (%), 95% credible interval (`ci_lower`, `ci_upper`), and
  `random_effect`.

- `estcoef`: Data frame of estimated regression coefficients.

- `hyperpar`: Data frame of hyperparameter posterior estimates.

- `random_effect_var`: Estimated area random effect variance.

- `random_effect_var_time`: Estimated temporal random effect variance
  (if temporal).

- `phi`: Estimated spatial variance proportion (for BYM2) or spatial
  autocorrelation.

- `rho_time`: Estimated temporal autocorrelation parameter (for AR1).

- `goodness`: Model fit metrics (DIC, WAIC, Marginal Log-Likelihood).

- `family`: Response family used.

- `spatial`: Spatial model type used.

- `temporal`: Temporal model type used.

- `st_interaction`: Spatio-temporal interaction type used.

- `fit`: Raw fitted model object.

- `call`: Matched function call.

## References

1.  Rao, J. N. K., and Molina, I. (2015). *Small Area Estimation*. John
    Wiley & Sons.

2.  Riebler, A., Sørbye, S. H., Simpson, D., and Rue, H. (2016). An
    intuitive Bayesian spatial model for disease mapping that accounts
    for scaling. *Statistical Methods in Medical Research*, 25(4),
    1145-1165.

3.  Marhuenda, Y., Molina, I., and Morales, D. (2013). Small area
    estimation with spatio-temporal Fay-Herriot models. *Computational
    Statistics & Data Analysis*, 58, 308-325.

4.  De Nicolò, S., and Gardini, A. (2024). The R Package tipsae: Tools
    for Mapping Proportions and Indicators on the Unit Interval.
    *Journal of Statistical Software*, 108(1), 1-36.

5.  Rue, H., Martino, S., and Chopin, N. (2009). Approximate Bayesian
    inference for latent Gaussian models by using integrated nested
    Laplace approximations. *Journal of the Royal Statistical Society:
    Series B*, 71(2), 319-392.

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  library(fastsae)

  # 1. Non-spatial Gaussian Fay-Herriot with INLA
  m_norm <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian"
  )

  # 2. Spatial BYM2 Gaussian Fay-Herriot with INLA
  m_bym2 <- hb_area(
    y ~ x1 + x2 + x3,
    data = mys,
    vardir = "vardir",
    family = "gaussian",
    spatial = "bym2",
    W = mys_proxmat
  )
}
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_area(formula = y ~ x1 + x2 + x3, data = mys, family = "gaussian", vardir =
#> "vardir")
#> 
#> ✔ Convergence: Yes (in - iterations)
#> Model: HB-GAUSSIAN (Non-spatial)
#> Random effect variance (sigma2_u): 1.666331 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.04292404  0.69337666  4.38855848  0.00001141  1.69553291
#> x1          -0.00117776  0.00880615 -0.13374311  0.89360572 -0.01878331
#> x2           0.05166286  0.05557978  0.92952611  0.35261650 -0.05677098
#> x3           0.03775957  0.05272564  0.71615189  0.47389757 -0.06709228
#>             ci_upper
#> (Intercept)   4.4274
#> x1            0.0160
#> x2            0.1621
#> x3            0.1405
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.346448    7.346448 0.7431500 0.5522719 10.11577 5.909856
#> 2      2 7.599650 6.505014    6.505014 0.8038211 0.6461283 12.35695 4.956225
#> 3      3 5.514137 5.053750    5.053750 0.7980968 0.6369584 15.79217 3.502536
#> 4      4 3.869326 4.302286    4.302286 0.7107171 0.5051189 16.51952 2.898599
#> 5      5 6.305063 6.322332    6.322332 0.8821486 0.7781861 13.95290 4.589879
#> 6      6 3.926807 4.084926    4.084926 0.5731446 0.3284948 14.03072 2.958638
#>   ci_upper random_effect    vardir
#> 1 8.822792    2.67523343 0.6618838
#> 2 8.107578    2.29154169 0.8374691
#> 3 6.634331    0.90446830 0.8822257
#> 4 5.687250   -1.16084235 0.6581716
#> 5 8.055282   -0.02614579 1.2788021
#> 6 5.206914   -0.72016021 0.3878004
#> ... and 36 more rows.
#> 
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_area(formula = y ~ x1 + x2 + x3, data = mys, family = "gaussian", spatial =
#> "bym2", W = mys_proxmat, vardir = "vardir")
#> 
#> ✔ Convergence: Yes (in - iterations)
#> Model: HB-GAUSSIAN (BYM2)
#> Random effect variance (sigma2_u): 1.666579 
#> Spatial autocorrelation (rho): 0.4816 
#> Spatial mixing fraction (phi): 0.4816 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  3.0436e+00  6.7545e-01  4.5060e+00  6.6051e-06  1.7297e+00
#> x1          -1.1809e-03  8.7848e-03 -1.3443e-01  8.9307e-01 -1.8693e-02
#> x2           5.1685e-02  5.5474e-02  9.3170e-01  3.5149e-01 -5.6614e-02
#> x3           3.7716e-02  5.2623e-02  7.1671e-01  4.7355e-01 -6.6717e-02
#>             ci_upper
#> (Intercept)   4.3885
#> x1            0.0159
#> x2            0.1617
#> x3            0.1403
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.348490    7.348490 0.7396469 0.5470775 10.06529 5.914594
#> 2      2 7.599650 6.506873    6.506873 0.8002646 0.6404234 12.29876 4.960632
#> 3      3 5.514137 5.054934    5.054934 0.7974642 0.6359491 15.77596 3.503539
#> 4      4 3.869326 4.301745    4.301745 0.7103568 0.5046068 16.51322 2.899893
#> 5      5 6.305063 6.322403    6.322403 0.8823432 0.7785295 13.95582 4.589801
#> 6      6 3.926807 4.084742    4.084742 0.5731888 0.3285454 14.03244 2.958605
#>   ci_upper random_effect    vardir
#> 1 8.814090    2.67927686 0.6618838
#> 2 8.097903    2.29516722 0.8374691
#> 3 6.632351    0.90612419 0.8822257
#> 4 5.686877   -1.16247193 0.6581716
#> 5 8.055458   -0.02605096 1.2788021
#> 6 5.207001   -0.72119734 0.3878004
#> ... and 36 more rows.
#> 
# }
```
