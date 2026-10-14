# Hierarchical Bayes for Unit-Level Small Area Estimation

Estimates small area parameters using unit-level Hierarchical Bayesian
models (Battese-Harter-Fuller model and generalized linear mixed models)
with Integrated Nested Laplace Approximations (INLA). Supports Gaussian,
Binomial (logistic regression for binary unit responses), and Poisson
likelihoods, with optional spatial random effect structures (BYM2,
Besag) and finite population adjustments.

## Usage

``` r
hb_unit(
  formula,
  unit_data,
  Xpop,
  domain_var,
  popsize_var = NULL,
  family = c("gaussian", "binomial", "poisson"),
  spatial = c("none", "bym2", "besag"),
  W = NULL,
  popnmean_xpop = NULL,
  engine = "inla",
  print_result = TRUE,
  ...
)
```

## Arguments

- formula:

  An object of class `formula` specifying the unit-level fixed-effects
  model (e.g., `y ~ x1 + x2`).

- unit_data:

  A `data.frame` containing the unit-level survey sample data.

- Xpop:

  A `data.frame` containing auxiliary population information for all
  domains (either domain population means or unit-level population
  records). Must contain the domain identifier specified by
  `domain_var`.

- domain_var:

  Character string specifying the column name for the domain identifier
  in both `unit_data` and `Xpop`.

- popsize_var:

  Optional character string specifying the column name for domain
  population sizes (\\N_d\\) in `Xpop`. If provided, finite population
  adjustment \\\hat{\bar{Y}}\_d = f_d \bar{y}\_{d,s} + (1 - f_d)
  \hat{\mu}\_{d,r}\\ is computed, where \\f_d = n_d / N_d\\. If `NULL`,
  the superpopulation expectation \\\hat{\bar{Y}}\_d = \bar{X}\_d^\top
  \hat{\beta} + \hat{u}\_d\\ is reported.

- family:

  Character string specifying the response likelihood. Options:

  - `"gaussian"`: Continuous response (Battese-Harter-Fuller model,
    default).

  - `"binomial"`: Binary / Bernoulli response (0 or 1, e.g. poverty
    indicator).

  - `"poisson"`: Count response (non-negative integer).

- spatial:

  Character string specifying the spatial random effect structure across
  domains:

  - `"none"`: Non-spatial independent and identically distributed (IID)
    domain effects (default).

  - `"bym2"`: Scaled Besag-York-Mollié 2 spatial model (Riebler et al.,
    2016).

  - `"besag"`: Intrinsic Conditional Autoregressive (ICAR) spatial
    model.

- W:

  Proximity or spatial adjacency matrix. Can be a square `matrix`,
  `Matrix`, or `spdep` `nb` or `listw` object. Dimensions must match the
  total number of unique domains in `Xpop`. Required when
  `spatial != "none"`.

- popnmean_xpop:

  Optional matrix or data frame of auxiliary population means per
  domain. If `NULL` (default), means are extracted directly from `Xpop`.

- engine:

  Estimation engine: `"inla"` (default) or an
  [`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
  object holding the approximation strategy, the priors (`prior_prec`,
  default `list(prior = "loggamma", param = c(0.01, 0.01))`; `prior_phi`
  for BYM2), `scale_model` and further
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) arguments.

- print_result:

  Logical. If `TRUE` (default), prints a summary of results.

- ...:

  Deprecated. `strategy`, `scale_model`, `prior_prec`, `prior_phi` and
  arguments of [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html)
  are now given to
  [`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md);
  passed here they are forwarded to it with a warning.

## Value

An object of class `c("fastsae_hb_unit", "fastsae_hb", "fastsae")`
containing:

- `df_hb`: Data frame with domain estimates, including:

  - `domain`: Domain identifier.

  - `hb`: Estimated domain mean or proportion.

  - `linear_pred`: Posterior mean of domain linear predictor.

  - `sd`: Posterior standard deviation (standard error).

  - `mse`: Posterior Mean Squared Error (\\\text{sd}^2\\).

  - `rse`: Relative Standard Error (%).

  - `ci_lower`: 2.5% quantile of posterior credible interval.

  - `ci_upper`: 97.5% quantile of posterior credible interval.

  - `samp_size`: Sample size (\\n_d\\) in the domain (0 for unsampled).

  - `pop_size`: Population size (\\N_d\\), if `popsize_var` is provided.

  - `estimated_total`: Estimated domain population total (if
    `popsize_var` is provided).

  - `sample_mean`: Direct sample mean (\\\bar{y}\_{d,s}\\).

  - `random_effect`: Posterior mean of domain random effect
    (\\\hat{u}\_d\\).

- `estcoef`: Data frame of estimated regression coefficients (posterior
  mean, sd, z-value, p-value, credible intervals).

- `hyperpar`: Data frame of hyperparameter posterior estimates.

- `random_effect_var`: Estimated domain random effect variance
  (\\\sigma_u^2\\).

- `residual_var`: Estimated residual variance (\\\sigma_e^2\\, for
  Gaussian).

- `phi`: Estimated spatial variance proportion (for BYM2).

- `goodness`: Model fit metrics (DIC, WAIC, Marginal Log-Likelihood).

- `family`: Response family used.

- `spatial`: Spatial model type used.

- `unsampled_domains`: Character vector of unsampled domain identifiers.

- `fit`: Raw fitted INLA model object.

- `engine`, `engine_spec`: Engine name (`"inla"`) and its settings.

- `call`: Matched function call.

## References

1.  Battese, G. E., Harter, R. M., & Fuller, W. A. (1988). An
    error-components model for prediction of county crop areas using
    survey and satellite data. *Journal of the American Statistical
    Association*, 83(401), 28-36.

2.  Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd
    ed.). John Wiley & Sons.

3.  Riebler, A., Sørbye, S. H., Simpson, D., & Rue, H. (2016). An
    intuitive Bayesian spatial model for disease mapping that accounts
    for scaling. *Statistical Methods in Medical Research*, 25(4),
    1145-1165.

4.  Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
    inference for latent Gaussian models by using integrated nested
    Laplace approximations. *Journal of the Royal Statistical Society:
    Series B*, 71(2), 319-392.

## Examples

``` r
# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
  library(fastsae)
  data(cornsoybean)
  data(cornsoybeanmeans)

  # Prepare auxiliary population means
  df_pop <- cornsoybeanmeans
  names(df_pop)[names(df_pop) == "MeanCornPixPerSeg"] <- "CornPix"
  names(df_pop)[names(df_pop) == "MeanSoyBeansPixPerSeg"] <- "SoyBeansPix"

  # Unit survey data
  df_sample <- cornsoybean
  names(df_sample)[names(df_sample) == "County"] <- "CountyIndex"

  # Fit unit-level Hierarchical Bayes model
  fit_hb_u <- hb_unit(
    formula = CornHec ~ CornPix + SoyBeansPix,
    unit_data = df_sample,
    Xpop = df_pop,
    domain_var = "CountyIndex",
    popsize_var = "PopnSegments"
  )
  print(fit_hb_u)
}
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_unit(formula = CornHec ~ CornPix + SoyBeansPix, unit_data = df_sample, Xpop
#> = df_pop, domain_var = "CountyIndex", popsize_var = "PopnSegments")
#> 
#> ✔ Convergence: Yes (in - iterations)
#> Model: HB-UNIT-GAUSSIAN (Non-spatial)
#> Engine: inla
#> Method: INLA (laplace)
#> Random effect variance (sigma2_u): 67.48867 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  1.7657e+01  3.1448e+01  5.6146e-01  5.7448e-01 -4.4280e+01
#> CornPix      3.6613e-01  6.5944e-02  5.5522e+00  2.8212e-08  2.3597e-01
#> SoyBeansPix -2.9102e-02  6.8518e-02 -4.2473e-01  6.7103e-01 -1.6411e-01
#>             ci_upper
#> (Intercept)  79.6880
#> CornPix       0.4959
#> SoyBeansPix   0.1060
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred    vardir       sd      mse      rse
#> 1      1 165.7600 122.2845    122.2048 292.00611 7.792395 60.72141 6.372348
#> 2      2  96.3200 123.3396    123.3875 292.00611 7.644423 58.43721 6.197867
#> 3      3  76.0800 113.2339    113.3281 292.00611 8.417644 70.85672 7.433857
#> 4      4 150.8900 115.2436    115.0746 146.00305 7.113667 50.60426 6.172724
#> 5      5 158.6233 135.5099    135.3863  97.33537 8.111452 65.79566 5.985875
#> 6      6 102.5233 108.1532    108.1830  97.33537 6.770167 45.83517 6.259793
#>    ci_lower ci_upper samp_size sample_mean random_effect pop_size
#> 1 107.01142 137.5576         1    165.7600      2.117395      545
#> 2 108.35650 138.3226         1     96.3200      1.485219      566
#> 3  96.73529 129.7325         1     76.0800     -4.671558      394
#> 4 101.30078 129.1864         2    150.8900     -2.522388      424
#> 5 119.61145 151.4083         3    158.6233      7.011414      564
#> 6  94.88369 121.4227         3    102.5233      3.613401      570
#>   estimated_total
#> 1        66645.06
#> 2        69810.20
#> 3        44614.14
#> 4        48863.27
#> 5        76427.58
#> 6        61647.33
#> ... and 6 more rows.
#> 
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_unit(formula = CornHec ~ CornPix + SoyBeansPix, unit_data = df_sample, Xpop
#> = df_pop, domain_var = "CountyIndex", popsize_var = "PopnSegments")
#> 
#> ✔ Convergence: Yes (in - iterations)
#> Model: HB-UNIT-GAUSSIAN (Non-spatial)
#> Engine: inla
#> Method: INLA (laplace)
#> Random effect variance (sigma2_u): 67.48867 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  1.7657e+01  3.1448e+01  5.6146e-01  5.7448e-01 -4.4280e+01
#> CornPix      3.6613e-01  6.5944e-02  5.5522e+00  2.8212e-08  2.3597e-01
#> SoyBeansPix -2.9102e-02  6.8518e-02 -4.2473e-01  6.7103e-01 -1.6411e-01
#>             ci_upper
#> (Intercept)  79.6880
#> CornPix       0.4959
#> SoyBeansPix   0.1060
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred    vardir       sd      mse      rse
#> 1      1 165.7600 122.2845    122.2048 292.00611 7.792395 60.72141 6.372348
#> 2      2  96.3200 123.3396    123.3875 292.00611 7.644423 58.43721 6.197867
#> 3      3  76.0800 113.2339    113.3281 292.00611 8.417644 70.85672 7.433857
#> 4      4 150.8900 115.2436    115.0746 146.00305 7.113667 50.60426 6.172724
#> 5      5 158.6233 135.5099    135.3863  97.33537 8.111452 65.79566 5.985875
#> 6      6 102.5233 108.1532    108.1830  97.33537 6.770167 45.83517 6.259793
#>    ci_lower ci_upper samp_size sample_mean random_effect pop_size
#> 1 107.01142 137.5576         1    165.7600      2.117395      545
#> 2 108.35650 138.3226         1     96.3200      1.485219      566
#> 3  96.73529 129.7325         1     76.0800     -4.671558      394
#> 4 101.30078 129.1864         2    150.8900     -2.522388      424
#> 5 119.61145 151.4083         3    158.6233      7.011414      564
#> 6  94.88369 121.4227         3    102.5233      3.613401      570
#>   estimated_total
#> 1        66645.06
#> 2        69810.20
#> 3        44614.14
#> 4        48863.27
#> 5        76427.58
#> 6        61647.33
#> ... and 6 more rows.
#> 
# }
```
