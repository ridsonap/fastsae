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
  strategy = c("laplace", "simplified.laplace"),
  scale_model = TRUE,
  prior_prec = list(prior = "loggamma", param = c(0.01, 0.01)),
  prior_phi = list(prior = "pc", param = c(0.5, 0.5)),
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

- strategy:

  INLA approximation strategy: `"laplace"` (default, highest accuracy)
  or `"simplified.laplace"` (faster).

- scale_model:

  Logical. If `TRUE` (default), scales the spatial graph so the marginal
  variance of the structured effect is approximately 1 (recommended for
  BYM2/Besag).

- prior_prec:

  List specifying the prior for domain random effect precision. Default
  is `list(prior = "loggamma", param = c(0.01, 0.01))`.

- prior_phi:

  List specifying the PC-prior for the spatial mixing parameter \\\phi\\
  in the BYM2 model. Default is
  `list(prior = "pc", param = c(0.5, 0.5))`.

- print_result:

  Logical. If `TRUE` (default), prints a summary of results.

- ...:

  Additional arguments passed to
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html).

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
#> Method: INLA (laplace)
#> Random effect variance (sigma2_u): 1.625993 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  1.8290e+01  3.1454e+01  5.8147e-01  5.6092e-01 -4.3687e+01
#> CornPix      3.6313e-01  6.5730e-02  5.5245e+00  3.3039e-08  2.3352e-01
#> SoyBeansPix -2.8599e-02  6.8534e-02 -4.1729e-01  6.7646e-01 -1.6358e-01
#>             ci_upper
#> (Intercept)  80.2739
#> CornPix       0.4926
#> SoyBeansPix   0.1065
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred    vardir       sd      mse      rse
#> 1      1 165.7600 120.7518    120.6691 290.94037 4.894334 23.95451 4.053219
#> 2      2  96.3200 122.0751    122.1207 290.94037 4.824955 23.28019 3.952447
#> 3      3  76.0800 116.2921    116.3944 290.94037 4.876114 23.77649 4.192989
#> 4      4 150.8900 117.0096    116.8490 145.47018 4.765913 22.71393 4.073095
#> 5      5 158.6233 131.1504    131.0035  96.98012 4.930084 24.30573 3.759107
#> 6      6 102.5233 105.8584    105.8761  96.98012 4.900004 24.01004 4.628828
#>    ci_lower ci_upper samp_size sample_mean random_effect pop_size
#> 1 111.15891 130.3447         1    165.7600     0.5692085      545
#> 2 112.61823 131.5321         1     96.3200     0.3799528      566
#> 3 106.73490 125.8493         1     76.0800    -1.2135491      394
#> 4 107.66843 126.3508         2    150.8900    -0.7011030      424
#> 5 121.48745 140.8134         3    158.6233     2.6425665      564
#> 6  96.25441 115.4624         3    102.5233     1.3097369      570
#>   estimated_total
#> 1        65809.73
#> 2        69094.53
#> 3        45819.08
#> 4        49612.08
#> 5        73968.83
#> 6        60339.30
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
#> Method: INLA (laplace)
#> Random effect variance (sigma2_u): 1.625993 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  1.8290e+01  3.1454e+01  5.8147e-01  5.6092e-01 -4.3687e+01
#> CornPix      3.6313e-01  6.5730e-02  5.5245e+00  3.3039e-08  2.3352e-01
#> SoyBeansPix -2.8599e-02  6.8534e-02 -4.1729e-01  6.7646e-01 -1.6358e-01
#>             ci_upper
#> (Intercept)  80.2739
#> CornPix       0.4926
#> SoyBeansPix   0.1065
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred    vardir       sd      mse      rse
#> 1      1 165.7600 120.7518    120.6691 290.94037 4.894334 23.95451 4.053219
#> 2      2  96.3200 122.0751    122.1207 290.94037 4.824955 23.28019 3.952447
#> 3      3  76.0800 116.2921    116.3944 290.94037 4.876114 23.77649 4.192989
#> 4      4 150.8900 117.0096    116.8490 145.47018 4.765913 22.71393 4.073095
#> 5      5 158.6233 131.1504    131.0035  96.98012 4.930084 24.30573 3.759107
#> 6      6 102.5233 105.8584    105.8761  96.98012 4.900004 24.01004 4.628828
#>    ci_lower ci_upper samp_size sample_mean random_effect pop_size
#> 1 111.15891 130.3447         1    165.7600     0.5692085      545
#> 2 112.61823 131.5321         1     96.3200     0.3799528      566
#> 3 106.73490 125.8493         1     76.0800    -1.2135491      394
#> 4 107.66843 126.3508         2    150.8900    -0.7011030      424
#> 5 121.48745 140.8134         3    158.6233     2.6425665      564
#> 6  96.25441 115.4624         3    102.5233     1.3097369      570
#>   estimated_total
#> 1        65809.73
#> 2        69094.53
#> 3        45819.08
#> 4        49612.08
#> 5        73968.83
#> 6        60339.30
#> ... and 6 more rows.
#> 
# }
```
