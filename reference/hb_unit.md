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
#> Random effect variance (sigma2_u): 64.90171 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  1.7745e+01  3.1594e+01  5.6165e-01  5.7436e-01 -4.4460e+01
#> CornPix      3.6560e-01  6.6219e-02  5.5210e+00  3.3702e-08  2.3477e-01
#> SoyBeansPix -2.8937e-02  6.8834e-02 -4.2039e-01  6.7420e-01 -1.6452e-01
#>             ci_upper
#> (Intercept)  80.0806
#> CornPix       0.4958
#> SoyBeansPix   0.1068
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred   vardir       sd      mse      rse
#> 1      1 165.7600 122.0304    121.9506 303.9354 7.438982 55.33845 6.096008
#> 2      2  96.3200 123.1504    123.1983 303.9354 7.275520 52.93320 5.907835
#> 3      3  76.0800 113.6687    113.7632 303.9354 8.075731 65.21742 7.104618
#> 4      4 150.8900 115.4670    115.2990 151.9677 6.806845 46.33314 5.895058
#> 5      5 158.6233 134.7432    134.6155 101.3118 8.108085 65.74104 6.017435
#> 6      6 102.5233 107.7689    107.7967 101.3118 6.604518 43.61965 6.128405
#>    ci_lower ci_upper samp_size sample_mean random_effect pop_size
#> 1 107.44998 136.6108         1    165.7600      1.687065      545
#> 2 108.89035 137.4104         1     96.3200      1.323445      566
#> 3  97.84032 129.4972         1     76.0800     -3.687659      394
#> 4 102.12556 128.8084         2    150.8900     -2.654040      424
#> 5 118.85136 150.6351         3    158.6233      5.831867      564
#> 6  94.82409 120.7138         3    102.5233      2.844545      570
#>   estimated_total
#> 1        66506.56
#> 2        69703.11
#> 3        44785.49
#> 4        48958.00
#> 5        75995.17
#> 6        61428.30
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
#> Random effect variance (sigma2_u): 64.90171 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  1.7745e+01  3.1594e+01  5.6165e-01  5.7436e-01 -4.4460e+01
#> CornPix      3.6560e-01  6.6219e-02  5.5210e+00  3.3702e-08  2.3477e-01
#> SoyBeansPix -2.8937e-02  6.8834e-02 -4.2039e-01  6.7420e-01 -1.6452e-01
#>             ci_upper
#> (Intercept)  80.0806
#> CornPix       0.4958
#> SoyBeansPix   0.1068
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred   vardir       sd      mse      rse
#> 1      1 165.7600 122.0304    121.9506 303.9354 7.438982 55.33845 6.096008
#> 2      2  96.3200 123.1504    123.1983 303.9354 7.275520 52.93320 5.907835
#> 3      3  76.0800 113.6687    113.7632 303.9354 8.075731 65.21742 7.104618
#> 4      4 150.8900 115.4670    115.2990 151.9677 6.806845 46.33314 5.895058
#> 5      5 158.6233 134.7432    134.6155 101.3118 8.108085 65.74104 6.017435
#> 6      6 102.5233 107.7689    107.7967 101.3118 6.604518 43.61965 6.128405
#>    ci_lower ci_upper samp_size sample_mean random_effect pop_size
#> 1 107.44998 136.6108         1    165.7600      1.687065      545
#> 2 108.89035 137.4104         1     96.3200      1.323445      566
#> 3  97.84032 129.4972         1     76.0800     -3.687659      394
#> 4 102.12556 128.8084         2    150.8900     -2.654040      424
#> 5 118.85136 150.6351         3    158.6233      5.831867      564
#> 6  94.82409 120.7138         3    102.5233      2.844545      570
#>   estimated_total
#> 1        66506.56
#> 2        69703.11
#> 3        44785.49
#> 4        48958.00
#> 5        75995.17
#> 6        61428.30
#> ... and 6 more rows.
#> 
# }
```
