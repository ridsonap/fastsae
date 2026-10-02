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
#> Random effect variance (sigma2_u): 32.6526 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  1.7648e+01  3.1139e+01  5.6673e-01  5.7089e-01 -4.3666e+01
#> CornPix      3.6612e-01  6.5299e-02  5.6069e+00  2.0599e-08  2.3729e-01
#> SoyBeansPix -2.9062e-02  6.7842e-02 -4.2837e-01  6.6838e-01 -1.6271e-01
#>             ci_upper
#> (Intercept)  79.0479
#> CornPix       0.4946
#> SoyBeansPix   0.1046
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred   vardir       sd      mse      rse
#> 1      1 165.7600 122.3343    122.2546 272.2962 7.678957 58.96638 6.277026
#> 2      2  96.3200 123.2725    123.3202 272.2962 7.527643 56.66541 6.106508
#> 3      3  76.0800 113.2017    113.2960 272.2962 8.307519 69.01487 7.338690
#> 4      4 150.8900 115.4470    115.2790 136.1481 6.987242 48.82156 6.052339
#> 5      5 158.6233 135.5072    135.3836  90.7654 8.085578 65.37657 5.966900
#> 6      6 102.5233 108.0618    108.0911  90.7654 6.682154 44.65118 6.183640
#>    ci_lower ci_upper samp_size sample_mean random_effect pop_size
#> 1 107.28355 137.3851         1    165.7600      2.117912      545
#> 2 108.51829 138.0267         1     96.3200      1.486738      566
#> 3  96.91893 129.4844         1     76.0800     -4.683639      394
#> 4 101.75199 129.1420         2    150.8900     -2.519906      424
#> 5 119.65945 151.3549         3    158.6233      6.985716      564
#> 6  94.96479 121.1588         3    102.5233      3.597851      570
#>   estimated_total
#> 1        66672.20
#> 2        69772.22
#> 3        44601.46
#> 4        48949.52
#> 5        76426.05
#> 6        61595.23
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
#> Random effect variance (sigma2_u): 32.6526 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  1.7648e+01  3.1139e+01  5.6673e-01  5.7089e-01 -4.3666e+01
#> CornPix      3.6612e-01  6.5299e-02  5.6069e+00  2.0599e-08  2.3729e-01
#> SoyBeansPix -2.9062e-02  6.7842e-02 -4.2837e-01  6.6838e-01 -1.6271e-01
#>             ci_upper
#> (Intercept)  79.0479
#> CornPix       0.4946
#> SoyBeansPix   0.1046
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred   vardir       sd      mse      rse
#> 1      1 165.7600 122.3343    122.2546 272.2962 7.678957 58.96638 6.277026
#> 2      2  96.3200 123.2725    123.3202 272.2962 7.527643 56.66541 6.106508
#> 3      3  76.0800 113.2017    113.2960 272.2962 8.307519 69.01487 7.338690
#> 4      4 150.8900 115.4470    115.2790 136.1481 6.987242 48.82156 6.052339
#> 5      5 158.6233 135.5072    135.3836  90.7654 8.085578 65.37657 5.966900
#> 6      6 102.5233 108.0618    108.0911  90.7654 6.682154 44.65118 6.183640
#>    ci_lower ci_upper samp_size sample_mean random_effect pop_size
#> 1 107.28355 137.3851         1    165.7600      2.117912      545
#> 2 108.51829 138.0267         1     96.3200      1.486738      566
#> 3  96.91893 129.4844         1     76.0800     -4.683639      394
#> 4 101.75199 129.1420         2    150.8900     -2.519906      424
#> 5 119.65945 151.3549         3    158.6233      6.985716      564
#> 6  94.96479 121.1588         3    102.5233      3.597851      570
#>   estimated_total
#> 1        66672.20
#> 2        69772.22
#> 3        44601.46
#> 4        48949.52
#> 5        76426.05
#> 6        61595.23
#> ... and 6 more rows.
#> 
# }
```
