# Hierarchical Bayes for Unit-Level Small Area Estimation

Estimates small area parameters using unit-level Hierarchical Bayesian
models (Battese-Harter-Fuller model and generalized linear mixed models)
with Integrated Nested Laplace Approximations (INLA) or by full MCMC
([`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md),
NumPyro NUTS). Supports Gaussian, Binomial (logistic regression for
binary unit responses), and Poisson likelihoods, with optional spatial
random effect structures (BYM2, Besag) and finite population
adjustments.

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

  The domain means of `"binomial"` and `"poisson"` are approximate and
  biased; see Details.

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
  object, or `"jax"` or a
  [`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
  object for MCMC (see the section Engines). An
  [`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
  holds the approximation strategy, the priors (`prior_prec`, default
  `list(prior = "loggamma", param = c(0.01, 0.01))`, vague on any
  response scale; `prior_phi` for BYM2), `scale_model` and further
  [`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) arguments.
  The Gaussian model is fitted in INLA's `"classic"` mode, the binomial
  and Poisson models in INLA's default `"compact"` mode (for them
  classic was much slower: 72 s against 2.6 s for 4000 binary units);
  `inla_engine(inla.mode = )` selects a mode.

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

- `random_effect_var`: Posterior mean of the domain random effect
  variance (\\\sigma_u^2\\). With few domains and the vague default
  prior its posterior is very skewed and INLA's summary of it is not
  accurate (on `cornsoybean`, 12 counties: about 70 against an exact
  posterior mean of 41 and median of 7, while the county means and SDs
  agree with the exact posterior). Do not read it as a REML estimate;
  the domain estimates in `df_hb` are the reliable output.

- `residual_var`: Estimated residual variance (\\\sigma_e^2\\, for
  Gaussian).

- `phi`: Estimated spatial variance proportion (for BYM2).

- `goodness`: Model fit metrics (DIC, WAIC, Marginal Log-Likelihood).

- `family`: Response family used.

- `spatial`: Spatial model type used.

- `unsampled_domains`: Character vector of unsampled domain identifiers.

- `fit`: Raw fitted INLA model object (the NumPyro output for
  [`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)).

- `engine`, `engine_spec`: Engine name (`"inla"` or `"jax"`) and its
  settings.

- `draws`, `rhat_max`, `device`, `elapsed_seconds`: With
  [`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md),
  the domains x draws matrix of the domain means (if
  `keep_draws = TRUE`), the largest split R-hat (`convergence` is
  `rhat_max < 1.05`), the device used and the sampling time.

- `call`: Matched function call.

## Details

**Non-Gaussian families.** The non-sampled part of a domain mean is
predicted at the mean covariates of the non-sampled units, \\g^{-1}(\bar
x\_{d,r}^\top\beta + u_d)\\. For the identity link of `"gaussian"` this
is the mean of \\x_j^\top\beta + u_d\\ over the units. For the logit and
log links of `"binomial"` and `"poisson"` it is not: the mean of
\\g^{-1}(x_j^\top\beta + u_d)\\ over the units differs from it whenever
the covariates vary within the domain (for the log link it is larger by
the factor \\E\[\exp((x_j - \bar x)^\top\beta)\] \> 1\\). In a
simulation with one normal covariate of within-domain variance 1, the
domain means of the sampled domains were biased by -8% (binomial, slope
1.2) and -16% (Poisson, slope 0.6), and the 95% credible intervals
covered the true means in 90% and 78% of the cases. Unbiased estimation
needs the covariates of every population unit, as in the empirical best
predictor of Molina and Rao (2010); a unit-level `Xpop` is currently
averaged per domain. A warning is given for these families.

## Engines

Both engines fit the same model from the sampled units and return the
same structure. With
[`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
the coefficients have normal priors (standard deviation 10, or 2.5 for
`"binomial"`) and the standard deviations of the domain effect and, for
`"gaussian"`, of the unit errors half-normal(1) priors, on the
standardised response scale for `"gaussian"`;
[`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
puts the vague log-gamma prior on the precisions. The domain means, SDs
and 2.5%/97.5% quantiles are computed from the MCMC draws; with
`popsize_var` each draw of the Gaussian domain mean includes the mean
error of the non-sampled units, \\N(0, \sigma_e^2 / (N_d - n_d))\\ for
that draw of \\\sigma_e^2\\ (the INLA engine adds its posterior mean to
the variance and gives a normal interval).

On simulated data (30 domains, three unsampled) a Stan fit of the same
model and priors gave domain means within 0.04 posterior SD and SDs
within 5% for the Gaussian (also with BYM2), Bernoulli and Poisson
models. With their own priors, `hbsae::fSAE(method = "HB")` (Gaussian,
numerical integration over the variance ratio, flat priors) and the INLA
engine gave domain means within 0.1 and 0.06 posterior SD. With few
domains the prior of the domain variance matters: on `cornsoybean` (12
counties) the county SDs were 1.1-1.4 times those of INLA, whose
log-gamma prior favours small variances, and 0.8-1.0 times those of
hbsae.

## References

1.  Molina, I., & Rao, J. N. K. (2010). Small area estimation of poverty
    indicators. *Canadian Journal of Statistics*, 38(3), 369-385.

2.  Battese, G. E., Harter, R. M., & Fuller, W. A. (1988). An
    error-components model for prediction of county crop areas using
    survey and satellite data. *Journal of the American Statistical
    Association*, 83(401), 28-36.

3.  Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd
    ed.). John Wiley & Sons.

4.  Riebler, A., Sørbye, S. H., Simpson, D., & Rue, H. (2016). An
    intuitive Bayesian spatial model for disease mapping that accounts
    for scaling. *Statistical Methods in Medical Research*, 25(4),
    1145-1165.

5.  Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
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
#> ✔ Convergence: Yes (INLA mode status 0)
#> Model: HB-UNIT-GAUSSIAN (Non-spatial)
#> Engine: inla
#> Method: INLA (laplace)
#> Random effect variance (sigma2_u): 90.54501 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  1.7975e+01  3.1524e+01  5.7020e-01  5.6854e-01 -4.4131e+01
#> CornPix      3.6432e-01  6.5991e-02  5.5208e+00  3.3749e-08  2.3409e-01
#> SoyBeansPix -2.8600e-02  6.8675e-02 -4.1645e-01  6.7708e-01 -1.6393e-01
#>             ci_upper
#> (Intercept)  80.1491
#> CornPix       0.4942
#> SoyBeansPix   0.1068
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred    vardir       sd      mse      rse
#> 1      1 165.7600 121.3256    121.2465 294.62568 6.358679 40.43280 5.241004
#> 2      2  96.3200 122.6207    122.6693 294.62568 6.192291 38.34447 5.049954
#> 3      3  76.0800 114.9777    115.0714 294.62568 6.930166 48.02720 6.027398
#> 4      4 150.8900 116.1428    115.9772 147.31284 5.921748 35.06710 5.098679
#> 5      5 158.6233 132.7227    132.5842  98.20856 7.283982 53.05640 5.488121
#> 6      6 102.5233 106.7526    106.7750  98.20856 5.954208 35.45259 5.577578
#>    ci_lower ci_upper samp_size sample_mean random_effect pop_size
#> 1 108.86258 133.7886         1    165.7600     1.1374532      545
#> 2 110.48386 134.7576         1     96.3200     0.8919957      566
#> 3 101.39462 128.5609         1     76.0800    -3.1312514      394
#> 4 104.53616 127.7494         2    150.8900    -1.3654053      424
#> 5 118.44609 146.9993         3    158.6233     5.1262952      564
#> 6  95.08232 118.4228         3    102.5233     2.4778564      570
#>   estimated_total
#> 1        66122.45
#> 2        69403.34
#> 3        45301.23
#> 4        49244.54
#> 5        74855.60
#> 6        60848.96
#> ... and 6 more rows.
#> 
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_unit(formula = CornHec ~ CornPix + SoyBeansPix, unit_data = df_sample, Xpop
#> = df_pop, domain_var = "CountyIndex", popsize_var = "PopnSegments")
#> 
#> ✔ Convergence: Yes (INLA mode status 0)
#> Model: HB-UNIT-GAUSSIAN (Non-spatial)
#> Engine: inla
#> Method: INLA (laplace)
#> Random effect variance (sigma2_u): 90.54501 
#> Spatial mixing fraction (phi): NA 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  1.7975e+01  3.1524e+01  5.7020e-01  5.6854e-01 -4.4131e+01
#> CornPix      3.6432e-01  6.5991e-02  5.5208e+00  3.3749e-08  2.3409e-01
#> SoyBeansPix -2.8600e-02  6.8675e-02 -4.1645e-01  6.7708e-01 -1.6393e-01
#>             ci_upper
#> (Intercept)  80.1491
#> CornPix       0.4942
#> SoyBeansPix   0.1068
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred    vardir       sd      mse      rse
#> 1      1 165.7600 121.3256    121.2465 294.62568 6.358679 40.43280 5.241004
#> 2      2  96.3200 122.6207    122.6693 294.62568 6.192291 38.34447 5.049954
#> 3      3  76.0800 114.9777    115.0714 294.62568 6.930166 48.02720 6.027398
#> 4      4 150.8900 116.1428    115.9772 147.31284 5.921748 35.06710 5.098679
#> 5      5 158.6233 132.7227    132.5842  98.20856 7.283982 53.05640 5.488121
#> 6      6 102.5233 106.7526    106.7750  98.20856 5.954208 35.45259 5.577578
#>    ci_lower ci_upper samp_size sample_mean random_effect pop_size
#> 1 108.86258 133.7886         1    165.7600     1.1374532      545
#> 2 110.48386 134.7576         1     96.3200     0.8919957      566
#> 3 101.39462 128.5609         1     76.0800    -3.1312514      394
#> 4 104.53616 127.7494         2    150.8900    -1.3654053      424
#> 5 118.44609 146.9993         3    158.6233     5.1262952      564
#> 6  95.08232 118.4228         3    102.5233     2.4778564      570
#>   estimated_total
#> 1        66122.45
#> 2        69403.34
#> 3        45301.23
#> 4        49244.54
#> 5        74855.60
#> 6        60848.96
#> ... and 6 more rows.
#> 
# }
```
