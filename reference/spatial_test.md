# Pre-Modeling Spatial Diagnostic and Autocorrelation Tests

Evaluates whether spatial autocorrelation is present in direct domain
estimates and regression residuals, conducts Lagrange Multiplier (LM)
tests (Anselin, 1988), and performs comparative model fit tests between
the standard Fay-Herriot
([`eblup_fh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_fh.md))
and Spatial Fay-Herriot
([`eblup_sfh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md))
models. Provides an explicit, automated recommendation on whether a
spatial model is statistically necessary.

## Usage

``` r
spatial_test(
  formula,
  vardir = NULL,
  data,
  W = NULL,
  alpha = 0.05,
  alternative = c("greater", "two.sided", "less")
)
```

## Arguments

- formula:

  Regression formula describing the model (e.g. `y ~ x1 + x2`).

- vardir:

  Optional vector, column name, or formula referencing sampling
  variances of the direct estimator. Required for full model comparison
  (Standard FH vs Spatial FH).

- data:

  A `data.frame` or `sf` polygon object containing the variables in
  `formula`.

- W:

  Proximity or spatial adjacency matrix. If `NULL` and `data` is an `sf`
  polygon object, it is automatically constructed via
  [`create_weights`](https://ridsonap.is-a.dev/fastsae/reference/create_weights.md).

- alpha:

  Significance level for hypothesis tests (default is 0.05).

- alternative:

  Alternative hypothesis for Moran's I (`"greater"`, `"two.sided"`,
  `"less"`).

## Value

An object of class `c("fastsae_spatial_test", "list")` containing:

- response_test:

  Global Moran's I test results on direct response variable \\y\\.

- residual_test:

  Moran's I, Anselin LM-Error, and Anselin LM-Lag test results on
  regression residuals.

- model_comparison:

  ML log-likelihoods, AIC, BIC, and the likelihood ratio test (LRT) of
  \\\rho = 0\\ comparing standard vs spatial FH, both fitted by ML (if
  `vardir` provided).

- spatial_recommended:

  Logical recommendation, from a single test at level `alpha` (see
  Details).

- decision_test:

  The test behind `spatial_recommended`: `"lrt"` or `"lm_error"`.

- reason:

  Character string explaining the recommendation.

## Details

Moran's \\I\\ uses the randomisation variance of Cliff and Ord (1981);
the LM-error and LM-lag statistics follow Anselin (1988). These tests
assume homoscedastic OLS errors, whereas direct estimates have
domain-specific sampling variances, so they are descriptive. The
recommendation therefore rests on one test, so that its false positive
rate is `alpha`: with `vardir`, the likelihood ratio test of \\\rho =
0\\ between
[`eblup_fh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_fh.md)
and
[`eblup_sfh`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md)
fitted by ML (which accounts for the sampling variances); without it,
the LM-error test.

## References

1.  Anselin, L. (1988). *Spatial Econometrics: Methods and Models*.
    Kluwer Academic Publishers.

2.  Cliff, A. D., & Ord, J. K. (1981). *Spatial Processes: Models &
    Applications*. Pion London.

3.  Petrucci, A., & Salvati, N. (2006). Small area estimation for
    spatial correlation in watershed erosion assessment. *Journal of
    Agricultural, Biological, and Environmental Statistics*, 11(2),
    169-182.

4.  Pratesi, M., & Salvati, N. (2008). Small area estimation: the EBLUP
    estimator based on spatially correlated random area effects.
    *Statistical Methods & Applications*, 17(1), 113-141.

## Examples

``` r
library(fastsae)
data(mys)
data(mys_proxmat)

# Run spatial diagnostics on simulated survey data
st <- spatial_test(
  formula = y ~ x1 + x2,
  vardir = "vardir",
  data = mys,
  W = mys_proxmat
)
print(st)
#> 
#> ── Pre-Modeling Spatial Diagnostic Tests (fastsae) ─────────────────────────────
#> 
#> ── 1. Direct Response Variable (ESDA) ──
#> 
#> • Global Moran's I: -0.2773 (Expected: -0.0323, z-stat: -2.049)
#> • p-value: "0.9798" [ns]
#> ℹ Direct estimates show no significant spatial clustering.
#> 
#> ── 2. Regression Residuals (Controlling for Auxiliary Variables) ──
#> 
#> • Residual Moran's I: -0.2577 (z-stat: -1.875, p-value: "0.9696") [ns]
#> • Anselin LM-Error Test: 3.967 (df = 1, p-value: "0.0464") [*]
#> • Anselin LM-Lag Test: 2.055 (df = 1, p-value: "0.1517")
#> 
#> ── 3. Model Comparison (Standard FH vs Spatial FH) ──
#> 
#> • Standard FH: sigma2_u = 2.2353, ML LogLik = -65.22, AIC = 138.44
#> • Spatial FH (SAR): sigma2_u = 0.8462, rho = -1, ML LogLik = -63.12, AIC =
#>   136.23
#> • Likelihood Ratio: Chisq = 4.2 (df = 1, p-value = "0.04042") [*]
#> • Information Delta: Delta AIC = -2.2, Delta BIC = -0.73
#> 
#> ── 4. Final Recommendation ──
#> 
#> Decision test: likelihood ratio test (FH vs Spatial FH) at alpha = 0.05.
#> ✔ SPATIAL MODEL RECOMMENDED: Significant spatial spillover detected.
#> ✔ The likelihood ratio test of rho = 0 rejects the standard Fay-Herriot model
#>   in favour of the Spatial Fay-Herriot model.
#> ℹ Recommended functions: `eblup_sfh()` (Frequentist) or `hb_area(..., spatial =
#>   'bym2')` (Bayesian).
```
