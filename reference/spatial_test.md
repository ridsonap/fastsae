# Pre-Modeling Spatial Diagnostic and Autocorrelation Tests

Evaluates whether spatial autocorrelation is present in direct domain
estimates and regression residuals, conducts Lagrange Multiplier (LM)
tests (Anselin, 1988), and performs comparative model fit tests between
the standard Fay-Herriot
([`eblup_fh`](https://ridsonap.github.io/fastsae/reference/eblup_fh.md))
and Spatial Fay-Herriot
([`eblup_sfh`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md))
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
  [`create_weights`](https://ridsonap.github.io/fastsae/reference/create_weights.md).

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

  Log-Likelihood, AIC, BIC, and Likelihood Ratio Test (LRT) comparing
  standard vs spatial FH (if `vardir` provided).

- recommendation:

  Character string with diagnosis and logical flag
  `spatial_recommended`.

## References

1.  Anselin, L. (1988). *Spatial Econometrics: Methods and Models*.
    Kluwer Academic Publishers.

2.  Cliff, A. D., & Ord, J. K. (1981). *Spatial Processes: Models &
    Applications*. Pion London.

3.  Petrucci, A., & Salvati, N. (2006). Small area estimation for
    spatial data in R. *Statistical Methods & Applications*, 14(3),
    391-407.

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
#> • Global Moran's I: -0.0404 (Expected: -0.0323, z-stat: -1.002)
#> • p-value: "0.8417" [ns]
#> ℹ Direct estimates show no significant spatial clustering.
#> 
#> ── 2. Regression Residuals (Controlling for Auxiliary Variables) ──
#> 
#> • Residual Moran's I: -0.0418 (z-stat: -1.165, p-value: "0.878") [ns]
#> • Anselin LM-Error Test: 0.835 (df = 1, p-value: "0.3608") [ns]
#> • Anselin LM-Lag Test: 0.826 (df = 1, p-value: "0.3633")
#> 
#> ── 3. Model Comparison (Standard FH vs Spatial FH) ──
#> 
#> • Standard FH: sigma2_u = 2.5692, REML LogLik = -65.29, AIC = 138.59
#> • Spatial FH (SAR): sigma2_u = 1.468, rho = -1, REML LogLik = -65.43, AIC =
#>   140.86
#> • Likelihood Ratio: Chisq = 0 (df = 1, p-value = "1") [ns]
#> • Information Delta: Delta AIC = 2.27, Delta BIC = 3.74
#> 
#> ── 4. Final Recommendation ──
#> 
#> ℹ STANDARD MODEL SUFFICIENT: No significant spatial spillover.
#> ✔ Neither the response variable nor regression residuals exhibit significant
#>   spatial autocorrelation. A standard model is sufficient.
#> ℹ Recommended function: `eblup_fh()` (avoids unnecessary g3 MSE inflation).
```
