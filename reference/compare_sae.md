# Compare Two Small Area Estimation Models

Provides a comprehensive evaluation and concordance comparison between
two fitted Small Area Estimation (SAE) models (e.g. comparing Classical
Fay-Herriot vs Spatial Fay-Herriot, Area-level vs Unit-level BHF, or
EBLUP vs EBP).

Computes empirical agreement and efficiency metrics including:

- Pearson linear correlation (\\r\\) and Spearman rank correlation
  (\\\rho\\).

- Mean Absolute Difference (MAE) and Root Mean Squared Difference
  (RMSD).

- Relative efficiency metrics: Mean MSE ratio (\\\text{MSE}\_1 /
  \text{MSE}\_2\\), and the proportion of domains where Model 2 achieves
  greater precision.

- Average Relative Standard Error (RSE %) across models and precision
  gains.

## Usage

``` r
compare_sae(model1, model2 = NULL, names = NULL, thresholds = c(20, 30), ...)

# Default S3 method
compare_sae(model1, model2 = NULL, names = NULL, thresholds = c(20, 30), ...)

# S3 method for class 'fastsae_comparison'
print(x, ...)

# S3 method for class 'fastsae_comparison'
summary(object, ...)

# S3 method for class 'fastsae_comparison'
plot(x, y = NULL, ...)

# S3 method for class 'fastsae_comparison'
autoplot(
  object,
  type = c("scatter", "comparison", "difference", "mse", "rse"),
  title = NULL,
  ...
)
```

## Arguments

- model1:

  A fitted `fastsae` model object, or a list containing two `fastsae`
  models.

- model2:

  A second fitted `fastsae` model object. Ignored if `model1` is a list.

- names:

  Optional character vector of length 2 specifying descriptive labels
  for the models. Default is inferred from model types or call
  arguments.

- thresholds:

  Numeric vector of length 2 defining the RSE (%) thresholds for
  reliability. Default is `c(20, 30)`.

- ...:

  Additional arguments.

- x:

  An object of class `fastsae_comparison` (for `print` and `plot`
  methods).

- object:

  An object of class `fastsae_comparison` (for `summary` and `autoplot`
  methods).

- y:

  Ignored argument for compatibility with the generic `plot` method.

- type:

  Character string indicating comparison plot type: `"scatter"`,
  `"comparison"`, `"difference"`, `"mse"`, or `"rse"`. Default is
  `"scatter"`.

- title:

  Optional character string specifying a custom plot title.

## Value

An S3 object of class `fastsae_comparison` containing:

- `metrics`: A data frame of overall concordance and efficiency
  statistics.

- `data`: A domain-level data frame with aligned estimates, MSE, and RSE
  values.

- `names`: Vector of the two model names.

- `thresholds`: The RSE thresholds used.

## See also

`autoplot.fastsae_comparison`,
[`map_sae`](https://ridsonap.github.io/fastsae/reference/map_sae.md)

## Examples

``` r
library(fastsae)
data(mys)

# Fit two models
fit_fh <- eblup_fh(y ~ x1 + x2, vardir = "vardir", data = mys)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_fh(formula = y ~ x1 + x2, vardir = "vardir", data = mys)
#> 
#> ✔ Convergence: Yes (in 7 iterations)
#> Model: Fay-Herriot (Area-level)
#> Method: eblup
#> Random effect variance (sigma2_u): 2.569182 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.0689476  0.7631917  4.0212018 0.0001
#> x1          -0.0041073  0.0090768 -0.4525074 0.6509
#> x2           0.0861173  0.0309576  2.7817847 0.0054
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    eblup    vardir random_effect       mse       rse
#> 1      1 8.359527 7.594807 0.6618838    2.96835204 0.5602338  9.855256
#> 2      2 7.599650 6.777056 0.8374691    2.52354833 0.6766501 12.137828
#> 3      3 5.514137 5.247512 0.8822257    0.77645304 0.7048881 15.999508
#> 4      4 3.869326 4.247072 0.6581716   -1.47453645 0.5556741 17.551750
#> 5      5 6.305063 6.353632 1.2788021   -0.09757891 0.9308379 15.185005
#> 6      6 3.926807 4.090117 0.3878004   -1.08192989 0.3497100 14.458337
#> ... and 36 more rows.
#> 
fit_sfh <- eblup_sfh(y ~ x1 + x2, vardir = "vardir", data = mys, W = mys_proxmat)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_sfh(formula = y ~ x1 + x2, vardir = "vardir", data = mys, W =
#> mys_proxmat)
#> 
#> ✔ Convergence: Yes (in 10 iterations)
#> Model: Spatial Fay-Herriot (Area-level SAR)
#> Method: REML
#> Random effect variance (sigma2_u): 1.467999 
#> Spatial autocorrelation (rho): -1 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  2.9786076  0.6242220  4.7717123 0.0000
#> x1          -0.0033919  0.0075283 -0.4505566 0.6523
#> x2           0.0857371  0.0261402  3.2799004 0.0010
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    vardir    eblup random_effect       mse       rse
#> 1      1 8.359527 0.6618838 7.253737    2.63803919 0.5206791  9.947713
#> 2      2 7.599650 0.8374691 6.415275    2.19528212 0.6143184 12.217475
#> 3      3 5.514137 0.8822257 5.121070    0.70679774 0.6380642 15.598092
#> 4      4 3.869326 0.6581716 4.399165   -1.28175155 0.5163702 16.334669
#> 5      5 6.305063 1.2788021 6.393072   -0.05118048 0.8064653 14.046990
#> 6      6 3.926807 0.3878004 4.159058   -0.98030040 0.3349916 13.916243
#> ... and 36 more rows.
#> 

# Compare models
comp <- compare_sae(fit_fh, fit_sfh, names = c("FH", "Spatial FH"))
print(comp)
#> 
#> ── Small Area Estimation Model Concordance & Comparison ────────────────────────
#> Model 1: FH
#> Model 2: Spatial FH
#> Matched Domains: 42
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#> Key Concordance & Efficiency Metrics:
#>                                                        Metric   Value
#>                                               Matched Domains      42
#>                                       Pearson Correlation (r)  0.9959
#>                                    Spearman Correlation (rho)  0.9948
#>                                Mean Absolute Difference (MAE)  0.1347
#>                           Root Mean Squared Difference (RMSD)  0.1668
#>                              Mean MSE Ratio (FH / Spatial FH)  1.2766
#>                            Median MSE Ratio (FH / Spatial FH)  1.1561
#>                        Domains where Spatial FH has lower MSE 41 / 42
#>  Percentage of domains where Spatial FH is more efficient (%)   97.6%
#>                                               Mean RSE % (FH)  25.59%
#>                                       Mean RSE % (Spatial FH)  22.16%
#>                                     Average RSE Reduction (%)   3.42%
#> 
summary(comp)
#> 
#> ── Small Area Estimation Model Concordance & Comparison ────────────────────────
#> Model 1: FH
#> Model 2: Spatial FH
#> Matched Domains: 42
#> ────────────────────────────────────────────────────────────────────────────────
#> 
#> Key Concordance & Efficiency Metrics:
#>                                                        Metric   Value
#>                                               Matched Domains      42
#>                                       Pearson Correlation (r)  0.9959
#>                                    Spearman Correlation (rho)  0.9948
#>                                Mean Absolute Difference (MAE)  0.1347
#>                           Root Mean Squared Difference (RMSD)  0.1668
#>                              Mean MSE Ratio (FH / Spatial FH)  1.2766
#>                            Median MSE Ratio (FH / Spatial FH)  1.1561
#>                        Domains where Spatial FH has lower MSE 41 / 42
#>  Percentage of domains where Spatial FH is more efficient (%)   97.6%
#>                                               Mean RSE % (FH)  25.59%
#>                                       Mean RSE % (Spatial FH)  22.16%
#>                                     Average RSE Reduction (%)   3.42%
#> 
#> Distribution of Domain Point Estimates:
#>          FH Spatial FH Difference
#> 0%   0.8555     0.9823    -0.4212
#> 25%  4.2455     4.2158    -0.1033
#> 50%  5.1879     4.9815    -0.0666
#> 75%  6.1655     6.2119     0.0704
#> 100% 9.0169     8.7804     0.2663
#> 
#> Distribution of Relative Standard Error (RSE %):
#>         FH Spatial FH RSE Diff (1 - 2)
#> 0%    9.86       9.95            -0.09
#> 25%  14.83      13.95             0.73
#> 50%  21.85      20.22             1.73
#> 75%  34.76      28.18             6.40
#> 100% 52.75      45.59             9.81
#> 

# Plot comparison
autoplot(comp, type = "scatter")

autoplot(comp, type = "difference")
```
