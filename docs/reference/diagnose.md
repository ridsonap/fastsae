# Diagnostic Evaluation of Small Area Estimation Models

Conducts a comprehensive diagnostic evaluation of small area estimation
models (HB, EBP, or EBLUP). Evaluates estimation precision, efficiency
gains over direct estimators, external calibration and bias diagnostics
(Brown et al., 2001), goodness-of-fit tests, and residual spatial
autocorrelation (Moran's I).

## Usage

``` r
diagnose(
  object,
  W = NULL,
  truth = NULL,
  rse_threshold = 25,
  alpha_level = 0.05
)
```

## Arguments

- object:

  A fitted model object of class `"fastsae"` (e.g., from
  [`hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md),
  [`eblup_fh`](https://ridsonap.github.io/fastsae/reference/eblup_fh.md),
  [`eblup_sfh`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md),
  or
  [`eblup_stfh`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md)).

- W:

  Optional spatial proximity or adjacency matrix of dimension \\D \times
  D\\. If `NULL` and the fitted model contains a spatial matrix (e.g.
  `eblup_sfh` or spatial `hb_area`), it is automatically extracted from
  `object`.

- truth:

  Optional numeric vector containing true parameter values \\\theta_d\\
  (useful for simulation studies).

- rse_threshold:

  Numeric. Threshold for reliable Relative Standard Error (default is
  25%).

- alpha_level:

  Numeric. Significance level for hypothesis tests (default is 0.05).

## Value

An object of class `c("fastsae_diagnose", "list")` containing:

- model_info:

  Basic model characteristics (model name, sample size, unsampled
  count).

- precision:

  Summary of direct vs SAE RSE, proportion of areas with RSE \<
  threshold, and MSE reduction ratios.

- brown_test:

  Results of the Brown et al. (2001) calibration test (\\H_0: \alpha =
  0, \beta = 1\\) and Chi-Square goodness-of-fit statistic \\W\\.

- spatial_test:

  Moran's I test for spatial autocorrelation in model residuals (if `W`
  is available).

- simulation_metrics:

  Relative bias (RB), Relative RMSE, and coverage rates (if `truth` is
  supplied).

- df_diag:

  Data frame combining direct estimates, SAE predictions, MSE, RSE, and
  residuals for diagnostics.

- status:

  Overall assessment summary string.

## References

1.  Brown, G., Chambers, R., Heady, P., & Heasman, D. (2001). Evaluation
    of small area estimation methods: An application to the British
    Labour Force Survey. *ONS Internal Report*.

2.  Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd
    ed.). John Wiley & Sons.

3.  Cliff, A. D., & Ord, J. K. (1981). *Spatial Processes: Models &
    Applications*. Pion London.

## Examples

``` r
library(fastsae)
data(mys)
data(mys_proxmat)

# 1. Fit Fay-Herriot model
fit_fh <- eblup_fh(y ~ x1 + x2, data = mys, vardir = ~vardir, domain = ~area)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_fh(formula = y ~ x1 + x2, vardir = ~vardir, domain = ~area, data = mys)
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
diag_fh <- diagnose(fit_fh)
print(diag_fh)
#> ── fastsae Small Area Estimation Diagnostic Report ─────────────────────────────
#> Model: "FH" | Domains: 42 (Sampled: 32, Unsampled: 10)
#> 
#> ── 1. Precision & Efficiency Gain ──
#> 
#> ! Domains with RSE < 25%: 52.4% (Caution: low precision)
#> ✔ Average RSE reduction: Direct "29.49%" -> SAE "25.59%" (Gain: 3.9%)
#> ✔ Variance reduction in 100% of areas (MSE ratio median: 1.33, max: 4.44)
#> 
#> ── 2. Brown et al. (2001) Calibration Tests ──
#> 
#> ! Bias Regression Test (H0: alpha = 0, beta = 1): F = 3.503, p-value = 0.0429 [Potential Systematic Bias]
#> Estimated parameters: alpha = -0.723, beta = 1.1587
#> ! Goodness-of-Fit Statistic W (Chi-Square): W = 4.57 (df = 32), p-value = 1 [Deviation from Survey Variance]
#> ────────────────────────────────────────────────────────────────────────────────
#> ! Final Assessment: CAUTION: Potential systematic bias detected between direct and model predictions.
#> ────────────────────────────────────────────────────────────────────────────────

# 2. Fit Spatial HB model (requires INLA)
if (requireNamespace("INLA", quietly = TRUE)) {
  fit_hb <- hb_area(y ~ x1 + x2, data = mys, vardir = "vardir")
  diag_hb <- diagnose(fit_hb)
  print(diag_hb)
}
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> hb_area(formula = y ~ x1 + x2, data = mys, vardir = "vardir")
#> 
#> ✔ Convergence: Yes (in - iterations)
#> Model: HB-GAUSSIAN (Non-spatial)
#> Random effect variance (sigma2_u): 1.68024 
#> 
#> Fixed Effects Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  2.9985e+00  6.9160e-01  4.3356e+00  1.4534e-05  1.6566e+00
#> x1          -3.6458e-03  8.0956e-03 -4.5035e-01  6.5246e-01 -1.9768e-02
#> x2           8.6096e-02  2.7872e-02  3.0890e+00  2.0085e-03  3.1167e-02
#>             ci_upper
#> (Intercept)   4.3809
#> x1            0.0122
#> x2            0.1410
#> 
#> HB Estimates (First 6 domains):
#>   domain        y       hb linear_pred        sd       mse      rse ci_lower
#> 1      1 8.359527 7.334884    7.334884 0.7431311 0.5522438 10.13146 5.897808
#> 2      2 7.599650 6.513222    6.513223 0.8022367 0.6435837 12.31705 4.966331
#> 3      3 5.514137 5.149918    5.149918 0.7842271 0.6150121 15.22795 3.622870
#> 4      4 3.869326 4.363556    4.363556 0.7084527 0.5019052 16.23567 2.962393
#> 5      5 6.305063 6.362973    6.362973 0.8810989 0.7763352 13.84728 4.631022
#> 6      6 3.926807 4.147447    4.147447 0.5686377 0.3233488 13.71055 3.028114
#>   ci_upper random_effect    vardir
#> 1 8.810730    2.71678122 0.6618838
#> 2 8.111717    2.28444112 0.8374691
#> 3 6.701018    0.72130815 0.8822257
#> 4 5.742045   -1.32531384 0.6581716
#> 5 8.092076   -0.08178645 1.2788021
#> 6 5.258777   -0.99655204 0.3878004
#> ... and 36 more rows.
#> 
#> ── fastsae Small Area Estimation Diagnostic Report ─────────────────────────────
#> Model: "HB-GAUSSIAN (Non-spatial)" | Domains: 42 (Sampled: 32, Unsampled: 10)
#> 
#> ── 1. Precision & Efficiency Gain ──
#> 
#> ! Domains with RSE < 25%: 61.9% (Caution: low precision)
#> ✔ Average RSE reduction: Direct "29.49%" -> SAE "22.95%" (Gain: 6.54%)
#> ✔ Variance reduction in 100% of areas (MSE ratio median: 1.55, max: 5.96)
#> 
#> ── 2. Brown et al. (2001) Calibration Tests ──
#> 
#> ✔ Bias Regression Test (H0: alpha = 0, beta = 1): F = 3, p-value = 0.0649 [Statistically Unbiased]
#> Estimated parameters: alpha = -0.805, beta = 1.1794
#> ! Goodness-of-Fit Statistic W (Chi-Square): W = 7.34 (df = 32), p-value = 1 [Deviation from Survey Variance]
#> 
#> ── 5. Bayesian Information Criteria & Predictive Diagnostics ──
#> 
#> ℹ WAIC: 117.61 (p_eff: 14.05) | DIC: 118.62 (p_eff: 20.03)
#> ✔ PIT Calibration Test vs Uniform(0,1): D = 0.084, p-value = 0.9651 [Well-calibrated predictive distribution]
#> ✔ Leave-One-Out CPO: No numerical approximation issues (min CPO = 0.0084)
#> ────────────────────────────────────────────────────────────────────────────────
#> ! Final Assessment: CAUTION: Model goodness-of-fit indicates notable deviation from survey variance.
#> ────────────────────────────────────────────────────────────────────────────────
```
