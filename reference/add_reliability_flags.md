# Add Official Statistical Reliability Flags

Categorizes Small Area Estimates into standard statistical publication
reliability tiers based on their Relative Standard Error (RSE / CV %),
following guidelines from official statistical agencies (such as
BPS-Statistics Indonesia, Eurostat, and US Census Bureau).

Default classification thresholds:

- **Reliable** (\\\text{RSE} \< 20\\\\): Suitable for unconditional
  official publication.

- **Use with Caution** (\\20\\ \le \text{RSE} \< 30\\\\): Usable with
  cautionary footnotes regarding sampling variability.

- **Unreliable** (\\\text{RSE} \ge 30\\\\): High sampling error;
  suppression or aggregation to higher geographic levels recommended.

## Usage

``` r
add_reliability_flags(object, rse_col = "rse", thresholds = c(20, 30))
```

## Arguments

- object:

  A `fastsae` model object or a data frame containing an RSE column.

- rse_col:

  Character string specifying the name of the RSE column. Default is
  `"rse"`.

- thresholds:

  Numeric vector of length 2 defining the RSE (%) cutoffs. Default is
  `c(20, 30)`.

## Value

The input object or data frame with an additional factor column
`reliability`.

## Examples

``` r
library(fastsae)
data(mys)

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
df_flagged <- add_reliability_flags(fit_fh)
table(df_flagged$reliability)
#> 
#>          Reliable (< 20%) Use with Caution (20-30%)        Unreliable (≥ 30%) 
#>                        19                         9                        14 
```
