# Export Small Area Estimation Results to Excel or CSV

Exports comprehensive estimation tables and model summaries to
multi-sheet Excel workbooks (`.xlsx`) or comma-separated files (`.csv`).
Formatted to meet official statistics dissemination standards.

The generated Excel workbook contains up to three dedicated sheets:

- **Estimates**: Domain IDs, Direct survey estimates, SAE model
  predictions, Standard Errors (SE), MSE, RSE (%), 95% confidence /
  credible intervals, and official Reliability Flags.

- **Model_Summary**: Model formula, method, convergence status,
  regression coefficients (\\\hat{\beta}\\, SE, p-values), variance
  components (\\\sigma_u^2\\, spatial \\\rho\\, temporal \\\rho_t\\,
  \\\phi\\), and goodness-of-fit metrics.

- **Benchmarked** (Optional): Pre- vs post-calibration values,
  adjustments, and percentage shifts when benchmarking calibration is
  applied.

## Usage

``` r
export_sae(
  object,
  file,
  benchmark = NULL,
  thresholds = c(20, 30),
  overwrite = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted `fastsae` object or a `fastsae_benchmark` object.

- file:

  Character string specifying the target file path (must end with
  `.xlsx` or `.csv`).

- benchmark:

  Optional `fastsae_benchmark` object to include calibration results
  alongside the model.

- thresholds:

  Numeric vector of length 2 defining the RSE (%) thresholds for
  reliability flags. Default is `c(20, 30)`.

- overwrite:

  Logical indicating whether to overwrite an existing file. Default is
  `TRUE`.

- ...:

  Additional arguments.

## Value

Invisibly returns a named list of data frames corresponding to the
exported sheets.

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
tmp_file <- tempfile(fileext = ".xlsx")
export_sae(fit_fh, file = tmp_file)
#> ✔ Successfully exported SAE results to Excel: /var/folders/j2/wt412qcx0g704rgp5p9y6l940000gn/T//RtmpY2XkmJ/file10275276929c0.xlsx
unlink(tmp_file)
```
