# Print a fastsae object

Print a fastsae object

## Usage

``` r
# S3 method for class 'fastsae'
print(x, ...)
```

## Arguments

- x:

  An object of class `fastsae`.

- ...:

  Additional arguments passed to print methods.

## Value

The original object invisibly.

## Examples

``` r
m <- eblup_fh(y ~ x1 + x2 + x3, vardir = "vardir", data = mys, print_result = FALSE)
print(m)
#> 
#> ── Fast Small Area Estimation (fastsae) ────────────────────────────────────────
#> Call:
#> eblup_fh(formula = y ~ x1 + x2 + x3, vardir = "vardir", data = mys,
#> print_result = FALSE)
#> 
#> ✔ Convergence: Yes (in 6 iterations)
#> Model: Fay-Herriot (Area-level)
#> Method: REML
#> Random effect variance (sigma2_u): 2.608103 
#> 
#> Fixed Effects Coefficients:
#>                   beta  std.error     zvalue pvalue
#> (Intercept)  3.1077510  0.7697687  4.0372527 0.0001
#> x1          -0.0019323  0.0098886 -0.1954019 0.8451
#> x2           0.0555184  0.0614129  0.9040187 0.3660
#> x3           0.0335344  0.0580013  0.5781663 0.5632
#> 
#> EBLUP Estimates (First 6 domains):
#>   domain        y    eblup    vardir random_effect       mse       rse
#> 1      1 8.359527 7.612738 0.6618838    2.94266461 0.5616180  9.844182
#> 2      2 7.599650 6.782316 0.8374691    2.54539588 0.6784873 12.144869
#> 3      3 5.514137 5.187060 0.8822257    0.96692722 0.7190378 16.347620
#> 4      4 3.869326 4.201545 0.6581716   -1.31646412 0.5619314 17.841554
#> 5      5 6.305063 6.323679 1.2788021   -0.03796677 0.9370273 15.307573
#> 6      6 3.926807 4.048590 0.3878004   -0.81904048 0.3548310 14.713194
#> ... and 36 more rows.
#> 
```
