# Summarize a fastsae object

Summarize a fastsae object

## Usage

``` r
# S3 method for class 'fastsae'
summary(object, ...)
```

## Arguments

- object:

  An object of class `fastsae`.

- ...:

  Additional arguments.

## Value

An object of class `summary.fastsae`.

## Examples

``` r
m <- eblup_fh(y ~ x1 + x2 + x3, vardir = "vardir", data = mys, print_result = FALSE)
s <- summary(m)
s$coefficients
#>                     beta   std.error stderr_beta     zvalue       pvalue
#> (Intercept)  3.107750953 0.769768746 0.769768746  4.0372527 0.0000540808
#> x1          -0.001932259 0.009888644 0.009888644 -0.1954019 0.8450783567
#> x2           0.055518370 0.061412854 0.061412854  0.9040187 0.3659854711
#> x3           0.033534410 0.058001322 0.058001322  0.5781663 0.5631518464
```
