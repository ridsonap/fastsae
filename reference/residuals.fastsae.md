# Extract residuals from a fastsae object

Extract residuals from a fastsae object

## Usage

``` r
# S3 method for class 'fastsae'
residuals(object, ...)
```

## Arguments

- object:

  An object of class `fastsae`.

- ...:

  Additional arguments.

## Value

Vector of residuals for sampled areas.

## Examples

``` r
m <- eblup_fh(y ~ x1 + x2 + x3, vardir = "vardir", data = mys, print_result = FALSE)
head(residuals(m))
#> [1]  0.74678890  0.81733383  0.32707611 -0.33221827 -0.01861583 -0.12178363
```
