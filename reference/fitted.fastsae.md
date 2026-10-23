# Extract fitted values (EBLUP or HB) from a fastsae object

Extract fitted values (EBLUP or HB) from a fastsae object

## Usage

``` r
# S3 method for class 'fastsae'
fitted(object, ...)
```

## Arguments

- object:

  An object of class `fastsae`.

- ...:

  Additional arguments.

## Value

Vector of fitted estimates.

## Examples

``` r
m <- eblup_fh(y ~ x1 + x2 + x3, vardir = "vardir", data = mys, print_result = FALSE)
head(fitted(m))
#> [1] 7.612738 6.782316 5.187060 4.201545 6.323679 4.048590
```
