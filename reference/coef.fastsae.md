# Extract coefficients from a fastsae object

Extract coefficients from a fastsae object

## Usage

``` r
# S3 method for class 'fastsae'
coef(object, ...)
```

## Arguments

- object:

  An object of class `fastsae`.

- ...:

  Additional arguments.

## Value

Named vector of estimated coefficients.

## Examples

``` r
m <- eblup_fh(y ~ x1 + x2 + x3, vardir = "vardir", data = mys, print_result = FALSE)
coef(m)
#>  (Intercept)           x1           x2           x3 
#>  3.107750953 -0.001932259  0.055518370  0.033534410 
```
