# Check whether the NumPyro / JAX engine is available

Checks that `hb_area(engine = "jax")` can run: the package reticulate
and a Python with the modules `jax` and `numpyro`. It looks for the
`r-numpyro-gpu` virtualenv (see
[`setup_numpyro_env`](https://ridsonap.github.io/fastsae/reference/setup_numpyro_env.md))
or any Python already configured for reticulate. The check never
installs or downloads anything: the managed ('uv') Python of reticulate
is not provisioned unless the `RETICULATE_USE_MANAGED_VENV` environment
variable is set.

## Usage

``` r
check_numpyro_available(device = "auto")
```

## Arguments

- device:

  Character: `"auto"`, `"metal"`, `"cuda"`, or `"cpu"`.

## Value

Logical indicating whether the environment is properly configured.

## See also

[`setup_numpyro_env`](https://ridsonap.github.io/fastsae/reference/setup_numpyro_env.md),
[`hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md)

## Examples

``` r
check_numpyro_available()
#> [1] FALSE
```
