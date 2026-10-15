# Check whether the NumPyro / JAX engine is available

Checks that `hb_area(engine = "jax")` can run: the package reticulate
and a Python with the modules `jax` and `numpyro`. It uses the Python
already initialised in the session, otherwise the one named by the
environment variables `RETICULATE_PYTHON` or `RETICULATE_PYTHON_ENV`,
otherwise the `r-numpyro-gpu` environment created by
[`setup_numpyro_env`](https://ridsonap.github.io/fastsae/reference/setup_numpyro_env.md)
(virtualenv, or conda when `conda` is on the `PATH`). When none of these
exists it returns `FALSE` at once, without searching the system for
other Python installations (which can take tens of seconds on some
platforms). To use another Python, set `RETICULATE_PYTHON`, or select it
with
[`reticulate::use_python()`](https://rstudio.github.io/reticulate/reference/use_python.html)
and initialise it (e.g.
[`reticulate::py_config()`](https://rstudio.github.io/reticulate/reference/py_config.html))
before the check. The check never installs or downloads anything: the
managed ('uv') Python of reticulate is not provisioned unless the
`RETICULATE_USE_MANAGED_VENV` environment variable is set.

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
