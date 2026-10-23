# Set up, or check, the Python environment for the NumPyro / JAX engine

Creates (or updates) a Python environment with `jax` and `numpyro` for
[`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
and activates it for reticulate. With `check_only = TRUE` it only checks
that the engine can run, and installs nothing.

## Usage

``` r
setup_numpyro_env(
  envname = "r-numpyro-gpu",
  method = "auto",
  device = "auto",
  check_only = FALSE
)
```

## Arguments

- envname:

  Name of the conda or virtualenv environment. Default
  `"r-numpyro-gpu"`.

- method:

  `"auto"` (virtualenv), `"virtualenv"`, or `"conda"`.

- device:

  `"auto"` (CUDA if `nvidia-smi` is found, otherwise CPU), `"metal"`
  (Apple Silicon), `"cuda"` (NVIDIA), or `"cpu"`. Versions are pinned to
  tested combinations: `jax==0.5.0` and `numpyro==0.19.0` for CPU/CUDA;
  `jax==0.4.34`, `jaxlib==0.4.34`, `jax-metal==0.1.1` for Metal. With
  `check_only = TRUE`, the device to check (`"auto"` uses the CPU on
  macOS).

- check_only:

  Logical; only check whether the engine can run (see Details), without
  installing anything. Default `FALSE`.

## Value

Invisible logical `TRUE` upon successful setup; with
`check_only = TRUE`, a logical telling whether the engine is available.

## Details

The check needs the package reticulate and a Python with the modules
`jax` and `numpyro`. It uses the Python already initialised in the
session, otherwise the one named by the environment variables
`RETICULATE_PYTHON` or `RETICULATE_PYTHON_ENV`, otherwise the
`r-numpyro-gpu` environment created by this function (virtualenv, or
conda when `conda` is on the `PATH`). When none of these exists it
returns `FALSE` at once, without searching the system for other Python
installations (which can take tens of seconds on some platforms). To use
another Python, set `RETICULATE_PYTHON`, or select it with
[`reticulate::use_python()`](https://rstudio.github.io/reticulate/reference/use_python.html)
and initialise it (e.g.
[`reticulate::py_config()`](https://rstudio.github.io/reticulate/reference/py_config.html))
before the check. The check never installs or downloads anything: the
managed ('uv') Python of reticulate is not provisioned unless the
`RETICULATE_USE_MANAGED_VENV` environment variable is set.

## See also

[`jax_engine`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md),
[`hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)

## Examples

``` r
setup_numpyro_env(check_only = TRUE)
#> [1] FALSE
if (FALSE) { # \dontrun{
setup_numpyro_env(device = "auto")
} # }
```
