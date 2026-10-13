# Set up the Python environment for the NumPyro / JAX engine

Creates (or updates) a Python environment with `jax` and `numpyro` for
`hb_area(engine = "jax")` and activates it for reticulate.

## Usage

``` r
setup_numpyro_env(envname = "r-numpyro-gpu", method = "auto", device = "auto")
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
  `jax==0.4.34`, `jaxlib==0.4.34`, `jax-metal==0.1.1` for Metal (NUTS
  still samples on CPU there).

## Value

Invisible logical `TRUE` upon successful setup.

## See also

[`check_numpyro_available`](https://ridsonap.github.io/fastsae/reference/check_numpyro_available.md),
[`hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md)

## Examples

``` r
if (FALSE) { # \dontrun{
setup_numpyro_env(device = "auto")
} # }
```
