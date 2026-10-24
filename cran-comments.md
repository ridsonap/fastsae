## Test environments
* local macOS 26.6.2 (aarch64), R 4.6.1, `--as-cran`
* win-builder: Windows, R-devel and R-release
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)
* GitHub Actions, CRAN-like: Ubuntu (release), `--as-cran`, without the non-CRAN package 'INLA'

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Possibly misspelled words in DESCRIPTION:
    BYM, Besag, ICAR, Phan, RW, Rao, Torabi, spatio

  These are author names (Besag, Phan, Rao, Torabi) and standard terms of the
  field: BYM (Besag-York-Mollie model), ICAR (intrinsic conditional
  autoregressive), RW (random walk, as in RW1/RW2) and "spatio-temporal".

The same check also reports, as INFO, that 'INLA' (Suggests) is not in a
mainstream repository; it is available from the `Additional_repositories` entry
(https://inla.r-inla-download.org/R/stable). 'INLA' is used conditionally
(`requireNamespace()`): the examples that fit INLA models run only inside
`if (requireNamespace("INLA"))`, the tests call `skip_if_not_installed("INLA")`,
and the vignette chunks that fit INLA models are evaluated only when 'INLA' is
installed.

## Optional Python backend
`hb_area(engine = jax_engine())` fits models by MCMC with the Python libraries
'JAX' and 'NumPyro' through 'reticulate' (Suggests). Nothing is installed or
downloaded unless the user calls `setup_numpyro_env()`;
`setup_numpyro_env(check_only = TRUE)` only inspects the existing Python configuration
(it disables the managed 'uv' environment of 'reticulate' during the check). It does not
search the system for Python installations when no Python is configured for
'reticulate': a previous win-builder run flagged its example at 23 s elapsed on Windows
(0.5 s CPU), spent in that search; it now returns `FALSE` in about 1 s. The corresponding
examples run only inside `if (setup_numpyro_env(check_only = TRUE))` and the tests are
skipped on CRAN.

## Submission summary
This is an update from version 0.1.0 to 1.0.0.

* Updated package title to "High-Performance Frequentist and Bayesian Small Area Estimation" to reflect the inclusion of both Frequentist and Bayesian methods.
* Added Hierarchical Bayes Small Area Estimation models via INLA (`hb_area`, `hb_unit`, `hb_twofold`), with engine objects `inla_engine()` and `jax_engine()` (optional MCMC backend in Python).
* Added Two-Fold subarea EBLUP models (`eblup_twofold`).
* Added model benchmarking (`benchmark_sae`), diagnostics (`diagnose`), comparison (`compare_sae`), mapping (`map_sae`), and export (`export_sae`).
* `eblup_stfh()` now predicts unsampled domains and missing periods (NA responses)
  instead of requiring a complete panel.
* Statistical corrections to the existing EBLUP models (see NEWS.md).

## Changes for users of 0.1.0
All functions of 0.1.0 are still exported; no reverse dependencies on CRAN.

* `eblup_sfh(mse_method = )` and `compute_mse = ` of `eblup_stfh()` and `eblup_bhf()` are
  deprecated in favour of `mse_type` and still work with a warning ('lifecycle').
* `eblup_bhf()` now computes the bootstrap MSE by default (`mse_type = "bootstrap"`);
  `mse_type = "none"` gives the 0.1.0 output.
* New arguments were inserted before existing ones (`transform` and `mse_type` in
  `eblup_fh()`, `mse_type` in `eblup_sfh()`), so calls that pass those later arguments by
  position change meaning; calls by name are unaffected.
