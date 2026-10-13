## Test environments
* local macOS 26.6 (aarch64), R 4.6.1
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)
* GitHub Actions, CRAN-like: Ubuntu (release), `--as-cran`, without the non-CRAN package 'INLA'
* win-builder: R-release (4.6.1)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Suggests or Enhances not in mainstream repositories: INLA
  Availability using Additional_repositories specification: INLA yes
  https://inla.r-inla-download.org/R/stable

  'INLA' is only in Suggests and used conditionally (`requireNamespace()`): the
  examples that fit INLA models run only inside `if (requireNamespace("INLA"))`, the
  tests call `skip_if_not_installed("INLA")`, and the vignette chunks that fit
  INLA models are evaluated only when 'INLA' is installed.

  Possibly misspelled words in DESCRIPTION:
    BYM, Besag, ICAR, Phan, RW, Rao, Torabi, spatio

  These are author names (Besag, Phan, Rao, Torabi) and standard terms of the
  field: BYM (Besag-York-Mollie model), ICAR (intrinsic conditional
  autoregressive), RW (random walk, as in RW1/RW2) and "spatio-temporal".

## Optional Python backend
`hb_area(engine = jax_engine())` fits models by MCMC with the Python libraries
'JAX' and 'NumPyro' through 'reticulate' (Suggests). Nothing is installed or
downloaded unless the user calls `setup_numpyro_env()`; `check_numpyro_available()`
only inspects the existing Python configuration (it disables the managed 'uv'
environment of 'reticulate' during the check). The corresponding examples run only
inside `if (check_numpyro_available())` and the tests are skipped on CRAN.

## Submission summary
This is an update from version 0.1.0 to 1.0.0.

* Updated package title to "High-Performance Frequentist and Bayesian Small Area Estimation in R" to reflect the inclusion of both Frequentist and Bayesian methods.
* Added Hierarchical Bayes Small Area Estimation models via INLA (`hb_area`, `hb_unit`, `hb_twofold`), with engine objects `inla_engine()` and `jax_engine()` (optional MCMC backend in Python).
* Added Two-Fold subarea EBLUP models (`eblup_twofold`).
* Added model benchmarking (`benchmark_sae`), diagnostics (`diagnose`), comparison (`compare_sae`), mapping (`map_sae`), and export (`export_sae`).
* Statistical corrections to the existing EBLUP models (see NEWS.md).
