# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What this is

`fastsae` is an R package (CRAN, GPL-3) for Small Area Estimation:
frequentist EBLUP models in compiled C++ (RcppArmadillo, RcppEigen,
OpenMP) and Bayesian hierarchical models via INLA or, optionally,
NumPyro/JAX MCMC through reticulate. Its main selling points are speed
and numerical agreement with reference packages: the EBLUPs agree with
`sae` (and `emdi`) to rounding error where both use the same method, and
the HB models are checked against exact posteriors, Stan and `tipsae`.
Changes to estimators should keep that agreement, and the tests compare
against those packages. The documented, deliberate differences from
`sae` (boundary `rho` and ML information in the SFH MSE; constrained
REML fallbacks) are listed in `NEWS.md` and the Rd files.

## Commands

Run from the repo root (see also `main.R`):

``` sh
Rscript -e 'Rcpp::compileAttributes()'          # after adding/renaming any // [[Rcpp::export]] in src/
Rscript -e 'devtools::document()'               # regenerate man/ from roxygen
Rscript -e 'devtools::test()'                   # full test suite
Rscript -e 'devtools::test(filter = "eblup_sfh")'   # single file: tests/testthat/test_eblup_sfh.R
Rscript -e 'devtools::check()'                  # R CMD check
R CMD INSTALL . --preclean                      # force a clean C++ rebuild
Rscript calc_cov.R                              # local coverage -> cov_res.rds (FLIBS path is Homebrew gcc 16 specific)
Rscript -e 'pkgdown::build_site()'              # site into docs/
```

- Bayesian tests and examples need **INLA**, which is not on CRAN:
  `install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)`.
  Without it, those tests are skipped (`skip_if_not_installed("INLA")`)
  and do not fail. The same goes for `sae`, `emdi`, `sf`, and `MASS`.
  The
  [`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
  tests need Python with `jax`/`numpyro`
  ([`setup_numpyro_env()`](https://ridsonap.is-a.dev/fastsae/reference/setup_numpyro_env.md))
  and are skipped on CRAN and without it.
- **Scripts in `benchmarks/`, `rjournal/` and the vignettes use the
  *installed* fastsae**, not `load_all()`. Run
  `R CMD INSTALL . --preclean` before benchmarks or renders, or the
  numbers come from an old build (this happened once: the installed
  build was days older than the source).
- **Indonesian documentation** (site only):
  `vignettes/articles/id/*.Rmd` translate the home page (`beranda.Rmd`,
  from `pkgdown/index.md`) and the 10 articles, plus `glosarium.Rmd`
  (term list). Technical terms stay in English (small area estimation,
  direct estimator, random effect, …); only the prose is Indonesian.
  Code, code comments and outputs are identical to the English article.
  **When an English article or `pkgdown/index.md` changes, update its
  Indonesian copy in the same PR.** `benchmarks/update_benchmark_docs.R`
  writes the generated blocks (data arrays, `BENCH-SUMMARY`,
  `BENCH-README`) into the Indonesian copies too, with Indonesian number
  format. The language switcher (`pkgdown/extra.js`) maps each English
  page to its `articles/id/` copy; the articles are listed in the
  `articles:` section of `_pkgdown.yml` (full vignette names,
  e.g. `articles/id/beranda`). They are excluded from the package build
  with the rest of `vignettes/articles/`.
- `pkgdown/index.md` is the site home page and is written separately
  from `README.Rmd`; `benchmarks/update_benchmark_docs.R` writes the
  benchmark tables into both, and into `vignettes/benchmarks.Rmd` and
  `vignettes/hb-area-inla.Rmd`.
- `README.md` is generated from `README.Rmd`. Edit the `.Rmd` and knit
  it.

## Package plumbing gotchas

- **`DESCRIPTION` has an explicit `Collate:` field.** A new `R/*.R` file
  must be added there or it won’t be loaded. Order matters: `utils.R`,
  `inla_utils.R`, `RcppExports.R`, and `methods.R` come before the model
  files.
- **`NAMESPACE` is maintained by hand.** It has no roxygen “Generated
  by” header, so `devtools::document()` will not rewrite it. Add new
  `export()` and `S3method()` lines manually.
- C++ functions are exported under dot-prefixed internal names via
  `// [[Rcpp::export(.name)]]`: `.eblup_core`, `.seblup_core`,
  `.seblup_pbmse`, `.seblup_npbmse`, `.eblup_stfh_core`, `.pbmse_stfh`,
  `.eblup_twofold_core`, `.bhf_fit_cpp`, `.pbmse_bhf_cpp`,
  `.eblup_bhf_cpp`. `R/RcppExports.R` and `src/RcppExports.cpp` are
  generated, so don’t edit them by hand.
- `src/Makevars` uses C++17 + OpenMP + LAPACK/BLAS. Stale `src/*.o` /
  `fastsae.so` files are in the tree locally, so use `--preclean` when
  header changes don’t seem to take effect.
- `.Rbuildignore` excludes `benchmarks/`, `rjournal/`, `paper/`,
  `playground/`, `docs/`, `data-raw/`, `vignettes/articles/`, `main.R`,
  `calc_cov.R`, and `inst/extdata`. Those directories are dev-only, and
  package code must not depend on them.
- Roxygen markdown is not enabled package-wide: a block that uses
  backticks or `*emphasis*` needs `@md`, or they appear literally in the
  Rd.

## Architecture

**Two engines, one result class.** Every model function returns a list
with S3 class `"fastsae"` (plus subclasses: `fastsae_unit`,
`fastsae_hb_area`, `fastsae_hb`, `fastsae_hb_unit`,
`fastsae_hb_twofold`). The shared methods in `R/methods.R` (`print`,
`summary`, `coef`, `fitted`, `residuals`, `plot`), as well as
`autoplot`, `compare_sae`, `benchmark_sae`, `diagnose`, `map_sae`, and
`export_sae`, all work on this common shape. They find estimates through
`.sae_table()` in `R/utils.R` (`df_hb %||% df_eblup`), so a new model
has to fill one of those data frames with the standard columns
(`domain, y, eblup|hb, vardir, mse, rse, ...`), plus `estcoef`,
`formula`, `call`, `data`, and `model`.

**Frequentist path (`eblup_*`)**: the R wrapper parses the
formula/`vardir`, builds `X` via `model.matrix`, calls a C++ core for
REML/ML Fisher scoring and MSE, then adds metadata and does
back-transformation in R. - `eblup_fh` → `src/eblup_fh.cpp` (supports
`transform = "log"`: conditional-mean back-transformation
`exp(mu + g1/2)` and delta-method MSE done in R, not the full Slud–Maiti
correction) - `eblup_sfh` (SAR spatial) → `src/eblup_sfh_core.cpp`, with
the sparse Cholesky machinery in `src/sar_engine.h` (no D x D inverse);
parametric and non-parametric bootstrap MSE live in
`eblup_sfh_pbmse.cpp` / `eblup_sfh_npbmse.cpp` (OpenMP, header
`eblup_sfh.h`) - `eblup_stfh` (spatio-temporal, Marhuenda et al.; REML
only) → `eblup_stfh_core.cpp` (+ analytical MSE g1+g2+2g3-g4,
block-structured, header `eblup_stfh_core.h`) and
`eblup_stfh_pbmse.cpp`; it sorts the panel into domain-major order
internally and predicts `NA` cells (unsampled domains, missing
periods) - `eblup_twofold` (sub-area level model of Torabi & Rao, 2014;
area means in `df_area`) → `eblup_twofold.cpp` - `eblup_bhf` (unit
level) fits REML/ML in C++ (`src/eblup_bhf_core.cpp`, likelihood
profiled in gamma = s2u/s2e) and predicts with `src/eblup_unit.cpp`; the
bootstrap MSE is driven by `.pbmse_unit()` in R with replicates in C++
(`.pbmse_bhf_cpp`). `lme4` is only in Suggests.

**Bayesian path (`hb_*`)**: `hb_area` (6 likelihood families; spatial
iid/`besag`/`bym`/`bym2`/`leroux`/`sar`/`slm`; temporal
`rw1`/`rw2`/`ar1`; Knorr-Held space-time interactions), `hb_unit` and
`hb_twofold` (Gaussian/binomial/Poisson) take an `engine` object from
`R/engine.R`:
[`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
(default; builds INLA formulas and calls
[`INLA::inla`](https://rdrr.io/pkg/INLA/man/inla.html)) or
[`jax_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
(NumPyro NUTS via reticulate; R side in `R/hb_*_jax.R` and
`R/python_env.R`, Python in `inst/python/hb_area_numpyro.py`). Both
engines share the default priors (half-normal sd(y) for Gaussian SDs,
half-normal(1) on the link scale otherwise; see
`vignettes/hb-priors.Rmd`). The helpers in `R/inla_utils.R` handle the
INLA install check, converting a weight matrix into an INLA graph, and
pulling posteriors out into `df_hb`. Gaussian models (and Leroux/SLM)
use INLA’s classic mode, others compact. The old `method = "laplace"`
lme4 route is removed (deprecated argument forwards to
[`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)).
`hb_area` temporarily sets INLA to single-threaded mode and restores the
setting on exit.

**Benchmarking comes in two forms. Don’t mix them up:** -
*Self-benchmarking* (`self_benchmark = TRUE` argument on the model
functions): `.setup_self_benchmark()` in `R/utils.R` augments the data
and formula before fitting, and `.compute_benchmark_summary()` reports
the calibration afterwards. - *Post-hoc benchmarking*
([`benchmark_sae()`](https://ridsonap.is-a.dev/fastsae/reference/benchmark_sae.md)
in `R/benchmark.R`; the `benchmark()` alias was removed) adjusts a
fitted object and returns class `fastsae_benchmark`.

**Supporting modules**: `create_weights` / `spatial_test` (spatial
weights from `sf` polygons, Moran/LM tests), `gvf_smooth` (generalized
variance function smoothing of `vardir`), `step_sae` (stepwise covariate
selection), `compare_sae` (Bland-Altman/CCC agreement, no Pearson r),
`add_reliability_flags` (BPS RSE tiers 25%/50%), and `sim_area_data` /
`sim_series_data` (simulation generators; `sim_spatial_weights` is
internal). Example datasets live in `data/` (`mys`, `mys_proxmat`,
`mys_panel`, `sim_area`, `sim_panel`, `cornsoybean*`) and are documented
in `R/data.R`.

## Tests

testthat edition 3, one file per module in `tests/testthat/`. Coverage
is close to 100% (`test_coverage_gaps.R` exists specifically to hit edge
branches), so new branches need tests. `helper-simulate_safe_sfh.R`
generates well-identified SFH data on a row-standardized grid. Use it
instead of `mys_proxmat` subsets for SFH bootstrap tests, because those
subsets push rho to ±1 and make the SVD fail. Some code comments and dev
notes (`playground/`, test helpers, `main.R`) are written in Indonesian.

## Non-package directories

- `benchmarks/`: speed and accuracy scripts against `sae`, `emdi`, and
  `tipsae`. `bench_frequentist.R` (with `bench_frequentist_child.R`)
  times fastsae, sae and emdi on the same data, each cell in its own
  process with peak RSS; `validation/agreement_reference.R` checks
  agreement with sae/emdi; `sim_mc_finite_pop_d500.R` is the
  design-based simulation; `validation/` also holds the prior and
  INLA-mode studies. Outputs go to `inst/extdata/*.rds` (read by the
  vignettes and the article) and `benchmarks/output/`. After re-running,
  `Rscript benchmarks/update_benchmark_docs.R` refreshes the tables in
  README, the site home and the vignettes.
- `rjournal/`: the R Journal article (`fastsae.Rmd`, `fastsae.bib`,
  `README_RJOURNAL.md` with the build steps and which script feeds which
  table). Verify every new reference (Crossref/DOI) before adding it:
  earlier versions of the bib and Rd files had wrong DOIs and
  non-existent references.
- `paper/`: old LaTeX build files only.
- `playground/`: standalone INLA tutorials, not part of the package.
