# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**fastsae** is an R package providing fast implementations of Small Area Estimation (SAE) methods using C++ (RcppArmadillo). It targets large-scale applications where standard R implementations are too slow—benchmarks show 100x–12,000x speedups over `sae` and `emdi` packages.

## Build & Development Commands

```bash
# Install dependencies and build package
R CMD build . && R CMD INSTALL fastsae_*.tar.gz

# Run full checks (including tests)
R CMD check --as-cran fastsae_*.tar.gz

# Run tests only
Rscript -e "testthat::test_local()"

# Run a specific test file
Rscript -e "testthat::test_file('tests/testthat/test_eblup_fh.R')"

# Compile C++ and regenerate Rcpp exports
Rscript -e "Rcpp::compileAttributes()"

# Run configure (generates src/Makevars for RcppParallel)
./configure

# Clean build artifacts
rm -f src/*.o src/*.so && R CMD build .
```

## Architecture

### R → C++ Boundary

The package follows a standard Rcpp pattern:
- **R wrapper functions** in `R/*.R` handle input validation, formula parsing, and result formatting
- **C++ core functions** in `src/*.cpp` perform the heavy computation
- **Exports** are auto-generated in `src/RcppExports.cpp` and `R/RcppExports.R` via `Rcpp::compileAttributes()`
- Internal C++ functions are prefixed with `.` (e.g., `.eblup_core`, `.seblup_pbmse`)

### Core Computation Modules

| Module | R File | C++ File | Purpose |
|--------|---------|----------|---------|
| Fay-Herriot | `eblup_fh.R` | `eblup_core.cpp` | Area-level EBLUP with analytical MSE |
| Spatial FH | `eblup_sfh.R` | `seblup_core.cpp` | Spatial EBLUP with SAR random effects |
| Bootstrap MSE | — | `seblup_pbmse.cpp`, `seblup_npbmse.cpp` | Parametric/non-parametric bootstrap |
| Spatio-Temporal FH | `eblup_stfh.R` | `eblup_stfh_core.cpp` | Space-time models with AR(1) |
| ST-FH Bootstrap MSE | `pbmse_stfh.R` | `eblup_stfh_pbmse.cpp` | Bootstrap MSE for ST-FH models |
| Unit-level BHF | `eblup_bhf.R` | `eblup_unit.cpp` | Battese-Harter-Fuller model |

### Key Design Patterns

**1. Fisher-Scoring Iteration**: The Fay-Herriot estimation uses a Fisher-scoring algorithm implemented in C++ for speed. The algorithm iterates to estimate the random effect variance (sigma²), then computes EBLUPs.

**2. Spatial Weight Matrix**: Spatial models require a row-standardized proximity matrix `W`. For spatial models, the matrix must cover ALL domains (including unsampled), not just sampled ones—the code handles unsampled areas by borrowing strength from spatial neighbors.

**3. Bootstrap Parallelization**: Parametric and non-parametric bootstrap MSE use OpenMP via RcppParallel. The `n_threads` parameter controls parallelism; `n_threads = 0` uses all available cores.

**4. Unsampled Domain Handling**:
   - Fay-Herriot (`eblup_fh`): Uses synthetic estimators for areas with `NA` response
   - Spatial FH (`eblup_sfh`): Uses spatial kriging (full-spatial synthetic) for unsampled areas
   - Spatio-Temporal FH (`eblup_stfh`): Does NOT support unsampled domains (complete panel required)

**5. Output Structure**: All fitting functions return a list with class `fastsae` containing:
- `df_eblup`: Data frame with estimates (eblup, mse, rse)
- `estcoef`: Regression coefficients with SE, z-value, p-value
- `random_effect_var`: Estimated variance of random effects
- `goodness`: Log-likelihood, AIC, BIC
- `n_iter`, `convergence`: Convergence information

### Data

Built-in datasets in `data/`:
- `mys`: Area-level data with 42 domains (10 unsampled)
- `mys_panel`: Panel version with space-time structure
- `mys_proxmat`: 42×42 spatial proximity matrix
- `cornsoybean`, `cornsoybeanmeans`: Unit-level BHF example data

## Testing Strategy

Tests compare against the `sae` package for correctness. Each model has a dedicated test file:
- `test_eblup_fh.R`: Basic Fay-Herriot
- `test_eblup_sfh.R`: Spatial Fay-Herriot
- `test_eblup_stfh.R`: Spatio-temporal Fay-Herriot
- `test_eblup_bhf.R`: Battese-Harter-Fuller unit-level
- `test_seblup_pbmse.R`, `test_seblup_npbmse.R`: Bootstrap MSE variants
- `test_pbmse_stfh.R`: Bootstrap MSE for ST-FH models

Tests use `skip_if_not_installed("sae")` since they depend on `sae` for reference values.

## Common Tasks

**Add a new estimation method**: Create a new C++ file in `src/`, add the exported function, run `Rcpp::compileAttributes()`, then wrap in an R function in `R/`.

**Modify MSE computation**: MSE formulas are in the C++ files. For parametric bootstrap MSE, see `seblup_pbmse.cpp`; for non-parametric, see `seblup_npbmse.cpp`.

**Change convergence tolerance**: Pass `precision` parameter to the C++ core functions. Default is `1e-4`.
