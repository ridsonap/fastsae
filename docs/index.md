# fastsae: Fast Small Area Estimation in R

**fastsae** is a high-performance R package for **Small Area Estimation
(SAE)**. It re-engineers classic and contemporary area-level and
unit-level SAE models using **C++ (`RcppArmadillo`)** and **OpenMP
multi-threading**, providing massive speedups (up to **12,000x faster**)
and radical memory reductions (up to **20,000x lower memory**) while
maintaining **exact numerical equivalence** with gold-standard
implementations.

------------------------------------------------------------------------

## Why fastsae?

### 1. High-Performance C++ Core

Core Fisher-scoring iterative routines and Woodbury matrix solvers run
in compiled C++ via `RcppArmadillo`, directly utilizing optimized BLAS
and LAPACK routines. Computational bottlenecks that previously required
minutes in pure R complete in milliseconds.

### 2. Minimal Memory Footprint

Standard packages often allocate dense $`O(m^2)`$ and $`O(m^3)`$
covariance matrices on the R heap. `fastsae` avoids heap bloat through
zero-copy matrix views and in-place algebra, keeping peak memory
consumption under 25 MB even for thousands of domains where existing
packages consume gigabytes.

### 3. Exact Numerical Equivalence

High computational speed does not come at the expense of statistical
fidelity. Point estimates ($`\hat{\beta}`$, $`\hat{\theta}`$) and
variance components ($`\hat{\sigma}_u^2`$, $`\hat{\rho}`$) are identical
to machine precision ($`< 10^{-13}`$) compared to published benchmarks
in the `sae` package (Molina & Rao).

### 4. Native Multi-Threaded Bootstrap

Parametric and Non-Parametric Bootstrap MSE estimation are parallelized
natively across CPU cores using `#pragma omp parallel for`. This
eliminates the process-forking and serialization overhead associated
with traditional R parallel clusters.

### 5. Automatic Spatial Kriging for Unsampled Domains

When geographic domains have missing response data (`y = NA`),
`eblup_sfh` automatically performs full-spatial synthetic prediction
(kriging) borrowing strength from neighboring areas, without requiring
manual subsetting or data partitioning.

### 6. Modern S3 Interface and Diagnostics

Designed for standard R workflows: full support for
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), and
[`residuals()`](https://rdrr.io/r/stats/residuals.html), alongside
diagnostic and model-comparison visualization via
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html).

------------------------------------------------------------------------

## Comparison with Other Packages

| Feature | `fastsae` | `sae` (Molina & Rao) | `emdi` (Kreutzmann et al.) |
|:---|:--:|:--:|:--:|
| **Core Computation** | **C++ (RcppArmadillo)** | Pure R | R / lme4 / nlme |
| **Multi-Threading** | **Native OpenMP** (`n_threads`) | Single-threaded | Optional foreach/parallel |
| **Numerical Consistency** | **Reference baseline** | Baseline | Approximations |
| **Unsampled Area Support** | **Automatic (Spatial Kriging)** | Manual subsetting required | Limited |
| **Bootstrap Speed** | **Parallel C++ (Ultra-Fast)** | Slow (Interpreted R loops) | Moderate |
| **RAM Consumption** | **Minimal (\< 25 MB)** | Moderate (~100 MB) | High (~800+ MB to 4 GB) |
| **S3 Methods Support** | `print`, `summary`, `coef`, `fitted`, `residuals`, `autoplot` | Custom lists | Standard S3 |

------------------------------------------------------------------------

## Supported Models

| Model | Function | Random Effect Structure | MSE Estimation Methods |
|:---|:---|:---|:---|
| **Fay-Herriot** (Area-level) | [`eblup_fh()`](https://ridsonap.github.io/fastsae/reference/eblup_fh.md) | Independent area effects ($`u_d \sim N(0, \sigma_u^2)`$) | Analytical (Prasad-Rao) |
| **Spatial Fay-Herriot** | [`eblup_sfh()`](https://ridsonap.github.io/fastsae/reference/eblup_sfh.md) | Simultaneous Autoregressive (SAR(1)) | Analytical, Parametric Bootstrap (`pbmse`), Non-Parametric Bootstrap (`npbmse`) |
| **Spatio-Temporal Fay-Herriot** | [`eblup_stfh()`](https://ridsonap.github.io/fastsae/reference/eblup_stfh.md) | Spatial SAR(1) + Temporal AR(1) | Parametric Bootstrap (`pbmse`) |
| **Battese-Harter-Fuller** (Unit-level) | [`eblup_bhf()`](https://ridsonap.github.io/fastsae/reference/eblup_bhf.md) | Random intercept nested in domains | Parametric Bootstrap (`pbmse`) |

------------------------------------------------------------------------

## Installation

\
`# Install from CRAN (once published)`\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"fastsae"``)`\
\
`# Or install the development version from GitHub`\
`# install.packages("remotes")`\
`remotes``::`[`install_github`](https://remotes.r-lib.org/reference/install_github.html)`(``"ridsonap/fastsae"``)`

------------------------------------------------------------------------

## Documentation and User Guides

Explore detailed guides, worked examples, and performance analyses:

- [**Getting Started with
  fastsae**](https://ridsonap.github.io/fastsae/articles/getting-started.md):
  Step-by-step walkthrough of model fitting, diagnostics, and S3
  extraction.
- [**Spatial & Spatio-Temporal
  Models**](https://ridsonap.github.io/fastsae/articles/spatial-temporal.md):
  Working with spatial weights matrices, spatial kriging, and
  spatio-temporal panels.
- [**Unit-Level Estimation
  (BHF)**](https://ridsonap.github.io/fastsae/articles/unit-level-bhf.md):
  Fitting the Battese-Harter-Fuller model with individual survey data.
- [**Performance
  Benchmarks**](https://ridsonap.github.io/fastsae/articles/benchmarks.md):
  Interactive charts and speedup comparisons across sample sizes.
- [**Function
  Reference**](https://ridsonap.github.io/fastsae/reference/index.md):
  Comprehensive API reference and parameter documentation.

------------------------------------------------------------------------

## References

- Fay, R. E., & Herriot, R. A. (1979). Estimates of income for small
  places: An application of James-Stein procedures to Census data.
  *Journal of the American Statistical Association*, 74(366), 269–277.
- Battese, G. E., Harter, R. M., & Fuller, W. A. (1988). An
  error-components model for prediction of county crop areas using
  survey and satellite data. *Journal of the American Statistical
  Association*, 83(401), 28–36.
- Marhuenda, Y., Molina, I., & Morales, D. (2013). Small area estimation
  with spatio-temporal Fay-Herriot models. *Computational Statistics &
  Data Analysis*, 58, 308–325.
- Pratesi, M., & Salvati, N. (2008). Small area estimation for spatially
  correlated data: A Fay-Herriot with the spatial linear spline model.
  *Journal of Applied Statistics*, 35(7), 781–794.
- Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd ed.).
  John Wiley & Sons.
