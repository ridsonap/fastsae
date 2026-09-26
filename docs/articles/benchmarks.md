# Performance & Scalability Benchmarks

## Executive Summary

Small Area Estimation often involves large administrative registries or
extensive spatial networks with hundreds or thousands of domains.
Traditional implementations in R often rely on interpreted loops, pure R
Fisher-scoring routines, or large memory allocations for $`D \times D`$
dense covariance matrices.

**fastsae** re-architects core estimation algorithms in compiled C++
using `RcppArmadillo` and native `OpenMP` multi-threading, delivering:

- **Up to 364x faster** than `sae` and **12,300x faster** than `emdi`
  for Fay-Herriot models at $`n = 1,000`$.
- **Up to 80x faster** than `sae` for Spatial Fay-Herriot models at
  $`n = 1,000`$.
- **Up to 31x faster** than `sae` for Spatio-Temporal models at
  $`n = 1,000`$.
- **Massive RAM reduction**: Peak memory footprint stays under **25 MB**
  where existing packages require hundreds of megabytes or several
  gigabytes.

------------------------------------------------------------------------

## Benchmark Results Table

The following benchmarks were conducted on simulated datasets across
domain sizes ranging from $`n = 30`$ to $`n = 1,000`$ (with 5 auxiliary
covariates):

| Metric | `fastsae` | `sae` (Molina & Rao) | `emdi` (Kreutzmann et al.) |
|:---|:--:|:--:|:--:|
| **Mean Time (EBLUP FH)** | **0.0015 s** | 0.291 s | 9.64 s |
| **Mean Time (Spatial FH)** | **0.165 s** | 12.60 s | 8.69 s |
| **Mean Time (Spatio-Temporal FH)** | **12.5 s** | 353.0 s | \- |
| **Peak Memory (EBLUP FH)** | **0.055 MB** | 16.3 MB | 824 MB |
| **Peak Memory (Spatial FH)** | **5.15 MB** | 408 MB | 824 MB |
| **Peak Memory (Spatio-Temporal FH)** | **0.289 MB** | 7,822 MB | \- |
| **Speedup at n = 1,000 (FH)** | **Baseline** | **~364x slower** | **~12,300x slower** |

------------------------------------------------------------------------

## Interactive Benchmark Explorer

Use the interactive controls below to compare execution time, RAM
consumption, and iterations per second across sample sizes:

Fay-Herriot (FH)

Spatial FH (SFH)

Spatio-Temporal (STFH)

Execution Time

Peak Memory

Throughput (iter/s)

Execution time (median, seconds) — FH · log scale

Speedup Factor (fastsae vs competitors) — FH

------------------------------------------------------------------------

## Architectural Insights: Why is fastsae so Fast?

1.  **Compiled C++ Linear Solvers**: Instead of interpreting nested
    loops in R, `fastsae` implements Fisher-scoring parameter search and
    Woodbury identity matrix inversions in Armadillo C++, directly
    leveraging optimized BLAS/LAPACK routines.

2.  **Zero-Copy Matrix Operations**: Memory allocations for large
    intermediate structures ($`V`$, $`V^{-1}`$, $`P`$) are avoided or
    reused across Fisher iterations rather than reallocated on the heap.

3.  **OpenMP Multi-Threaded Bootstrap**: Bootstrap resampling runs
    natively in parallel across CPU cores using
    `#pragma omp parallel for`, avoiding the serialization overhead of R
    worker processes
    ([`parallel::makeCluster`](https://rdrr.io/r/parallel/makeCluster.html)
    / `foreach`).

4.  **Woodbury Identity for Panel Data**: In `eblup_stfh`, inversion of
    the block $`(DT \times DT)`$ covariance matrix $`V`$ is reduced to
    operations on individual $`D \times D`$ and $`T \times T`$ blocks
    via Kronecker and Woodbury decomposition, transforming an
    $`O((DT)^3)`$ bottleneck into scalable $`O(D^3 + T^3)`$ steps.

5.  **Integrated Nested Laplace Approximations (INLA)**: For complex
    hierarchical generalized and spatio-temporal models (`ebp_area`),
    `fastsae` utilizes INLA to compute analytical posterior margins
    directly from sparse Gaussian Markov Random Fields, bypassing Monte
    Carlo sampling entirely.

------------------------------------------------------------------------

## Bayesian Spatio-Temporal Benchmark: `fastsae` (INLA) vs `tipsae` (Stan MCMC)

To assess Bayesian small area estimation performance,
[`fastsae::ebp_area()`](https://ridsonap.github.io/fastsae/reference/ebp_area.md)
was benchmarked against
[`tipsae::fit_sae()`](https://rdrr.io/pkg/tipsae/man/fit_sae.html), the
state-of-the-art Stan MCMC package for spatio-temporal Beta small area
models.

The test used the official `tipsae` panel dataset: **`emilia`** (38
health districts in Emilia-Romagna over 5 years, $`N = 190`$ domains
$`\times`$ years) with spatial polygon contiguity matrix $`W`$ from
**`emilia_shp`**.

| Metric | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Stan MCMC) | [`fastsae::ebp_area`](https://ridsonap.github.io/fastsae/reference/ebp_area.md) (INLA) | Advantage |
|:---|:---|:---|:---|
| **Model Specification** | Besag ICAR + Domain RW(1) | `spatial = "besag"`, `temporal = "rw1"`, `st_interaction = "domain-specific"` | **Exact structural equivalence** |
| **Computational Engine** | Hamiltonian Monte Carlo (NUTS Stan) | Integrated Nested Laplace Approximation | Analytical & deterministic |
| **Execution Time** | **15.8 s** (1 chain, 200 iter) / ~120 s (4 chains) | **1.93 s** (`simplified.laplace`) | **~9x to 50x+ faster** |
| **Pearson Correlation ($`r`$)** | Baseline | **0.9893** | **Near-identical point estimates** |
| **Mean Absolute Error (MAE)** | Baseline | **0.00304** | **Negligible numerical error** |
| **Convergence Overhead** | Requires $`\hat{R} < 1.05`$ checks, warmup, and tuning | None (closed-form Laplace expansions) | Instant convergence |

### Interactive Beta SAE Benchmark Explorer ($`n = 30`$ to $`n = 1,000`$)

Use the interactive controls below to compare runtime, memory footprint,
and speedup factors between
[`fastsae::ebp_area`](https://ridsonap.github.io/fastsae/reference/ebp_area.md)
and [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html):

Execution Time

Peak Memory

Throughput (iter/s)

**Area-Level Beta SAE**: INLA vs Stan MCMC

Execution time (median, seconds) — Beta SAE · log scale

Speedup Factor & Efficiency Multiplier (fastsae vs tipsae) — Beta SAE

#### Comparison of Area Estimates & MSE ($`n = 1,000`$)

The table below shows the first 10 domain point estimates and
corresponding MSE from
[`fastsae::ebp_area`](https://ridsonap.github.io/fastsae/reference/ebp_area.md)
and [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html)
across $`n = 1,000`$ areas (overall correlation $`r = 0.9993`$):

| area | estimasi fastsae | mse fastsae | estimasi tipsae | mse tipsae |
|:-----|:-----------------|:------------|:----------------|:-----------|
| 1    | 0.54531          | 0.00019     | 0.53319         | 0.00009    |
| 2    | 0.62443          | 0.00026     | 0.59912         | 0.00012    |
| 3    | 0.52038          | 0.00042     | 0.51852         | 0.00020    |
| 4    | 0.57258          | 0.00029     | 0.55830         | 0.00013    |
| 5    | 0.58710          | 0.00041     | 0.57171         | 0.00020    |
| 6    | 0.64299          | 0.00030     | 0.61513         | 0.00013    |
| 7    | 0.25747          | 0.00044     | 0.29868         | 0.00017    |
| 8    | 0.36546          | 0.00040     | 0.39219         | 0.00016    |
| 9    | 0.48903          | 0.00022     | 0.48951         | 0.00011    |
| 10   | 0.62417          | 0.00040     | 0.59155         | 0.00017    |
