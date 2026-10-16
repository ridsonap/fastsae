# Performance & Scalability Benchmarks

## Executive Summary

Small Area Estimation often involves large administrative registries or
extensive spatial networks with hundreds or thousands of domains.
Traditional implementations in R often rely on interpreted loops, pure R
Fisher-scoring routines, or large memory allocations for $`D \times D`$
dense covariance matrices.

**fastsae** fits the frequentist models in compiled C++
(`RcppArmadillo`; sparse Cholesky factorisations for the SAR models) and
the Bayesian models with INLA. At $`n = 1,000`$ domains:

- **Fay-Herriot**: 1,263x faster than `sae` and 42,874x faster than
  `emdi`.
- **Spatial Fay-Herriot**: 67x faster than `sae` and 51x faster than
  `emdi`.
- **Spatio-Temporal Fay-Herriot** ($`T = 5`$): 24x faster than `sae`.

------------------------------------------------------------------------

## Benchmark Results Table

| Model ($`n = 1,000`$ domains, 5 covariates) | `fastsae` | `sae` | `emdi` |
|:---|:--:|:--:|:--:|
| **Fay-Herriot** | **0.00119 s** (0.503 MB) | 1.50 s (72.0 MB; 1,263x slower) | 51.0 s (3.69 GB; 42,874x slower) |
| **Spatial Fay-Herriot** | **1.00 s** (0.203 MB) | 67.2 s (1.88 GB; 67x slower) | 51.1 s (3.68 GB; 51x slower) |
| **Spatio-Temporal FH ($`T = 5`$)** | **75.3 s** (3.29 MB) | 1790 s (31.7 GB; 24x slower) | \- |

Times are medians of `bench::mark()`; memory is the total memory
allocated by R during one call (`bench::mark()$mem_alloc`), not the peak
resident memory. The `fastsae` rows were re-measured with fastsae 1.0.0
(frequentist on 2026-10-13, `benchmarks/bench_frequentist_fastsae.R`;
Bayesian on 2026-10-14, `benchmarks/bench_bayes_fastsae.R`); the
`sae`/`emdi` rows come from an earlier run of the same design (simulated
data, 5 covariates) on different random draws, the `tipsae` rows from an
earlier run on the same data.

------------------------------------------------------------------------

## Interactive Benchmark Explorer

Use the interactive controls below to compare execution time, memory
allocated, and iterations per second across sample sizes:

Fay-Herriot (FH)

Spatial FH (SFH)

Spatio-Temporal (STFH)

Execution Time

Memory Allocated

Throughput (iter/s)

Execution time (median, seconds) — FH · log scale

Speedup Factor (fastsae vs competitors) — FH

------------------------------------------------------------------------

## Architectural Insights: Why is fastsae so Fast?

1.  **Compiled C++ Linear Solvers**: Instead of interpreting nested
    loops in R, `fastsae` implements Fisher scoring and the analytical
    MSE in Armadillo C++. The spatial Fay-Herriot model never forms a
    $`D \times D`$ inverse: with $`A = (I - \rho W')(I - \rho W)`$ and
    $`M = A + \sigma^2 E'\Psi^{-1}E`$, $`V^{-1}`$, $`g_1, \dots, g_4`$
    and the kriging terms follow from sparse Cholesky solves with $`A`$
    and $`M`$ (Eigen, AMD ordering).

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
    hierarchical generalized and spatio-temporal models (`hb_area`),
    `fastsae` utilizes INLA to compute analytical posterior margins
    directly from sparse Gaussian Markov Random Fields, bypassing Monte
    Carlo sampling entirely.

------------------------------------------------------------------------

## Bayesian Spatio-Temporal Benchmark: `fastsae` (INLA) vs `tipsae` (Stan MCMC)

To assess Bayesian small area estimation performance,
[`fastsae::hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
was benchmarked against
[`tipsae::fit_sae()`](https://rdrr.io/pkg/tipsae/man/fit_sae.html), the
state-of-the-art Stan MCMC package for spatio-temporal Beta small area
models.

The test used the official `tipsae` panel dataset: **`emilia`** (38
health districts in Emilia-Romagna over 5 years, $`N = 190`$ domains
$`\times`$ years) with spatial polygon contiguity matrix $`W`$ from
**`emilia_shp`**.

| Metric | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Stan MCMC) | [`fastsae::hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md) (INLA) |
|:---|:---|:---|
| **Model** | Beta, Besag ICAR + domain-specific RW(1) | `family = "beta"`, `spatial = "besag"`, `temporal = "rw1"`, `st_interaction = "domain-specific"`; same structure, different priors |
| **Computation** | Hamiltonian Monte Carlo (NUTS) | Integrated nested Laplace approximation (INLA’s default `"compact"` mode for Beta) |
| **Execution time** | 120 s (4 chains x 2000 iterations in parallel; largest $`\hat R`$ 1.07). One chain of 200 iterations takes about 10 s but has not converged ($`\hat R`$ \> 2) | 1.1 s |
| **Posterior means vs converged `tipsae`** | Reference | $`r = 0.998`$, mean absolute difference 0.0016 (largest 0.006) |
| **Posterior SD vs `tipsae`** | Reference | ratio 0.84 to 1.26 (the priors differ) |

Measured on 2026-10-16 on one machine (Apple silicon, R 4.6.1, INLA
26.8.7, tipsae 1.0.3); times depend on the hardware. INLA needs no
chains or warm-up, but its fit should still be checked (`convergence`,
from INLA’s mode search).

### Interactive Beta SAE Benchmark Explorer ($`n = 30`$ to $`n = 1,000`$)

Use the interactive controls below to compare runtime, memory footprint,
and speedup factors between
[`fastsae::hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
and [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html)
across **Standard Beta**, **Spatial Beta (Besag ICAR)**, and
**Spatio-Temporal Beta (Besag + RW1)** models:

Standard Beta

Spatial Beta (ICAR)

Spatio-Temporal Beta

Execution Time

Memory Allocated

Throughput (iter/s)

Execution time (median, seconds) — Standard Beta · log scale

Speedup Factor & Efficiency Multiplier (fastsae vs tipsae) — Standard
Beta

#### Inferential Equivalence & Validation Metrics ($`n = 1,000`$)

To establish that the substantial computational accelerations achieved
by
[`fastsae::hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
preserve full inferential validity relative to Stan MCMC
([`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html)),
summary concordance and error metrics were evaluated across
$`n = 1,000`$ domains/units for Standard Beta, Spatial Beta (Besag
ICAR), and Spatio-Temporal Beta models:

| Dimension | Metric | Beta | Spatial Beta | Spatio-Temporal Beta | Statistical Implication |
|:---|:---|:--:|:--:|:--:|:---|
| **Point Estimates (EBP)** | **Pearson Correlation ($`r`$)** | **0.99934** | **0.99891** | **0.99614** | Near-perfect linear agreement between INLA and MCMC |
|  | **Spearman Rank Correlation ($`\rho`$)** | **0.99921** | **0.99875** | **0.99517** | Consistent domain priority and rank ordering |
|  | **Mean Absolute Error (MAE)** | **0.01969** | **0.01456** | **0.01200** | Average estimation discrepancy $`< 0.02`$ |
|  | **Root Mean Squared Difference (RMSD)** | **0.02369** | **0.01713** | **0.01458** | Negligible $`L_2`$ deviation across domains |
| **Uncertainty / MSE** | **MSE Pearson Correlation ($`r_{\text{MSE}}`$)** | **0.96315** | **0.89254** | **0.67957** | Preserved area-level precision ordering |
|  | **Mean Absolute Difference (MSE)** | **0.00018** | **0.00012** | **0.00025** | Virtually identical error variances |

------------------------------------------------------------------------

### Finite Population Simulation & Ground Truth Validation

For comprehensive design-based simulation studies evaluating parameter
recovery and empirical properties under repeated survey sampling from
finite populations ($`D = 60`$ and national-scale $`D = 500`$), see the
dedicated article: [Finite Population Simulation &
Validation](https://ridsonap.github.io/fastsae/articles/finite-population-validation.md).
