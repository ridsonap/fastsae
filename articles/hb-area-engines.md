# hb_area() Engines: INLA vs JAX on CPU and Apple GPU (Metal)

[`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
can be estimated by two engines:

- **[`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)**
  (default): integrated nested Laplace approximation with
  [INLA](https://www.r-inla.org).
- **[`jax_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)**:
  full MCMC with NumPyro’s No-U-Turn sampler on
  [JAX](https://jax.readthedocs.io). It runs on the CPU, on an NVIDIA
  GPU (`device = "cuda"`) or on the GPU of Apple silicon
  (`device = "metal"`).

This page compares the speed of the engines and how closely their
estimates agree, on one Apple silicon machine.

## Setup

| Item | Value |
|:---|:---|
| Machine | Darwin arm64 Apple M4 |
| R | R version 4.6.1 (2026-06-24) |
| fastsae | 1.0.0 |
| INLA | 26.8.7 |
| JAX / jax-metal / NumPyro | 0.5.0 / 0.1.1 / 0.19.0 |
| MCMC settings | jax_engine() defaults: 2 chains x 1000 draws, 500 warmup |
| Data | `sim_area_data(D = n, n_unsampled = 0, seed = 42 + n)`; `y ~ x1 + x2` |
| Timing | Each fit run 2 times on the same data; the first run includes JAX compilation |

Each engine runs in a separate R process, because JAX fixes its device
at the first fit of a session. The code is in
[`benchmarks/bench_hb_area_engines.R`](https://github.com/ridsonap/fastsae/blob/main/benchmarks/bench_hb_area_engines.R).

## Computation time

![Warm computation time (log scales). Lower is
faster.](hb-area-engines_files/figure-html/time-plot-1.png)

Warm computation time (log scales). Lower is faster.

| Model | D | INLA | JAX CPU | JAX Metal | JAX Metal, first run | Metal vs CPU |
|:---|---:|---:|---:|---:|---:|---:|
| Beta (proportions) | 100 | 1.3 s | 5.7 s | 119 s | 124 s | 20.8x slower |
| Beta (proportions) | 1,000 | 1.7 s | 104 s | 149 s | 149 s | 1.4x slower |
| Beta (proportions) | 5,000 | 3.2 s | 340 s | 236 s | 256 s | **1.4x faster** |
| Beta (proportions) | 10,000 | 4.8 s | 659 s | 314 s | 328 s | **2.1x faster** |
| Fay-Herriot (Gaussian) | 100 | 1.4 s | 5.3 s | 122 s | 127 s | 22.8x slower |
| Fay-Herriot (Gaussian) | 1,000 | 1.5 s | 37 s | 159 s | 160 s | 4.3x slower |
| Fay-Herriot (Gaussian) | 5,000 | 1.7 s | 207 s | 237 s | 239 s | 1.1x slower |
| Fay-Herriot (Gaussian) | 10,000 | 1.9 s | 606 s | 332 s | 332 s | **1.8x faster** |

The results:

- **INLA is the fastest engine at every size.** It does not sample, so
  its time grows slowly with $`D`$. For example, the Gaussian model at
  $`D = 10,000`$ takes 1.9 s.
- **Metal has a large fixed cost per fit but grows slowly with $`D`$.**
  NUTS builds each trajectory in a loop whose length depends on the
  data, and on Metal every pass through that loop waits for the GPU and
  the CPU to synchronise. For the Gaussian model, going from $`D = 100`$
  to $`D = 10,000`$ multiplies the Metal time by 2.7 and the CPU time by
  113. 
- **CPU or GPU?** The CPU is faster for small and medium numbers of
  domains. The GPU catches up as $`D`$ grows: for the Gaussian model
  Metal is faster at $`D = 10,000`$ and above, and for the Beta model
  Metal is faster at $`D = 5,000`$ and above (of the sizes tested).
  `jax_engine(device = "auto")` therefore keeps the CPU on macOS, and
  `device = "metal"` pays off only for large models.

## Agreement of the estimates

The JAX fits are compared with INLA on the same data. `cor` is the
correlation of the posterior means, `max |diff|` is the largest absolute
difference between them, and `sd ratio` is the median ratio of the
posterior standard deviations (JAX / INLA).

| Model                  | Engine    | D      | cor     | max \|diff\| | sd ratio |
|:-----------------------|:----------|:-------|:--------|-------------:|:---------|
| Beta (proportions)     | JAX CPU   | 100    | 0.99994 |      0.00430 | 1.032    |
| Beta (proportions)     | JAX CPU   | 1,000  | 1.00000 |      0.00100 | 0.956    |
| Beta (proportions)     | JAX CPU   | 5,000  | 1.00000 |      0.00054 | 0.961    |
| Beta (proportions)     | JAX CPU   | 10,000 | 1.00000 |      0.00043 | 0.967    |
| Beta (proportions)     | JAX Metal | 100    | 0.99991 |      0.00340 | 1.032    |
| Beta (proportions)     | JAX Metal | 1,000  | 0.99999 |      0.00140 | 0.950    |
| Beta (proportions)     | JAX Metal | 5,000  | 1.00000 |      0.00049 | 0.946    |
| Beta (proportions)     | JAX Metal | 10,000 | 1.00000 |      0.00044 | 0.943    |
| Fay-Herriot (Gaussian) | JAX CPU   | 100    | 0.99999 |      0.02400 | 1.004    |
| Fay-Herriot (Gaussian) | JAX CPU   | 1,000  | 0.99999 |      0.01600 | 1.000    |
| Fay-Herriot (Gaussian) | JAX CPU   | 5,000  | 0.99999 |      0.02400 | 1.000    |
| Fay-Herriot (Gaussian) | JAX CPU   | 10,000 | 0.99999 |      0.01500 | 0.999    |
| Fay-Herriot (Gaussian) | JAX Metal | 100    | 0.99999 |      0.02000 | 1.007    |
| Fay-Herriot (Gaussian) | JAX Metal | 1,000  | 0.99999 |      0.01500 | 1.002    |
| Fay-Herriot (Gaussian) | JAX Metal | 5,000  | 0.99999 |      0.01600 | 0.999    |
| Fay-Herriot (Gaussian) | JAX Metal | 10,000 | 0.99999 |      0.01800 | 0.999    |

The posterior means of both MCMC engines agree with INLA, with
correlations of at least 0.9999. The posterior standard deviations agree
closely for the Gaussian model (ratio 1.00 to 1.01). For the Beta model
they differ from INLA’s by a few percent (ratio 0.94 to 1.03, smaller
from $`D = 1,000`$ on), on the CPU and on Metal alike. This difference
comes from the method and the priors, which are not identical in the two
engines, and not from the device. The GPU evaluates floating-point
operations in a different order from the CPU, so the Metal and CPU
chains differ slightly, but they give the same estimates. The largest
split R-hat over all JAX fits is 1.009.

## Using the Apple GPU

``` r

# One-time setup of the Python environment with jax-metal
fastsae::setup_numpyro_env(device = "metal")

# In a new R session: request Metal before any other JAX fit
fit <- hb_area(y ~ x1 + x2, vardir = "vardir", data = df,
               engine = jax_engine(device = "metal"))
fit$device   # "metal:METAL:0"
```

jax-metal does not implement two operations used by NumPyro’s NUTS
sampler: the bit count `popcnt`, and some `lax.cond` statements with
constant outputs. When `device = "metal"`, fastsae replaces both with
equivalent operations that Metal supports, and samples the chains
sequentially. If a model still cannot be compiled on Metal, the fit
falls back to the CPU with a warning.

In short: for most surveys (hundreds to a few thousand domains),
[`inla_engine()`](https://ridsonap.github.io/fastsae/reference/hb_engine.md)
is the fastest choice, and `jax_engine(device = "cpu")` is the faster
way to get full MCMC. `jax_engine(device = "metal")` becomes worthwhile
for very large numbers of domains.
