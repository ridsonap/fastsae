# Bayesian Area-Level Small Area Estimation with INLA

## Introduction

Classical Small Area Estimation (SAE) models—such as the Fay-Herriot
(1979) model—traditionally assume a Gaussian distribution for direct
survey estimators with known sampling error variances. However,
real-world official statistics frequently encounter non-Gaussian
indicators: - **Poverty rates, unemployment rates, and prevalence
ratios**: strictly bounded in the open interval $`(0, 1)`$. - **Disease
counts, crime incidents, and rare events**: discrete counts subject to
Poisson or Negative Binomial overdispersion with population exposures. -
**Binary survey aggregates**: proportions derived from small sample
sizes (Binomial trials). - **Skewed economic expenditures**: strictly
positive continuous outcomes (Gamma).

Furthermore, when data are collected across geographical areas and
repeated over time (panel/longitudinal surveys), outcomes exhibit both
spatial dependency and temporal autocorrelation.

The **fastsae** package introduces
[`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md),
a unified hierarchical Bayesian framework powered by **Integrated Nested
Laplace Approximations (INLA)**. It provides fast, deterministic, and
highly accurate analytical posterior approximations without the
computational bottlenecks, convergence diagnostics, or chain-tuning
associated with Markov Chain Monte Carlo (MCMC).

------------------------------------------------------------------------

## Supported Distribution Families

[`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
supports 6 probability distributions for the direct response:

| Family | Typical SAE Use Case | Required Sampling Input | Model Scale Parameter |
|:---|:---|:---|:---|
| **`"gaussian"`** | Continuous indicators (e.g., mean income, mean years of schooling) | `vardir` (sampling variance $`D_d`$) | Inverted sampling variance $`\text{scale} = 1/D_d`$ |
| **`"beta"`** | Rates and proportions strictly in $`(0, 1)`$ | `vardir` or `trials` | Sampling dispersion via Janicki (2020): $`\phi_d = \frac{y_d(1 - y_d)}{V_d} - 1`$ |
| **`"poisson"`** | Counts and rare event rates | `exposure` (expected count / population offset) | Log-link offset $`\log(E_d)`$ |
| **`"nbinomial"`** | Overdispersed counts and clustered events | `exposure` | Log-link offset with negative binomial dispersion |
| **`"binomial"`** | Aggregated binary proportions | `trials` (domain sample size $`n_d`$) | Number of trials per domain |
| **`"gamma"`** | Positive, right-skewed indicators (e.g., expenditures) | `vardir` | Shape-rate parameterization with known sampling variance |

------------------------------------------------------------------------

## Mathematical Formulation

The hierarchical Bayesian model underlying
[`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
extends the classical Fay-Herriot framework to handle non-Gaussian
responses and structured random effects. The model comprises three
levels:

### Level 1: Sampling Model (Likelihood)

For domain $`d`$ with direct estimator $`y_d`$ and known sampling
variance $`D_d`$:

``` math
\pi(y_d \mid \eta_d) = \text{ExpFamily}(y_d \mid \eta_d, D_d)
```

The link function connects the linear predictor $`\eta_d`$ to the mean
$`\mu_d = E(y_d)`$:

| Family | Link Function | Linear Predictor |
|:---|:---|:---|
| Gaussian | Identity: $`\eta_d = \mu_d`$ | $`\eta_d = x_d^\top \beta + u_d`$ |
| Beta | Logit: $`\eta_d = \log(\mu_d / (1-\mu_d))`$ | $`\eta_d = x_d^\top \beta + u_d`$ |
| Poisson | Log: $`\eta_d = \log(\mu_d)`$ | $`\eta_d = x_d^\top \beta + u_d + \log(E_d)`$ |
| Binomial | Logit: $`\eta_d = \log(\mu_d / (1-\mu_d))`$ | $`\eta_d = x_d^\top \beta + u_d`$ |
| Negative Binomial | Log: $`\eta_d = \log(\mu_d)`$ | $`\eta_d = x_d^\top \beta + u_d + \log(E_d)`$ |
| Gamma | Log: $`\eta_d = \log(\mu_d)`$ | $`\eta_d = x_d^\top \beta + u_d`$ |

### Level 2: Prior on Fixed Effects

``` math
\beta \sim N(0, \tau_\beta^2 I_p)
```

A weakly informative normal prior with large variance $`\tau_\beta^2`$
is placed on regression coefficients.

### Level 3: Prior on Random Effects

``` math
u \sim N(0, \Sigma)
```

The random effects $`u = (u_1, \ldots, u_m)^\top`$ follow a zero-mean
Gaussian Markov Random Field (GMRF) with precision matrix
$`\Sigma^{-1}`$ encoding the chosen spatial or temporal structure (BYM2,
Besag ICAR, RW1, etc.).

### INLA Approximation

Rather than MCMC sampling, INLA approximates the posterior marginals:

``` math
\pi(\theta \mid y) \approx \tilde{\pi}(\theta \mid y) = \frac{\pi(\theta) \pi(y \mid \theta)}{\pi(y \mid \theta)}
```

and

``` math
\pi(u_d \mid y) \approx \tilde{\pi}(u_d \mid y) = \sum_k \tilde{\pi}(u_d \mid \theta_k, y) \tilde{\pi}(\theta_k \mid y) \Delta_k
```

where the integrals are evaluated using **Laplace approximations** and
**numerical integration**. Three computational strategies are available:

- **`strategy = "gaussian"`**: Fastest, best for well-behaved problems
- **`strategy = "laplace"`**: Default, good accuracy/speed balance
- **`strategy = "simplified.laplace"`**: Slightly faster Laplace, good
  for large models

------------------------------------------------------------------------

## Spatial and Spatio-Temporal Structures

### Spatial Priors

1.  **BYM2 (`spatial = "bym2"`)**: The scaled Besag-York-Mollié model
    (Simpson et al., 2017) decomposes area random effects into a spatial
    structured component and an unstructured IID component:
    ``` math
    u_d = \frac{1}{\sqrt{\tau_u}} \left( \sqrt{1 - \phi} \, v_d + \sqrt{\phi} \, u_* d \right)
    ```
    where $`\phi \in [0, 1]`$ measures the proportion of variance
    explained by spatial structure, regulated by Penalized Complexity
    (PC) priors.
2.  **Besag ICAR (`spatial = "besag"`)**: Intrinsic Conditional
    Autoregressive spatial model.
3.  **Spatial Lag Model (`spatial = "slm"`)**: Simultaneous spatial
    autoregression $`(I - \rho W)^{-1}`$.

### Temporal Models

When data are observed longitudinally across $`T \ge 2`$ time periods,
specify `time`: - **`temporal = "rw1"`**: First-order Random Walk
($`\Delta t_\tau \sim N(0, \sigma_t^2)`$). - **`temporal = "rw2"`**:
Second-order Random Walk for smooth non-linear temporal trends. -
**`temporal = "ar1"`**: First-order Autoregressive process with
autocorrelation $`|\rho_t| < 1`$. - **`temporal = "iid"`**: Unstructured
time random effects.

### Space-Time Interaction Structures (`st_interaction`)

- **`"none"` (Additive)**: Separate spatial main effect and temporal
  trend:
  ``` math
  \eta_{dt} = x_{dt}^\top \beta + s_d + \gamma_t
  ```
- **`"domain-specific"`** (matching `tipsae`): Domain-specific random
  walks / AR(1) processes with shared precision:
  ``` math
  \eta_{dt} = x_{dt}^\top \beta + s_d + t_{d, t}, \quad t_{d, \cdot} \sim \text{RW1}(\sigma_t^2)
  ```
- **`"separable"`** (matching classical `eblup_stfh`): Spatial random
  field evolving dynamically across time with AR(1) temporal
  correlation.
- **`"type1"` to `"type4"`**: Full Knorr-Held (2000) classifications of
  space-time interactions.

------------------------------------------------------------------------

## Case Study: Spatio-Temporal Beta SAE on Emilia-Romagna Poverty

We demonstrate
[`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
using the Italian poverty dataset `emilia` from the `tipsae` package (38
health districts in Emilia-Romagna across 5 years, $`N = 190`$). The
goal is to estimate the Head Count Ratio (`hcr`), which represents the
proportion of households below the poverty line.

### 1. Load Data and Construct Spatial Adjacency Matrix

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`fastsae`](https://ridsonap.github.io/fastsae/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(``tipsae``)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`spdep`](https://github.com/r-spatial/spdep/)`)`\
\
[`data`](https://rdrr.io/r/utils/data.html)`(``"emilia"``)`\
[`data`](https://rdrr.io/r/utils/data.html)`(``"emilia_shp"``)`\
\
`# Construct binary spatial adjacency matrix W from district polygons`\
`nb`` ``<-`` `[`poly2nb`](https://r-spatial.github.io/spdep/reference/poly2nb.html)`(``emilia_shp``)`\
`W`` ``<-`` `[`nb2mat`](https://r-spatial.github.io/spdep/reference/nb2mat.html)`(``nb``, style ``=`` ``"B"``, zero.policy ``=`` ``TRUE``)`\
[`rownames`](https://rdrr.io/r/base/colnames.html)`(``W``)`` ``<-`` `[`colnames`](https://rdrr.io/r/base/colnames.html)`(``W``)`` ``<-`` `[`as.character`](https://rdrr.io/r/base/character.html)`(``emilia_shp``$``NAME_DISTRICT``)`\
\
[`head`](https://rdrr.io/r/utils/head.html)`(``emilia``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"id"``, ``"year"``, ``"hcr"``, ``"vars"``, ``"x"``)``]``)`\
`#>                    id year    hcr         vars       x`\
`#> 1 CASALECCHIO DI RENO 2014 0.0404 9.090478e-05 -0.2624`\
`#> 2   CITTA' DI BOLOGNA 2014 0.0825 6.404001e-05 -0.0008`\
`#> 3               IMOLA 2014 0.1033 3.120275e-04 -0.0522`\
`#> 4         PIANURA EST 2014 0.0633 1.025764e-04 -0.4007`\
`#> 5       PIANURA OVEST 2014 0.0625 1.562500e-04 -0.2277`\
`#> 6      PORRETTA TERME 2014 0.1276 6.643609e-04 -0.4434`

### 2. Fit Spatio-Temporal Beta Model

We model the poverty proportion using a Beta likelihood with a Besag
ICAR spatial random effect, a domain-specific RW(1) temporal dynamic,
and known sampling variances `vars`:

\
`fit_beta_st`` ``<-`` `[`hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md)`(`\
`  formula ``=`` ``hcr`` ``~`` ``x``,`\
`  data ``=`` ``emilia``,`\
`  domain ``=`` ``"id"``,`\
`  time ``=`` ``"year"``,`\
`  vardir ``=`` ``"vars"``,`\
`  family ``=`` ``"beta"``,`\
`  spatial ``=`` ``"besag"``,`\
`  temporal ``=`` ``"rw1"``,`\
`  st_interaction ``=`` ``"domain-specific"``,`\
`  W ``=`` ``W``,`\
`  print_result ``=`` ``FALSE`\
`)`\
\
[`summary`](https://rdrr.io/r/base/summary.html)`(``fit_beta_st``)`\
`#> `\
`#> Variance Components:`\
`#> sigma2_u: 0.054876 `\
`#> sigma2_t (temporal): 0.005803 `\
`#> `\
`#> Coefficients:`\
`#>                    beta   std.error      zvalue      pvalue    ci_lower`\
`#> (Intercept) -2.2517e+00  1.5232e-02 -1.4783e+02  0.0000e+00 -2.2815e+00`\
`#> x            4.0613e-01  6.4372e-02  6.3090e+00  2.8081e-10  2.7947e-01`\
`#>             ci_upper`\
`#> (Intercept)  -2.2218`\
`#> x             0.5321`\
`#> `\
`#> Hyperparameters:`\
`#>                                  mean       sd 0.025quant  0.5quant 0.975quant`\
`#> Precision for ..domain_id..  18.22296  6.62633   8.634528  17.10218    34.3497`\
`#> Precision for ..time_id..   172.31261 78.43491  72.152963 155.44977   373.0499`\
`#>                                  mode`\
`#> Precision for ..domain_id..  15.06372`\
`#> Precision for ..time_id..   127.08334`\
`#> `\
`#> Goodness of Fit:`\
`#>                                                   DIC `\
`#>                                            -976.34160 `\
`#>                                                    pD `\
`#>                                              56.02359 `\
`#>                                                  WAIC `\
`#>                                            -984.78963 `\
`#>                                                 pWAIC `\
`#>                                              38.70809 `\
`#> Marginal_LogLik.log marginal-likelihood (integration) `\
`#>                                             439.32543 `\
`#> `\
`#> HB Summary Statistics:`\
`#>        hb           linear_pred           sd                mse           `\
`#>  Min.   :0.05270   Min.   :-2.894   Min.   :0.005185   Min.   :2.688e-05  `\
`#>  1st Qu.:0.08271   1st Qu.:-2.413   1st Qu.:0.007702   1st Qu.:5.932e-05  `\
`#>  Median :0.09635   Median :-2.244   Median :0.009371   Median :8.781e-05  `\
`#>  Mean   :0.09802   Mean   :-2.251   Mean   :0.009566   Mean   :9.767e-05  `\
`#>  3rd Qu.:0.11011   3rd Qu.:-2.095   3rd Qu.:0.011051   3rd Qu.:1.221e-04  `\
`#>  Max.   :0.15657   Max.   :-1.687   Max.   :0.016768   Max.   :2.812e-04  `\
`#>       rse        `\
`#>  Min.   : 6.370  `\
`#>  1st Qu.: 8.725  `\
`#>  Median : 9.813  `\
`#>  Mean   : 9.827  `\
`#>  3rd Qu.:10.979  `\
`#>  Max.   :13.670`

### 3. Inspect Posterior Predictions and Uncertainty

The output object contains full domain-period predictions, posterior
standard errors, and 95% Credible Intervals:

\
[`head`](https://rdrr.io/r/utils/head.html)`(``fit_beta_st``$``df_hb``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"domain"``, ``"time"``, ``"y"``, ``"hb"``, ``"sd"``, ``"rse"``, ``"ci_lower"``, ``"ci_upper"``)``]``)`\
`#>                domain time      y         hb          sd       rse   ci_lower`\
`#> 1 CASALECCHIO DI RENO 2014 0.0404 0.05554648 0.006127591 11.031465 0.04408737`\
`#> 2   CITTA' DI BOLOGNA 2014 0.0825 0.08225896 0.006131555  7.453966 0.07068304`\
`#> 3               IMOLA 2014 0.1033 0.09343257 0.010214708 10.932706 0.07511426`\
`#> 4         PIANURA EST 2014 0.0633 0.06509829 0.006294779  9.669652 0.05342028`\
`#> 5       PIANURA OVEST 2014 0.0625 0.07112657 0.007870648 11.065693 0.05671277`\
`#> 6      PORRETTA TERME 2014 0.1276 0.09608722 0.012916432 13.442403 0.07314329`\
`#>     ci_upper`\
`#> 1 0.06810573`\
`#> 2 0.09473234`\
`#> 3 0.11517430`\
`#> 4 0.07814089`\
`#> 5 0.08759939`\
`#> 6 0.12376134`

------------------------------------------------------------------------

## Comparison and Validation Against `tipsae` (Stan MCMC)

In `tipsae::fit_sae(..., spatial_error = TRUE, temporal_error = TRUE)`,
the model is estimated using Hamiltonian Monte Carlo (HMC) in Stan. In
[`fastsae::hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md),
the identical structural model is estimated analytically via INLA.

| Characteristic | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Stan MCMC) | [`fastsae::hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md) (INLA) | Advantage |
|:---|:---|:---|:---|
| **Model Specification** | Besag ICAR + Domain RW(1) | `spatial = "besag"`, `temporal = "rw1"`, `st_interaction = "domain-specific"` | Exact structural match |
| **Estimation Method** | MCMC (NUTS) | Integrated Nested Laplace Approximations | Deterministic, no burn-in needed |
| **Execution Time** | ~15.8 s (1 chain, 200 iter) / ~120 s (4 chains) | **1.93 s** (`simplified.laplace`) | **~9x to 50x+ faster** |
| **Pearson Correlation ($`r`$)** | Baseline | **0.9893** | Near-identical estimates |
| **Mean Absolute Error (MAE)** | Baseline | **0.00304** | Negligible numerical discrepancy |

### Interactive Scalability Explorer ($`n = 30`$ to $`n = 1,000`$)

Use the interactive controls below to compare runtime, memory footprint,
and speedup factors between
[`fastsae::hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
and [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html)
across both **Standard Beta** and **Spatial Beta (Besag ICAR)** models:

Standard Beta

Spatial Beta (ICAR)

Spatio-Temporal Beta

Execution Time

Peak Memory

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

## Cross-Sectional Count Data: Poisson & Negative Binomial

For count data such as disease incidence or crime events,
[`hb_area()`](https://ridsonap.github.io/fastsae/reference/hb_area.md)
seamlessly supports Poisson and Negative Binomial distributions with
population offsets:

\
`# Simulate 42-domain count data with known spatial structure`\
[`data`](https://rdrr.io/r/utils/data.html)`(``"sim_area"``, package ``=`` ``"fastsae"``)`\
[`data`](https://rdrr.io/r/utils/data.html)`(``"mys_proxmat"``, package ``=`` ``"fastsae"``)`\
\
`# Fit Poisson Spatial Model with population exposure offset`\
`fit_pois`` ``<-`` `[`hb_area`](https://ridsonap.github.io/fastsae/reference/hb_area.md)`(`\
`  formula ``=`` ``y_poisson`` ``~`` ``x1`` ``+`` ``x2``,`\
`  data ``=`` ``sim_area``,`\
`  exposure ``=`` ``"exposure"``,`\
`  family ``=`` ``"poisson"``,`\
`  spatial ``=`` ``"bym2"``,`\
`  W ``=`` ``mys_proxmat``,`\
`  print_result ``=`` ``FALSE`\
`)`\
\
[`summary`](https://rdrr.io/r/base/summary.html)`(``fit_pois``)`\
`#> `\
`#> Variance Components:`\
`#> sigma2_u: 0.058754 `\
`#> rho (spatial): 0.4756 `\
`#> phi (spatial fraction): 0.4756 `\
`#> `\
`#> Coefficients:`\
`#>                    beta   std.error      zvalue      pvalue    ci_lower`\
`#> (Intercept)  1.4536e-01  1.3524e-01  1.0748e+00  2.8246e-01 -1.2125e-01`\
`#> x1           3.8369e-01  4.9803e-02  7.7040e+00  1.3186e-14  2.8563e-01`\
`#> x2          -2.9735e-01  3.3473e-02 -8.8834e+00  6.4852e-19 -3.6363e-01`\
`#>             ci_upper`\
`#> (Intercept)   0.4120`\
`#> x1            0.4820`\
`#> x2           -0.2317`\
`#> `\
`#> Hyperparameters:`\
`#>                                   mean       sd 0.025quant   0.5quant`\
`#> Precision for ..domain_id.. 17.0201186 4.670462 9.50752941 16.4615886`\
`#> Phi for ..domain_id..        0.4755961 0.261137 0.05377607  0.4645501`\
`#>                             0.975quant      mode`\
`#> Precision for ..domain_id..  27.745714 15.439052`\
`#> Phi for ..domain_id..         0.936717  0.229371`\
`#> `\
`#> Goodness of Fit:`\
`#>                                                   DIC `\
`#>                                             317.61486 `\
`#>                                                    pD `\
`#>                                              31.75464 `\
`#>                                                  WAIC `\
`#>                                             309.63790 `\
`#>                                                 pWAIC `\
`#>                                              17.13688 `\
`#> Marginal_LogLik.log marginal-likelihood (integration) `\
`#>                                            -169.74605 `\
`#> `\
`#> HB Summary Statistics:`\
`#>        hb          linear_pred            sd               mse           `\
`#>  Min.   :0.3080   Min.   :-1.1892   Min.   :0.02948   Min.   :0.0008689  `\
`#>  1st Qu.:0.7659   1st Qu.:-0.2776   1st Qu.:0.05787   1st Qu.:0.0033496  `\
`#>  Median :1.3851   Median : 0.3240   Median :0.07672   Median :0.0058864  `\
`#>  Mean   :1.5312   Mean   : 0.2393   Mean   :0.16093   Mean   :0.0665241  `\
`#>  3rd Qu.:2.2717   3rd Qu.: 0.7917   3rd Qu.:0.12182   3rd Qu.:0.0148529  `\
`#>  Max.   :3.5986   Max.   : 1.2802   Max.   :0.75386   Max.   :0.5683065  `\
`#>       rse        `\
`#>  Min.   : 2.790  `\
`#>  1st Qu.: 5.156  `\
`#>  Median : 6.256  `\
`#>  Mean   :10.339  `\
`#>  3rd Qu.: 9.882  `\
`#>  Max.   :29.442`

------------------------------------------------------------------------

## References

- De Nicolò, S., & Gardini, A. (2022). tipsae: Mapping proportions and
  rates using spatio-temporal Beta small area models. *Journal of
  Statistical Software*.
- Janicki, H. (2020). Properties of the Beta-logistic model for small
  area estimation. *Survey Methodology*, 46(1), 89–112.
- Knorr-Held, L. (2000). Bayesian modelling of inhomogeneous spatial and
  temporal variation in rates. *Statistics in Medicine*, 19(17-18),
  2555–2567.
- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
  inference for latent Gaussian models by using integrated nested
  Laplace approximations. *Journal of the Royal Statistical Society:
  Series B*, 71(2), 319–392.
- Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H.
  (2017). Penalising model component complexity: A principled, practical
  approach to constructing priors. *Statistical Science*, 32(1), 1–28.
