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
`#> x            4.0612e-01  6.4372e-02  6.3090e+00  2.8084e-10  2.7947e-01`\
`#>             ci_upper`\
`#> (Intercept)  -2.2218`\
`#> x             0.5321`\
`#> `\
`#> Hyperparameters:`\
`#>                                  mean        sd 0.025quant  0.5quant 0.975quant`\
`#> Precision for ..domain_id..  18.22307  6.626532   8.634569  17.10219    34.3505`\
`#> Precision for ..time_id..   172.31221 78.434145  72.154753 155.44927   373.0483`\
`#>                                  mode`\
`#> Precision for ..domain_id..  15.06357`\
`#> Precision for ..time_id..   127.08259`\
`#> `\
`#> Goodness of Fit:`\
`#>                                                   DIC `\
`#>                                            -976.34166 `\
`#>                                                    pD `\
`#>                                              56.02386 `\
`#>                                                  WAIC `\
`#>                                            -984.78987 `\
`#>                                                 pWAIC `\
`#>                                              38.70815 `\
`#> Marginal_LogLik.log marginal-likelihood (integration) `\
`#>                                             439.32542 `\
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
`#> 1 CASALECCHIO DI RENO 2014 0.0404 0.05554640 0.006127610 11.031517 0.04408725`\
`#> 2   CITTA' DI BOLOGNA 2014 0.0825 0.08225895 0.006131567  7.453981 0.07068300`\
`#> 3               IMOLA 2014 0.1033 0.09343262 0.010214743 10.932738 0.07511425`\
`#> 4         PIANURA EST 2014 0.0633 0.06509827 0.006294797  9.669684 0.05342022`\
`#> 5       PIANURA OVEST 2014 0.0625 0.07112654 0.007870670 11.065728 0.05671270`\
`#> 6      PORRETTA TERME 2014 0.1276 0.09608738 0.012916482 13.442434 0.07314336`\
`#>     ci_upper`\
`#> 1 0.06810568`\
`#> 2 0.09473235`\
`#> 3 0.11517443`\
`#> 4 0.07814090`\
`#> 5 0.08759940`\
`#> 6 0.12376161`

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
`#> x1           3.8369e-01  4.9804e-02  7.7040e+00  1.3187e-14  2.8563e-01`\
`#> x2          -2.9735e-01  3.3473e-02 -8.8834e+00  6.4859e-19 -3.6363e-01`\
`#>             ci_upper`\
`#> (Intercept)   0.4120`\
`#> x1            0.4820`\
`#> x2           -0.2317`\
`#> `\
`#> Hyperparameters:`\
`#>                                   mean        sd 0.025quant   0.5quant`\
`#> Precision for ..domain_id.. 17.0200956 4.6704720 9.50753133 16.4615515`\
`#> Phi for ..domain_id..        0.4755909 0.2611461 0.05377634  0.4645349`\
`#>                             0.975quant       mode`\
`#> Precision for ..domain_id.. 27.7457544 15.4389788`\
`#> Phi for ..domain_id..        0.9367355  0.2292976`\
`#> `\
`#> Goodness of Fit:`\
`#>                                                   DIC `\
`#>                                             317.61486 `\
`#>                                                    pD `\
`#>                                              31.75465 `\
`#>                                                  WAIC `\
`#>                                             309.63789 `\
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
`#>  Mean   :1.5312   Mean   : 0.2393   Mean   :0.16093   Mean   :0.0665243  `\
`#>  3rd Qu.:2.2717   3rd Qu.: 0.7917   3rd Qu.:0.12182   3rd Qu.:0.0148529  `\
`#>  Max.   :3.5986   Max.   : 1.2802   Max.   :0.75386   Max.   :0.5683086  `\
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
