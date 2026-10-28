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
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md), a
unified hierarchical Bayesian framework powered by **Integrated Nested
Laplace Approximations (INLA)**. It provides fast deterministic
approximations of the posterior marginals, without the chains, warm-up
and tuning of Markov chain Monte Carlo (MCMC); against long Stan runs
with the same model and priors its posterior means and SDs agreed
closely (see the Engines section of
[`?hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)).

------------------------------------------------------------------------

## Supported Distribution Families

[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
supports 6 probability distributions for the direct response:

| Family | Typical SAE Use Case | Required Sampling Input | Model Scale Parameter |
|:---|:---|:---|:---|
| **`"gaussian"`** | Continuous indicators (e.g., mean income, mean years of schooling) | `vardir` (sampling variance $`D_d`$) | Inverted sampling variance $`\text{scale} = 1/D_d`$ |
| **`"beta"`** | Rates and proportions strictly in $`(0, 1)`$ | `vardir` or `trials` | Precision $`\phi_d = \hat\mu_d(1 - \hat\mu_d)/V_d - 1`$ from a first fit (Janicki, 2020), or $`n_d - 1`$ from `trials` |
| **`"poisson"`** | Counts and rare event rates | `exposure` (expected count / population offset) | Log-link offset $`\log(E_d)`$ |
| **`"nbinomial"`** | Overdispersed counts and clustered events | `exposure`; optional `vardir` | Log-link offset; size $`\hat\mu_d^2/(V_d - \hat\mu_d)`$ from a first fit when `vardir` is given |
| **`"binomial"`** | Aggregated binary proportions | `trials` (domain sample size $`n_d`$, or the effective sample size $`n_d/\mathrm{deff}_d`$ for survey-weighted proportions; not rounded) | Number of trials per domain |
| **`"gamma"`** | Positive, right-skewed indicators (e.g., expenditures) | `vardir` | Precision (shape) $`y_d^2/V_d`$ |

------------------------------------------------------------------------

## Mathematical Formulation

The hierarchical Bayesian model underlying
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
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
\beta_0 \propto 1, \qquad \beta_j \sim N(0, 1000), \quad j = 1, \ldots, p
```

The intercept has INLA’s flat prior and the other coefficients a normal
prior with variance 1000, which is nearly flat on the scale of most
covariates.

### Level 3: Prior on Random Effects

``` math
u \mid \sigma_u \sim N(0, \sigma_u^2 R^{-}), \qquad \sigma_u \sim \text{half-normal}(s)
```

The random effects $`u = (u_1, \ldots, u_m)^\top`$ follow a zero-mean
Gaussian Markov random field whose structure matrix $`R`$ encodes the
chosen spatial or temporal dependence ($`R = I`$ for iid effects; Besag,
BYM2, RW1 and so on otherwise). The scale of the half-normal prior on
their standard deviation is $`s = \mathrm{sd}(y)`$, the standard
deviation of the direct estimates, for Gaussian models and $`s = 1`$ on
the logit or log scale for the other families. The further
hyperparameters (the BYM2 mixing parameter, autocorrelations, dispersion
parameters) and the reasons for these defaults are listed in
[`vignette("hb-priors")`](https://ridsonap.is-a.dev/fastsae/articles/hb-priors.md);
[`inla_engine()`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)
replaces any of them.

### INLA Approximation

Rather than MCMC sampling, INLA approximates the posterior marginals of
the hyperparameters $`\theta`$ and of the latent field
$`x = (\beta, u)`$ (Rue et al., 2009):

``` math
\tilde{\pi}(\theta \mid y) \propto \left. \frac{\pi(x, \theta, y)}{\tilde{\pi}_G(x \mid \theta, y)} \right|_{x = x^*(\theta)}
```

where $`\tilde{\pi}_G`$ is the Gaussian approximation of the full
conditional of $`x`$ and $`x^*(\theta)`$ its mode, and

``` math
\tilde{\pi}(x_d \mid y) = \sum_k \tilde{\pi}(x_d \mid \theta_k, y) \, \tilde{\pi}(\theta_k \mid y) \, \Delta_k
```

where the integrals are evaluated using **Laplace approximations** and
**numerical integration**. Three computational strategies are available
through `hb_area(..., engine = inla_engine(strategy = ...))`, which also
takes the priors and any further
[`INLA::inla()`](https://rdrr.io/pkg/INLA/man/inla.html) arguments:

- **`strategy = "gaussian"`**: Fastest, best for well-behaved problems
- **`strategy = "laplace"`**: Default, good accuracy/speed balance
- **`strategy = "simplified.laplace"`**: Slightly faster Laplace, good
  for large models

The same models can be fitted by MCMC on CPU or GPU with
`engine = jax_engine()` (NumPyro NUTS on JAX; see
[`?jax_engine`](https://ridsonap.is-a.dev/fastsae/reference/hb_engine.md)).
An MCMC fit is only usable when its chains have converged. The fit
reports the rank-normalized R-hat and the effective sample size (ESS) of
every parameter (Vehtari et al., 2021) and warns when R-hat is at least
1.05, the ESS is below 100 per chain, or there are divergent
transitions; the trace and rank plots show which parameter is the
problem:

``` r

fit_jax <- hb_area(y ~ x1 + x2 + x3,
  data = mys, vardir = "vardir",
  engine = jax_engine(chains = 4, samples = 1000)
)
fit_jax$rhat_max # below 1.01 is recommended
fit_jax$ess_min # above 100 per chain
head(fit_jax$mcmc_diagnostics) # R-hat and ESS of each parameter
autoplot(fit_jax, type = "rank") # chains near the dashed line: mixed
autoplot(fit_jax, type = "trace", pars = c("sigma2_u", "hb[5]"))
```

------------------------------------------------------------------------

## Spatial and Spatio-Temporal Structures

### Spatial Priors

1.  **BYM2 (`spatial = "bym2"`)**: The scaled Besag-York-Mollié model
    (Riebler et al., 2016) decomposes area random effects into a spatial
    structured component and an unstructured IID component:
    ``` math
    u_d = \frac{1}{\sqrt{\tau_u}} \left( \sqrt{1 - \phi} \, v_d + \sqrt{\phi} \, u^*_d \right)
    ```
    where $`\phi \in [0, 1]`$ measures the proportion of variance
    explained by spatial structure and has the penalised complexity (PC)
    prior $`P(\phi < 0.5) = 0.5`$ (Simpson et al., 2017).
2.  **Besag ICAR (`spatial = "besag"`)**: Intrinsic Conditional
    Autoregressive spatial model.
3.  **SAR area effect (`spatial = "sar"`)**:
    $`u = (I - \rho W)^{-1}\varepsilon`$, the random effect of the
    spatial Fay-Herriot model of
    [`eblup_sfh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md).
4.  **Spatial Lag Model (`spatial = "slm"`)**:
    $`\eta = (I - \rho W)^{-1}(X\beta + \varepsilon)`$; the covariates
    are filtered too, so this is not the model of
    [`eblup_sfh()`](https://ridsonap.is-a.dev/fastsae/reference/eblup_sfh.md).

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
- **`"type1"` to `"type4"`**: Knorr-Held (2000) interactions
  $`\delta_{dt}`$ added to the main effects,
  $`\eta_{dt} = x_{dt}^\top \beta + s_d + \gamma_t + \delta_{dt}`$:
  - `"type1"`: $`\delta_{dt}`$ independent (unstructured space
    $`\times`$ unstructured time);
  - `"type2"`: an independent `temporal` process for each domain
    (unstructured space $`\times`$ structured time);
  - `"type3"`: an independent ICAR field for each period (structured
    space $`\times`$ unstructured time);
  - `"type4"`: an ICAR field evolving with the `temporal` structure
    (structured space $`\times`$ structured time).

  `"type3"` and `"type4"` need a structured `spatial` effect (`"bym2"`,
  `"bym"`, `"besag"` or `"leroux"`).

------------------------------------------------------------------------

## Case Study: Spatio-Temporal Beta SAE on Emilia-Romagna Poverty

We demonstrate
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
using the Italian poverty dataset `emilia` from the `tipsae` package (38
health districts in Emilia-Romagna across 5 years, $`N = 190`$). The
goal is to estimate the Head Count Ratio (`hcr`), which represents the
proportion of households below the poverty line.

### 1. Load Data and Construct Spatial Adjacency Matrix

``` r

library(fastsae)
library(tipsae)
library(spdep)

data("emilia")
data("emilia_shp")

# Construct binary spatial adjacency matrix W from district polygons
nb <- poly2nb(emilia_shp)
W <- nb2mat(nb, style = "B", zero.policy = TRUE)
rownames(W) <- colnames(W) <- as.character(emilia_shp$NAME_DISTRICT)

head(emilia[, c("id", "year", "hcr", "vars", "x")])
#>                    id year    hcr         vars       x
#> 1 CASALECCHIO DI RENO 2014 0.0404 9.090478e-05 -0.2624
#> 2   CITTA' DI BOLOGNA 2014 0.0825 6.404001e-05 -0.0008
#> 3               IMOLA 2014 0.1033 3.120275e-04 -0.0522
#> 4         PIANURA EST 2014 0.0633 1.025764e-04 -0.4007
#> 5       PIANURA OVEST 2014 0.0625 1.562500e-04 -0.2277
#> 6      PORRETTA TERME 2014 0.1276 6.643609e-04 -0.4434
```

### 2. Fit Spatio-Temporal Beta Model

We model the poverty proportion using a Beta likelihood with a Besag
ICAR spatial random effect, a domain-specific RW(1) temporal dynamic,
and known sampling variances `vars`:

``` r

fit_beta_st <- hb_area(
  formula = hcr ~ x,
  data = emilia,
  domain = "id",
  time = "year",
  vardir = "vars",
  family = "beta",
  spatial = "besag",
  temporal = "rw1",
  st_interaction = "domain-specific",
  W = W,
  print_result = FALSE
)

summary(fit_beta_st)
#> 
#> Variance Components:
#> sigma2_u: 0.07407 
#> sigma2_t (temporal): 0.008528 
#> 
#> Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept) -2.2696e+00  1.5370e-02 -1.4766e+02  0.0000e+00 -2.2998e+00
#> x            4.1691e-01  6.6161e-02  6.3015e+00  2.9478e-10  2.8685e-01
#>             ci_upper
#> (Intercept)  -2.2395
#> x             0.5465
#> 
#> Hyperparameters:
#>                                  mean        sd 0.025quant  0.5quant 0.975quant
#> Precision for ..domain_id..  15.25291  5.508698    7.22469  14.33838   28.60708
#> Precision for ..time_id..   136.59835 57.259684   60.48608 125.04665  281.22075
#>                                  mode
#> Precision for ..domain_id..  12.67088
#> Precision for ..time_id..   104.92244
#> 
#> Goodness of Fit:
#>                                                   DIC 
#>                                            -971.94545 
#>                                                    pD 
#>                                              60.69152 
#>                                                  WAIC 
#>                                            -979.78139 
#>                                                 pWAIC 
#>                                              42.46189 
#> Marginal_LogLik.log marginal-likelihood (integration) 
#>                                             432.98607 
#> 
#> HB Summary Statistics:
#>        hb           linear_pred           sd                mse           
#>  Min.   :0.04905   Min.   :-2.970   Min.   :0.004907   Min.   :2.408e-05  
#>  1st Qu.:0.08142   1st Qu.:-2.428   1st Qu.:0.007907   1st Qu.:6.253e-05  
#>  Median :0.09444   Median :-2.266   Median :0.009723   Median :9.454e-05  
#>  Mean   :0.09673   Mean   :-2.268   Mean   :0.009911   Mean   :1.057e-04  
#>  3rd Qu.:0.10925   3rd Qu.:-2.103   3rd Qu.:0.011574   3rd Qu.:1.339e-04  
#>  Max.   :0.15774   Max.   :-1.678   Max.   :0.018014   Max.   :3.245e-04  
#>       rse        
#>  Min.   : 6.545  
#>  1st Qu.: 9.167  
#>  Median :10.224  
#>  Mean   :10.312  
#>  3rd Qu.:11.594  
#>  Max.   :14.895
```

### 3. Inspect Posterior Predictions and Uncertainty

The output object contains full domain-period predictions, posterior
standard errors, and 95% Credible Intervals:

``` r

head(fit_beta_st$df_hb[, c("domain", "time", "y", "hb", "sd", "rse", "ci_lower", "ci_upper")])
#>                domain time      y         hb          sd       rse   ci_lower
#> 1 CASALECCHIO DI RENO 2014 0.0404 0.05102821 0.005596864 10.968176 0.04062425
#> 2   CITTA' DI BOLOGNA 2014 0.0825 0.08182497 0.006281355  7.676576 0.07003561
#> 3               IMOLA 2014 0.1033 0.09310380 0.010890483 11.697140 0.07366716
#> 4         PIANURA EST 2014 0.0633 0.06365559 0.006408113 10.066851 0.05185443
#> 5       PIANURA OVEST 2014 0.0625 0.06937248 0.007927322 11.427186 0.05493779
#> 6      PORRETTA TERME 2014 0.1276 0.09446590 0.014070634 14.894934 0.06968697
#>     ci_upper
#> 1 0.06256393
#> 2 0.09467502
#> 3 0.11635703
#> 4 0.07700998
#> 5 0.08602893
#> 6 0.12479148
```

------------------------------------------------------------------------

## Comparison and Validation Against `tipsae` (Stan MCMC)

In `tipsae::fit_sae(..., spatial_error = TRUE, temporal_error = TRUE)`,
the model is estimated by Hamiltonian Monte Carlo (HMC) in Stan.
[`fastsae::hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
fits a model with the same structure by INLA; the priors of the two
packages differ.

| Characteristic | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Stan MCMC) | [`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (INLA) |
|:---|:---|:---|
| **Model** | Beta, Besag ICAR + domain-specific RW(1) | `spatial = "besag"`, `temporal = "rw1"`, `st_interaction = "domain-specific"` |
| **Estimation** | MCMC (NUTS) | Integrated nested Laplace approximation (INLA’s default `"compact"` mode for Beta) |
| **Execution time** | 22 min (4 chains x 10,000 iterations, largest $`\hat R`$ 1.007); 4 chains x 2,000 iterations (120 s) have not converged | 1.1 s |
| **Posterior means vs converged `tipsae`** | Reference | $`r = 0.999`$, mean absolute difference 0.0017 |

Measured on 2026-10-16 (`tipsae` reference refitted on 2026-10-17) on
one machine (Apple silicon); times depend on the hardware. The accuracy
of [`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
was also checked against Stan with identical models and priors; see the
Engines section of
[`?hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md).

### Interactive Scalability Explorer ($`n = 30`$ to $`n = 1,000`$)

Use the interactive controls below to compare runtime, memory footprint,
and speedup factors between
[`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
and [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html)
across both **Standard Beta** and **Spatial Beta (Besag ICAR)** models:

Standard Beta

Spatial Beta (ICAR)

Spatio-Temporal Beta

Execution Time

Memory Allocated

Throughput (iter/s)

Execution time (median, seconds) — Standard Beta · log scale

Speedup Factor & Efficiency Multiplier (fastsae vs tipsae) — Standard
Beta

#### Agreement with Stan MCMC (Beta model, $`n = 1,000`$)

To check that the speed of
[`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
does not come at the cost of its estimates, both packages were fitted
with their defaults to the $`n = 1,000`$ data set of the Beta benchmark
(independent area effects, sampling variances from the direct estimates;
`benchmarks/bench_beta_concordance.R`). The two packages do not fit
identical models (their priors differ, and fastsae fixes the precision
at a first-stage estimate while tipsae makes it a function of the model
mean), so the comparison is between what a user obtains from each
package:

| Metric | [`fastsae::hb_area`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md) (INLA) | [`tipsae::fit_sae`](https://rdrr.io/pkg/tipsae/man/fit_sae.html) (Stan) |
|:---|:--:|:--:|
| RMSE against the true area means | 0.0678 | 0.0684 |
| Coverage of the 95% intervals | 93.3% | 91.1% |
| Mean width of the 95% intervals | 0.252 | 0.259 |
| Time | 4.3 s | 184 s |

The posterior means of the two packages differed by 0.0096 on average
(absolute difference), and the posterior SDs of fastsae were 0.99 times
those of tipsae in the median area (90% of the areas between 0.88 and
1.17). The largest $`\hat{R}`$ of the tipsae fit was 1.008.

------------------------------------------------------------------------

## Cross-Sectional Count Data: Poisson & Negative Binomial

For count data such as disease incidence or crime events,
[`hb_area()`](https://ridsonap.is-a.dev/fastsae/reference/hb_area.md)
seamlessly supports Poisson and Negative Binomial distributions with
population offsets:

``` r

# Simulate 42-domain count data with known spatial structure
data("sim_area", package = "fastsae")

# The BYM2 model uses a neighbourhood graph: the 4 nearest centroids of each domain
W_knn <- create_weights(cbind(sim_area$x_coord, sim_area$y_coord), method = "knn", k = 4)

# Fit Poisson Spatial Model with population exposure offset
fit_pois <- hb_area(
  formula = y_poisson ~ x1 + x2,
  data = sim_area,
  exposure = "exposure",
  family = "poisson",
  spatial = "bym2",
  W = W_knn,
  print_result = FALSE
)

summary(fit_pois)
#> 
#> Variance Components:
#> sigma2_u: 0.069376 
#> rho (spatial): 0.1595 
#> phi (spatial fraction): 0.1595 
#> 
#> Coefficients:
#>                    beta   std.error      zvalue      pvalue    ci_lower
#> (Intercept)  1.4161e-01  1.3846e-01  1.0228e+00  3.0643e-01 -1.3129e-01
#> x1           3.8646e-01  5.0100e-02  7.7137e+00  1.2219e-14  2.8774e-01
#> x2          -2.9875e-01  3.3665e-02 -8.8741e+00  7.0474e-19 -3.6547e-01
#>             ci_upper
#> (Intercept)   0.4149
#> x1            0.4853
#> x2           -0.2327
#> 
#> Hyperparameters:
#>                                   mean        sd  0.025quant   0.5quant
#> Precision for ..domain_id.. 15.6290050 4.4994512 8.488591703 15.0604006
#> Phi for ..domain_id..        0.1595282 0.1562824 0.007218553  0.1055143
#>                             0.975quant        mode
#> Precision for ..domain_id.. 26.0450791 14.01277938
#> Phi for ..domain_id..        0.5934515  0.01699368
#> 
#> Goodness of Fit:
#>                                                   DIC 
#>                                             317.74616 
#>                                                    pD 
#>                                              31.88699 
#>                                                  WAIC 
#>                                             309.67486 
#>                                                 pWAIC 
#>                                              17.15907 
#> Marginal_LogLik.log marginal-likelihood (integration) 
#>                                            -185.87986 
#> 
#> HB Summary Statistics:
#>        hb          linear_pred            sd               mse           
#>  Min.   :0.3061   Min.   :-1.1959   Min.   :0.02962   Min.   :0.0008776  
#>  1st Qu.:0.7648   1st Qu.:-0.2791   1st Qu.:0.05813   1st Qu.:0.0033800  
#>  Median :1.3850   Median : 0.3240   Median :0.07685   Median :0.0059062  
#>  Mean   :1.5294   Mean   : 0.2381   Mean   :0.16304   Mean   :0.0690514  
#>  3rd Qu.:2.3016   3rd Qu.: 0.8023   3rd Qu.:0.12189   3rd Qu.:0.0148681  
#>  Max.   :3.6016   Max.   : 1.2810   Max.   :0.76548   Max.   :0.5859662  
#>       rse        
#>  Min.   : 2.791  
#>  1st Qu.: 5.164  
#>  Median : 6.265  
#>  Mean   :10.478  
#>  3rd Qu.: 9.915  
#>  Max.   :30.137
```

------------------------------------------------------------------------

## References

- De Nicolò, S., & Gardini, A. (2024). The R package tipsae: Tools for
  mapping proportions and indicators on the unit interval. *Journal of
  Statistical Software*, 108(1), 1–36.
- Janicki, R. (2020). Properties of the beta regression model for small
  area estimation of proportions and application to estimation of
  poverty rates. *Communications in Statistics - Theory and Methods*,
  49(9), 2264–2284.
- Knorr-Held, L. (2000). Bayesian modelling of inhomogeneous spatial and
  temporal variation in rates. *Statistics in Medicine*, 19(17-18),
  2555–2567.
- Riebler, A., Sørbye, S. H., Simpson, D., & Rue, H. (2016). An
  intuitive Bayesian spatial model for disease mapping that accounts for
  scaling. *Statistical Methods in Medical Research*, 25(4), 1145–1165.
- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
  inference for latent Gaussian models by using integrated nested
  Laplace approximations. *Journal of the Royal Statistical Society:
  Series B*, 71(2), 319–392.
- Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H.
  (2017). Penalising model component complexity: A principled, practical
  approach to constructing priors. *Statistical Science*, 32(1), 1–28.
- Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
  (2021). Rank-normalization, folding, and localization: An improved
  R-hat for assessing convergence of MCMC. *Bayesian Analysis*, 16(2),
  667–718.
