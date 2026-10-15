# Catatan teknis: lapisan Bayes/INLA (`hb_area`, `hb_unit`, `hb_twofold`, `inla_utils`)

Sumber: kode yang benar-benar dibaca di `/Volumes/work/_MainR/fastsae-cov/R/`. Semua nomor baris = baris file asli.
Klaim semantik INLA diverifikasi dari dok INLA terpasang
(`/Library/Frameworks/R.framework/Versions/4.6/Resources/library/INLA/documentation/*.pdf`, diekstrak `gs -sDEVICE=txtwrite`)
dan dari eksperimen kecil `INLA::inla()` (ditandai **[E]**). Tidak ada `\doi{}` di roxygen keempat file (grep `doi` di `R/` = 0 hit).

## 0. Fakta lintas-fungsi (terverifikasi)

| Fakta | Bukti |
|---|---|
| `control.fixed` tidak pernah di-set paket → default INLA berlaku: `mean=0, mean.intercept=0, prec=0.001, prec.intercept=0 (flat), expand.factor.strategy="model.matrix", compute=TRUE, correlation.matrix=FALSE, cdf=NULL, quantiles=NULL` | `INLA::inla.set.control.fixed.default()` **[E]**; tidak ada `control.fixed` di 4 file |
| Gaussian INLA: `Var(y)=1/(s·τ)`, `s` = argumen `scale` | `documentation/likelihood/gaussian.pdf`; **[E]** `scale=4`, `prec` fixed `initial=0` → sd(μ̂)=0.25=1/√(4·4) |
| Hyper likelihood yang `fixed=TRUE` hilang dari `summary.hyperpar` | **[E]** `family=gaussian` + `scale` + `prec` fixed → baris hanya `"Precision for id"` |
| `summary.fitted.values` = **inverse-link(η)** (bukan E[y] bila ada offset/Ntrials) | **[E]** binomial `y=5,6,4,7`, `Ntrials=10` → fitted 0.55 (=p), bukan 5.5; poisson `E=100` → fitted 0.055 (=rate), bukan 5.5; nbinomial `E` juga menerima `E` (fitted 0.055 vs 5.53 tanpa E) |
| `pc.prec` param `c(u,α)`: **Prob(σ > u) = α**, σ=1/√τ | `documentation/prior/pc.prec.pdf` |
| BYM2 phi `pc` param `c(u,α)`: **Prob(ϕ ≤ u) = α**; default INLA BYM2 = `pc.prec(1,0.01)` + `pc(0.5,0.5)`, `constr=TRUE` | `documentation/latent/bym2.pdf`; `inla.models()` **[E]** |
| `generic1`: **Q = τ(I − (β/λ_max)·C)**, β∈[0,1), θ2 = logit(β), prior default `gaussian(0,0.1)`, `constr=FALSE`, hanya butuh `Cmatrix` (bukan `Rmatrix`/`graph`) | `documentation/latent/generic1.pdf` |
| `slm`: x=(I−ρW)⁻¹(Xβ+ε); ρ* =(ρ−ρmin)/(ρmax−ρmin), θ2=logit(ρ*); β ada di **ekor** `summary.random[[1]]`; re-scale ρ: `rho.min + x*(rho.max-rho.min)`; prior default `normal(0,10)` | `documentation/latent/slm.pdf` |
| `constr` default INLA: bym2/bym/besag/rw1/rw2 **TRUE**; iid/ar1/generic1/slm **FALSE** | `inla.models()` **[E]** (pakai sendiri tidak pernah mengoper `constr`) |
| Skala indeks `f(index, group)`: posisi = `(group−1)·n_index + index` | **[E]** uji `f(time_id, group=domain_id, control.group=iid)` → 4 baris pertama = domain 1 dan monoton sesuai waktu ⇒ **`inla_utils.R:239` benar** |
| **INLA melarang dua `f()` memakai covariate/index yang sama** | **[E]** error: "The covariate `..time_id..' are used in the following f()-terms: 1, 2 / Only one is allowed…" |
| `inla()` tidak punya argumen `...` (42 formal) → daftar argumen ber-nama ganda bikin error saat `do.call` | **[E]** `formal argument "control.compute" matched by multiple actual arguments` |

---

## 1. `R/hb_area.R` (813 baris)

### 1.1 Signature + export
`export(hb_area)` (`NAMESPACE:10`). Definisi `hb_area` **R/hb_area.R:148-314**:

```r
hb_area <- function(formula, data, domain = NULL, time = NULL,
  family = c("gaussian","binomial","poisson","nbinomial","beta","gamma"),
  spatial = c("none","bym2","bym","besag","generic1","slm"),
  temporal = c("none","rw1","rw2","ar1","iid"),
  st_interaction = c("none","domain-specific","separable","type1","type2","type3","type4"),
  W = NULL, vardir = NULL, trials = NULL, exposure = NULL,
  method = c("inla","laplace"),
  strategy = c("laplace","simplified.laplace","gaussian"),
  link = NULL, scale_model = TRUE,
  prior_prec  = list(prior = "pc.prec", param = c(1, 0.01)),
  prior_phi   = list(prior = "pc",     param = c(0.5, 0.5)),
  prior_rho = NULL, prior_prec_time = NULL, prior_rho_time = NULL,
  print_result = TRUE, ...)
```
- `match.arg(tolower(.), choices=...)`: 174-179 (family, spatial, temporal, st_interaction, method, strategy).
- Return: class `c("fastsae_hb_area","fastsae")` (dibuat `inla_utils.R:394`; jalur lme4 juga `hb_area.R:811`), **tanpa** elemen `fastsae_hb` (berbeda dari hb_unit/twofold).

### 1.2 Routing (pseudocode bernomor)
1. `call_matched <- match.call()` (173); validasi `data` data.frame (181-183).
2. ID: `domain` → `.get_variable` (`R/utils.R:19-38`: string kolom | formula satu sisi | vektor berpanjang `nrow(data)`); `domain=NULL & time=NULL` → `domain_vec=1..n` (188-193); `time` wajib ada `domain` (192); `temporal!="none"` wajib `time` (202-204) dan ≥2 periode (209-211).
3. `mf <- model.frame(formula, data, na.action=na.pass)`; `y <- model.response` (214-215); `vardir/trials/exposure` diekstrak via `.get_variable` (218-220).
4. Binomial (224-240): `trials` wajib (226); jika y semua di [0,1] dan ada pecahan → `y <- as.integer(round(y*trials_vec))` dan kolom respons di `data` ditimpa (231-234); selain itu `y <- round(y)` (236-239).
5. Gamma: `y > 0` wajib (242-247). `spatial!="none"` → `W_obj <- .convert_spatial_weights(W, n_domains=n_unique_domains, ...)` (250-253).
6. **Jalur lme4** hanya bila `method=="laplace" && spatial=="none" && temporal=="none" && family %in% c("binomial","poisson")` → `.fit_glmm_laplace` (258-267). Selain itu masuk INLA; bila `method=="laplace"` → `strategy` dipaksa `"laplace"` (269-276).
7. `.fit_inla_area(...)` (280-304) → `out$W <- W` (307) → `print(out)` bila `print_result` (309-311).

### 1.3 `.fit_inla_area` (318-704): rumus INLA
Index: `..domain_id..` (363-366), `..obs_id..=1..n_obs` (367), `..time_id..` (372-377). Pemetaan family → `inla_family` = nama sama (380-388). `fixed_str` = `paste(term.labels," + ")` atau `"1"` (391-393).

**Prior waktu (399-405), kode persis:**
```r
if (is.null(prior_prec_time)) {
  prior_prec_time <- list(prior = "pc.prec", param = c(1, 0.01))
}
hyper_time <- list(prec = prior_prec_time)
if (temporal == "ar1" && !is.null(prior_rho_time)) {
  hyper_time$rho <- prior_rho_time
}
```

**Bagian spasial** — string persis (semua di-`paste0` lalu `as.formula`):

| spatial | string `f(...)` | baris |
|---|---|---|
| none | `f(..domain_id.., model = 'iid', hyper = list(prec = prior_prec))` | 411-413 (hanya bila `temporal=="none"` atau `st_interaction %in% none/domain-specific/type1/type2`) |
| bym2 | `f(..domain_id.., model = 'bym2', graph = W_obj$graph, scale.model = <scale_model>, hyper = list(prec = prior_prec, phi = prior_phi))` | 414-421 (skip bila `st_interaction=="separable"`) |
| bym | `f(..domain_id.., model = 'bym', graph = W_obj$graph, scale.model = <>, hyper = list(prec.unstruct = prior_prec, prec.spatial = prior_prec))` | 422-429 |
| besag | `f(..domain_id.., model = 'besag', graph = W_obj$graph, scale.model = <>, hyper = list(prec = prior_prec))` | 430-437 |
| generic1 | `f(..domain_id.., model = 'generic1', Cmatrix = C_mat, hyper = hyper_generic1)`; `R_mat <- diag(rowSums(adj)) - adj; C_mat <- diag(nrow(adj)) - R_mat` (439-441) | 438-450 |
| slm | `f(..domain_id.., model = 'slm', args.slm = args_slm, hyper = hyper_slm)` — **tanpa guard separable** | 451-488 |

- `hyper_generic1` (442-447): `list(prec = prior_prec)`; `beta <- prior_rho` bila ada; else `beta <- prior_phi` bila `prior_phi$prior != "pc"`.
- `args_slm` (471-477): `rho.min/rho.max` dari akar eigen `W_std` (453-463, `W_std = W_mat * (1/rowSums)` = row-standardisasi; clamp `[-0.999, 0.999]`), `W = Matrix(W_std,sparse)`, `X = model.matrix` (465-467), `Q.beta = Diagonal(p_cov, 1e-4)` (469).
- `hyper_slm` (478-485): `list(prec = prior_prec)`; `rho <- prior_rho`; else `prior_phi` bila `prior != "pc"`; **else default `list(prior = "logitbeta", param = c(1, 1))`** (484).

**Bagian temporal/spatio-temporal** (491-567; `scale_model_time` hanya untuk `rw1`/`rw2`, 492):

| `st_interaction` | string tambahan | baris | padanan Knorr-Held |
|---|---|---|---|
| none | `f(..time_id.., model = '<temporal>'<,scale.model>, hyper = hyper_time)` | 494-498 | efek utama aditif |
| domain-specific | `f(..time_id.., model = '<temporal>', group = ..domain_id.., control.group = list(model = 'iid')<,scale.model>, hyper = hyper_time)` | 499-504 | (tipsae): RW/AR per domain, var bersama |
| separable | grup spasial: `f(..domain_id.., model='<bym2\|besag\|bym\|generic1>', graph=…, group = ..time_id.., control.group = list(model = '<temporal>'), scale.model=<>, hyper=…)`; fallback `model='iid'` (527-531) | 505-532 | Q_space ⊗ Q_time |
| type1 | efek waktu + `f(..obs_id.., model = 'iid', hyper = list(prec = prior_prec))` | 533-537 | **KH I** (ruang tak-berstruktur × waktu tak-berstruktur) |
| type2 | efek waktu + `f(..time_id.., model='<temporal>', group=..domain_id.., control.group=list(model='iid')<…>, hyper=hyper_time)` | 538-545 | **KH II** |
| type3 | efek waktu + `f(..domain_id.., model='<spatial>', graph=…, group=..time_id.., control.group=list(model='iid'), scale.model=<>, hyper=list(prec=prior_prec))` bila `spatial!="none"` | 546-555 | **KH III** |
| type4 | efek waktu + `f(..domain_id.., model='<spatial>', graph=…, group=..time_id.., control.group=list(model='<<temporal>>'), scale.model=<>, hyper=list(prec=prior_prec))` bila `spatial!="none"` | 556-566 | **KH IV** |

Perakitan (569-575): `spatial=="slm"` → `paste(y, "~ -1 +", rand_terms)`; selain itu `paste(y, "~", fixed_str, "+", rand_terms)`.

### 1.4 Kontrol INLA (578-609)
```r
control_compute  <- list(dic = TRUE, waic = TRUE, cpo = TRUE, mlik = TRUE, config = TRUE)  # 578-584
control_predictor <- list(compute = TRUE, link = 1)                                         # 587
control_inla     <- list(strategy = strategy)                                               # 590
inla_args <- list(formula=…, family=…, data=inla_data, control.compute=control_compute,
                  control.predictor=control_predictor, control.inla=control_inla, ...)       # 593-601
```
`num.threads`: hanya dipaksa `"1:1"` bila env CRAN (`_R_CHECK_LIMIT_CORES_`/`_R_CHECK_PACKAGE_NAME_`) dan user belum mengisi (604-609); opsi global `num.threads="1:1"` + `mc.cores=1` juga hanya saat CRAN (344-357). Di luar CRAN: **default `inla.getOption("num.threads")`**. `control.fixed` tidak di-set.

### 1.5 Likelihood per family (vardir/trials/exposure)
- **gaussian + `vardir`** (612-620): validasi `vardir>0` (613-615); `inla_args$scale <- 1 / vardir`; `inla_args$control.family <- list(hyper = list(prec = list(initial = 0, fixed = TRUE)))` ⇒ Var = 1/(s·τ) = vardir **[E]**.
- **gaussian tanpa `vardir`**: tidak ada `scale`/`control.family` → τ residual diestimasi (default `loggamma(1,5e-05)`).
- **beta** (623-642): dengan `vardir` → `phi_dir <- (y*(1-y)/vardir) - 1; phi_dir[phi_dir<1] <- 1; scale <- phi_dir`, `control.family = list(hyper = list(theta = list(initial = 0, fixed = TRUE)))`; dengan `trials` → `scale <- pmax(trials-1,1)` (634-641). Sesuai `likelihood/beta.pdf`: φ = s·exp(θ), θ fixed 0 ⇒ φ = s.
- **gamma + `vardir`** (645-657): `s_gamma = y^2/vardir` (default 1 utk NA/nonpositif) sebagai `scale`, `control.family = list(hyper = list(theta = list(initial = 0, fixed = TRUE)))` ⇒ Var ≈ µ²/s.
- **binomial** (660-662): `inla_args$Ntrials <- trials` (link default logit; `trials` wajib, 226).
- **poisson/nbinomial** (665-667): `inla_args$E <- exposure` (link default log; `E` diterima keduanya **[E]**).
- **nbinomial** tanpa exposure: tanpa `E`.
- Link: `link` argumen **tidak pernah dipakai** (lihat §5-D6); `control.predictor$link = 1` = link default family.

### 1.6 EBP & ketidakpastian
- `fit <- tryCatch(do.call(INLA::inla, inla_args), …)` (670-678) → `.extract_inla_results(...)` (681-701).
- Titik prediksi (`inla_utils.R:274-284`): `hb_est = summary.fitted.values$mean` (inverse-link, skala respons; utk area-level = rata-rata area), `hb_sd = $sd`, `mse = sd^2`, `rse = 100·sd/|est|`, CI = `0.025quant/0.975quant`, `linear_pred = summary.linear.predictor$mean`. **Tanpa transformasi tambahan ke skala rata-rata area** (memang sudah skala area). Kolom tambahan: `vardir`/`precision` (beta/gamma) 321-331, `trials`/`estimated_total = hb*trials` 334-341, `exposure`/`rate`/`estimated_count = hb*exposure` 344-348.
- **Domain unsampled**: `hb_area` tidak pernah membuat baris prediksi. Baris yang hilang dari `data` tidak dievaluasi; baris dengan `y = NA` tetap dihitung (`na.pass` 214, `seq_len(n_obs)` di 275-276). Tidak ada flag `unsampled_domains` (ada di hb_unit).
- **Ketidakpastian**: murni analitik `mse = sd²` (280). **Tidak ada `inla.posterior.sample` di `hb_area`**; agregasi area tidak dilakukan.

### 1.7 `.fit_glmm_laplace` (708-813)
1. `formula_lmer = y ~ <fixed> + (1 | ..domain_id..)` (724).
2. gaussian → `lme4::lmer(..., REML=FALSE)` (730); binomial+trials → `glmer(cbind(successes, failures) ~ …, family=binomial(logit), nAGQ=1)` (732-736); binomial tanpa trials → `glmer(y ~ …)` (738); poisson+exposure → `offset(log(exposure))`, `family=poisson(log)`, `nAGQ=1` (740-744); poisson tanpa exposure (746).
3. `estcoef`: β, se, z, p=2·pnorm(|z|,lower=F), CI ±1.96·se (752-761).
4. `sigma2_u = VarCorr(fit)[[1]]` (764-765); prediksi `predict(type="response"/"link")` (768-769); `ranef` (770).
5. `df_hb`: `sd/mse/rse/ci_lower/ci_upper = NA_real_` (777-781); `goodness = c(AIC, BIC, LogLik)` (786-790); `hyperpar = data.frame(Parameter="sigma2_u", Estimate=…)` (797); `model = "HB-<FAMILY> (Laplace GLMM)"` (804).

### 1.8 `@references` (110-121, persis)
```
\enumerate{
  \item Rao, J. N. K., and Molina, I. (2015). \emph{Small Area Estimation}. John Wiley & Sons.
  \item Riebler, A., Sørbye, S. H., Simpson, D., and Rue, H. (2016). An intuitive Bayesian spatial model
    for disease mapping that accounts for scaling. \emph{Statistical Methods in Medical Research}, 25(4), 1145-1165.
  \item Marhuenda, Y., Molina, I., and Morales, D. (2013). Small area estimation with spatio-temporal Fay-Herriot models.
    \emph{Computational Statistics & Data Analysis}, 58, 308-325.
  \item De Nicolò, S., and Gardini, A. (2024). The R Package tipsae: Tools for Mapping Proportions and Indicators on the Unit Interval.
    \emph{Journal of Statistical Software}, 108(1), 1-36.
  \item Rue, H., Martino, S., and Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian
    models by using integrated nested Laplace approximations. \emph{Journal of the Royal Statistical Society: Series B}, 71(2), 319-392.
}
```

---

## 2. `R/hb_unit.R` (481 baris)

### 2.1 Signature (124-140) + `match.arg` (144-146)
```r
hb_unit <- function(formula, unit_data, Xpop, domain_var, popsize_var = NULL,
  family = c("gaussian","binomial","poisson"),
  spatial = c("none","bym2","besag"),
  W = NULL, popnmean_xpop = NULL,
  strategy = c("laplace","simplified.laplace"),
  scale_model = TRUE,
  prior_prec = list(prior = "loggamma", param = c(0.01, 0.01)),
  prior_phi  = list(prior = "pc", param = c(0.5, 0.5)),
  print_result = TRUE, ...)
```
Return: `c("fastsae_hb_unit","fastsae_hb","fastsae")` (471), `invisible(result)` (480).

### 2.2 Pseudocode bernomor
1. `.check_inla_installed()` (142); `match.arg` family/spatial/strategy (144-146); validasi formula/data/Xpop/domain_var/popsize_var (148-169).
2. Sampel: `model.frame(na.action=na.omit)` (172), `dom_s` diselaraskan dgn NA (174-176), `y_s` (178), wajib ≥1 prediktor (182-184).
3. Populasi: `all_domains = unique(Xpop[[domain_var]])` (187); sampel tak ada di Xpop → abort (191-196); `unsampled_domains = setdiff(all_domains, unique(dom_s))` (199).
4. Bila `nrow(Xpop) > n_domains` → agregasi `mean` per domain dgn `aggregate(cbind(<terms>) ~ <domain_var>)` (202-208), lalu `match(all_domains, …)` (214).
5. `meanxpop`: `model.matrix(reformulate(term_labels), Xpop_domain)` bila `popnmean_xpop=NULL` (217-219); bila diberi & kolomnya = p−1 → tambah kolom `(Intercept)` (221-224).
6. W → `.convert_spatial_weights` → `adj_graph` (228-232).
7. Data INLA: baris sampel (`df_s`, y nyata) + **1 baris prediksi per domain dengan `y = NA`** (244-249) → `combined_df`; `..dom_id..` = 1..n_domains (236-247).
8. Formula: `y ~ <fixed_str> + rand_str`, `rand_str` (253-263):
   - none: `f(..dom_id.., model = 'iid', hyper = list(prec = prior_prec))`
   - bym2: `f(..dom_id.., model = 'bym2', graph = adj_graph, scale.model = <>, hyper = list(prec = prior_prec, phi = prior_phi))`
   - besag: `f(..dom_id.., model = 'besag', graph = adj_graph, scale.model = <>, hyper = list(prec = prior_prec))`
9. Thread: `inla.setOption("num.threads","1:1")` + `on.exit` restore (268-272); args (274-282): `control.predictor=list(compute=TRUE, link=1)`, `control.inla=list(strategy=strategy)`, `control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE, mlik=TRUE)` (**tanpa `config`**), `num.threads = 1`; `...` **menimpa** nama yang sudah ada: `inla_args[names(extra_args)] <- extra_args` (283-284).
10. `do.call(INLA::inla, inla_args)` (287-295).

### 2.3 Likelihood
`family` dikirim apa adanya ke `inla()` (277) → link default: gaussian identitas, binomial logit, poisson log.
- **Tidak ada `Ntrials`** → binomial = Bernoulli per unit; **tidak ada `E`** → poisson tanpa exposure/offset.
- **Tidak ada `scale`/`control.family` sama sekuali** → tidak ada jalan untuk `vardir` diketahui di level unit.
- Gaussian: `sigma2_e` dibaca dari `summary.hyperpar` baris `"Gaussian"` (386-390) lalu dipakai utk `vardir = σ̂²_e/n_d` (401-409); binomial → `p̄(1−p̄)/n`; poisson → `ȳ/n` (403-407).

### 2.4 EBP, domain unsampled, ketidakpastian
- `pred_indices <- (n_sample + 1):nrow(combined_df)` (298); `hb = summary.fitted.values$mean`, `sd`, `0.025quant`, `0.975quant`; `linear_pred = summary.linear.predictor$mean` (299-306) → **superpopulation** tanpa transformasi skala-populasi.
- Fin populasi bila `popsize_var` (315-343): `f_d = samp_size/pop_size` (dipangkas ke [0,1], 317-319); `hb_adj = f_d·ȳ_s + (1−f_d)·hb_pred` (323); `sd_adj = (1−f_d)·sd` (326) + floor eps (328); CI = `hb_adj ± 1.96·sd_adj` (330-331); clamp binomial [0,1], poisson ≥0 (333-338). Domain unsampled → `f_d = 0` ⇒ `hb_adj = hb_pred` (322).
- Tanpa `popsize_var`: pakai kuantil INLA apa adanya (344-350).
- `mse = sd²`, `rse = 100·sd/|hb|` (412-413). **Tanpa posterior sampling.**
- `random_effect = summary.random$..dom_id..$mean[1:n_domains]` (353-354).
- `unsampled_domains` dilaporkan + `cli_alert_info` bila print (461, 475-477).

### 2.5 Hyperpar & goodness
`grep("..dom_id..", rownames)` → `σ²_u = 1/prec[1]` (378-382); `grep("Gaussian")` → `σ²_e` (386-390); `grep("Phi", ignore.case)` → `phi` (394-397). `goodness` = **list** `dic, p_eff_dic, waic, p_eff_waic, mlik = fit$mlik[1,1]` (438-444). `hyperpar = fit$summary.hyperpar` apa adanya (454).

### 2.6 `@references` (83-94, persis)
```
\enumerate{
  \item Battese, G. E., Harter, R. M., & Fuller, W. A. (1988). An error-components model
    for prediction of county crop areas using survey and satellite data.
    \emph{Journal of the American Statistical Association}, 83(401), 28-36.
  \item Rao, J. N. K., & Molina, I. (2015). \emph{Small Area Estimation} (2nd ed.). John Wiley & Sons.
  \item Riebler, A., S\enc{ø}{o}rbye, S. H., Simpson, D., & Rue, H. (2016). An intuitive Bayesian spatial model
    for disease mapping that accounts for scaling. \emph{Statistical Methods in Medical Research}, 25(4), 1145-1165.
  \item Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian inference for latent Gaussian
    models by using integrated nested Laplace approximations.
    \emph{Journal of the Royal Statistical Society: Series B}, 71(2), 319-392.
}
```

---

## 3. `R/hb_twofold.R` (514 baris) + alias `hb_tfh`

### 3.1 Signature (139-160) + alias
```r
hb_twofold <- function(formula, vardir = NULL, domain, subarea = NULL, data,
  family = c("gaussian","binomial","poisson"),
  spatial = c("none","bym2","besag"), W = NULL, weight = NULL,
  trials = NULL, exposure = NULL,
  strategy = c("laplace","simplified.laplace"), scale_model = TRUE,
  prior_prec_area   = list(prior = "pc.prec", param = c(1, 0.01)),
  prior_prec_subarea= list(prior = "pc.prec", param = c(1, 0.01)),
  prior_phi = list(prior = "pc", param = c(0.5, 0.5)),
  compute_area = TRUE, n_samples = 200, print_result = TRUE, ...)
```
`hb_tfh <- hb_twofold` (513-514); `@aliases hb_tfh` (115); `NAMESPACE:12-13` mengekspor keduanya. Return class `c("fastsae_hb_twofold","fastsae_hb","fastsae")` (504), `invisible` (510).

### 3.2 Pseudocode bernomor
1. `.check_inla_installed()` (162); `match.arg` family/spatial/strategy (164-166); `domain` wajib (178-180); `subarea` NULL → `1..n` (183-187); `area_id = match(domain_vec, unique_domains)` (193).
2. `y` & `term_labels` dari `model.frame(na.pass)` (197-199); `vardir/weight/trials/exposure` via `.get_variable` (202-205).
3. Validasi likelihood (207-223): gaussian **wajib `vardir`** → `scale_vec = 1/vardir`, baris NA/≤0/inf diganti `1` (213-215); binomial wajib `trials`; poisson `scale_vec=NULL`.
4. W → `.convert_spatial_weights` (226-230). Data: `..area_id..`, `..subarea_id..=1..n_subareas` (232-236).
5. Formula (242-258): `y ~ <fixed> + <rand_area> + f(..subarea_id.., model = 'iid', hyper = list(prec = prior_prec_subarea))`, dengan `rand_area`:
   - none: `f(..area_id.., model = 'iid', hyper = list(prec = prior_prec_area))`
   - bym2: `f(..area_id.., model = 'bym2', graph = adj_graph, scale.model = <>, hyper = list(prec = prior_prec_area, phi = prior_phi))`
   - besag: `f(..area_id.., model = 'besag', graph = adj_graph, scale.model = <>, hyper = list(prec = prior_prec_area))`
6. Thread `num.threads="1:1"` global (261-265); args (268-276): `control.predictor=list(compute=TRUE, link=1)`, `control.inla=list(strategy=strategy)`, `control.compute=list(dic,waic,cpo,mlik, config = compute_area)`, `num.threads = 1`.
7. Likelihood (278-287): gaussian → `scale = scale_vec` + `control.family = list(hyper = list(prec = list(initial = 0, fixed = TRUE)))`; binomial → `Ntrials = as.integer(round(trials_vec))`; poisson → `E = exposure` (opsional). `...` menimpa (289-290).
8. `do.call(INLA::inla, …)` (293-301).

### 3.3 Sub-area & agregasi area
- Prediksi dari `summary.fitted.values[1..n_subareas]` + `summary.linear.predictor` (304-311); `mse = sd²`, `rse` (321-322); RE area/subarea diambil dengan pengindeksan `area_id`/`subarea_id` (314-318); `df_hb` (324-339) + alias `df_subarea` (486).
- **Bobot** (345-356): per area dinormalisasi; `weight=NULL` → `1/length(idx)`; NA/negatif → 0, bila `sum_w==0` → `1/length(idx)`.
- **Titik area**: `area_hb = Σ norm_weight·hb_pred` per domain (359).
- **Ketidakpastian area** (367-408): `ps <- INLA::inla.posterior.sample(n = n_samples, fit)` (368); baris diambil dari `rownames(ps[[1]]$latent)` dengan kunci `"Predictor:1".."Predictor:<n_subareas>"` (369-371); bila ada `NA` → gagal (373); transformasi balik link manual: binomial `1/(1+exp(-x))` (378), poisson `exp(x)` (380), gaussian tanpa transformasi; `area_draws = colSums(w_i * samples)` (387) → `sd`, `quantile 0.025/0.975` (388-390).
- **Fallback** (398-408) bila sampling gagal: `area_sd = sqrt(sum((w_i·hb_sd)^2))` (asumsi independensi) dan CI `±1.96·sd` (405-406). `mse_area = sd²`, `rse_area` (410-411), `df_area` (413-423); bila `compute_area=FALSE` → `df_area = NULL` (342).
- Hyperpar (444-466): `grep("..area_id..")` → `σ²_v = 1/prec`; `grep("..subarea_id..")` → `σ²_u`; `grep("Phi")` → `phi`; `random_effect_var = c(sigma2_v=, sigma2_u=)` (468). `goodness` = list (471-477).

### 3.4 `@references` (105-112, persis)
```
\enumerate{
  \item Torabi, M., & Rao, J. N. K. (2014). On small area estimation under a sub-area
    level model. \emph{Journal of Multivariate Analysis}, 127, 36-55.
  \item Rao, J. N. K., & Molina, I. (2015). \emph{Small Area Estimation} (2nd ed.). John Wiley & Sons.
  \item Riebler, A., S\enc{ø}{o}rbye, S. H., Simpson, D., & Rue, H. (2016). An intuitive Bayesian spatial model
    for disease mapping that accounts for scaling. \emph{Statistical Methods in Medical Research}, 25(4), 1145-1165.
}
```

---

## 4. `R/inla_utils.R` (397 baris) — semua `@noRd`, tidak diekspor

### 4.1 `.check_inla_installed()` (8-15)
`requireNamespace("INLA")`; kalau tidak ada → `cli_abort` + instruksi repo INLA.

### 4.2 `.convert_spatial_weights(W, n_domains, spatial="bym2", domain_names=NULL)` (27-88)
1. `W=NULL` → abort (28-30).
2. `listw` → `spdep::listw2mat` (33-38); `nb` → `spdep::nb2mat(style="B", zero.policy=TRUE)` (39-44); `matrix`/`Matrix` → `as.matrix` (45-46); **`inla.graph` → `return(list(graph=W, W_mat=NULL))`** tanpa validasi dimensi (47-48); lainnya abort (49-51).
3. Bila ada `rownames` & semua nama domain cocok → subset urut `domain_names` (54-59).
4. Cek dimensi = `n_domains × n_domains` (62-67).
5. `adj_mat <- (W_mat > 0 | t(W_mat) > 0) * 1; diag(adj_mat) <- 0` (72-73) — **graf biner simetris, bobot dibuang**.
6. `degree==0` → `cli_warn` (isolated nodes) (76-82); return `list(graph=adj_mat, W_mat=W_mat)` (84-87).

### 4.3 `.extract_inla_results(...)` (93-396)
Argumen: `fit, data, y, domain, time, family, spatial, temporal, st_interaction, vardir, trials, exposure, X_mat, rho_range, domain_id, time_id, unique_domains, unique_times, call`.
1. **Koefisien** (118-144): `spatial=="slm"` → baris `(n_domains+1):(n_domains+p_cov)` dari `fit$summary.random[[1]]` dengan `row.names = colnames(X_mat)` (119-132, sesuai dok slm **[E]**); selain itu `fit$summary.fixed`. Kolom: `beta, std.error=sd, zvalue=mean/sd (NA bila sd≤0), pvalue=2·pnorm(|z|,lower=FALSE), ci_lower=0.025quant, ci_upper=0.975quant`.
2. **Hyperpar** (146-206): `hyperpar = as.data.frame(fit$summary.hyperpar)` (147, bisa NULL).
   - `random_effect_var`: baris `grep("Precision for (\\.\\.domain_id\\.\\.|domain)")` → `[1]` → `1/mean` bila `mean>0` (158-168); fallback `grep("Precision for")` **hanya bila `spatial=="none" && temporal=="none"`** (159-162).
   - `random_effect_var_time`: `grep("Precision for (\\.\\.time_id\\.\\.|time)")` → `1/mean` (171-177).
   - `phi`/`rho`: `grep("Phi for")` → `phi=rho=mean[1]`; else `grep("Beta for")` (generic1/Leroux) → sama; else `grep("Rho for (\\.\\.domain_id\\.\\.|domain)")` (slm) → **re-scale linear** `rho = rho_range[1] + raw*(rho_range[2]-rho_range[1])` (196; sama dgn resep dok INLA **[E]**); `rho_time = grep("Rho for (\\.\\.time_id\\.\\.|time)")` (203-205).
3. **Random effect per observasi** (208-272): spasial `summary.random$..domain_id..$mean[1:n_domains]` lalu `[domain_id]` (217-227); temporal bila `st_interaction=="domain-specific"` **atau** `length(mean)==n_domains*n_times` → `(domain_id-1)*n_times + time_id` (233-241, **benar [E]**); selain itu `[1:n_times][time_id]` (242-248); total = spasial+temporal (253-258); fallback `summary.random[[1]]` (260-267); bila `length != n_obs` → `rep(NA_real_, n_obs)` (270-272).
4. **Fitted** (274-284): `summary.fitted.values[1:n_obs]` & `summary.linear.predictor[1:n_obs]`; `hb=mean`, `sd`, `mse=sd²`, `rse`, CI kuantil 0.025/0.975, `linear_pred=mean`.
5. **`df_hb`** (286-318): dengan kolom `time` bila `time` tidak NULL + `random_effect_spatial`/`random_effect_temporal` bila tersedia (302-303).
6. Kolom turunan (320-348): `vardir` (+ `precision` beta `pmax(y(1-y)/vardir-1,1)`, gamma `precision=y²/vardir`, `cv_dir=sqrt(vardir)/y`); `trials` (+ `estimated_total=hb*trials`, beta `precision=pmax(trials-1,1)`); `exposure` (+ `rate=hb`, `estimated_count=hb*exposure`).
7. `goodness` = **named numeric vector** `DIC, pD, WAIC, pWAIC, Marginal_LogLik=fit$mlik[1,1]` (351-357).
8. Output (371-395): `df_hb, hb, df_eblup, estcoef, hyperpar, random_effect_var, random_effect_var_time, phi, rho, rho_time, goodness, family, spatial, temporal, st_interaction, level="area", model="HB-<FAM> (<SPATIAL> + <TEMPORAL> + ST:<INTERACTION>)", convergence=TRUE, fit, call`; class `c("fastsae_hb_area","fastsae")`.

---

## 5. DISCREPANCIES (kode vs dokumentasi/konvensi)

- **D1 — `st_interaction="type2"` selalu gagal.** Dua `f()` memakai `..time_id..` (hb_area.R:539-541 dan 542-545). INLA: *"The covariate `..time_id..' are used in the following f()-terms: 1, 2 / Only one is allowed"* **[E]** ⇒ error, bukan model KH II.
- **D2 — `type3`/`type4` dgn `spatial != "none"` selalu gagal.** Term spasial ditambah dua kali: blok spasial (416/424/432/449/487, guard hanya `!= "separable"`) dan blok interaksi (551-554 / 561-564). INLA error sama (covariate `..domain_id..` dipakai f()-term 1 & 3) **[E]**.
- **D3 — `type3`/`type4` dgn `spatial="none"`** → interaksi tidak dibuat (guard 550, 560) **dan** tidak ada term `..domain_id..` (guard 411 gagal) ⇒ model = efek waktu saja, `st_interaction` diabaikan diam-diam (juga: `temporal=="none"` ⇒ seluruh blok 491-567 dilewati).
- **D4 — `hyper` salah untuk `type3`/`type4`** (552-553, 562-563 selalu `hyper = list(prec = prior_prec)`): `bym` → INLA error `Unknown keyword in 'hyper' 'prec'` **[E]**; `generic1` → `graph=` (bukan `Cmatrix`) → error *"For generic models the Cmatrix has to be provided"* **[E]**; `slm` → `graph=` tanpa `args.slm` → error **[E]**; `bym2` → `prior_phi` **diam-diam diabaikan**, phi jatuh ke default INLA `pc(0.5,0.5)` (uji hyper parsial: baris `Phi for ..domain_id..` tetap muncul) **[E]**. (Saat ini ketiganya tertimpa D2 lebih dulu.)
- **D5 — `spatial="slm"` tidak punya guard `!= "separable"`** (487 vs 415/423/431/448) ⇒ slm + separable menambah 2 term `..domain_id..` (slm + grouped iid, 528-531) → error mekanisme D2. Selain itu slm + separable menurunkan spasial ke `iid`.
- **D6 — argumen `link` mati.** Dok (51-52) menjanjikan pilihan link; `link` diteruskan (295, 333) tetapi **tidak pernah dipakai** di `.fit_inla_area` (satu-satunya `link` di 587 adalah `control.predictor$link = 1`). ⇒ selalu link default family.
- **D7 — `prior_phi` dipakai ganda** (443-447 generic1 `beta`, 481-483 slm `rho`) hanya bila `prior_phi$prior != "pc"`. Akibat default (`prior="pc"`): `beta` generic1 memakai default INLA `gaussian(0,0.1)` pada skala logit(β) (bukan prior paket, dok generic1), sementara slm tetap dapat `logitbeta(1,1)` (484). Guard juga memblokir user yang mengisi `prior_phi = list(prior="pc.cor0", …)` utk beta/rho.
- **D8 — prior default tidak konsisten antar fungsi** utk komponen yang sama: `hb_area` `pc.prec(1,0.01)` (165) & `hb_unit` **`loggamma(0.01,0.01)`** (136) — dok hb_unit (45-46) memang menulis loggamma, tetapi paket lain memakai PC prior; `hb_twofold` `pc.prec(1,0.01)` (153-154).
- **D9 — `random_effect_var` utk `spatial="bym"` mengambil komponen yang salah.** `inla_utils.R:158` mencocokkan dua baris — urutan INLA (uji **[E]**): `Precision for the Gaussian observations || Precision for id (iid component) || Precision for id (spatial component)` — `[1]` = **iid component**, jadi `random_effect_var` = σ² tak-berstruktur, bukan varian spasial (baris spasial = `[2]`).
- **D10 — fallback regex hyperpar berisiko salah sasaran.** `inla_utils.R:161` `grep("Precision for")` mengembalikan baris pertama; utk gaussian tanpa `vardir` baris pertama = `"Precision for the Gaussian observations"` **[E]** ⇒ akan melaporkan σ²_error sbg `random_effect_var`. Saat ini tidak tercapai karena index selalu `..domain_id..` (baris 412), tetapi rawan bila pemanggilan berubah.
- **D11 — `scale.model` tidak seragam.** Diset untuk bym2/bym/besag (417/425/433) dan rw1/rw2 (492), **tidak** untuk generic1 (449), slm (487), iid, ar1 (492) — sementara dok (53-54) menyebutnya "recommended for BYM2/Besag".
- **D12 — `constr` tidak pernah dioper** → mengikuti default INLA (bym2/bym/besag/rw1/rw2 TRUE; iid/ar1/generic1/slm FALSE) **[E]**; tidak didokumentasikan di roxygen.
- **D13 — `...` bisa bikin error argumen ganda.** hb_area menyusun `inla_args <- list(..., ...)` (593-601) tanpa menimpa; `inla()` tidak punya `...` **[E]** ⇒ user yang mengoper `control.compute`/`control.family`/`num.threads` lewat `...` mendapat `formal argument matched by multiple actual arguments`. hb_unit (284) & hb_twofold (290) menimpa (`inla_args[names(extra_args)] <- extra_args`) — perilaku berbeda.
- **D14 — `method="laplace"` sering diam-diam jatuh ke INLA**: gerbang 258 hanya utk non-spatial non-temporal binomial/poisson; gaussian/spatial/temporal → INLA dengan `strategy` dipaksa `"laplace"` (269-276).
- **D15 — jalur lme4 tidak menghasilkan ketidakpastian**: `sd/mse/rse/ci_* = NA` (777-781), `goodness = AIC/BIC/LogLik` (786-790) dan `hyperpar = data.frame(Parameter, Estimate)` (797) — format berbeda dari jalur INLA (`goodness` named vector, `hyperpar` = summary.hyperpar); `df_hb` juga tanpa kolom `time`.
- **D16 — kelas tidak konsisten**: `hb_area` → `c("fastsae_hb_area","fastsae")` (394/811); `hb_unit`/`hb_twofold` menambah `fastsae_hb` (471/504) → `hb_area` tidak punya kelas `fastsae_hb`.
- **D17 — domain unsampled**: hb_unit sengaja menambah baris `y=NA` per domain (244-249) + pesan (476); hb_area/hb_twofold **tidak** membuat baris ⇒ area yang tidak ada di `data` tidak pernah diprediksi; area dengan `y=NA` tetap dievaluasi (275-276 / 304).
- **D18 — tipe CI berubah-ubah di `hb_twofold`**: jalur sampling memakai kuantil posterior (389-390), fallback memakai normal `±1.96` (405-406); demikian pula `hb_unit` + `popsize_var` mengganti kuantil INLA dgn `±1.96·sd_adj` (330-331).
- **D19 — `hb_unit` + `popsize_var`**: kolom `linear_pred` tetap η model tanpa FPC sementara `hb` sudah dicampur `f_d·ȳ_s + (1−f_d)·pred` (323 vs 419) ⇒ dua kolom tidak lagi saling invers; `sd_adj` diskalakan `(1−f_d)` (326) — aproksimasi, bukan varian postifik-pengamatan.
- **D20 — `.convert_spatial_weights` membakar bobot ke graf biner** (72-73) tetapi `W_mat` mentah (dipakai slm, hb_area 452-455) dikembalikan apa adanya; jalur `inla.graph` (47-48) melewati cek dimensi dan mengembalikan `W_mat=NULL` → `spatial="slm"` dgn `W` berupa `inla.graph` akan gagal di `rowSums(W_mat)` (453).
- **D21 — petunjuk index**: `st_interaction="domain-specific"` (239) dan ekstraksi koefisien slm (123) **terverifikasi benar** (bukan discrepancy; dicatat agar tidak dianggap bug). `slm` rho re-scale (196) juga cocok dgn resep dok INLA **[E]**.
