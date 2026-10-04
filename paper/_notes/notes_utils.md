# Catatan Teknis — Lapisan Utilitas `fastsae` (benchmark, diagnose, compare, export, map, methods, sim)

Semua klaim di bawah berasal dari pembacaan langsung kode. Nomor baris = `R/<file>:<baris>`.
Aturan keras dipatuhi: tidak ada file di luar `paper/_notes/` yang diubah.

---

## 1. `R/benchmark.R` — kalibrasi benchmark

### 1.1 API & argumen
- Generic `benchmark_sae(object, ...)` (`R/benchmark.R:120-122`); alias `benchmark()` → `benchmark_sae()` (`:264-266`), plus method `benchmark.fastsae`/`benchmark.default` (`:269-276`).
- `benchmark_sae.fastsae(object, target=NULL, weight=NULL, method=c("ratio","difference","optimal","logit"), group=NULL, type=c("mean","total"), national_target=NULL, outer_target=NULL, ...)` (`:126-136`).
- `benchmark_sae.default` sama, `object` harus numerik (`:212-231`); `mse = NULL` dikirim ke worker (`:256`), `model_family="generic"` (`:258`).
- `nat_target <- outer_target %||% national_target` (`:139`, `:225`) → `outer_target` adalah alias eksplisit (`:59`).
- Ekstraksi dari model: `df_est <- object$df_hb %||% object$df_ebp %||% object$df_eblup` (`:142`); `domain_vec <- df_est$domain %||% df_est$area %||% seq_len(n)` (`:147`); `y_hat <- df_est$hb %||% ebp %||% eblup %||% est` (`:148`); `mse_vec <- df_est$mse %||% (sd^2)` (`:153`).
- `weight`: string kolom (dicari di `object$data` lalu `df_est`, `:157-164`), vektor numerik panjang-n (`:165-169`), atau `NULL` → bobot seragam `rep(1, n)` (`:170-173`; default seragam juga di `:235`).
- `group`: string kolom (`:178-185`) atau vektor (`:186-191`); jika `NULL`, worker memakai `group = "All"` (`:311`).

### 1.2 Metode dan rumus persis yang dikodekan (semua di `.benchmark_worker`)
Persetujuan agregat per grup: `w_calc <- if (type=="mean") w_g/sum(w_g) else w_g` (`:451`), `agg_initial <- sum(w_calc*y_g)` (`:452`). Karena itu `type="mean"` = rata-rata tertimbang ternormalisasi, `type="total"` = total populasi (`:52-54`, `:451`).

| method | rumus di kode | baris |
|---|---|---|
| `"ratio"` | `ratio_adj <- T_g/agg_initial`; `y_bm <- y_g*ratio_adj`. Error bila `|agg_initial| < .Machine$double.eps` | `:454-460` |
| `"difference"` | `diff_val <- T_g - agg_initial`; `y_bm <- y_g + diff_val/sum(w_calc)` | `:462-465` |
| `"optimal"` | fallback ke difference + `cli_warn` bila `mse` NA atau semua `<=0`; selain itu `inv_prec <- pmax(mse_g, 1e-8)`; `lambda <- diff_val/sum(inv_prec)`; `y_bm <- y_g + (inv_prec/w_calc)*lambda` | `:467-480` |
| `"logit"` | wajib `0<y<1` (`:484-489`) dan (untuk `type="mean"`) `0<T_g<1` (`:490-492`); `logit_y <- qlogis(y_g)`; akar `f_obj(alpha) <- sum(w_calc*plogis(logit_y+alpha)) - T_g`; bracket `[-20,20]`, jika tanda sama → `[-50,50]` (`:499-506`); `uniroot(..., tol=1e-9)` (`:508`); `y_bm <- plogis(logit_y + alpha_opt)` (`:510`) |

Kepatuhan constraint terverifikasi secara algebra untuk keempatnya: Σ`w_calc`·`y_bm` = `T_g` (mis. ratio: `agg*(T/agg)=T`).

### 1.3 Dua tahap hierarkis (target nasional)
Aktivasi: `is_hierarchical <- !is.null(national_target) && n_groups > 1` (`:318`).

**Agregat per grup** (`:328-341`): `W_g=sum(w_g)`; `w_norm=w_g/W_g`; `group_orig[g]=sum(w_norm*y_hat)` (SELALU ternormalisasi, lihat ketidakcocokan D-8); `group_mse[g]=sum(w_norm^2 * pmax(mse,1e-8,na.rm=TRUE))` bila mse tidak semua NA, else `NA` (`:335-340`).

**Parsing `target`** (`:349-393`): data.frame harus punya kolom `group`,`target` (`:351-361`); vektor numerik: skalar & 1 grup (`:363-364`), skalar & >1 grup non-hierarkis → di-recycle (`:365-366`), named vector (`:367-373`), panjang == n_groups (`:374-375`), selain itu error (`:377`,`:380`). Jika `target=NULL` & hierarkis → `target_map <- group_orig` (`:383-386`); jika `national_target` ada tapi hanya 1 grup → target grup = `national_target` (`:387-390`); jika keduanya tak ada → error (`:391-393`).

**STAGE 1 — harmonisasi grup ke target nasional** (`:395-438`):
- `W_share <- group_weights/sum(group_weights)` (`:399`)
- `nat_current_agg <- if (type=="mean") sum(W_share*target_map) else sum(target_map)` (`:402`); `diff_nat <- national_target - nat_current_agg` (`:403`); hanya dieksekusi bila `|diff_nat| > 1e-7` (`:405`).
- ratio: `target_map <- target_map * (national_target/nat_current_agg)` (`:407-408`)
- difference: tambah `diff_nat` (mean) atau `diff_nat/n_groups` (total) (`:410-411`)
- optimal: `inv_prec_g <- pmax(group_mse,1e-8)`; bila ada NA → semua diset 1 (`:413-414`); `lambda_g <- diff_nat/sum(inv_prec_g)`; `target_map <- target_map + (inv_prec_g/(W_share atau 1))*lambda_g` (`:415-416`)
- logit: `logit_t <- qlogis(target_map)`; akar `sum(W_share*plogis(logit_t+a)) - national_target` via `uniroot(interval=c(-50,50), tol=1e-9)` (`:418-423`)
- `stage1_df` (kolom: `group, n_domains, weight_total, weight_share_pct=round(W_share*100,3), original_group_estimate, initial_target, calibrated_target, adjustment`) (`:427-437`)

**STAGE 2 — kalibrasi level kecil dalam tiap grup** (`:440-512`): `df_work$target <- target_map[df_work$group]` (`:441`), lalu 4 metode di §1.2.

### 1.4 Objek return
`data.frame` dengan kelas `c("fastsae_benchmark","data.frame")` (`:571`), kolom: `domain`, (`group` bila ada, `:519-521`), `weight`, `original`, `benchmarked`, `adjustment = benchmarked-original` (`:525`), `rel_adjustment_pct = adjustment/original*100` atau NA bila `|original| <= eps` (`:526-527`), `target`, dan `national_target` (berulang) bila hierarkis (`:529-531`).
Atribut (`:562-569`): `method`, `type`, `hierarchical`, `group_var`, `stage1_summary`, `verification` (list per grup: `target, original_sum, benchmarked_sum, discrepancy=|bm_sum-target|`, `:541-546`), `national_verification` (strukur sama, `:554-559`), `call`.

### 1.5 Metode terkait
- `print.fastsae_benchmark` (`:577-637`): judul hierarkis/non-hierarkis (`:581-584`), label metode (`:586-592`), status `[CONSISTENT]` bila `discrepancy < 1e-5` selain itu `[APPROX]` (`:602-605`, `:618-619`), hanya 6 grup pertama bila >8 dan `all_groups=FALSE` (`:614`), `head(...,6)` domain + jumlah sisanya (`:629-633`).
- `summary.fastsae_benchmark` (`:641-660`): print + `stage1_summary` head-10 (`:646-651`) + `summary(adjustment)` & `summary(rel_adjustment_pct)` (`:654-658`).
- `autoplot.fastsae_benchmark` (`:664-700`): scatter `original` vs `benchmarked`, garis 1:1 putus-putus (`:672-673`), warna per `group` bila >1 grup (`:675-680`), legend disembunyikan bila >10 grup (`:695-697`).
- `plot.fastsae_benchmark` (`:704-710`): `sf_geom` di `...` → `map_sae(x, ...)`, else `autoplot`.
- `@references` (disalin persis, `:79-93`):
  1. Rao, J. N. K., and Molina, I. (2015). *Small Area Estimation* (2nd ed.). John Wiley & Sons. Chapter 10: "Benchmarking and Other Practical Issues", pp. 297-315.
  2. You, Y., and Rao, J. N. K. (2002). A pseudo-empirical best linear unbiased prediction approach to small area estimation using survey weights. *The Canadian Journal of Statistics*, 30(3), 431-439.
  3. Steorts, R. C., Hall, P., and Ghosh, M. (2014). General benchmarking under quadratic loss with applications to small area estimation. *Journal of Survey Statistics and Methodology*, 2(2), 173-193.
  4. Datta, G. S., Ghosh, M., Steorts, R., and Maples, J. (2011). Bayesian benchmarking with applications to small area estimation. *Test*, 20(3), 574-588.
  5. Berg, E., and Fuller, W. A. (2014). Small area prediction of proportions with a constrained multinomial logit model. *Journal of Survey Statistics and Methodology*, 2(3), 256-283.

---

## 2. `R/diagnose.R` — diagnostik

### 2.1 API
`diagnose(object, W=NULL, truth=NULL, rse_threshold=25, alpha_level=0.05)` (`R/diagnose.R:53-57`) — fungsi biasa, **bukan** generic S3. Wajib kelas `"fastsae"` (`:58-60`).

### 2.2 Ekstraksi & kuantitas dasar
- df dari `df_hb %||% df_ebp %||% df_eblup` (`:63-71`); kolom estimasi `hb`/`ebp`/`eblup` (`:73`).
- `N_sampled = sum(!is.na(direct_y))`, `N_unsampled = N_total - N_sampled` (`:84-87`).
- `direct_rse = sqrt(vardir)/|y| * 100`, hanya utk area ter-sampel dengan `vardir>0` dan `y!=0` (`:90-92`).
- `eff_ratio = vardir/mse` (rasio efisiensi; >1 berarti SAE lebih presisi) (`:95-97`).
- `residuals = direct_y - sae_pred` (`:100`); `std_residuals = residuals/sqrt(vardir)` (`:101-103`).
- Precision summary (`:114-130`): `prop_reliable = mean(rse_clean < rse_threshold)*100` (`:109`), `prop_gain = mean(eff_ratio>1)*100` (`:112`), mean/median direct & SAE RSE, `eff_ratio_summary` (Min,Q1,Median,Mean,Q3,Max, `:121-128`).

### 2.3 Tes Brown et al. (2001) — rumus persis (`:133-191`)
Syarat: minimal 5 pasang (y, ŷ) ter-sampel (`:138`).
- **(a) Regresi bias**: `fit_ols <- lm(y_sub ~ pred_sub)` → `alpha=coefs[1]`, `beta=coefs[2]`, `vcov=vcov(fit_ols)` (`:144-146`). **OLS polos tanpa bobot.**
- **Wald H0: (alpha,beta)=(0,1)**: `diff_vec <- coefs - c(0,1)` (`:149`); `f_stat <- t(diff_vec) %*% solve(vcov) %*% diff_vec / 2` (`:152`, dibagi 2); `df1=2`, `df2=n_sub-2` (`:153-154`); `p_val <- pf(f_stat, 2, n_sub-2, lower.tail=FALSE)` (`:155`); `is_unbiased <- (p_val >= alpha_level)` (`:164`). Bila `solve(vcov)` gagal, tes dilewati (`:150-151`).
- **(b) Statistik GOF chi-kuat W**: `denom <- vardir + mse` (area valid bila keduanya non-NA dan `(vd+mse)>0`, minimal 3, `:169-172`); **`W <- sum(((y - pred)^2)/denom)`** (`:173-174`); `df_w <- sum(valid_w)` (`:175`); `p_value <- pchisq(W, df_w, lower.tail=FALSE)` (`:176`); penerimaan H0 memakai **interval dua sisi**: `crit_low=qchisq(alpha/2, df)`, `crit_high=qchisq(1-alpha/2, df)`, `is_good_fit <- (W >= crit_low && W <= crit_high)` (`:178-180`).

### 2.4 Moran's I pada residual — rumus persis (`:193-246`)
Matriks: argumen `W` → `object$W` → NULL (`:197-203`). Hanya dijalankan bila `nrow(W_mat)==N_total` (`:207`), lalu **dipotong ke area ter-sampel**: `res_clean <- residuals[sampled_mask]`, `W_sub <- W_mat[sampled_mask, sampled_mask]` (`:208-209`); syarat `n_spat>=5` dan semua residual non-NA (`:212`), dan `s0=sum(W_sub)>0` (`:216`).
- `z <- res - mean(res)` (`:214`)
- **`I <- (n/s0) * (sum(W * outer(z,z)) / sum(z^2))`** (`:217-219`)
- `E[I] <- -1/(n-1)` (`:222`)
- `s1 <- 0.5*sum((W + t(W))^2)` (`:223`); `s2 <- sum((rowSums(W)+colSums(W))^2)` (`:224`); `kurt <- n*sum(z^4)/(sum(z^2)^2)` (`:225`)
- **`Var(I) <- [ n*((n^2-3n+3)*s1 - n*s2 + 3*s0^2) - kurt*((n^2-n)*s1 - 2n*s2 + 6*s0^2) ] / [(n-1)(n-2)(n-3)*s0^2] - E[I]^2`** (`:227-229`); `var_I <- max(1e-8, var_I)` (`:231`)
- `z_stat <- (I - E[I])/sqrt(var_I)`; `p <- 2*pnorm(-|z_stat|)`; `no_residual_autocorrelation <- (p >= alpha_level)` (`:232-241`)
Catatan: residual yang dipakai adalah residual mentah direct−pred (bukan standardized, `:100`,`:208`), tanpa koreksi randomisasi.

### 2.5 Metrik simulasi (bila `truth` diberikan, `:251-271`)
- `rb <- (sae_pred - truth)/truth * 100` (`:255`)
- `rrmse <- sqrt((sae_pred-truth)^2)/|truth| * 100` (= |error relatif| per domain) (`:256`)
- `coverage_rate <- mean(truth>=ci_lower & truth<=ci_upper)*100` bila kolom `ci_lower`/`ci_upper` ada, else NA (`:259-262`)
- Keluaran: `mean_rb, mean_abs_rb=mean(|rb|), mean_rrmse=mean(rrmse), coverage_rate` (`:264-269`)

### 2.6 Diagnostik Bayes (opsional, `:276-332`)
Deteksi `is_bayesian` = kelas `fastsae_hb_area`/`fastsae_ebp_area` atau ada `WAIC`/`DIC` di `object$goodness` (`:281`). Ambil `WAIC, pWAIC, DIC, pD, Marginal_LogLik` (`:284-288`); CPO/PIT dari `object$fit$cpo` (`:277`, `:293-310`); `failures = sum(cpo$failure>0)` (`:303`); **tes KS** `ks.test(pit_clean,"punif",0,1)` dengan `is_calibrated <- (p >= alpha_level)`, syarat ≥5 nilai PIT di (0,1) (`:312-320`).

### 2.7 Threshold & status (`:334-349`)
- `is_unbiased` (tes Wald ada & lolos), `is_fit` (GOF ada & lolos) — keduanya dianggap lolos bila tesnya tidak tersedia (`:337-338`).
- `has_high_precision <- is.na(prop_reliable) || prop_reliable >= 80` (`:339`).
- Urutan status (`:341-349`): semua lolos → **`"PASS: Model is well-calibrated, statistically unbiased, and reliable."`**; `!is_unbiased` → **`"CAUTION: Potential systematic bias detected between direct and model predictions."`**; `!is_fit` → **`"CAUTION: Model goodness-of-fit indicates notable deviation from survey variance."`**; selain itu → **`"WARNING: High proportions of domains with RSE >= threshold."`**
- Default: `rse_threshold = 25` (`:56`), `alpha_level = 0.05` (`:57`). Ambang lain di `print`: `prop_reliable>=80` (High precision, `:405-409`), `prop_gain>=80` (sukses efisiensi, `:421-424`), `|mean_rb|<5` (low bias, `:470-474`), `coverage_rate>=90` (nominal terpenuhi, `:477-481`).

### 2.8 Return & print
Return kelas `c("fastsae_diagnose","list")` (`:386`) berisi: `model_info` (`:371-376`), `precision` (`:377`), `brown_test` (`:378`), `spatial_test` (`:379`), `bayesian_metrics` (`:380`), `simulation_metrics` (`:381`), `df_diag` (`:382`; kolom `domain, direct_y, sae_pred, vardir, mse, direct_rse, sae_rse, eff_ratio, residuals, std_residuals, is_sampled` + `cpo`,`pit` bila ada, `:351-368`), `status` (`:383`).
`print.fastsae_diagnose` (`:391-525`) mencetak 5 bagian (Precision & Efficiency, Brown Calibration, Moran's I, Simulation, Bayesian) + Final Assessment (ikon sukses/warn/danger dari prefiks PASS/CAUTION/WARNING, `:515-521`).
`@references` (disalin persis, `:27-33`):
1. Brown, G., Chambers, R., Heady, P., & Heasman, D. (2001). Evaluation of small area estimation methods: An application to the British Labour Force Survey. *ONS Internal Report*.
2. Rao, J. N. K., & Molina, I. (2015). *Small Area Estimation* (2nd ed.). John Wiley & Sons.
3. Cliff, A. D., & Ord, J. K. (1981). *Spatial Processes: Models & Applications*. Pion London.

---

## 3. `R/compare_sae.R` — perbandingan dua model

- Generic `compare_sae(model1, model2=NULL, names=NULL, thresholds=c(20,30), ...)` (`:62-64`); **hanya ada method `.default`** (terdaftar `NAMESPACE:33`) — `model1` boleh list 2 model (`:70-76`).
- Penamaan otomatis via `.get_model_short_name` (peta `FH/SFH/ST/BHF/TFH/TWOFOLD/EBP` → label manusia, `:447-464`); bila sama ditambah "(1)"/"(2)" (`:86-89`).
- Merge kunci: `domain`, atau `c(domain,subarea)` bila kolom `subarea` ada di keduanya, atau `c(domain,time)` (deteksi `time|period|year` di `.extract_domain_data`, `:413-416`), `merge(..., all=FALSE)` — inner join, error bila 0 overlap (`:99-109`).
- **Yang dibandingkan**: estimasi (kolom `estimate_1/2`), `mse_1/2`, `rse_1/2` (fallback `rse = sqrt(mse)/|est|*100`, `:403`). Tidak ada argumen `truth` → **tidak ada RMSE/bias vs kebenaran**. Yang ada adalah *RMSD* antar dua model: `diff = est2-est1`, `abs_diff`, `rmsd = sqrt(mean(diff^2))`, `mae = mean(|diff|)` (`:112-135`).
- Metrik (`:132-176`): Pearson `r` (`:132`), Spearman `rho` (`:133`), MAE, RMSD, `mse_ratio = mse_1/mse_2` (`:118`), mean & median rasio MSE, jumlah & % domain dengan `mse_2 < mse_1` (`:139-140`), mean RSE masing-masing model, `rse_gain = mean_rse1 - mean_rse2` (`:144`). Semua disusun sebagai `data.frame(Metric, Value)` dengan format string (`:146-176`).
- **Kelompok/grup**: TIDAK ada argumen pengelompokan; satu-satunya "kelompok" adalah kunci merge (domain/time/subarea) dan label model. `thresholds` hanya dipakai untuk garis horizontal pada plot tipe `"rse"` (`:343-348`) — tidak dipakai perhitungan metrik.
- Return: `structure(list(metrics, data, names, thresholds, model1, model2), class="fastsae_comparison")` (`:178-188`).
- `print` (`:193-205`), `summary` (`:209-234`; kuantil estimasi & selisih, kuantil RSE bila ada), `plot` → `ggplot2::autoplot` (`:238-240`).
- `autoplot.fastsae_comparison(object, type=c("scatter","comparison","difference","mse","rse"), title=NULL)` (`:244-249`):
  - `scatter`: titik warna `diff` dengan `scale_color_gradient2`, garis 1:1, `coord_fixed`, subjudul berisi `r` (`:255-276`).
  - `difference`: stem/lollipop per domain terurut menurut `diff`, garis 0 (`:278-306`).
  - `mse`: bar dodged per domain per model (`:308-331`).
  - `rse`: line+point per model + 2 garis ambang `thresholds` (`:333-359`).
  - `comparison` (default): titik+garis terdodge per model (`:361-386`).

---

## 4. `R/export_sae.R` — flag reliabilitas & ekspor

### 4.1 `add_reliability_flags(object, rse_col="rse", thresholds=c(20,30))` (`:33`)
- Validasi `thresholds` panjang-2 menaik (`:34-36`). Ambil df dari model atau `as.data.frame(object)` (`:38-39`).
- Bila kolom RSE tak ada, hitung dari `mse`: `sqrt(mse)/|est|*100` (`:45-53`).
- Tier via `cut(breaks=c(-Inf,20,30,Inf), right=FALSE)` → label `"Reliable (< 20%)"`, `"Use with Caution (20-30%)"`, `"Unreliable (≥ 30%)"` (`:55-66`) — konsisten dengan doc `:14-16`.
- Return: **`df` saja** (bukan objek model yang dimodifikasi); penugasan `object$df_hb <- df` di `:69-71` tidak dipertahankan (R copy-on-modify) (`:68-75`).

### 4.2 `export_sae(object, file, benchmark=NULL, thresholds=c(20,30), overwrite=TRUE, ...)` (`:119-126`)
- Ekstensi hanya `xlsx`/`csv` (`:131-134`); `overwrite=FALSE` + file ada → error (`:136-138`).
- Sheet: `Estimates` (`:144`), `Model_Summary` (`:147`), `Benchmarked` (`:150-153`) dan **`Benchmarked_Groups`** dari `attr(,"stage1_summary")` (`:153-156`) → bisa 4 sheet.
- `.build_estimates_sheet` (`:192-255`): utk objek benchmark → kolom `Domain, Original_Estimate, Benchmarked_Estimate, Adjustment, Adjustment_Pct` (`:196-203`); utk model → `Domain, Direct_Estimate, Model_Estimate, Standard_Error=sqrt(mse), MSE, RSE_Pct, CI_Lower_95=est-1.96*se, CI_Upper_95=est+1.96*se (fallback bila tak ada kolom CI, :226-227), Reliability_Flag` (cut yang sama, `:234-239`), dengan pembulatan (`:241-252`).
- `.build_model_summary_sheet` (`:259-301`): baris `Parameter/Value` — Package, Model Type, Family, Method, Call, Convergence, n_domains (`:269-275`), komponen ragam `sigma2_u/rho/rho_t/phi` (`:278-281`), `goodness` apa pun (`:284-288`), koefisien `Beta: <term>` = `beta (SE: …, p: …)` dari `estcoef` kolom `beta/std.error/pvalue` (`:291-298`).
- Penulisan: xlsx → `writexl::write_xlsx` dulu, fallback `openxlsx::write.xlsx`, fallback terakhir CSV hanya sheet Estimates (`:160-172`); csv → Estimates + `<nama>_summary.csv` (`:173-181`). Sheet `Benchmarked` **tidak** ditulis pada mode csv.
- Return: `invisible(sheets)` — named list data frame (`:183`).

---

## 5. `R/map_sae.R` — peta choropleth

- Generic `map_sae(object, ...)` (`:93-95`); methods: `.fastsae` (`:99-174`), `.fastsae_benchmark` (`:178-265`), `.list` (`:269-304`), `.default` (data frame, `:308-369`). Wajib paket `sf` (`:113-115`).
- Argumen `.fastsae`: `model2, sf_geom, key, type=c("estimate","rse","reliability","comparison","difference"), indicator (alias type, :117), palette, thresholds=c(20,30), facet_scales="fixed", title, subtitle, ...` (`:100-111`).
- Geometri diambil dari `sf_geom` atau `object$data/df_hb/df_ebp/df_eblup` yang berkelas `sf` (`:139-150`).
- Pencocokan kunci `.match_spatial_keys` (`:507-577`): (i) `key` string tunggal = kolom di `sf_geom` (`:515-520`); (ii) `key` named = pasangan `c(model_col=sf_col)` (`:521-533`); (iii) otomatis: nama persis → daftar kandidat `c("domain","area","id","code","kd_kab","kode","kabupaten","provinsi","district","region")` (case-insensitive, `:543-548`) → **overlap set terbesar dan diadopsi bila `>= 40%` jumlah domain** (`:550-565`); gagal → error (`:569-574`).
- `.merge_model_with_sf` (`:581-593`): kunci sementara `..merge_key..`, `merge(all.x=TRUE)`.
- `.generate_sae_map` (`:597-702`) — kelas output **`ggplot`**:
  - `"estimate"`: `geom_sf` + `scale_fill_viridis_c(option=pal|"viridis")` (`:607-619`).
  - `"rse"`: viridis `option="magma", direction=-1` (`:621-633`).
  - `"reliability"`: skala manual 3 warna `#198754` (hijau), `#ffc107` (amber), `#dc3545` (merah) sesuai label tier (`:635-651`).
  - `"comparison"`: facet `~Source` Direct Survey vs SAE; bila `direct` semua NA → warning dan dipaksa ke `"estimate"` (`:653-675`).
  - `"difference"`: `diff_val = estimate - direct` dengan `scale_fill_gradient2(low="#2b83ba", mid="#ffffbf", high="#d7191c")`; error bila `direct` tak ada (`:677-695`).
- `.extract_model_data_for_map` (`:454-503`): kolom `domain/direct/estimate/mse/rse`; flag reliability **hard-coded** `breaks=c(-Inf,20,30,Inf)` (`:482-487`); label model dari `object$model`/`family` (`:489-499`).
- Dua model via `model2`/list → `.map_two_models` (`:377-450`): `"difference"` = est2 − est1 dengan gradient2 (`:414-428`); `"comparison"` = facet `~Model` viridis (`:430-449`); pencocokan model2 memakai `rownames(df2) <- domain` lalu indeks karakter (`:409-411`) — domain model2 yang tidak ada di hasil merge jadi NA.
- Objek benchmark (`:178-265`): `type=c("comparison","estimate","difference")`, `comparison` = facet Original vs Benchmarked (`:215-233`), `difference` = isi `adjustment` (`:235-249`), `estimate` = nilai `benchmarked` (`:251-264`); `sf_geom` wajib (`:208-210`).
- `map_sae.default` (data frame): kolom domain kandidat `c("domain","area","id","code","kd_kab")` (`:327-333`), kolom estimasi `c("hb","ebp","eblup","estimate","est","y_hat","pred")` (`:335-340`), RSE dihitung dari `mse` bila perlu (`:342-348`).

---

## 6. `R/methods.R` + `R/autoplot.R` — metode S3

Helper `%||%` didefinisikan di `R/methods.R:6`.

**Dari `methods.R`:**
| method | baris | perilaku singkat |
|---|---|---|
| `print.fastsae` | `:15-109` | `cli_h1`, Call, Convergence (+`n_iter`), tipe model dari `x$model` (`:34-43`), temporal/ST, method, varian ragam, `rho/rho_time/phi` (`:58-78`), `printCoefmat(estcoef)` (`:80-87`), `head(est_df,6)` + sisa baris (`:89-105`) |
| `summary.fastsae` | `:118-146` | membangun list `summary.fastsae` (20 komponen, `:119-143`) |
| `print.summary.fastsae` | `:155-256` | Call, konvergensi, komponen ragam/`estvarcomp`, koefisien, hyperpar, goodness of fit, `summary()` kolom estimasi terpilih (`:243-252`) |
| `coef.fastsae` | `:265-273` | vektor bernama `estcoef$beta`, fallback `fit$beta` |
| `fitted.fastsae` | `:282-293` | `df_hb$hb` / `df_ebp$ebp` / `df_eblup$eblup` |
| `residuals.fastsae` | `:302-316` | `y - hb/ebp/eblup`, fallback `residuals(fit$lme)` |
| `plot.fastsae` | `:332-344` | `y` = model fastsae + (`sf_geom`/data sf) → `map_sae(x, model2=y)`; `y` model tanpa sf → `ggplot2::autoplot(list(M1,M2))`; tanpa `y` + sf → `map_sae(x)`; else `ggplot2::autoplot(x)` |

**Dari `autoplot.R`:**
- Re-export `ggplot2::autoplot` (`:9`), `@import ggplot2`, `@importFrom rlang .data` (`:1-3`).
- `autoplot.fastsae(object, type=c("comparison","ribbon","mse","estimates","scatter","rse","map"))` (`:46-63`): `map`→`map_sae`; `comparison` & `ribbon` → `.autoplot_single_comparison` (poin + ribbon CI; CI dibuat `est ± 1.96*sqrt(mse)` bila belum ada, `:93-99`); `mse` → bar MSE per domain (`:436-464`, error bila mse tak ada `:439-444`); `estimates` → scatter `y` vs prediksi + garis 1:1 + garis rata-rata putus-putus, wajib kolom `y` (`:128-176`); `scatter` single → warning lalu fallback ke estimates (`:179-182`); `rse` → RSE SAE vs RSE direct + 2 ambang default `c(20,30)` + anotasi teks (`:185-244`).
- `autoplot.list(object, type=c("comparison","mse","scatter","map"))` (`:67-81`): bila bukan semua fastsae → `NextMethod()` (`:68-71`); helper multi: `.autoplot_multi_comparison` (errorbar CI, butuh `names(x)`, `:251-304`), `.autoplot_multi_mse` (dodge bila domain identik, else facet, `:307-368`), `.autoplot_multi_scatter` (tepat 2 model, r di judul, `:371-429`).
- `autoplot.fastsae_diagnose(object, type=c("all","calibration","rse","residuals","qq","pit","cpo"))` (`:486-630`): `calibration` = direct vs SAE + 1:1 + `geom_smooth(lm)` (`:493-507`); `rse` = direct RSE vs SAE RSE + garis ambang dari `object$precision$rse_threshold` (`:509-524`); `residuals` = std residual vs fitted + garis 0, ±2 (`:526-540`); `qq` = `stat_qq`/`stat_qq_line` (`:542-557`); `pit` = histogram density + garis 1 (`:559-577`); `cpo` = stem per domain (`:579-600`); `all` = kalibrasi dengan warna = SAE RSE (`:602-628`).
- Metode S3 lain di file utilitas: `autoplot.fastsae_benchmark` (`benchmark.R:664`), `plot/summary/print.fastsae_benchmark` (`:704/:641/:577`), `print/summary/autoplot/plot.fastsae_comparison` (`compare_sae.R:193/:209/:244/:238`), `print.fastsae_diagnose` (`diagnose.R:391`), `print.fastsae_sim_data` (`sim_area_data.R:204`), `print.fastsae_sim_series` (`sim_series_data.R:315`). Daftar registrasi resmi: `NAMESPACE:25-54`.

---

## 7. Simulasi, utils, data

### 7.1 `sim_spatial_weights(D=40, type=c("knn","grid","ring"), style=c("B","W"), k=4, coords=NULL, seed=NULL)` (`R/sim_spatial_weights.R:47-52`)
- `set.seed` hanya bila `seed` non-NULL (`:53-55`); `D>=2` (`:58-60`).
- **knn** (`:67-93`): koordinat `matrix(runif(2D), ncol=2)` bila tak diberikan (`:69`, dokumentasi `Uniform(0,1)^2`, `:25`); `k_eff <- min(max(1,k), D-1)` (`:81`); `dmat=dist(coords)`; tetangga terdekat kecuali diri (`:85-90`); **simetrisasi mutual** `A <- pmax(A, t(A))`, `diag(A)=0` (`:92-93`).
- **grid** (`:95-107`): `r=floor(sqrt(D))`, `c=ceiling(D/r)` (`:96-97`), `expand.grid` lalu ambil `D` baris pertama (`:99-100`), Rook via jarak Manhattan == 1 (`:104-106`).
- **ring** (`:109-121`): koordinat lingkaran `cos/sin` (`:110-111`), `i` terhubung `i-1` dan `i+1` dengan wraparound (`:115-120`).
- **style**: `W` → `A/rowSums` dengan baris-0 diset penyebut 1 (`:124-127`); `B` → matriks biner (`:129`).
- Return matriks `D×D` dengan `dimnames` "1".."D" dan atribut `coords`, `type`, `style` (`:132-135`).

### 7.2 `sim_area_data(D=42, W=NULL, spatial_type=c("knn","grid","ring"), rho=0.5, phi=0.6, sigma_u=0.5, beta=c(1.0,0.8,-0.5), n_unsampled=6, seed=NULL)` (`R/sim_area_data.R:47-55`)
DGP persis:
1. `W` dibuat `sim_spatial_weights(..., style="B")` bila NULL; `W_std = W/rowSums` (baris-0 → penyebut 1) (`:63-97`). Fallback koordinat: MDS `cmdscale` dari `1 - W/max(W+1e-6)` lalu `runif` (`:73-83`).
2. Koveriat: `x1 ~ N(2,1)`, `x2 ~ U(0,5)` (`:100-101`); `eta_fix = b0 + b1*x1 + b2*x2` (`:102`).
3. Parameter dipotong: `rho ∈ [-0.95,0.95]`, `phi ∈ [0,1]`, `sigma_u ≥ 0.01` (`:105-107`).
4. Efek spasial SAR: `v_raw = solve(I - rho*W_std, rnorm(D))`, lalu **distandardisasi** (mean 0, sd 1) (`:109-113`); `u_iid` = normal distandardisasi (`:115-117`).
5. **BYM2-style**: `u_total = sigma_u * (sqrt(phi)*v + sqrt(1-phi)*u_iid)` (`:119`); `eta = eta_fix + u_total` (`:120`).
6. Respon (semua bentuk persis, `:122-154`):
   - Gaussian: `vardir ~ U(0.04,0.16)`, `e~N(0,sqrt(vardir))`, `y_gaussian = eta + e` (`:124-126`).
   - Poisson: `exposure = round(U(80,400))`; `lambda_rate = exp(clip(eta/2 - 0.5, -4, 4))`; `y ~ Pois(exposure*lambda_rate)` (`:129-131`).
   - Binomial: `trials = round(U(50,250))`; `prob = plogis(eta-1)`; `y ~ Binom(trials, prob)` (`:134-136`).
   - Negatif Binomial: `theta=3`; `mu = exp(clip(eta/2,-3,4))*15`; `rate ~ Gamma(shape=3, scale=mu/3)`; `y ~ Pois(rate)` (`:139-142`).
   - Beta: `mu = clip(plogis(eta-1), 0.02, 0.98)`, `kappa=25`, `y ~ Beta(mu*25,(1-mu)*25)`, hasil dipotong ke `[0.001,0.999]` (`:145-149`).
   - Gamma: `shape=5`, `mu = exp(clip(eta/2,-3,3))`, `y ~ Gamma(shape=5, scale=mu/5)` (`:152-154`).
7. Unsampled: `n_unsampled` dipotong ke `[0, D-1]`, `sample()` indeks, keenam respon diset `NA` (`:157-166`).
8. Return kelas `c("fastsae_sim_data","list")` berisi `data` (15 kolom, termasuk `truth_gaussian = eta` dan kolom **`domain`**, `:168-185`), `W`, `W_std`, `coords`, `u_spatial=v`, `u_iid`, `u_total`, `phi`, `rho` (`:187-200`).

### 7.3 `sim_series_data(D=42, T=5, time_start=2022, W=NULL, spatial_type=..., rho_s=0.5, rho_t=0.6, sigma_s=0.4, sigma_t=0.3, trend=c("linear","random_walk","ar1","none"), trend_slope=0.05, beta=c(1.0,0.8,-0.5), n_unsampled=4, prop_intermittent=0.05, sort_order=c("domain-major","time-major"), seed=NULL)` (`R/sim_series_data.R:66-81`)
1. `set.seed` bila `seed` (`:82-84`); `T>=2` (`:86-87`); `time_vec = time_start + (1:T) - 1` (`:89-90`).
2. `W`/`W_std` seperti §7.1–7.2 (`:97-131`).
3. Efek spasial `u1`: SAR `solve(I - rho_s*W_std, rnorm(D))`, distandardisasi, **dikalikan `sigma_s`** (`:135-141`).
4. Efek spasial-waktu `u2` (AR(1) per domain, matriks D×T): `u2[d,1] ~ N(0, sigma_t)`; `u2[d,t] = rho_t*u2[d,t-1] + N(0, innov_sd)` dengan `innov_sd = sigma_t*sqrt(max(0.01, 1-rho_t^2))` → varian stasioner = `sigma_t^2` (`:144-154`).
5. Tren deterministik (`:157-167`): `linear = trend_slope*(t-1)`; `random_walk = cumsum(c(0, rnorm(T-1, sd=trend_slope)))`; `ar1`: `rw[t]=0.7*rw[t-1]+rnorm(sd=trend_slope)`; `none = 0`.
6. Koveriat dinamis: `x1_base~N(2,0.8)`, `x2_base~U(1,4)` (`:171-172`); `x1 = x1_base + N(0,0.25)`, `x2 = x2_base + U(-0.4,0.4)` (`:191-192`).
7. `eta = b0+b1*x1+b2*x2+trend_vec[t] + u1[d] + u2[d,t]` (`:195-198`).
8. Respon sama persis formula §7.2 di N = D·T baris, dengan pengecualian: `vardir=round(U(0.04,0.16),5)` (`:202`), `exposure = round(U(90,420)*(1+0.02*(t-1)))` (`:207`), `trials = round(U(60,280))` (`:212`), NB/Beta/Gamma identik (`:217-232`).
9. Pengurutan baris (`:174-181`): `domain-major` → `expand.grid(time_idx, d_idx)` (time tercepat → semua periode utk domain 1 dulu, sesuai syarat `eblup_stfh`); `time-major` → `expand.grid(d_idx, time_idx)`.
10. Missingness: `n_unsampled` dipotong `[0,D-1]`, domain persisten diset `NA` di semua respon (`:236-248`); `prop_intermittent` dipotong `[0,0.3]`, `round(n*prop)` baris acak dari area ter-sampel diset `NA` (`:251-264`).
11. Return kelas `c("fastsae_sim_series","list")`: `data` (15 kolom: `area, year, x1, x2, y_gaussian, vardir, y_poisson, exposure, y_binomial, trials, y_beta, y_nbinomial, y_gamma, x_coord, y_coord`, `:270-287`), `W`, `W_std`, `coords`, `u_spatial=u1`, `u_temporal=u2`, `parameters` (D,T,rho_s,rho_t,sigma_s,sigma_t,beta,trend,trend_slope,sort_order) (`:289-308`).

### 7.4 `R/utils.R` (38 baris)
- `utils::globalVariables(c("density",".data"))` (`:5`).
- `.get_variable(data, variable)` (`:19-38`): menerima string kolom (`:20-25`), formula satu-sisi (`:26-32`), atau vektor panjang-n (`:33-34`); selain itu `cli_abort` (`:24`,`:31`,`:36`). `@noRd` (tidak diekspor).

### 7.5 `R/data.R` (dokumentasi dataset)
- `mys` (`:1-18`): 42 baris × 9 var, 10 domain non-sampel; kolom `area,y,vardir,rse,x1,x2,x3,n,weight`; source = simulasi berbasis Susenas/BPS (`:17`).
- `mys_panel` (`:20-38`): 462 baris × 10 var, 42 wilayah × 11 tahun (2016–2026).
- `mys_proxmat` (`:40-44`): matriks 42×42 **row-standardized**, source "Simulated based on contiguous administrative boundaries".
- `cornsoybean` (`:47-75`) & `cornsoybeanmeans` (`:78-106`): dataset asli dari paket `sae`; `@source` Battese, Harter & Fuller (1988), *JASA* 83, 28–36 (`:68-71`, `:100-103`).
- `sim_area` (`:109-145`): 42×14, kolom `area,...` (tanpa `truth_gaussian`); source `sim_area_data` + struktur `mys_proxmat` (`:138-139`).
- `sim_panel` (`:148-187`): 210×15, 2022–2026, domain-major; source `sim_series_data` + `mys_proxmat` (`:180-181`).
- Verifikasi terhadap `.rda` yang dikirim: `data/sim_area.rda` = 42×14 dengan kolom `area` (TANPA `truth_gaussian`), `data/sim_panel.rda` = 210×15 dengan `area,year,...` — lihat D-11.

---

## 8. `@references` / DOI
- **Tidak ada satu pun DOI atau URL di blok `@references` seluruh `R/*.R`** (grep `doi|DOI|http` hanya menemukan URL repositori INLA di `R/inla_utils.R:12`).
- Blok `@references` hanya ada di 10 file; yang termasuk lingkup catatan ini: `benchmark.R:79-93` (5 entri, disalin persis di §1.5) dan `diagnose.R:27-33` (3 entri, disalin persis di §2.8). File lain dalam lingkup (compare/export/map/methods/autoplot/sim/utils/data) **tidak punya `@references`**; `data.R` hanya punya `@source` (`:17`,`:37`,`:43`,`:68-71`,`:100-103`,`:138-139`,`:180-181`).

---

## 9. DISCREPANCIES vs rujukan publik / vs dokumentasi internal

| ID | Isu | Lokasi kode |
|---|---|---|
| D-1 | Tes Brown: regresi memakai **OLS tanpa bobot** (`lm(y_sub ~ pred_sub)`); formulasi Brown et al. (2001) yang lazim dikutip memakai WLS dengan bobot ∝ 1/(v_d + mse_d). Status: **TIDAK PASTI dari repo ini** (tidak ada teks rujukan di repo), tapi berbeda dari rumus standar yang umum dikutip. | `R/diagnose.R:144` |
| D-2 | Statistik Wald dibagi 2 dan dibandingkan dengan **F(2, n−2)**, bukan χ²(2) asimptotik klasik dari `c'V⁻¹c`. Kedua versi beredar; angka p akan berbeda. | `R/diagnose.R:152-155` |
| D-3 | Keputusan GOF memakai penerimaan **interval dua sisi** `[qchisq(α/2), qchisq(1−α/2)]`, sedangkan `p_value` yang dicetak adalah uji ekor-atas tunggal. Akibatnya `is_good_fit` bisa FALSE sementara `p_value ≥ α` (W terlalu kecil), sehingga `status` boleh jadi CAUTION padahal p-value tampak baik. | `R/diagnose.R:176-180` vs `:338`, cetak `:447-449` |
| D-4 | Metrik yang dinamai "Mean Relative RMSE" sebenarnya `sqrt((pred−truth)^2)/|truth|` per domain lalu diambil rerata = **MAPE-like (mean absolute relative error)**, bukan akar kuadrat rerata kuadrat error relatif. | `R/diagnose.R:256`, `:267`, label cetak `:475` |
| D-5 | Coverage dicetak sebagai "95% CI Coverage" tanpa memverifikasi level CI (memakai `ci_lower/ci_upper` apa adanya); lolos bila ≥90%. | `R/diagnose.R:259-262`, `:477-481` |
| D-6 | Metode `"optimal"`: penyesuaian `δ_d = (mse_d / w_d)·λ` dengan `λ = diff/Σmse_d`. Solusi KKT standar dari min Σ δ_d²/mse_d s.t. Σ w_dδ_d = diff adalah `δ_d ∝ w_d·mse_d`. Constraint tetap terpenuhi (kode diverifikasi benar), tetapi **alokasi beda ketika bobot tidak seragam**; dengan bobot seragam keduanya sama-sama ∝ mse. Variabel juga dinamai `inv_prec` padahal isinya `mse`. | `R/benchmark.R:476-479` (dan varian tahap-1 `:413-416`) |
| D-7 | Tahap-1 `logit` selalu memakai bentuk `Σ W_share·plogis(...)` (bertipe mean) meski `type="total"` memakai `Σ target_map` untuk `nat_current_agg`/`diff_nat` → ketidakcocokan skala pada `type="total"`. | `R/benchmark.R:402` vs `:420` |
| D-8 | `group_orig` selalu dihitung dengan bobot ternormalisasi (skala mean) meski `type="total"`; saat `target=NULL` + hierarkis, target awal berskala mean sementara STAGE 2 `type="total"` mengagregasi tanpa normalisasi → skala target tidak konsisten. | `R/benchmark.R:333-334`, `:383-386` vs `:451-452` |
| D-9 | `add_reliability_flags` dokumen: "input object atau data frame dengan kolom `reliability`"; kode selalu mengembalikan `df` saja (penugasan ke `object$df_*` tidak efektif). | doc `R/export_sae.R:23` vs `:68-75` |
| D-10 | `export_sae` dokumen "hingga tiga sheet"; kode bisa menghasilkan **empat** (`Benchmarked_Groups`), dan mode CSV hanya menulis `Estimates` + `_summary.csv` (sheet Benchmarked tidak pernah ditulis di CSV). | doc `R/export_sae.R:89-99` vs `:153-156`, `:173-181` |
| D-11 | `sim_area_data()` menghasilkan kolom **`domain`** + **`truth_gaussian`** (15 kolom), sedangkan dataset terlampir `sim_area` & dokumentasinya (`data.R:120-136`) memakai **`area`** dan 14 kolom tanpa truth (diverifikasi langsung dari `data/sim_area.rda`). | `R/sim_area_data.R:169-185` vs `R/data.R:120-136` |
| D-12 | `map_sae` argumen `thresholds=c(20,30)` diteruskan ke `.generate_sae_map` tetapi **tidak pernah dipakai**; pemotongan reliability hard-coded 20/30 di extractor. Parameter `thresholds` efektif inert. | `R/map_sae.R:107`,`:169`,`:602` vs `:482-487` |
| D-13 | `diagnose()` bukan generic S3 (fungsi biasa) dan **tidak ada** `summary`/`plot` method untuk `fastsae_diagnose` (hanya `print` + `autoplot`); `@return` doc tidak menyebut komponen `bayesian_metrics` yang dikembalikan. | `R/diagnose.R:53`, `:16-25` vs `:380`; `NAMESPACE:43,47` |
| D-14 | `compare_sae` mendaftarkan **hanya** method `.default` (tidak ada `compare_sae.fastsae`); `thresholds` disimpan tapi tidak memengaruhi metrik apa pun. | `NAMESPACE:33`; `R/compare_sae.R:183`, pemakaian hanya `:343-348` |
| D-15 | `autoplot.list` memanggil `NextMethod()` bila elemen bukan fastsae — bisa gagal dispatch bila tidak ada method lanjutan; `.autoplot_multi_comparison` bergantung pada `names(x)` (list tanpa nama → data kosong). | `R/autoplot.R:68-71`, `:257` |
| D-16 | Pemetaan dua model (`map_sae` list/`model2`) menyusun `data.frame(matched_sf[, c("domain","geometry")], ...)` — bila kolom geometri `sf` bernama selain `"geometry"` (atau nama domain bukan `"domain"` setelah merge) ekspresi ini akan gagal. | `R/map_sae.R:217-220`, `:432-435` |

---

## 10. Pseudocode bernomor

### 10.1 `.benchmark_worker(domain, y_hat, target, weight, method, group, group_var, type, national_target, mse, call, model_family)` — `R/benchmark.R:280-573`
1. `n <- length(y_hat)`; abort bila ada `weight <= 0` (`:294-298`).
2. Susun `df_work` (domain, y_hat, weight, mse atau NA) dan kolom `group` (atau `"All"`) (`:301-312`).
3. `unique_groups`, `n_groups`; `is_hierarchical <- !is.null(national_target) && n_groups > 1` (`:314-318`).
4. Untuk tiap grup: `W_g = Σw`, `w_norm = w/W_g`, `group_orig = Σ w_norm·y_hat`, `group_mse = Σ w_norm²·max(mse,1e-8)` atau NA (`:328-341`).
5. Isi `target_map` (data.frame/named vector/vektor/skalar), atau `group_orig` bila hierarkis tanpa target, atau `national_target` bila 1 grup, atau abort (`:349-393`).
6. **STAGE 1** (bila hierarkis): `W_share = W_g/ΣW_g`; `agg = Σ W_share·target` (mean) atau `Σ target` (total); `diff = national_target − agg`; bila `|diff| > 1e-7` sesuaikan `target_map` dengan metode terpilih (ratio × konstanta; difference + diff atau diff/G; optimal + (mse_g/W_share)·λ; logit lewat `uniroot` pada [-50,50], tol 1e-9) (`:397-425`); simpan `stage1_df` (`:427-437`).
7. **STAGE 2** per grup: `w_calc = w/Σw` (mean) atau `w` (total); `agg = Σ w_calc·y`; hitung `y_bm` sesuai metode (lihat §1.2; logit lewat `uniroot` bracket ±20 → ±50, tol 1e-9) (`:444-512`).
8. Susun `res_df` (domain, group?, weight, original, benchmarked, adjustment, rel_adjustment_pct, target, national_target?) (`:515-531`).
9. Verifikasi per grup `|Σ w_calc·y_bm − T_g|` dan (bila hierarkis) verifikasi nasional (`:534-560`).
10. Set atribut (`method, type, hierarchical, group_var, stage1_summary, verification, national_verification, call`), kelas `c("fastsae_benchmark","data.frame")`, kembalikan (`:562-572`).

### 10.2 Tes Brown et al. (2001) — `R/diagnose.R:133-191`
1. Ambil area ter-sampel dengan `direct_y` dan `sae_pred` non-NA; bila jumlahnya `< 5`, tes dilewati (`:137-138`).
2. Regresi OLS `direct_y ~ sae_pred`; ambil `coefs` dan `vcov` (`:144-146`).
3. Wald: `diff = (α,β) − (0,1)`; `F = diff'·vcov⁻¹·diff / 2`; `df1=2`, `df2 = n−2`; `p = P(F_{2,n−2} > F)`; `is_unbiased = (p ≥ alpha_level)` (`:149-165`). Jika `vcov` tidak bisa diinversikan, lewati Wald (`:150-151`).
4. GOF: untuk area dengan `vardir+mse > 0` (≥ 3): `W = Σ (y − ŷ)²/(vardir + mse)`; `df = jumlah area valid`; `p = P(χ²_df > W)`; terima bila `qchisq(α/2,df) ≤ W ≤ qchisq(1−α/2,df)` (`:169-189`).

### 10.3 Moran's I pada residual — `R/diagnose.R:193-246`
1. Pilih matriks W: argumen → `object$W` → selesai (`:197-203`).
2. Bila `nrow(W) == N_total`: potong ke area ter-sampel `res = y − ŷ`, `W_sub = W[sampled, sampled]`; bila `n < 5` atau ada NA → selesai (`:207-212`).
3. `z = res − mean(res)`; `S0 = Σ W_sub`; bila `S0 ≤ 0` → selesai (`:214-216`).
4. `I = (n/S0) · [ Σ_i Σ_j w_ij z_i z_j / Σ z² ]` (`:217-219`).
5. `E[I] = −1/(n−1)`; `S1 = ½Σ_ij (w_ij+w_ji)²`; `S2 = Σ_i (Σ_j w_ij + Σ_j w_ji)²`; `K = n·Σz⁴/(Σz²)²` (`:222-225`).
6. `Var(I) = { n[(n²−3n+3)S1 − nS2 + 3S0²] − K[(n²−n)S1 − 2nS2 + 6S0²] } / [(n−1)(n−2)(n−3)S0²] − E[I]²`; batasi `Var ≥ 1e-8` (`:227-231`).
7. `z = (I − E[I])/√Var(I)`; `p = 2·Φ(−|z|)`; `no_residual_autocorrelation = (p ≥ alpha_level)` (`:232-242`).

---
*Status: semua angka/rumus diambil langsung dari kode; item yang tidak dapat dipastikan dari repo ditandai TIDAK PASTI (terutama D-1).*
