# Catatan teknis: model frekuentis fastsae (FH, BHF, two-fold)

Sumber (semua baris = baris aktual; tidak ada file di luar `paper/_notes/` yang disentuh):
`R/eblup_fh.R`, `src/eblup_fh.cpp`, `R/eblup_bhf.R`, `src/eblup_unit.cpp`,
`R/eblup_twofold.R`, `src/eblup_tfh.cpp`, `R/RcppExports.R:4-33`, `R/utils.R:19-38`.
Pembanding `sae` 1.3 (`mseFH`, `eblupBHF`, `pbmseBHF`) dibaca via `Rscript` (read-only).
Helper `.get_variable(data, variable)` (`R/utils.R:19-38`): string kolom / formula satu sisi /
vektor apa adanya; error bila tidak cocok.

## 1. Fay-Herriot area-level: `eblup_fh()` + `.eblup_core`

### 1.1 Signature & return
- `R/eblup_fh.R:49-58`
  ```r
  eblup_fh(formula, vardir, domain = NULL, data, method = c("REML","ML"),
           maxiter = 100, precision = 1e-4, print_result = TRUE)
  ```
  `method` di-`match.arg` baris 59 (default **REML**); `domain` default `1:nrow(data)` (60-64);
  `vardir`/`domain` bisa nama kolom/formula/vektor (63, 68); `model.frame(na.pass)` (67) agar
  baris `y=NA` tetap ada; `X` dari `model.matrix` (75); NA pada X = error (78-80);
  panjang `vardir` = nrow(mf) (70-72).
- Panggilan C++: `.eblup_core(Xall, yall, vardirall, method, maxiter, precision)` di
  `R/eblup_fh.R:82-89` (wrapper `R/RcppExports.R:4-5`; fungsi `eblup_core`
  `src/eblup_fh.cpp:8-15`, export `[[Rcpp::export(.eblup_core)]]` baris 8).
- Return (list, **kelas S3 `"fastsae"`**, `R/eblup_fh.R:102`): elemen C++
  (`src/eblup_fh.cpp:252-261`) `random_effect_var, estcoef, df_eblup, goodness, n_iter,
  convergence, method="eblup", level="area"`; plus R: `formula` (93), `model="FH"` (94),
  `call` (100), `data` (101). `df_eblup` diurutkan (97-98) → `domain, y, eblup, vardir,
  random_effect, mse, rse`. `estcoef` = `beta, std.error, stderr_beta, zvalue, pvalue`
  (`src/eblup_fh.cpp:235-241`; dua kolom SE identik), rownames = kolom X (`R/eblup_fh.R:92`).
  `goodness` = `loglikelihood, AIC, BIC` (`src/eblup_fh.cpp:176-180`).

### 1.2 Model sebagaimana dikodekan
- $y_d = x_d'\beta + u_d + e_d$, $u_d\sim N(0,\sigma^2)$, $e_d\sim N(0,\psi_d)$ (psi = `vardir`),
  prior flat untuk $\beta$ (implisit via GLS); tidak ada blok komentar model di file ini, yang
  ada $V_d=\psi_d+\sigma^2$, $w_d=1/(\psi_d+\sigma^2)$ (`src/eblup_fh.cpp:89`).
- $\hat\beta = (X'V^{-1}X)^{-1}X'V^{-1}y$ (135-143), $B_d^{code} = \psi_d/(\psi_d+\hat\sigma^2)$
  (161-162), $\hat u_d = \frac{\hat\sigma^2}{\hat\sigma^2+\psi_d}(y_d-x_d'\hat\beta)$ (164),
  $\hat\theta_d = x_d'\hat\beta+\hat u_d$ (165).
- **Domain unsampled** = baris `y = NA` (`find_nonfinite`, 33-47): hanya `Xns` disimpan,
  estimasi dari sampel saja (40-42); prediksi $x_{ns}'\hat\beta$ (208),
  $\mathrm{MSE}=\hat\sigma^2+x_{ns}'Qx_{ns}$ (209), `random_effect=NA` (216).
- `vardir` harus **> 0 untuk area tersampel** (57-59) else `stop`.
- Loglik yang dilaporkan = loglik **ML** $-\frac12\sum_d[\log(2\pi(\hat\sigma^2+\psi_d))+r_d^2/(\hat\sigma^2+\psi_d)]$
  (170-172), AIC $=-2\ell+2(p+1)$, BIC $=-2\ell+(p+1)\log m$ ($m$ = area tersampel, 173-174)
  — **tetap ML walau `method="REML"`**.

### 1.3 Estimasi varians (Fisher scoring)
- **Starting**: $\hat\sigma^{2(0)}=\mathrm{median}(\psi_d)$, dipangkas ≥0 (64-65).
- Loop `while ((diff > precision) && (k < maxiter))` (**88**); default `precision=1e-4`,
  `maxiter=100` (C++ 14-15, R 55-56).
- Per iterasi: $W=\mathrm{diag}(w_d)$ (89), $X'WX$ (92-93), $Q=(X'WX)^{-1}$ via `inv_sympd`
  + fallback `solve` (95-96), $Py = Wy - XQX'Wy$ (98-99).
- **Persamaan skor** (101-108): ML $s = -\frac12\sum_d w_d + \frac12(Py)'(Py)$ (103);
  REML $s = -\frac12\mathrm{tr}P + \frac12(Py)'(Py)$ dengan
  $\mathrm{tr}P = \sum_d w_d - \mathrm{tr}(Q\,X'W^2X)$ (105-107).
- **Informasi Fisher** (110-118): ML $I=\frac12\sum_d w_d^2$ (112); REML $I=\frac12\mathrm{tr}(P^2)$
  dengan $\mathrm{tr}(P^2)=\sum w^2 - 2\sum_d w_d^2(x_d'Qx_d)w_d + \mathrm{tr}(K^2)$,
  $K=Q\,X'W^2X$ (114-117; `M_diag` $=w_d^2 x_d'Qx_d$).
- **Update** $\hat\sigma^{2\leftarrow}=\hat\sigma^2+s/I$ (122); **batas**: $I\le0$ →
  `numeric_limits::min()` (120), update non-finite/negatif → **dipaksa 0** (123, lagi di 130).
- **Konvergensi** (125): $\texttt{diff}=|(\sigma^{2(new)}-\sigma^2)/\max(\sigma^2,10^{-12})|$;
  loop berhenti bila `diff<=precision` (88); `convergence=(diff<=precision)` (258),
  `n_iter=k` (257); R memperingati bila gagal (`R/eblup_fh.R:104-107`).
- Pass akhir dengan $\sigma^2$ final (135-144): $\beta$, `se=sqrt(diag(Q))`, z, p dua sisi (149-154).

### 1.4 MSE (analitik; `src/eblup_fh.cpp:183-202`)
Dengan $w_i=1/(\psi_i+\hat\sigma^2)$ (135), $\sum_d w_d^2$ (185),
$\mathrm{Var}A=2/\sum_d w_d^2$ (186), $h_d=x_d'Qx_d$ (189):
- $g_1 = \psi_d(1-B_d)$ (191); $g_2 = B_d^2 h_d$ (192); $g_3 = B_d^2\,\mathrm{Var}A/(\hat\sigma^2+\psi_d)$ (193).
- REML: $\mathrm{MSE}=g_1+g_2+2g_3$ (201); ML: $\mathrm{MSE}=g_1+g_2+2g_3-b\,B_d^2$ (199) dengan
  $b=-\mathrm{tr}(Q\,X'W^2X)/\max(\sum w_d^2,10^{-30})$ (197-198).
- Area unsampled: $\hat\sigma^2+x'Qx$ (209) — tanpa $g_3$.
- $\mathrm{RSE}(\%)=100\sqrt{\mathrm{MSE}}/|\hat\theta|$ (227-230), NaN bila $\hat\theta=0$ (229).
- **Tidak ada bootstrap** di FH; hanya rumus di atas.

### 1.5 OpenMP
**Tidak ada.** `src/eblup_fh.cpp` tanpa `#include <omp.h>`/`#pragma omp`. `#pragma omp` di
paket hanya di `eblup_tfh.cpp:481`, `eblup_sfh_pbmse.cpp:108`, `eblup_sfh_npbmse.cpp:159`,
`eblup_stfh_pbmse.cpp:352`.

### 1.6 `@references` (disalin persis, `R/eblup_fh.R:5-8`; tidak ada `@doi`)
```
#' @references
#' \enumerate{
#'  \item Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley & Sons.
#' }
```

### 1.7 DISCREPANCIES (FH)
1. **Parametrisasi $B_d$**: kode memakai $B_d=\psi_d/(\psi_d+\sigma^2)$ (162) = komplemen
   $B_d^{RM}=\sigma^2/(\sigma^2+\psi_d)$ yang lazim di Rao–Molina; identik dengan
   `sae::mseFH` (`Bd <- vardir/(A+vardir)`).
2. $g_1$ ditulis $\psi_d(1-B_d)=\sigma^2\psi_d/(\sigma^2+\psi_d)$ (191) — ekuatip dengan
   $(B_d^{RM})^2\psi_d + (1-B_d^{RM})^2\sigma^2$; bentuk $(1-B_d)^2(v_d+\sigma^2)$ yang
   kadang muncul di literatur **bukan** bentuk yang dikodekan.
3. `goodness` selalu memakai loglik ML (170-174) walau `method="REML"`; AIC/BIC memakai
   $m$ = area tersampel. Bandingkan two-fold yang menghitung loglik REML (bagian 3).
4. Koreksi bias $-bB_d^2$ **hanya untuk ML** (196-199); REML tanpa koreksi orde-2 (201).
   (Ini juga yang dilakukan `sae::mseFH`.)
5. `estcoef` memuat dua kolom identik `std.error` dan `stderr_beta` (237-238).
6. `vardir` harus `> 0` ketat (57-59); `sae::mseFH` tidak memeriksa ini.
7. Area unsampled didefinisikan **hanya** dari `y = NA` (33), bukan dari indeks area; tidak
   ada argumen `domain` yang dipakai untuk pemisahan sampel/tak-sampel (`domain` hanya label
   di `R/eblup_fh.R:97`).

### 1.8 Pseudocode `eblup_fh` / `.eblup_core`
1. `R/eblup_fh.R:59-80` — `match.arg(method)`, resolve `domain`, `model.frame(na.pass)`,
   ambil `y`, `X`, `vardir`; validasi panjang & `anyNA(X)`.
2. `.eblup_core(...)` dipanggil (`R/eblup_fh.R:82`).
3. `src/eblup_fh.cpp:33-47` — pisahkan baris `y` finite (sampel) / nonfinite (unsampled);
   simpan `Xns`. Baris 53-59 — cek `m>0`, method, `vardir>0`.
4. Baris 64-65 — $\sigma^{2(0)} = \max(\mathrm{median}(\psi),0)$.
5. Baris 88-128 — Fisher scoring: $W,Q,Py$; skor $s$; informasi $I$;
   $\sigma^{2(new)}=\sigma^2+s/I$; clamp ke 0 bila <0/non-finite (123); `diff` relatif (125);
   `k++`; stop bila `diff<=precision` atau `k>=maxiter`.
6. Baris 135-154 — rehitung $Q,\hat\beta$, SE, z, p dengan $\sigma^2$ final.
7. Baris 159-180 — $\hat\theta_d = x_d'\hat\beta+\hat u_d$; loglik ML, AIC, BIC.
8. Baris 183-202 — MSE $g_1,g_2,g_3$; REML → $g_1+g_2+2g_3$; ML → $-bB_d^2$.
9. Baris 207-222 — cabang unsampled: $x'\hat\beta$, $\mathrm{MSE}=\hat\sigma^2+x'Qx$,
   `random_effect=NA`; gabung ke vektor penuh.
10. Baris 227-261 — RSE, tabel koefisien, `df_eblup`, list output.
11. `R/eblup_fh.R:92-111` — metadata, sisipkan `domain`, `class<-"fastsae"`, cek
    `convergence`, `print()`.

## 2. BHF unit-level: `eblup_bhf()` + `.eblup_bhf_cpp` + `.pbmse_unit`

### 2.1 Signature & return
- `R/eblup_bhf.R:44-57`
  ```r
  eblup_bhf(formula, unit_data, Xpop, domain_var, popsize_var,
            method = c("REML","ML"), popnmean_xpop = NULL, B = 100,
            compute_mse = FALSE, n_threads = 1, seed = -1, print_result = TRUE)
  ```
- `model.frame(na.action=na.omit)` (61); `dom` dipotong mengikuti baris yang di-omit (63-65);
  `selectdom = unique(dom)` lalu **digabung dengan domain `Xpop`**
  (`unique(c(selectdom, pop_dom))`, 88) sehingga domain populasi tanpa sampel ikut masuk;
  `Xpop` diurutkan sebaris dengan `selectdom` (89).
- **Variance components & $\beta$ diestimasi di R dengan `lme4::lmer`** (`y ~ x1+x2+...+(1|dom)`,
  dibangun dari `term.labels`, 73-82, `REML=(method=="REML")`): $\hat\sigma_u^2$ =
  `VarCorr(fit)$dom[1,1]` (107), $\hat\sigma_e^2$ = `attr(VarCorr,"sc")^2` (108),
  $\hat\beta$ = `fixef` (83), $\hat u_d$ = `ranef(fit)$dom` (84).
- `meanxpop`: `popnmean_xpop=NULL` → `model.matrix(~x1+x2, Xpop)` (92-94); bila diberikan,
  intercept disisipkan bila kolom kurang (97-99); `ncol` wajib = `ncol(Xs)` (102-104).
- Panggilan C++ `.eblup_bhf_cpp(selectdom, dom, Xs, meanxpop, ys, popnsize, betaest, upred)`
  di `R/eblup_bhf.R:111-120` (wrapper `R/RcppExports.R:32-33`; fungsi `eblup_bhf_cpp`
  `src/eblup_unit.cpp:8-18`).
- Return (list, **kelas `c("fastsae","fastsae_unit")`**, 205): `df_eblup`, `eblup`
  (duplikat `df_eblup`, 186), `estcoef`, `random_effect_var`,
  `fit = list(method, random_effect_var, sigma2_e, beta, random_effect, lme)`, `formula`,
  `model="BHF"`, `level="unit"`, `convergence = TRUE` (200), `data`, `call`. **Tidak ada
  `goodness`/`n_iter`.** `df_eblup` (160-168): `domain, eblup, samp_size` + `mse, rse`
  (atau `mse=rse=NA`). `estcoef` (171-181): `beta, std.error, stderr_beta, tvalue, zvalue,
  pvalue` dengan `zvalue = tvalue` (t-value lmer) dan `pvalue = 2*pnorm(|t|)`.

### 2.2 Prediktor (C++ `src/eblup_unit.cpp:57-101`)
Untuk domain $d$ dengan $m_d$ unit sampel, $f_d = m_d/N_d$ (baris 64), $\bar x_{s,d}$ (67-72),
$\bar y_{s,d}$ (75-77), $\bar x_{pop,d}$ = baris `meanxpop[i,]`:
$$\hat\theta_d = f_d\,\bar y_{s,d} + (\bar x_{pop,d} - f_d\,\bar x_{s,d})'\hat\beta + (1-f_d)\hat u_d$$
(baris 80-83, 91). $\hat u_d$ diambil dari map `upred` berbasis rowname data.frame (38-45).
Domain **tanpa unit sampel** (94-100): prediktor sintetis $\bar x_{pop,d}'\hat\beta$, `samp_size=0`,
didorong ke `warn_domains` (99) lalu `cli_warn` di R (128-132).
Validasi: `ncol(meanxpop)==p`, `nrow(meanxpop)==I`, `length(popnsize)==I` (23-31) — jika domain
sampel tidak ada di `Xpop` → **error** (bukan NA).
Seluruh perhitungan berat/gabungan ada di C++; tidak ada iterasi di C++ (beta/v sudah dari lmer).

### 2.3 Bootstrap MSE parametrik `.pbmse_unit` (`R/eblup_bhf.R:215-369`)
- Signature (215-225): `(.pbmse_unit(formula, unit_data, Xpop, domain_var, popsize_var,
  method=c("REML","ML"), B=100, n_threads=1, seed=-1))`, `@noRd` (214).
- Persiapan identik dengan `eblup_bhf` (229-267) → refit lmer, ambil $\hat\beta,\hat u,
  \hat\sigma_u^2,\hat\sigma_e^2$.
- EBLUP awal (270-282) — dipakai untuk `eblup` di return, **tidak dipakai** oleh `eblup_bhf`.
- Seed: `if (seed >= 0) set.seed(seed)` (285) — default `-1` = tanpa seed.
- **Jenis bootstrap: parametrik.** Tiap replikat $b$ (loop `for (b in seq_len(B))`, 293):
  1. $u^*_{d}\sim N(0,\hat\sigma_u^2)$ semua domain (295), $e^*_{j}\sim N(0,\hat\sigma_e^2)$
     semua unit sampel (298); $y^*_{j} = x_j'\hat\beta + u^*_{d(j)} + e^*_j$ (301-305).
  2. Refit `lmer` (308-316); **gagal → `next`** (317); lalu $\hat\beta^*,\hat u^*$ (319-320) dan
     EBLUP bootstrap via `.eblup_bhf_cpp` dengan `ys=y_boot` (323-332).
  3. "Nilai kebenaran" populasi (335-354): $\mu_d = \bar x_{pop,d}'\hat\beta$ (342);
     bila $n_d>0,\ r_d=N_d-n_d>0$: $\bar e^*_{s,d}n_d/N_d + \mathcal N(0,\hat\sigma_e^2/r_d)\,r_d/N_d$
     (343-347); bila $r_d=0$ (full enumeration): hanya $\bar e^*_{s,d}$ (348-350);
     bila $n_d=0$: $\mu_d + u_d^* + \mathcal N(0,\hat\sigma_e^2/N_d)$ (351-352).
  4. Akumulasi $(\hat\theta^*_d-\theta^{*,true}_d)^2$ (356-357).
- $\mathrm{MSE}_d = \frac1B\sum_b(\cdot)^2$ (**dibagi `B` di 360**, termasuk bila ada replikat
  gagal); return `list(eblup, mse, method, B)` (363-368). Digabung dengan `merge(..., sort=FALSE)`
  (157); `rse = 100*sqrt(mse)/|eblup|` dengan guard `|eblup| < .Machine$double.eps → NA` (161-164).
- **Paralelisme: tidak ada.** `n_threads` diterima (54, 148, 223) tapi **tidak pernah dipakai**
  (tak ada `foreach`/`parallel`/OpenMP) → bootstrap = loop `for` serial di R.

### 2.4 OpenMP
Tidak ada di `src/eblup_unit.cpp` (tanpa `omp.h`/`#pragma omp`).

### 2.5 `@references` (disalin persis, `R/eblup_bhf.R:22-25`; tidak ada `@doi`)
```
#' @references
#' Battese, G. E., Harter, R. M., and Fuller, W. A. (1988). An error-components
#' model for prediction of county crop areas using survey and satellite data.
#' *Journal of the American Statistical Association*, 83(401), 28-36.
```

### 2.6 DISCREPANCIES (BHF)
1. **`n_threads` didokumentasikan tapi tak dipakai** (`@param` di `R/eblup_bhf.R:16`;
   argumen 54 → 148 → 223 tanpa pemakaian) → MSE bootstrap selalu serial.
2. **Default `B = 100`** (52, 222) vs `sae::pbmseBHF` default `B = 200`.
3. **`seed = -1` default = tanpa seed** (55, 285); `sae::pbmseBHF` tidak punya argumen seed.
4. **Replikat gagal tetap ikut dibagi `B`** (`next` 317 + `mse <- mse/B` 360) → MSE ter-bias
   ke bawah saat `lmer` gagal; `sae::pbmseBHF` hanya menaikkan `b` bila refit sukses.
5. **`popnmean_xpop` hanya untuk estimasi titik** (92-100) dan **tidak diteruskan** ke
   `.pbmse_unit` (140-150); bootstrap selalu membangun `meanxpop` dari `Xpop` (246-250) →
   titik estimasi dan MSE bisa memakai matriks rata-rata populasi berbeda.
6. Domain sampel **tidak ada di `Xpop`** → `stop` dari C++ (`src/eblup_unit.cpp:26-28`);
   `sae::eblupBHF` memberi `NA` + warning.
7. `convergence = TRUE` di-hardcode (200); kegagalan `lmer` melempar error (82, tanpa
   `tryCatch`) — tidak ada status konvergensi di model ini.
8. `pvalue` dibuat dari t-value lmer dengan distribusi normal (179); `lmer` sendiri tidak
   menyediakan p-value.
9. Prediktor = $f_d\bar y_s + (\bar x_{pop}-f_d\bar x_s)'\hat\beta + (1-f_d)\hat u_d$
   (`src/eblup_unit.cpp:91`) — **identik dengan `sae::eblupBHF`**; ekuatip dengan
   $f_d\bar y_s + (1-f_d)(\bar x_{unsamp,d}'\hat\beta+\hat u_d)$,
   $\bar x_{unsamp}=(\bar x_{pop}-f_d\bar x_s)/(1-f_d)$, **bukan** bentuk alternatif
   $f_d\bar y + (1-f_d)(\bar x_{pop}'\hat\beta+\hat u_d)$ (beda $f_d(\bar x_{pop}-\bar x_s)'\hat\beta$).
10. `Xpop` diurutkan dua kali (89 dan 243) — idempoten, tapi urutan bergantung pada
    `unique(c(selectdom,pop_dom))`.
11. Tidak ada AIC/BIC/loglik di output BHF (berbeda dengan FH & two-fold).

### 2.7 Pseudocode `eblup_bhf` (1-10) dan `.pbmse_unit` (b-1..b-6)
1. `44-66` — `match.arg(method)`; `model.frame(na.omit)`; `dom` dipotong untuk baris NA;
   `selectdom = unique(dom)`.
2. `69-70` — `ys`, `Xs = model.matrix`. 3. `73-84` — formula `y ~ terms + (1|dom)`,
   `lme4::lmer(REML=...)`, `fixef`, `ranef`; `107-108` → $\hat\sigma_u^2,\hat\sigma_e^2$.
4. `87-90` — `selectdom = unique(c(selectdom, domains(Xpop)))`; urutkan `Xpop`; `popnsize`.
5. `92-104` — bangun/validasi `meanxpop` (intercept disisipkan bila perlu).
6. `111-120` → `.eblup_bhf_cpp`: loop per domain `src/eblup_unit.cpp:57-101`: $f_d$,
   $\bar x_{s,d}$, $\bar y_{s,d}$, $\hat u_d$ → prediktor; cabang tanpa-sampel → sintetis +
   `warn_domains`.
7. `122-132` — `data.frame(domain, eblup, samp_size)`; warning `warn_domains`.
8. `135-158` — bila `compute_mse`: `.pbmse_unit(...)` lalu `merge` MSE per domain.
9. `160-181` — kolom `mse`/`rse` (atau NA); `estcoef` dari `summary(fit)$coefficients`.
10. `184-211` — rakit output, `class = c("fastsae","fastsae_unit")`, `print()`.

b-1. `229-250` — model frame, `selectdom`, urutkan `Xpop`, `popnsize`, `meanxpop`
   (intercept disisipkan bila kurang). b-2. `253-267` — refit lmer awal →
   $\hat\beta,\hat u,\hat\sigma_u^2,\hat\sigma_e^2$. b-3. `270-282` — EBLUP awal (hanya untuk return).
b-4. `285-291` — `set.seed` bila `seed>=0`; `mse=numeric(I)`; indeks baris per domain.
b-5. `293-358` — untuk $b=1..B$: (a) gambar $u^*,e^*$; (b) $y^*=X\hat\beta+u^*+e^*$;
   (c) refit lmer, gagal → `next`; (d) `.eblup_bhf_cpp` dengan $y^*$;
   (e) `truemean_boot` (cabang $n_d>0,r_d>0$ / $r_d=0$ / $n_d=0$); (f) akumulasi selisih kuadrat.
b-6. `360-368` — `mse <- mse/B`; return `list(eblup, mse, method, B)`; lalu `eblup_bhf 157-168`
   `merge(..., sort=FALSE)` dan `rse = 100*sqrt(mse)/|eblup|`.

## 3. Two-fold subarea: `eblup_twofold()` (alias `eblup_tfh`) + `.eblup_tfh_core`

### 3.1 Signature & return
- `R/eblup_twofold.R:65-78`
  ```r
  eblup_twofold(formula, vardir, domain = NULL, subarea = NULL, data,
                method = c("REML","ML"), mse = c("analytical","bootstrap"),
                B = 200, seed = NULL, maxiter = 100, precision = 1e-4,
                print_result = TRUE)
  ```
  `method` default REML, `mse` default **analytical** (79-80).
  `domain` **wajib** (81-83); `subarea` default `seq_len(nrow(data))` (85-89).
  `area_idx <- as.integer(factor(domain)) - 1L` (107) → indeks **0-based** (level faktor =
  urut unik). `set.seed(seed)` bila `!is.null(seed)` (109).
- Panggilan C++ `R/eblup_twofold.R:111-121`: `.eblup_tfh_core(Xall, yall, vardirall, area,
  method, mse_type, B, maxiter, precision)` (wrapper `R/RcppExports.R:28-29`; fungsi
  `eblup_tfh_core` `src/eblup_tfh.cpp:273-277`).
- Return: list **kelas `"fastsae"`** (`R/eblup_twofold.R:137`); elemen C++
  (`src/eblup_tfh.cpp:550-553`) `random_effect_var` = vektor bernama `sigma2_v` (area),
  `sigma2_u` (subarea) (544), `estcoef`, `df_eblup`, `goodness` (`loglikelihood, AIC, BIC`;
  AIC $=-2\ell+2(p+2)$, BIC $=-2\ell+(p+2)\log m_{fit}$, 545-548), `n_iter`, `convergence`,
  `method="eblup"`, `level="subarea"`; plus `formula` (124), `model="TWOFOLD"` (126),
  `call` (135), `data` (136). `df_eblup` (129-133): `domain, subarea, y, eblup, vardir,
  random_effect_area, random_effect_subarea, mse, rse`.
- **Alias**: `eblup_tfh <- eblup_twofold` (`R/eblup_twofold.R:150-151`), keduanya diekspor
  (`NAMESPACE:8-9`), `@aliases eblup_tfh` (44) — objek fungsi yang sama.

### 3.2 Model (komentar `src/eblup_tfh.cpp:4-8`)
- Level 1 (sampling): $y_{ij} = \theta_{ij} + e_{ij}$, $e_{ij}\sim N(0,\psi_{ij})$ (`psi` = `vardir`).
- Level 2: $\theta_{ij} = x_{ij}'\beta + v_i + u_{ij}$, $v_i\sim N(0,\sigma_v^2)$ (area),
  $u_{ij}\sim N(0,\sigma_u^2)$ (subarea) — **notasi kertas**: `s2v` = area, `s2u` = subarea (baris 8, 63).
- Per area (blok): $V_d = \sigma_v^2\mathbf 1\mathbf 1' + \mathrm{diag}(\sigma_u^2+\psi_d)$
  (baris 377), $S_d = \sigma_v^2\mathbf 1\mathbf 1' + \sigma_u^2 I$ (378).
- Woodbury: `dd` $=1/(\sigma_u^2+\psi_d)$ (98/345), $s=\sum \text{dd}$, $c=\sigma_v^2/(1+\sigma_v^2 s)$
  (100/347, =0 bila $\sigma_v^2=0$), $V_d^{-1}a = \text{dd}\odot a - c\,\text{dd}(\text{dd}'a)$
  (105, 108, 350).
- Domain **unsampled** = `y = NA` (290-292); subarea tak-sampel di area tersampel:
  $\hat\theta = x'\hat\beta + \hat v_i$ dan $\hat u = 0$ (358-363); area tanpa subarea tersampel:
  default sintetis $\hat\theta = x'\hat\beta$ (333-334).
- Validasi (281-311): method/mse_type sah; `B>=1` bila bootstrap; panjang `area`/`y`/`vardir`;
  `area>=0`; `Ns>0`; $N_s > p$ (301); `vardir>0` (302); minimal **2 area tersampel** (311).

### 3.3 Estimasi varians — Fisher scoring (`tfh_pass` 80-169, `tfh_fit` 173-269)
- **Starting**: $\sigma_v^{2(0)}=\sigma_u^{2(0)}=0.5\,\mathrm{median}(\psi)$ (181-182);
  loop `while ((diff > precision) && (k < maxiter))` (**186**), `precision=1e-4`, `maxiter=100`.
- Per iterasi `tfh_pass` (80-169): akumul $X'V^{-1}X$ (110) & $X'V^{-1}y$ (111),
  $Q=(X'V^{-1}X)^{-1}$ (114-115), $\hat\beta=QX'V^{-1}y$ (116), $r=V^{-1}(y-X\hat\beta)$ (125),
  $q_d=\mathbf 1'V^{-1}\mathbf 1-(X'V^{-1}\mathbf 1)'Q(X'V^{-1}\mathbf 1)$ (128,134),
  $\mathrm{tr}P=\mathrm{tr}V-\mathrm{tr}(QX'V^{-2}X)$ (127,135), $\mathrm{tr}V^2$ (144), $\mathrm{tr}P^2$ (145).
- **Skor & informasi, orde $(v,u)$** — REML (131-151):
  $s_v=-\frac12q_d+\frac12(\mathbf 1'r)^2$, $s_u=-\frac12\mathrm{tr}P+\frac12r'r$,
  $I_{vv}=\frac12q_d^2$, $I_{uu}=\frac12\mathrm{tr}(P^2)$,
  $I_{vu}=\frac12[\mathbf 1'V^{-1}\mathbf 1-2\,\mathbf z'Q\overline x+\overline x'QMQ\overline x]$ (148-151).
  ML (152-158): $s_v=-\frac12\mathbf 1'V^{-1}\mathbf 1+\frac12(\mathbf 1'r)^2$,
  $s_u=-\frac12\mathrm{tr}V+\frac12r'r$, $I_{vv}=\frac12(\mathbf 1'V^{-1}\mathbf 1)^2$,
  $I_{uu}=\frac12\mathrm{tr}(V^{-2})$, $I_{vu}=\frac12\mathbf z'\mathbf z$, $\mathbf z=V^{-1}\mathbf 1$ (129).
- **Update** $\Delta=\mathcal I^{-1}s$ (193-199), fallback $\Delta_k=s_k/I_{kk}$ bila determinan/
  `inv_sympd` gagal (201-204); $\sigma^2\leftarrow\sigma^2+\Delta$ (206-207);
  **batas**: non-finite → nilai lama (208-209), negatif → **dipaksa 0** (210-211).
- **Konvergensi** (213-215): $\texttt{diff}=\max_k|(\sigma^{2(new)}_k-\sigma^2_k)/\max(\sigma^2_k,10^{-12})|$;
  `conv=(diff<=precision)` (230), `niter=k` (229).
- **Pass final** pada varian final (223) → $\beta$, $Q$, $\mathcal I^{-1}$ sesuai ML/REML
  (232-237); komentar revisi 26-40 menyebut ini perbaikan ($\beta/Q$ dulu dari iterate terakhir,
  informasi selalu ML).
- **Loglik** (239-267): core $=\sum_d[\log|V_d|+r_d'V_d^{-1}r_d]$ (240-253);
  ML $\ell=-\frac12(n\log2\pi+\text{core})$ (256-257); REML
  $\ell=-\frac12[(n-p)\log2\pi+\text{core}+\log|X'V^{-1}X|-\log|X'X|]$ (259-265).
- **EBLUP**: $\hat v_i=\sigma_v^2\mathbf 1'V^{-1}r$, $\hat u_{ij}=\sigma_u^2V^{-1}r$ (351-352),
  $\hat\theta=x\hat\beta+\hat v+\hat u$ (353).

### 3.4 MSE analitik (`src/eblup_tfh.cpp:367-469`)
Subarea **tersampel** (per area, 369-411), dengan $V=V_d$, $S=S_d$, $Q$:
- $g_1 = \mathrm{diag}(S - SV^{-1}S)$ (383-386).
- $g_2 = \mathrm{diag}(D Q D')$, $D=(I-SV^{-1})X_d$ (388-392).
- $g_3$: turunan prediktor (komentar 15-17):
  $M_u = V^{-1} - SV^{-2}$ (395), $M_v = \mathbf 1(\mathbf V^{-1}\mathbf 1)' - SV^{-1}\mathbf 1(\mathbf V^{-1}\mathbf 1)'$ (396-398),
  $\mathbb E[(\partial_{\sigma_u^2}\tilde\mu)^2] = \mathrm{diag}(M_uVM_u')$ (399-406),
  lalu $g_3 = \mathbb E_vv\,\mathrm{var}_v + \mathbb E_uu\,\mathrm{var}_u + 2\mathbb E_{vu}\mathrm{cov}_{vu}$ (407)
  dengan $\mathrm{var}_v,\mathrm{var}_u,\mathrm{cov}_{vu}$ dari $(\mathcal I^{-1})$ final (324).
- $\mathrm{MSE} = g_1+g_2+g_3$ (408); nilai negatif dalam $(-10^{-8},0)$ dipangkas ke 0 (409).
Subarea **tak-sampel** (415-468):
- Area **tanpa** subarea tersampel: $x'Qx + \sigma_v^2 + \sigma_u^2$ (421-427).
- Area tersampel: $g_1^* = \sigma_v^2 - \sigma_v^4(\mathbf 1'V^{-1}\mathbf 1) + \sigma_u^2$ (445);
  turunan $\mathbf a_v = \mathbf w - \sigma_v^2(\mathbf 1'V^{-1}\mathbf 1)\mathbf w$, $\mathbf a_u = -\sigma_v^2V^{-2}\mathbf 1$ (448-449);
  $g_3^* = \mathbb E_{vv}\mathrm{var}_v+\mathbb E_{uu}\mathrm{var}_u+2\mathbb E_{vu}\mathrm{cov}_{vu}$ (451-456);
  $g_2^* = (x - \sigma_v^2 X'V^{-1}\mathbf 1)'Q(x-\sigma_v^2 X'V^{-1}\mathbf 1)$ (461-463);
  $\mathrm{MSE}^* = g_1^*+g_2^*+g_3^*$ (464), clamp (465).

### 3.5 MSE bootstrap parametrik (`src/eblup_tfh.cpp:470-518`)
- Default `B = 200` (C++ 277 / R 73); seed di R (`R/eblup_twofold.R:109`).
- **Semua gambaran dibuat lebih dulu, serial**: $v^*\in\mathbb R^{m\times B}$ (475),
  $u^*\in\mathbb R^{N\times B}$ (476), $e^*\in\mathbb R^{N_s\times B}$ (477-478).
- **OpenMP**: `#pragma omp parallel for schedule(static)` **baris 481**, loop atas replikat
  $b = 0..B-1$ (482-516). Tiap thread menulis kolomnya sendiri dari `SqDiff` (479).
  Tidak ada `omp_set_num_threads` → jumlah thread lewat `OMP_NUM_THREADS` (argumen `seed`
  tidak memengaruhi jumlah thread; `eblup_twofold` tidak punya `n_threads`).
- Per replikat: $\theta^* = X\hat\beta + v^*_{i} + u^*_{j}$ (486-487), $y^* = \theta^*_{sampel}+e^*$ (488),
  refit `tfh_fit(..., want_loglik=false)` (490-491), prediktor bootstrap dengan prediktor yang
  **sama** seperti estimasi titik (termasuk cabang tak-sampel, 495-513),
  `(eb_star - theta_star)^2` (515); $\mathrm{MSE} = \mathrm{rata}(\cdot)$ per baris (517).
- Reprodukibilitas: semua RNG dipanggil **sebelum** pragma paralel → `set.seed` cukup,
  tidak bergantung pada jumlah thread.

### 3.6 OpenMP
Hanya `src/eblup_tfh.cpp:42-44` (`#ifdef _OPENMP` + `#include <omp.h>`) dan
`#pragma omp parallel for schedule(static)` di **baris 481**, loop replikat bootstrap (482-516).
Estimasi komponen varian, EBLUP, dan MSE analitik **serial**.

### 3.7 `@references` (disalin persis, `R/eblup_twofold.R:7-12`; tidak ada `@doi`)
```
#' @references
#' \enumerate{
#'  \item Torabi, M., & Rao, J. N. K. (2014). On small area estimation under a sub-area
#'    level model. \emph{Journal of Multivariate Analysis}, 127, 36--55.
#'  \item Rao, J. N. K., & Molina, I. (2015). \emph{Small Area Estimation}. John Wiley & Sons.
#' }
```

### 3.8 DISCREPANCIES (two-fold)
1. **Bukan rumus tertutup Torabi–Rao**: `src/eblup_tfh.cpp:19-24` menyatakan eksplisit bahwa
   eq. (3.4) Torabi & Rao (2014) "found to contain a spurious [A][B] term and a structurally
   incorrect var_u coefficient K (wrong shape… paper s4 also has wrong sign)"; pendekatan
   matriks untuk $g_3$ divalidasi dengan turunan beda-pemisahan (korelasi 0.97/0.96/0.74),
   koreksi $\sigma_u^{-2}$ "rejected as dimensionally wrong".
2. **MSE = $g_1+g_2+g_3$** (408), bukan $g_1+g_2+2g_3$; faktor 2 sudah ada di suku silang
   $2\mathbb E_{vu}\mathrm{cov}_{vu}$ (407) — jangan disamakan dengan penulisan PR untuk FH.
3. $\mathrm{var}_v,\mathrm{var}_u,\mathrm{cov}_{vu}$ dari **invers matriks informasi Fisher**
   (235, 324), bukan bentuk tertutup $\mathrm{Var}(\hat\sigma^2)$ seperti FH (`VarA=2/Σw²`,
   `eblup_fh.cpp:186`); g1 juga beda bentuk: $\psi_d(1-B_d)$ vs $\mathrm{diag}(S-SV^{-1}S)$.
4. Boundary: $\sigma_v^2,\sigma_u^2$ dipangkas ke 0 (210-211) tanpa proyeksi konstrain;
   fallback skor diagonal bila informasi singular (201-204).
5. Bootstrap menggambar $u^*$ untuk **semua** $N$ subarea (476) padahal tak-sampel tidak
   pernah masuk $y^*$; $e^*$ hanya untuk $N_s$ baris sampel (477).
6. MSE bootstrap memakai prediktor yang **sama** dengan estimasi titik, termasuk $\hat v_i$
   untuk subarea tak-sampel (komentar 509; komentar revisi 34-35 menyebut ini perbaikan bug).
7. `level="subarea"` (553) vs `"area"` (FH 260) vs `"unit"` (BHF); penamaan
   `random_effect_area`/`random_effect_subarea` mengikuti notasi kertas (543-544).
8. `set.seed(seed)` dieksekusi untuk **semua** jalur (109), bukan hanya `mse="bootstrap"`.
9. Negatif kecil hanya dipangkas bila $>-10^{-8}$ (409, 465); lebih negatif dibiarkan
   (RSE jadi `NaN` via `sqrt`, 521-525).
10. AIC/BIC memakai $m_{fit}$ = jumlah area **dengan subarea tersampel** (306-310, 548),
    berbeda basis dengan FH yang memakai $m$ = jumlah baris sampel (174).

### 3.9 Pseudocode `eblup_twofold` / `.eblup_tfh_core`
1. `R/eblup_twofold.R:79-107` — `match.arg(method, mse)`; wajib `domain`; resolve `subarea`;
   `model.frame(na.pass)`; cek panjang `vardir`; `anyNA(X)`; `area = factor(domain)-1`.
2. `109-121` — `set.seed(seed)` bila tidak NULL; panggil `.eblup_tfh_core(...)`.
3. `src/eblup_tfh.cpp:281-318` — validasi argumen/ukuran; `Ns<=p` → stop; `m_fit<2` → stop;
   pisahkan sampel/tak-sampel; `sidx[d]` (baris sampel per area), `nsg[d]` (baris tak-sampel).
4. `321` → `tfh_fit`: **(a)** start 181-182; **(b)** loop 186-219: `tfh_pass` →
   $\beta,Q$, skor $s$, informasi $\mathcal I$; $\Delta=\mathcal I^{-1}s$ (fallback diagonal
   201-204); clamp ≥0/non-finite (206-211); `diff` relatif maksimum (213-215); `k++`;
   **(c)** pass final 223-237 → $\beta,Q,\mathcal I^{-1}$; **(d)** loglik 239-267.
5. `327-364` — prediktor: default sintetis $x'\hat\beta$ untuk baris `NA` (333-334); per area
   tersampel $\hat v_i,\hat u_{ij},\hat\theta$ (351-353); subarea tak-sampel di area tersampel
   $x'\hat\beta+\hat v_i$, $\hat u=0$ (358-363).
6. `367-411` — **analitik**: per area bangun $V_d,S_d$; $g_1,g_2,g_3$; clamp kecil (409).
7. `415-468` — MSE subarea tak-sampel: area tanpa sampel → $x'Qx+\sigma_v^2+\sigma_u^2$;
   area bersampel → $g_1^*,g_2^*,g_3^*$.
8. `470-518` — **bootstrap**: gambar $v^*,u^*,e^*$ serial (475-478); `#pragma omp parallel
   for` (481) atas $b$: $\theta^*,y^*$, refit `tfh_fit(want_loglik=false)`, prediktor dengan
   cabang tak-sampel, selisih kuadrat; `mse_all = mean(SqDiff,1)` (517).
9. `521-553` — RSE (%), tabel koefisien (z/p normal dari $Q$), `df_eblup`,
   `random_effect_var` (`sigma2_v`,`sigma2_u`), `goodness`, list output.
10. `R/eblup_twofold.R:124-147` — metadata, `domain`/`subarea`, `class<-"fastsae"`,
    cek `convergence`, `print()`.
