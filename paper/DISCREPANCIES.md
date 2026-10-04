# DISCREPANCIES — perbedaan implementasi vs rumus/dokumentasi publik

Laporan teknis `paper/`. Setiap butir ditelusuri ke kode; nomor baris mengacu pada
repo `main` (commit `0252055`). Butir yang ditandai **TIDAK PASTI** tidak dapat
dipastikan dari repo (tidak ada rujukan di dalam repo) dan ditulis apa adanya.
Kode pendukung lengkap: `paper/_notes/*.md`.

Pembacaan: **kode adalah sumber kebenaran**; laporan selalu menggambarkan apa yang
dilakukan kode, lalu memberi *implementation note* yang menunjuk ke butir di sini.

---

## A. Fay–Herriot area-level (`eblup_fh`)

| ID | Butir | Lokasi |
|---|---|---|
| FH-1 | Parametrisasi $B_d$ memakai $B_d=\psi_d/(\psi_d+\sigma^2)$, yaitu komplemen dari $B_d=\sigma^2/(\sigma^2+\psi_d)$ yang lazim di Rao–Molina; identik dengan `sae::mseFH` (`Bd <- vardir/(A+vardir)`). | `src/eblup_fh.cpp:162` |
| FH-2 | $g_1=\psi_d(1-B_d)=\sigma^2\psi_d/(\sigma^2+\psi_d)$ — ekuatip dengan $(B_d^{RM})^2\psi_d+(1-B_d^{RM})^2\sigma^2$; bentuk $(1-B_d)^2(v_d+\sigma^2)$ yang kadang muncul di literatur **bukan** yang dikodekan. | `src/eblup_fh.cpp:191` |
| FH-3 | `goodness` selalu memakai loglik **ML** walau `method="REML"`; AIC/BIC memakai $m$ = jumlah area tersampel. (Two-fold berbeda: loglik REML dihitung benar.) | `src/eblup_fh.cpp:170-174` |
| FH-4 | Koreksi bias orde-kedua $-bB_d^2$ **hanya untuk ML**; REML tanpa koreksi — sama seperti `sae::mseFH`. | `src/eblup_fh.cpp:196-201` |
| FH-5 | `estcoef` memuat dua kolom identik `std.error` dan `stderr_beta`. | `src/eblup_fh.cpp:237-238` |
| FH-6 | `vardir` harus `> 0` ketat untuk area tersampel; `sae::mseFH` tidak memeriksa ini. | `src/eblup_fh.cpp:57-59` |
| FH-7 | Area *unsampled* didefinisikan **hanya** dari `y = NA`, bukan dari argumen `domain` (domain hanya label). | `src/eblup_fh.cpp:33`, `R/eblup_fh.R:97` |

## B. Battese–Harter–Fuller (`eblup_bhf`)

| ID | Butir | Lokasi |
|---|---|---|
| BHF-1 | `n_threads` didokumentasikan tetapi **tidak pernah dipakai** (diteruskan sampai ke `.pbmse_unit` tanpa pemakaian) → bootstrap BHF selalu serial. | `R/eblup_bhf.R:16,54,148,223` |
| BHF-2 | Default `B = 100`; `sae::pbmseBHF` memakai `B = 200`. | `R/eblup_bhf.R:52,222` |
| BHF-3 | Default `seed = -1` = tanpa seed; `sae::pbmseBHF` tidak punya argumen seed. | `R/eblup_bhf.R:55,285` |
| BHF-4 | Replikat bootstrap yang gagal tetap ikut membagi $B$ (`next` lalu `mse <- mse/B`) → MSE ter-bias ke bawah bila `lmer` gagal; `sae::pbmseBHF` hanya menaikkan penghitung bila refit sukses. | `R/eblup_bhf.R:317,360` |
| BHF-5 | `popnmean_xpop` hanya dipakai untuk estimasi titik dan **tidak diteruskan** ke `.pbmse_unit` → titik estimasi dan bootstrap bisa memakai matriks rata-rata populasi yang berbeda. | `R/eblup_bhf.R:92-100` vs `:140-150,246-250` |
| BHF-6 | Domain sampel yang tidak ada di `Xpop` → **error** dari C++; `sae::eblupBHF` memberi `NA` + warning. | `src/eblup_unit.cpp:26-28` |
| BHF-7 | `convergence = TRUE` di-hardcode; kegagalan `lmer` melempar error (tanpa `tryCatch`) — tidak ada status konvergensi. | `R/eblup_bhf.R:200,82` |
| BHF-8 | `pvalue` dibuat dari t-value `lmer` dengan distribusi normal; `lmer` sendiri tidak menyediakan p-value. | `R/eblup_bhf.R:179` |
| BHF-9 | Prediktor $f_d\bar y_s+(\bar x_{pop}-f_d\bar x_s)'\hat\beta+(1-f_d)\hat u_d$ **identik** dengan `sae::eblupBHF`; ekuatip dengan $f_d\bar y_s+(1-f_d)(\bar x_{unsamp}'\hat\beta+\hat u_d)$ — **bukan** bentuk alternatif $f_d\bar y+(1-f_d)(\bar x_{pop}'\hat\beta+\hat u_d)$. | `src/eblup_unit.cpp:91` |
| BHF-10 | Tidak ada AIC/BIC/loglik pada output BHF (berbeda dengan FH dan two-fold). | `R/eblup_bhf.R:160-211` |

## C. Two-fold subarea (`eblup_twofold`)

| ID | Butir | Lokasi |
|---|---|---|
| TFH-1 | **Bukan rumus tertutup Torabi–Rao.** Kode menyatakan eksplisit bahwa eq. (3.4) Torabi & Rao (2014) mengandung istilah $[A][B]$ yang semu dan koefisien varian $u$ yang salah bentuk (dan tanda salah di lampiran); pendekatan matriks dipakai untuk $g_3$, divalidasi dengan turunan beda pemisah (korelasi 0,97/0,96/0,74). | `src/eblup_twofold.cpp:19-24` |
| TFH-2 | MSE = $g_1+g_2+g_3$ (faktor 2 sudah ada di suku silang $2\mathrm{cov}_{vu}$), bukan $g_1+g_2+2g_3$ seperti penulisan Prasad–Rao untuk FH. | `src/eblup_twofold.cpp:407-408` |
| TFH-3 | $\mathrm{var}_v,\mathrm{var}_u,\mathrm{cov}_{vu}$ berasal dari **invers matriks informasi Fisher**, bukan bentuk tertutup $\mathrm{Var}(\hat\sigma^2)$ seperti FH (`VarA=2/\sum w^2`). | `src/eblup_twofold.cpp:235,324` vs `eblup_fh.cpp:186` |
| TFH-4 | Boundary: $\sigma_v^2,\sigma_u^2$ dipangkas ke 0 tanpa proyeksi konstrain; fallback skor diagonal bila informasi singular. | `src/eblup_twofold.cpp:201-211` |
| TFH-5 | Bootstrap menggambar $u^*$ untuk **semua** $N$ subarea padahal tak-sampel tidak pernah masuk $y^*$; $e^*$ hanya untuk baris sampel. | `src/eblup_twofold.cpp:476-477` |
| TFH-6 | MSE bootstrap memakai prediktor yang sama dengan estimasi titik, termasuk $\hat v_i$ untuk subarea tak-sampel (komentar kode menyebut ini perbaikan bug sebelumnya). | `src/eblup_twofold.cpp:34-35,509` |
| TFH-7 | `set.seed(seed)` dieksekusi untuk **semua** jalur, bukan hanya `mse="bootstrap"`. | `R/eblup_twofold.R:108` |
| TFH-8 | Negatif kecil hanya dipangkas bila $>-10^{-8}$; lebih negatif dibiarkan → RSE `NaN`. | `src/eblup_twofold.cpp:409,465` |
| TFH-9 | Basis AIC/BIC = jumlah area **dengan subarea tersampel**, berbeda dengan FH ($m$ = baris sampel). | `src/eblup_twofold.cpp:306-310` vs `eblup_fh.cpp:174` |

## D. Fay–Herriot spasial (`eblup_sfh`)

| ID | Butir | Lokasi |
|---|---|---|
| SFH-1 | Skor **ML** memakai kuadratik $y'PDPy$ dengan $P$ = proyektor REML, sedangkan informasi ML memakai $V^{-1}$ — ML tidak konsisten (kuadratik REML, informasi ML). | `src/eblup_sfh_core.cpp:330-348` |
| SFH-2 | `std.error` $\sigma_u^2,\rho$ dihitung dari $I_{dev}$ berbentuk **REML** meskipun `method="ML"`. | `src/eblup_sfh_core.cpp:502-510,647-648` |
| SFH-3 | `loglikelihood/AIC/BIC` selalu berbentuk ML (`log|V|`, $r'V^{-1}r$) juga untuk fit REML; tanpa koreksi $\log|X'V^{-1}X|$. | `src/eblup_sfh_core.cpp:465-470` |
| SFH-4 | Boundary $\rho$: selama iterasi di-clamp $\pm0.999$, tetapi pada hasil **akhir** dipetakan ke $\pm1.0$ persis, tempat $A$ singular → `inv_sympd` gagal → `pinv`. Jalur bootstrap mempertahankan $\pm0.999$. Dua jalur dengan perilaku boundary berbeda. | `src/eblup_sfh_core.cpp:373-374,389-399` |
| SFH-5 | $\sigma_u^2$ tidak pernah di-clamp nonnegatif selama iterasi; hanya $\sigma^2=\max(0,\cdot)$ di akhir. | `src/eblup_sfh_core.cpp:376,394` |
| SFH-6 | Dok menyebut kolom `tvalue`; kode mengembalikan `zvalue`. | `R/eblup_sfh.R:39` vs `eblup_sfh_core.cpp:629-635` |
| SFH-7 | Dok menyebut "analytical MSE approximation" untuk area tak-tersampel; aktual **hanya $g_1+g_2$**, tanpa $g_3/g_4$ dan tanpa koreksi ML. | `R/eblup_sfh.R:24-26,59-60` vs `eblup_sfh_core.cpp:596-599` |
| SFH-8 | Inkonsistensi sub-blok $W$: fit memakai $W_{ss}$, kriging memakai sub-blok dari $A_{full}^{-1}$; $[A_{full}^{-1}]_{ss}\neq[(I-\rho W_{ss})'(I-\rho W_{ss})]^{-1}$. | `src/eblup_sfh_core.cpp:257` vs `:580-587` |
| SFH-9 | `mse_method="npbmse"` + `method="ML"` baru dicek di C++ (bukan di R). | `R/eblup_sfh.R:167` vs `eblup_sfh_npbmse.cpp:32-34` |
| SFH-10 | Rumus bias-corrected $2(g_1+g_2)-\hat g_1-\hat g_2+\hat g_3$ memakai $\hat g_3$ dari bobot fit awal, bukan dari replikat; kesesuaian dengan rumus publik **TIDAK PASTI** (tidak ada rujukan yang mengikat rumus ini di repo). | `src/eblup_sfh_pbmse.cpp:127-128,164`; `eblup_sfh_npbmse.cpp:184-185,222` |

## E. Fay–Herriot spasio-temporal (`eblup_stfh`)

| ID | Butir | Lokasi |
|---|---|---|
| STFH-1 | Struktur spasial = efek spasial **konstan-waktu** $kron(\Omega_1,1_T)$ — bukan salah satu interaksi Knorr–Held I–IV; $u_2$ = AR(1) per-domain. Konsisten dengan `sae::eblupSTFH`, tetapi harus ditulis eksplisit di laporan. | `src/eblup_stfh_core.cpp:281,436-440` |
| STFH-2 | **Tidak ada MSE analitik**: `mse`/`rse` = `NA` kecuali `compute_mse=TRUE` (bootstrap) — bertolak belakang dengan SFH. | `R/eblup_stfh.R:41-42,57,269-272` |
| STFH-3 | Bootstrap STFH **tidak mengestimasi ulang komponen varian tiap replikat** (fixed-$\theta$), sedangkan `sae::pbmseSTFH` memanggil `eblupSTFH()` ulang tiap replikat; juga tanpa versi bias-corrected. | `src/eblup_stfh_pbmse.cpp:91,368-374,386` |
| STFH-4 | `method` di-hardcode `"REML"` — tidak ada pilihan ML (berbeda dengan `eblup_sfh`). | `R/eblup_stfh.R:225` |
| STFH-5 | `mse_pb` dibagi `n_valid` (bukan `B`) dan `res$B` diset ke `n_valid` → bisa `< B` yang diminta; SFH memaksa tepat $B$ replikat valid. | `src/eblup_stfh_pbmse.cpp:410-420`, `R/eblup_stfh.R:266` |
| STFH-6 | Kolom bernama `tvalue` tetapi p-value dihitung dari distribusi normal. | `src/eblup_stfh_core.cpp:483-486` |
| STFH-7 | $F^{-1}$ untuk `std.error` berasal dari iterasi terakhir $\theta_k$, sedangkan estimasi yang dilaporkan $\theta_{k+1}$. | `src/eblup_stfh_core.cpp:319-323` vs `:496-512` |
| STFH-8 | Cabang `param_invalid` mengembalikan `convergence = TRUE` dengan `estcoef = NULL` (praktis tak tercapai, tetapi struktur return inkonsisten); nama `goodness` berbeda dari SFH (`loglike` vs `loglikelihood`). | `src/eblup_stfh_core.cpp:346-374,478` |

## F. Lapisan Bayes/INLA (`hb_area`, `hb_unit`, `hb_twofold`)

Semua butir ini **terverifikasi lewat eksperimen kecil `INLA::inla()`** oleh catatan
`notes_bayes.md` (ditandai [E] di sana).

| ID | Butir | Lokasi |
|---|---|---|
| HB-1 | `st_interaction="type2"` **selalu gagal**: dua `f()` memakai kovariat `..time_id..` yang sama → INLA menolak ("Only one is allowed"). | `R/hb_area.R:539-545` |
| HB-2 | `type3`/`type4` dengan `spatial != "none"` **selalu gagal**: term spasial ditambah dua kali (blok spasial + blok interaksi) → INLA menolak kovariat ganda `..domain_id..`. | `R/hb_area.R:416-449` vs `:551-564` |
| HB-3 | `type3`/`type4` dengan `spatial="none"` → interaksi tidak dibuat **dan** tidak ada term `..domain_id..` ⇒ model tinggal efek waktu, `st_interaction` diabaikan diam-diam. | `R/hb_area.R:411,550,560` |
| HB-4 | `hyper` salah untuk `type3/type4`: selalu `list(prec = prior_prec)` → `bym` error "Unknown keyword 'prec'", `generic1` dipanggil dengan `graph=` (bukan `Cmatrix`), `slm` tanpa `args.slm`, `bym2` → `prior_phi` diam-diam diabaikan (jatuh ke default INLA `pc(0.5,0.5)`). | `R/hb_area.R:552-553,562-563` |
| HB-5 | `spatial="slm"` tidak punya guard `!= "separable"` → slm + separable menambah dua term spasial (dan menurunkan spasial ke `iid`). | `R/hb_area.R:487` vs `:415-448` |
| HB-6 | Argumen `link` **mati**: diteruskan tetapi tidak pernah dipakai; `control.predictor$link = 1` selalu. | `R/hb_area.R:51-52,295,333` vs `:587` |
| HB-7 | `prior_phi` dipakai ganda (juga untuk `beta` generic1 dan `rho` slm), hanya bila `prior_phi$prior != "pc"`; default paket (`"pc"`) → `beta` generic1 memakai default INLA `gaussian(0,0.1)` pada skala logit($\beta$). | `R/hb_area.R:443-447,481-484` |
| HB-8 | Prior default tidak konsisten antar fungsi: `hb_area` `pc.prec(1,0.01)`, `hb_unit` **`loggamma(0.01,0.01)`**, `hb_twofold` `pc.prec(1,0.01)`. | `R/hb_area.R:165`, `R/hb_unit.R:136`, `R/hb_twofold.R:152-153` |
| HB-9 | `random_effect_var` untuk `spatial="bym"` mengambil komponen **iid**, bukan komponen spasial (urutan baris `summary.hyperpar` INLA). | `R/inla_utils.R:158` |
| HB-10 | Fallback regex `grep("Precision for")` mengembalikan baris pertama; untuk Gaussian tanpa `vardir` baris pertama = "Precision for the Gaussian observations" → akan melaporkan $\sigma^2_{error}$ sebagai `random_effect_var` (saat ini tidak tercapai, tetapi rapuh). | `R/inla_utils.R:161` |
| HB-11 | `scale.model` diset untuk bym2/bym/besag dan rw1/rw2, **tidak** untuk generic1, slm, iid, ar1 — sementara dok menyebut "recommended for BYM2/Besag". | `R/hb_area.R:417,425,433,449,487,492` |
| HB-12 | `constr` tidak pernah dioper → mengikuti default INLA (bym2/bym/besag/rw1/rw2 `TRUE`); tidak didokumentasikan di roxygen. | `R/hb_area.R` (tidak ada `constr`) |
| HB-13 | `hb_area` menyusun `list(..., ...)` tanpa menimpa → `...` berisi `control.compute` dsb. memicu error argumen ganda dari `inla()`; `hb_unit`/`hb_twofold` menimpa. | `R/hb_area.R:593-601` vs `hb_unit.R:284`, `hb_twofold.R:289` |
| HB-14 | `method="laplace"` sering diam-diam jatuh ke INLA; gerbang `glmm` hanya untuk binomial/poisson tanpa struktur spasial/temporal, dan `strategy` dipaksa `"laplace"` pada jalur INLA-nya. | `R/hb_area.R:258-276` |
| HB-15 | Jalur `lme4` (`glmm`) tidak menghasilkan ketidakpastian (`sd/mse/rse/ci_* = NA`) dan memakai format `hyperpar`/`goodness` yang berbeda dari jalur INLA. | `R/hb_area.R:777-797` |
| HB-16 | Kelas tidak konsisten: `hb_area` → `c("fastsae_hb_area","fastsae")` (tanpa `fastsae_hb`), `hb_unit`/`hb_twofold` → `c(..., "fastsae_hb", "fastsae")`. | `R/hb_area.R:394,811`; `hb_unit.R:471`; `hb_twofold.R:503` |
| HB-17 | Domain unsampled: `hb_unit` menambah baris `y=NA` per domain; `hb_area`/`hb_twofold` **tidak** membuat baris → area yang tidak ada di `data` tidak pernah diprediksi (area dengan `y=NA` tetap dievaluasi). | `R/hb_unit.R:244-249`; `hb_area.R:275-276`; `hb_twofold.R:303` |
| HB-18 | Tipe interval berubah-ubah: jalur sampling memakai kuantil posterior, fallback memakai normal $\pm1.96$. | `R/hb_twofold.R:389-406`; `hb_unit.R:330-331` |
| HB-19 | `hb_unit` + `popsize_var`: kolom `linear_pred` tetap η model tanpa FPC sementara `hb` sudah dicampur $f_d\bar y_s+(1-f_d)\hat\theta$ ⇒ dua kolom tidak lagi saling invers; `sd_adj` diskalakan $(1-f_d)$ (aproksimasi). | `R/hb_unit.R:323-331,419` |
| HB-20 | `.convert_spatial_weights` membakar bobot ke graf biner untuk `inla.graph`, tetapi `W_mat` mentah dipakai apa adanya untuk `slm`; jalur `inla.graph` mengembalikan `W_mat=NULL` → `spatial="slm"` dengan `W` berupa graf akan gagal di `rowSums(W_mat)`. | `R/inla_utils.R:47-48,72-73`; `hb_area.R:452-455` |
| HB-21 | `control.fixed` tidak pernah di-set paket → default INLA berlaku: `mean=0`, `mean.intercept=0`, `prec=0.001`, `prec.intercept=0` (flat), `expand.factor.strategy="model.matrix"`. Ini adalah prior tetap yang dipakai laporan. | tidak ada di `R/*.R`; diverifikasi `INLA::inla.set.control.fixed.default()` |

## G. Utilitas (diagnostik, benchmarking, ekspor, peta, simulasi)

| ID | Butir | Lokasi |
|---|---|---|
| UT-1 | Tes Brown: regresi **OLS tanpa bobot**, sementara formulasi yang lazim dikutip memakai WLS dengan bobot $\propto 1/(v_d+\mathrm{mse}_d)$. **TIDAK PASTI dari repo** (tidak ada teks rujukan dalam repo). | `R/diagnose.R:144` |
| UT-2 | Statistik Wald dibagi 2 lalu dibandingkan dengan $F(2,n-2)$, bukan $\chi^2(2)$ asimptotik $c'V^{-1}c$. Kedua versi beredar; angka p berbeda. | `R/diagnose.R:152-155` |
| UT-3 | Keputusan GOF memakai penerimaan interval dua sisi `[qchisq(α/2), qchisq(1−α/2)]`, sedangkan `p_value` yang dicetak adalah uji ekor-atas tunggal → `is_good_fit` bisa FALSE walau `p_value ≥ α`. | `R/diagnose.R:176-180` vs `:338,447-449` |
| UT-4 | Metrik "Mean Relative RMSE" sebenarnya rerata galat relatif absolut (MAPE-like), bukan akar kuadrat rerata kuadrat. | `R/diagnose.R:256,267,475` |
| UT-5 | "95% CI Coverage" dicetak tanpa memverifikasi level CI; lolos bila ≥90%. | `R/diagnose.R:259-262,477-481` |
| UT-6 | Metode `optimal`: penyesuaian $\delta_d=(\mathrm{mse}_d/w_d)\lambda$ dengan $\lambda=\Delta/\sum\mathrm{mse}_d$. Solusi KKT standar dari $\min\sum\delta_d^2/\mathrm{mse}_d$ s.t. $\sum w_d\delta_d=\Delta$ adalah $\delta_d\propto w_d\,\mathrm{mse}_d$. Constraint tetap terpenuhi (kode diverifikasi benar) tetapi **alokasi berbeda** bila bobot tidak seragam; dengan bobot seragam keduanya $\propto$ mse. Variabel juga dinamai `inv_prec` padahal isinya mse. | `R/benchmark.R:476-479,413-416` |
| UT-7 | Tahap-1 `logit` selalu memakai bentuk $\sum W_{share}\,\mathrm{logit}^{-1}(\cdot)$ (bertipe mean) meski `type="total"` memakai $\sum$ target → ketidakcocokan skala. | `R/benchmark.R:402` vs `:420` |
| UT-8 | `group_orig` dihitung dengan bobot ternormalisasi (skala mean) meski `type="total"`; pada hierarki berskala total target bisa tidak konsisten. | `R/benchmark.R:333-334,383-386` vs `:451-452` |
| UT-9 | `add_reliability_flags` didokumentasikan menerima "objek atau data frame", tetapi selalu mengembalikan `df` saja (penugasan ke `object$df_*` tidak efektif). | doc `R/export_sae.R:23` vs `:68-75` |
| UT-10 | `export_sae` didokumentasikan "hingga tiga sheet"; kode bisa menghasilkan **empat** (`Benchmarked_Groups`), dan mode CSV hanya menulis `Estimates` + `_summary.csv`. | doc `R/export_sae.R:89-99` vs `:153-156,173-181` |
| UT-11 | `sim_area_data()` menghasilkan kolom `domain` + `truth_gaussian` (15 kolom), sedangkan dataset terlampir `sim_area` memakai `area` dan 14 kolom tanpa truth. | `R/sim_area_data.R:169-185` vs `R/data.R:120-136` |
| UT-12 | `map_sae`: argumen `thresholds=c(20,30)` diteruskan tetapi **tidak pernah dipakai**; pemotongan reliability hard-coded 20/30. | `R/map_sae.R:107,169,602` vs `:482-487` |
| UT-13 | `diagnose()` bukan generic S3; tidak ada `summary`/`plot` untuk `fastsae_diagnose` (hanya `print` + `autoplot`); `@return` tidak menyebut `bayesian_metrics`. | `R/diagnose.R:53,16-25` vs `:380` |
| UT-14 | `compare_sae` hanya mendaftarkan method `.default`; `thresholds` disimpan tetapi tidak memengaruhi metrik. | `NAMESPACE:31`; `R/compare_sae.R:183,343-348` |
| UT-15 | `autoplot.list` memanggil `NextMethod()` bila elemen bukan fastsae (bisa gagal dispatch); `.autoplot_multi_comparison` bergantung pada `names(x)`. | `R/autoplot.R:68-71,257` |
| UT-16 | Pemetaan dua model menyusun `data.frame(matched_sf[, c("domain","geometry")], ...)` — akan gagal bila kolom geometri `sf` bukan bernama `"geometry"`. | `R/map_sae.R:217-220,432-435` |

## H. Dokumentasi vs data benchmark (README, vignette, artikel)

Angka dipercaya dari `inst/extdata/*.rds`; klaim dokumentasi yang bertolak belakang dicatat di sini.
Perhitungan: `paper/_notes/bench_methodology.md` §6, §8.

| ID | Klaim | Data | Status |
|---|---|---|---|
| DOC-1 | `README.Rmd:118` — tabel ditandai "Benchmark across $n=1{,}000$ areas" | angka tabel = rata-rata lintas $n=30\ldots1000$; nilai $n=1000$ sesungguhnya berbeda (mis. FH 0,00413357905927 / 1,50363106854 / 51,0223606815 s) | **BERLAWANAN** (label $n$ salah) |
| DOC-2 | `README.Rmd:122` — "~190×" dan "~6,400×" | rasio mean = **215,7768×** / **7150,544×**; rasio pada $n=1000$ = **363,7601×** / **12343,39×** | **BERLAWANAN** dengan keduanya |
| DOC-3 | `README.Rmd:124` — fastsae "< 10 MB" | fastsae >10 MB pada `seblup` $n=1000$ **23,154800415**, `beta` 133,58505249, `spatial_beta` 161,025146484, `spatio_temporal` 109,008224487, `hb` 57,7934341431 | **BERLAWANAN**, kecuali dibatasi pada FH (maks 0,185699462891) |
| DOC-4 | `README.Rmd:124` — pembanding "~400 MB" / "~800 MB" | mean `eblup` sae = 16,334924062093; mean `seblup` sae = 471,57688013713; emdi 835,88609568278 | campur dua file; untuk FH tidak cocok |
| DOC-5 | `README.Rmd:34-35` — "roughly 50 to 6,000 times faster" | rentang teramati 2,639× … 12343,39× | **BERLAWANAN** di kedua ujung |
| DOC-6 | `README.Rmd:87` — "10–100× speedups" | 1,887× … 544,5762×; pada Bayesian `hbsae` justru lebih cepat (rasio 0,0162–0,2127) | **BERLAWANAN** |
| DOC-7 | `vignettes/benchmarks.Rmd:39` — emdi spasial "8,69 s" | 9,186901320597 s | **BERLAWANAN** |
| DOC-8 | Label "**Peak** Memory" (`README.Rmd:124`, `vignettes/benchmarks.Rmd:41-43`, `paper/article.tex:284`) | nilai yang dilaporkan = **mean** `bench::mark()$mem_alloc`, bukan puncak per-n | **SALAH LABEL** |
| DOC-9 | `paper/article.tex:100,480` — "73× to 544×" | nilai serendah 1,887× | di luar rentang teramati |
| DOC-10 | `vignettes/benchmarks.Rmd:407-414,634-639` — kolom Spatial/ST untuk $n\ge100$ | tidak ada data pendukung di `inst/extdata` | **TIDAK DIDUKUNG DATA** |
| DOC-11 | Reproduksibilitas | 4 file (`eblup_benchmark.rds`, `seblup_benchmark.rds`, `steblup_benchmark.rds`, `beta_estimates_n1000.rds`) **tidak punya skrip pembuat** di repo; `main.R:158-160` membacanya dari `inst/exdata/` (salah ketik) yang tidak ada | hasil tidak dapat direproduksi dari repo |

**Konsekuensi untuk laporan:** semua tabel/figur dihasilkan langsung dari
`inst/extdata/*.rds` dengan skrip `scripts/*.R`; klaim README/vignette yang
bertolak belakang ditulis sebagai temuan di bagian 3 (*caveat*), bukan dipakai
sebagai angka laporan.
