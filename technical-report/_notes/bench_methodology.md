# catatan: metodologi benchmark (`bench_methodology.md`)

Dibuat oleh subagen (baca-saja) untuk laporan teknis. **Tidak ada file di luar
`paper/_notes/` yang diubah**; tidak ada benchmark dijalankan ulang. Semua angka
di bawah dihitung langsung dari file `.rds` (R, `digits` tinggi) atau dikutip
persis dari sumber dengan nomor baris. File pasangan: `bench_data_dump.md`
(dump mentah 9 file `inst/extdata/*.rds`).

Sumber yang dibaca:
`benchmarks/*.R` (11 skrip), `benchmarks/output/` (isi), `inst/extdata/` (9 file),
`vignettes/benchmarks.Rmd`, `README.Rmd`, `paper/article.tex`, `main.R`,
`paper/_notes/README_for_agents.md`.

---

## 1. Inventaris skrip di `benchmarks/`

| # | Skrip | Comparator (dari kode) | Apa yang diukur | n (default) | Iterasi | Seed | Paralel | Menulis |
|---|-------|------------------------|-----------------|-------------|---------|------|---------|---------|
| 1 | `bench_beta_tipsae.R` | `fastsae::hb_area` vs `tipsae::fit_sae` (`:9-10`) | `bench::mark(memory=TRUE, time_unit="s")` (`:60,84-85`) | `c(250,500)` (`:28`) | `iter_arg=20L` (`:22`) | data `42+n` (`:54`), tipsae `seed=42` (`:79`) | – | `inst/exdata` + `inst/extdata`/`beta_benchmark.rds` (`:119`; dir `:35`) |
| 2 | `bench_spatial_beta.R` | fastsae vs `tipsae` + `sf`,`spdep` (`:9-12`) | `bench::mark` (`:93,121-122`) | `c(30,50,100)` (`:30`) | `20L` (`:24`) | `42+n` (`:87`), tipsae `seed=42` (`:116`) | – | `spatial_beta_benchmark.rds` ke 2 lokasi (`:43-51`) |
| 3 | `bench_spatio_temporal_beta.R` | fastsae vs `tipsae` (`:9-10`) | `bench::mark` (`:108,140-141`) | `c(30,…,1000)`, `T_arg=5L` (`:30,36`) | `5L` (`:24`) | `42+n` (`:101`), tipsae `seed=42` (`:135`) | – | `spatio_temporal_beta_benchmark.rds` (`:50-56`) |
| 4 | `bench_hb_normal.R` | fastsae vs `hbsae::fSAE.Area` vs `fastsaehb` (`:10-12,34`) | `bench::mark` (`:57,81-82`) | `c(30,50,100,250)` (`:30`) | `10L` (`:24`) | `42+n` (`:51`) | – | `hb_normal_benchmark.rds` **dan** `hb_benchmark.rds` (`:107-108`) |
| 5 | `update_hb_benchmark.R` | fastsae vs `hbsae` (baris `fastsaehb` dihapus, `:23-25`) | `bench::mark` (`:43,60-61`) | `c(500,1000)` (`:30`) | `10` (`:58`) | `42+n` (`:37`) | – | kedua file `hb_*.rds` (`:94-95`) |
| 6 | `run_hb_area_full_benchmark.R` | fastsae vs `fastsaehb` vs `tipsae (Stan HMC)`; `library(rstan)` (`:2-10,97,181`) | `Sys.time()` diff (`:33,43,46,54,105,120,…`) | tidak ada loop n (skenario tunggal) | – | `set.seed(123/42)` (`:32,104,201`), data `seed=123/42` | – | **tidak menulis file apa pun** (hanya konsol) |
| 7 | `run_multipackage_hb_benchmark.R` | fastsae vs `hbsae`, `CARBayes`, `fastsaehb`, `tipsae` (`:2-15`) | `Sys.time()` diff + `bench::mark(iterations=5, memory=TRUE)` (`:258-265`) | – | 5 (`:263`) | `set.seed(123/42)` (`:104,248`) | – | **tidak menulis file apa pun** |
| 8 | `compare_ebp_tipsae_truth.R` | fastsae vs `tipsae` vs direct (akurasi vs `theta_true`) | `Sys.time()` diff (`:102-129,156-187`) | `D=60` (`:29`) | – | `seed_val=42` (`:30-31`) | – | `benchmarks/output/comparison_beta_truth_results.rds` (`:464`; `:25`) |
| 9 | `sim_finite_population_beta.R` | fastsae vs `tipsae` (finite population) | `Sys.time()` diff (`:171-247`) | – | – | `set.seed(2024L)` (`:27`), tipsae `seed=42` (`:199,245`) | – | `benchmarks/output/finite_population_simulation_results.rds` (`:430`) |
| 10 | `sim_mc_finite_population.R` | direct vs fastsae vs tipsae, MC | `Sys.time()` per rep + total (`:154-191,214-223`) | `D=60` (`:44`), `R=100` (`:20`) | R=100 | `set.seed(2024L)` (`:43`), `rep_seed=5000+r` (`:119-120`) | `mclapply` `cores=min(6, detectCores()-2)` (`:26,216`) | `output/mc_finite_population_results_R100.rds` (`:413`) |
| 11 | `sim_mc_finite_pop_d500.R` | idem, skala nasional | idem (`:154-191,214-223`) | `D=500` (`:44`), `R=25` (`:20`) | R=25 | `set.seed(999L)` (`:43`), `rep_seed=8000+r` (`:121-122`) | `mclapply` `cores=min(5, detectCores()-2)` (`:26,216`) | `output/mc_finite_population_results_D500_R25.rds` (`:400`) |

Catatan skrip:
- Skrip 1–3 menulis ke `inst/exdata` (typo, dibuat sendiri `:35/:38/:45`)
  **dan** `inst/extdata`; `inst/exdata` tidak ada di repo (hanya `inst/extdata`,
  `inst/spelling`) → bacaan `main.R:158-160` (`inst/exdata/…`) akan gagal.
- Skrip 4 `:38` memakai `out_dirs <- c("inst/extdata","inst/extdata")` (duplikat).
- Semua skrip 1–5 adalah *incremental*: baris lama per-`n` dihapus lalu diganti
  (mis. `bench_beta_tipsae.R:101-107`) → file bisa berisi campuran run.

## 2. Script → file `inst/extdata` (jejak pembuat)

| File `inst/extdata` | Pembuat | Bukti |
|---|---|---|
| `beta_benchmark.rds` | `bench_beta_tipsae.R` | `:119` |
| `spatial_beta_benchmark.rds` | `bench_spatial_beta.R` | `:49` |
| `spatio_temporal_beta_benchmark.rds` | `bench_spatio_temporal_beta.R` | `:56` |
| `hb_normal_benchmark.rds` | `bench_hb_normal.R` **+** `update_hb_benchmark.R` | `:107` / `:94` |
| `hb_benchmark.rds` | sama (tulis identik) | `:108` / `:95` |
| `eblup_benchmark.rds` | **tidak ada skrip** | hanya dibaca `main.R:158` |
| `seblup_benchmark.rds` | **tidak ada skrip** | hanya dibaca `main.R:159` |
| `steblup_benchmark.rds` | **tidak ada skrip** | hanya dibaca `main.R:160` |
| `beta_estimates_n1000.rds` | **tidak ada skrip** | tidak direferensikan di kode mana pun |

Verifikasi: `grep -rn "eblup_benchmark\|seblup_benchmark\|steblup_benchmark\|beta_estimates_n1000"`
seluruh repo hanya menemukan `main.R:158-160` dan catatan `_notes`;
`git log -S "eblup_benchmark" --all` hanya menunjuk pembaca (`README.Rmd`,
`main.R`, `docs/steblup.html`) → **skrip pembuat 4 file itu tidak pernah ada di
riwayat git yang tersedia** (`TIDAK PASTI` apakah dijalankan di luar repo).

## 3. `benchmarks/output/` (isi)

4 file RDS + 13 PNG (tanggal file 27 Sep):

- `comparison_beta_truth_results.rds` (12 226 B): `param_recovery`, `metrics_cs`,
  `metrics_sp`, `data_cs`/`data_sp` (60×14), `runtimes$cs = 1.997771025,
  2.080733061, 1.041527300`; `runtimes$sp = 1.13122797, 20.56911492, 18.18299712`.
- `finite_population_simulation_results.rds` (12 031 B): `eval_cs`, `eval_sp`
  (RMSE/MAE/coverage/efficiency/runtime, 60 domain), `data_cs`/`data_sp` (60×18).
- `mc_finite_population_results_R100.rds`: `R = 100`, `D = 60`,
  `mc_elapsed_sec = 269.121217`.
- `mc_finite_population_results_D500_R25.rds`: `R = 25`, `D = 500`,
  `mc_elapsed_sec = 225.4484899`.
- PNG: `plot1_estimates_vs_truth.png`, `plot2_method_agreement.png`,
  `plot3_credible_interval_coverage.png`, `plot4_parameter_recovery.png`,
  `plot_finite_pop_{coverage,errors,estimates_vs_true}.png`,
  `plot_mc_empirical_{coverage,rmse,rrmse_box}.png`,
  `plot_d500_mc_{agreement,empirical_rmse,rrmse_box}.png`.

Tidak ada hasil waktu/memori FH di `output/` — yang itu di `inst/extdata/`.

## 4. Metadata: tercatat vs tidak tercatat

| Item | Status |
|---|---|
| Hardware | Hanya `paper/article.tex:277`: "on an Apple Silicon processor under R 4.4.2". Model CPU/jumlah inti/RAM/OS: **tidak tercatat** |
| R version | `R 4.4.2` hanya di `article.tex:277`; tidak ada `sessionInfo()` di skrip/vignette |
| Threads | **Tidak ada** pengaturan `OMP_NUM_THREADS`/`omp_set_num_threads` di skrip mana pun. `README.Rmd:126` hanya: "Timings depend on hardware and the number of threads" |
| Tanggal run, seed, versi paket pembanding | **tidak tercatat** di dalam file `.rds` (tidak ada atribut/seed/kolom tanggal) |
| Fungsi ukur waktu | `bench::mark()` → kolom `Min (s)` = `res$min`, `Median (s)` = `res$median`, `Iterations/sec` = `res$itr/sec`, `n_iter` = `res$n_itr` (`bench_beta_tipsae.R:92-98`); skrip 6,7,8,9,10,11 memakai `difftime(Sys.time(), …)` |
| Fungsi ukur memori | `Memory (MB)` = `as.numeric(res$mem_alloc)/1024^2` (`bench_beta_tipsae.R:96`) = **total alokasi selama run**, bukan peak RSS. Label "Peak" di README/vignette/article adalah istilah yang tidak persis |
| n_iter aktual vs permintaan | `beta_benchmark.rds`: 3/16/17/20; `eblup`: 16/19/20; `seblup`: 17/19/20; `steblup`: 20 (n=1000 → 3); `spatial_beta`: 5/20; `spatio_temporal`: 5; `hb*`: 10. Beberapa run besar tampaknya memakai `--iter=3/5` (**tidak tercatat di file mana argumen itu**) |

## 5. Angka dasar hasil (dihitung dari file)

Rasio = `Median (s)` pembanding / `Median (s)` fastsae.

**Rata-rata lintas n (n = 30, 50, 100, 250, 500, 1000):**

| File | fastsae (s) | pembanding (s) | rasio | `Memory (MB)` mean fastsae / pembanding |
|---|---|---|---|---|
| `eblup_benchmark.rds` | 0.00134774859665 | sae 0.29081284627803 / emdi 9.63713592576581 | 215.7768× / 7150.544× | 0.055046081543 / 16.334924062093 / 824.466003417969 |
| `seblup_benchmark.rds` | 0.167441420422 | sae 12.571459845169 / emdi 9.186901320597 | 75.07975× / 54.86636× | 5.14562861125 / 471.57688013713 / 835.88609568278 |
| `steblup_benchmark.rds` | 12.450401489 | sae 352.791667070 | 28.33577× | 0.288579305013 / 7821.687708536783 |
| `beta_benchmark.rds` | 1.36922060452 | tipsae 48.66849496249 | 35.54467× | 81.0358517965 / 5860.8719825745 |
| `spatial_beta_benchmark.rds` | 1.18091287303 | tipsae 166.75466978117 | 141.2083× | 70.6704216003 / 7366.1817372640 |
| `spatio_temporal_beta_benchmark.rds` (hanya n=30,50) | 1.55486241356 | tipsae 724.37321689632 | 465.8761× | 105.869155884 / 3913.604370117 |
| `hb_benchmark.rds` = `hb_normal_benchmark.rds` | 1.2636070713440 | **hbsae 0.0956100388624** | **0.07566437× (hbsae lebih cepat)** | 23.7776565552 / 67.0100059509 |

**Baris n = 1000:**

| File | fastsae | pembanding | rasio |
|---|---|---|---|
| `eblup` | 0.00413357905927 s / 0.185699462891 MB | sae 1.50363106854 s / 71.976524353 MB; emdi 51.0223606815 s / 3773.59893799 MB | 363.7601× / 12343.39× |
| `seblup` | 0.832157648023 s / 23.154800415 MB | sae 67.151753158 s / 1920.00576019 MB; emdi 51.0593093895 s / 3772.56869507 MB | 80.69595× / 61.35774× |
| `steblup` | 58.100478563 s / 0.893669128418 MB | sae 1786.49415005 s / 32478.3168335 MB | 30.74836× |
| `beta` | 1.63966884301 s / 133.58505249 MB | tipsae 187.161585684 s / 18335.2039948 MB | 114.146× |
| `spatial_beta` | 1.72812134121 s / 161.025146484 MB | tipsae 404.341828743 s / 23038.6522675 MB | 233.9777× |
| `hb` | 1.27669157891 s / 25.0034255981 MB | hbsae 0.271584635484 s / 220.352233887 MB | 0.2127253× |
| `spatio_temporal_beta` | **tidak ada baris n=1000** (hanya n=30, 50) | – | – |

**Rentang rasio per-n (min → max):** eblup-sae 2.639× (n=30) → 363.7601×;
eblup-emdi 18.5× → 12343.39×; seblup-sae 5.583× → 80.69595×;
seblup-emdi 11.21× (n=100) → 61.35774×; steblup-sae 3.227× → 30.74836×;
beta-tipsae **1.887× (n=30)** → 114.146×; spatial_beta-tipsae 41.91× (n=30),
72.84× (n=100) → 233.9777×; spatio_temporal-tipsae 374.9× (n=30), **544.5762× (n=50)**;
hb **hbsae selalu lebih cepat** (0.0162×–0.2127×).

`beta_estimates_n1000.rds` (1000 baris, 5 kolom) → Pearson 0.9993429,
Spearman 0.9992081, MAE 0.0196916, RMSD 0.0236901, r(MSE) 0.9631480,
MAD(MSE) 0.0001801.

---

## 6. Klaim vs data (kontradiksi eksplisit)

### 6.1 `README.Rmd`

| Klaim (baris) | Data file | Status |
|---|---|---|
| `:34-35` "roughly **50 to 6,000** times faster than `sae` and `emdi`" | rentang teramati 2.639× … 12343.39× (sae 2.639–363.76; emdi 18.5–12343.39) | **BERLAWANAN** di kedua ujung |
| `:87` "**10–100×** speedups for large problems" | Bayesian: 1.887× (beta n=30) … 544.5762×; hbsae 0.0162–0.2127 (lebih cepat dari fastsae) | **BERLAWANAN** |
| `:118` header "Benchmark across **n = 1,000** areas" | angka di tabel = rata-rata lintas n=30…1000 (0.00134774859665; 0.29081284627803; 9.63713592576581; 0.167441420422; 12.571459845169; 9.186901320597). Nilai n=1000 sesungguhnya: 0.00413357905927 / 1.50363106854 / 51.0223606815 dan 0.832157648023 / 67.151753158 / 51.0593093895 | **BERLAWANAN** (label n salah) |
| `:122` FH "0.0015 s / 0.291 s (~190x) / 9.64 s (~6,400x)" | mean 0.00134774859665 / 0.29081284627803 / 9.63713592576581 → rasio **215.7768× / 7150.544×**; n=1000 → **363.7601× / 12343.39×** | **BERLAWANAN** (190x/6400x tidak cocok dengan salah satu pun) |
| `:123` Spatial "0.165 s / 12.60 s (~76x) / 8.69 s (~53x)" | mean 0.167441420422 / 12.571459845169 / **9.186901320597** → 75.07975× / 54.86636×; n=1000 → 80.69595× / 61.35774× | 76x/53x ≈ mean (dalam pembulatan); **emdi 8.69 s BERLAWANAN** (9.186901320597); jika `:118` dibaca literal (n=1000) → 80.70×/61.36× ≠ 76x/53x |
| `:124` "**< 10 MB**" (fastsae) | fastsae >10 MB pada: `seblup` n=1000 **23.154800415**, `beta` 133.58505249, `spatial_beta` 161.025146484, `spatio_temporal` 109.008224487, `hb` 57.7934341431 | **BERLAWANAN** kecuali dibatasi EBLUP saja (maks 0.185699462891) |
| `:124` sae "**~400 MB**" | mean `eblup` 16.334924062093; mean `seblup` 471.57688013713 (maks 1920.00576019) | **BERLAWANAN** untuk FH (16.3), hanya kasar untuk SFH; tampak campur dua file |
| `:124` emdi "**~800 MB**" | mean 824.466003417969 (eblup) / 835.88609568278 (seblup); puncak 3773.59893799 / 3772.56869507 | mean ≈ didukung; sebagai "peak" **BERLAWANAN** |
| `:126` "Timings depend on hardware and the number of threads" | satu-satunya keterangan hardware; jumlah thread tidak pernah diatur/dicatat | didukung sebagai disclaimer, tapi **spek hardware tidak tercatat** |

### 6.2 `vignettes/benchmarks.Rmd`

| Klaim (baris) | Data file | Status |
|---|---|---|
| `:25` "Up to **364x** … **12,300x** … at n=1000" | 363.7601× / 12343.39× | **DIDUKUNG** |
| `:26` "Up to **80x** (Spatial, n=1000)" | 80.69595× | **DIDUKUNG** |
| `:27` "Up to **31x** (Spatio-Temporal, n=1000)" | 30.74836× | **DIDUKUNG** |
| `:28` "Peak memory footprint stays under **25 MB**" | keluarga FH/SFH/STFH ✓ (maks 0.185699462891 / 23.154800415 / 0.893669128418); keluarga Bayesian ✗ (beta 133.58505249, spatial_beta 161.025146484, spatio_temporal 109.008224487, hb 57.7934341431) | **sebagian BERLAWANAN**; juga nilai di tabel adalah mean, bukan "peak" |
| `:38` EBLUP "0.0015 / 0.291 / 9.64" | 0.00134774859665 / 0.29081284627803 / 9.63713592576581 | 0.291 & 9.64 ✓; **0.0015 sedikit meleset** dari mean |
| `:39` Spatial "0.165 / 12.60 / **8.69**" | 0.167441420422 / 12.571459845169 / **9.186901320597** | **emdi 8.69 BERLAWANAN**; dua lainnya ≈ |
| `:40` ST "12.5 / 353.0" | 12.450401489 / 352.791667070 (`steblup`) | **DIDUKUNG** |
| `:41` "Peak Memory (EBLUP) 0.055 / 16.3 / 824" | mean 0.055046081543 / 16.334924062093 / 824.466003417969; puncak 0.185699462891 / 71.976524353 / 3773.59893799 | angka ✓ sebagai **mean**, label "Peak" **BERLAWANAN** |
| `:42` "Peak Memory (Spatial) 5.15 / **408** / **824**" | mean `seblup` 5.14562861125 / **471.57688013713** / **835.88609568278** | 5.15 ✓; **408 tidak ada di file mana pun**; **824 = mean emdi `eblup` (824.466003417969)** → kemungkinan salah salin |
| `:43` "Peak Memory (ST) 0.289 / 7,822" | mean `steblup` 0.288579305013 / 7821.687708536783 | ✓ sebagai mean; label "Peak" ✗ (puncak 0.893669128418 / 32478.3168335) |
| `:44` "~364x / ~12,300x slower" | 363.7601× / 12343.39× | **DIDUKUNG** |
| `:89-143` `RAW_DATA` grafik (median/IPS/mem per n) | dicocokkan per baris, mis. n=1000: FH 0.00413/0.186/176; SFH 0.832/23.15/1.19; STFH 58.10/0.894 & 1786.49/32478.3 = nilai file | **DIDUKUNG** |
| `:378-414` `BETA_DATA.ST`, baris **n = 100** (`:412-413`: 2.626624 / 1840.6304 s) | `spatio_temporal_beta_benchmark.rds` hanya n=30,50 | **TIDAK DIDUKUNG** (tidak ada data mentahnya) |
| `:628-640` metrik **kolom Beta** (0.99934 / 0.99921 / 0.01969 / 0.02369 / 0.96315 / 0.00018) | persis dihitung dari `beta_estimates_n1000.rds` (0.9993429 / 0.9992081 / 0.0196916 / 0.0236901 / 0.9631480 / 0.0001801) | **DIDUKUNG** |
| `:634-639` kolom **Spatial Beta** (0.99891, 0.01456, 0.01713, 0.89254, 0.00012) dan **Spatio-Temporal** (0.99614, 0.00120, 0.01458, 0.67957, 0.00025) | tidak ada file estimasi spatial/ST di repo | **TIDAK DAPAT DIVERIFIKASI** (`TIDAK PASTI`) |
| `:37` vs `README.Rmd:120` label paket: "sae (Molina & Rao)" vs "sae (Molina & Marhuenda)" | – | inkonsistensi sitasi (bukan angka) |

### 6.3 `paper/article.tex`

| Klaim (baris) | Data file | Status |
|---|---|---|
| `:27` "up to 364× FH, 80× SFH, 544× spatio-temporal Beta; peak RAM reductions up to 99.9%" | 363.7601 / 80.69595 / 544.5762; pengurangan `eblup` n=1000 vs emdi = 1 − 0.185699462891/3773.59893799 = 99.9951% | **DIDUKUNG** |
| `:99` "up to 364× sae / 12,300× emdi at D=1,000; … gigabytes → under 25 megabytes" | 363.7601 / 12343.39; fastsae D=1000: FH 0.185699462891, SFH 23.154800415, STEBLUP 0.893669128418 MB; emdi 3773.59893799 MB, sae STEBLUP 32478.3168335 MB | **DIDUKUNG** (untuk keluarga FH) |
| `:100` dan `:480` "**73× to 544×** speedups over Stan-based MCMC" | 544.5762 ✓ (spatio n=50); 72.84 (spatial n=100) ≈73 ✓; tetapi nilai teramati juga 1.887×, 2.835×, 6.473× (beta n=30/50/100) dan 41.91× (spatial n=30) yang **di bawah 73** | **SEBAGIAN BERLAWANAN**: "73–544" bukan rentang hasil, hanya seleksi; `tipsae` disebut Stan HMC di `run_hb_area_full_benchmark.R:97,181` tetapi skrip itu tidak menulis file |
| `:277` "bench package on an Apple Silicon processor under R 4.4.2" | satu-satunya pernyataan hardware; tidak ada jumlah inti/OS/RAM/tanggal | didukung sebagai klaim, **detail tidak tercatat** |
| `:288-305` tabel `tab:bench_freq` (median & "Peak Memory") | seluruh baris yang dicocokkan identik dengan `eblup/seblup/steblup_benchmark.rds` (mis. D=1000: 0.00413/1.50363/51.0224; 0.83216/67.1518/51.0593; 58.1005/1786.49; mem 0.186/71.977/3773.60; 23.155/1920.01/3772.57; 0.894/32478.3) | **DIDUKUNG**, tetapi 4 file itu **tidak punya skrip pembuat** (bagian 2) dan judul kolom "Peak Memory" = `mem_alloc` |
| `:307,310` "364× FH, 80× SFH, 31× ST-FH at D=1,000"; "4 ms vs 1.5 s (364×) … 51 s (12,343×)"; "58 s, <1 MB" | 363.7601 / 80.69595 / 30.74836; 0.00413357905927 vs 1.50363106854 vs 51.0223606815; 58.100478563 s / 0.893669128418 MB | **DIDUKUNG** |
| `:337` "1.73 s vs 404 s (234×)"; "908 s vs 1.67 s (544×)"; "r=0.9993, ρ=0.9992, MAD=0.0197" | 1.72812134121 / 404.341828743 → 233.9777×; 908.173061557 / 1.66766950605 → 544.5762×; 0.9993429 / 0.9992081 / 0.0196916 | **DIDUKUNG** |
| `:478` "RAM under 25 MB at D=1,000 (vs 32 GB in sae)" | 23.154800415 MB vs 32478.3168335 MB (≈31.7 GiB / 32.5 GB desimal) | **DIDUKUNG** |
| implikasi "fastsae tercepat" | `hb_benchmark.rds`: hbsae 4.7×–62× lebih cepat dari fastsae (rasio 0.0162–0.2127) | **BERLAWANAN** dengan bacaan "selalu tercepat" (README/vignette tidak menyebut hbsae) |

## 7. Apakah hasil mentah n = 1000 ada?

- **Ya** (baris n=1000 ada di file): `eblup_benchmark.rds`,
  `seblup_benchmark.rds`, `steblup_benchmark.rds`, `beta_benchmark.rds`,
  `spatial_beta_benchmark.rds`, `hb_benchmark.rds`, `hb_normal_benchmark.rds`,
  dan `beta_estimates_n1000.rds` (estimasi, bukan waktu).
- **Tidak ada**: `spatio_temporal_beta_benchmark.rds` (hanya n=30, 50) →
  grafik `BETA_DATA.ST` di `vignettes/benchmarks.Rmd:407-414` (baris n=100: `:412-413`) untuk n=100/250/500/1000
  **tidak punya data mentah di repo**.
- Tidak ada hasil mentah n=1000 untuk klaim "73×" (angka itu dari spatial n=100)
  dan untuk metrik korelasi Spatial/ST (`vignettes/benchmarks.Rmd:634-639`).

## 8. Ringkasan kontradiksi (prioritas perbaikan)

1. `README.Rmd:118` menandai tabel "n = 1,000" tetapi angkanya rata-rata n=30…1000.
2. `README.Rmd:122` "~190x / ~6,400x" vs data 215.7768× / 7150.544× (mean) atau 363.7601× / 12343.39× (n=1000).
3. `README.Rmd:124` "< 10 MB" dibantah `seblup` n=1000 = 23.154800415 MB (dan beta 133.58505249, spatial 161.025146484, hb 57.7934341431); "~400 MB" tak cocok untuk FH (16.334924062093).
4. `README.Rmd:34-35` "50 to 6,000×" vs rentang 2.639×–12343.39×; `README.Rmd:87` "10–100×" vs 1.887×–544.5762× (dan hbsae lebih cepat).
5. `vignettes/benchmarks.Rmd:39` emdi spatial "8.69 s" vs 9.186901320597 s; `:42` "408 MB" & "824 MB" vs 471.57688013713 & 835.88609568278 (824 = mean emdi file EBLUP → salah salin).
6. Label "**Peak** Memory" di `README.Rmd:124`, `vignettes/benchmarks.Rmd:41-43`, `paper/article.tex:284` untuk nilai yang sebenarnya **mean** dari `bench::mark()$mem_alloc` (bukan puncak; puncak per-n jauh lebih besar).
7. `paper/article.tex:100,480` "73× to 544×" bukan rentang teramati (nilai serendah 1.887×).
8. 4 file (`eblup/seblup/steblup/beta_estimates_n1000`) tidak punya skrip pembuat → hasil tidak dapat direproduksi dari repo; `main.R:158-160` membacanya dari `inst/exdata/` yang tidak ada.
9. `vignettes/benchmarks.Rmd:407-414` (ST n≥100) dan `:634-639` kolom Spatial/ST tidak punya data pendukung.

## 9. Item `TIDAK PASTI`

- Asal-usul "8.69 s" dan "408 MB" (`vignettes/benchmarks.Rmd:39,42`): tidak ada
  agregat (mean/median/min/maks) dari file mana pun yang menghasilkannya →
  kemungkinan run lama; **TIDAK PASTI**.
- Kapan, di mesin apa, dengan versi paket apa `eblup/seblup/steblup` dan
  `beta_estimates_n1000` dihasilkan: **TIDAK PASTI** (tidak tercatat, tidak ada skrip).
- Jumlah thread OpenMP yang dipakai saat benchmark: **TIDAK PASTI** (tidak diatur/dicatat).
- Apakah `hb_normal_benchmark.rds` dan `hb_benchmark.rds` memang sengaja identik
  (byte-identik, MD5 `dc15e7c1…3106`): keduanya ditulis oleh skrip yang sama,
  tetapi alasan pembuatannya ganda **TIDAK PASTI**.
- Apakah tabel/vignette dihitung dari versi file yang lebih lama (angka 8.69/408):
  **TIDAK PASTI**.
