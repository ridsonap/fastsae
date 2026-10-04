# Catatan teknis: model spasial & spasio-temporal frekuentis (`fastsae`)

Sumber (semua dibaca langsung, nomor baris = baris file):
`R/eblup_sfh.R`, `R/eblup_stfh.R`, `R/RcppExports.R`,
`src/eblup_sfh.h`, `src/eblup_sfh_core.cpp`, `src/eblup_sfh_pbmse.cpp`, `src/eblup_sfh_npbmse.cpp`,
`src/eblup_stfh_core.h`, `src/eblup_stfh_core.cpp`, `src/eblup_stfh_pbmse.cpp`, `src/Makevars`.
Pembanding `sae` 1.3 dibaca via `deparse(sae::eblupSTFH)` / `deparse(sae::pbmseSTFH)` (read-only).

---

## 1. `eblup_sfh()` — wrapper R (spatial Fay-Herriot)

### 1.1 Signature (`R/eblup_sfh.R:85-99`)
```r
eblup_sfh(formula, vardir, domain = NULL, data,
          method = c("REML","ML"), mse_method = c("analytical","pbmse","npbmse"),
          W = NULL, B = 100, n_threads = 1, seed = -1,
          maxiter = 100, precision = 1e-4, print_result = TRUE)
```
- `method` → `match.arg(toupper(...))` (`:100`), default **REML**; `mse_method` default **analytical** (`:101`).
- `domain` default `1:nrow(data)` (`:102-106`).
- Validasi: `vardir` length (`:113`), `anyNA(X)` → abort (`:121`), `vardir>0` hanya area tersampel (`:126-129`),
  `W` wajib & bujur sangkar `n_total x n_total` (**termasuk area tak-tersampel**, `:131-141`).
- `seed = -1` = jangan sentuh RNG R (dok `:30-31`; dieksekusi di C++ `src/eblup_sfh_pbmse.cpp:29`).
- `n_threads <= 0` = semua core (dok `:28-29`; dieksekusi `eblup_sfh_pbmse.cpp:73-75`).

### 1.2 Alur (`R/eblup_sfh.R`)
- `y`/`X` dari `model.frame(..., na.action=na.pass)` (`:109,117-118`) → `y` boleh berisi `NA`.
- `idx_s <- !is.na(y)` (`:144`).
- **pbmse/npbmse**: subset ke area tersampel `Xs,ys,vardirs,Ws = W[idx_s,idx_s]` (`:155-158`) lalu
  `.seblup_pbmse()` (`:161`) / `.seblup_npbmse()` (`:167`).
  Jika ada area tak-tersampel: fit ulang `.seblup_core(...)` penuh (`:179-182`), timpa baris tersampel
  dengan hasil bootstrap (`:187-189`), kolom bootstrap diisi `NA` untuk area tak-tersampel (`:195-200`).
- **analytical**: `.seblup_core()` langsung dengan `W` penuh (`:205-213`).
- Metadata: `estcoef` rownames=`colnames(X)`, `formula`, `model="SFH"`, `domain` jadi kolom pertama
  `df_eblup`, `call`, `data`, `W`, `class="fastsae"` (`:217-228`). Gagal konvergen → alert + return
  tanpa print (`:230-233`).

### 1.3 `@references` (disalin persis, `R/eblup_sfh.R:5-8`)
```
#' @references
#' \enumerate{
#'  \item Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley & Sons.
#' }
```
**Tidak ada DOI** di roxygen (grep `doi|@references` di `R/` — `R/eblup_sfh.R:5` satu-satunya untuk file ini).

---

## 2. `seblup_core()` — inti SFH (`src/eblup_sfh_core.cpp`)

### 2.1 Signature & return
- Export Rcpp: `[[Rcpp::export(.seblup_core)]]` (`:209`), definisi `:210-219`:
  `seblup_core(Xall, yall, vardirall, Wall, method="REML", maxiter=100, precision=1e-4, only_core=false)`.
  Wrapper R: `R/RcppExports.R:8-10`.
- Forward declaration & struct `SeblupFitArma{theta,g1d,g2d,sigma2,rho_fix,n_iter,converged}` di `src/eblup_sfh.h:9-17,32-41`.
- **`only_core=TRUE`** (dipakai bootstrap) → list ringan 13 elemen (`:431-445`):
  `convergence, sigma2_u, rho, beta, Xbeta, theta, GVi, g1d, g2d, Q, XtVi, Vi, G`.
- **Mode penuh** (`:655-667`): `estcoef` (data.frame `beta, std.error, stderr_beta, zvalue, pvalue` `:629-635`),
  `random_effect_var=sigma2`, `rho`, `estvarcomp` (parameter/estimate/std.error: `sigma2_u`,`rho` `:649-653`),
  `goodness` (`loglikelihood,AIC,BIC` `:472-476`), `df_eblup`
  (`y, vardir, eblup, random_effect, mse, rse` `:637-644`), `model="SFH"`, `level="area"`, `n_iter`, `convergence`, `method`.

### 2.2 Model persis seperti dikodekan
- **y**: `yall` (bisa `NA`); pemisahan area tersampel/tak-tersampel `:240-254` (`find_nonfinite`/`find_finite`).
- **Matriks W**: user-supplied, row-standardized (dok `R/eblup_sfh.R:16-21`); fungsi tidak membangun W.
  Rook contiguity hanya tersedia lewat pembuat data `sim_spatial_weights(D, type=c("knn","grid","ring"), style=c("B","W"), k=4, coords=NULL, seed=NULL)` (`R/sim_spatial_weights.R:47-52`; `type="grid"` = Rook `:8,30-32`; `style="W"` = row-standardized `:10-12`).
- **Efek spasial = SAR / spatial error (bukan CAR proper)**, dikodekan lewat *precision*:
  ```
  A      = (I - rho*W') (I - rho*W)                 // :311 (juga :63, :150, :397, :580)
  derSigma = A^{-1}                                 // :312-313 (inv_sympd, fallback pinv)
  G      = sigma2 * A^{-1}                          // :401 (:154 di seblup_fit_core_arma)
  V      = G + diag(vardir)                         // :318, :402
  ```
  ⇒ `Cov(u) = σ²[(I−ρW)'(I−ρW)]^{-1}` = `σ²(I−ρW)^{-1}(I−ρW')^{-1}` ⇒ `u = ρWu + ε` (simultaneous/conditional AR).
  Turunan: `dA/dρ = -WpWt + 2ρ·WtW` (`:315`, `WtW=W'W`, `WpWt=W+W'` di `:283-284`),
  `dV/dρ = -σ² A^{-1} (dA/dρ) A^{-1}` (`:316`) — persis `−σ²·derSigma·derRho·derSigma`.
- **X/Z**: tidak ada Z terpisah; `X` = `model.matrix` (`R/eblup_sfh.R:118`). `m = nrow(X)` (semua area), `p = ncol(X)` (`:259-260`).
- **EBLUP**: `eblup = Xβ + G V^{-1} (y − Xβ)` (`:462`, dan `:423` untuk only_core).
- **β**: `β = (X'V^{-1}X)^{-1} X'V^{-1} y` (`:407-414`); `Q = (X'V^{-1}X)^{-1}`, `P = V^{-1} − V^{-1}XQ X'V^{-1}` (`:496-497`, `:330`).

### 2.3 Matriks identitas / strategi memori
- **Tidak ada Woodbury, tidak ada invers blok** di SFH. Inversi penuh `inv_sympd(Vi, V)` per iterasi
  (`:319`, fallback `pinv` `:320`) dan `inv_sympd(derSigma, A)` (`:312`) → **O(m³) per iterasi**.
- Bukti alokasi **O(m²)**: `mat A(m,m), Vi(m,m), derSigma(m,m); mat derRho(m,m), derVRho(m,m), V(m,m); mat P(m,m)` (`:294-298`).
  `log_det(V)` penuh untuk loglik (`:465-466`) → O(m³) sekali di akhir.
- Prekomputasi di luar loop: `Wt, I, WtW, WpWt` (`:279-284`).
- `seblup_fit_core_arma()` adalah salinan loop yang sama tanpa objek R/RNG (`:11-16`, `:17-195`) — aman OpenMP.

### 2.4 Fisher scoring
- **Yang diiterasikan**: 2 parameter `theta = (sigma2_u, rho)` (`:291-293`).
- **Starting values**: `sigma2_u(0) = median(vardir)`, `rho(0) = 0.5` (`:288-289`; `:39-40`).
- **Loop**: `while (diff > precision && k < maxiter)` (`:307`), `diff = precision+1` awal (`:300`).
- **Skor** (REML, `:350-351`) / (ML, `:342-343`):
  - REML: `s_a = −½ tr(P·D_a) + ½ y' P D_a P y`
  - ML : `s_a = −½ tr(V^{-1} D_a) + ½ y' P D_a P y` ← **lihat D1 di §5**
- **Matriks informasi** (`:345-348` ML; `:353-356` REML): `I_ab = ½ tr(M_a M_b)` dengan `M_a = V^{-1}D_a` (ML) / `P D_a` (REML).
- **Update**: `theta_{k+1} = theta_k + I^{-1} s` (`:365-371`, `solve(step, Idev, s)`).
- **Boundary**:
  - `rho` di-clamp `(-0.999, 0.999)` tiap iterasi (`:373-374`); di `seblup_fit_core_arma` `:126-127`.
  - `sigma2` **tidak** di-clamp saat iterasi; hanya `if (sigma2<0) sigma2=0` di finalisasi (`:394`, `:147`).
  - Final `rho_fix`: **dipetakan ke ±1.0 persis** jika == ±0.999 (`:389-391`) — beda dengan `seblup_fit_core_arma` yang mempertahankan ±0.999 (`:142-144`). Lihat D4.
  - Gagal `solve`/non-finite → `k = maxiter; break` (`:365-369`, `:380-383`) ⇒ `convergence=FALSE` (`:416`).
- **Konvergensi**: `diff = max(|Δθ| / (|θ_k|+1e-12))` relatif (`:379`); sukses = `k < maxiter && diff <= precision` (`:416`).
  `seblup_fit_core_arma` menambah cek finiteness + `g1d>=0`, batas `1e6` (`:183-184`) yang **tidak** ada di mode penuh.

### 2.5 MSE analitik (g1…g4) — semua di `src/eblup_sfh_core.cpp`
- **g1, g2** (`:481-485`, duplikat di only_core `:425-429`, dan di `seblup_fit_core_arma:175-180`):
  ```cpp
  Ga = G - GVi*G;      // G - G V^{-1} G
  R  = X - GVi*X;
  g1d = Ga.diag();
  g2d = sum((R*Q) % R, 1);   // diag(R Q R')
  ```
- **g3** (`:487`, `:512-526`): matriks turunan `Amat = dV/dρ = -σ² A^{-1}(dA/dρ)A^{-1}` (`:491-494`),
  `l1 = V^{-1}D - σ² V^{-1}D V^{-1}D`, `l2 = V^{-1}Amat - σ² V^{-1}Amat V^{-1}D` (`:513-517`),
  lalu per-area `i`:
  ```cpp
  L(2,m): L.row(0)=l1.t().row(i); L.row(1)=l2.t().row(i);
  g3d(i) = trace(L * V * L.t() * Idevi);            // :521-526
  ```
  `Idev` = informasi REML 2×2 (`:502-506`), `Idevi = Idev^{-1}` (`:508-510`).
- **g4** (`:528-541`): `psi=diag(vardir)`; `D12aux=-der3`, `D22aux = 2σ² der3·rhosigma − 2σ² derSigma·WtW·derSigma`;
  `D = (psiVi·D12aux·Vipsi)(Idevi01+Idevi10) + (psiVi·D22aux·Vipsi)·Idevi11`; `g4d(i)=0.5*D(i,i)`.
- **Gabungan**: `mse = g1d + g2d + 2.0*g3d - g4d` (`:543`).
- **Koreksi bias ML** (`:546-571`, hanya jika `method=="ML"`): `h_a = -tr(Q X'V^{-1} D_a V^{-1}X)`, `bML = Idevi·h/2`,
  turunan g1 `dg1_dA`, `dg1_dp`, lalu `mse -= bML'·gradg1` (`:570`). **Hanya g1 dikoreksi.**
- **Area tak-tersampel** (`:576-604`): `A_full = [(I−ρW')'(I−ρW')]^{-1}` dari `Wall` penuh (`:580-583`),
  `KrigW = G_rs V^{-1}` (`:589`), `u_ns = KrigW·resid` (`:590`), `eblup_ns = X_ns β + u_ns` (`:591`),
  **MSE = g1_ns + g2_ns saja** (`:596-599`, `mse_ns = g1d_ns + g2d_ns`) — tanpa g3/g4.
- **RSE**: `100·sqrt(mse)/|eblup|`, `NA` jika `eblup==0` (`:621-624`).

### 2.6 Loglik/AIC/BIC
`loglike = −0.5(m log 2π + log|V| + r'V^{-1}r)` (`:465-468`); `AIC = −2ll+2(p+2)`, `BIC = −2ll+(p+2)log m` (`:469-470`).

---

## 3. Bootstrap SFH

### 3.1 Parametrik — `src/eblup_sfh_pbmse.cpp`
- Signature `:16-28`: `seblup_pbmse(X, y, vardir, W, method="REML", maxiter=100, precision=1e-4, B=100, n_threads=0, max_attempts_factor=5, seed=-1)`; export `[[Rcpp::export(.seblup_pbmse)]]` `:15`; wrapper `R/RcppExports.R:16-18`. `max_attempts_factor` **tidak** diekspos `eblup_sfh`.
- Seed: `if (seed>=0) set_r_seed(seed)` (`:29`); `set_r_seed` memanggil `base::set.seed` (`eblup_sfh_core.cpp:203-207`).
- Fit awal `only_core=TRUE` (`:34-38`); warning jika tak konvergen (`:40-42`).
- **RNG dibangkitkan SEKUENSIAL di luar region paralel** (`:90-98`): `U_boot ~ N(0,σ²)` (`rnorm`, `:94`), `E_boot(i,b) ~ N(0,vardir_i)` (`R::rnorm`, `:96`).
- **Generate `u*`**: `V_boot = solve(I − ρW, U_boot)` (`:57,100`) ⇒ `u* = (I−ρW)^{-1}e*`, kovarians `σ²A^{-1}` (eksak).
- **Paralel**: `#pragma omp parallel for schedule(dynamic)` `:108`; per-replicate:
  `theta_boot = Xβ̂ + v_boot` (`:111`), `y* = theta_boot + e*` (`:112`), **refit penuh Fisher-scoring** (`:114-116`),
  valid jika `converged && sigma2>=0 && theta/g1d/g2d finite` (`:118-124`).
  Proxy g3: `β_sblup = Q X'V^{-1} y*`, `θ_sblup = Xβ_sblup + G V^{-1}(y* − Xβ_sblup)` dengan **bobot dari fit awal** (`:127-128`).
- **Skema batch** sampai tepat B replikat valid (`:77-153`), anggaran `B*max_attempts_factor` (`:71,81-88`), stop jika habis (`:84-85`), warning bila 0 valid di satu round (`:147-152`).
- **Rumus** (`:158-164`):
  ```cpp
  mse_pb  = mean_b (θ̂*_b − θ*_b)^2                 // sum_mse/B
  mse_pbbc = 2*(g1+g2) − g1_pb − g2_pb + g3_pb      // g3_pb = mean (θ̂*_b − θ_sblup_b)^2
  ```
- Refit ulang `only_core=FALSE` untuk output lengkap (`:166-170`), kolom ditambah `mse_pb`, `mse_pbbc` (`:172-175`), plus `n_rounds`, `n_attempted_total` (`:177-178`).

### 3.2 Nonparametrik — `src/eblup_sfh_npbmse.cpp`
- Signature `:16-28` identik pbmse; export `:15`.
- **Hanya REML**: `stop` bila `method != "REML"` (`:32-34`).
- Residual terstandar (`:70-110`): `P = V^{-1} − V^{-1}X Q X'V^{-1}` (`:77`),
  `Ve = diag(vardir) P diag(vardir)` (`:73,79`), `Vu = (I−ρW) G P G (I−ρW)'` (`:65,80`),
  simetrisasi (`:83-84`), **eigen** `eig_sym` lalu ambil kolom `p..m-1` (buang p eigen kecil pertama, asumsi urut naik) → `Vei05`, `Vui05` (`:87-97`),
  `ustim = Vui05 (I−ρW)(G V^{-1} r)`, `estim = Vei05 (r − G V^{-1} r)` (`:99-100`),
  standardisasi mean/sd empiris, `u_std` diskalakan `sqrt(σ²)` (`:102-110`).
- **Resampling indeks SEKUENSIAL** (`:143-149`, `R::unif_rand`), **paralel** `#pragma omp parallel for schedule(dynamic)` `:159`.
- Per-replicate (`:161-192`): `e_boot = sqrt(vardir) %*% e_samp` (`:166`), `v_boot = (I−ρW)^{-1}u_boot` (`:167`),
  `θ* = Xβ̂+v_boot` (`:168`), refit penuh (`:171-173`), proxy sblup bobot-awal (`:184-185`).
- **Rumus** (`:216-222`): `mse_npb = mean(θ̂*−θ*)²`; `mse_npbbc = 2(g1+g2) − g1_npb − g2_npb + g3_npb`.
- Output tambahan `mse_npb`, `mse_npbbc` (`:230-233`), `n_rounds`, `n_attempted_total` (`:235-236`).

### 3.3 OpenMP
| File | Baris |
|---|---|
| `src/eblup_sfh_pbmse.cpp` | `:6-8` include, `:13` `// [[Rcpp::plugins(openmp)]]`, `:73-75` `omp_set_num_threads`, **`:108` `#pragma omp parallel for schedule(dynamic)`** |
| `src/eblup_sfh_npbmse.cpp` | `:6-8`, `:13`, `:125-127`, **`:159` `#pragma omp parallel for schedule(dynamic)`** |
| `src/eblup_stfh_pbmse.cpp` | `:6-8`, `:13`, `:345-347`, **`:352` `#pragma omp parallel for schedule(dynamic)`** |
| `src/Makevars` / `Makevars.win` | `:3` `PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)`; `:4` `PKG_LIBS` + LAPACK/BLAS/FLIBS |

Tidak ada pragma di `eblup_sfh_core.cpp` maupun `eblup_stfh_core.cpp` (hanya `eblup_tfh.cpp:481` di file lain).

---

## 4. `eblup_stfh()` — spatio-temporal Fay-Herriot

### 4.1 Signature (`R/eblup_stfh.R:103-122`)
```r
eblup_stfh(formula, vardir, data, domain, time, W,
           model = c("ST","S"), maxiter = 100, precision = 1e-4,
           sigma21_start = NULL, rho1_start = 0.5,
           sigma22_start = NULL, rho2_start = 0.5,
           compute_mse = FALSE, B = 100, n_threads = 1, seed = -1, print_result = TRUE)
```
- **Panel lengkap wajib**: `NA` pada respons → `cli_abort` (`:129-139`); dok `:64-67`.
- Validasi `vardir>0` (`:152`), dimensi `nrow(X) == D*T` (`:163-171`), `W` `D x D` tanpa NA (`:174-181`),
  `sigma*_start >= 0` (`:184-189`), `rho1_start ∈ (-1,1)` (`:190-192`), `rho2_start ∈ (-1,1)` **hanya bila `model=="ST"`** (`:193-195`).
- Panggil `.eblup_stfh_core(...)` (`:198-212`); `sigma*_start = NULL` dikirim sebagai **`-1`** (penanda "pakai default", `:208,210`).
- Metadata: rownames `estcoef`, `stderr_beta = std.error`, `zvalue = tvalue` (`:216-220`),
  `formula`, `level="area"`, **`method <- "REML"` (hard-coded, `:225`)**, `call`, `data`, `class="fastsae"` (`:223-230`).
- PBMSE opsional (`:241-266`): `.pbmse_stfh(...)` → `mse_pb`, `mse`, `rse` (`:263-265`), `res$B`.
  Tanpa MSE: `mse_pb/mse/rse = NA`, `B = NA_integer_` (`:269-272`).
- Urutan kolom `df_eblup`: `domain, time, y, eblup, vardir, mse, rse, mse_pb, random_effect_u1, random_effect_u2` (`:276-284`).

### 4.2 `@references` (disalin persis, `R/eblup_stfh.R:10-16`)
```
#' @references
#' \enumerate{
#'  \item Marhuenda, Y., Molina, I., & Morales, D. (2013). Small area estimation
#'    with spatio-temporal Fay-Herriot models. Computational Statistics & Data
#'    Analysis, 58, 308-325.
#'  \item Rao, J. N., & Molina, I. (2015). Small area estimation. John Wiley & Sons.
#' }
```
**Tanpa DOI.** Dokumen eksplisit menyatakan reimplementasi `sae::eblupSTFH()` (`:3-8`).

### 4.3 `eblup_stfh_core()` (`src/eblup_stfh_core.cpp`)
- Export `:160`, signature `:161-175`:
  `eblup_stfh_core(Xall,yall,vardirall,proxmat,D,Tt,model="ST",maxiter=100,precision=1e-4,
   sigma21_start=-1, rho1_start=0.5, sigma22_start=-1, rho2_start=0.5)`.
  Default start: `sigma2* = 0.5*median(vardir)` jika `<0` (`:201-203`); `nparam = 4` (ST) / `3` (S) (`:207-208`).
- **Return sukses** (`:527-535`): `estcoef` (`beta,std.error,tvalue,pvalue` `:488-493`),
  `estvarcomp` (`estimate,std.error` untuk `sigma21,rho1,sigma22[,rho2]` `:514-517`; nama `:219-221`),
  `goodness` (**`loglike,AIC,BIC`** `:477-479`), `df_eblup` (`eblup, random_effect_u1, random_effect_u2` `:520-524`),
  `model`, `convergence`, `n_iter`.
- **Return gagal** (`:236-246`): semua komponen `NULL`, `convergence`, `n_iter`.

#### Model persis seperti dikodekan
Kerangka `y_dt = x_dt'β + u1_d + u2_dt + e_dt`, `e_dt ~ N(0, vardir_dt)`, penulisan **domain-major** (dok `R/eblup_stfh.R:20-23`).
- **Spatial effect (bukan Knorr-Held type I–IV)**: `u1` **konstan sepanjang waktu dalam satu domain**,
  `Vu1 = σ²₁ · Ω₁` dengan `Ω₁ = [(I − ρ₁W)'(I − ρ₁W)]^{-1}` (`:74-80`);
  di ekspansi matriks: `Va[0] = kron(Ω₁, OnesT)` (`:281`, `OnesT` = matriks `T×T` semua satu, `:217`)
  ⇒ identik dengan `Z₁ Ω₁ Z₁'` pada `sae::eblupSTFH` (`Va[[1]] <- Z1 %*% Omega1rho1_k %*% tZ1`).
- **Temporal**: AR(1) per-domain `Ω₂(i,j) = ρ₂^{|i−j|} / (1−ρ₂²)` (`build_omega2` `:23-46`; diagonal `:43`),
  turunan `∂Ω₂/∂ρ₂` (`:36-39`, `:45`).
  Blok `A_d = σ²₂ Ω₂ + diag(vardir_d)` (`:108-109`); **model `"S"`**: `A_d = (σ²₂ + vardir_d)` diagonal, `Ω₂=I` (`:100-106`).
- **Matriks waktu/spasi**: `invA` blok-diagonal `M×M` (`:90-123`), `Z₁` = indikator domain (implisit lewat `build_invAZ1` `:128-136`).
- `nparam` turunan V (`:280-289`):
  `Va[0]=kron(Ω₁,1_T)` (∂/∂σ²₁), `Va[1]=kron(−σ²₁Ω₁(dA₁/dρ₁)Ω₁, 1_T)` (∂/∂ρ₁, `:276-277`),
  `Va[2]=kron(I_D,Ω₂)` atau `I_M` untuk model `"S"` (∂/∂σ²₂), `Va[3]=kron(I_D, σ²₂ ∂Ω₂/∂ρ₂)` (∂/∂ρ₂).

#### Inversi: Woodbury (bukti kode)
`src/eblup_stfh_core.cpp:139-153`:
```cpp
Cmat_out = mp.invVu1 + diagmat(mp.diagC);   // D x D  (= Ω1^{-1}/σ21 + diag(Z1' A^{-1} Z1))
inv_sympd(Cinv, Cmat_out);
invV_out = mp.invA - invAZ1 * Cinv * invAZ1.t();   // :151  Woodbury
```
- `invA` diisi per blok domain, inversi `T×T` per domain (`:108-121`) ⇒ **O(D·T³) + O(D³)**, bukan O(M³).
- `diagC(d) = accu(invAd)` (`:120`) = diagonal `Z₁'A^{-1}Z₁` (sah karena `invA` blok-diagonal per domain).
- `logdetV = logdetA + logdetVu1 + logdetC` (`:456-471`, `logdetVu1 = D log σ²₁ − logdet A₁` `:463`)
  = determinant lemma → **tanpa inversi/log-det M×M**.
- **Tetap O(M²) memori** di tempat lain: `invA` dan `invV` disimpan `M×M` (`:90, :151`), plus `Va[nparam]` dan `PV[nparam]` matriks `M×M` penuh (`:280-297`), `P` `M×M` (`:272`).
  Bukan `O(m^2)` vs `O(m^3)` untuk memori — keduanya O(·²); yang diturunkan adalah **biaya inversi**.

#### Fisher scoring
- Loop `while (diff > precision && k < maxiter)` (`:248`); **REML saja** (`P = invV − tXinvV' Q tXinvV` `:272`).
- Skor/informasi (`:310-316`):
  `S_a = −½ tr(P Va_a) + ½ y'P Va_a P y`; `F_ab = ½ tr(P Va_a P Va_b)` via `accu(PV[i] % PV[j].t())` (`:303`).
- Update `θ_{k+1} = θ_k + F^{-1} S` (`:319-323`).
- **Starting**: dari argumen (`:224-227`); default `0.5·median(vardir)`, `0.5` (`:201-203`).
- **Boundary**: `ρ₁,ρ₂` di-clamp `±0.999` (`:326-331`); `σ²₁,σ²₂` **tidak** di-clamp selama iterasi, hanya `max(·,0)` di final (`:346-347`).
- **Konvergensi**: ganti `0` dengan `1e-4` lalu `diff = max|(θ_k − θ_{k+1})/θ_k|` (`:334-336`).
- `k >= maxiter && diff >= precision` → `make_fail_result(false)` (`:339-341`); kegagalan inversi di tengah → return gagal (`:259,263,270,321`).
- **SE variance components** = `sqrt(diag(F^{-1}))` (`:496-512`), dengan guard `diagFinv<0` → return gagal (`:497-511`).
- **SE β** = `sqrt(diag(Q))`, `pvalue` dari **normal** `R::pnorm` (`:482-486`).
- EBLUP (`:430-453`): `u1_est = Vu1·(Z₁'V^{-1}r)` (`:435`), ekspansi konstan-waktu (`:436-440`),
  `u2` per blok `σ²₂Ω₂·(V^{-1}r)_d` (`:446-450`) / `σ²₂ V^{-1}r` untuk `"S"` (`:444`), `eblup = Xβ + u1 + u2` (`:453`).
- AIC/BIC: `−2ll + 2(p+nparam)`, `−2ll + log(M)(p+nparam)` (`:473-475`).

### 4.4 PBMSE STFH (`src/eblup_stfh_pbmse.cpp`)
- Signature `:233-246`: `pbmse_stfh(Xall,yall,vardirall,proxmat,D,Tt,model="ST",maxiter=100,precision=1e-4,B=100,n_threads=0,seed=-1)`; export `:232`.
- Seed `set.seed` (`:18-24`, `:247`).
- Fit awal memanggil `eblup_stfh_core` dengan start `0.5·med, 0.5, 0.5·med, 0.5` (`:257-262`); stop bila gagal (`:264-266`).
- Cholesky `Ω₁` (`:290-298`).
- **RNG SEKUENSIAL** (`:303-337`), **sebelum** pragma:
  - `u1 ~ N(0, σ²₁Ω₁)` via `sqrt(σ²₁)·L·z` (`:59-70`).
  - `u2` AR(1): `u2_0 ~ N(0, σ²₂/(1−ρ₂²))`, `u2_t = ρ₂ u2_{t−1} + N(0,σ²₂)` (`:75-87`) — ekuivalen dengan `sae`.
  - model `"S"`: `u2 ~ N(0,σ²₂)` i.i.d. per domain×waktu (`:322-331`).
  - `ε ~ N(0, vardir_i)` (`:334-336`).
- **Paralel** `#pragma omp parallel for schedule(dynamic)` (`:352`): `y_boot = Xβ̂ + u1* + u2* + ε*` (`:366`),
  lalu `fit_stfh_eblup_only()` (`:93-227`, dipanggil `:369-374`) — **variance components TIDAK diestimasi ulang**
  (komentar `:91` "Uses FIXED theta from initial fit - no re-estimation" dan `:368`); β dihitung ulang (`:187-191`),
  Woodbury dipakai lagi (`:160-184`). Gagal → isi `theta_est`, `valid_count=0` (`:376-382`).
- **MSE** (`:386-412`): `mudt_b = Xβ̂ + u1_b + u2_b` (`:404`, = θ* replikat by construction),
  `mse_pb = mean_b (eblup*_b − mudt_b)²` dibagi `n_valid` (bukan B) (`:405-412`).
  Komentar `:386` menyatakan "matching sae::pbmseSTFH formula".
- Return (`:417-430`): `eblup, mse_pb, B(=n_valid), B_total, sigma21, rho1, sigma22, rho2, beta, n_iter, convergence`.
- **Tidak ada** MSE analitik dan **tidak ada** nonparametrik bootstrap untuk STFH
  (tidak ada file `src/eblup_stfh_npbmse.cpp`; `R/eblup_stfh.R:269-272` mengisi `mse/rse = NA`).

---

## 5. DISCREPANCIES (implementasi vs rumus publik / vs dokumentasi)

| # | Baris kode | Masalah |
|---|---|---|
| **D1** | `eblup_sfh_core.cpp:330-331` (`P`), `:342-343` | **Skor ML** memakai `y'P D P y` dengan `P` = matriks REML (`V^{-1}−V^{-1}XQX'V^{-1}`), bukan `y'V^{-1}D V^{-1}y` seperti `sae`/rumus ML. Informasi ML (`:345-348`) benar memakai `V^{-1}`. Jadi ML tidak konsisten: kuadratik = REML, trace/informasi = ML. |
| **D2** | `:502-510`, `:647-648` | `std.error` σ²_u dan ρ dihitung dari `Idev` **bentuk REML** (`tr(PD·PD)`) **meskipun `method="ML"`**; `Idev` dihitung ulang setelah loop tanpa cabang ML. |
| **D3** | `:465-470` | `loglikelihood/AIC/BIC` selalu **bentuk ML** (`log|V|`, `r'V^{-1}r`) juga untuk fit REML; tidak ada koreksi `log|X'V^{-1}X|`. |
| **D4** | `:373-374` vs `:389-391` vs `:142-144` | Boundary ρ: selama iterasi di-clamp `±0.999`, tetapi `seblup_core` **final** memetakan `ρ=±0.999` → **`±1.0` persis**, di mana `A=(I−ρW')'(I−ρW)` singular → `inv_sympd` gagal → `pinv` (`:397-399`). `seblup_fit_core_arma` (jalur bootstrap) mempertahankan `±0.999`. Dua jalur estimasi berbeda perilaku di boundary. |
| **D5** | `:376` (loop), `:394` | `σ²_u` tidak pernah di-clamp nonnegatif saat iterasi (bisa negatif → `G` indefinited sementara); hanya `σ²=max(0,·)` di final. |
| **D6** | `R/eblup_sfh.R:39` vs `eblup_sfh_core.cpp:629-635` | Dok menyebut kolom **`tvalue`**; kode mengembalikan **`zvalue`** (dan `stderr_beta` duplikat `std.error`). |
| **D7** | `R/eblup_sfh.R:24-26,59-60` vs `eblup_sfh_core.cpp:596-599` | Dok menyebut "analytical MSE approximation" untuk area tak-tersampel; aktual **hanya `g1+g2`**, tanpa `g3`/`g4` (dan tanpa koreksi ML). |
| **D8** | `:257` vs `:580-587` | Inkonsistensi sub-blok W: fit area tersampel memakai `Wall(idx_s,idx_s)` (kovarians `Ω` dari `W_ss`), tetapi kriging memakai sub-blok dari `A_full^{-1}` (`G_full` dari `Wall` penuh). `[A_full^{-1}]_{ss} ≠ [(I−ρW_{ss})'(I−ρW_{ss})]^{-1}`. |
| **D9** | `R/eblup_sfh.R:167` vs `eblup_sfh_npbmse.cpp:32-34` | `mse_method="npbmse"` + `method="ML"` tidak dicek di R; baru `stop` di C++. |
| **D10** | `:183-184` vs `:183`(mode penuh tanpa cek) | Replikat bootstrap dinilai valid dengan kriteria lebih ketat (`g1d>=0`, batas `1e6`) daripada titik estimasi utama; `seblup_core` mode penuh tidak punya cek finiteness yang sama. |
| **D11** | `eblup_sfh_pbmse.cpp:164`, `eblup_sfh_npbmse.cpp:222` | Rumus bias-corrected `2(g1+g2) − ĝ1 − ĝ2 + ĝ3`, dengan `ĝ3 = mean(θ̂* − θ_sblup*)²` di mana `θ_sblup*` memakai **bobot `Q·X'V^{-1}` dan `G·V^{-1}` dari fit awal**, bukan dari replikat (`:127-128`, `:184-185`). Kesesuaian dengan rumus publik: **TIDAK PASTI** (tidak ada referensi yang mengikat rumus ini). |
| **D12** | `eblup_stfh_pbmse.cpp:91,368-374` vs `sae::pbmseSTFH` | **STFH bootstrap tidak mengestimasi ulang komponen varian tiap replikat** (fixed-θ), sedangkan `sae::pbmseSTFH` memanggil `eblupSTFH()` ulang tiap replikat (kode `sae` 1.3, baris deparse 102-105). Juga **tanpa** versi bias-corrected, berbeda dengan SFH. Klaim "matching sae formula" (`:386`) hanya berlaku untuk bentuk MSE `mean(θ̂*−θ*)²`, bukan untuk skema bootstrap-nya. |
| **D13** | `R/eblup_stfh.R:225` | Tidak ada pilihan ML; `method` di-hardcode `"REML"` — berbeda dari `eblup_sfh` yang punya `method=c("REML","ML")`. |
| **D14** | `R/eblup_stfh.R:193-195` vs `eblup_stfh_core.cpp:205` | Validasi `rho2_start ∈ (-1,1)` di R **hanya untuk `model=="ST"`**, tetapi C++ memvalidasi **tanpa syarat** → `model="S"` dengan `rho2_start=2` lolos R lalu error di C++. |
| **D15** | `eblup_sfh_core.cpp:473` vs `eblup_stfh_core.cpp:478` | Nama `goodness` berbeda: SFH `loglikelihood`, STFH `loglike`. |
| **D16** | `eblup_stfh_core.cpp:346-347,360-374` | Cabang `param_invalid` mengembalikan **`convergence = TRUE` dengan `estcoef = NULL`**; praktis tak tercapai karena `σ²` sudah di-clamp `max(0,·)` dan `ρ` sudah di-clamp `±0.999`, tapi struktur return-nya inkonsisten. |
| **D17** | `:319-323` vs `:496-512` | `F^{-1}` yang dipakai untuk `std.error` berasal dari iterasi terakhir dengan `θ_k`, sedangkan estimasi yang dilaporkan adalah `θ_{k+1}` → SE dievaluasi di titik sebelum update terakhir. |
| **D18** | `R/eblup_stfh.R:41-42,57,269-272` | STFH **tidak punya MSE analitik**: `mse`/`rse` = `NA` kecuali `compute_mse=TRUE`. Bertolak belakang dengan SFH yang selalu punya `mse` analitik. |
| **D19** | `eblup_stfh_core.cpp:281,436-440` | Struktur spasial = efek spasial **konstan-waktu** (`kron(Ω₁,1_T)`), bukan tipe interaksi Knorr-Held I–IV; `u2` = AR(1) per-domain. Konsisten dengan `sae::eblupSTFH` tetapi perlu ditulis eksplisit di makalah. |
| **D20** | `eblup_stfh_core.cpp:483-486` | Kolom bernama `tvalue` tetapi `pvalue` dihitung dari **distribusi normal** (`R::pnorm`), bukan t. SFH memakai `R::pnorm5` (`:456`) dan menamainya `zvalue`. |
| **D21** | `eblup_stfh_pbmse.cpp:410-412,420` + `R/eblup_stfh.R:266` | `mse_pb` dibagi `n_valid` (bukan `B`) dan `res$B` di-set ke `n_valid`, jadi bisa `< B` yang diminta; SFH sebaliknya memaksa tepat B replikat valid (batch, `eblup_sfh_pbmse.cpp:77-88`). |

---

## 6. Pseudocode bernomor (mengikuti alur kode)

### 6.1 `seblup_core` (mode penuh) — `src/eblup_sfh_core.cpp`
```
 1  validasi Wall square vs nrow(Xall)                         [:223-225]
 2  adaNA = yall.has_nan(); split idx_ns/idx_s; W = Wall[idx_s,idx_s]   [:240-257]
 3  siapkan Wt, I, WtW, WpWt (luar loop)                        [:279-284]
 4  sigma2_u[0]=median(vardir); rho[0]=0.5                     [:286-289]
 5  k=0; diff=precision+1                                      [:300-301]
 6  while diff>precision && k<maxiter:                         [:307]
 7    k++; A=(I-rho*W')'(I-rho*W); derSigma=A^{-1} (pinv fallback)  [:308-313]
 8    derRho=-WpWt+2*rho*WtW; derVRho=-sigma2*derSigma*derRho*derSigma [:315-316]
 9    V=sigma2*derSigma+diag(vardir); Vi=V^{-1}                 [:318-320]
10    XtVi=X'Vi; Q=(XtViX)^{-1}; P=Vi-XtVi' Q XtVi; Py=P y     [:323-331]
11    PD=P*derSigma; PR=P*derVRho                               [:334-335]
12    if ML: s=-.5tr(Vi D)+.5 y'PDPy ; I=.5tr(ViDa ViDb)        [:338-348]
13    else : s=-.5tr(P D)+.5 y'PDPy ; I=.5tr(PDa PDb)           [:349-357]
14    step = I^{-1} s; gagal -> k=maxiter, break                [:365-369]
15    rho di-clamp ke (-0.999,0.999)                            [:371-374]
16    diff = max|Δθ|/(|θ|+1e-12); non-finite -> k=maxiter,break [:376-383]
17  rho_fix: 0.999->1.0, -0.999->-1.0 ; sigma2=max(0,·)         [:389-394]
18  hitung ulang A^{-1}, G=sigma2*A^{-1}, V, Vi, Q, beta        [:396-414]
19  is_converged = (k<maxiter && diff<=precision)               [:416]
20  if only_core: return {convergence,sigma2,rho,beta,Xbeta,theta,GVi,g1d,g2d,Q,XtVi,Vi,G} [:419-446]
21  stderr/s.e./p (normal), eblup = Xb + GVi (y-Xb)             [:451-462]
22  loglike/AIC/BIC (bentuk ML)                                 [:464-476]
23  g1=diag(G-GVi G); g2=diag(R Q R'); g3=trace(L V L' Idevi); g4=0.5 D_ii [:481-541]
24  mse = g1+g2+2 g3 - g4 ; jika ML: mse -= bML' grad(g1)       [:543-571]
25  jika adaNA: kriging G_full, u_ns, mse_ns = g1_ns+g2_ns       [:576-604]
26  rse=100 sqrt(mse)/|eblup| ; susun List output               [:621-669]
```

### 6.2 `.seblup_pbmse` — `src/eblup_sfh_pbmse.cpp`
```
 1  if seed>=0 -> base::set.seed(seed)                          [:29]
 2  fit0 = seblup_core(only_core=TRUE); warning bila tak konv.  [:34-42]
 3  siapkan sigma2, rho, beta, Q, XtVi, GVi, g1, g2, I-rho*W    [:46-58]
 4  while n_valid < B (batch; anggaran B*max_attempts_factor)   [:77,71]
 5    bangkitkan U_boot (N(0,sigma2)), E_boot (N(0,vardir)) SEKUENSIAL [:91-98]
 6    V_boot = (I-rho*W)^{-1} U_boot                             [:100]
 7    #pragma omp parallel for schedule(dynamic)                [:108]
 8      theta*=X beta+v*; y*=theta*+e*                           [:111-112]
 9      refit Fisher-scoring penuh -> res (sigma2,rho,theta,g1,g2) [:114-116]
10      validasi res; θ_sblup = X(QX'V^{-1}y*)+GVi(y*-Xβ_sblup)  [:118-128]
11      simpan (θ̂*−θ*)², (θ̂*−θ_sblup)², g1*, g2*                [:130-134]
12    agregasi per round; n_valid+=round                          [:137-145]
13  mse_pb = mean(θ̂*−θ*)² ; mse_pbbc = 2(g1+g2)−g1*−g2*+g3*      [:158-164]
14  refit seblup_core(only_core=FALSE); tambah mse_pb/mse_pbbc    [:166-178]
```
`.seblup_npbmse` identik pada langkah 1-4 & 12-14, dengan langkah tambahan
**residual terstandar + eigen-dekomposisi (buang p eigen terkecil) + resampling indeks** (`eblup_sfh_npbmse.cpp:70-149`)
dan `stop` bila `method != "REML"` (`:32-34`).

### 6.3 `eblup_stfh_core` — `src/eblup_stfh_core.cpp`
```
 1  validasi model/D/T/W/NA (NA -> stop)                        [:176-194]
 2  sigma2*_start default 0.5*median(vardir); rho start validated [:201-205]
 3  θ0 = (sigma21, rho1, sigma22[, rho2]); nparam=4|3           [:223-227]
 4  while diff>precision && k<maxiter:                          [:248]
 5    Ω1=[(I-rho1 W)'(I-rho1 W)]^{-1}; Vu1=sigma21*Ω1           [:74-80]
 6    invA: per blok d, A_d = sigma22*Ω2 + diag(vardir_d) (T×T), atau 1/(sigma22+vardir) untuk "S" [:96-123]
 7    invV = invA - (invA Z1) C^{-1} (invA Z1)'   (Woodbury, C D×D) [:139-153]
 8    Q=(X'invV X)^{-1}; P=invV-X'invV Q invV'X; Py=P y         [:265-274]
 9    Va = {kron(Ω1,1_T), kron(-s21 Ω1 Ω1'ρ Ω1,1_T), kron(I,Ω2)|I_M, kron(I,s22 Ω2')} [:276-289]
10    PV_i=P Va_i; trPV_i; trPVPV_ij=accu(PV_i ⊙ PV_j')         [:292-307]
11    S_a=-.5 trPV_a+.5 Py'Va_a Py ; F_ab=.5 trPVPV_ab          [:310-316]
12    θ_{k+1}=θ_k+F^{-1}S ; rho1,rho2 di-clamp ±0.999           [:319-331]
13    diff = max|(θ_k−θ_{k+1})/θ_k| dengan 0 diganti 1e-4        [:334-336]
14  gagal konvergen -> return {semua NULL, convergence=FALSE}    [:339-341]
15  σ²=max(0,·); rebuild pieces; beta=Q X'invV y; resid          [:346-427]
16  u1 = Vu1 (Z1' invV resid), ekspasi konstan-waktu; u2 per blok; eblup = Xb+u1+u2 [:430-453]
17  logdetV = logdetA + [D log σ21 − logdet A1] + logdetC ; loglike/AIC/BIC [:456-479]
18  se(beta)=sqrt(diag Q); se(theta)=sqrt(diag F^{-1})           [:482-512]
19  return {estcoef, estvarcomp, goodness, df_eblup, model, convergence=TRUE, n_iter} [:527-535]
```

### 6.4 `.pbmse_stfh` — `src/eblup_stfh_pbmse.cpp`
```
 1  if seed>=0 -> set.seed                                       [:247]
 2  init = eblup_stfh_core(start default 0.5med/0.5)             [:257-266]
 3  ekstrak sigma21, rho1, sigma22, rho2, beta, theta_est=X*beta [:269-280]
 4  Cholesky Ω1                                                  [:290-298]
 5  SEKUENSIAL: u1~N(0,s21 Ω1); u2 AR(1) per domain; eps~N(0,vardir) [:303-337]
 6  #pragma omp parallel for schedule(dynamic)                  [:352]
 7    y* = X*beta + u1* (ekspansi) + u2* + eps*                  [:358-366]
 8    fit_stfh_eblup_only(y*) : beta dihitung ulang, σ²/ρ TETAP  [:369-374, :93-227]
 9    gagal -> Eblup*=theta_est, valid=0                         [:376-382]
10  mse_pb = Σ_valid (eblup* − (X*beta+u1*+u2*))² / n_valid      [:391-412]
11  return {eblup, mse_pb, B=n_valid, B_total, σ²/ρ, beta, n_iter, convergence} [:417-430]
```
