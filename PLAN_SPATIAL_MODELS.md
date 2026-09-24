# Plan: Add Additional Spatial Models to `ebp_area()`

## Overview

Tambahkan model spasial baru pada fungsi `ebp_area()`:
1. **Besag Proper** (`besagproper`) - Proper CAR model

> **Note**: SLM model tidak ditambahkan karena INLA memiliki issue dengan cara pass parameter `rho.min`/`rho.max` dalam formula string. Bisa ditambahkan di masa depan dengan penelitian lebih lanjut.

---

## Status: ✅ COMPLETED

### Changes Made

#### 1. Modified `R/ebp_area.R`

- Updated dokumentasi untuk menambahkan `besagproper`
- Updated `match.arg()` untuk spatial
- Updated validasi W matrix
- Updated random_part switch statement

#### 2. Modified `R/ebp_area_utils.R`

- Updated `_convert_inla_to_fastsae()` untuk ekstraksi variance dari model `besagproper`

#### 3. Added tests in `tests/testthat/test_ebp_area.R`

- Tests untuk `besagproper` model
- Validation tests untuk semua model spasial
- Random effect variance tests

---

## Model Spasial yang Tersedia Sekarang

| Model | INLA Name | Keterangan |
|-------|-----------|------------|
| **none** | `iid` | No spatial random effects |
| **bym** | `bym` | Besag-York-Mollie model |
| **bym2** | `bym2` | Proper BYM2 model |
| **besagproper** | `besagproper` | Besag proper model (proper CAR) |
