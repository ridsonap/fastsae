# The R Journal Submission Bundle: fastsae

This directory contains the complete submission bundle for **The R Journal**, prepared in strict accordance with the official submission guidelines at [https://journal.r-project.org/submissions.html](https://journal.r-project.org/submissions.html).

---

## 1. Article Metadata

- **Title**: *fastsae: High-Performance Frequentist and Bayesian Small Area Estimation in R*
- **Category**: Add-on Package / Software Article
- **Target Journal**: *The R Journal*
- **Authors**:
  1. **Ridson Al Farizal P.** (Department of Computational Statistics, Politeknik Statistika STIS, Jakarta, Indonesia; ORCID: [0000-0003-0617-0214](https://orcid.org/0000-0003-0617-0214); Email: `alfrzlp@gmail.com`)
  2. **Azka Ubaidillah** (Department of Statistics, Politeknik Statistika STIS, Jakarta, Indonesia; ORCID: [0000-0002-3597-0459](https://orcid.org/0000-0002-3597-0459); Email: `azka@stis.ac.id`)
- **Keywords**: small area estimation, Fay-Herriot, empirical best prediction, INLA, OpenMP, Rcpp, hierarchical benchmarking, spatial statistics, R
- **Page Count**: 15 pages (standard R Journal layout)

---

## 2. Directory Contents

| File | Description |
|---|---|
| `fastsae.Rmd` | Modern R Markdown manuscript designed for `rjtools` (`rjtools::rjournal_web_article` and `rjtools::rjournal_pdf_article`). Contains reproducible knitr chunks. |
| `RJwrapper.tex` | Master LaTeX wrapper for classic R Journal compilation. |
| `fastsae.tex` | LaTeX article body with complete sections, equations, tables, and R Journal typography. |
| `RJournal.sty` | Official R Journal LaTeX style package. |
| `fastsae.bib` | Comprehensive BibTeX database (45+ entries) with complete DOIs, volume, and issue numbers. |
| `RJwrapper.pdf` | Pre-compiled 15-page publication-quality PDF (0 errors, 0 undefined citations). |
| `figures/` | High-resolution vector figures in PDF format: |
|   * `fig-time-fh.pdf` | Runtime scaling across domain sample sizes (D) for frequentist FH models. |
|   * `fig-time-beta.pdf` | Runtime scaling across domain sample sizes (D) for Bayesian Beta models. |
|   * `fig-mem.pdf` | Peak RAM memory footprint across domain sizes (D) on log-log scale. |
|   * `fig-ratio.pdf` | Speedup ratios (364x FH, 80x SFH, 544x Spatio-temporal Beta). |
|   * `fig-equiv.pdf` | Numerical equivalence between `fastsae` (INLA) and `tipsae` (Stan HMC) (r = 0.9993). |
|   * `fig-diagnostics.pdf` | Calibration and precision diagnostic panel generated from `autoplot()`. |

---

## 3. How to Build and Verify

### Option A: Using the Modern `rjtools` Workflow (Recommended by The R Journal)
```r
# Install rjtools if not already installed
install.packages("rjtools")

# Check and validate the submission
rjtools::check_article("paper/rjournal/fastsae.Rmd")

# Render HTML and PDF versions
rmarkdown::render("paper/rjournal/fastsae.Rmd")
```

### Option B: Using the Traditional LaTeX Workflow
```bash
cd paper/rjournal
pdflatex RJwrapper.tex
bibtex RJwrapper
pdflatex RJwrapper.tex
pdflatex RJwrapper.tex
```

---

## 4. Compliance Checklist (R Journal Guidelines)

- [x] **Clear Motivation & Gap**: Identifies the triple operational barrier in official SAE: computational bottlenecks in bootstrap, MCMC friction in non-Gaussian likelihoods, and lack of two-stage hierarchical benchmarking.
- [x] **Comprehensive Comparison Table**: Table 1 explicitly contrasts `fastsae` against `sae`, `emdi`, `tipsae`, `mcmcsae`, and `SUMMER` across 12 capability dimensions.
- [x] **Rigorous Mathematical Methodology**: Exact equations for FH, SFH, STFH, BHF, INLA Latent Gaussian Models with exact survey dispersion matching (Beta, Gamma, Binomial, Poisson, NegBin), two-stage benchmarking, and RSE reliability classification.
- [x] **Reproducible Worked Examples**: Section 5 provides end-to-end executable code demonstrating data preparation, model fitting, model comparison (`compare_sae`), hierarchical calibration (`benchmark_sae`), reliability tagging (`add_reliability_flags`), and diagnostic plotting (`autoplot`).
- [x] **Complete Benchmark & Simulation**: Rigorous runtime/memory benchmarks (D from 30 to 1000) and design-based Monte Carlo simulations under known ground truth (D=500, R=25, N=250,000) demonstrating 1.92x efficiency gain over direct survey estimator.
- [x] **Institutional Dissemination Impact**: Table 5 provides a synthesis of empirical evidence and concrete policy implications ("So What").
- [x] **Clean Bibliography**: All 45+ citations include valid DOIs or URLs.
- [x] **Compilation Verified**: Compiled to `RJwrapper.pdf` (15 pages) with zero LaTeX errors and zero missing citations.

