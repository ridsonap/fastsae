# R Journal submission: fastsae

Article: *fastsae: Frequentist and Bayesian Small Area Estimation in R*, by Ridson Al Farizal P.
and Azka Ubaidillah (Politeknik Statistika STIS). Submission guidelines:
<https://journal.r-project.org/submissions.html>.

## Files

| File | Content |
|---|---|
| `fastsae.Rmd` | The article (`rjtools::rjournal_article`); renders to HTML and PDF |
| `fastsae.bib` | References |
| `motivation-letter/` | Motivation letter |
| `RJournal.sty`, `RJwrapper.tex` | Style and wrapper used by the PDF output |

The article computes its tables and figures from the stored benchmark results, so it
renders in a few minutes:

| Result | Script | Stored in |
|---|---|---|
| Agreement with sae and emdi | `benchmarks/validation/agreement_reference.R` | `inst/extdata/agreement_reference.rds` |
| Frequentist timings (fastsae, sae, emdi on the same data) | `benchmarks/bench_frequentist.R` (calls `bench_frequentist_child.R`) | `inst/extdata/{eblup,seblup,steblup}_benchmark.rds`, `frequentist_benchmark_meta.rds` |
| Beta timings | `benchmarks/bench_beta_tipsae.R`, `bench_spatial_beta.R`, `bench_spatio_temporal_beta.R` | `inst/extdata/{beta,spatial_beta,spatio_temporal_beta}_benchmark.rds` |
| fastsae vs tipsae estimates | `benchmarks/bench_beta_concordance.R` | `inst/extdata/beta_concordance_n1000.rds` |
| Estimation engines | `benchmarks/bench_hb_area_engines.R` | `inst/extdata/hb_engine_benchmark.rds` |
| Gaussian HB, fastsae vs hbsae | `benchmarks/bench_hb_normal.R` | `inst/extdata/hb_normal_benchmark.rds` |
| Design-based simulation | `benchmarks/sim_mc_finite_pop_d500.R` | `benchmarks/output/mc_design_based_D500_R200.rds` |

The frequentist timings run every cell in its own R process; run them on an otherwise idle
machine (about two hours, most of it `sae::eblupSTFH()` at 500 and 1,000 areas). The
design-based simulation takes about an hour with 5 workers of 2 cores each.

## Build and check

Run from this directory (`rjtools::initial_check_article()` reads stdin through pandoc, so
give it an empty stdin when run non-interactively):

```sh
Rscript -e 'rmarkdown::render("fastsae.Rmd", output_format = "rjtools::rjournal_web_article")'
Rscript -e 'rmarkdown::render("fastsae.Rmd", output_format = "rjtools::rjournal_pdf_article")'
Rscript -e 'rjtools::initial_check_article(".")' < /dev/null
```

The article needs fastsae, INLA, tipsae (for the example data), sf, spdep, ggplot2,
patchwork and kableExtra.

`initial_check_article()` of rjtools 1.0.18.1 stops in `check_spelling()`, which looks for a
`\section*{References}` that the rjtools template does not produce; run the checks after it
one by one. It also flags the lower-case package name in the title and the proper noun in
"Example: poverty in Emilia-Romagna"; both are intended.

## Before submitting

- fastsae 1.0.0 must be on CRAN (a requirement of the journal); CRAN has 0.1.0, which lacks
  the HB functions the article describes.
- Set `date` in the YAML header to the submission date (`check_date()` compares it with the
  current date).
- When the engine benchmark ended, other processes (Finder, Spotlight) used about one CPU
  core; re-run `benchmarks/bench_hb_area_engines.R` on an idle machine for final numbers.
