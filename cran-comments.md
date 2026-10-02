## Test environments
* local macOS, R 4.4 / 4.5
* win-builder (devel and release)
* GitHub Actions (macOS, Windows, Ubuntu devel/release/oldrel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Submission summary
This is an update from version 0.1.0 to 1.0.0.

* Updated package title to "Frequentist and Bayesian Small Area Estimation in R" to reflect the inclusion of both Frequentist and Bayesian methods.
* Added Hierarchical Bayes Small Area Estimation models via INLA (`hb_area`, `hb_unit`, `hb_twofold`).
* Added Two-Fold subarea EBLUP models (`eblup_twofold`).
* Added model benchmarking (`benchmark_sae`), diagnostics (`diagnose`), comparison (`compare_sae`), mapping (`map_sae`), and export (`export_sae`).
