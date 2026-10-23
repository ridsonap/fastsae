---
output: pdf_document
fontsize: 12pt
---

\thispagestyle{empty}
\today

Editor   
The R Journal  
\bigskip

Dear Editors,
\bigskip

Please consider our article titled "fastsae: Frequentist and Bayesian Small Area Estimation in R" for publication in The R Journal.

The article introduces the R package `fastsae`, which brings the main area-level and unit-level small area estimation models used in official statistics together in one interface. Its frequentist part implements the Fay-Herriot model, its spatial and spatio-temporal extensions, and the unit-level and two-fold nested error models with analytical and bootstrap mean squared errors in compiled code that scales to thousands of areas. Its hierarchical Bayes part fits area-level, unit-level and two-fold models with six likelihoods and spatial, temporal and space-time effects, either by integrated nested Laplace approximations or by Hamiltonian Monte Carlo in JAX through the same call, with the dispersion of each area fixed by the sampling variance of its direct estimate and with documented default priors. Benchmarking, reliability classification, model comparison, maps and export complete the workflow.

The article compares the package with sae, emdi and tipsae in runtime and in the accuracy and coverage of the estimates, reports a design-based simulation, and illustrates a complete workflow with poverty data. The benchmark scripts and their stored results are in the package repository, and the article renders from the stored results in a few minutes; the timings of sae and emdi come from an earlier run, as stated in the article.

We believe the article will interest researchers and practitioners in official statistics and survey methodology who need small area estimates in production.

\bigskip
\bigskip

Regards,

Ridson Al Farizal P.  
Department of Computational Statistics  
Politeknik Statistika STIS  
Jl. Otto Iskandardinata No. 64C  
Jakarta Timur 13330, Indonesia  
alfrzlp@gmail.com

Azka Ubaidillah  
Department of Statistics  
Politeknik Statistika STIS  
Jl. Otto Iskandardinata No. 64C  
Jakarta Timur 13330, Indonesia  
azka@stis.ac.id
