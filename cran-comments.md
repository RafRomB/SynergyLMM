## Resubmission

This submission fixes CRAN check additional issues found when re-building of vignette outputs that were occurring on platforms with no long double support (`noLD`) and systems using OpenBLAS. I have made the following changes:

* Removed lines 231-235 in vignettes/SynergyLMM.Rmd to avoid convergence problems in model fitting.

* Added 'lme4' package to Suggests in DESCRIPTION since 'performance' and 'clubSandwich' packages have it as suggested package and it is required for some of the functions used in this package.

## R CMD check results

0 errors | 0 warnings | 1 note
