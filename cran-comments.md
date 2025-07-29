## CRAN submission

This is a minor update to the SynergyLMM package (1.1.0)

## Changes

- Gompertz growth model is now available for fitting data using non-linear mixed 
effect models.
- `lmmSynergy()` now allows for multiple testing correction and p-value 
adjustment.
- `lmmModel_estimates()` now also reports the standard error of the 
fixed effect coefficients.
- `lmmSynergy()` now also returns a data frame with the estimated coefficients, 
calculated with `lmmModel_estimates()`, for each time point.
- `CookDistance()` now allows to choose between calculating Cook's distances 
based on changes of the fitted values, or changes of the fixed effects.
- Correction of combination index values in `lmmSynergy()` when `method = "RA"`.
- Updated documentation for clarity and consistency.
- Updated vignette with updates and precomputed results to reduce building time.


## R CMD check results

0 errors | 0 warnings | 1 note
