# SynergyLMM 1.1.1

Patch update to the SynergyLMM package (1.1.1) so it does not break with the new ggplot2 v4.0.0 update.

## Fixes

Replaced fragile test for correct ggplot output in plot_SynergyLMM() function.

# SynergyLMM 1.1.0

## New features:

* Gompertz growth model is now available for fitting data using non-linear mixed 
effect models.
* `lmmSynergy()` now allows for multiple testing correction and p-value 
adjustment.

## Breaking changes:

* `lmmModel_estimates()` now also reports the standard error of the 
fixed effect coefficients.
* `lmmSynergy()` now also returns a data frame with the estimated coefficients, 
calculated with `lmmModel_estimates()`, for each time point.
* `CookDistance()` now allows to choose between calculating Cook's distances 
based on changes of the fitted values, or changes of the fixed effects.

## Minor improvements and fixes:

* Correction of combination index values in `lmmSynergy()` when `method = "RA"` 
and negative values in the denominator appear. A message is prompted
and CI values are capped at 100 to maintain interpretability. 

# SynergyLMM 1.0.1

Correcting minor errors in vignette.

# SynergyLMM 1.0.0

## Initial Release
We are excited to announce the release of **SynergyLMM**, an R package providing a comprehensive statistical framework for designing and analyzing _in vivo_ drug combination experiments.
