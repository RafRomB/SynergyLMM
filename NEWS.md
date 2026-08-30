# SynergyLMM (development version)

# SynergyLMM 1.2.0

The 'SynergyLMM' Shiny application is now shipped with the package.

## New features

* New `runSynergyLMMApp()` launches the interactive 'SynergyLMM' application
locally. The app is bundled under `inst/shiny` and covers the full workflow:
data upload, model fitting, synergy evaluation, model diagnostics, and the
_post hoc_ and _a priori_ power analyses. The packages it needs are listed in
`Suggests`, so they are not a hard dependency; the function reports every
missing one in a single message. The app is still available online at
<https://synergylmm.uiocloud.no/>.
* `PostHocPwr()` gains an optional `progress` argument, a function of no
arguments called once per simulation. It allows front-ends to report progress
during long runs and is ignored when `NULL` (the default).
* `APrioriPwr()`, `PwrSampleSize()` and `PwrTime()` now attach the plot they
draw to their returned object as an attribute named `"plot"`. It can be
retrieved with `attr(result, "plot")` and passed to, for example,
`ggplot2::ggsave()`. The returned values themselves are unchanged.

## Fixes

* The application previously ran its own copy of the analysis code, which had
drifted from the package: it fitted the _a priori_ power models with an
intercept (`~ Time:Treatment` rather than `~ 0 + Time:Treatment`), used the
old `subject` identifier column, and did not support the variance-function
argument. The app now calls the package functions directly, so the two can no
longer diverge. Reported power is unaffected by this correction, as the
synergy contrasts are zero-sum and therefore invariant to the intercept
parametrisation: `APrioriPwr()`, `PwrSampleSize()` and `PwrTime()` agree with
the previous application output to within rounding.

# SynergyLMM 1.1.4

Patch update fixing a test failure reported by the CRAN additional checks with OpenBLAS.

## Fixes

* Tests no longer assert simulated outcomes. Whether a p-value is exactly 0
depends on the draws from `MASS::mvrnorm()`, which are not reproducible across
BLAS/LAPACK implementations. The wording of the corresponding `lmmSynergy()`
warning is now covered by a deterministic unit test.
* The reported resolution of simulated p-values is now computed as `1/nsim`.
It is unchanged for any `nsim` that is a power of ten, including all defaults.
* `plot_lmmSynergy()` now flags p-values approximated to 0 whenever at least one
is present, and draws them in the colour announced in the legend.
* Replaced the use of `.data$` inside `dplyr::select()`, deprecated since
tidyselect 1.2.0, with column names.

# SynergyLMM 1.1.3

Patch update to update the package maintainer's email address.

# SynergyLMM 1.1.2

Patch update to the SynergyLMM package (1.1.2) to add the citation to the published article and some minor updates.

## Fixes

* Added citation to published article: https://doi.org/10.1038/s41467-025-65218-9.
* `APrioriPwr()`, `PwrSampleSize()`, and `PwrTime()` now allow to control the output of exemplary data plot.
* Minor corrections in `lmmModel_estimates()` documentation and naming of columns in the resulting data frame.

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
