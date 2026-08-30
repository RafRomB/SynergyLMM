### A Priori Power Analysis

Prospective, _a priori_, statistical power analysis allows to assess how the statistical power of the analysis varies with variations of some parameters such as the sample size, or drug combination effect size, while keeping the other parameters constant.

_A priori_ power analysis are only available for **exponential growth** models, and **Bliss** and **HSA** synergy reference models.

#### Sample Size Power Analysis

SynergyLMM web-app uses the estimated parameters from the input dataset and fitted model to calculate the _a priori_ estatistical power with different sample sizes.

#### Time Power Analysis


The number of measurement points per subject is another experimental factor that can influence the experimental design and statistical power. This depends on the duration of tumor growth follow-up and the frequency of the measurements during that period.

If `max` option is selected, the user can explore the _a priori_ power analysis with experimental durations. The estimated growth rates for the control and treatment groups are obtained from the input data and fitted model. The user can then specify different experimental variables, such as the sample size per group, the minimum and maximum time of follow-up, or the interval between each tumor measurement.

If `freq` option is selected, a similar analysis can be performed, but in this case varying the time interval between each tumor measurement, keeping constant the experiment duration.

#### Variability Power Analysis

A key determinant of statistical power is the effect size of the drug combination effect, given by the magnitude of the estimated coefficient for the growth rate after combination treatment, relative to that observed in the monotherapy and control groups. Treatments with larger effect sizes are easier to detect, even with smaller sample sizes,  resulting in higher statistical power.
Another critical factor is the variance components of the model, including residual variance (within-subject variance) and random effects variance (between-subject variance). High estimated coefficients for these variances reduce precision, thereby diminishing statistical power.

This section allows the user to change any of these factors. By default, the estimated values from the input datased and fitted model are used, but the user can vary these values to simulate any drug combination experiment, including varying the estimated growth rate for the control, Drug A, Drug B, and Combination groups, or the variance for the residuals and random effects.

**Note: Currently, _a priori_ power analysis is only available for 2-drugs combination studies. If the input dataset corresponds to a 3-drugs combination experiment, only the estimated coefficients for Drug A and Drug B will be used (Drug C will be ignored).**