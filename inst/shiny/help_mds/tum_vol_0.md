## Null Tumor Measurements

This parameter allows to control the behavior regarding measurements in which the tumor measurement is 0, and therefore the logarithmic transformation is not possible.

If set to `ignore` (the default behavior), these measurements will be excluded from the analysis. 

If set to `transform`, one unit will be added to all measurements (e.g., measurement $X$ will be transformed in $X'=X+1$), before the $\text{log}$ transformation. 