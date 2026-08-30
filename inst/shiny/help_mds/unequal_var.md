## Define Unequal Variances

Checking this option allows to define a different variance for each subject, time point, or treatment. Otherwise, the model assumes that the variance is the same for all subjects. 

Allowing for different variance per time point, subject, or treatment may improve the model diagnostics, and therefore, ensure reliable results.

If the model diagnostics indicate some violations of the distributional assumptions, it is recommended trying to define the variance structure to improve the model. 

You can use the diagnostics plots to identify if there is some heteroscedasticity of the residuals by sample, time point, or treatment, that may need to be specified in the variance structure (see examples in **Tutorial tab**).