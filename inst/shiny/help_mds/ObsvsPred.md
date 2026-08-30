### Indices of Model Performance

Some indices of the model performance are provided:
- Akaike’s Information Criterion (AIC)
- Small sample AIC (AICc)
- Bayesian Information Criterion (BIC)
- [Nakagawa’s $R^2$](https://doi.org/10.1098/rsif.2017.0213)
- Root Mean Squared Error (RMSE) of the residuals
- Standard Deviation (sigma) of the residuals

Smaller values of AIC, AICc, BIC, RMSE, and sigma of the residuals indicate a better-fitting model.

$R^2$ values range from 0 to 1. An $R^2$ value closer to 1 suggests that the model explains a larger proportion of the variance in the data, indicating a better fit.

_Note: $R^2$ values are only available for models fitted using the exponential growth model._

### Plots of Observed vs Predicted Values

In order to assess the adequacy of the model to the data and the different assumptions, such as the exponential tumor growth, plots of the observed versus predicted values of $log$(relative tumor volume) vs Time for each subject are provided. 

The actual measurements are shown as blue dots, the dashed orange line indicates the regression line for each subject, while the continous blue line indicates the marginal, treatment-specific, regression line for each treatment group.