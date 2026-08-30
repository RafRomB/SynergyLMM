## Influential Diagnostics

#### Cook's Distances

The identification of influential subjects is based on the calculation of Cook's distances. The Cook's distances
can be calculated based on the change in fitted values or fixed effects.

Use Cook's distances **based on the change in fitted values** to assess the influence of each subject in the model fit. You can use Cook's distance **based on the change in fixed effect values** if your main concern is bias in the estimated parameters.

- **Cook's distances based on the change in fitted values** 

The Cook's distances
are calculated as the normalized change in fitted response values due to the removal of a subject from the model.
Firts, a leave-one-subject-out model is fitted, removing individually each subject to fit the model. Then, the Cook's
distance for subject $i$, $(D_i)$, is calculated as:
$$D_i=\frac{\sum_{j=1}^n\Bigl(\hat{y}_{j}-\hat{y}_{j_{(-i)}}\Bigl)^2}{p\cdot MSE}$$

where $\hat{y}_j$ is the $j^{th}$ fitted response value using the complete model, and $\hat{y}_{j_{(-i)}}$ is the 
$j^{th}$ fitted response value obtained using the model where subject $i$ has been removed.

The denominator of the expression is equal to the number of the fixed-effects coefficients and the Cook distance is normalized by the
mean square error $(MSE)$ of the model. 

- **Cook's distances based on the change in fixed effects values**

The identification of the subjects with a greater influence in the estimated fixed effects is based on the calculation of Cook's distances, as
described in Gałecki and Burzykowsk (2013). To compute the Cook's distance for the fixed effect estimates (i.e., the contribution to each subject to the coefficients of its treatment group), 
first a matrix containing the leave-one-subject-out estimates or the fixed effects is calculated. Then, the Cook's distances are calculated according to:

$$D_i \equiv  \frac{(\hat{\beta} - \hat{\beta}_{(-i)})[\widehat{Var(\hat{\beta})}]^{-1}(\hat{\beta} - \hat{\beta}_{(-i)})}{p}$$

where $\beta$ represents the vector of fixed effects and $\hat{\beta}_{(-i)}$ is the estimate of the parameter vector \eqn{\beta} obtained by fitting the model to the data with the $i$-th subject excluded. The denominator of 
the expression is equal to the number of the fixed-effects coefficients.


By default, the threshold to highlight influential subjects is set to three times the mean of the Cook's distances values of all individuals. This can be modified in the _Advanced Options_ menu.

![Cook's Distances](cooksD.png){width="100%"}
_Example Cook's distances plot highliting three subjects with a high influence in the coefficient for the group._

