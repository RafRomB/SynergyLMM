# Observed vs predicted ----

#' @title Observed vs predicted values and performance of the model
#' 
#' @description
#' `ObsvsPred` allows the user to have a straight forward idea about how the model is fitting the data, providing
#' plots of the predicted regression lines versus the actual data points.
#' 
#' 
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param nrow Number of rows of the layout to organize the observed vs predicted plots.
#' @param ncol Number of columns of the layout to organize the observed vs predicted plots.
#' @param ... Additional arguments to be passed to [performance::model_performance()].
#' @details
#'  The function provides visual and quantitative information about the performance of the model:
#' - A layout of the observed and predicted values of \eqn{log}(relative tumor volume) vs Time for each SampleID (i.e. subject), 
#' with the actual measurements, the regression line for each SampleID, and the marginal, treatment-specific, 
#' regression line for each treatment group.
#' - Performance metrics of the model obtain calculated using [performance::model_performance()]. The maximum likelihood-based Akaike's Information Criterion (AIC), 
#' small sample AIC (AICc), and Bayesian Information Criterion, and the Nakagawa's r-squared 
#' root mean squared error (RMSE) of the residuals, and the standard deviation of the residuals (sigma) are provided.
#' 
#' @returns A layout of the observed vs predicted values for each SampleID and model performance metrics.
#' 
#' @references
#' - Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' - Lüdecke et al., (2021). _performance: An R Package for Assessment, Comparison and Testing of Statistical Models_. Journal of Open Source Software, 6(60), 3139. https://doi.org/10.21105/joss.03139
#' - Sakamoto, Y., M. Ishiguro, and G. Kitagawa. 1984. _Akaike Information Criterion Statistics_. Mathematics and Its Applications. Reidel.
#' - Nakagawa, Shinichi, and Holger Schielzeth. 2013. _A General and Simple Method for Obtaining r2 from Generalized Linear Mixed-effects Models_. Methods in Ecology and Evolution 4 (February): 133–42. https://doi.org/10.1111/j.2041-210x.2012.00261.x.
#' - Johnson, Paul C. D. 2014. _Extension of Nakagawa & Schielzeth’s r 2 GLMM to Random Slopes Models_. Methods in Ecology and Evolution 5 (September): 944–46. https://doi.org/10.1111/2041-210X.12225.
#' - Nakagawa, Shinichi, Paul C. D. Johnson, and Holger Schielzeth. 2017. _The Coefficient of Determination r2 and Intra-Class Correlation Coefficient from Generalized Linear Mixed-Effects Models Revisited and Expanded_. Journal of The Royal Society Interface 14 (September): 20170213. https://doi.org/10.1098/rsif.2017.0213.
#' @examples
#' # Load the example data
#' data(grwth_data)
#' # Fit the model
#' lmm <- lmmModel(
#'   data = grwth_data,
#'   sample_id = "subject",
#'   time = "Time",
#'   treatment = "Treatment",
#'   tumor_vol = "TumorVolume",
#'   trt_control = "Control",
#'   drug_a = "DrugA",
#'   drug_b = "DrugB",
#'   combination = "Combination"
#'   )
#'# Obtain Observed vs Predicted plots, and model performance metrics
#' ObsvsPred(model = lmm, nrow = 4, ncol = 8)
#' @export
ObsvsPred <- function(model, nrow = 4, ncol = 5, ...) {
  print(performance::model_performance(model, metrics = c("AIC", "AICc", "BIC", "R2", "RMSE", "SIGMA")), ...)
  plot_ObsvsPred(model, nrow, ncol)
}
