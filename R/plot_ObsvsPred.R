#' @import graphics
NULL

#' @title Plots of Observed vs Predicted Values
#' @description
#' Visualization of observed vs predicted values by a fitted linear mixed model of tumor growth data.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param nrow Number of rows of the layout to organize the observed vs predicted plots.
#' @param ncol Number of columns of the layout to organize the observed vs predicted plots.
#' @returns A layout (arranged in `nrow` rows and `ncol` columns) of the observed and predicted values of \eqn{log}(relative tumor volume) vs Time for each SampleID (i.e. subject), 
#' with the actual measurements, the regression line for each SampleID, and the marginal, treatment-specific, 
#' regression line for each treatment group.
#' @examples
#' #' data(grwth_data)
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
#'   combination = "Combination",
#'   show_plot = FALSE
#'   )
#' # Obtain the plots
#' plot_ObsvsPred(lmm, nrow = 4, ncol = 8)    
#' 
#' @export
plot_ObsvsPred <- function(model, nrow = 4, ncol = 5){
  TV.df <- model$data
  aug.Pred <- nlme::augPred(model, primary = ~Time, level = 0:1, length.out = 2, minimum = 0)
  plot(aug.Pred, layout = c(ncol, nrow, 1), lty = c(1,2),
       key = list(lines = list(lty = c(1,2), col = c("slateblue", "orange")),
                  text = list(c("Marginal", "Subject-specific")),
                  columns = 2,
                  space="top"), 
       pch = 20, lwd = 1.5, main = "Observed and Predicted Values by Time",
       par.strip.text=list(col="black", cex=1),
       xlab = list("Time", cex = 1.2), ylab = list("log RTV", cex = 1.2))
}
