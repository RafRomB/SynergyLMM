#' @import graphics
#' @importFrom stats fitted qqnorm resid
NULL



#' @title Plots for residuals diagnostics
#' @description
#' Visualization of residuals diagnostics for a fitted linear mixed model of tumor growth data.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @returns A list with different plots for evaluating the normality and homoscedasticity of the normalized residuals
#' (standardized residuals pre-multiplied by the inverse square-root factor of the estimated error correlation matrix, see [nlme::residuals.lme]), including:
#' - A normal Q-Q plot of the normalized residuals of the model.
#' - A normal Q-Q plot of the normalized residuals of the model by Time.
#' - A normal Q-Q plot of the normalized residuals of the model by Treatment.
#' - A dotplot of normalized residuals vs fitted values.
#' - A dotplot of the normalized residuals by Time and Treatment.
#' @examples
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
#'   combination = "Combination",
#'   show_plot = FALSE
#'   )
#' # Generate plots 
#' plot_residDiagnostics(lmm)
#' # Access to specific plots
#' plot_residDiagnostics(lmm)[[1]]
#' plot_residDiagnostics(lmm)[[2]]
#' @export
plot_residDiagnostics <- function(model){
  # Individual Plots
  p1 <- qqnorm(model, ~resid(., type = "normalized"),pch = 20, main = "Q-Q Plot of Normalized Residuals",
               xlab = "Normalized residuals", abline = c(0,1))
  p2 <- qqnorm(model, ~resid(., type = "normalized")|Time,pch = 20, main = "Q-Q Plot of Normalized Residuals by Time",
               par.strip.text=list(col="black", cex=1),
               xlab = "Normalized residuals", abline = c(0,1))
  p3 <- qqnorm(model, ~resid(., type = "normalized")|Treatment, pch = 20, main = "Q-Q Plot of Normalized Residuals by Treatment",
               par.strip.text=list(col="black", cex=1),
               xlab = "Normalized residuals", abline = c(0,1))
  p4 <- plot(model, resid(., type = "normalized") ~ fitted(.), main = "Normalized Residuals vs Fitted Values", 
             pch = 20, ylab = "Normalized residuals", abline = 0)
  p5 <- plot(model, resid(., type = "normalized") ~Time|Treatment, id = 0.05, pch=20,
             adj = -0.03, cex = 1, main = "Normalized Residuals per Time and Treatment", 
             par.strip.text=list(col="black", cex=1), ylab = "Normalized residuals", abline = 0)
  # Arranged plots
  p6 <- cowplot::plot_grid(p1,p2,p3,p4,p5, nrow = 3, ncol = 2, align = "hv")
  return(list(p1,p2,p3,p4,p5,p6))
}