#' @importFrom ggplot2 aes element_text ggplot labs stat_qq stat_qq_line theme xlab ylab
#' @import graphics
#' @importFrom stats fitted qqnorm resid residuals
NULL


#' @title Plots for random effects diagnostics
#' @description
#' Visualization of random effects diagnostics for a fitted linear mixed model of tumor growth data.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @returns A list with different plots for evaluating the normality and homoscedasticity of the random effects, including:
#' - A normal Q-Q plot of the random effects of the model.
#' - A normal Q-Q plot of the residuals by sample.
#' - Boxplots of the "raw" residuals (observed - fitted) by sample.
#' - Scatter plots of the normalized residuals (standardized residuals pre-multiplied by the inverse square-root factor of the estimated error correlation matrix, see [nlme::residuals.lme])
#' vs fitted values by sample. Observations with absolute standardized (normalized) residuals greater than the \eqn{1-0.05/2} quantile of the standard normal distribution 
#' are identified in the plots labelled with the time point corresponding to the observation.
#' 
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
#' plot_ranefDiagnostics(lmm)
#' # Access to specific plots
#' plot_ranefDiagnostics(lmm)$plots[[1]]
#' plot_ranefDiagnostics(lmm)$plots[[2]]
#' @export
#' 

plot_ranefDiagnostics <- function(model) {
  
  ranef_mod <- nlme::ranef(model)
  ranef_names <- colnames(ranef_mod)
  p1 <- list()
  for (i in 1:ncol(ranef_mod)){
    p <- ggplot(ranef_mod, aes(sample = ranef_mod[,ranef_names[i]])) + stat_qq(col = "gray20") + 
      stat_qq_line() +
      labs(title = paste("Normal Q-Q Plot of", ranef_names[i],"Random Effects")) + 
      xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + cowplot::theme_cowplot() +
      theme(plot.title = element_text(size = 10, hjust = 0.5), axis.title = element_text(size = 12))
    p1[[i]] <- p
  }
  
  
  p2 <- qqnorm(model, ~resid(., type = "normalized")|SampleID, pch=20, cex = 0.5, col = "gray20",
               main = list("Normal Q-Q Plot of Normalized Residuals by Sample", cex = 0.8), 
               par.strip.text=list(col="black", cex=0.8), xlab = "Normalized Residuals", abline = c(0,1))
  p3 <- plot(model, SampleID ~ resid(., type = "response"), abline = 0, main = list("Raw Residuals by Sample", cex = 0.8),
             xlab = "Residuals")
  p4 <- plot(model, residuals(., type = "normalized") ~ fitted(.)|SampleID, id = 0.05, adj = -0.03, pch = 20, col = "slateblue4", cex=0.75,
             main = list("Normalized Residuals vs Fitted Values by Sample", cex =0.8),par.strip.text=list(col="black", cex=0.8), idLabels = ~Time,
             abline = 0, ylab = "Normalized residuals")
  p <- append(p1, list(p2,p3,p4))
  
  # Arranged plots
  p5 <- cowplot::plot_grid(plotlist = p, ncol = 2)
  return(list(plots = p, grid = p5))
}

