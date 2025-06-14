# Residual Diagnostics ----

#' @title Diagnostics of residuals of the linear mixed model
#' 
#' @description
#' `residDiagnostics` provides several plots as well as statistical test for the examination
#' of the normality and homoscedasticity of the residuals of the input model.
#' 
#' @param model  An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param pvalue Threshold for the p-value of outlier observations based on their normalized residuals.
#' @param norm_test String indicating the function for testing the normality of the normalized residuals. A collection of functions from
#'  \code{fBasics::\link[fBasics:normalTest]{normalTest}} is available. We recommend using one of "shapiroTest", "dagoTest",
#'  or "adTest" for performing Shapiro - Wilk, 
#' D'Agostino, or Anderson - Darling normality test, respectively.
#' @param verbose Logical indicating if the normality and homoscedasticity tests results, and the list of potential
#' outlier observations should be printed to the console.
#' 
#' @details
#' One of the assumption of the model fit by [`lmmModel()`] is that the residuals are normally distributed.
#' For the evaluation of this assumption, `residDiagnostics` provides Q-Q plots of the normalized residuals
#' (standardized residuals pre-multiplied by the inverse square-root factor of the estimated error correlation matrix, see [nlme::residuals.lme]),
#' together with statistical assessment of their 
#' normality using Shapiro-Wilk, D'Agostini and Anderson-Darling normality tests. Additionally, Q-Q plots of the normalized residuals by time point and 
#' treatment group are provided to be able to detect time points or treatment groups which could be notably different from the others and be 
#' affecting the adequacy of the model. 
#' 
#' Scatter plots of the normalized residuals versus fitted values and normalized residuals 
#' per time and per treatment are also provided to give information about variability of the residuals and possible outlier observations. These plots are accompanied
#' by Levene and Fligner-Killend homogeneity of variance test results.
#' 
#' Observations with absolute standardized (normalized) residuals greater than the \eqn{1-0.05/2} quantile of the standard normal distribution 
#' are identified and reported as potential outlier observations.
#' 
#' @returns A list with different elements for the diagnostics of the residuals are produced:
#' - `plots`: Different plots for evaluating the normality and homocedasticity of the residuals.
#' - `outliers`: Data frame with the identified outliers based on the Pearson residuals and the value of `pval`. The column `resid.p` contains the
#' value of the Pearson residuals for each observation.
#' - `Normality`: Results from the test of the normality of the normalized residuals of the model.
#' - `Levene.test`: List with the Levene homoscedasticity test results of the normalized residuals by Time and Treatment.
#' - `Fligner.test`: List with the Fligner-Killeen homoscedasticity test results of the normalized residuals by Time and Treatment.
#' 
#' @references
#' - Pinheiro JC, Bates DM (2000). _Mixed-Effects Models in S and S-PLUS_. Springer, New York. \doi{doi:10.1007/b98882}.
#' - Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' 
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
#'   
#' # Residuals diagnostics
#' resid_diag <- residDiagnostics(model = lmm, pvalue = 0.05)
#' 
#' # Access outliers data frame
#' resid_diag$Outliers
#' 
#' # Access individual plots
#' resid_diag$Plots[1]
#' resid_diag$Plots[2]
#' 
#' # Access results of normality tests
#' resid_diag$Normality
#' 
#' # Access to homoscedasticity test results
#' 
#' resid_diag$Levene.test
#' 
#' resid_diag$Fligner.test
#'
#' @export
residDiagnostics <- function(model, 
                             pvalue=0.05,
                             norm_test = "shapiroTest",
                             verbose = TRUE) {
  # Plots
  resid_plot <- plot_residDiagnostics(model)
  plot(resid_plot[[6]])
  
  # Normality test
  
  norm_res <- resid(model, type = "normalized")
  
  # Normality test
  
  norm_test <- getFromNamespace(norm_test, "fBasics")
  Normality <- norm_test(norm_res)
  
  # Homocedasticity test
  
  # Normalized Residuals by Time
  
  model$data$normalized_resid <- norm_res
  
  levene <- list()
  fligner <- list()
  
  levene$Time <- car::leveneTest(normalized_resid ~ as.factor(Time), data = model$data)
  
  fligner$Time <- fligner.test(normalized_resid ~ as.factor(Time), data = model$data)
  
  levene$Treatment <- car::leveneTest(normalized_resid ~ Treatment, data = model$data)
  
  fligner$Treatment <- fligner.test(normalized_resid ~ Treatment, data = model$data)
  
  # List of outliers
  outliers.idx <- within(model$data, {
    resid.std <- resid(model, type = "normalized") # Standardized resids.
    idx <- abs(resid.std) > -qnorm(pvalue / 2) # Indicator vector
  })
  
  outliers <- subset(outliers.idx, idx) # Data with outliers
  outliers <- outliers[,-c(9:10)]
  
  if (verbose) {
    writeLines("\nNormalized Residuals Normality Test")
    print(Normality)
    
    writeLines("\nNormalized Residuals Levene Homoscedasticity Test by Time")
    print(levene$Time)
    
    writeLines("\nNormalized Residuals Fligner-Killeen Homoscedasticity Test by Time")
    print(fligner$Time)
    
    writeLines("\nNormalized Residuals Levene Homoscedasticity Test by Treatment")
    print(levene$Treatment)
    
    writeLines("\nNormalized Residuals Fligner-Killeen Homoscedasticity Test by Treatment")
    print(fligner$Treatment)
    
    writeLines("\nOutlier observations")
    print(outliers)
  }
  
  
  return(invisible(
    list(
      Plots = resid_plot,
      Outliers = outliers,
      Normality = Normality, 
      Levene.test = levene,
      Fligner.test = fligner
    )
  ))
}