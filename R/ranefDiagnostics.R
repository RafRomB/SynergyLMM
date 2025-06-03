#' @importFrom utils getFromNamespace
#' @importFrom fBasics normalTest
# Random Effects Diagnostics ----
#' @title Diagnostics of random effects of the linear mixed model
#' 
#' @description
#' `ranefDiagnostics` provides several plots as well as statistical test for the examination
#' of the normality of the random effects of the input model.
#' 
#' @param model An object of class "lme" or "nlme" representing the mixed-effects model fitted by [`lmmModel()`].
#' @param norm_test String indicating the function for testing the normality of the random effects. A collection of functions from
#'  \code{fBasics::\link[fBasics:normalTest]{normalTest}} is available. We recommend using one of "shapiroTest", "dagoTest",
#'  or "adTest" for performing Shapiro - Wilk, 
#' D'Agostino, or Anderson - Darling normality test, respectively.
#' @param verbose Logical indicating if the normality and homoscedasticity tests results should be printed to the console.
#' @details
#' One of the assumptions of the model obtained with [`lmmModel()`] (as in any other linear mixed model) is that
#' the random effects are normally distributed:
#' 
#' \deqn{b_i = N(0,\psi)}
#' 
#' For the evaluation of this assumption, `ranefDiagnostics` provides Q-Q plots of random effects, 
#' together with statistical assessment of their normality using Shapiro-Wilk, D'Agostini and Anderson-Darling normality tests. 
#' Additionally, Q-Q plots of the normalized residuals
#' (standardized residuals pre-multiplied by the inverse square-root factor of the estimated error correlation matrix, see [nlme::residuals.lme])
#' by sample are provided to allow for the identification of subjects 
#' which could be notably different from the others and be affecting the adequacy of the model. 
#' Additionally, boxplots of the "raw" residuals (observed - fitted) by sample and scatter plots of the normalized residuals versus fitted values by sample 
#' are provided to give information about variability of the residuals by subject and possible outlier observations. Observations with absolute standardized (normalized) 
#' residuals greater than the \eqn{1-0.05/2} quantile of the standard normal distribution 
#' are identified in the scatter plots labelled with the time point corresponding to the observation.
#' 
#' @returns A list with different elements for the diagnostics of the random effects are produced:
#' - `plots`: Different plots for evaluating the normality and homoscedasticity of the random effects are produced.
#' - `Normality`: Results from the test of the normality of the random effects.
#' - `Levene.test`: results from Levene homoscedasticity test ([car::leveneTest()]) of the normalized residuals by SampleID (i.e., by subject).
#' - `Fligner.test`: results from Fligner homoscedasticity test ([stats::fligner.test()]) of the normalized residuals by SampleID (i.e., by subject).
#' 
#' @references
#' - Pinheiro JC, Bates DM (2000). _Mixed-Effects Models in S and S-PLUS_. Springer, New York. \doi{doi:10.1007/b98882}.
#' - Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' @seealso [plot_ranefDiagnostics()]
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
#' # Run random effects diagnostics
#' ranef_diag <- ranefDiagnostics(lmm)
#' 
#' #Access to individual plots
#' 
#' ranef_diag$Plots[1]
#' ranef_diag$Plots[2]
#' 
#' # Access to normality tests
#' 
#' ranef_diag$Normality
#' 
#' # Access to homoscedasticity tests of residuals by subject
#' 
#' ranef_diag$Levene.test
#' 
#' ranef_diag$Fligner.test
#' 
#' @export

ranefDiagnostics <- function(model,
                             norm_test = "shapiroTest",
                             verbose = TRUE) {
  # Plots
  ranef_plot <- plot_ranefDiagnostics(model)
  plot(ranef_plot$grid)
  ranef_plot <- ranef_plot$plots
  
  # Normality test
  
  norm_test <- getFromNamespace(norm_test, "fBasics")
  
  ranef_mod <- nlme::ranef(model)
  ranef_names <- colnames(ranef_mod)
  Normality <- list()
  
  for (i in 1:ncol(ranef_mod)) {
    Normality[[i]] <- norm_test(ranef_mod[,ranef_names[i]],
                                 description = paste("Normality Test of",
                                 ranef_names[i],"random effects"))
  }
  
  names(Normality) <- ranef_names
  
  # Homoscedasticity test
  
  # Normalized Residuals
  
  norm_res <- residuals(model, type = "normalized")
  norm_res <- data.frame(
    normalized_resid = norm_res,
    SampleID = names(norm_res),
    stringsAsFactors = T
  )
  colnames(norm_res) <- c("normalized_resid", "SampleID")
  
  levene <- car::leveneTest(normalized_resid ~ SampleID, data = norm_res)
  
  fligner <- fligner.test(normalized_resid ~ SampleID, data = norm_res)
  
  if (verbose) {
    writeLines("\nNormality Test of Random Effects")
    print(Normality)
    
    writeLines("\nNormalized Residuals Levene Homoscedasticity Test by Sample")
    print(levene)
    
    writeLines("\nNormalized Residuals Fligner-Killeen Homoscedasticity Test by Sample")
    print(fligner)
  }
  
  return(invisible(
    list(
      Plots = ranef_plot,
      Normality = Normality,
      Levene.test = levene,
      Fligner.test = fligner
    )
  ))
}

