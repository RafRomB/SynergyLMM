# Random Effects Diagnostics ----
#' @title Diagnostics of random effects of the linear mixed model
#' 
#' @description
#' `ranefDiagnostics` provides several plots as well as statistical test for the examination
#' of the normality of the random effects of the input model.
#' 
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
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
#' - `Normality`: List with the results from 3 different test of the normality of the random effects: Shapiro - Wilk normality test, 
#' D'Agostino normality test and Anderson - Darling normality test.
#' - `Levene.test`: results from Levene homocedasticity test of the normalized residuals by SampleID (i.e., by subject).
#' - `Fligner.test`: results from Fligner homocedasticity test of the normalized residuals by SampleID (i.e., by subject).
#' 
#' @references
#' - Pinheiro JC, Bates DM (2000). _Mixed-Effects Models in S and S-PLUS_. Springer, New York. doi:10.1007/b98882 <https://doi.org/10.1007/b98882>.
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
#' ranef_diag$Normality$Shapiro.test
#' 
#' # Access to homoscedasticity tests of residuals by subject
#' 
#' ranef_diag$Levene.test
#' 
#' ranef_diag$Fligner.test
#' 
#' @export

ranefDiagnostics <- function(model){
  # Plots
  ranef_plot <- plot_ranefDiagnostics(model)
  print(ranef_plot[[5]])
  
  # Normality test
  print(
    ranef_shapiro <- fBasics::shapiroTest(nlme::ranef(model)$Time, 
                                          title = "Shapiro - Wilk Normality Test of random effects")
  )
  
  if (length(nlme::ranef(model)$Time) < 20) {
    ranef_DAgostino <- warning("Sample size must be at least 20 for D'Agostino Normality Test")
  } else {
    print(
      ranef_DAgostino <- fBasics::dagoTest(nlme::ranef(model)$Time, 
                                           title = "D'Agostino Normality Test of random effects")
    )
  }
  print(
    ranef_ad <- fBasics::adTest(nlme::ranef(model)$Time, 
                                title = "Anderson - Darling Normality Test of random effects")
  )
  
  Normality <- list(
    Shapiro.test = ranef_shapiro,
    DAgostino.test = ranef_DAgostino,
    Anderson.Darling.test = ranef_ad
  )
  
  # Homocedasticity test
  
  # Normalized Residuals
  
  norm_res <- residuals(model, type = "normalized")
  norm_res <- data.frame(
    normalized_resid = norm_res,
    SampleID = names(norm_res),
    stringsAsFactors = T
  )
  colnames(norm_res) <- c("normalized_resid", "SampleID")
  
  writeLines("\nNormalized Residuals Levene Homoscedasticity Test by Sample")
  print(levene <- car::leveneTest(normalized_resid ~ SampleID, data = norm_res))
  
  writeLines("\nNormalized Residuals Fligner-Killeen Homoscedasticity Test by Sample")
  print(fligner <- fligner.test(normalized_resid ~ SampleID, data = norm_res))
  
  return(invisible(
    list(
      Plots = ranef_plot,
      Normality = Normality,
      Levene.test = levene,
      Fligner.test = fligner
    )
  ))
}



