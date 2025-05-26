#' @title Get estimates from a linear mixed model of tumor growth data
#' @description
#' `lmmModel_estimates` allows the user to easily extract some of the interesting model estimates for further use in other functions, 
#' such as for power calculation.
#' 
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @details
#' The model estimates provided by `lmmModel_estimates` include:
#' - Fixed effect coefficients: \eqn{\hat{\beta}_{Control}}, \eqn{\hat{\beta}_A}, \eqn{\hat{\beta}_B}, ( \eqn{\hat{\beta}_{C}}), \eqn{\hat{\beta}_{Combination}}, 
#' which represent the estimated specific growth rates for the Control, Drug A, Drug B, (Drug C, if present), and Combination groups, respectively.
#' - The standard deviation (sd) corresponding to each of the fixed effect coefficients. 
#' - Standard deviation of the random effects (between-subject variance). Column `sd_ranef`.
#' - Standard deviation of the residuals (within-subject variance). Column `sd_resid`.
#' 
#' @returns A data frame with the estimated values for the coefficients of the tumor growth for each treatment, their standard error,
#' the standard deviation of the random effects, and the standard deviation of the residuals of the model.
#' These values can be useful for the power analysis of the model using [`APrioriPwr()`].
#' @examples
#' data("grwth_data")
#' # Fit example model
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
#' # Get the estimates
#' lmmModel_estimates(lmm)
#' @export

lmmModel_estimates <- function(model){
  
  # Get table with coefficients of fixed effects and std errors
  tTable <- summary(model)$tTable[,c("Value", "Std.Error")]
  # Build table with fixed effects and residuals estimates
  dt <- data.frame(matrix(t(tTable), nrow = 1, byrow = TRUE), sqrt(model$modelStruct$reStruct[[1]][1]), model$sigma)
  
  # Name columns according to treatments
  trt_names <- names(model$coefficients$fixed)
  trt_names <- sub("Time:Treatment", replacement = "", trt_names)
  trt_names <- trt_names[-length(trt_names)]
  
  trt_names <- as.vector(rbind(trt_names, paste0("sd_", trt_names)))
  
  colnames(dt) <- c(trt_names,"Combination", "sd_Combination","sd_ranef", "sd_resid")
  
  rownames(dt) <- "estimate"
  return(dt)
}
