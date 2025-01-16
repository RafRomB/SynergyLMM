#' @title Get estimates from a linear mixed model of tumor growth data
#' @description
#' `lmmModel_estimates` allows the user to easily extract some of the interesting model estimates for further use in other functions, 
#' such as for power calculation.
#' 
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @details
#' The model estimates provided by `lmmModel_estimates` include:
#' - Fixed effect coefficients: \eqn{\hat{\beta}_C}, \eqn{\hat{\beta}_A}, \eqn{\hat{\beta}_B}, \eqn{\hat{\beta}_{AB}}, 
#' which represent the estimated specific growth rates for the Control, Drug A, Drug B and Combination groups, respectively.
#' These are shown in columns `control`, `drug_a`, `drug_b`, and `combination`, respectively.
#' - Standard error of the random effects (between-subject variance). Column `sd_ranef`.
#' - Standard error of the residuals (within-subject variance). Column `sd_resid`.
#' 
#' @returns A data frame with the estimated values for the coefficients of the tumor growth for each treatment,
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
  dt <- data.frame(t(model$coefficients$fixed), sqrt(model$modelStruct$reStruct[[1]][1]), model$sigma)
  trt_names <- names(model$coefficients$fixed)
  trt_names <- sub("Time:Treatment", replacement = "", trt_names)
  if (ncol(dt) == 7) {
    colnames(dt) <- c(trt_names[1:4],"Combination", "sd_ranef", "sd_resid")
  } else {
    colnames(dt) <- c(trt_names[1:3], "Combination", "sd_ranef", "sd_resid")
  }
  rownames(dt) <- "estimate"
  return(dt)
}
