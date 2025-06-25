#' @title Get estimates from a linear mixed model of tumor growth data
#' @description
#' `lmmModel_estimates` allows the user to easily extract some of the interesting model estimates for further use in other functions, 
#' such as for power calculation.
#' 
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param robust If TRUE, sandwich-based robust estimators
#' of the standard error of the regression coefficient estimates by [clubSandwich::conf_int()] are provided. Sandwich-based robust estimators
#' are only available for exponential growth models ('explme').
#' @param type Character string specifying which small-sample adjustment should be used when `robust = True`. Available options are "CR0", "CR1", "CR1p", "CR1S", "CR2", or "CR3". 
#' See "Details" section of [clubSandwich::vcovCR()] for further information.
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

lmmModel_estimates <- function(model,
                               robust = FALSE,
                               type = "CR2") {
  
  UseMethod("lmmModel_estimates")
}

## lmmModel_estimates.explme method
#' This is the method for \code{lmmModel_estimates()} generic function to calculate synergy
#' using the exponential tumor growth model.
#' @rdname lmmModel_estimates
#' @method lmmModel_estimates explme
#' @export
#' 
lmmModel_estimates.explme <- function(model,
                                      robust = FALSE,
                                      type = "CR2"){
  
  # Get table with coefficients of fixed effects and std errors
  
  # Using sandwich robust estimators
  if (robust){
    tTable <- as.data.frame(clubSandwich::conf_int(model, vcov = clubSandwich::vcovCR(model, type = type)))
    tTable <- tTable[,c("beta", "SE")]
  } else {
    # If not, obtain them directly from the model
    tTable <- summary(model)$tTable[,c("Value", "Std.Error")]
  }
  

  # Build table with fixed effects and residuals estimates
  dt <- data.frame(matrix(t(tTable), nrow = 1, byrow = TRUE), # Fixed effects and their std. dev.
                   sqrt(diag(nlme::pdMatrix(model$modelStruct$reStruct[[1]]))), # std. dev. of random effects
                   model$sigma, row.names = NULL) # std.dev of residuals
  
  # Name columns according to treatments
  trt_names <- names(model$coefficients$fixed)
  trt_names <- sub("Time:Treatment", replacement = "", trt_names)
  trt_names <- trt_names[-length(trt_names)]
  
  trt_names <- as.vector(rbind(trt_names, paste0("sd_", trt_names)))
  
  colnames(dt) <- c(trt_names,"Combination", "sd_Combination","sd_ranef", "sd_resid")
  
  return(dt)
}


## lmmModel_estimates.gompertzlme method
#' This is the method for \code{lmmModel_estimates()} generic function to calculate synergy
#' using the Gompertz tumor growth model.
#' @rdname lmmModel_estimates
#' @method lmmModel_estimates gompertzlme
#' @export
#' 
lmmModel_estimates.gompertzlme <- function(model,
                                           robust = FALSE,
                                           type = "CR2"){
  
  if(robust) {
    warning("Sandwich-based robust estimators are only available for exponential growth models.")
  }
  
  # Get table with coefficients of fixed effects and std errors
  
  tTable <- summary(model)$tTable[, c("Value", "Std.Error")]
  
  # Build table with fixed effects and residuals estimates
  dt <- data.frame(matrix(t(tTable), nrow = 1, byrow = TRUE), # Fixed effects and their std. dev.
                   matrix(sqrt(diag(nlme::pdMatrix(model$modelStruct$reStruct[[1]]))), nrow = 1, byrow = TRUE), # std. dev. of random effects
                   model$sigma, row.names = NULL) # std.dev of residuals
  
  # Name columns according to treatments
  trt_names <- names(model$coefficients$fixed)
  trt_names <- sub("Treatment", replacement = "", trt_names)
  trt_names_r0 <- trt_names[grepl("r0", trt_names)]
  trt_names_r0[length(trt_names_r0)] <- "r0.Combination"
  
  trt_names_rho <- trt_names[grepl("rho", trt_names)]
  trt_names_rho[length(trt_names_rho)] <- "rho.Combination"
  
  
  trt_names_r0 <- as.vector(rbind(trt_names_r0, paste0("sd_", trt_names_r0)))
  trt_names_rho <- as.vector(rbind(trt_names_rho, paste0("sd_", trt_names_rho)))
  trt_names <- c(trt_names_r0, trt_names_rho)
  
  colnames(dt) <- c(trt_names,"sd_r0_ranef", "sd_rho_ranef","sd_resid")
  
  return(dt)
}


