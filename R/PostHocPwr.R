#' @importFrom rlang .data
NULL
## Power Evaluation using Simulations 

#' @title Post hoc power calculation based on simulations of the synergy evaluation using LMM.
#' @description
#' `PostHocPwr` allows for the _post hoc_ power analysis of the synergy hypothesis testing for Bliss and HSA refence
#' models for a given tumor growth data fitted model.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param nsim Number of simulations to perform.
#' @param method String indicating the method for synergy calculation. Possible methods are "Bliss" and "HSA",
#' corresponding to Bliss and highest single agent, respectively.
#' @param pvalue Threshold for the p-value of synergy calculation to be considered statistically significant.
#' @param time Time point for which to calculate the statistical power. If not specified, the last time point is used
#' by default.
#' @param ... Additional parameters to be passed to [nlmeU::simulateY]:
#' @details
#' The _post hoc_ power calculation relies on simulation of the dependent variable, using [nlmeU::simulateY].
#' 1. For a given fitted model of the tumor growth data, `nsim` simulations of the dependent variable (\eqn{\log (RTV)})
#' are done, based on the marginal distribution implied by the fitted model.
#' 2. The model is then fitted to the new values of the dependant variable.
#' 3. For each simulation, the new estimates from each model are then used for the synergy hypothesis testing as
#' explained in [lmmSynergy], and the p-values stored.
#' 4. The power is returned as the proportion of simulations resulting in a significant synergy hypothesis testing
#' (p-value < `pvalue`).
#' 
#' When `time` is specified, `PostHocPwr` refits the model using the data from the `time_start` time point defined in [lmmModel()] until `time`, and report the
#' statistical power for that model. If `time` is not specified, the model fitted using all data points is used for statistical power calculation.
#' 
#' @returns Returns a numeric value of the power for the synergy calculation for the model using the method specified in `method`. 
#' The power is expressed as the proportion of simulations that provides a p-value below the threshold specified in `pvalue`.
#' @references 
#' Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' @seealso [APrioriPwr()], [PwrSampleSize()], [PwrTime()].
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
#'   combination = "Combination"
#'   )
#'  PostHocPwr(lmm, nsim = 10) # 10 simulations for shorter computing time
#'  # Using a seed to obtain reproducible results
#'  PostHocPwr(lmm, seed = 123, nsim = 10)
#'  # Calculating the power for an specific day
#'  PostHocPwr(lmm, nsim = 10, time = 6)
#' 
#' @export

PostHocPwr <- function(model,
                       nsim = 1000,
                       method = "Bliss",
                       pvalue = 0.05,
                       time = NA,
                       ...) {
  
  tryCatch({
    if (!"explme" %in% class(model)) {
      stop(
        "Post hoc power calculation is only available for exponential growth models ('explme')."
      )
    }
    
  }, error = function(e) {
    stop("Error in model class: ", e$message)
  })
  
  # Validate method input
  valid_methods <- c("Bliss", "HSA")
  if (!method %in% valid_methods) {
    stop("Invalid 'method' provided. Choose from 'Bliss' or 'HSA'.")
  }
  
  fixef_betas <- nlme::fixef(model)
  
  if (method == "Bliss") {
    if (length(fixef_betas) == 5) {
      contrast <- "b5 = b2 + b3 + b4 - 2*b1"
    } else {
      contrast <- "b4 = b2 + b3 - b1"
    }
  }
  
  if (method == "HSA") {
    if (length(fixef_betas) == 5) {
      if (which.min(fixef_betas[2:4]) == 1) {
        contrast <- "b5 = b2"
      } else if (which.min(fixef_betas[2:4]) == 2){
        contrast <- "b5 = b3"
      } else {
        contrast <- "b5 = b4"
      }
    } else {
      if (which.min(fixef_betas[2:3]) == 1) {
        contrast <- "b4 = b2"
      } else{
        contrast <- "b4 = b3"
      } 
    }
  }
  
  if(!is.na(time)){
    data <- model$data %>% dplyr::filter(.data$Time <= time)
    model <- update(model, data = data)
  }
  
  simA <- nlmeU::simulateY(model, nsim = nsim, ...) # Simulation
  dt <- model$data # working copy
  simfmA <- apply(simA, 2, function(y) {
    dt$logRTV <- y
    auxFit <- update(model, data = dt)
    marginaleffects::hypotheses(auxFit, hypothesis = contrast)
  })
  FstateE <-
    sapply(simfmA, function(x)
      x$p.value)
  
  powerE <- sum(FstateE < pvalue) / nsim
  return(powerE)
}

