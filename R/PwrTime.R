#' @importFrom ggplot2 aes geom_hline geom_line ggplot labs scale_x_continuous xlab
#' @importFrom rlang .data
NULL

## Power with varying times of follow-up or frequency of measurements

#' @title _A Priori_ Synergy Power Analysis Based on Time
#' @description
#'  _A priori_ power calculation for a hypothetical two-drugs combination study of synergy
#' depending on the time of follow-up or the frequency of measurements.
#' @param npg Number of mouse per group.
#' @param time A list in which each element is a vector with the times at which the tumor volume measurements have been performed.
#' If `type` is set to "max", each vector in the list should represent measurements taken at the same interval and differ in the final
#' time of follow-up. If `type` is set to "freq", each vector in the list should have the same final time of follow-up and
#' differ in the intervals at which the measurements have been taken. 
#' @param type String indicating whether to calculate the power depending on the time of follow-up ("max"), or the frequency
#' of measurements ("freq").
#' @param grwrControl Coefficient for Control treatment group tumor growth rate.
#' @param grwrA Coefficient for Drug A treatment group tumor growth rate.
#' @param grwrB Coefficient for Drug B treatment group tumor growth rate.
#' @param grwrComb Coefficient for Combination (Drug A + Drug B) treatment group tumor growth rate.
#' @param sd_ranef Random effects standard deviation for the model.
#' @param sgma Residuals standard deviation for the model.
#' @param method String indicating the method for synergy calculation. Possible methods are "Bliss" and "HSA",
#' corresponding to Bliss and highest single agent, respectively.
#' @param vF An optional [nlme::varFunc] object or one-sided formula describing the within-group heteroscedasticity
#' structure. If given as a formula, it is used as the argument to [nlme::varFixed], corresponding to fixed variance weights. 
#' See the documentation on [nlme::varClasses] for a description of the available [nlme::varFunc] classes. Defaults to NULL, corresponding to 
#' homoscedastic within-group errors.
#' @param plot_exmpDt Logical indicating if a plot representing the hypothetical data should also be returned.
#' @param ... Additional parameters to be passed to [nlmeU::Pwr.lme] method.
#' @details
#' `PwrTime` allows the user to define an hypothetical drug combination study, customizing several 
#' experimental parameters, such as the sample size, time of measurements, or drug effect,
#' for the power evaluation of synergy for Bliss and HSA reference models. The power calculation is
#' based on F-tests of the fixed effects of the model as previously described (Helms, R. W. (1992), 
#' Verbeke and Molenberghs (2009), Gałecki and Burzykowski (2013)). 
#' 
#' The focus the power analysis with `PwrTime` is on the **time** at which the measurements are done. The function allows
#' for the evaluation of how the statistical power changes when the time of follow-up varies while the frequency
#' of measurements is keep constant. It also allows to how the statistical power changes when the time of follow-up is
#' kept constant, but the frequency of measurements varies.
#' 
#' For other _a priori_ power analysis see also [`APrioriPwr()`] and [`PwrSampleSize()`].
#' 
#' - `npg`, `grwrControl`, `grwrA`, `grwrB`, `grwrComb`, `sd_ranef` and `sgma` are parameters referring to
#' the initial exemplary data set and corresponding fitted model. These values can be obtained from a fitted model, using [`lmmModel_estimates()`],
#' or be defined by the user.
#' -  `time` is a list in which each element is a vector with the times at which the tumor volume measurements have been performed, and for
#' which the statistical power is going to be evaluated, keeping the rest of parameters fixed.
#' 
#' @returns The functions returns a plot showing the values of the power calculation depending on the values assigned to 
#' `Time`. If `type` is set to "max", the plot shows how the power varies depending on the maximum time of follow-up. 
#' If `type` is set to "freq", the plot shows how the power varies depending on how frequently the measurements have
#' been performed.
#' If `plot_exmpDt = TRUE`, a plot representing the hypothetical data, with the regression lines for each
#' treatment group according to `grwrControl`, `grwrA`, `grwrB` and `grwrComb` values is also plotted. The values 
#' assigned to `sd_ranef` and `sgma` are also shown.
#' 
#' The function also returns the data frame with the power for the analysis for each value specified in ` Time`.
#' 
#' @references
#' - Helms, R. W. (1992). _Intentionally incomplete longitudinal designs: I. Methodology and comparison of some full span designs_. Statistics in Medicine, 11(14–15), 1889–1913. https://doi.org/10.1002/sim.4780111411
#' - Verbeke, G. & Molenberghs, G. (2000). _Linear Mixed Models for Longitudinal Data_. Springer New York. https://doi.org/10.1007/978-1-4419-0300-6
#' - Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' @seealso [PostHocPwr], [APrioriPwr()], [PwrSampleSize()].
#' @examples
#' # Power analysis maintaining the frequency of measurements 
#' # and varying the time of follow-up ('type = "max"')
#' PwrTime(time = list(seq(0, 9, 3), 
#'                     seq(0, 12, 3), 
#'                     seq(0, 15, 3), 
#'                     seq(0, 21, 3), 
#'                     seq(0, 30, 3)), 
#'                     type = "max")
#' 
#' # Power analysis maintaining the time of follow-up 
#' # and varying the frequency of measurements ('type = "freq"')
#' PwrTime(time = list(seq(0, 10, 1), 
#'                     seq(0, 10, 2), 
#'                     seq(0, 10, 5), 
#'                     seq(0, 10, 10)), 
#'                     type = "freq")

#' @export 

PwrTime <- function(npg = 5,
                    time = list(seq(0, 9, 3), seq(0, 21, 3), seq(0, 30, 3)),
                    type = "max",
                    grwrControl = 0.08,
                    grwrA = 0.07,
                    grwrB = 0.06,
                    grwrComb = 0.03,
                    sd_ranef = 0.01,
                    sgma = 0.1 ,
                    method = "Bliss",
                    vF = NULL,
                    plot_exmpDt = FALSE,
                    ...) {
  
  
  # Validate method input
  valid_methods <- c("Bliss", "HSA")
  if (!method %in% valid_methods) {
    stop("Invalid 'method' provided. Choose from 'Bliss' or 'HSA'.")
  }
  
  # Validate type input
  valid_type <- c("max", "freq")
  if (!type %in% valid_type) {
    stop(paste(type, ": Invalid 'type' provided. Choose from 'max' or 'freq'.", sep = ""))
  }
  
  # Validate time input
  Time <- time # List with the times for the meassurements
  
  if(type == "max"){
    m <- c()
    for (i in Time) {
      m <- c(m, max(i))
    }
    m <- unique(m)
    if(length(m) < length(Time)){
      warning(
        "Your list 'time' has several vectors with the same maximum time of follow-up.\nConsider using 'type = freq' to evaluate the effect of frequency of masurements on power calculation."
      )
    }
  }
  
  ## Constructing an exemplary dataset
  
  time_vector <- c()
  Pwr_vector <- c()
  
  for(d in Time){ # Vector with times
    
    SampleID <- 1:(4*npg) # Subjects' ids
    Treatment <- gl(4, npg, labels = c("Control", "DrugA", "DrugB", "Combination")) # Treatment for each subject
    dts <- data.frame(SampleID, Treatment) # Subject-level data
    
    dtL <- list(Time = d, SampleID = SampleID)
    dtLong <- expand.grid(dtL) # Long format
    mrgDt <- merge(dtLong, dts, sort = FALSE) # Merged
    
    # Define betas according to doubling time for each condition (doubling time = log(2)/beta)
    
    # Controls growth rate
    C <- grwrControl
    
    # Treatment A growth rate
    A <- grwrA
    
    # Treatment B growth rate
    B <- grwrB
    
    # Combination growth rate
    AB <- grwrComb
    
    exmpDt <- within(mrgDt, {
      m0 <- C*Time
      mA <- A*Time
      mB <- B*Time
      mAB <- AB*Time
    })
    
    exmpDt$mA[exmpDt$Treatment == "Control"] <- exmpDt$m0[exmpDt$Treatment == "Control"]
    exmpDt$mA[exmpDt$Treatment == "DrugA"] <- exmpDt$mA[exmpDt$Treatment == "DrugA"]
    exmpDt$mA[exmpDt$Treatment == "DrugB"] <- exmpDt$mB[exmpDt$Treatment == "DrugB"]
    exmpDt$mA[exmpDt$Treatment == "Combination"] <- exmpDt$mAB[exmpDt$Treatment == "Combination"]
    
    exmpDt$mAB <- NULL
    exmpDt$mB <- NULL
    
    # Build lme object
    
    ## Objects of class pdMat
    
    sgma <- sgma
    D <- log(sd_ranef)
    
    pd1 <- pdDiag(D, form = ~0+Time)
    
    cntrl <- lmeControl(maxIter = 0, msMaxIter = 0, niterEM = 0, returnObject = TRUE, opt = "optim")
    
    fmA <- do.call(nlme::lme, c(
      list(
        fixed = mA ~  0 + Time:Treatment,
        random = list(SampleID = pd1),
        data = exmpDt,
        control = cntrl,
        weights = vF
      )
    ))
    
    # Use of Pwr() function for a priori power calculations
    
    # Ploting power curve
    
    if(method == "Bliss"){
      dtF <- Pwr(fmA, sigma = sgma, L = c("Time:TreatmentControl" = 1,"Time:TreatmentDrugA" = -1,
                                          "Time:TreatmentDrugB" = -1,"Time:TreatmentCombination" = 1), ...)
    }
    if(method == "HSA"){
      if(which.min(c(grwrA, grwrB)) == 1){
        dtF <- Pwr(fmA, sigma = sgma, L = c("Time:TreatmentDrugA" = -1,"Time:TreatmentCombination" = 1), ...)
      } else{
        dtF <- Pwr(fmA, sigma = sgma, L = c("Time:TreatmentDrugB" = -1,"Time:TreatmentCombination" = 1), ...)
      }
    }
    
    if(type == "max"){
      time_vector <- c(time_vector, max(d))
      title <- "Power depending on\ntime of follow-up for"
      x.lab <- "Maximum time of follow-up"
    }
    if(type == "freq"){
      time_vector <- c(time_vector, length(d))
      title <- "Power depending on\nfrequency of measurements for"
      x.lab <- "Number of measurements"
    }
    Pwr_vector <- c(Pwr_vector, dtF$Power)
  }
  
  p1 <- .plot_exmpDt(exmpDt, grwrControl = C, grwrA = A, grwrB = B, grwrComb = AB, sd_ranef = sd_ranef, sgma = sgma)
  
  npg_Pwr <- data.frame(Time = time_vector, Power = Pwr_vector)
  
  p2 <- npg_Pwr %>% ggplot(aes(x = .data$Time, y = .data$Power)) + geom_line() + cowplot::theme_cowplot() + xlab(x.lab) + 
    labs(title = paste(title, method)) + scale_x_continuous(breaks = time_vector) +
    geom_hline(yintercept = 0.8, lty = "dashed")
  
  if (plot_exmpDt == TRUE) {
    plot(plot_grid(p1,p2, ncol = 2))
  } else {
    plot(p2)
  }
  
  return(npg_Pwr)
}


