#' @importFrom ggplot2 aes geom_hline geom_line ggplot labs scale_x_continuous xlab 
#' @importFrom rlang .data
NULL

## Power with varying number of subjects 

#' @title _A Priori_ Synergy Power Analysis Based on Sample Size
#' @description
#'  _A priori_ power calculation for a hypothetical two-drugs combination study of synergy evaluation using linear-mixed models
#' depending on the sample size per group.
#' @param npg A vector with the sample size (number of subjects) per group to calculate the power of 
#' the synergy analysis.
#' @param time Vector with the times at which the tumor volume measurements have been performed.
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
#' `PwrSampleSize` allows the user to define an hypothetical drug combination study, customizing several 
#' experimental parameters, such as the sample size, time of measurements, or drug effect,
#' for the power evaluation of synergy for Bliss and HSA reference models. The power calculation is
#' based on F-tests of the fixed effects of the model as previously described (Helms, R. W. (1992), 
#' Verbeke and Molenberghs (2009), Gałecki and Burzykowski (2013)). 
#' 
#' The focus the power analysis with `PwrSampleSize` is on the **sample size per group**. The function allows
#' for the evaluation of how the statistical power changes when the sample size per group varies while the
#' other parameters are kept constant. For other _a priori_ power analysis see also [`APrioriPwr()`] and [`PwrTime()`].
#' 
#' - `time`, `grwrControl`, `grwrA`, `grwrB`, `grwrComb`, `sd_ranef` and `sgma` are parameters referring to
#' the initial exemplary data set and corresponding fitted model. These values can be obtained from a fitted model, using [`lmmModel_estimates()`],
#' or be defined by the user.
#' -  `npg` is a vector indicating the different sample sizes for which the statistical power is going to be evaluated, keeping the 
#' rest of parameters fixed.
#' @returns The functions returns a plot showing the values of the power calculation depending on the values assigned to 
#' `npg`.
#' If `plot_exmpDt = TRUE`, a plot representing the hypothetical data, with the regression lines for each
#' treatment group according to `grwrControl`, `grwrA`, `grwrB` and `grwrComb` values is also plotted. The values 
#' assigned to `sd_ranef` and `sgma` are also shown.
#' 
#' The function also returns the data frame with the power for the analysis for each sample size
#' specified in `npg`.
#' @references
#' - Helms, R. W. (1992). _Intentionally incomplete longitudinal designs: I. Methodology and comparison of some full span designs_. Statistics in Medicine, 11(14–15), 1889–1913. https://doi.org/10.1002/sim.4780111411
#' - Verbeke, G. & Molenberghs, G. (2000). _Linear Mixed Models for Longitudinal Data_. Springer New York. https://doi.org/10.1007/978-1-4419-0300-6
#' - Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' @seealso [PostHocPwr], [APrioriPwr()], [PwrTime()].
#' @examples
#' PwrSampleSize(npg = 1:20)
#' 
#' @export 


PwrSampleSize <- function(npg = c(5, 8, 10),
                          time = c(0, 3, 5, 10),
                          grwrControl = 0.08,
                          grwrA = 0.07,
                          grwrB = 0.06,
                          grwrComb = 0.03,
                          sd_ranef = 0.01,
                          sgma = 0.1,
                          method = "Bliss",
                          vF = NULL,
                          plot_exmpDt = FALSE,
                          ...) {
  
  ## Constructing an exemplary dataset
  
  # Validate method input
  valid_methods <- c("Bliss", "HSA")
  if (!method %in% valid_methods) {
    stop("Invalid 'method' provided. Choose from 'Bliss' or 'HSA'.")
  }
  
  Time <- time # Vector with times for tumor volume measurements
  
  npg_vector <- c()
  Pwr_vector <- c()
  
  for (n in npg) {
    # No of subjects per group
    SampleID <- 1:(4 * n) # Subjects' ids
    Treatment <- gl(4, n, labels = c("Control", "DrugA", "DrugB", "Combination")) # Treatment for each subject
    dts <- data.frame(SampleID, Treatment) # Subject-level data
    
    dtL <- list(Time = Time, SampleID = SampleID)
    dtLong <- expand.grid(dtL) # Long format
    mrgDt <- merge(dtLong, dts, sort = FALSE) # Merged
    
    # Controls growth rate
    C <- grwrControl
    
    # Treatment A growth rate
    A <- grwrA
    
    # Treatment B growth rate
    B <- grwrB
    
    # Combination growth rate
    AB <- grwrComb
    
    exmpDt <- within(mrgDt, {
      m0 <- C * Time
      mA <- A * Time
      mB <- B * Time
      mAB <- AB * Time
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
    
    pd1 <- pdDiag(D, form = ~ 0 + Time)
    
    cntrl <- lmeControl(
      maxIter = 0,
      msMaxIter = 0,
      niterEM = 0,
      returnObject = TRUE,
      opt = "optim"
    )
    
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
    
    if (method == "Bliss") {
      dtF <- Pwr(
        fmA,
        sigma = sgma,
        L = c(
          "Time:TreatmentControl" = 1,
          "Time:TreatmentDrugA" = -1,
          "Time:TreatmentDrugB" = -1,
          "Time:TreatmentCombination" = 1
        ),
        ...
      )
    }
    if (method == "HSA") {
      if (which.min(c(grwrA, grwrB)) == 1) {
        dtF <- Pwr(
          fmA,
          sigma = sgma,
          L = c(
            "Time:TreatmentDrugA" = -1,
            "Time:TreatmentCombination" = 1
          ),
          ...
        )
      } else{
        dtF <- Pwr(
          fmA,
          sigma = sgma,
          L = c(
            "Time:TreatmentDrugB" = -1,
            "Time:TreatmentCombination" = 1
          ),
          ...
        )
      }
    }
    
    npg_vector <- c(npg_vector, n)
    Pwr_vector <- c(Pwr_vector, dtF$Power)
  }
  
  # Ploting exemplary data
  p1 <- .plot_exmpDt(
    exmpDt,
    grwrControl = C,
    grwrA = A,
    grwrB = B,
    grwrComb = AB,
    sd_ranef = sd_ranef,
    sgma = sgma
  )
  
  npg_Pwr <- data.frame(N = npg_vector, Power = Pwr_vector)
  
  p2 <- npg_Pwr %>% ggplot(aes(x = .data$N, y = .data$Power)) + geom_line() + cowplot::theme_cowplot() + xlab("N per group") +
    labs(title = paste("Power depending on\nnumber of subjects per group for", method)) + scale_x_continuous(breaks = npg) +
    geom_hline(yintercept = 0.8, lty = "dashed")
  
  if (plot_exmpDt == TRUE) {
    plot(plot_grid(p1, p2, ncol = 2))
  } else {
    plot(p2)
  }
  
  return(npg_Pwr)
}
