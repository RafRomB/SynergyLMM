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
#'   drug_ab = "Combination"
#'   )
#'  PostHocPwr(lmm, nsim = 100) # 100 simulations for shorter computing time
#'  # Using a seed to obtain reproducible results
#'  PostHocPwr(lmm, seed = 123, nsim = 100)
#' 
#' @export

PostHocPwr <- function(model,
                       nsim = 1000,
                       method = "Bliss",
                       pvalue = 0.05,
                       ...) {
  # Validate method input
  valid_methods <- c("Bliss", "HSA")
  if (!method %in% valid_methods) {
    stop("Invalid 'method' provided. Choose from 'Bliss' or 'HSA'.")
  }
  
  if (method == "Bliss") {
    contrast <- "b4 = b2 + b3 - b1"
  }
  if (method == "HSA") {
    fixef_betas <- nlme::fixef(model)[2:3]
    if (which.min(fixef_betas) == 1) {
      contrast <- "b4 = b2"
    } else{
      contrast <- "b4 = b3"
    }
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


## A priori Power Calculations

#' @title _A Priori_ Synergy Power Analysis Based on Variability and Drug Effect 
#' @description
#' _A priori_ power calculation for a hypothetical study of synergy evaluation using linear-mixed models
#' with varying drug combination effect and/or experimental variability.
#' @param npg Number of subjects per group.
#' @param time Vector with the times at which the tumor volume measurements have been performed.
#' @param grwrControl Coefficient for Control treatment group tumor growth rate.
#' @param grwrA Coefficient for Drug A treatment group tumor growth rate.
#' @param grwrB Coefficient for Drug B treatment group tumor growth rate.
#' @param grwrComb Coefficient for Combination (Drug A + Drug B) treatment group tumor growth rate.
#' @param sd_ranef Random effects standard deviation (between-subject variance) for the model.
#' @param sgma Residuals standard deviation (within-subject variance) for the model.
#' @param sd_eval A vector with values representing the standard deviations of random effects,
#' through which the power for synergy calculation will be evaluated.
#' @param sgma_eval A vector with values representing the standard deviations of the residuals,
#' through which the power for synergy calculation will be evaluated.
#' @param grwrComb_eval A vector with values representing the coefficients for Combination treatment group tumor growth rate,
#' through which the power for synergy calculation will be evaluated.
#' @param method String indicating the method for synergy calculation. Possible methods are "Bliss" and "HSA",
#' corresponding to Bliss and highest single agent, respectively.
#' @param ... Additional parameters to be passed to [nlmeU::Pwr.lme] method.
#' @details
#' `APrioriPwr` allows for total customization of an hypothetical drug combination study and allows the user
#' to define several experimental parameters, such as the sample size, time of measurements, or drug effect,
#' for the power evaluation of synergy for Bliss and HSA reference models. The power calculation is
#' based on F-tests of the fixed effects of the model as previously described (Helms, R. W. (1992), 
#' Verbeke and Molenberghs (2009), Gałecki and Burzykowski (2013)). 
#' 
#' The focus the power analysis with `APrioriPwr` is on the **drug combination effect** and the **variability** in the
#' experiment. For other _a priori_ power analysis see also [`PwrSampleSize()`] and [`PwrTime()`].
#' 
#' - `npg`, `time`, `grwrControl`, `grwrA`, `grwrB`, `grwrComb`, `sd_ranef` and `sgma` are parameters referring to
#' the initial exemplary data set and corresponding fitted model. These values can be obtained from a fitted model, using [`lmmModel_estimates()`],
#' or be defined by the user.
#' - `sd_eval`, `sgma_eval`, and `grwrComb_eval` are the different values that will be modified in the initial
#' exemplary data set to fit the corresponding model and calculate the power.
#' 
#' @returns The functions returns several plots:
#' - A plot representing the hypothetical data, with the regression lines for each
#' treatment group according to `grwrControl`, `grwrA`, `grwrB` and `grwrComb` values. The values 
#' assigned to `sd_ranef` and `sgma` are also shown.
#' - A plot showing the values of the power calculation depending on the values assigned to 
#' `sd_eval` and `sgma_eval`,
#' - A plot showing the values of the power calculation depending on the values assigned to
#' `grwrComb_eval`.
#' 
#' If saved to a variable, the function saves the exemplary data frame built for the hypothetical study.
#' @import ggplot2
#' @importFrom nlme lme lmeControl pdDiag
#' @importFrom cowplot theme_cowplot plot_grid
#' @references
#' - Helms, R. W. (1992). _Intentionally incomplete longitudinal designs: I. Methodology and comparison of some full span designs_. Statistics in Medicine, 11(14–15), 1889–1913. https://doi.org/10.1002/sim.4780111411
#' - Verbeke, G. & Molenberghs, G. (2000). _Linear Mixed Models for Longitudinal Data_. Springer New York. https://doi.org/10.1007/978-1-4419-0300-6
#' - Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' @seealso [PostHocPwr],[PwrSampleSize()], [PwrTime()].
#' @examples
#' APrioriPwr(
#' sd_eval = seq(0.01, 0.2, 0.01),
#' sgma_eval = seq(0.01, 0.2, 0.01),
#' grwrComb_eval = seq(0.01, 0.1, 0.005)
#' )
#' 
#' @export


APrioriPwr <- function(npg = 5,
                       time = c(0, 3, 5, 10),
                       grwrControl = 0.08,
                       grwrA = 0.07,
                       grwrB = 0.06,
                       grwrComb = 0.03,
                       sd_ranef = 0.01,
                       sgma = 0.1,
                       sd_eval = NULL,
                       sgma_eval = NULL,
                       grwrComb_eval = NULL,
                       method = "Bliss",
                       ...) {
  if (is.null(sd_eval) & is.null(sgma_eval) & is.null(grwrComb_eval)) {
    stop(
      "One of the following, 'sd_eval' and 'sgma_eval', or 'grwrComb_eval', arguments must be specified"
    )
  }
  if (!is.null(sd_eval) | !is.null(sgma_eval)) {
    stopifnot("Both, 'sd_eval' and 'sgma_eval', must be specified" = c(!is.null(sd_eval), !is.null(sgma_eval)))
  }
  
  # Validate method input
  valid_methods <- c("Bliss", "HSA")
  if (!method %in% valid_methods) {
    stop("Invalid 'method' provided. Choose from 'Bliss' or 'HSA'.")
  }
  
  ## Constructing an exemplary dataset
  
  npg <- npg # No of subjects per group
  Time <- time # Vector with times of tumor volume measurements
  
  subject <- 1:(4 * npg) # Subjects' ids
  Treatment <- gl(4, npg, labels = c("Control", "DrugA", "DrugB", "Combination")) # Treatment for each subject
  dts <- data.frame(subject, Treatment) # Subject-level data
  
  dtL <- list(Time = Time, subject = subject)
  dtLong <- expand.grid(dtL) # Long format
  mrgDt <- merge(dtLong, dts, sort = FALSE) # Merged
  
  # Control growth rate
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
  
  fmA <- lme(
    mA ~ Time:Treatment,
    random = list(subject = pd1),
    data = exmpDt,
    control = cntrl
  )
  
  fmB <- fmA # Save copy of the model to use with different
             # values of grwrComb
  
  # Use of Pwr() function for a priori power calculations
  
  # Power for different values of variance
  
  if (!is.null(sd_eval) & !is.null(sgma_eval)) {
    # Ploting power curve
    
    if (method == "Bliss") {
      print(Pwr(
        fmA,
        sigma = sgma,
        L = c(
          "Time:TreatmentControl" = 1,
          "Time:TreatmentDrugA" = -1,
          "Time:TreatmentDrugB" = -1,
          "Time:TreatmentCombination" = 1
        )
      ), ...)
    }
    if (method == "HSA") {
      if (which.min(c(grwrA, grwrB)) == 1) {
        print(Pwr(
          fmA,
          sigma = sgma,
          L = c(
            "Time:TreatmentDrugA" = -1,
            "Time:TreatmentCombination" = 1
          )
        ), ...)
      } else{
        print(Pwr(
          fmA,
          sigma = sgma,
          L = c(
            "Time:TreatmentDrugB" = -1,
            "Time:TreatmentCombination" = 1
          )
        ), ...)
      }
    }
    
    power.df <- data.frame(matrix(ncol = 3, nrow = 0), row.names = NULL)
    colnames(power.df) <- c("Power", "SD", "sigma")
    
    idx <- 1
    
    for (i in 1:length(sd_eval)) {
      for (j in 1:length(sgma_eval)) {
        D <- log(sd_eval[i])
        pd1 <- pdDiag(D, form = ~ 0 + Time)
        fmA <- lme(
          mA ~ 0 + Time:Treatment,
          random = list(subject = pd1),
          data = exmpDt,
          control = cntrl
        )
        if (method == "Bliss") {
          dtF <- Pwr(
            fmA,
            sigma = sgma_eval[j],
            L = c(
              "Time:TreatmentControl" = 1,
              "Time:TreatmentDrugA" = -1,
              "Time:TreatmentDrugB" =
                -1,
              "Time:TreatmentCombination" = 1
            ),
            ...
          )
        }
        if (method == "HSA") {
          if (which.min(c(grwrA, grwrB)) == 1) {
            dtF <- Pwr(
              fmA,
              sigma = sgma_eval[j],
              L = c(
                "Time:TreatmentDrugA" = -1,
                "Time:TreatmentCombination" = 1
              ),
              ...
            )
          } else{
            dtF <- Pwr(
              fmA,
              sigma = sgma_eval[j],
              L = c(
                "Time:TreatmentDrugB" = -1,
                "Time:TreatmentCombination" = 1
              ),
              ...
            )
          }
        }
        power.df[idx, "SD"] <- sd_eval[i]
        power.df[idx, "sigma"] <- sgma_eval[j]
        power.df[idx, "Power"] <- dtF$Power
        idx <- idx + 1
      }
    }
    
    p2 <- power.df %>% ggplot(aes(
      x = .data$SD,
      y = .data$sigma,
      z = .data$Power
    )) + geom_raster(aes(fill = .data$Power)) +
      scale_fill_continuous(type = "viridis") + cowplot::theme_cowplot() + labs(title = paste("Power for", method, sep = " ")) +
      xlab("SD for random effects") + ylab("Sigma")
    print(plot_grid(p1, p2, ncol = 2))
  }
  
  if (!is.null(grwrComb_eval)) {
    # Ploting power curve
    
    dif <- grwrComb_eval
    dim(dif) <- c(length(dif), 1)
    
    colnames(dif) <- "Time:TreatmentCombination"
    if (method == "Bliss") {
      print(Pwr(
        fmB,
        sigma = sgma,
        L = c(
          "Time:TreatmentControl" = 1,
          "Time:TreatmentDrugA" = -1,
          "Time:TreatmentDrugB" = -1,
          "Time:TreatmentCombination" = 1
        )
      ))
      dtF <- Pwr(
        fmB,
        sigma = sgma,
        L = c(
          "Time:TreatmentControl" = 1,
          "Time:TreatmentDrugA" = -1,
          "Time:TreatmentDrugB" = -1,
          "Time:TreatmentCombination" = 1
        ),
        altB = dif
      )
    }
    if (method == "HSA") {
      if (which.min(c(grwrA, grwrB)) == 1) {
        print(Pwr(
          fmB,
          sigma = sgma,
          L = c(
            "Time:TreatmentDrugA" = -1,
            "Time:TreatmentCombination" = 1
          )
        ))
        dtF <- Pwr(
          fmB,
          sigma = sgma,
          L = c(
            "Time:TreatmentDrugA" = -1,
            "Time:TreatmentCombination" = 1
          ),
          altB = dif
        )
      } else{
        print(Pwr(
          fmB,
          sigma = sgma,
          L = c(
            "Time:TreatmentDrugB" = -1,
            "Time:TreatmentCombination" = 1
          )
        ))
        dtF <- Pwr(
          fmB,
          sigma = sgma,
          L = c(
            "Time:TreatmentDrugB" = -1,
            "Time:TreatmentCombination" = 1
          ),
          altB = dif
        )
      }
    }
    p3 <- dtF %>% ggplot(aes(
      x = .data$`Time:TreatmentCombination`,
      y = .data$Power
    )) + geom_line() + cowplot::theme_cowplot() +
      labs(title = paste(
        "Power across growth rate\nvalues for combination treatment for ",
        method
      )) + xlab("Growth rate (logRTV/Times)") +
      geom_hline(yintercept = 0.8, lty = "dashed")
    print(plot_grid(p1, p3, ncol = 2))
  }
  if (!is.null(sd_eval) &
      !is.null(sgma_eval) & !is.null(grwrComb_eval)) {
    print(plot_grid(p1, p2, p3, ncol = 3))
  }
  return(invisible(exmpDt))
}

## Power with varying number of subjects 

#' @title _A Priori_ Synergy Power Analysis Based on Sample Size
#' @description
#'  _A priori_ power calculation for a hypothetical study of synergy evaluation using linear-mixed models
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
#' @returns The functions returns two plots:
#' - A plot representing the hypothetical data, with the regression lines for each
#' treatment group according to `grwrControl`, `grwrA`, `grwrB` and `grwrComb` values. The values 
#' assigned to `sd_ranef` and `sgma` are also shown.
#' - A plot showing the values of the power calculation depending on the values assigned to 
#' `npg`.
#' 
#' If saved to a variable, the fuction returns the data frame with the power for the analysis for each sample size
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
    subject <- 1:(4 * n) # Subjects' ids
    Treatment <- gl(4, n, labels = c("Control", "DrugA", "DrugB", "Combination")) # Treatment for each subject
    dts <- data.frame(subject, Treatment) # Subject-level data
    
    dtL <- list(Time = Time, subject = subject)
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
    
    fmA <- lme(
      mA ~ Time:Treatment,
      random = list(subject = pd1),
      data = exmpDt,
      control = cntrl
    )
    
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
  
  print(plot_grid(p1, p2, ncol = 2))
  return(invisible(npg_Pwr))
}

## Power with varying times of follow-up or frequency of measurements

#' @title _A Priori_ Synergy Power Analysis Based on Time
#' @description
#'  _A priori_ power calculation for a hypothetical study of synergy evaluation using linear-mixed models
#' depending on depending on the time of follow-up or the frequency of measurements.
#' @param npg Number of mouse per group.
#' @param time A list in which each element is a vector with the times at which the tumor volume measurements have been performed.
#' If `type` is set to "max", each vector in the list should have the measurements taken at the same interval and differ in the final
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
#' @returns The functions returns two plots:
#' - A plot representing the hypothetical data, with the regression lines for each
#' treatment group according to `grwrControl`, `grwrA`, `grwrB` and `grwrComb` values. The values 
#' assigned to `sd_ranef` and `sgma` are also shown.
#' - A plot showing the values of the power calculation depending on the values assigned to 
#' `Time`. If `type` is set to "max", the plot shows how the power varies depending on the maximum time of follow-up. 
#' If `type` is set to "freq", the plot shows how the power varies depending on how frequently the measurements have
#' been performed.
#' 
#' If saved to a variable, the function returns the power for the analysis for each value specified in ` Time`.
#' 
#' #' @references
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
    
    subject <- 1:(4*npg) # Subjects' ids
    Treatment <- gl(4, npg, labels = c("Control", "DrugA", "DrugB", "Combination")) # Treatment for each subject
    dts <- data.frame(subject, Treatment) # Subject-level data
    
    dtL <- list(Time = d, subject = subject)
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
    
    fmA <- lme(mA ~ Time:Treatment, random = list(subject = pd1), data = exmpDt, control = cntrl)
    
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
  
  print(plot_grid(p1,p2, ncol = 2))
  return(invisible(npg_Pwr))
}


