## Power Evaluation using Simulations 

#' @title Post hoc power calculation based on simulations of the synergy evaluation using LMM.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param nsim Number of simulations to perform.
#' @param method String indicating the method for synergy calculation. Possible methods are "Bliss", "HSA" and "RA",
#' corresponding to Bliss, highest single agent and response additivity, respectively.
#' @param pvalue Threshold for the p-value of synergy calculation to be considered statistically significant.
#' @param ... Additional parameters to be passed to [nlmeU::simulateY].
#' @returns Returns a numeric value of the power for the synergy calculation for the model using the method specified in `method`. 
#' The power is expressed as the proportion of simulations that provides a p-value below the threshold specified in `pvalue`. 
#' @export

PostHocPwr <- function(model,
                       nsim = 1000,
                       method = "Bliss",
                       pvalue = 0.05,
                       ...) {
  # Validate method input
  valid_methods <- c("Bliss", "HSA", "RA")
  if (!method %in% valid_methods) {
    stop("Invalid 'method' provided. Choose from 'Bliss', 'HSA', or 'RA'.")
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
  if (method == "RA") {
    contrast <- "b4 = log(exp(b2) + exp(b3) - exp(b1))"
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

#' @title A priori power calculation for a hypothetical study of synergy evaluation using LMM.
#' @param npg Number of mouse per group.
#' @param time Vector with the times at which the tumor volume measurements have been performed.
#' @param grwrControl Coefficient for Control treatment group tumor growth rate.
#' @param grwrA Coefficient for Drug A treatment group tumor growth rate.
#' @param grwrB Coefficient for Drug B treatment group tumor growth rate.
#' @param grwrComb Coefficient for Combination (Drug A + Drug B) treatment group tumor growth rate.
#' @param sd_ranef Random effects standard deviation for the model.
#' @param sgma Residuals standard deviation for the model.
#' @param sd_eval A vector with values representing the standard deviations of random effects,
#' through which the power for synergy calculation will be evaluated.
#' @param sgma_eval A vector with values representing the standard deviations of the residuals,
#' through which the power for synergy calculation will be evaluated.
#' @param grwrComb_eval A vector with values representing the coefficients for Combination treatment group tumor growth rate,
#' through which the power for synergy calculation will be evaluated.
#' @param method String indicating the method for synergy calculation. Possible methods are "Bliss" and "HSA",
#' corresponding to Bliss and highest single agent, respectively.
#' @param ... Additional parameters to be passed to [nlmeU::Pwr.lme] method.
#' @returns The functions returns several plots:
#' - A plot representing the hypothetical data, with the regression lines for each
#' treatment group according to `grwrControl`, `grwrA`, `grwrB` and `grwrComb` values. The values 
#' assigned to `sd_ranef` and `sgma` are also shown.
#' - A plot showing the values of the power calculation depending on the values assigned to 
#' `sd_eval` and `sgma_eval`,
#' - A plot showing the values of the power calculation depending on the values assigned to
#' `grwrComb_eval`.
#' If saved to a variable, the function saves the exemplary data frame built for the hypothetical study.
#' @import ggplot2
#' @importFrom nlme lme lmeControl pdDiag
#' @importFrom cowplot theme_cowplot plot_grid
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
  p1 <- plot_exmpDt(
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

#' @title A priori power calculation for a hypothetical study of synergy evaluation using LMM
#' depending on the sample size per group
#' @param npg A vector with the sample size (number of animals) per group to calculate the power of 
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
#' @returns The functions returns two plots:
#' - A plot representing the hypothetical data, with the regression lines for each
#' treatment group according to `grwrControl`, `grwrA`, `grwrB` and `grwrComb` values. The values 
#' assigned to `sd_ranef` and `sgma` are also shown.
#' - A plot showing the values of the power calculation depending on the values assigned to 
#' `npg`.
#' 
#' If save to a variable, the fuction returns the data frame with the power for the analysis for each sample size
#' especified in `npg`.
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
  p1 <- plot_exmpDt(
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
    labs(title = paste("Power depending on\nnumber of animals for", method)) + scale_x_continuous(breaks = npg) +
    geom_hline(yintercept = 0.8, lty = "dashed")
  
  print(plot_grid(p1, p2, ncol = 2))
  return(invisible(npg_Pwr))
}

## Power with varying times of follow-up or frequency of measurements

#' @title A priori power calculation for a hypothetical study of synergy evaluation using LMM
#' depending on the time of follow-up or the frequency of measurements
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
#' @export 

PwrTime <- function(npg = 5, time = list(seq(0, 9, 3), seq(0, 21, 3), seq(0, 30, 3)), type = "max", 
                        grwrControl, grwrA, grwrB, grwrComb, sd_ranef, sgma, method = "Bliss", ...){
  
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
      title <- "Power depending on\nfrequency of measurements"
      x.lab <- "Number of measurements"
    }
    Pwr_vector <- c(Pwr_vector, dtF$Power)
  }
  
  p1 <- plot_exmpDt(exmpDt, grwrControl = C, grwrA = A, grwrB = B, grwrComb = AB, sd_ranef = sd_ranef, sgma = sgma)
  
  npg_Pwr <- data.frame(Time = time_vector, Power = Pwr_vector)
  
  p2 <- npg_Pwr %>% ggplot(aes(x = .data$Time, y = .data$Power)) + geom_line() + cowplot::theme_cowplot() + xlab(x.lab) + 
    labs(title = paste(title, method)) + scale_x_continuous(breaks = time_vector) +
    geom_hline(yintercept = 0.8, lty = "dashed")
  
  print(plot_grid(p1,p2, ncol = 2))
  return(invisible(npg_Pwr))
}


