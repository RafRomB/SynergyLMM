#' @importFrom ggplot2 aes annotate coord_cartesian geom_line geom_point ggplot labs scale_color_manual scale_fill_continuous xlab
#' @importFrom rlang .data 
NULL

#' @title Helper function to plot exemplary data for power calculation
#' @description
#' `plot_exmpDt` plots the regression lines of any exemplary data produced for a priori power
#' calculation.
#' @param exmpDt Data frame with exemplary data obtained with [APrioriPwr()].
#' @param grwrControl Value for the label of the coefficient for Control treatment group tumor growth rate.
#' @param grwrA Value for the label of the coefficient for Drug A treatment group tumor growth rate.
#' @param grwrB Value for the label of the coefficient for Drug B treatment group tumor growth rate.
#' @param grwrComb Value for the label of the coefficient for Combination (Drug A + Drug B) treatment group tumor growth rate.
#' @param sd_ranef Value for the label of the random effects standard deviation of the model.
#' @param sgma Value for the label of the residuals standard deviation of the model.
#' @returns A ggplot2 plot (see [ggplot2::ggplot()] for more details) showing the regression lines corresponding 
#' to the fixed effects for each treatment of the exemplary data for power calculations.
#' @keywords internal
#' @noRd
.plot_exmpDt <- function(exmpDt, grwrControl = NULL, grwrA = NULL, grwrB=NULL, grwrComb=NULL, sd_ranef=NULL, sgma=NULL){
  # Ploting exemplary data
  
  selDt <- with(exmpDt,{
    lvls <- levels(Treatment)
    i <- match(lvls, Treatment)
    subj <- subject[i]
    subset(exmpDt, subject %in% subj)
  })
  
  selDt %>% ggplot(aes(x = .data$Time, y = .data$mA)) + geom_line(aes(colour = .data$Treatment), lwd = 2) + 
    labs(title = "Exemplary Data") + ylab("log (RTV)") + xlab("Time since start of treatment") + cowplot::theme_cowplot() +
    annotate(geom = "text", x = max(selDt$Time)+3, y = max(selDt$mA), label = paste("GR Control=",grwrControl), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Time)+3, y = max(selDt$mA)*0.95, label = paste("GR Drug A=",grwrA), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Time)+3, y = max(selDt$mA)*0.9, label = paste("GR Drug B=",grwrB), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Time)+3, y = max(selDt$mA)*0.85, label = paste("GR Combination=", grwrComb), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Time)+3, y = max(selDt$mA)*0.8, label = paste("SD=",sd_ranef), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Time)+3, y = max(selDt$mA)*0.75, label = paste("Sigma=",sgma), hjust = -0.05, size = 4) +
    coord_cartesian(xlim = c(0, max(selDt$Time)+5), clip = "off") +
    scale_color_manual(values = c("#3c3c3b", "#d50c52", "#00a49c", "#601580"))
}

## A priori Power Calculations

#' @title _A Priori_ Synergy Power Analysis Based on Variability and Drug Effect 
#' @description
#' _A priori_ power calculation for a hypothetical two-drugs combination study of synergy evaluation using linear-mixed models
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
#' `sd_eval` and `sgma_eval`. The power result corresponding to the values assigned to `sd_ranef` and
#' `sgma` is shown with a red dot.
#' - A plot showing the values of the power calculation depending on the values assigned to
#' `grwrComb_eval`. The vertical dashed line indicates the value of `grwrComb`. The horizontal
#' line indicates the power of 0.80.
#' 
#' The statistical power for the fitted model for the initial data set according to the values given by
#' `npg`, `time`, `grwrControl`, `grwrA`, `grwrB`, `grwrComb`, `sd_ranef` and `sgma` is also shown in the console.
#' 
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
      pwr.result <- Pwr(
        fmA,
        sigma = sgma,
        L = c(
          "Time:TreatmentControl" = 1,
          "Time:TreatmentDrugA" = -1,
          "Time:TreatmentDrugB" = -1,
          "Time:TreatmentCombination" = 1
        )
        ,
        ...
      )
    }
    if (method == "HSA") {
      if (which.min(c(grwrA, grwrB)) == 1) {
        pwr.result <- Pwr(
          fmA,
          sigma = sgma,
          L = c(
            "Time:TreatmentDrugA" = -1,
            "Time:TreatmentCombination" = 1
          ),
          ...
        )
      } else {
        pwr.result <- Pwr(
          fmA,
          sigma = sgma,
          L = c(
            "Time:TreatmentDrugB" = -1,
            "Time:TreatmentCombination" = 1
          )
          ,
          ...
        )
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
    )) + ggplot2::geom_raster(aes(fill = .data$Power)) +
      scale_fill_continuous(type = "viridis") + cowplot::theme_cowplot() + labs(title = paste("Power for", method, sep = " ")) +
      xlab("SD for random effects") + ylab("SD for residuals") + geom_point(x = sd_ranef, y = sgma, shape = 18, size = 5, color = "firebrick3")
    print(plot_grid(p1, p2, ncol = 2))
  }
  
  if (!is.null(grwrComb_eval)) {
    # Ploting power curve
    
    dif <- grwrComb_eval
    dim(dif) <- c(length(dif), 1)
    
    colnames(dif) <- "Time:TreatmentCombination"
    if (method == "Bliss") {
      pwr.result <- Pwr(
        fmB,
        sigma = sgma,
        L = c(
          "Time:TreatmentControl" = 1,
          "Time:TreatmentDrugA" = -1,
          "Time:TreatmentDrugB" = -1,
          "Time:TreatmentCombination" = 1
        )
        ,
        ...
      )
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
        pwr.result <- Pwr(
          fmB,
          sigma = sgma,
          L = c(
            "Time:TreatmentDrugA" = -1,
            "Time:TreatmentCombination" = 1
          )
          ,
          ...
        )
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
        pwr.result <- Pwr(
          fmB,
          sigma = sgma,
          L = c(
            "Time:TreatmentDrugB" = -1,
            "Time:TreatmentCombination" = 1
          )
          ,
          ...
        )
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
      ggplot2::geom_hline(yintercept = 0.8, lty = "dashed") + ggplot2::geom_vline(xintercept = grwrComb, lty=3) 
    print(plot_grid(p1, p3, ncol = 2))
  }
  if (!is.null(sd_eval) &
      !is.null(sgma_eval) & !is.null(grwrComb_eval)) {
    print(plot_grid(p1, p2, p3, ncol = 3))
  }
  return(pwr.result)
}