

## 4. Power Evaluation using Simulations ----

Pwr_simulate <- function(model, nsim=1000, method = "Bliss", pvalue = 0.05){
  
  if(method == "Bliss"){
    contrast <- "b4 = b2 + b3 - b1"
  }
  if(method == "HSA" ){
    fixef_betas <- fixef(model)[2:3]
    if(which.min(fixef_betas) == 1){
      contrast <- "b4 = b2"
    } else{
      contrast <- "b4 = b3"
    }
  }
  if(method == "RA"){
    contrast <- "b4 = log(exp(b2) + exp(b3) - exp(b1))"
  }
  
  simA <- simulateY(model, nsim = nsim) # Simulation
  dt <- model$data # working copy
  simfmA <- apply(simA,
                  2,
                  function(y){
                    dt$logRTV <- y
                    auxFit <- update(model, data = dt)
                    hypotheses(auxFit, hypothesis = contrast)
                  })
  FstateE <- 
    sapply(simfmA, function(x) x$p.value)
  
  powerE <- sum(FstateE < pvalue)/nsim
  return(powerE)
}

## 5. A priori Power Calculations ----

#' @import ggplot2
#' @importFrom nlme lme lmeControl pdDiag
#' @importFrom cowplot theme_cowplot plot_grid
#' @export

Pwr_lmm <- function(npg = 5, Day = c(0,3,5,10), grwrControl, grwrA, grwrB, grwrComb, sd_ranef, sgma, sd_eval = NULL, 
                    sgma_eval = NULL, grwrComb_eval = NULL, method = "Bliss"){
  
  if(is.null(sd_eval) & is.null(sgma_eval) & is.null(grwrComb_eval)){
    stop("One of the following, 'sd_eval' and 'sgma_eval', or 'grwrComb_eval', arguments must be specified")
  }
  if(!is.null(sd_eval) | !is.null(sgma_eval)){
    stopifnot("Both, 'sd_eval' and 'sgma_eval', must be specified" = c(!is.null(sd_eval), !is.null(sgma_eval)))
  }
  
  ## Constructing an exemplary dataset
  
  npg <- npg # No of subjects per group
  subject <- 1:(4*npg) # Subjects' ids
  Treatment <- gl(4, npg, labels = c("Control", "TreatA", "TreatB", "Combination")) # Treatment for each subject
  dts <- data.frame(subject, Treatment) # Subject-level data
  
  dtL <- list(Day = Day, subject = subject)
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
    m0 <- C*Day
    mA <- A*Day
    mB <- B*Day
    mAB <- AB*Day
  })
  
  exmpDt$mA[exmpDt$Treatment == "Control"] <- exmpDt$m0[exmpDt$Treatment == "Control"]
  exmpDt$mA[exmpDt$Treatment == "TreatA"] <- exmpDt$mA[exmpDt$Treatment == "TreatA"]
  exmpDt$mA[exmpDt$Treatment == "TreatB"] <- exmpDt$mB[exmpDt$Treatment == "TreatB"]
  exmpDt$mA[exmpDt$Treatment == "Combination"] <- exmpDt$mAB[exmpDt$Treatment == "Combination"]
  
  exmpDt$mAB <- NULL
  exmpDt$mB <- NULL
  
  # Ploting exemplary data
  selDt <- with(exmpDt,{
    lvls <- levels(Treatment)
    i <- match(lvls, Treatment)
    subj <- subject[i]
    subset(exmpDt, subject %in% subj)
  })
  
  p1 <- selDt %>% ggplot(aes(x = Day, y = mA)) + geom_line(aes(colour = Treatment), lwd = 2) + 
    labs(title = "Exemplary Data") + ylab("logRTV") + xlab("Days") + theme_cowplot() +
    annotate(geom = "text", x = max(selDt$Day)+3, y = max(selDt$mA), label = paste("GR Control=",C), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Day)+3, y = max(selDt$mA)*0.95, label = paste("GR Drug A=",A), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Day)+3, y = max(selDt$mA)*0.9, label = paste("GR Drug B=",B), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Day)+3, y = max(selDt$mA)*0.85, label = paste("GR Combination=",AB), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Day)+3, y = max(selDt$mA)*0.8, label = paste("SD=",sd_ranef), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Day)+3, y = max(selDt$mA)*0.75, label = paste("Sigma=",sgma), hjust = -0.05, size = 4) +
    coord_cartesian(xlim = c(0, max(selDt$Day)+5), clip = "off")
  
  # Build lme object
  
  ## Objects of class pdMat
  
  sgma <- sgma
  D <- log(sd_ranef)
  
  pd1 <- pdDiag(D, form = ~0+Day)
  
  cntrl <- lmeControl(maxIter = 0, msMaxIter = 0, niterEM = 0, returnObject = TRUE, opt = "optim")
  
  fmA <- lme(mA ~ Day:Treatment, random = list(subject = pd1), data = exmpDt, control = cntrl)
  
  # Use of Pwr() function for a priori power calculations
  
  # Power for different values of variance
  
  if(!is.null(sd_eval) & !is.null(sgma_eval)){
    
    # Ploting power curve
    
    if(method == "Bliss"){
      print(Pwr(fmA, sigma = sgma, L = c("Day:TreatmentControl" = 1,"Day:TreatmentTreatA" = -1,"Day:TreatmentTreatB" = -1,"Day:TreatmentCombination" = 1)))
    }
    if(method == "HSA"){
      if(which.min(c(grwrA, grwrB)) == 1){
        print(Pwr(fmA, sigma = sgma, L = c("Day:TreatmentTreatA" = -1,"Day:TreatmentCombination" = 1)))
      } else{
        print(Pwr(fmA, sigma = sgma, L = c("Day:TreatmentTreatB" = -1,"Day:TreatmentCombination" = 1)))
      }
    }
    
    power.df <- data.frame(matrix(ncol = 3, nrow = 0), row.names = NULL)
    colnames(power.df) <- c("Power", "SD", "sigma")
    
    idx <- 1
    
    for (i in 1:length(sd_eval)) {
      for(j in 1:length(sgma_eval)){
        D <- log(sd_eval[i])
        pd1 <- pdDiag(D, form = ~0+Day)
        fmA <- lme(mA ~ 0+Day:Treatment, random = list(subject = pd1), data = exmpDt, control = cntrl)
        if(method == "Bliss"){
          dtF <- Pwr(fmA,sigma = sgma_eval[j], L = c("Day:TreatmentControl"=1, "Day:TreatmentTreatA"=-1,
                                                     "Day:TreatmentTreatB"=-1, "Day:TreatmentCombination"=1))
        }
        if(method == "HSA"){
          if(which.min(c(grwrA, grwrB)) == 1){
            dtF <- Pwr(fmA, sigma = sgma_eval[j], L = c("Day:TreatmentTreatA" = -1,"Day:TreatmentCombination" = 1))
          } else{
            dtF <- Pwr(fmA, sigma = sgma_eval[j], L = c("Day:TreatmentTreatB" = -1,"Day:TreatmentCombination" = 1))
          }
        }
        power.df[idx, "SD"] <- sd_eval[i]
        power.df[idx, "sigma"] <- sgma_eval[j]
        power.df[idx, "Power"] <- dtF$Power
        idx <- idx + 1
      }
    }
    
    p2 <- power.df %>% ggplot(aes(x = SD, y = sigma, z = Power)) + geom_raster(aes(fill = Power)) + 
      scale_fill_continuous(type = "viridis") + theme_cowplot() + labs(title = paste("Power for", method, sep = " ")) +
      xlab("SD for random effects") + ylab("Sigma")
    print(plot_grid(p1,p2, ncol = 2))
  }
  
  if(!is.null(grwrComb_eval)){
    
    # Ploting power curve
    
    dif <- grwrComb_eval
    dim(dif) <- c(length(dif), 1)
    
    colnames(dif) <- "Day:TreatmentCombination"
    if(method == "Bliss"){
      print(Pwr(fmA, sigma = sgma, L = c("Day:TreatmentControl" = 1,"Day:TreatmentTreatA" = -1,"Day:TreatmentTreatB" = -1,"Day:TreatmentCombination" = 1)))
      dtF <- Pwr(fmA, sigma = sgma, L = c("Day:TreatmentControl" = 1,"Day:TreatmentTreatA" = -1,"Day:TreatmentTreatB" = -1,"Day:TreatmentCombination" = 1), altB = dif)
    }
    if(method == "HSA"){
      if(which.min(c(grwrA, grwrB)) == 1){
        print(Pwr(fmA, sigma = sgma, L = c("Day:TreatmentTreatA" = -1,"Day:TreatmentCombination" = 1)))
        dtF <- Pwr(fmA, sigma = sgma, L = c("Day:TreatmentTreatA" = -1,"Day:TreatmentCombination" = 1), altB = dif)
      } else{
        print(Pwr(fmA, sigma = sgma, L = c("Day:TreatmentTreatB" = -1,"Day:TreatmentCombination" = 1)))
        dtF <- Pwr(fmA, sigma = sgma, L = c("Day:TreatmentTreatB" = -1,"Day:TreatmentCombination" = 1), altB = dif)
      }
    }
    p3 <- dtF %>% ggplot(aes(x = `Day:TreatmentCombination`, y = Power)) + geom_line() + theme_cowplot() +
      labs(title = paste("Power across growth rate\nvalues for combination treatment for ", method)) + xlab("Growth rate (logRTV/Days)") +
      geom_hline(yintercept = 0.8, lty = "dashed")
    print(plot_grid(p1,p3, ncol = 2))
  }
  if(!is.null(sd_eval) & !is.null(sgma_eval) & !is.null(grwrComb_eval)){
    print(plot_grid(p1, p2, p3, ncol = 3))
  }
}

