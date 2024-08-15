sigmaTolme <- function(object, value){ 
  ### Use this function only with Pwr(), because it  corrupts lme.object
  funNm <- "sigmaTolme"
  
  .functionName <- "simgmaTolme"                                # Function name
  .traceR <- if (is.null(options()$traceR)) function(...){} else options()$.traceR 
  
  .traceR(1, lbl = "-> sigmaTolme STARTS")
  sigma0 <- object$sigma 
  val <- value * value
  sc  <- sqrt(val)/sigma0  
  object$sigma <- sqrt(val)
  #resids <- object$residuals
  #resids <- resids * sc  
  #std <- attr(resids,"std")*sc
  ## attr(object$residuals,"std") <- std
  
  attr(object$fixDF, "varFixFact") <- 
    sc*attr(object$fixDF, "varFixFact") # Rescaled for anova
  .traceR(5, vcov(object)*sc*sc, funNm)
  
  object$varFix  <- object$varFix*sc*sc  # vcov rescaled  
  .traceR(1, lbl = "sigmaTolme ENDS <-")
  object
}


## -> Pwr function
#' Calculates power based on a model fit
#'
#' This function is generic; method functions can be written to handle specific classes of objects.
#' @export
#' @param object an object containing the results returned by a model fitting function (e.g., \code{lme}).
#' @param \dots some methods for this generic function may require additional arguments.
#' @return numeric scalar value.
#' @author Andrzej Galecki and Tomasz Burzykowski
#' @seealso \code{\link{Pwr.lme}}
#' @examples
#'  \dontrun{
#'   Pwr (fm1)
#' }
Pwr <-  function(object, ...) UseMethod("Pwr")


## -> Pwr.lme function
#' Performs power calculations
#'
#' This is method for \code{Pwr()} generic function. It works fine for an
#' example given in the book. It may require additional testing, especially for
#' post-hoc power analysis
#'
#' @method Pwr lme
#' @param object an object containing \code{lme} fit, which provides information needed for power calculations
#' @param \dots some additional arguments may be required.
#' @param type an optional character string specifying the type of sum of squares to be used in F-tests 
#'  needed for power calculations. Syntax is the same as for \code{anova.lme()} in \code{nlme} package.
#' @param Terms an optional integer or character vector specifying which terms
#'  in the model should be jointly tested to be zero using a Wald F-test. See
#'  \code{anova.lme} in \code{nlme} package for details.
#' @param L an optional numeric vector or array specifying linear combinations
#'  of the coefficients in the model that should be tested to be zero. See
#'  \code{anova.lme} in \code{nlme} package for details.
#' @param verbose an optional logical value. See \code{anova.lme} in nlme package for details.
#' @param sigma numeric scalar value.
#' @param ddf numeric scalar value. Argument can be used to redefine default number of denominator degrees of freedom
#' @param alpha numeric scalar value. By default 0.05.
#' @param altB matrix/vector containing alternative values for beta parameters
#' @param tol numeric scalar value.
#' @return a data frame inheriting from class Pwr.lme
#' @exportS3Method Pwr lme
#' @author Andrzej Galecki and Tomasz Burzykowski
#' @seealso \code{\link{anova.lme}}
#' 
Pwr.lme <- function (object, ...,
                     type = c("sequential", "marginal"), 
                     Terms, L, verbose = FALSE, sigma, ddf= numeric(0), alpha=0.05,
                     altB = NULL, tol = 1e-10, dt1efs = NULL) 
{
  ## Arguments: 
  ##  1. object: one object only
  #   2. adjustSigma set to FALSE. (Try to explore adjustSigma argument if missing(sigma) 
  #   3. Alternative altB name
  
  .functionName <- "Pwr.lme"                                # Function name
  .traceR <- if (is.null(options()$traceR)) function(...){} else options()$.traceR    
  
  .traceR(1, lbl = "-> Pwr.lme starts")
  if (!inherits(object, "lme")) {
    stop("Object must inherit from class \"lme\" ")
  }
  
  Lmiss <- missing(L)
  Tmiss <- missing(Terms)    
  Lx <- !Tmiss || !Lmiss     # Contrasts matrix L can be created 
  fixefs <- object$coefficients$fixed
  fixefNms <- names(fixefs)
  if (!Lx && !missing(sigma)) stop("L or Terms arguments need to be specified with non-missing sigma argument") 
  ### optx <- if (object$sigma < tol) FALSE else TRUE
  
  .traceR(10) 
  if (!Tmiss && Lmiss){          # IF 1  (Check this part)
    .traceR(101, lbl = "IF1", store = FALSE)
    ##  Based on Terms argument L matrix is created (with colnames)
    ##  Colnames assigned
    
    cLnms <- fixefNms
    assign <- attr(object$fixDF, "assign")
    nTerms <- length(assign)
    nX <- length(unlist(assign))
    L <- diag(nX)[unlist(assign[Terms]), , drop = FALSE]
    colnames(L) <- cLnms 
    cLnms <- NULL
  }
  
  if (!missing(sigma)) object <- sigmaTolme(object,  value=sigma) # nlmeU::: 
  
  x <- anova(object, adjustSigma=FALSE,
             test=TRUE,type=type, L=L, verbose=verbose)  # ANOVA results stored in x 
  .traceR(15)  
  
  
  rt <- attr(x,"rt")   ### Check whether rt is needed
  
  ndf <- x[["numDF"]]
  
  ddf2 <- if (length(ddf)) ddf else x[["denDF"]] #  ddf2 is needed
  rankL <- ndf
  
  
  .traceR(20, lbl = "Before if Lx")
  
  if (Lx){ 
    .traceR(201, lbl = "if Lx", store = FALSE) 
    ## Terms or L present: cLnms created for later use
    dimL <-  if (Lmiss) NULL else dim(L)             # Extract clNms from L argument
    condt <- !Lmiss && is.null(dimL)
    ## cLnms over-written
    cLnms <- if (condt) names(L) else colnames(L)  
    
  }
  
  .traceR(202, lbl = "ifLxBefore")
  if (!is.null(altB) && Lx){ # START
    .traceR(210, lbl = "IF 4", store = FALSE)
    
    if (!Lmiss && is.null(dimL))  {
      .traceR(211, lbl = "IF 3", store = FALSE)
      #print("if here 2")
      dim(L) <- c(1,length(L))  # L is matrix
      colnames(L) <- cLnms
    }
    .traceR(215, lbl = "altB")
    ### Go through altB
    altBdt <- as.data.frame(altB)
    altBnrow <- nrow(altBdt)
    altBnms  <- names(altBdt)
    fixefx <- object$coefficients$fixed
    dt1 <- data.frame(matrix(rep(fixefx, altBnrow), nrow= altBnrow, byrow = TRUE))
    names(dt1) <- fixefNms
    #print(names(dt1))
    dt1[, altBnms] <- altBdt
    #print(names(dt1))
    
    ##  Trimming dimensions
    vcovb <- object$varFix
    if (length(cLnms) < length(fixefNms)){
      vcovb <- vcovb[cLnms,cLnms]
      dt1efs  <- as.matrix(dt1[cLnms])
    }
    if(length(cLnms) == length(fixefNms)){
      dt1efs  <- as.matrix(dt1)
    }
    .traceR(216, lbl = "altB")
    
    
    Fstat <- function(bx) {
      b1 <- L %*% bx
      res0 <- t(b1) %*% solve( L %*% vcovb %*% t(L)) %*% b1
      # Return Lcontrasts i.e. b1 vector together with res0
      res0/rankL
    }
    
    FstatAll <- apply(dt1efs, 1, Fstat)
    .traceR(217, lbl = "FstatAll")
    #dtAll1  <- data.frame(dt1,  numDF = ndf, denDF = ddf2, 
    #          Fvalue = FstatAll)
    
    dtAll2 <- within(dt1, {
      numDF <- ndf
      denDF <- ddf2
      Fvalue <- FstatAll
      Fcrit <- qf(1-alpha, numDF, denDF)
      nc   <- Fvalue * numDF
      Power <- 1- pf(Fcrit, numDF,denDF,nc)
    })
    .traceR(218, lbl = "EXITdtAll2")
    return(dtAll2)
  }  # END
  .traceR(21, lbl = "ifLx After")
  Fcrit <- qf(1-alpha, ndf, ddf2)
  
  ncx <- x[["F-value"]]*ndf  # Rescaling not needed
  
  Power <- 1- pf(Fcrit, ndf, ddf2, ncx) 
  
  F0val <- x[["F-value"]] 
  .traceR(25, lbl=  "F0val")
  
  ret <- data.frame(x$numDF, ddf2, F0val, ncx, Power) #Fcrit omitted
  
  
  varNames <- c("numDF", "denDF","F-value", "nc", "Power")  # Fcrit omitted
  vcovb <- object$varFix 
  
  .traceR(30, lbl=  "if attr(x,L) before")
  if (!is.null(axL <- attr(x, "L"))) {  # L mtx specified
    .traceR(301, lbl= "IF5", store = FALSE)      
    if (!Tmiss) lab <- paste("Power calculations for effect(s):", 
                             paste(Terms, collapse = ", "), "\n",
                             " represented by linear combination: \n")
    
    if (!Lmiss)  lab <- "Power calculations for a linear combination: \n"
    dimL <- dim(L)
    if (is.null(dimL)) names(L) <- cLnms
    attr(ret,"L") <- L
    names(ret) <- varNames
  } else {   #  L not-specified. Effects tested one by one.
    .traceR(302, lbl= "ELSE5", store = FALSE)       
    # lab <- paste("Power calculations for effect(s):", Terms,"\n")  
    dimnames(ret) <- list(rownames(x),varNames)  ## ???
  }
  .traceR(30, lbl=  "if attr(x,L) after")  
  
  ### Modify lab.
  ### if alpha= ne 0.05 then paste(lab, alpha=)
  ### if ddf is NULL then warning ddf2 incorrect
  if (!is.null(attr(x, "label"))) {  
    attr(ret,"label") <- lab 
  } else {
    attr(ret,"label") <- "Power calculations: \n" 
  }
  attr(ret, "coefficients") <- object$coefficients$fixed
  attr(ret,"varFixed")   <- vcovb
  
  attr(ret,"alpha") <- alpha
  attr(ret,"rt") <- rt
  
  class(ret)  <- c("Pwr.lme","data.frame")
  .traceR(1, lbl="Pwr.lme ENDS <-")
  ret
}

#' @exportS3Method print Pwr
print.Pwr <- function (x, verbose = attr(x, "verbose"), ...) 
{
  if ((rt <- attr(x, "rt")) == 1) {
    if (!is.null(lab <- attr(x, "label"))) {
      cat(lab)
      if (!is.null(L <- attr(x, "L"))) {
        print(zapsmall(L))
      }
    }
    ## cat("??print \n")
    pval <- format(round(x[, "Power"], 4))
    pval[as.double(pval) == 0] <- "<.0001"
    x[, "F-value"] <- format(zapsmall(x[, "F-value"]))
    x[, "nc"] <- format(zapsmall(x[, "nc"]))
    ## x[, "Fcrit"] <- format(zapsmall(round(x[, "Fcrit"], 3)))
    x[, "Power"] <- pval
    print(as.data.frame(x))
  }
}



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

