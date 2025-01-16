#' @import nlme
#' @import lattice
NULL

# Influential Diagnostics ----

## Leave-one-out 

#' @title Helper function for creating the object `lmeUall` containing fitted models
#' with leave-one-out data
#' @param cx Subject to remove from the data to build the model
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`],
#' and fitted using maximum likelihood.
#' @returns A list with the leave-one-out model fits
#' @keywords internal
#' @noRd

.lmeU <- function(cx, model){
  SampleID <- NULL
  dfU <- subset(model$data, SampleID != cx) ## LOO data
  update(model, data = dfU)
}

# logLik1 function for varStruct with LOO data

#' @title Modified [nlmeU::logLik1] helper function to deal with varIdent structure for leave-one-out data
#' @param modfit An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`],
#' with a variance structure defined by [nlme::varIdent], fitted using maximum likelihood and using leave-one-out data.
#' @param dt1 A data frame with data for one subject, for whom the log-likelihood function is to be evaluated
#' @param dtInit An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`],
#' with a variance structure defined by [nlme::varIdent], fitted using maximum likelihood, with the complete data.
#' @param var_name Name of the variable for the weights of the model in the case that a variance structure has been specified using [nlme::varIdent()].
#' @returns Numeric scalar value representing contribution of a given subject to the overall 
#' log-likelihood returned by `logLik()` function applied to the "lme" object defined by `modfit` argument.
#' @keywords internal
#' @noRd

.logLik1.varIdent_loo <- function(modfit, dt1, dtInit, var_name){
  var_name <- as.character(unique(dt1[,var_name]))
  m <- modfit$modelStruct # Model structure
  sigma <- modfit$sigma # sigma
  D <- as.matrix(m$reStruct[[1]]) 
  D <- D  * sigma^2 # Matrix D 
  
  vecR <- sigma/(nlme::varWeights(dtInit$modelStruct$varStruct)) # AugDiagonal of R_i
  if(length(var_name) == 1){
    vecR <- vecR[names(vecR) == var_name]
  } else{
    vecR <- vecR[var_name]
  }  
  
  # Ensure right dimension of vecR
  vecR <- vecR[seq_len(nrow(dt1))]
  
  vecR2 <- vecR^2
  R <- diag(vecR2, nrow=length(vecR)) # R_i matrix     
  Z <- model.matrix(m$reStruc, data=dt1) # Z_i matrix
  V <- Z %*% D %*% t(Z) + R # V_i matrix
  predY <- predict(modfit, dt1, level=0) # Predict fixed
  
  dvName <- as.character(formula(modfit)[[2]])
  
  r <- dt1[[dvName]] - predY # Residuals
  n <- nrow(dt1) # No. of obs for subject
  lLik <- n*log(2*pi) + log(det(V)) + 
    t(r) %*% solve(V) %*% r
  return(-0.5 * as.numeric(lLik))
}

# Calculation of the likelihood displacement
#' @title Helper function for the calculation of the likelihood displacement for every subject
#' @param cx Subject that has been remove from the data to build the model with leave-one-out data.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`],
#' fitted using maximum likelihood.
#' @param lmeUall A list with the leave-one-out model fits obtained by [`.lmeU()`].
#' @param var_name Name of the variable for the weights of the model in the case that a variance structure has been specified using [nlme::varIdent()].
#' @returns Numeric value indicating the displacement in the log-likelihood due to removal of subject `cx`.
#' @keywords internal
#' @noRd
#' 
.lLik <- function(cx, model, lmeUall, var_name){
  SampleID <- NULL
  lmeU <- lmeUall[[cx]] # LOO fit extracted
  lLikU <- logLik(lmeU, REML = FALSE) # LOO log-likelihood
  df.s <- subset(model$data, SampleID == cx) # Data for subject cx...
  if(is.null(model$modelStruct$varStruct)){
    lLik.s <- nlmeU::logLik1(lmeU, df.s) # ... and log-likelihood  
  } else{
    if(is.null(var_name)){
      stop("`var_name` cannot be NULL if a variance estructure has been specified in the model")
    }
    lLik.s <- .logLik1.varIdent_loo(lmeU, df.s, model, var_name) # ... and log-likelihood
  }
  return(lLikU + lLik.s) # "Displaced" log-likelihood
} 

#' @title Likelihood displacements for the model
#' @description
#' `logLikSubjectDisplacements` allows the user to evaluate the log-likelihood displacement for each subject, 
#' indicating the influence of every subject to the model.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param disp_thrh Numeric value indicating the threshold of log-likelihood displacement. If not specified, the threshold is set to the 90% percentile of the log-likelihood
#' displacement values.
#' @param label_angle Numeric value indicating the angle for the label of subjects with a log-likelihood displacement greater than `disp_thrh`.
#' @param var_name Name of the variable for the weights of the model in the case that a variance structure has been specified using [nlme::varIdent()].
#' (See examples in [`lmmModel()`]).
#' @param ... Extra arguments, if any, for [lattice::panel.xyplot].
#' @details
#' The evaluation of the log-likelihood displacement is based in the analysis proposed in Verbeke and Molenberghs (2009) and Gałecki and Burzykowski (2013).
#' First, a list of models fitted to leave-one-subject-out datasets are obtained. Then, for each model, the maximum likelihood estimate obtained by fitting the 
#' model to all data and the maximum likelihood estimate obtained by fitting the model to the data with the \eqn{i}-th subject removed are obtained and used for the 
#' log-likelihood displacement calculation. The likelihood displacement, \eqn{LDi} , is defined as twice the difference between the log-likelihood computed at a 
#' maximum and displaced values of estimated parameters (Verbeke and Molenberghs (2009), Gałecki and Burzykowsk (2013)):
#' 
#' \deqn{LD_i \equiv 2 \times \Bigr[\ell_\textrm{Full}(\widehat{\Theta};\textrm{y})-\ell_\textrm{Full}(\widehat{\Theta}_{(-i)};\textrm{y})\Bigr]}
#' 
#' where \eqn{\widehat{\Theta}} is the maximum-likelihood estimate of \eqn{\Theta} obtained by fitting the model to all data, while \eqn{\widehat{\Theta}_{-i}} is
#' the maximum-likelihood estimate obtained by fitting the model to the data with the \eqn{i}-subject excluded.
#' 
#' @returns Returns a plot of the log-likelihood displacement values for each subject, indicating those subjects
#' whose contribution is greater than `disp_thrh`.
#' 
#' @references 
#' - Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' - Verbeke, G. & Molenberghs, G. (2000). _Linear Mixed Models for Longitudinal Data_. Springer New York. https://doi.org/10.1007/978-1-4419-0300-6
#' @examples
#' # Load the example data
#' data(grwth_data)
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
#' # Obtain log-likelihood displacement for each subject
#' logLikSubjectDisplacements(model = lmm)
#' # Modifying the threshold for log-likelihood displacement
#' logLikSubjectDisplacements(model = lmm, disp_thrh = 1)
#' 
#' # Calculating the log-likelihood contribution in a model with a variance structure specified
#' lmm_var <- lmmModel(
#'   data = grwth_data,
#'   sample_id = "subject",
#'   time = "Time",
#'   treatment = "Treatment",
#'   tumor_vol = "TumorVolume",
#'   trt_control = "Control",
#'   drug_a = "DrugA",
#'   drug_b = "DrugB",
#'   combination = "Combination",
#'   weights = nlme::varIdent(form = ~ 1|SampleID)
#'   ) 
#' # Calculate the log-likelihood contribution
#' logLikSubjectDisplacements(model = lmm, var_name = "SampleID")
#' 
#' @export
logLikSubjectDisplacements <- function(model,
                                       disp_thrh = NA,
                                       label_angle = 0,
                                       var_name = NULL,
                                       ...) {
  
  model <- update(model, method = "ML")
  
  # Fitting the model to "leave-one-out" data
  
  subject.c <- levels(model$data$SampleID)
  
  lmeUall <- lapply(subject.c, .lmeU, model = model)
  names(lmeUall) <- subject.c
  
  # Likelihood Displacement for Model
  
  lLikUall <- sapply(
    subject.c,
    .lLik,
    model = model,
    lmeUall = lmeUall,
    var_name = var_name
  ) #... for all subjects
  dif.2Lik <- 2 * (logLik(model) - lLikUall) # Vector for LDi
  
  # Plot of the likelihood displacements with an indication of outlying values
  
  if(is.na(disp_thrh)){
    disp_thrh <- round(quantile(dif.2Lik, probs = 0.9),3)
  }
  
  outL <- dif.2Lik > disp_thrh # Outlying LDi's
  
  if(sum(outL) == 0){
    writeLines(paste("No subject with a log-likelihood displacement greater than", disp_thrh))
  }
  
  print(paste(
    "Outliers with Log Likelihood displacement greater than:",
    disp_thrh
  ))
  print(dif.2Lik[outL])
  
  subject.f <- factor(subject.c, levels = subject.c)
  myPanel <- function(x, y, ...) {
    x1 <- as.numeric(x)
    panel.xyplot(x1, y, ...)
    ltext( # Label outlying LDi's
      x1[outL],
      y[outL],
      subject.c[outL],
      cex = 1,
      adj = c(0.5, 1),
      srt = label_angle
    ) 
    panel.xyplot(
      x = subject.f[outL],
      y = dif.2Lik[outL],
      type = "p",
      pch = 20
    )
    panel.abline(h = disp_thrh, lty = 2)
  }
  dtp <- dotplot(
    dif.2Lik ~ subject.f,
    panel = myPanel,
    type = "h",
    ylab = "Log Likelihood-displacement",
    xlab = "Subject",
    main = "Log Likelihood-displacement values vs Subjects rank"
  )
  lxlims <- length(dtp$x.limits)
  print(update(
    dtp,
    xlim = rep("", lxlims),
    grid = "h",
    key = list(
      lines = list(lty = 2, lwd = 0.5),
      text = list(c(
        paste("logLik displacement threshold:", as.character(disp_thrh))
      ), cex = 1.2),
      space = "top"
    )
  ))
  return(invisible(dif.2Lik))
}

