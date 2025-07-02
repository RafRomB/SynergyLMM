
## Cook's distance for beta estimates
#' @title Helper function for the calculation of Cook's distance of beta estimates
#' @param betaU A vector with beta estimates calculated for leave-one-out models.
#' @param beta0 A vector with beta estimates calculated for the complete data model.
#' @param vb.inv A matrix with the inverse variance-covariance matrix of beta estimates of the complete data model.
#' @returns A numeric value of the numerator of Cook's distance, as described in Gałecki, A., & Burzykowski, T. (2013).
#' @keywords internal
#' @noRd
.CookDfun <- function(betaU, beta0, vb.inv){
  dbetaU <- betaU - beta0 # beta(-i) - beta
  CookD.value <- t(dbetaU) %*% vb.inv %*% dbetaU # Compute Cook distance (Gałecki, A., & Burzykowski, T. (2013)
}

## Leave-one-out prediction

#' @title Helper function for calculating Cook's distance with leave-one-out data
#' @param cx Subject to remove from the data to build the model
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`],
#' and fitted using maximum likelihood.
#' @param maxIter Limit of maximum number of iterations for the optimization algorithm. Default to 1000. 
#' @returns Sum of squared errors via the leave-one-out model fits
#' @keywords internal
#' @noRd

.lmeU2 <- function(cx, model, maxIter){
  SampleID <- NULL
  dfU <- subset(model$data, SampleID != cx) ## LOO data
  maxIter.tmp <- 50
  repeat{
    tmp.update <- suppressWarnings(
      try(update(model, data = dfU, control = list(maxIter = maxIter.tmp)), silent = T)
    )
    if (any(class(tmp.update) == "try-error")) {
      maxIter.tmp <- 2 * maxIter.tmp
      if (maxIter.tmp > maxIter) {
        warning(paste0("Maximum number of iterations (", maxIter, ") exceeded for SampleID = ", cx,
                       ". Cook distance not calculated for SampleID = ", cx))
        return(0)
      }
    } else {
      break
    }
  }
  sum((predict(model, newdata = model$data, level = 0) - 
    predict(tmp.update, newdata = model$data, level = 0))^2)
}

#' @title Cook's distance for individual subjects
#' @description
#' `CookDistance` allows the user to identify those subjects with a greater influence in the predicted values or in the estimation of the
#' fixed effects for the treatment group, based in the calculation of Cook's distances.
#' 
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param type Type of Cook's distance to calculated. Possible options are `fitted`, to calculte Cook's distances
#' based on the change in fitted values, or `fixef` to calculate Cook's distances based on the change in the fixed effects.
#' See Details section for more information.
#' @param cook_thr Numeric value indicating the threshold for the Cook's distance. If not specified, the threshold is set to three times the mean of the Cook's
#' distance values.
#' @param label_angle Numeric value indicating the angle for the label of subjects with a Cook's distance greater than `cook_thr`.
#' @param maxIter Limit of maximum number of iterations for the optimization algorithm. Default to 1000. 
#' @param verbose Logical indicating if the subjects with a Cook's distance greater than `cook_thr` should be printed to the console.
#' @details
#' The identification of influential subjects is based on the calculation of Cook's distances. The Cook's distances
#' can be calculated based on the change in fitted values or fixed effects.
#' - **Cook's distances based on the change in fitted values** 
#' 
#' When `type = "fitted"`, the Cook's distances
#' are calculated as the normalized change in fitted response values due to the removal of a subject from the model.
#' Firts, a leave-one-subject-out model is fitted, removing individually each subject to fit the model. Then, the Cook's
#' distance for subject \eqn{i}, (\eqn{D_i}), is calculated as:
#' \deqn{D_i=\frac{\sum_{j=1}^n\Bigl(\hat{y}_{j}-\hat{y}_{j_{(-i)}}\Bigl)^2}{rank(X)\cdot MSE}}
#' 
#' where \eqn{\hat{y}_j} is the \eqn{j^{th}} fitted response value using the complete model, and \eqn{\hat{y}_{j_{(-i)}}} is the 
#' \eqn{j^{th}} fitted response value obtained using the model where subject \eqn{i} has been removed.
#' 
#' The denominator of the expression is equal to the number of the fixed-effects coefficients, which, under the assumption that
#' the design matrix is of full rank, is equivalent to the rank of the design matrix, and the Cook distance is normalized by the
#' mean square error (\eqn{MSE}) of the model. 
#' 
#' - **Cook's distances based on the change in fixed effects values**
#' 
#' The identification of the subjects with a greater influence in the estimated fixed effects is based on the calculation of Cook's distances, as
#' described in Gałecki and Burzykowsk (2013). To compute the Cook's distance for the fixed effect estimates (i.e., the contribution to each subject to the coefficients of its treatment group), 
#' first a matrix containing the leave-one-subject-out estimates or the fixed effects is calculated. Then, the Cook's distances are calculated according to:
#' 
#' \deqn{D_i \equiv  \frac{(\hat{\beta} - \hat{\beta}_{(-i)})[\widehat{Var(\hat{\beta})}]^{-1}(\hat{\beta} - \hat{\beta}_{(-i)})}{rank(X)}}
#' 
#' where \eqn{\beta} represents the vector of fixed effects and \eqn{\hat{\beta}_{(-i)}} is the estimate of the parameter vector \eqn{\beta} obtained by fitting the model to the data with the \eqn{i}-th subject excluded. The denominator of 
#' the expression is equal to the number of the fixed-effects coefficients, which, under the assumption that the design matrix is of full rank, is equivalent to the rank of the design matrix.
#'  
#' 
#' @returns A plot of the Cook's distance value for each subject, indicating those subjects
#' whose Cook's distance is greater than `cook_thr`. 
#' 
#' If saved to a variable, the function returns a vector with the Cook's distances for each subject.
#' 
#' @references 
#' - Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' @examples
#' #' # Load the example data
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
#' # Calulate Cook's distances for each subject
#' CookDistance(model = lmm)
#' # Change the Cook's distance threshold
#' CookDistance(model = lmm, cook_thr = 0.15)
#' 
#' @export
CookDistance <- function(model,
                         type = "fitted",
                         cook_thr = NA,
                         label_angle = 0,
                         maxIter = 1000,
                         verbose = TRUE) {
  subject.c <- levels(model$data$SampleID)
  
  beta0 <- fixef(model) # estimated betas
  n.fixeff <- length(beta0) # Number of fixed effects
  rankX <- n.fixeff # Rank of matrix X
  
  if (type == "fitted"){
    CookD.num <- sapply(subject.c , .lmeU2, model = model, maxIter = maxIter)
    
    CookD <- CookD.num / rankX / mean(model$residuals^2) # Cook's distance Di
    
  } else if (type == "fixef") {
    
    lmeUall <- lapply(subject.c, .lmeU, model = model, maxIter = maxIter)
    names(lmeUall) <- subject.c
    
    idx <- c()
    
    for (i in lmeUall){
      idx <- c(idx, (typeof(i)))
    }
    
    betaUall <- sapply(lmeUall[idx == "list"], fixef) # Matrix with betas(-i) estimates
    
    vcovb <- vcov(model) # Variance-covariance matrix for betas
    vb.inv <- solve(vcovb) # Inverse of of the var-cov matrix of betas
    
    CookD.num <- apply(betaUall, 2, .CookDfun, beta0 = beta0, vb.inv = vb.inv)
    CookD <- CookD.num / rankX # Cook's distance Di
  } else {
    stop(paste("type = ", type, " is not a valid option for Cook's distance calculation."))
  }
  
  if(is.na(cook_thr)){
    cook_thr <- round(3*mean(CookD[CookD != 0]),3)
  }
  
  outD <- CookD > cook_thr # Outlying Di's
  
  if(sum(outD) == 0){
      message(paste("No subject with a Cook's distance greater than:", cook_thr))
  } else {
    if (verbose){
      print(paste("Subjects with Cook's distance greater than:", cook_thr))
      print(CookD[outD]) 
    }
  }
  
  subject.x <- 1:length(subject.c)
  plot(
    CookD,
    ylab = "Cook's Distance",
    type = "h",
    xlab = "Subject",
    main = "Cook's Distance vs Subject"
  )
  
  if (sum(outD) > 0) {
    points(subject.x[outD], CookD[outD], pch = 20)
    text(
      x = subject.x[outD],
      y = CookD[outD],
      labels = subject.c[outD],
      cex = 1,
      adj = c(0.5, 0.1),
      srt = label_angle,
      xpd = TRUE
    )
  }
  
  abline(h = cook_thr, lty = "dashed")
  legend(
    "topright",
    inset = c(0, -0.1),
    legend = c(as.character(cook_thr)),
    lty = 2,
    title = "Cook's distance threshold",
    bty = "n",
    cex = 0.8,
    xpd = TRUE
  )
  return(invisible(CookD))
}

