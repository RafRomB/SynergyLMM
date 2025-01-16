
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

#' @title Cook's distance for the coefficient estimates
#' @description
#' `CookDistance` allows the user to identify those subjects with a greater influence in the estimation of the
#' \eqn{\beta} (tumor growth rate) for the treatment group, based in the calculation of Cook's distances.
#' 
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param cook_thr Numeric value indicating the threshold for the Cook's distance. If not specified, the threshold is set to the 90% percentile of the Cook's
#' distance values.
#' @param label_angle Numeric value indicating the angle for the label of subjects with a Cook's distance greater than `cook_thr`.
#' @details
#' The identification of the subjects with a greater influence in each estimated \eqn{\beta} representing the tumor growth is based on the calculation of Cook's distances, as
#' described in Gałecki and Burzykowsk (2013). To compute the Cook's distance for the \eqn{\beta} estimates (i.e., the contribution to each subject to the coefficient of its treatment group), 
#' first a matrix containing the leave-one-subject-out estimates or \eqn{\beta} is calculated. Then, the Cook's distances are calculated according to:
#' 
#' \deqn{D_i \equiv  \frac{(\hat{\beta} - \hat{\beta}_{(-i)})[\widehat{Var(\hat{\beta})}]^{-1}(\hat{\beta} - \hat{\beta}_{(-i)})}{rank(X)}}
#' 
#' where \eqn{\hat{\beta}_{(-i)}} is the estimate of the parameter vector \eqn{\beta} obtained by fitting the model to the data with the \eqn{i}-th subject excluded. The denominator of 
#' the expression is equal to the number of the fixed-effects coefficients, which, under the assumption that the design matrix is of full rank, is equivalent to the rank of the design matrix.
#' 
#' @returns Plot of the Cook's distance value for each subject, indicating those subjects
#' whose Cook's distance is greater than `cook_thr`.
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
                         cook_thr = NA,
                         label_angle = 0) {
  subject.c <- levels(model$data$SampleID)
  lmeUall <- lapply(subject.c, .lmeU, model = model)
  names(lmeUall) <- subject.c
  
  betaUall <- sapply(lmeUall, fixef) # Matrix with betas(-i) estimates
  
  vcovb <- vcov(model) # Variance-covariance matrix for betas
  vb.inv <- solve(vcovb) # Inverse of of the var-cov matrix of betas
  
  beta0 <- fixef(model) # estimated betas
  
  CookD.num <- apply(betaUall, 2, .CookDfun, beta0 = beta0, vb.inv = vb.inv)
  
  n.fixeff <- length(beta0) # Number of fixed effects
  rankX <- n.fixeff # Rank of matrix X
  CookD <- CookD.num / rankX # Cook's distance Di
  
  if(is.na(cook_thr)){
    cook_thr <- round(quantile(CookD, probs = 0.9),3)
  }
  
  outD <- CookD > cook_thr # Outlying Di's
  
  if(sum(outD) == 0){
    writeLines(paste("No subject with a Cook's distance greater than", cook_thr))
  }
  
  print(paste("Subjects with Cook's distance greater than:", cook_thr))
  print(CookD[outD])
  
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

