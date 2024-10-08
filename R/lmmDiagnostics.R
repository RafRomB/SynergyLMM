# Random Effects Diagnostics ----
#' @title Diagnostics of random effects of the linear mixed model
#' 
#' @description
#' `ranefDiagnostics` provides several plots as well as statistical test for the examination
#' of the normality of the random effects of the input model.
#' 
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @import ggplot2
#' @details
#' One of the assumptions of the model obtained with [`lmmModel()`] (as in any other linear mixed model) is that
#' the random effects are normally distributed:
#' 
#' \deqn{b_i = N(0,\psi)}
#' 
#' For the evaluation of this assumption, `ranefDiagnostics` provides Q-Q plots of the normalized residuals, 
#' together with statistical assessment of their normality using Shapiro-Wilk, D'Agostini and Anderson-Darling normality tests. 
#' Additionally, Q-Q plots of the normalized residuals by time point and treatment group are provided to be able to detect time 
#' points or treatment groups which could be notably different from the others and be affecting the adequacy of the model. 
#' Additionally, scatter plots of the conditional Pearson residuals versus fitted values and Pearson residuals per time and per 
#' treatment are provided to give information about variability of the residuals and possible outlier observations. 
#' Observations with absolute standardized (normalized) residuals greater than the \eqn{1-0.05/2} quantile of the standard normal distribution 
#' are identified and reported as potential outlier observations.
#' 
#' @returns A list with different elements for the diagnostics of the random effects are produced:
#' - `plots`: Different plots for evaluating the normality and homoscedasticity of the random effects are produced.
#' - `Normality`: List with the results from 3 different test of the normality of the random effects: Shapiro - Wilk normality test, 
#' D'Agostino normality test and Anderson - Darling normality test.
#' - `Levene.test`: results from Levene homocedasticity test of the conditional, marginal and normalized residuals by SampleID (i.e. by subject).
#' - `Fligner.test`: results from Fligner homocedasticity test of the conditional, marginal and normalized residuals by SampleID (i.e. by subject).
#' 
#' @references
#' - Pinheiro JC, Bates DM (2000). _Mixed-Effects Models in S and S-PLUS_. Springer, New York. doi:10.1007/b98882 <https://doi.org/10.1007/b98882>.
#' - Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' 
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
#'   drug_ab = "Combination"
#'   )
#' # Run random effects diagnostics
#' ranef_diag <- ranefDiagnostics(lmm)
#' 
#' #Access to individual plots
#' 
#' ranef_diag$Plots[1]
#' ranef_diag$Plots[2]
#' 
#' # Access to normality tests
#' 
#' ranef_diag$Normality
#' 
#' ranef_diag$Normality$Shapiro.test
#' 
#' # Access to homoscedasticity tests of residuals by subject
#' 
#' ranef_diag$Levene.test
#' 
#' ranef_diag$Fligner.test
#' 
#' @export

ranefDiagnostics <- function(model){
  # Plots
  ranef_plot <- plot_ranefDiagnostics(model)
  print(ranef_plot[[5]])
  
  # Normality test
  print(
    ranef_shapiro <- fBasics::shapiroTest(nlme::ranef(model)$Time, 
                                          title = "Shapiro - Wilk Normality Test of random effects")
  )
  
  if (length(nlme::ranef(model)$Time) < 20) {
    ranef_DAgostino <- warning("Sample size must be at least 20 for D'Agostino Normality Test")
  } else {
    print(
      ranef_DAgostino <- fBasics::dagoTest(nlme::ranef(model)$Time, 
                                           title = "D'Agostino Normality Test of random effects")
    )
  }
  print(
    ranef_ad <- fBasics::adTest(nlme::ranef(model)$Time, 
                                title = "Anderson - Darling Normality Test of random effects")
  )
  
  Normality <- list(
    Shapiro.test = ranef_shapiro,
    DAgostino.test = ranef_DAgostino,
    Anderson.Darling.test = ranef_ad
  )
  
  # Homocedasticity test
  
  # Conditional Residuals
  cond_res <- residuals(model, type = "response")
  cond_res <- data.frame(
    conditional_resid = cond_res,
    SampleID = names(cond_res),
    stringsAsFactors = T
  )
  colnames(cond_res) <- c("conditional_resid", "SampleID")
  
  # Marginal Residuals
  mar_res <- residuals(model, type = "response", level = 0)
  mar_res <- data.frame(
    marginal_resid = mar_res,
    SampleID = names(mar_res),
    stringsAsFactors = T
  )
  colnames(mar_res) <- c("marginal_resid", "SampleID")
  
  # Normalized Residuals
  
  norm_res <- residuals(model, type = "normalized")
  norm_res <- data.frame(
    normalized_resid = norm_res,
    SampleID = names(norm_res),
    stringsAsFactors = T
  )
  colnames(norm_res) <- c("normalized_resid", "SampleID")
  
  levene <- list()
  writeLines("Conditional Residuals")
  print(levene$conditional_resid <- car::leveneTest(conditional_resid ~ SampleID, data = cond_res))
  writeLines("Marginal Residuals")
  print(levene$marginal_resid <- car::leveneTest(marginal_resid ~ SampleID, data = mar_res))
  writeLines("Normalized Residuals")
  print(levene$normalized_resid <- car::leveneTest(normalized_resid ~ SampleID, data = norm_res))
  
  fligner <- list()
  
  writeLines("Conditional Residuals")
  print(fligner$conditional_resid <- fligner.test(conditional_resid ~ SampleID, data = cond_res))
  writeLines("Marginal Residuals")
  print(fligner$marginal_resid <- fligner.test(marginal_resid ~ SampleID, data = mar_res))
  writeLines("Normalized Residuals")
  print(fligner$normalized_resid <- fligner.test(normalized_resid ~ SampleID, data = norm_res))
  
  return(invisible(
    list(
      Plots = ranef_plot,
      Normality = Normality,
      Levene.test = levene,
      Fligner.test = fligner
    )
  ))
}

# Residual Diagnostics ----

#' @title Diagnostics of residuals of the linear mixed model
#' 
#' @description
#' `residDiagnostics` provides several plots as well as statistical test for the examination
#' of the normality of the residuals of the input model.
#' 
#' @param model  An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param pvalue Threshold for the p-value of outlier observations based on their Pearson residuals.
#' 
#' @details
#' One of the assumption of the model fit by [`lmmModel()`] is that the residuals are normally distributed:
#' 
#' \deqn{\varepsilon_{i} \sim \mathcal{N}(0, \sigma^2)}
#' 
#' For the evaluation of this assumption, `residDiagnostics` provides Q-Q plots of the normalized residuals, together with statistical assessment of their 
#' normality using Shapiro-Wilk, D'Agostini and Anderson-Darling normality tests. Additionally, Q-Q plots of the normalized residuals by time point and 
#' treatment group are provided to be able to detect time points or treatment groups which could be notably different from the others and be 
#' affecting the adequacy of the model. Additionally, scatter plots of the conditional Pearson residuals versus fitted values and Pearson residuals 
#' per time and per treatment are provided to give information about variability of the residuals and possible outlier observations. 
#' Observations with absolute standardized (normalized) residuals greater than the \eqn{1-0.05/2} quantile of the standard normal distribution 
#' are identified and reported as potential outlier observations.
#' 
#' @returns A list with different elements for the diagnostics of the residuals are produced:
#' - `plots`: Different plots for evaluating the normality and homocedasticity of the residuals.
#' - `outliers`: Data frame with the identified outliers based on the Pearson residuals and the value of `pval`. The column `resid.p` contains the
#' value of the Pearson residuals for each observation.
#' - `Normality`: List with the results from 3 different test of the normality of the normalized residuals of the model: Shapiro - Wilk normality test, 
#' D'Agostino normality test and Anderson - Darling normality test.
#' 
#' @references
#' - Pinheiro JC, Bates DM (2000). _Mixed-Effects Models in S and S-PLUS_. Springer, New York. doi:10.1007/b98882 <https://doi.org/10.1007/b98882>.
#' - Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' 
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
#'   drug_ab = "Combination"
#'   )
#' # Residuals diagnostics
#' resid_diag <- residDiagnostics(model = lmm, pvalue = 0.05)
#' # Access outliers data frame
#' resid_diag$outliers
#' # Access individual plots
#' resid_diag$plots[1]
#' resid_diag$plots[2]
#' 
#' # Access results of normality tests
#' resid_diag$Normality
#' resid_diag$Normality$Shapiro.test
#'
#' @export
residDiagnostics <- function(model, pvalue=0.05){
  # Plots
  resid_plot <- plot_residDiagnostics(model)
  print(resid_plot[[6]])
  
  # Normality test
  
  norm_res <- resid(model, type = "normalized")
  
  print(
    res_shapiro <- fBasics::shapiroTest(norm_res, description = "Shapiro - Wilk Normality Test of normalized residuals")
  )
  print(
    res_DAgostino <- fBasics::dagoTest(norm_res, description = "D'Agostino Normality Test of normalized residuals")
  )
  print(
    res_ad <- fBasics::adTest(norm_res, description = "Anderson - Darling Normality Test of normalized residuals")
  )
  
  Normality <- list(
    Shapiro.test = res_shapiro,
    DAgostino.test = res_DAgostino,
    Anderson.Darling.test = res_ad
  )
  
  
  # List of outliers
  outliers.idx <- within(model$data, {
    resid.p <- resid(model, type = "pearson") # Pearson resids.
    idx <- abs(resid.p) > -qnorm(pvalue / 2) # Indicator vector
  })
  print("Outlier observations")
  outliers <- subset(outliers.idx, idx) # Data with outliers
  outliers <- outliers[,-8]
  print(outliers)
  
  return(invisible(
    list(
      plots = resid_plot,
      outliers = outliers,
      Normality = Normality
    )
  ))
}

# Observed vs predicted ----

#' @title Observed vs predicted values and performance of the model
#' 
#' @description
#' `ObsvsPred` allows the user to have a straight forward idea about how the model is fitting the data, providing
#' plots of the predicted regression lines versus the actual data points.
#' 
#' 
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param nrow Number of rows of the layout to organize the observed vs predicted plots.
#' @param ncol Number of columns of the layout to organize the observed vs predicted plots.
#' @param ... Additional arguments to be passed to [performance::model_performance()].
#' @details
#'  The function provides visual and quantitative information about the performance of the model:
#' - A layout of the observed and predicted values of \eqn{log}(relative tumor volume) vs Time for each SampleID (i.e. subject), 
#' with the actual measurements, the regression line for each SampleID, and the marginal, treatment-specific, 
#' regression line for each treatment group.
#' - Performance metrics of the model obtain calculated using [performance::model_performance()]. The maximum likelihood-based Akaike's Information Criterion (AIC), 
#' small sample AIC (AICc), and Bayesian Information Criterion, and the Nakagawa's r-squared 
#' root mean squared error (RMSE) of the residuals, and the standard deviation of the residuals (sigma) are provided.
#' 
#' @returns A layout of the observed vs predicted values for each SampleID and model performance metrics.
#' 
#' @references
#' - Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' - Lüdecke et al., (2021). _performance: An R Package for Assessment, Comparison and Testing of Statistical Models_. Journal of Open Source Software, 6(60), 3139. https://doi.org/10.21105/joss.03139
#' - Sakamoto, Y., M. Ishiguro, and G. Kitagawa. 1984. _Akaike Information Criterion Statistics_. Mathematics and Its Applications. Reidel.
#' - Nakagawa, Shinichi, and Holger Schielzeth. 2013. _A General and Simple Method for Obtaining r2 from Generalized Linear Mixed-effects Models_. Methods in Ecology and Evolution 4 (February): 133–42. https://doi.org/10.1111/j.2041-210x.2012.00261.x.
#' - Johnson, Paul C. D. 2014. _Extension of Nakagawa & Schielzeth’s r 2 GLMM to Random Slopes Models_. Methods in Ecology and Evolution 5 (September): 944–46. https://doi.org/10.1111/2041-210X.12225.
#' - Nakagawa, Shinichi, Paul C. D. Johnson, and Holger Schielzeth. 2017. _The Coefficient of Determination r2 and Intra-Class Correlation Coefficient from Generalized Linear Mixed-Effects Models Revisited and Expanded_. Journal of The Royal Society Interface 14 (September): 20170213. https://doi.org/10.1098/rsif.2017.0213.
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
#'   drug_ab = "Combination"
#'   )
#'# Obtain Observed vs Predicted plots, and model performance metrics
#' ObsvsPred(model = lmm, nrow = 4, ncol = 8)
#' @export
ObsvsPred <- function(model, nrow = 4, ncol = 5, ...) {
  print(performance::model_performance(model, metrics = c("AIC", "AICc", "BIC", "R2", "RMSE", "SIGMA")), ...)
  plot_ObsvsPred(model, nrow, ncol)
}

# Influential Diagnostics ----

#' @title Modified [nlmeU::logLik1] helper function to deal with varIdent structure
#' @param modfit An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`],
#' with a variance structure defined by [nlme::varIdent], and fitted using maximum likelihood.
#' @param dt1 A data frame with data for one subject, for whom the log-likelihood function is to be evaluated
#' @param var_name Name of the variable for the weights of the model in the case that a variance structure has been specified using [nlme::varIdent()].
#' @returns Numeric scalar value representing contribution of a given subject to the overall 
#' log-likelihood returned by `logLik()` function applied to the "lme" object defined by `modfit` argument.
#' @keywords internal
#' @noRd

.logLik1.varIdent <- function(modfit, dt1, var_name){
  var_name <- as.character(unique(dt1[,var_name]))
  m <- modfit$modelStruct # Model structure
  sigma <- modfit$sigma # sigma
  D <- as.matrix(m$reStruct[[1]]) # "subject"
  D <- D  * sigma^2 # Matrix D 
  
  vecR <- sigma/(nlme::varWeights(modfit$modelStruct$varStruct)) # AugDiagonal of R_i
  if(length(var_name) == 1){
    vecR <- vecR[names(vecR) == var_name]
  } else{
    vecR <- vecR[var_name]
  }
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

## Individual contributions to logLik 

#' @title Contributions of the individual subjects to the log-likelihood for the model.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param lLik_thrh Threshold of log-likelihood contribution per-observation per-subject.
#' @param label_angle Angle for the label of subjects with a log-likelihood contribution smaller than `lLik_thrh`.
#' @param var_name Name of the variable for the weights of the model in the case that a variance structure has been specified using [nlme::varIdent()].
#' @returns Returns a plot of the per-observation individual-subject log-likelihood contibutions to the model, indicating those subjects
#' whose contribution is smaller than `lLik_thrh`. 
#' @export
logLikSubjectContributions <- function(model,
                        lLik_thrh = 0,
                        label_angle = 0,
                        var_name = NULL) {
  
  # Update model to use maximum-likelihood estimation
  
  model <- update(model, method = "ML")
  
  # Fixed-effect estimates and their variance-covariance matrix
  
  beta0 <- nlme::fixef(model) # estimated betas
  
  vcovb <- vcov(model) # Variance-covariance matrix for betas
  
  # Contribution of Individual Subjects to the Log-Likelihood
  
  if(is.null(model$modelStruct$varStruct)) {
    lLik.i <- by(
      model$data,
      model$data$SampleID,
      FUN = function(dfi) # Function to calculate log likelihood for subject i
        nlmeU::logLik1(model, dfi)
    ) 
  } else {
    if (is.null(var_name)) {
      stop("`var_name` cannot be NULL if a variance estructure has been specified in the model")
    }
    lLik.i <- by(
      model$data,
      model$data$SampleID,
      FUN = function(dfi) # Function to calculate log likelihood for subject i with a variance structure.
        .logLik1.varIdent(model, dfi, var_name = var_name)
    ) 
  }
  
  lLik.i <- as.vector(lLik.i) # Coerce array to a vector
  
  # Plot of individual contributions to the log-likelihood (traditional graphics)
  
  nx <- by(model$data, model$data$SampleID, nrow) # ni
  lLik.n <- lLik.i / as.vector(nx) # logLiki/ni
  outL <- lLik.n < lLik_thrh # TRUE for value < lLik_thrh
  lLik.n[outL] # logLiki/ni < lLik_thrh
  
  if (sum(outL) == 0) {
    writeLines(paste("No subject with a log-likelihood contribution below", lLik_thrh))
  }
  
  subject.c <- levels(model$data$SampleID)
  subject.x <- 1:length(subject.c)
  plot(
    lLik.n ~ subject.x,
    type = "h",
    ylab = "Contribution to log Likelihood",
    xlab = "Subject",
    main = "Per-observation subject\nlog-likelihood contributions"
  )
  points(subject.x[outL], lLik.n[outL], type = "p", pch = 20)
  text(
    x = subject.x[outL],
    y = lLik.n[outL],
    labels = subject.c[outL],
    cex = 0.8,
    adj = c(0.5, 1),
    srt = label_angle,
    xpd = TRUE
  )
  abline(h = lLik_thrh, lty = 2)
  legend(
    "topright",
    inset = c(-0.02, -0.1),
    legend = c(as.character(lLik_thrh)),
    lty = 2,
    title = "logLik contribution\nthreshold",
    bty = "n",
    cex = 0.66,
    xpd = TRUE
  )
  names(lLik.n) <- subject.c
  print(lLik.n)
}

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
  D <- as.matrix(m$reStruct[[1]]) # "subject"
  D <- D  * sigma^2 # Matrix D 
  
  vecR <- sigma/(nlme::varWeights(dtInit$modelStruct$varStruct)) # AugDiagonal of R_i
  if(length(var_name) == 1){
    vecR <- vecR[names(vecR) == var_name]
  } else{
    vecR <- vecR[var_name]
  }  
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
#' - Molenberghs, G., & Verbeke, G. (2000). _Linear Mixed Models for Longitudinal Data_. Springer New York. https://doi.org/10.1007/978-1-4419-0300-6
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
#'   drug_ab = "Combination"
#'   ) 
#' # Obtain log-likelihood displacement for each subject
#' logLikSubjectDisplacements(model = lmm)
#' # Modifying the threshold for log-likelihood displacement
#' logLikSubjectDisplacements(model = lmm, disp_thrh = 1.5)
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
#'   drug_ab = "Combination",
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
      cex = 0.66,
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
      ), cex = 0.66),
      space = "top"
    )
  ))
  return(invisible(dif.2Lik))
}


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
#'   drug_ab = "Combination"
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
      cex = 0.6,
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
    title = "Cook's distance\n threshold",
    bty = "n",
    cex = 0.66,
    xpd = TRUE
  )
  return(invisible(CookD))
}



