#' @title Diagnostics of random effects of the linear mixed model.
#' 
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`LMM_Model()`].
#' @import ggplot2
#' @import stats
#' @import graphics
#' @returns A list with different elements for the diagnostics of the random effects are produced:
#' - `plots`: Different plots for evaluating the normality and homoscedasticity of the random effects are produced.
#' - `Normality`: List with the results from 3 different test of the normality of the random effects: Shapiro - Wilk normality test, 
#' D'Agostino normality test and Anderson - Darling normality test.
#' - `Levene.test`: results from Levene homocedasticity test of the conditional, marginal and normalized residuals by mouse.
#' - `Fligner.test`: results from Fligner homocedasticity test of the conditional, marginal and normalized residuals by mouse.
#' @export

ranef_diag <- function(model){
  # Plots
  p1 <- ggplot(nlme::ranef(model), aes(sample = nlme::ranef(model)$Day)) + stat_qq(col = "gray20") + stat_qq_line() +
    labs(title = "Normal Q-Q Plot of Random Effects") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + cowplot::theme_cowplot()+
    theme(plot.title = element_text(size = 10, hjust = 0.5), axis.title = element_text(size = 12))
  p2 <- qqnorm(model, ~resid(.)|Mouse, pch=20, cex = 0.5, col = "gray20",
               main = list("Normal Q-Q Plot by Mouse", cex = 0.8), par.strip.text=list(col="black", cex=.66))
  p3 <- plot(model, Mouse ~ resid(., type = "response"), abline = 0, main = list("Conditional Residuals by Mouse", cex = 0.8))
  p4 <- plot(model, residuals(., type = "pearson") ~ fitted(.)|Mouse, id = 0.05, adj = -0.03, pch = 20, col = "slateblue4", cex=0.5,
             main = list("Pearson Residuals vs Fitted Values by Mouse", cex =0.8),par.strip.text=list(col="black", cex=.66))
  cowplot::plot_grid(p1,p2,p3,p4, ncol = 2)
  ranef_plot <- cowplot::plot_grid(p1,p2,p3,p4, ncol = 2)
  
  # Normality test
  ranef_shapiro <- fBasics::shapiroTest(nlme::ranef(model)$Day, description = "Shapiro - Wilk Normality Test of random effects")
  
  if(length(nlme::ranef(model)$Day)<20){
    ranef_DAgostino <- print("Sample size must be at least 20 for D'Agostino Normality Test")
  } else {
    ranef_DAgostino <- fBasics::dagoTest(nlme::ranef(model)$Day, description = "D'Agostino Normality Test of random effects") 
  }
  ranef_ad <- fBasics::adTest(nlme::ranef(model)$Day, description = "Anderson - Darling Normality Test of random effects")
  
  Normality <- list(Shapiro.test = ranef_shapiro, DAgostino.test = ranef_DAgostino, Anderson.Darling.test = ranef_ad)
  
  # Homocedasticity test
  
  # Conditional Residuals
  cond_res <- residuals(model, type = "response")
  cond_res <- data.frame(Conditional.Resid = cond_res, Mouse = names(cond_res), stringsAsFactors = T)
  colnames(cond_res) <- c("Conditional.Resid", "Mouse")
  
  # Marginal Residuals
  mar_res <- residuals(model, type = "response", level = 0)
  mar_res <- data.frame(Marginal.Resid = mar_res, Mouse = names(mar_res), stringsAsFactors = T)
  colnames(mar_res) <- c("Marginal.Resid", "Mouse")
  
  # Normalized Residuals
  
  norm_res <- residuals(model, type = "normalized")
  norm_res <- data.frame(Normalized.Resid = norm_res, Mouse = names(norm_res), stringsAsFactors = T)
  colnames(norm_res) <- c("Normalized.Resid", "Mouse")
  
  levene <- list()
  
  levene$Conditional.Resid <- car::leveneTest(Conditional.Resid ~ Mouse, data = cond_res)
  levene$Marginal.Resid <- car::leveneTest(Marginal.Resid ~ Mouse, data = mar_res)
  levene$Normalized.Resid <- car::leveneTest(Normalized.Resid ~ Mouse, data = norm_res)
  
  fligner <- list()
  
  fligner$Conditional.Resid <- fligner.test(Conditional.Resid ~ Mouse, data = cond_res)
  fligner$Marginal.Resid <- fligner.test(Marginal.Resid ~ Mouse, data = mar_res)
  fligner$Normalized.Resid <- fligner.test(Normalized.Resid ~ Mouse, data = norm_res)
  
  return(list(plots = ranef_plot, Normality = Normality, Levene.test = levene, Fligner.test = fligner))
}


#' @title Diagnostics of residuals of the linear mixed model.
#' 
#' @param model  An object of class "lme" representing the linear mixed-effects model fitted by [`LMM_Model()`].
#' @param pval Threshold for the p-value of outlier observations based on their Pearson residuals.
#' @returns A list with different elements for the diagnostics of the residuals are produced:
#' - `plots`: Different plots for evaluating the normality and homocedasticity of the residuals.
#' - `outliers`: Data frame with the identified outliers based on the Pearson residuals and the value of `pval`.
#' - `Normality`: List with the results from 3 different test of the normality of the normalized residuals of the model: Shapiro - Wilk normality test, 
#' D'Agostino normality test and Anderson - Darling normality test.
#' @export
resid_diag <- function(model, pval=0.05){
  # Plots
  p1 <- qqnorm(model, ~resid(., type = "normalized"),pch = 20, main = "Q-Q Plot of Normalized Residuals")
  p2 <- qqnorm(model, ~resid(., type = "normalized")|Day,pch = 20, main = "Q-Q Plot of Normalized Residuals by Day",
               par.strip.text=list(col="black", cex=.5))
  p3 <- qqnorm(model, ~resid(., type = "normalized")|Treatment, pch = 20, main = "Q-Q Plot of Normalized Residuals by Treatment",
               par.strip.text=list(col="black", cex=.5))
  p4 <- plot(model,main = "Pearson Residuals vs Fitted Values", pch = 20)
  p5 <- plot(model, resid(., type = "pearson") ~Day|Treatment, id = 0.05, pch=20,
             adj = -0.03, cex = 0.6, main = "Pearson Residuals per Time and Treatment", 
             par.strip.text=list(col="black", cex=.5))
  
  cowplot::plot_grid(p1,p2,p3,p4,p5, nrow = 3, ncol = 2, align = "hv")
  resid_plot <-  cowplot::plot_grid(p1,p2,p3,p4,p5, nrow = 3, ncol = 2, align = "hv")
  
  # Normality test
  
  norm_res <- resid(model, type = "normalized")
  
  res_shapiro <- fBasics::shapiroTest(norm_res, description = "Shapiro - Wilk Normality Test of normalized residuals")
  res_DAgostino <- fBasics::dagoTest(norm_res, description = "D'Agostino Normality Test of normalized residuals")
  res_ad <- fBasics::adTest(norm_res, description = "Anderson - Darling Normality Test of normalized residuals")
  
  Normality <- list(Shapiro.test = res_shapiro, DAgostino.test = res_DAgostino, Anderson.Darling.test = res_ad)
  
  
  # List of outliers
  id <- pval # Argument for qnorm
  outliers.idx <- within(model$data, 
                         {
                           resid.p <- resid(model, type = "pearson") # Pearson resids.
                           idx <- abs(resid.p) > -qnorm(id/2) # Indicator vector
                         })
  outliers <- subset(outliers.idx, idx) # Data with outliers
  
  return(list(plots = resid_plot, outliers = outliers, Normality = Normality))
}

## Observed vs predicted ----

#' @title Observed vs predicted values and performance
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`LMM_Model()`].
#' @param nrow Number of rows of the layout to organize the observed vs predicted plots.
#' @param ncol Number of columns of the layout to organize the observed vs predicted plots.
#' @param ... Additional arguments to be passed to [performance::model_performance()].
#' @details
#'  The function provides visual and quantitative information about the performance of the model:
#' - A layout of the observed and predicted values of log(relative tumor volume) vs day for each mouse, 
#' with the actual measurements, the regression line for each mouse, and the marginal, treatment-specific, 
#' regression line for each treatment group.
#' - Performance metrics of the model obtain calculated using [performance::model_performance()].
#' @returns A layout of the observed vs predicted values for each mouse and model performance metrics.
#' @export
ObsvsPred_plot <- function(model, nrow, ncol, ...){
  print(performance::model_performance(model, metrics = c("AIC", "AICc","BIC", "R2", "RMSE", "SIGMA")), ...)
  TV.df <- model$data
  aug.Pred <- nlme::augPred(model, primary = ~Day, level = 0:1, length.out = 2)
  plot(aug.Pred, layout = c(ncol, nrow, 1), lty = c(1,2),
       key = list(lines = list(lty = c(1,2), col = c("slateblue", "orange")),
                  text = list(c("Marginal", "Mouse-specific")),
                  columns = 2,
                  space="top"), 
       pch = 20, lwd = 1.5, main = "Observed and Predicted Values by Day",
       par.strip.text=list(col="black", cex=.5))
}

# Influential Diagnostics

# Modified logLik1 function to deal with varIdent structure

logLik1.varIdent <- function(modfit, dt1, varName){
  varName <- as.character(unique(dt1[,varName]))
  m <- modfit$modelStruct # Model structure
  sigma <- modfit$sigma # sigma
  D <- as.matrix(m$reStruct[[1]]) # "subject"
  D <- D  * sigma^2 # Matrix D 
  
  vecR <- sigma/(nlme::varWeights(modfit$modelStruct$varStruct)) # AugDiagonal of R_i
  if(length(varName) == 1){
    vecR <- vecR[names(vecR) == varName]
  } else{
    vecR <- vecR[varName]
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
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`LMM_Model()`].
#' @param lLik_thrh Threshold of log-likelihood contribution.
#' @param label_angle Angle for the label of subjects with a log-likelihood contribution smaller than `lLik_thrh`.
#' @param varName Name of the variable for the weights of the model in the case that a variance structure has been specified using [nlme::varIdent()].
#' @returns Returns a plot of the per-observation individual-subject log-likelihood contibutions to the model, indicating those subjects
#' whose contribution is smaller than `lLik_thrh`. 
#' @export
logLik_cont <- function(model, lLik_thrh, label_angle = 0, varName = NULL){
  
  # Fixed-effect estimates and their variance-covariance matrix
  
  beta0 <- nlme::fixef(model) # estimated betas
  
  vcovb <- vcov(model) # Variance-covariance matrix for betas
  
  # Contribution of Individual Subjects to the Log-Likelihood
  
  if(is.null(model$modelStruct$varStruct)){
    lLik.i <- by(model$data, model$data$Mouse, FUN = function(dfi) nlmeU::logLik1(model, dfi)) # Function to calculate log likelihood for subject i
  } else {
    if(is.null(varName)){
      stop("`varName` cannot be NULL is a variance estructure has been specified in the model")
    }
    lLik.i <- by(model$data, model$data$Mouse, FUN = function(dfi) logLik1.varIdent(model, dfi, varName = varName)) # Function to calculate log likelihood for subject i
  }
  
  lLik.i <- as.vector(lLik.i) # Coerce array to a vector
  
  # Plot of individual contributions to the log-likelihood (traditional graphics)
  
  nx <- by(model$data, model$data$Mouse, nrow) # ni
  lLik.n <- lLik.i/as.vector(nx) # logLiki/ni
  outL <- lLik.n < lLik_thrh # TRUE for value < lLik_thrh
  lLik.n[outL] # logLiki/ni < lLik_thrh
  subject.c <- levels(model$data$Mouse) 
  subject.x <- 1:length(subject.c)
  plot(lLik.n~subject.x, type = "h", ylab = "Contribution to log Likelihood",
       xlab = "Subject", main = "Per-observation subject\nlog-likelihood contributions")
  points(subject.x[outL], lLik.n[outL], type = "p", pch = 20)
  text(x = subject.x[outL], y = lLik.n[outL], labels = subject.c[outL], cex = 0.8, adj = c(0.5,0.1), srt = label_angle, xpd = TRUE)
  abline(h = lLik_thrh, lty = 2)
  legend("topright", inset = c(-0.02,-0.1), legend = c(as.character(lLik_thrh)),lty = 2,
         title = "logLik contribution\nthreshold", bty = "n", cex = 0.66,xpd = TRUE)
}

## Leave-one-out 

# Creating object lmeUall containing fitted models

lmeU <- function(cx, model){
  Mouse <- NULL
  dfU <- subset(model$data, Mouse != cx) ## LOO data
  update(model, data = dfU)
}

# logLik1 function for varStruct with LOO data

logLik1.varIdent_loo <- function(modfit, dt1, dtInit, varName){
  varName <- as.character(unique(dt1[,varName]))
  m <- modfit$modelStruct # Model structure
  sigma <- modfit$sigma # sigma
  D <- as.matrix(m$reStruct[[1]]) # "subject"
  D <- D  * sigma^2 # Matrix D 
  
  vecR <- sigma/(nlme::varWeights(dtInit$modelStruct$varStruct)) # AugDiagonal of R_i
  if(length(varName) == 1){
    vecR <- vecR[names(vecR) == varName]
  } else{
    vecR <- vecR[varName]
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

lLik <- function(cx, model, lmeUall, varName){
  Mouse <- NULL
  lmeU <- lmeUall[[cx]] # LOO fit extracted
  lLikU <- logLik(lmeU, REML = FALSE) # LOO log-likelihood
  df.s <- subset(model$data, Mouse == cx) # Data for subject cx...
  if(is.null(model$modelStruct$varStruct)){
    lLik.s <- nlmeU::logLik1(lmeU, df.s) # ... and log-likelihood  
  } else{
    if(is.null(varName)){
      stop("`varName` cannot be NULL is a variance estructure has been specified in the model")
    }
    lLik.s <- logLik1.varIdent_loo(lmeU, df.s, model, varName) # ... and log-likelihood
  }
  return(lLikU + lLik.s) # "Displaced" log-likelihood
} 

#' @title Likelihood displacement for the model.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`LMM_Model()`].
#' @param disp_thrh Threshold of log-likelihood displacement.
#' @param label_angle Angle for the label of subjects with a log-likelihood displacement greater than `disp_thrh`.
#' @param varName Name of the variable for the weights of the model in the case that a variance structure has been specified using [nlme::varIdent()].
#' @param ... Extra arguments, if any, for [lattice::panel.xyplot]. Usually passed on as 
#' graphical parameters to low level plotting functions, or to the panel functions performing smoothing, if applicable.
#' @import nlme
#' @import lattice
#' @returns Returns a plot of the log-likelihood displacement values for each subject, indicating those subjects
#' whose contribution is greater than `disp_thrh`.
#' @export
loo_logLik_disp <- function(model, disp_thrh, label_angle = 0, varName = NULL, ...){
  
  # Fitting the model to "leave-one-out" data
  
  subject.c <- levels(model$data$Mouse)
  
  lmeUall <- lapply(subject.c, lmeU, model = model)
  names(lmeUall) <- subject.c
  
  # Likelihood Displacement for Model
  
  lLikUall <- sapply(subject.c, lLik, model = model, lmeUall = lmeUall, varName = varName) #... for all subjects
  dif.2Lik <- 2*(logLik(model) - lLikUall) # Vector for LDi

  # Plot of the likelihood displacements with an indication of outlying values
  
  disp_thrh <- disp_thrh
  
  outL <- dif.2Lik > disp_thrh # Outlying LDi's
  print(paste("Outliers with Log Likelihood displacement greater than:", disp_thrh))
  print(dif.2Lik[outL])
  
  subject.f <- factor(subject.c, levels = subject.c)
  myPanel <- function(x, y, ...){
    x1 <- as.numeric(x)
    panel.xyplot(x1, y, ...)
    ltext(x1[outL], y[outL], subject.c[outL], cex = 0.66, adj = c(0.5,0.1), srt = label_angle) # Label outlying LDi's
    panel.xyplot(x = subject.f[outL], y = dif.2Lik[outL], type = "p", pch = 20)
    panel.abline(h = disp_thrh, lty = 2)
  }
  dtp <- dotplot(dif.2Lik~subject.f, panel = myPanel, type = "h", ylab = "Log Likelihood-displacement",
                 xlab = "Subject", main = "Log Likelihood-displacement values vs Subjects rank")
  lxlims <- length(dtp$x.limits)
  update(dtp, xlim = rep("", lxlims), grid = "h", 
         key = list(lines = list(lty = 2, lwd = 0.5),
                    text = list(c(paste("logLik displacement threshold:", 
                                        as.character(disp_thrh))), cex = 0.66), 
                    space = "top"))
}


## Cook's distance for beta estimates

CookDfun <- function(betaU, beta0, vb.inv){
  dbetaU <- betaU - beta0 # beta(-i) - beta
  CookD.value <- t(dbetaU) %*% vb.inv %*% dbetaU # Compute Cook distance
}

#' @title Cook's distance for the coefficient estimates.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`LMM_Model()`].
#' @param cook_thr Threshold for the Cook's distance.
#' @param label_angle Angle for the label of subjects with a Cook's distance greater than `cook_thr`.
#' @import nlme
#' @returns Plot of the Cook's distance value for each subject, indicating those subjects
#' whose Cook's distance is greater than `cook_thr`.
#' @export
Cooks_dist <- function(model, cook_thr, label_angle = 0){
  subject.c <- levels(model$data$Mouse)
  lmeUall <- lapply(subject.c, lmeU, model = model)
  names(lmeUall) <- subject.c
  
  betaUall <- sapply(lmeUall, fixef) # Matrix with betas(-i) estimates
  
  vcovb <- vcov(model) # Variance-covariance matrix for betas
  vb.inv <- solve(vcovb) # Inverse of of the var-cov matrix of betas
  
  beta0 <- fixef(model) # estimated betas
  
  CookD.num <- apply(betaUall, 2, CookDfun, beta0 = beta0, vb.inv = vb.inv)
  
  n.fixeff <- length(beta0) # Number of fixed effects
  rankX <- n.fixeff # Rank of matrix X
  CookD <- CookD.num/rankX # Cook's distance Di
  
  outD <- CookD > cook_thr # Outlying Di's
  
  print(paste("Subjects with Cook's distance greater than:", cook_thr))
  print(CookD[outD])
  
  subject.x <- 1:length(subject.c)
  plot(CookD, ylab = "Cook's Distance", type = "h", xlab = "Subject", main = "Cook's Distance vs Subject")
  points(subject.x[outD], CookD[outD], pch = 20)
  text(x = subject.x[outD], y = CookD[outD], labels = subject.c[outD], cex = 0.6, adj = c(0.5,0.1), srt = label_angle, xpd = TRUE)
  abline(h = cook_thr, lty = "dashed")
  legend("topright", inset = c(0,-0.1), legend = c(as.character(cook_thr)),lty = 2,
         title = "Cook's distance\n threshold", bty = "n", cex = 0.66,xpd = TRUE)
}



