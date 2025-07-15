#' @importFrom rlang .data
NULL

## Synergy Calculation using LMM

#' @importFrom marginaleffects hypotheses
#' @import stats
#' 
#' @title Synergy calculation using linear-mixed and non-linear mixed-effect models
#' @description
#' `lmmSynergy` allows for the calculation of synergy using 3 different references models: Bliss independence, highest single agent and
#' response additivity. The calculation of synergy is based on hypothesis testing on the coefficient estimates from the model fitted by
#' [`lmmModel()`].
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param method String indicating the method for synergy calculation. Possible methods are "Bliss", "HSA" and "RA",
#' corresponding to Bliss, highest single agent, and response additivity, respectively.
#' @param min_time Minimun time for which to start calculating synergy.
#' @param conf_level Numeric value between 0 and 1. Confidence level to use to build a confidence interval and obtain p-values. The default value is 0.95.
#' @param padj String indicating the correction method for adjusting p-values. Possible options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' More details in [stats::p.adjust()].
#' @param robust If TRUE, uncertainty is estimated using sandwich-based robust estimators
#' of the variance-covariance matrix of the regression coefficient estimates provided by [clubSandwich::vcovCR.lme]. Sandwich-based robust estimators
#' are only available for exponential growth models ('explme').
#' @param type Character string specifying which small-sample adjustment should be used if '`robust = TRUE`', with available options "CR0", "CR1", "CR1p", "CR1S", "CR2", or "CR3". 
#' See "Details" section of [clubSandwich::vcovCR()] for further information.
#' @param nsim Number of random sampling to calculate the synergy for Response Additivity model, and 
#' for synergy assessment in Gompertz growth models.
#' @param set_seed Logical indicating if the seed for those methods based on simulations (RA synergy and Gompertz growth models) should be fixed for reproducible results.
#' The seed can be also set before running `lmmSynergy()` using `set.seed()` function.
#' @param show_plot Logical indicating if a plot with the results of the synergy calculation should be generated.
#' @param ... Additional arguments to be passed to [marginaleffects::hypotheses()].
#' @details
#' `lmmSynergy` uses the statistical description provided by Demidenko and Miller (2019) for the calculation of synergy. It is based on hypothesis testing 
#' on the coefficients estimates from the model fitted by [`lmmModel()`].
#' 
#' ## Exponential Growth Model
#' 
#' Estimated coefficients \eqn{\hat{\beta}_C}, \eqn{\hat{\beta}_A}, \eqn{\hat{\beta}_B}, \eqn{\hat{\beta}_{AB}}, 
#' which represent the estimated specific growth rates for the Control, Drug A, Drug B and Combination groups, respectively, are to used to calculate synergy.
#' 
#' **Bliss Independence**
#' 
#' For Bliss model, `lmmSynergy` test the following null hypothesis:
#' 
#' _Two-drugs combination experiment_:
#' \deqn{H_0: \beta_{combination} = \beta_A + \beta_B - \beta_{control}}
#' 
#'  _Three-drugs combination experiment_:
#' \deqn{H_0: \beta_{combination} = \beta_A + \beta_B + \beta_C - 2\beta_{control}}
#' 
#' **Highes Single Agent (HSA)**
#' 
#' For the HSA model, `lmmSynergy` test the following null hypothesis:
#' 
#' _Two-drugs combination experiment_:
#' 
#' \deqn{H_0: \beta_{combination} = \min(\beta_A, \beta_B)}
#' 
#' _Three-drugs combination experiment_:
#' 
#' \deqn{H_0: \beta_{combination} = \min(\beta_A, \beta_B, \beta_C)}
#' 
#' **Response Additivity (RA)**
#' 
#' For the RA model, `lmmSynergy` test the following null hypothesis:
#' 
#' _Two-drugs combination experiment_:
#' \deqn{H_0: e^{\beta_{combination}t} = e^{\beta_At}+e^{\beta_Bt}-e^{\beta_{control}t}}
#' 
#' _Three-drugs combination experiment_:
#' \deqn{H_0: e^{\beta_{combination}t} = e^{\beta_At}+e^{\beta_Bt}+e^{\beta_Ct}-2e^{\beta_{control}t}}
#' 
#' For **Bliss** and **HSA** models, `lmmSynergy` uses [marginaleffects::hypotheses()] to conduct hypothesis tests on the estimated coefficients of the model.
#' 
#' In the case of the **RA** model, the null hypothesis is tested comparing the area under the curve (i.e. cumulative effect from the beginning of a treatment to 
#' a time point of interest) obtained from each side of the equation for the null hypothesis, based on `ra_sim` random samplings from the
#' distribution of the coefficients.
#' 
#' ## Gompertz Growth Model
#' 
#' Estimated coefficients \eqn{r_{0_{T_i}}} and \eqn{\rho_{T_i}} are to used to calculate synergy.
#' - \eqn{r_{0_{T_i}}} is the fixed effect for the initial growth rate for treatment group \eqn{T_i}.
#' - \eqn{\rho_{T_i}} is the fixed effect for the constant accounting for the reduction in the tumor growth rate for treatment group \eqn{T_i}.
#' 
#' The following expressions are simplified with these symbols:
#' 
#' \deqn{\gamma = \frac{r_{0_{Control}}}{\rho_{Control}} \cdot (1-e^{-\rho_{Control}\cdot t})}
#' \deqn{A = \frac{r_{0_{DrugA}}}{\rho_{DrugA}} \cdot (1-e^{-\rho_{DrugA}\cdot t})}
#' \deqn{B = \frac{r_{0_{DrugB}}}{\rho_{DrugB}} \cdot (1-e^{-\rho_{DrugB}\cdot t})}
#' \deqn{C = \frac{r_{0_{DrugC}}}{\rho_{DrugC}} \cdot (1-e^{-\rho_{DrugC}\cdot t})}
#' \deqn{\Lambda = \frac{r_{0_{Combination}}}{\rho_{Combination}} \cdot (1-e^{-\rho_{Combination}\cdot t})}
#' 
#' **Bliss Independence**
#' 
#' For Bliss model, `lmmSynergy` test the following null hypothesis:
#' 
#' _Two-drugs combination experiment_:
#' \deqn{H_0: \Lambda = A + B - \gamma}
#' 
#'  _Three-drugs combination experiment_:
#' \deqn{H_0: \Lambda = A + B + C- 2\gamma}
#' 
#' **Highes Single Agent (HSA)**
#' 
#' For the HSA model, `lmmSynergy` test the following null hypothesis:
#' 
#' _Two-drugs combination experiment_:
#' 
#' \deqn{H_0: \Lambda = \min(A,B)}
#' 
#' _Three-drugs combination experiment_:

#' \deqn{H_0: \Lambda = \min(A,B,C)}
#' 
#' **Response Additivity (RA)**
#' 
#' For the RA model, `lmmSynergy` test the following null hypothesis:
#' 
#' _Two-drugs combination experiment_:
#' \deqn{H_0: e^{\Lambda} = e^{A}+e^{B}-e^{\gamma}}
#' 
#' _Three-drugs combination experiment_:
#' \deqn{H_0: e^{\Lambda} = e^{A}+e^{B}+e^{C}-2e^{\gamma}}
#' 
#' For the **Gompertz models**, the null hypothesis is tested comparing the area under the curve (i.e. cumulative effect from the beginning of a treatment to 
#' a time point of interest) obtained from each side of the equations for the null hypothesis, based on `nsim` random samplings from the
#' distribution of the coefficients.
#' 
#' ## Combination Index and Synergy Score
#' 
#' The results obtained by `lmmSynergy` include the synergy score (SS) and combination index (CI) for the model, for each time point, together with their confidence interval,
#' and the corresponding p-value. The values of SS and CI provided by `lmmSynergy` follow previous definitions of these metrics so they have the same interpretation:
#' 
#' - The SS has been defined as the excess response due to drug interaction compared to the reference model (Ianevski et al. (2017), Ianevski, Giri, and Aittokallio (2022), Mao and Guo (2023)).
#' Following this definition, a \eqn{SS>0}, \eqn{SS=0}, and \eqn{SS<0}, represent synergistic, additive and antagonistic effects, respectively.
#' 
#' - According to the common definition of the CI, a \eqn{CI<1}, \eqn{CI=1}, and \eqn{CI>1} represent synergistic, additive and antagonistic effects, respectively (Yadav et al. (2015), Demidenko and Miller (2019), 
#' Mao and Guo (2023)), and provides information about the observed drug combination effect versus the expected additive effect provided by the reference synergy model. 
#' A drug combination effect larger than the expected (\eqn{CI < 1}) would indicate synergism, a drug combination effect equal to the expected (\eqn{CI = 1}) would indicate additivity, 
#' and a lower drug combination effect than the expected (\eqn{CI > 1}) would indicate antagonism.
#' 
#' As mentioned above, the results include the synergy results for **each day**. This means that `lmmSynergy` refits the model using the data from `time_start` defined in [lmmModel()] until 
#' each time point, providing the synergy results for each of these models and for that specific time point. 
#' 
#' **Uncertainty estimation using robust estimators**
#' 
#' If `robust = TRUE`, `lmmSynergy` deals with possible model misspecifications, allowing for cluster-robust variance estimation using  [clubSandwich::vcovCR.lme].
#' When using `robust = TRUE`, setting `type = "CR2"` is recommended. See more details in [clubSandwich::vcovCR()]. This option is only available for 'explme' exponential tumor growth models.
#' 
#' _Note_: When a variance structure has been specified in the model it is recommended to use always `robust = TRUE` to get a better estimation. 
#' 
#' @returns The function returns a list with two elements:
#' - `Constrasts`: List with the outputs of the linear test for the synergy null hypothesis obtained by [marginaleffects::hypotheses()] for each time.
#' See [marginaleffects::hypotheses()] for more details.
#' - `Synergy`: Data frame with the synergy results, indicating the model of synergy ("Bliss", "HSA" or "RA"), the metric (combination index and synergy score),
#' the value of the metric estimate (with upper and lower confidence interval bounds) and the p-value, for each time.
#' - `Estimates`: Data frame with the estimates from each model at each time point, obtained with [`lmmModel_estimates()`] function.
#' For 'explme' exponential tumor growth models, if `robust=TRUE`, sandwich-based robust estimators for the standard errors of the estimated coefficients are reported.
#' 
#' If `show_plot = TRUE`, a plot with the synergy results obtained with [plot_lmmSynergy()] is also shown.
#' @references 
#' - Demidenko, Eugene, and Todd W. Miller. 2019. _Statistical Determination of Synergy Based on Bliss Definition of Drugs Independence._ PLoS ONE 14 (November). https://doi.org/10.1371/journal.pone.0224137.
#' - Yadav, Bhagwan, Krister Wennerberg, Tero Aittokallio, and Jing Tang. 2015. _Searching for Drug Synergy in Complex Dose–Response Landscapes Using an Interaction Potency Model._ Computational and Structural Biotechnology Journal 13: 504–13. https://doi.org/10.1016/j.csbj.2015.09.001.
#' - Ianevski, Aleksandr, Liye He, Tero Aittokallio, and Jing Tang. 2017. _SynergyFinder: A Web Application for Analyzing Drug Combination Dose–Response Matrix Data._ Bioinformatics 33 (August): 2413–15. https://doi.org/10.1093/bioinformatics/btx162.
#' - Ianevski, Aleksandr, Anil K Giri, and Tero Aittokallio. 2022. _SynergyFinder 3.0: An Interactive Analysis and Consensus Interpretation of Multi-Drug Synergies Across Multiple Samples._ Nucleic Acids Research 50 (July): W739–43. https://doi.org/10.1093/nar/gkac382.
#' - Mao, Binchen, and Sheng Guo. 2023. _Statistical Assessment of Drug Synergy from in Vivo Combination Studies Using Mouse Tumor Models._ Cancer Research Communications 3 (October): 2146–57. https://doi.org/10.1158/2767-9764.CRC-23-0243.
#' - Vincent Arel-Bundock, Noah Greifer, and Andrew Heiss. Forthcoming. How to Interpret Statistical Models Using marginaleffects in R and Python. Journal of Statistical Software. https://marginaleffects.com
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
#' # Most simple use with default values
#' syn <- lmmSynergy(lmm)
#' # Accessing to synergy results data frame
#' syn$Synergy
#' # Selecting different reference models:
#' ## Bliss
#' lmmSynergy(lmm, method = "Bliss")
#' ## HSA
#' lmmSynergy(lmm, method = "HSA")
#' ## RA
#' lmmSynergy(lmm, method = "RA", ra_sim = 1000)
#' 
#' # Only calculate synergy from Time 12 onwards
#' lmmSynergy(lmm, min_time = 12)
#' 
#' # Using robust standard errors
#' lmmSynergy(lmm, method = "Bliss", robust = TRUE, type = "CR2")
#' 

#' @export

lmmSynergy <- function(model,
                       method = "Bliss",
                       min_time = 0,
                       conf_level = 0.95,
                       padj = "none",
                       robust = FALSE,
                       type = "CR2",
                       nsim = 1000,
                       set_seed = TRUE,
                       show_plot = TRUE,
                       ...) {
  UseMethod("lmmSynergy")
}


## lmmSynergy.explme method
#' This is the method for \code{lmmSynergy()} generic function to calculate synergy
#' using the exponential tumor growth model.
#' @rdname lmmSynergy
#' @method lmmSynergy explme
#' @export
lmmSynergy.explme <- function(model,
                              method = "Bliss",
                              min_time = 0,
                              conf_level = 0.95,
                              padj = "none",
                              robust = FALSE,
                              type = "CR2",
                              nsim = 1000,
                              set_seed = TRUE,
                              show_plot = TRUE,
                              ...) {
  
  # Validate method input
  valid_methods <- c("Bliss", "HSA", "RA")
  if (!method %in% valid_methods) {
    stop("Invalid 'method' provided. Choose from 'Bliss', 'HSA', or 'RA'.")
  }
  
  estimates <- data.frame()
  ss <- data.frame()
  ci <- data.frame()
  Contrasts <- list()
  
  fixef_betas <- nlme::fixef(model)
  
  if(method == "RA") {
    
    # Function to compute LHS AUC
    lhs_auc <- function(beta_AB) {
      integrate(function(t)
        exp(beta_AB * t), t1, t2)$value
    }
    
    # Function to compute RHS AUC
    # 3 drugs
    rhs_auc3 <- function(beta_A, beta_B, beta_Z, beta_C) {
      integrate(function(t) {
        (exp(beta_A * t) + exp(beta_B * t) + exp(beta_Z * t) - 2*exp(beta_C * t))
      }, t1, t2)$value
    }
    # 2 drugs
    rhs_auc <- function(beta_A, beta_B, beta_C) {
      integrate(function(t) {
        (exp(beta_A * t) + exp(beta_B * t) - exp(beta_C * t))
      }, t1, t2)$value
    }
    
    # Define initial time point for RA calculation
    t1 <- min(model$dt1$Time)
    
    # Times to calculate synergy
    times <- unique(model$data$Time)
    times <- times[times >= min_time]
    times <- times[order(times)]
    
    Contrasts <- NULL
    
    for (d in times) {
      data <- model$data %>% dplyr::filter(.data$Time <= d)
      model_time <- update(model, data = data)
      model_time_estimates <- lmmModel_estimates.explme(model_time, robust = robust, type = type)
      model_time_estimates$Time <- d
      rownames(model_time_estimates) <- paste0("estimates_Time_",d)
      estimates <- rbind(estimates, model_time_estimates)
      
      # Define final time point for each model
      t2 <- d
      
      if (robust) {
        cluster_robust_vcov <- clubSandwich::vcovCR(model_time, type = type) # Cluster-robust variance-covariance matrix
        betas <- fixef(model_time) # model betas estimates
        
        if (set_seed) {
          set.seed(123)
        }
        
        betas_mvnorm <- MASS::mvrnorm(n = nsim, mu = betas, Sigma = cluster_robust_vcov) # Simulate from the multivariate normal distribution
        
        b1 <- betas_mvnorm[,1]
        b2 <- betas_mvnorm[,2]
        b3 <- betas_mvnorm[,3]
        if (length(fixef_betas) == 5){
          b4 <- betas_mvnorm[,4]
          b5 <- betas_mvnorm[,5]
          # Compute AUC for each simulation sample
          lhs_aucs <- sapply(b5, lhs_auc)
          rhs_aucs <- mapply(rhs_auc3, beta_A = b2, beta_B = b3, beta_Z = b4, beta_C = b1)
        } else {
          b4 <- betas_mvnorm[,4]
          # Compute AUC for each simulation sample
          lhs_aucs <- sapply(b4, lhs_auc)
          rhs_aucs <- mapply(rhs_auc, beta_A = b2, beta_B = b3, beta_C = b1)
        }
        
        # Substitute negative values by 0
        rhs_aucs[which(rhs_aucs<0)] <- 0
        
        # Compute the difference in AUCs
        delta_aucs <- lhs_aucs - rhs_aucs
        
        ratio_aucs <- lhs_aucs/rhs_aucs
        
        if (sum(ratio_aucs == Inf) > 0) {
          message(paste0("Some Combination Index values were infinite (+Inf) due to zero denominators for time: ", d, ".
                   These values have been capped at a maximum of 100 to preserve interpretability."))
          ratio_aucs[which(ratio_aucs == Inf)] <- 100
        }
        
        # Estimates
        
        delta <- median(delta_aucs)
        ratio <- median(ratio_aucs)
        sd_delta <- sd(delta_aucs)
        
        # 95% Confidence interval
        ci_delta <- quantile(delta_aucs, c((1-conf_level)/2, 1-(1-conf_level)/2))
        
        ci_ratio <- quantile(ratio_aucs, c((1-conf_level)/2, 1-(1-conf_level)/2))
        
        # p-value (two-tailed test)
        p_delta <-  2 * min(mean(delta_aucs <= 0), mean(delta_aucs >= 0))
        
        p_ratio <-  2 * min(mean(ratio_aucs <= 1), mean(ratio_aucs >= 1))

      } else {
        
        vcov_mtx <- vcov(model_time) # Variance-Covariance Matrix for the Fitted Model Object
        betas <- fixef(model_time) # model betas estimates
        
        if (set_seed) {
          set.seed(123)
        }
        
        betas_mvnorm <- MASS::mvrnorm(n = nsim, mu = betas, Sigma = vcov_mtx) # Simulate from the multivariate normal distribution
        
        b1 <- betas_mvnorm[,1]
        b2 <- betas_mvnorm[,2]
        b3 <- betas_mvnorm[,3]
        if (length(fixef_betas) == 5){
          b4 <- betas_mvnorm[,4]
          b5 <- betas_mvnorm[,5]
          # Compute AUC for each simulation sample
          lhs_aucs <- sapply(b5, lhs_auc)
          rhs_aucs <- mapply(rhs_auc3, beta_A = b2, beta_B = b3, beta_Z = b4, beta_C = b1)
        } else {
          b4 <- betas_mvnorm[,4]
          # Compute AUC for each simulation sample
          lhs_aucs <- sapply(b4, lhs_auc)
          rhs_aucs <- mapply(rhs_auc, beta_A = b2, beta_B = b3, beta_C = b1)
        }
        
        # Substitute negative values by 0
        rhs_aucs[which(rhs_aucs<0)] <- 0
        
        # Compute the difference in AUCs
        delta_aucs <- lhs_aucs - rhs_aucs
        
        ratio_aucs <- lhs_aucs/rhs_aucs
        
        if (sum(ratio_aucs == Inf) > 0) {
          message(paste0("Some Combination Index values were infinite (+Inf) due to zero denominators for time: ", d, ".
                   These values have been capped at a maximum of 100 to preserve interpretability."))
          ratio_aucs[which(ratio_aucs == Inf)] <- 100
        }
        
        # Estimates
        
        delta <- median(delta_aucs)
        ratio <- median(ratio_aucs)
        
        sd_delta <- sd(delta_aucs)
        
        # 95% Confidence interval
        ci_delta <- quantile(delta_aucs, c((1-conf_level)/2, 1-(1-conf_level)/2))
        
        ci_ratio <- quantile(ratio_aucs, c((1-conf_level)/2, 1-(1-conf_level)/2))
        
        # p-value (two-tailed test)
        p_delta <-  2 * min(mean(delta_aucs <= 0), mean(delta_aucs >= 0))
        
        p_ratio <-  2 * min(mean(ratio_aucs <= 1), mean(ratio_aucs >= 1))
        
      }
      
      ss <- rbind(ss,
                  data.frame(
                    method,
                    "SS",
                    - delta/sd_delta,
                    - ci_delta[2]/sd_delta,
                    - ci_delta[1]/sd_delta,
                    p_delta,
                    d
                  ))
      
      ci <- rbind(ci, data.frame(
        method,
        "CI",
        ratio,
        ci_ratio[1],
        ci_ratio[2],
        p_ratio,
        d
      ))
    }
  } else {
    
    if (method == "Bliss") {
      if (length (fixef_betas) == 5) {
        contrast <- "b5 = b2 + b3 + b4 - 2*b1"
      } else {
        contrast <- "b4 = b2 + b3 - b1"
      }
    }
    if (method == "HSA") {
      if (length(fixef_betas) == 5) {
        if (which.min(fixef_betas[2:4]) == 1) {
          contrast <- "b5 = b2"
        } else if (which.min(fixef_betas[2:4]) == 2){
          contrast <- "b5 = b3"
        } else {
          contrast <- "b5 = b4"
        }
      } else {
        if (which.min(fixef_betas[2:3]) == 1) {
          contrast <- "b4 = b2"
        } else{
          contrast <- "b4 = b3"
        } 
      }
    }
    
    times <- unique(model$data$Time)
    times <- times[times >= min_time]
    times <- times[order(times)]
    
    i <- 1
    for (d in times) {
      data <- model$data %>% dplyr::filter(.data$Time <= d)
      model_time <- update(model, data = data)
      model_time_estimates <- lmmModel_estimates.explme(model_time, robust = robust, type = type)
      model_time_estimates$Time <- d
      estimates <- rbind(estimates, model_time_estimates)
      
      if (robust) {
        Test <- hypotheses(
          model_time,
          hypothesis = contrast,
          vcov = clubSandwich::vcovCR(model_time, type = type),
          conf_level = conf_level,
          ...
        )
      } else {
        Test <- hypotheses(model_time, hypothesis = contrast, conf_level = conf_level, ...)
      }
      
      ss <- rbind(ss,
                  data.frame(
                    method,
                    "SS",
                    - (Test$estimate)/(Test$std.error),
                    - (Test$conf.high)/(Test$std.error),
                    - (Test$conf.low)/(Test$std.error),
                    Test$p.value,
                    d
                  ))
      
      ci <- rbind(ci, data.frame(
        method,
        "CI",
        exp(Test$estimate*d),
        exp(Test$conf.low*d),
        exp(Test$conf.high*d),
        Test$p.value,
        d
      ))
      
      Contrasts[[i]] <- Test
      i <- i + 1
    }
    names(Contrasts) <- paste("Time", times, sep = "")
  }
  
  colnames(ci) <- c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Time")
  if (padj != "none") {
    ci$padj <- stats::p.adjust(p = ci$pval, method = padj)
    ci <- ci[,c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj", "Time")]
  }
  
  
  colnames(ss) <- c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Time")
  if (padj != "none") {
    ss$padj <- stats::p.adjust(p = ss$pval, method = padj)
    ss <- ss[,c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj", "Time")]
  }
  
  df <- rbind(ci, ss)
  
  rownames(df) <- NULL
  if (sum(df$pval == 0) > 0) {
    ndec <- strsplit(format(nsim, scientific = T), split = "\\+")[[1]][2]
    apx_p <- paste("p<1e-",ndec, sep = "")
    warning(paste("p-values below", apx_p, "are approximated to 0."),
            " If you used method = 'RA' consider increasing 'nsim' value for",
            " more precise p-values.")
  }
  result <- list(Contrasts = Contrasts, Synergy = df, Estimates = estimates, nsim = nsim)
  if(show_plot) {
    plot(plot_lmmSynergy(result)$CI_SS)
  }
  attr(result, "SynergyLMM") <- "lmmSynergy"
  return(result)
}


## lmmSynergy.gompertzlme method
#' This is the method for \code{lmmSynergy()} generic function to calculate synergy
#' using the Gompertz tumor growth model.
#' @rdname lmmSynergy
#' @method lmmSynergy gompertzlme
#' @export
lmmSynergy.gompertzlme <- function(model,
                                   method = "Bliss",
                                   min_time = 0,
                                   conf_level = 0.95,
                                   padj = "none",
                                   robust = FALSE,
                                   type = "CR2",
                                   nsim = 10000,
                                   set_seed = TRUE,
                                   show_plot = TRUE,
                              ...) {
  
  # Validate method input
  valid_methods <- c("Bliss", "HSA", "RA")
  if (!method %in% valid_methods) {
    stop("Invalid 'method' provided. Choose from 'Bliss', 'HSA', or 'RA'.")
  }
  
  if(robust) {
    warning("Sandwich-based robust estimators are only available for exponential growth models, 
    and the results are from standard estimation techniques.")
  }
  
  ss <- data.frame()
  ci <- data.frame()

  fixef_betas <- nlme::fixef(model)
  
  # Define initial time point for calculation
  t1 <- min(model$dt1$Time)
  
  # Times to calculate synergy
  times <- unique(model$data$Time)
  times <- times[times >= min_time]
  times <- times[order(times)]
  
  data <- model$data
  model_time <- model
  estimates <- data.frame(lmmModel_estimates(model_time), row.names = NULL)
  
  vcov_mtx <- vcov(model_time) # Variance-Covariance Matrix for the Fitted Model Object
  betas <- fixef(model_time) # model betas estimates
  
  if (set_seed) {
    set.seed(123)
  }
  
  betas_mvnorm <- MASS::mvrnorm(n = nsim, mu = betas, Sigma = vcov_mtx) # Simulate from the multivariate normal distribution
  
  if (length(unique(data$Treatment)) == 4) { # Check number of treatments: for 2-drug combination
    r0Ctrl <- betas_mvnorm[,1]
    r0A <- betas_mvnorm[,2]
    r0B <- betas_mvnorm[,3]
    r0AB <- betas_mvnorm[,4]
    rhoCtrl <- betas_mvnorm[,5]
    rhoA <- betas_mvnorm[,6]
    rhoB <- betas_mvnorm[,7]
    rhoAB <- betas_mvnorm[,8]
    
  } else if (length(unique(data$Treatment)) == 5) { # Check number of treatments: for 3-drug combination
    r0Ctrl <- betas_mvnorm[,1]
    r0A <- betas_mvnorm[,2]
    r0B <- betas_mvnorm[,3]
    r0Z <- betas_mvnorm[,4]
    r0AB <- betas_mvnorm[,5]
    rhoCtrl <- betas_mvnorm[,6]
    rhoA <- betas_mvnorm[,7]
    rhoB <- betas_mvnorm[,8]
    rhoZ <- betas_mvnorm[,9]
    rhoAB <- betas_mvnorm[,10]
    
  } else {
    stop("Incorrect number of treatments. Synergy analysis can only be performed for 2 or 3 drug combinations.")
  }
  
  
  # LHS AUC computation
  
  lhs_auc <- function(r0_AB, rho_AB) { # r0_AB; rho_AB: coefficients for combination group
    integrate(function(t) {
      (r0_AB/rho_AB)*(1-exp(-rho_AB*t))
    }, t1, t2)$value
  }
  
  # RHS AUC computation using Bliss
  
  # Function to compute RHS AUC
  # 3 drugs
  rhs_auc3_bliss <- function(r0_A, rho_A, r0_B, rho_B, r0_Z, rho_Z, r0_C, rho_C) { # A, B, Z, C: coefficients for first, second, third, and control groups, respectively
    integrate(function(t) {
      (r0_A/rho_A)*(1-exp(-rho_A*t)) + (r0_B/rho_B)*(1-exp(-rho_B*t)) + (r0_Z/rho_Z)*(1-exp(-rho_Z*t)) - 2*((r0_C/rho_C)*(1-exp(-rho_C*t)))
    }, t1, t2)$value
  }
  # 2 drugs
  rhs_auc_bliss <- function(r0_A, rho_A, r0_B, rho_B, r0_C, rho_C) { # A, B, C: coefficients for first, second, and control groups, respectively
    integrate(function(t) {
      (r0_A/rho_A)*(1-exp(-rho_A*t)) + (r0_B/rho_B)*(1-exp(-rho_B*t)) - (r0_C/rho_C)*(1-exp(-rho_C*t))
    }, t1, t2)$value
  }
  
  # RHS AUC computation using HSA
  
  rhs_auc_hsa <- function(r0, rho) {
    integrate(function(t) {
      (r0/rho) * (1 - exp(-rho * t))
    }, t1, t2)$value
  }
  
  # RHS AUC computation using RA
  
  # 3 drugs
  rhs_auc3_ra <- function(r0_A, rho_A, r0_B, rho_B, r0_Z, rho_Z, r0_C, rho_C) { # A, B, Z, C: coefficients for first, second, third, and control groups, respectively
    integrate(function(t) {
      exp((r0_A/rho_A)*(1-exp(-rho_A*t))) + exp((r0_B/rho_B)*(1-exp(-rho_B*t))) + exp((r0_Z/rho_Z)*(1-exp(-rho_Z*t))) - 2*exp(((r0_C/rho_C)*(1-exp(-rho_C*t))))
    }, t1, t2)$value
  }
  # 2 drugs
  rhs_auc_ra <- function(r0_A, rho_A, r0_B, rho_B, r0_C, rho_C) { # A, B, C: coefficients for first, second, and control groups, respectively
    integrate(function(t) {
      exp((r0_A/rho_A)*(1-exp(-rho_A*t))) + exp((r0_B/rho_B)*(1-exp(-rho_B*t))) - exp((r0_C/rho_C)*(1-exp(-rho_C*t)))
    }, t1, t2)$value
  }
  
  
  if(method == "Bliss") {
    
    for (d in times) {
      
      # Define final time point for each model
      t2 <- d
      
      lhs_aucs <- mapply(lhs_auc, r0_AB = r0AB, rho_AB = rhoAB)
      
      if (length(unique(data$Treatment)) == 4) { # Check number of treatments: for 2-drug combination
        
        rhs_aucs <- mapply(rhs_auc_bliss, r0_A = r0A, rho_A = rhoA, r0_B = r0B, rho_B = rhoB, r0_C = r0Ctrl, rho_C = rhoCtrl)
        
      } else if (length(unique(data$Treatment)) == 5) { # Check number of treatments: for 3-drug combination
        
        rhs_aucs <- mapply(rhs_auc3_bliss, r0_A = r0A, rho_A = rhoA, r0_B = r0B, rho_B = rhoB, r0_Z = r0Z, rho_Z = rhoZ, r0_C = r0Ctrl, rho_C = rhoCtrl)
      } else {
        stop("Incorrect number of treatments. Synergy analysis can only be performed for 2 or 3 drug combinations.")
      }
      
      # Substitute negative values by 0
      rhs_aucs[which(rhs_aucs<0)] <- 0
      
      # Compute the difference in AUCs
      delta_aucs <- lhs_aucs - rhs_aucs
      
      ratio_aucs <- exp(delta_aucs)
      
      # Estimates
      
      delta <- median(delta_aucs)
      ratio <- median(ratio_aucs)
      
      sd_delta <- sd(delta_aucs)
      
      # 95% Confidence interval
      ci_delta <- quantile(delta_aucs, c((1-conf_level)/2, 1-(1-conf_level)/2))
      
      ci_ratio <- quantile(ratio_aucs, c((1-conf_level)/2, 1-(1-conf_level)/2))
      
      # p-value (two-tailed test)
      p_delta <-  2 * min(mean(delta_aucs <= 0), mean(delta_aucs >= 0))
      
      p_ratio <-  2 * min(mean(ratio_aucs <= 1), mean(ratio_aucs >= 1))
      
      ss <- rbind(ss,
                  data.frame(
                    method,
                    "SS",
                    - delta/sd_delta,
                    - ci_delta[2]/sd_delta,
                    - ci_delta[1]/sd_delta,
                    p_delta,
                    d
                  ))
      
      ci <- rbind(ci, data.frame(
        method,
        "CI",
        ratio,
        ci_ratio[1],
        ci_ratio[2],
        p_ratio,
        d
      ))
    }
  } else if (method == "HSA") {
    
    for (d in times) {
      
      # Define final time point for each model
      t2 <- d
      
      lhs_aucs <- mapply(lhs_auc, r0_AB = r0AB, rho_AB = rhoAB)
      rhs_aucsA <- mapply(rhs_auc_hsa, r0 = r0A, rho = rhoA)
      rhs_aucsB <- mapply(rhs_auc_hsa, r0 = r0B, rho = rhoB)
      
      if (length(unique(data$Treatment)) == 4) { # Check number of treatments: for 2-drug combination
        
        rhs_aucs <- list(rhs_aucsA, rhs_aucsB)
        
        delta_aucs_AB <- rhs_aucsA - rhs_aucsB
        
        if (2 * min(mean(delta_aucs_AB <= 0), mean(delta_aucs_AB >= 0)) > 0.05){
          drugs <- c("'drug_a'", "'drug_b'")
          warning(paste("No significant difference between single drugs effects.",
                        "Treatment assigned to", drugs[which.min(c(median(rhs_aucsA), median(rhs_aucsB)))], 
                        "is used for the HSA effect comparison."))
        }
        
        rhs_aucs <- rhs_aucs[[which.min(c(median(rhs_aucsA), median(rhs_aucsB)))]]
        
      } else if (length(unique(data$Treatment)) == 5) { # Check number of treatments: for 3-drug combination
        
        rhs_aucsZ <- mapply(rhs_auc_hsa, r0 = r0Z, rho = rhoZ)
        rhs_aucs <- list(rhs_aucsA, rhs_aucsB, rhs_aucsZ)
        
        delta_aucs_AB <- rhs_aucsA - rhs_aucsB
        delta_aucs_AZ <- rhs_aucsA - rhs_aucsZ
        delta_aucs_BZ <- rhs_aucsB - rhs_aucsZ
        
        if (2 * min(mean(delta_aucs_AB <= 0), mean(delta_aucs_AB >= 0)) > 0.05 &
            2 * min(mean(delta_aucs_AZ <= 0), mean(delta_aucs_AZ >= 0)) > 0.05 &
            2 * min(mean(delta_aucs_BZ <= 0), mean(delta_aucs_BZ >= 0)) > 0.05) {
          drugs <- c("'drug_a'", "'drug_b'", "'drug_c'")
          warning(paste("No significant difference between single drugs effects.",
                        "Treatment assigned to", drugs[which.min(c(median(rhs_aucsA), median(rhs_aucsB), median(rhs_aucsZ)))], 
                        "is used for the HSA effect comparison."))
        }
        rhs_aucs <- rhs_aucs[[which.min(c(median(rhs_aucsA), median(rhs_aucsB), median(rhs_aucsZ)))]]
      } else {
        stop("Incorrect number of treatments. Synergy analysis can only be performed for 2 or 3 drug combinations.")
      }
      
      # Substitute negative values by 0
      rhs_aucs[which(rhs_aucs<0)] <- 0
      
      # Compute the difference in AUCs
      delta_aucs <- lhs_aucs - rhs_aucs
      
      ratio_aucs <- exp(delta_aucs)
      
      # Estimates
      
      delta <- median(delta_aucs)
      ratio <- median(ratio_aucs)
      
      sd_delta <- sd(delta_aucs)
      
      # 95% Confidence interval
      ci_delta <- quantile(delta_aucs, c((1-conf_level)/2, 1-(1-conf_level)/2))
      
      ci_ratio <- quantile(ratio_aucs, c((1-conf_level)/2, 1-(1-conf_level)/2))
      
      # p-value (two-tailed test)
      p_delta <-  2 * min(mean(delta_aucs <= 0), mean(delta_aucs >= 0))
      
      p_ratio <-  2 * min(mean(ratio_aucs <= 1), mean(ratio_aucs >= 1))
      
      ss <- rbind(ss,
                  data.frame(
                    method,
                    "SS",
                    - delta/sd_delta,
                    - ci_delta[2]/sd_delta,
                    - ci_delta[1]/sd_delta,
                    p_delta,
                    d
                  ))
      
      ci <- rbind(ci, data.frame(
        method,
        "CI",
        ratio,
        ci_ratio[1],
        ci_ratio[2],
        p_ratio,
        d
      ))
    }
  } else if (method == "RA") {
    
    for (d in times) {
      
      # Define final time point for each model
      t2 <- d
      
      lhs_aucs <- mapply(lhs_auc, r0_AB = r0AB, rho_AB = rhoAB)
      
      if (length(unique(data$Treatment)) == 4) { # Check number of treatments: for 2-drug combination
        rhs_aucs <- mapply(rhs_auc_ra, r0_A = r0A, rho_A = rhoA, r0_B = r0B, rho_B = rhoB, r0_C = r0Ctrl, rho_C = rhoCtrl)
        
      } else if (length(unique(data$Treatment)) == 5) { # Check number of treatments: for 3-drug combination
        rhs_aucs <- mapply(rhs_auc3_ra, r0_A = r0A, rho_A = rhoA, r0_B = r0B, rho_B = rhoB, r0_Z = r0Z, rho_Z = rhoZ, r0_C = r0Ctrl, rho_C = rhoCtrl)
        
      } else {
        stop("Incorrect number of treatments. Synergy analysis can only be performed for 2 or 3 drug combinations.")
      }
      
      # Substitute negative values by 0
      rhs_aucs[which(rhs_aucs<0)] <- 0
      
      # Compute the difference in AUCs
      delta_aucs <- lhs_aucs - rhs_aucs
      
      ratio_aucs <- lhs_aucs/rhs_aucs
      
      if (sum(ratio_aucs == Inf) > 0) {
        message(paste0("Some Combination Index values were infinite (+Inf) due to zero denominators for time: ", d, ".
                   These values have been capped at a maximum of 100 to preserve interpretability."))
        ratio_aucs[which(ratio_aucs == Inf)] <- 100
      }
      
      # Estimates
      
      delta <- median(delta_aucs)
      ratio <- median(ratio_aucs)
      
      sd_delta <- sd(delta_aucs)
      
      # 95% Confidence interval
      ci_delta <- quantile(delta_aucs, c((1-conf_level)/2, 1-(1-conf_level)/2))
      
      ci_ratio <- quantile(ratio_aucs, c((1-conf_level)/2, 1-(1-conf_level)/2))
      
      # p-value (two-tailed test)
      p_delta <-  2 * min(mean(delta_aucs <= 0), mean(delta_aucs >= 0))
      
      p_ratio <-  2 * min(mean(ratio_aucs <= 1), mean(ratio_aucs >= 1))
      
      ss <- rbind(ss,
                  data.frame(
                    method,
                    "SS",
                    - delta/sd_delta,
                    - ci_delta[2]/sd_delta,
                    - ci_delta[1]/sd_delta,
                    p_delta,
                    d
                  ))
      
      ci <- rbind(ci, data.frame(
        method,
        "CI",
        ratio,
        ci_ratio[1],
        ci_ratio[2],
        p_ratio,
        d
      ))
    }
  }
  
  colnames(ci) <- c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Time")
  if (padj != "none") {
    ci$padj <- stats::p.adjust(p = ci$pval, method = padj)
    ci <- ci[,c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj", "Time")]
  }
  
  
  colnames(ss) <- c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Time")
  if (padj != "none") {
    ss$padj <- stats::p.adjust(p = ss$pval, method = padj)
    ss <- ss[,c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj", "Time")]
  }
  
  df <- rbind(ci, ss)
  rownames(df) <- NULL
  if (sum(df$pval == 0) > 0) {
    ndec <- strsplit(format(nsim, scientific = T), split = "\\+")[[1]][2]
    apx_p <- paste("p<1e-",ndec, sep = "")
    warning(paste("p-values below", apx_p, "are approximated to 0."),
            " If you used a Gompertz model, consider increasing 'nsim' value for",
            " more precise p-values.")
  }
  result <- list(Synergy = df, Estimates = estimates, nsim = nsim)
  if(show_plot) {
    plot(plot_lmmSynergy(result)$CI_SS)
  }
  attr(result, "SynergyLMM") <- "lmmSynergy"
  return(result)
}
