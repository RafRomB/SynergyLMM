
## Synergy Calculation using LMM

#' @importFrom marginaleffects hypotheses
#' 
#' @title Synergy calculation using linear-mixed models.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param method String indicating the method for synergy calculation. Possible methods are "Bliss", "HSA" and "RA",
#' corresponding to Bliss, highest single agent and response additivity, respectively.
#' @param min_time Minimun time from which to start calculating synergy.
#' @param robustSE If TRUE, uncertainty is estimated using robust standard errors 
#' using a sandwich estimate of the variance-covariance matrix of the regression coefficient estimates provided by [clubSandwich::vcovCR()].
#' @param type Character string specifying which small-sample adjustment should be used, with available options "CR0", "CR1", "CR1p", "CR1S", "CR2", or "CR3". 
#' See "Details" section of [clubSandwich::vcovCR()] for further information.
#' @param ra_nsim Number of simulations to calculate the synergy for Response Additivity model.
#' @param norm_test If `method` is set to "RA", string indicating the test for checking the 
#' normality of the hypothesis expression of synergy (see "Warning#3" in "Description" section of [marginaleffects::hypotheses()]).
#' Possible values are "shapiroTest", "dagoTest" or "adTest", for Shapiro-Wilk, D'Agostino and Anderson-Darling normality tests, respectively.
#' @param t_ci Time for calculation of combination index value. The value of the CI represents the proportion of tumor cell survival at time `t_ci` 
#' in the drug combination group compared to the expected tumor cell survival according to the reference model. If not specified, `t_ci` is set to the
#' last time point of follow-up.
#' @param show_plot Logical indicating if a plot with the results of the synergy calculation should be generated.
#' @param ... Additional arguments to be passed to [marginaleffects::hypotheses()].
#' @returns The function returns a list with two elements:
#' - `Constrasts`: List with the outputs of the (non)-linear test for the synergy null hypothesis obtained by [marginaleffects::hypotheses()] for each time.
#' See [marginaleffects::hypotheses()] for more details.
#' - `Synergy`: Data frame with the synergy results, indicating the model of synergy ("Bliss", "HSA" or "RA"), the metric (combination index and synergy score),
#' the value of the metric estimate (with upper and lower confidence intervals) and the p-value, for each time.
#' @export

lmmSynergy <- function(model,
                    method = "Bliss",
                    min_time = 0,
                    robustSE = FALSE,
                    type = "CR2",
                    ra_nsim = 1000,
                    norm_test = "shapiroTest",
                    t_ci = NULL,
                    show_plot = TRUE,
                    ...) {
  
  # Validate method input
  valid_methods <- c("Bliss", "HSA", "RA")
  if (!method %in% valid_methods) {
    stop("Invalid 'method' provided. Choose from 'Bliss', 'HSA', or 'RA'.")
  }
  
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
    
    i <- 1
    for (d in times) {
      data <- model$data %>% dplyr::filter(.data$Time <= d)
      model_time <- update(model, data = data)
      
      # Define final time point for each model
      t2 <- d
      
      if (robustSE) {
        cluster_robust_vcov <- clubSandwich::vcovCR(model_time, type = type)
        summ <- clubSandwich::coef_test(model_time, vcov = cluster_robust_vcov)
        
        b1 <- rnorm(n = ra_nsim, mean = summ$beta[1], sd = summ$SE[1])
        b2 <- rnorm(n = ra_nsim, mean = summ$beta[2], sd = summ$SE[2])
        b3 <- rnorm(n = ra_nsim, mean = summ$beta[3], sd = summ$SE[3])
        b4 <- rnorm(n = ra_nsim, mean = summ$beta[4], sd = summ$SE[4])
        
        # Compute AUC for each simulation sample
        lhs_aucs <- sapply(b4, lhs_auc)
        rhs_aucs <- mapply(rhs_auc, beta_A = b2, beta_B = b3, beta_C = b1)
        
        # Compute the difference in AUCs
        delta_aucs <- lhs_aucs - rhs_aucs
        
        ratio_aucs <- lhs_aucs/rhs_aucs
        
        # Estimates
        
        delta <- median(delta_aucs)
        ratio <- median(ratio_aucs)
        median_rhs <- median(rhs_aucs)
        
        # 95% Confidence interval
        ci_delta <- quantile(delta_aucs, c(0.025, 0.975))
        
        ci_ratio <- quantile(ratio_aucs, c(0.025, 0.975))
        
        # p-value (two-tailed test)
        p_delta <-  2 * min(mean(delta_aucs <= 0), mean(delta_aucs >= 0))
        
        p_ratio <-  2 * min(mean(ratio_aucs <= 1), mean(ratio_aucs >= 1))

      } else {
        
        summ <- summary(model_time)
        
        b1 <- summ$tTable[1, c("Value", "Std.Error")]
        b2 <- summ$tTable[2, c("Value", "Std.Error")]
        b3 <- summ$tTable[3, c("Value", "Std.Error")]
        b4 <- summ$tTable[4, c("Value", "Std.Error")]
        
        b1 <- rnorm(n = 1000, mean = b1["Value"], sd = b1["Std.Error"])
        b2 <- rnorm(n = 1000, mean = b2["Value"], sd = b2["Std.Error"])
        b3 <- rnorm(n = 1000, mean = b3["Value"], sd = b3["Std.Error"])
        b4 <- rnorm(n = 1000, mean = b4["Value"], sd = b4["Std.Error"])
        
        # Compute AUC for each simulation sample
        lhs_aucs <- sapply(b4, lhs_auc)
        rhs_aucs <- mapply(rhs_auc, beta_A = b2, beta_B = b3, beta_C = b1)
        
        # Compute the difference in AUCs
        delta_aucs <- lhs_aucs - rhs_aucs
        
        ratio_aucs <- lhs_aucs/rhs_aucs
        
        # Estimates
        
        delta <- median(delta_aucs)
        ratio <- median(ratio_aucs)
        
        median_rhs <- median(rhs_aucs)
        
        
        # 95% Confidence interval
        ci_delta <- quantile(delta_aucs, c(0.025, 0.975))
        
        ci_ratio <- quantile(ratio_aucs, c(0.025, 0.975))
        
        # p-value (two-tailed test)
        p_delta <-  2 * min(mean(delta_aucs <= 0), mean(delta_aucs >= 0))
        
        p_ratio <-  2 * min(mean(ratio_aucs <= 1), mean(ratio_aucs >= 1))
        
      }
      
      ss <- rbind(ss,
                  data.frame(
                    method,
                    "SS",
                    -100 * delta/median_rhs,
                    -100 * ci_delta[2]/median_rhs,
                    -100 * ci_delta[1]/median_rhs,
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
      i <- i + 1
    }
  } else {
    
    if (method == "Bliss") {
      contrast <- "b4 = b2 + b3 - b1"
    }
    if (method == "HSA") {
      if (which.min(fixef_betas[2:3]) == 1) {
        contrast <- "b4 = b2"
      } else{
        contrast <- "b4 = b3"
      }
    }
    
    times <- unique(model$data$Time)
    times <- times[times >= min_time]
    times <- times[order(times)]
    
    if(is.null(t_ci)){
      t_ci <- max(times)
    }
    i <- 1
    for (d in times) {
      data <- model$data %>% dplyr::filter(.data$Time <= d)
      model_time <- update(model, data = data)
      if (robustSE) {
        Test <- hypotheses(
          model_time,
          hypothesis = contrast,
          vcov = clubSandwich::vcovCR(model, type = type),
          ...
        )
      } else {
        Test <- hypotheses(model_time, hypothesis = contrast, ...)
      }
      
      ss <- rbind(ss,
                  data.frame(
                    method,
                    "SS",
                    -100 * (Test$estimate)/(fixef_betas[2]+fixef_betas[3]-fixef_betas[1]),
                    -100 * (Test$conf.high)/(fixef_betas[2]+fixef_betas[3]-fixef_betas[1]),
                    -100 * (Test$conf.low)/(fixef_betas[2]+fixef_betas[3]-fixef_betas[1]),
                    Test$p.value,
                    d
                  ))
      
      ci <- rbind(ci, data.frame(
        method,
        "CI",
        exp(Test$estimate*t_ci),
        exp(Test$conf.low*t_ci),
        exp(Test$conf.high*t_ci),
        Test$p.value,
        d
      ))
      
      Contrasts[[i]] <- Test
      i <- i + 1
    }
    names(Contrasts) <- paste("Time", times, sep = "")
  }
  
  colnames(ss) <- c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Time")
  colnames(ci) <- c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Time")
  df <- rbind(ss, ci)
  rownames(df) <- NULL
  result <- list(Contrasts = Contrasts, Synergy = df)
  if(show_plot) {
    print(plot_lmmSynergy(result))
  }
  return(result)
}