
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
        cluster_robust_vcov <- clubSandwich::vcovCR(model_time, type = type) # Cluster-robust variance-covariance matrix
        betas <- fixef(model_time) # model betas estimates
        
        betas_mvnorm <- MASS::mvrnorm(n = ra_nsim, mu = betas, Sigma = cluster_robust_vcov) # Simulate from the multivariate normal distribution
        
        b1 <- betas_mvnorm[,1]
        b2 <- betas_mvnorm[,2]
        b3 <- betas_mvnorm[,3]
        b4 <- betas_mvnorm[,4]
        
        # Compute AUC for each simulation sample
        lhs_aucs <- sapply(b4, lhs_auc)
        rhs_aucs <- mapply(rhs_auc, beta_A = b2, beta_B = b3, beta_C = b1)
        
        # Compute the difference in AUCs
        delta_aucs <- lhs_aucs - rhs_aucs
        
        ratio_aucs <- lhs_aucs/rhs_aucs
        
        # Estimates
        
        delta <- median(delta_aucs)
        ratio <- median(ratio_aucs)
        sd_delta <- sd(delta_aucs)
        
        # 95% Confidence interval
        ci_delta <- quantile(delta_aucs, c(0.025, 0.975))
        
        ci_ratio <- quantile(ratio_aucs, c(0.025, 0.975))
        
        # p-value (two-tailed test)
        p_delta <-  2 * min(mean(delta_aucs <= 0), mean(delta_aucs >= 0))
        
        p_ratio <-  2 * min(mean(ratio_aucs <= 1), mean(ratio_aucs >= 1))

      } else {
        
        vcov_mtx <- vcov(model_time) # Variance-Covariance Matrix for the Fitted Model Object
        betas <- fixef(model_time) # model betas estimates
        
        betas_mvnorm <- MASS::mvrnorm(n = ra_nsim, mu = betas, Sigma = vcov_mtx) # Simulate from the multivariate normal distribution
        
        b1 <- betas_mvnorm[,1]
        b2 <- betas_mvnorm[,2]
        b3 <- betas_mvnorm[,3]
        b4 <- betas_mvnorm[,4]
        
        # Compute AUC for each simulation sample
        lhs_aucs <- sapply(b4, lhs_auc)
        rhs_aucs <- mapply(rhs_auc, beta_A = b2, beta_B = b3, beta_C = b1)
        
        # Compute the difference in AUCs
        delta_aucs <- lhs_aucs - rhs_aucs
        
        ratio_aucs <- lhs_aucs/rhs_aucs
        
        # Estimates
        
        delta <- median(delta_aucs)
        ratio <- median(ratio_aucs)
        
        sd_delta <- sd(delta_aucs)
        
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
                    - (Test$estimate)/(Test$std.error),
                    - (Test$conf.high)/(Test$std.error),
                    - (Test$conf.low)/(Test$std.error),
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