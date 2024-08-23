
## Synergy Calculation using LMM

#' @importFrom marginaleffects hypotheses
#' 
#' @title Synergy calculation using linear-mixed models.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`LMM_Model()`].
#' @param method String indicating the method for synergy calculation. Possible methods are "Bliss", "HSA" and "RA",
#' corresponding to Bliss, highest single agent and response additivity, respectively.
#' @param min_day Minimun day from which to start calculating synergy.
#' @param robustSE If TRUE, uncertainty is estimated using robust standard errors 
#' using a sandwich estimate of the variance-covariance matrix of the regression coefficient estimates provided by [clubSandwich::vcovCR()].
#' @param type Character string specifying which small-sample adjustment should be used, with available options "CR0", "CR1", "CR1p", "CR1S", "CR2", or "CR3". 
#' See "Details" section of [clubSandwich::vcovCR()] for further information.
#' @param norm_test If `method` is set to "RA", string indicating the test for checking the 
#' normality of the hypothesis expression of synergy (see "Warning#3" in "Description" section of [marginaleffects::hypotheses()]).
#' Possible values are "shapiroTest", "dagoTest" or "adTest", for Shapiro-Wilk, D'Agostino and Anderson-Darling normality tests, respectively.
#' @param ... Additional arguments to be passed to [marginaleffects::hypotheses()].
#' @returns The function returns a list with two elements:
#' - `Constrasts`: List with the outputs of the (non)-linear test for the synergy null hypothesis obtained by [marginaleffects::hypotheses()] for each day.
#' See [marginaleffects::hypotheses()] for more details.
#' - `Synergy`: Data frame with the synergy results, indicating the model of synergy ("Bliss", "HSA" or "RA"), the metric (combination index and synergy score),
#' the value of the metric estimate (with upper and lower confidence intervals) and the p-value, for each day.
#' @export

LMM_Syn <- function(model,
                    method = "Bliss",
                    min_day = 0,
                    robustSE = FALSE,
                    type = "CR2",
                    norm_test = "shapiroTest",
                    ...) {
  
  ci <- data.frame()
  ss <- data.frame()
  Contrasts <- list()
  if (method == "Bliss") {
    contrast <- "b4 = b2 + b3 - b1"
  }
  if (method == "HSA") {
    fixef_betas <- nlme::fixef(model)[2:3]
    if (which.min(fixef_betas) == 1) {
      contrast <- "b4 = b2"
    } else{
      contrast <- "b4 = b3"
    }
  }
  if(method == "RA") {
    contrast <- "b4 = log(exp(b2) + exp(b3) - exp(b1))"
    
    summ <- summary(model)
    
    b1 <- summ$tTable[1, c("Value", "Std.Error")]
    b2 <- summ$tTable[2, c("Value", "Std.Error")]
    b3 <- summ$tTable[3, c("Value", "Std.Error")]
    
    b1 <- rnorm(n = 1000, mean = b1["Value"], sd = b1["Std.Error"])
    b2 <- rnorm(n = 1000, mean = b2["Value"], sd = b2["Std.Error"])
    b3 <- rnorm(n = 1000, mean = b3["Value"], sd = b3["Std.Error"])
    b4 <- log(exp(b2) + exp(b3) - exp(b1))
    par(mfrow = c(1, 2))
    hist(log(exp(b2) + exp(b3) - exp(b1)), main = "Histogram for RA\nexpression values")
    qqnorm(log(exp(b2) + exp(b3) - exp(b1)))
    qqline(b4)
    par(mfrow = c(1, 1))
    if (norm_test == "shapiroTest") {
      print(fBasics::shapiroTest(b4))
    }
    if (norm_test == "dagoTest") {
      print(fBasics::dagoTest(b4))
    }
    if (norm_test == "adTest") {
      print(fBasics::adTest(b4))
    }
  }
  
  days <- unique(model$data$Day)
  days <- days[days >= min_day]
  days <- days[order(days)]
  i <- 1
  for (d in days) {
    data <- model$data %>% dplyr::filter(.data$Day <= d)
    model_day <- update(model, data = data)
    if (robustSE) {
      Test <- hypotheses(
        model_day,
        hypothesis = contrast,
        vcov = clubSandwich::vcovCR(model, type = type),
        ...
      )
    } else {
      Test <- hypotheses(model_day, hypothesis = contrast, ...)
    }
    ci <- rbind(ci, data.frame(
      method,
      "CI",
      exp(Test$estimate),
      exp(Test$conf.low),
      exp(Test$conf.high),
      Test$p.value,
      d
    ))
    ss <- rbind(ss,
                data.frame(
                  method,
                  "SS",
                  -100 * (Test$estimate),
                  -100 * (Test$conf.low),
                  -100 * (Test$conf.high),
                  Test$p.value,
                  d
                ))
    Contrasts[[i]] <- Test
    i <- i + 1
  }
  names(Contrasts) <- paste("Day", days, sep = "")
  colnames(ci) <- c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Day")
  colnames(ss) <- c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Day")
  df <- rbind(ci, ss)
  result <- list(Contrasts = Contrasts, Synergy = df)
  print(Plot_LMM_Syn(result))
  return(result)
}