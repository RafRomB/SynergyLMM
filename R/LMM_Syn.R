
## Synergy Calculation using LMM

#' @importFrom marginaleffects hypotheses
#' @export

LMM_Syn <- function(model, method = "Bliss", min = 0, robustSE = FALSE, type = "CR2", test = "shapiroTest", ...){
  
  ci <- data.frame()
  ss <- data.frame()
  Contrasts <- list()
  if(method == "Bliss"){
    contrast <- "b4 = b2 + b3 - b1"
  }
  if(method == "HSA" ){
    fixef_betas <- nlme::fixef(model)[2:3]
    if(which.min(fixef_betas) == 1){
      contrast <- "b4 = b2"
    } else{
      contrast <- "b4 = b3"
    }
  }
  if(method == "RA"){
    contrast <- "b4 = log(exp(b2) + exp(b3) - exp(b1))"
    
    summ <- summary(model)
    
    b1 <- summ$tTable[1, c("Value", "Std.Error")]
    b2 <- summ$tTable[2, c("Value", "Std.Error")]
    b3 <- summ$tTable[3, c("Value", "Std.Error")]
    
    b1 <- rnorm(n = 1000, mean = b1["Value"], sd = b1["Std.Error"])
    b2 <- rnorm(n = 1000, mean = b2["Value"], sd = b2["Std.Error"])
    b3 <- rnorm(n = 1000, mean = b3["Value"], sd = b3["Std.Error"])
    b4 <- log(exp(b2)+exp(b3)-exp(b1))
    par(mfrow = c(1,2))
    hist(log(exp(b2)+exp(b3)-exp(b1)), main = "Histogram for RA\nexpression values")
    qqnorm(log(exp(b2)+exp(b3)-exp(b1)))
    qqline(b4)
    par(mfrow = c(1,1))
    if(test == "shapiroTest"){
      print(fBasics::shapiroTest(b4))
    }
    if(test == "dagoTest"){
      print(fBasics::dagoTest(b4))
    }
    if(test == "adTest"){
      print(fBasics::adTest(b4))
    }
  }
  
  days <- unique(model$data$Day)
  days <- days[days>=min]
  days <- days[order(days)]
  i <- 1
  for(d in days){
    data <- model$data %>% dplyr::filter(Day <= d)
    model_day <- update(model, data = data)
    if(robustSE){
      Test <- hypotheses(model_day, hypothesis = contrast, vcov = clubSandwich::vcovCR(model, type = type), ...)
    } else {
      Test <- hypotheses(model_day, hypothesis = contrast, ...)
    }
    ci <- rbind(ci, data.frame(method, "CI", exp(Test$estimate), exp(Test$conf.low), exp(Test$conf.high), Test$p.value, d))
    ss <- rbind(ss, data.frame(method, "SS", -100*(Test$estimate), -100*(Test$conf.low), -100*(Test$conf.high), Test$p.value, d))
    Contrasts[[i]] <- Test
    i <- i+1
  }
  names(Contrasts) <- paste("Day",days, sep = "")
  colnames(ci) <- c("Model","Metric","Estimate", "lwr", "upr", "pval", "Day")
  colnames(ss) <- c("Model","Metric","Estimate", "lwr", "upr", "pval", "Day")
  df <- rbind(ci, ss)
  result <- list(Contrasts = Contrasts, Synergy = df)
  return(result)
}