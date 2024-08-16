## Power Evaluation using Simulations 

#' @export

Pwr_simulate <- function(model, nsim=1000, method = "Bliss", pvalue = 0.05){
  
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
  }
  
  simA <- nlmeU::simulateY(model, nsim = nsim) # Simulation
  dt <- model$data # working copy
  simfmA <- apply(simA,
                  2,
                  function(y){
                    dt$logRTV <- y
                    auxFit <- update(model, data = dt)
                    marginaleffects::hypotheses(auxFit, hypothesis = contrast)
                  })
  FstateE <- 
    sapply(simfmA, function(x) x$p.value)
  
  powerE <- sum(FstateE < pvalue)/nsim
  return(powerE)
}
