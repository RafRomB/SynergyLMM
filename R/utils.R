#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL

#' @title Helper function to simulate tumor growth data for a two-drug combination experiment.
#' @param npg Number of samples per group.
#' @param timepoints Vector with the time points at which the tumor volume measurements have been performed.
#' @param initial_volume Initial volume for the tumor growth.
#' @param grwrControl Coefficient for Control treatment group tumor growth rate.
#' @param grwrA Coefficient for Drug A treatment group tumor growth rate.
#' @param grwrB Coefficient for Drug B treatment group tumor growth rate.
#' @param grwrComb Coefficient for Combination (Drug A + Drug B) treatment group tumor growth rate.
#' @param sd Variability for the tumor growth.
#' @details 
#' The function simulates the tumor growth following exponential kinetics,
#' given by
#' \deqn{TV(t) = TV_0 \cdot e^{\beta_i \cdot t}}
#' 
#' where \eqn{TV_0} is given by `initial_volume`, \eqn{t} is given by `timepoints` 
#' and \eqn{\beta_i} are the coefficients given by `grwrControl`, `grwrA`, `grwrB`, and `grwrComb`.
#' 
#' The variability is simulated using the `sd` argument to add random noise from a normal distribution \eqn{N(1, SD)},
#' with \eqn{SD} = `sd`.
#' @returns A data frame with tumor growth data in long format.
#' @examples
#' # This code generates the 'grwth_data' example dataset:
#' set.seed(123)
#' grwth_data <- simulateTumorGrowth(npg = 8, timepoints = seq(0,30,3), 
#'                                    initial_volume = 200, grwrControl = 0.08,
#'                                    grwrA = 0.07, grwrB = 0.065, 
#'                                    grwrComb = 0.04, sd = 0.2)
#' 
#' @export
simulateTumorGrowth <- function(npg = 5,
                                 timepoints = c(0, 3, 5, 10),
                                 initial_volume = 100,
                                 grwrControl = 0.08,
                                 grwrA = 0.07,
                                 grwrB = 0.06,
                                 grwrComb = 0.04,
                                 sd = 0.1) {
  
  subject <- 1:(4 * npg) # Subjects' ids
  Treatment <- gl(4, npg, labels = c("Control", "DrugA", "DrugB", "Combination")) # Treatment for each subject
  dts <- data.frame(subject, Treatment) # Subject-level data
  
  dtL <- list(Time = timepoints, subject = subject)
  dtLong <- expand.grid(dtL) # Long format
  mrgDt <- merge(dtLong, dts, sort = FALSE) # Merged
  
  mrgDt$TumorVolume <- NULL
  
  # Simulate exponential growth for each subject
  mrgDt$TumorVolume[mrgDt$Treatment == "Control"] <- initial_volume * exp(grwrControl * timepoints)
  mrgDt$TumorVolume[mrgDt$Treatment == "DrugA"] <- initial_volume * exp(grwrA * timepoints)
  mrgDt$TumorVolume[mrgDt$Treatment == "DrugB"] <- initial_volume * exp(grwrB * timepoints)
  mrgDt$TumorVolume[mrgDt$Treatment == "Combination"] <- initial_volume * exp(grwrComb * timepoints)
  
  # Add random noise to simulate variability
  mrgDt$TumorVolume <- mrgDt$TumorVolume * rnorm(nrow(mrgDt), mean = 1, sd = sd)
  mrgDt <- as.data.frame(mrgDt)
  return(mrgDt)
}
