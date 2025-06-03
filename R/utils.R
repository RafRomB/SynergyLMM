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


#' @title Gompertz log model function to calculate starting values
#' @description
#' Define the Gompertz log model function
#' @param Time Time variable in the data
#' @param r0 First parameter of the model, representing the initial tumor growth rate
#' @param rho Second parameter of the mode, representing a constant that corrects the growth rate with time
#' @returns logRTV based on the Gompertz log model for the given parameters
#' @keywords internal
#' @noRd

gompertzLog_fun <- function(Time, r0, rho) {
  (r0 / rho) * (1 - exp(-rho * Time))
}

#' @title Definition of self-starting funtion for start value of Gompertz log model
#' @description
#' Definition of the self-starting function to estimate initial values for Gompertz log model
#' @param mCall Matched call to the `selfStart` model
#' @param LHS Expression on the left-hand side of the model formula in the call to [stats::nls]
#' @param data `data.frame` in which to find the variable named in the other two arguments
#' @keywords internal
#' @noRd
gompertzLog_init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["Time"]], LHS, data)
  Time <- xy$x
  logRTV <- xy$y
  
  # Simple estimation:
  # Use initial slope for r0 and asymptote for r0 / rho => estimate rho
  r0_start <- (logRTV[2] - logRTV[1]) / (Time[2] - Time[1])  # crude initial slope
  asymptote <- max(logRTV)  # assume saturation
  rho_start <- r0_start / asymptote  # rearranged from model: asymptote = r0 / rho
  
  value <- c(r0 = r0_start, rho = rho_start)
  names(value) <- mCall[c("r0", "rho")]
  value
}

#' @title Create selfStart model
#' @description
#' Self-starting Gompertz log model
#' @keywords internal
#' @noRd

SSgompertzLog <- selfStart(
  model = gompertzLog_fun,
  initial = gompertzLog_init,
  parameters = c("r0", "rho")
)


