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

#' @title Definition of self-starting function for start value of Gompertz log model
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



#' @title Warn when simulated p-values are approximated to zero
#' @description
#' Emits the warning used by [lmmSynergy()] when one or more of the reported
#' p-values is exactly zero. Simulation-based p-values cannot be smaller than
#' the resolution allowed by `nsim`, so a zero simply means that every simulated
#' draw fell on the same side of the null.
#' @param pvals Numeric vector of p-values.
#' @param nsim Number of simulations used to obtain `pvals`.
#' @param gompertz Logical. Whether the p-values come from a Gompertz model fit,
#' which changes the advice given in the warning.
#' @return Invisibly, `TRUE` if a warning was emitted and `FALSE` otherwise.
#' @keywords internal
#' @noRd

warn_zero_pval <- function(pvals, nsim, gompertz = FALSE) {
  if (!any(pvals == 0, na.rm = TRUE)) {
    return(invisible(FALSE))
  }
  apx_p <- approx_pval_label(nsim)
  advice <- if (gompertz) {
    " If you used a Gompertz model, consider increasing 'nsim' value for more precise p-values."
  } else {
    " If you used method = 'RA' consider increasing 'nsim' value for more precise p-values."
  }
  warning(paste("p-values below", apx_p, "are approximated to 0."), advice,
          call. = FALSE)
  invisible(TRUE)
}

#' @title Label for the resolution of simulated p-values
#' @description
#' Builds the `"p<1e-03"` style label describing the smallest p-value that
#' `nsim` simulations can resolve.
#' @param nsim Number of simulations.
#' @return A character string.
#' @keywords internal
#' @noRd

approx_pval_label <- function(nsim) {
  paste0("p<", format(1 / nsim, scientific = TRUE, digits = 1))
}
