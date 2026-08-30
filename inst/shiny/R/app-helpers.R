# Helpers for the SynergyLMM Shiny app.
#
# Shiny auto-sources every file in this directory before ui.R/server.R, so these are available
# unqualified in server.R. They are *presentation and plumbing only*: every number shown by the app
# is computed by the installed SynergyLMM package, so the app can never drift from it.


#' Fit the model with the optimizer controls exposed in the UI.
#'
#' `SynergyLMM::lmmModel()` forwards `...` to `nlme::lme()` / `nlme::nlme()`, so the control object
#' and the variance function are simply passed through.
app_lmmModel <- function(...,
                         grwth_model = "exp",
                         maxIter = 50,
                         msMaxIter = 50,
                         niterEM = 25,
                         msMaxEval = 200,
                         pnlsMaxIter = 7,
                         opt = "nlminb",
                         weights = NULL) {

  if (grwth_model == "exp") {
    cntrl <- nlme::lmeControl(
      maxIter = maxIter,
      msMaxIter = msMaxIter,
      niterEM = niterEM,
      msMaxEval = msMaxEval,
      opt = opt
    )
  } else {
    cntrl <- nlme::nlmeControl(
      maxIter = maxIter,
      msMaxIter = msMaxIter,
      niterEM = niterEM,
      pnlsMaxIter = pnlsMaxIter,
      opt = opt
    )
  }

  SynergyLMM::lmmModel(
    ...,
    grwth_model = grwth_model,
    control = cntrl,
    weights = weights
  )
}


#' Rename an estimates table to the column names the app displays.
#'
#' The package returns machine-friendly names (`se_*`, `sd_*`); the app shows human-readable ones.
#' This is a rename only - the values are the package's.
app_estimate_names <- function(est) {
  nms <- names(est)
  nms <- sub("^se_", "SD ", nms)
  nms[nms == "sd_ranef"] <- "SD Random Effects"
  nms[nms == "sd_resid"] <- "SD Residuals"
  nms[nms == "sd_r0_ranef"] <- "SD r0 Random Effects"
  nms[nms == "sd_rho_ranef"] <- "SD rho Random Effects"
  names(est) <- nms
  est
}


#' Model estimates with the column names the app displays.
app_lmmModel_estimates <- function(model, ...) {
  app_estimate_names(SynergyLMM::lmmModel_estimates(model, ...))
}


#' Retrieve the plot drawn by the three power-analysis functions.
#'
#' They return a data frame and attach the assembled plot as the "plot" attribute, which is what the
#' download handlers need to hand to `ggplot2::ggsave()`.
app_pwr_plot <- function(x) {
  p <- attr(x, "plot")
  if (is.null(p)) {
    stop("No plot is available for this power analysis result.", call. = FALSE)
  }
  p
}
