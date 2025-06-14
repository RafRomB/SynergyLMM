% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_ranefDiagnostics.R
\name{plot_ranefDiagnostics}
\alias{plot_ranefDiagnostics}
\title{Plots for random effects diagnostics}
\usage{
plot_ranefDiagnostics(model)
}
\arguments{
\item{model}{An object of class "lme" representing the linear mixed-effects model fitted by \code{\link[=lmmModel]{lmmModel()}}.}
}
\value{
A list with different plots for evaluating the normality and homoscedasticity of the random effects, including:
\itemize{
\item A normal Q-Q plot of the random effects of the model.
\item A normal Q-Q plot of the residuals by sample.
\item Boxplots of the "raw" residuals (observed - fitted) by sample.
\item Scatter plots of the normalized residuals (standardized residuals pre-multiplied by the inverse square-root factor of the estimated error correlation matrix, see \link[nlme:residuals.lme]{nlme::residuals.lme})
vs fitted values by sample. Observations with absolute standardized (normalized) residuals greater than the \eqn{1-0.05/2} quantile of the standard normal distribution
are identified in the plots labelled with the time point corresponding to the observation.
}
}
\description{
Visualization of random effects diagnostics for a fitted linear mixed model of tumor growth data.
}
\examples{
data(grwth_data)
# Fit the model
lmm <- lmmModel(
  data = grwth_data,
  sample_id = "subject",
  time = "Time",
  treatment = "Treatment",
  tumor_vol = "TumorVolume",
  trt_control = "Control",
  drug_a = "DrugA",
  drug_b = "DrugB",
  combination = "Combination",
  show_plot = FALSE
  )
# Generate plots 
plot_ranefDiagnostics(lmm)
# Access to specific plots
plot_ranefDiagnostics(lmm)$plots[[1]]
plot_ranefDiagnostics(lmm)$plots[[2]]
}
