% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_ObsvsPred.R
\name{plot_ObsvsPred}
\alias{plot_ObsvsPred}
\title{Plots of Observed vs Predicted Values}
\usage{
plot_ObsvsPred(model, nrow = 4, ncol = 5)
}
\arguments{
\item{model}{An object of class "lme" representing the linear mixed-effects model fitted by \code{\link[=lmmModel]{lmmModel()}}.}

\item{nrow}{Number of rows of the layout to organize the observed vs predicted plots.}

\item{ncol}{Number of columns of the layout to organize the observed vs predicted plots.}
}
\value{
A layout (arranged in \code{nrow} rows and \code{ncol} columns) of the observed and predicted values of \eqn{log}(relative tumor volume) vs Time for each SampleID (i.e. subject),
with the actual measurements, the regression line for each SampleID, and the marginal, treatment-specific,
regression line for each treatment group.
}
\description{
Visualization of observed vs predicted values by a fitted linear mixed model of tumor growth data.
}
\examples{
#' data(grwth_data)
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
# Obtain the plots
plot_ObsvsPred(lmm, nrow = 4, ncol = 8)    

}
