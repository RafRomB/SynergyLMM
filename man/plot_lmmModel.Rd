% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_lmmModel.R
\name{plot_lmmModel}
\alias{plot_lmmModel}
\title{Plotting of tumor growth data from a fitted model}
\usage{
plot_lmmModel(
  model,
  trt_control = "Control",
  drug_a = "Drug_A",
  drug_b = "Drug_B",
  drug_c = NA,
  combination = "Combination"
)
}
\arguments{
\item{model}{An object of class "lme" representing the linear mixed-effects model fitted by \code{\link[=lmmModel]{lmmModel()}}.}

\item{trt_control}{String indicating the name assigned to the 'Control' group.}

\item{drug_a}{String indicating the name assigned to the 'Drug A' group.}

\item{drug_b}{String indicating the name assigned to the 'Drug B' group.}

\item{drug_c}{String indicating the name assigned to the 'Drug C' group (if present).}

\item{combination}{String indicating the name assigned to the Combination ('Drug A' + 'Drug B', or 'Drug A' + 'Drug B' + 'Drug C') group.}
}
\value{
A ggplot2 plot (see \code{\link[ggplot2:ggplot]{ggplot2::ggplot()}} for more details) showing the tumor growth data represented as log(relative tumor volume) versus time since treatment initiation.
The regression lines corresponding to the fixed effects for each treatment group are also plotted.
}
\description{
Vizualization of tumor growth data and linear mixed model fitted regression line for the fixed effects. This functions returns a \href{https://ggplot2.tidyverse.org/}{ggplot2} plot, allowing for
further personalization.
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
# Default plot
plot_lmmModel(lmm,
trt_control = "Control",
drug_a = "DrugA",
drug_b = "DrugB",
combination = "Combination"
)
# Adding ggplot2 elements
plot_lmmModel(lmm,
trt_control = "Control",
drug_a = "DrugA",
drug_b = "DrugB",
combination = "Combination"
) + ggplot2::labs(title = "Example Plot") + ggplot2::theme(legend.position = "top")

}
