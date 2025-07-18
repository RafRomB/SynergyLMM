% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmmModel.R
\name{lmmModel}
\alias{lmmModel}
\title{Linear Mixed Effect Model for Tumor Growth}
\usage{
lmmModel(
  data,
  grwth_model = "exp",
  sample_id = "SampleID",
  time = "Time",
  treatment = "Treatment",
  tumor_vol = "TV",
  trt_control = "Control",
  drug_a = "Drug_A",
  drug_b = "Drug_B",
  drug_c = NA,
  combination = "Combination",
  time_start = NULL,
  time_end = NULL,
  min_observations = 1,
  show_plot = TRUE,
  tum_vol_0 = "ignore",
  start_values = c(0.05, 0.01),
  ...
)
}
\arguments{
\item{data}{A data frame with the tumor growth data, in long format.
It must contain at least the following columns: mice IDs, time of follow-up (numeric number), treatment and tumor volume (numeric number).}

\item{grwth_model}{Tumor growth model to use. Possible options are "exp", for exponential tumor growth model, or "gompertz" for Gompertz tumor growth model.}

\item{sample_id}{String indicating the name of the column with the mice IDs.}

\item{time}{String indicating the name of the column with the time of follow-up.}

\item{treatment}{String indicating the name of the column with the treatment corresponding to each sample.}

\item{tumor_vol}{String indicating the name of the column with the tumor volume (or any other measurement representing the tumor growth).}

\item{trt_control}{String indicating the name assigned to the 'Control' group.}

\item{drug_a}{String indicating the name assigned to the 'Drug A' group.}

\item{drug_b}{String indicating the name assigned to the 'Drug B' group.}

\item{drug_c}{String indicating the name assigned to the 'Drug C' group (if present).}

\item{combination}{String indicating the name assigned to the Combination ('Drug A' + 'Drug B', or 'Drug A' + 'Drug B' + 'Drug C') group.}

\item{time_start}{Numeric value indicating the time point at which the treatment started. If not
specified, the minimum value in the \code{time} column is used as the starting time point.}

\item{time_end}{Numeric value indicating the last time point to be included in the analysis. If not
specified, the maximum value in the \code{time} column is used as the final time point.}

\item{min_observations}{Minimum number of observation for each sample to be included in the analysis.}

\item{show_plot}{Logical indicating if a plot for the log of the relative tumor volume (RTV) vs Time for each sample,
and the model calculated marginal slope for each treatment, should be produced.}

\item{tum_vol_0}{Select the behavior of the function regarding measurements in which the tumor measurement is 0, and therefore the logarithmic
transformation is not possible. Possible options are 'ignore' (default), to ignore these measurements, or 'transform', to add 1 unit to all
measurements before the log transformation.}

\item{start_values}{Numeric vector of length 2 with initial estimates for the fixed effects of the Gompertz model (\eqn{r_0, \rho}). It can be set
to "selfStart", in which case initial values will be derived from a Gompertz model derived from a call to \code{stats::\link[stats:nls]{nls}}}

\item{...}{Additional arguments to be passed to \code{nlme::\link[nlme:lme]{lme}} or \code{nlme::\link[nlme:nlme]{nlme}}.}
}
\value{
An object of class "lme" (see \code{nlme::\link[nlme:lme]{lme}} for details) or "nlme" (see \code{nlme::\link[nlme:nlme]{nlme}} for details) representing the linear mixed-effects model fit.
If \code{show_plot = TRUE}, the plot
of the tumor growth data obtained with \code{\link[=plot_lmmModel]{plot_lmmModel()}} is also shown.
}
\description{
\code{lmmModel()} fits a mixed effect model from a tumor growth dataset. The input data frame must be in long format and include at least the following columns: column with the sample ids,
column with the time at which each measurement has been done, a column indicating the treatment group, and a column with the tumor measurement (e.g. tumor volume).
}
\details{
\code{lmmModel()} fits a mixed effect model describing the tumor growth. Currently, two growth kinetics models are available: exponential growth, and Gompertz growth.

\strong{Exponential}

\code{lmmModel()} will fit a linear mixed-effect model (LMM) assuming that the tumor growth follows an exponential kinetics. Any departure from this assumption can be tested using the diagnostics functions \code{\link[=ranefDiagnostics]{ranefDiagnostics()}},
\code{\link[=residDiagnostics]{residDiagnostics()}}, and \code{\link[=ObsvsPred]{ObsvsPred()}}.

The model formula for the LMM following the \strong{exponential tumor growth} is:
\deqn{\log RTV_{i}(t) = \beta_{T_i} \cdot t + b_i \cdot t + \varepsilon_{i} (t).}
\itemize{
\item \eqn{\log RTV_{i}(t)} denotes the value of the logarithm of the relative tumor volume measured for subject \eqn{i} at time \eqn{t}.
\item \eqn{\beta_{T_i}} represents the fixed effects at time \eqn{t} for each  treatment \eqn{T_i}, where \eqn{T_i \in \{Control, DrugA, DrugB, Combination\}} in the case of
two-drugs combination experiments, or \eqn{T_i \in \{Control, DrugA, DrugB, DrugC, Combination\}} in the case of three-drugs combination experiments, and indicates the
tumor-specific growth rate for each treatment group.
\item \eqn{b_i \cdot t} corresponds to the subject-specific random slope that takes into account the longitudinal nature of the data, where \eqn{b_i \sim \mathcal{N}(0,\sigma^2_b)} is the random effect for subject \eqn{i}.
\item \eqn{\varepsilon_{i}(t) \sim \mathcal{N}(0,\sigma^2)} is the residual error term.
}

\strong{Gompertz}

\code{lmmModel()} will fit a non-linear mixed effect (NLME) model assuming the tumor growth follows a Gompertz growth kinetics. Any departure from this assumption can be tested using the diagnostics functions \code{\link[=ranefDiagnostics]{ranefDiagnostics()}},
\code{\link[=residDiagnostics]{residDiagnostics()}}, and \code{\link[=ObsvsPred]{ObsvsPred()}}.

The model formula for the non-linear mixed-effect model following the \strong{Gompertz tumor growth} is:
\deqn{\log RTV_{i}(t) = \frac{r_{0_{T_i}}+b_{0_i}}{\rho_{T_i}+b_{1_i}}(1-e^{-(\rho_{T_i}+b_{1_i})\cdot t})+\varepsilon_{i}(t).}
\itemize{
\item \eqn{\log RTV_{i}(t)} denotes the value of the logarithm of the relative tumor volume measured for subject \eqn{i} at time \eqn{t}.
\item \eqn{T_i \in \{Control, DrugA, DrugB, Combination\}} in the case of
two-drugs combination experiments, or \eqn{T_i \in \{Control, DrugA, DrugB, DrugC, Combination\}} in the case of three-drugs combination experiments,
indicates the treatment groupt of individual \eqn{i}.
\item \eqn{r_{0_{T_i}}} is the fixed effect for the initial growth rate for treatment group \eqn{T_i}.
\item \eqn{\rho_{T_i}} is the fixed effect for the constant accounting for the reduction in the tumor growth rate for treatment group \eqn{T_i}.
\item \eqn{b_{0_i} \sim \mathcal{N}(0,\sigma^2_{{r_0}})} is the random effect on \eqn{r_0} for individual \eqn{i}.
\item \eqn{b_{1_i} \sim \mathcal{N}(0,\sigma^2_{{\rho}})} is the random effect on \eqn{\rho} for individual \eqn{i}.
\item \eqn{\varepsilon_{i}(t) \sim \mathcal{N}(0,\sigma^2)} is the residual error term.
}

The implementation of the linear mixed model in \code{lmmModel()} is done using \code{nlme::\link[nlme:lme]{lme}} for the exponential model, or \code{nlme::\link[nlme:nlme]{nlme}}
for the Gompertz model. Both also allow for the
specification of within-group correlations structures and/or unequal variances. These, and additional parameters,
can be passed to the \code{nlme::\link[nlme:lme]{lme}} or \code{nlme::\link[nlme:nlme]{nlme}} functions through the \code{...} argument for fitting the model (see examples below).
}
\examples{
data("grwth_data")
# Most simple model
lmmModel(
 data = grwth_data,
 sample_id = "subject",
 time = "Time",
 treatment = "Treatment",
 tumor_vol = "TumorVolume",
 trt_control = "Control",
 drug_a = "DrugA",
 drug_b = "DrugB",
 combination = "Combination"
 )
# Changing the last time point of follow-up
lmmModel(
 data = grwth_data,
 sample_id = "subject",
 time = "Time",
 treatment = "Treatment",
 tumor_vol = "TumorVolume",
 trt_control = "Control",
 drug_a = "DrugA",
 drug_b = "DrugB",
 combination = "Combination",
 time_end = 21
 )
# Adding additional parameters for model fitting
lmmModel(
 data = grwth_data,
 sample_id = "subject",
 time = "Time",
 treatment = "Treatment",
 tumor_vol = "TumorVolume",
 trt_control = "Control",
 drug_a = "DrugA",
 drug_b = "DrugB",
 combination = "Combination",
 # Adding variance function to represent a different variance per subject
 weights = nlme::varIdent(form = ~1|SampleID),
 # Specifiying control values for lme Fit (useful when convergence problems appear)
 control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 100, msMaxEval = 1000)
 )

# Fit a model specifying a different variance per Time
lmmModel(
 data = grwth_data,
 sample_id = "subject",
 time = "Time",
 treatment = "Treatment",
 tumor_vol = "TumorVolume",
 trt_control = "Control",
 drug_a = "DrugA",
 drug_b = "DrugB",
 combination = "Combination",
 # Adding variance function to represent a different variance per Time
 weights = nlme::varIdent(form = ~1|Time)
 )

# Fit a model using the Gompertz model

lmmModel(
data = grwth_data, 
grwth_model = "gompertz", # Selecting Gompertz model
sample_id = "subject",
time = "Time",
treatment = "Treatment",
tumor_vol = "TumorVolume",
trt_control = "Control",
drug_a = "DrugA",
drug_b = "DrugB",
combination = "Combination"
)


}
\references{
\itemize{
\item Pinheiro JC, Bates DM (2000). \emph{Mixed-Effects Models in S and S-PLUS}. Springer, New York. doi:10.1007/b98882 \doi{doi:10.1007/b98882}.
\item Pinheiro J, Bates D, R Core Team (2024). \emph{nlme: Linear and Nonlinear Mixed Effects Models}. R package version 3.1-166, \url{https://CRAN.R-project.org/package=nlme}.
\item Andrzej Galecki & Tomasz Burzykowski (2013) \emph{Linear Mixed-Effects Models Using R: A Step-by-Step Approach} First Edition. Springer, New York. ISBN 978-1-4614-3899-1
}
}
