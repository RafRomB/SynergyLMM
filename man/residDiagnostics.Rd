% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residDiagnostics.R
\name{residDiagnostics}
\alias{residDiagnostics}
\title{Diagnostics of residuals of the linear mixed model}
\usage{
residDiagnostics(
  model,
  pvalue = 0.05,
  norm_test = "shapiroTest",
  verbose = TRUE
)
}
\arguments{
\item{model}{An object of class "lme" representing the linear mixed-effects model fitted by \code{\link[=lmmModel]{lmmModel()}}.}

\item{pvalue}{Threshold for the p-value of outlier observations based on their normalized residuals.}

\item{norm_test}{String indicating the function for testing the normality of the normalized residuals. A collection of functions from
\code{fBasics::\link[fBasics:normalTest]{normalTest}} is available. We recommend using one of "shapiroTest", "dagoTest",
or "adTest" for performing Shapiro - Wilk,
D'Agostino, or Anderson - Darling normality test, respectively.}

\item{verbose}{Logical indicating if the normality and homoscedasticity tests results, and the list of potential
outlier observations should be printed to the console.}
}
\value{
A list with different elements for the diagnostics of the residuals are produced:
\itemize{
\item \code{plots}: Different plots for evaluating the normality and homocedasticity of the residuals.
\item \code{outliers}: Data frame with the identified outliers based on the Pearson residuals and the value of \code{pval}. The column \code{resid.p} contains the
value of the Pearson residuals for each observation.
\item \code{Normality}: Results from the test of the normality of the normalized residuals of the model.
\item \code{Levene.test}: List with the Levene homoscedasticity test results of the normalized residuals by Time and Treatment.
\item \code{Fligner.test}: List with the Fligner-Killeen homoscedasticity test results of the normalized residuals by Time and Treatment.
}
}
\description{
\code{residDiagnostics} provides several plots as well as statistical test for the examination
of the normality and homoscedasticity of the residuals of the input model.
}
\details{
One of the assumption of the model fit by \code{\link[=lmmModel]{lmmModel()}} is that the residuals are normally distributed.
For the evaluation of this assumption, \code{residDiagnostics} provides Q-Q plots of the normalized residuals
(standardized residuals pre-multiplied by the inverse square-root factor of the estimated error correlation matrix, see \link[nlme:residuals.lme]{nlme::residuals.lme}),
together with statistical assessment of their
normality using Shapiro-Wilk, D'Agostini and Anderson-Darling normality tests. Additionally, Q-Q plots of the normalized residuals by time point and
treatment group are provided to be able to detect time points or treatment groups which could be notably different from the others and be
affecting the adequacy of the model.

Scatter plots of the normalized residuals versus fitted values and normalized residuals
per time and per treatment are also provided to give information about variability of the residuals and possible outlier observations. These plots are accompanied
by Levene and Fligner-Killend homogeneity of variance test results.

Observations with absolute standardized (normalized) residuals greater than the \eqn{1-0.05/2} quantile of the standard normal distribution
are identified and reported as potential outlier observations.
}
\examples{
# Load the example data
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
  combination = "Combination"
  )
  
# Residuals diagnostics
resid_diag <- residDiagnostics(model = lmm, pvalue = 0.05)

# Access outliers data frame
resid_diag$Outliers

# Access individual plots
resid_diag$Plots[1]
resid_diag$Plots[2]

# Access results of normality tests
resid_diag$Normality

# Access to homoscedasticity test results

resid_diag$Levene.test

resid_diag$Fligner.test

}
\references{
\itemize{
\item Pinheiro JC, Bates DM (2000). \emph{Mixed-Effects Models in S and S-PLUS}. Springer, New York. \doi{doi:10.1007/b98882}.
\item Andrzej Galecki & Tomasz Burzykowski (2013) \emph{Linear Mixed-Effects Models Using R: A Step-by-Step Approach} First Edition. Springer, New York. ISBN 978-1-4614-3899-1
}
}
