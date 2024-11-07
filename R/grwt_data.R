#' Example Tumor Growth Data
#' 
#' A long format data frame with example tumor growth data, generated with [`simulateTumorGrowth()`], 
#' representing a simulated _in vivo_ two-drugs combination experiment, involving 8 subjects per group, with 30 days of follow-up, taking measurement
#' every 3 days.
#' 
#' @format
#' A data frame with 352 rows and 4 columns:
#' \describe{
#'  \item{subject}{Subject IDs.}
#'  \item{Time}{Time for each measurement in an arbitrary unit, for example, days.}
#'  \item{Treatment}{Treatment group for each subject (Control, DrugA, DrugB, and Combination).}
#'  \item{TumorVolume}{Measurement representing the tumor volume in an arbitrary unit, for example, mm^3.}
#' }
#' @usage data(grwth_data)
"grwth_data"