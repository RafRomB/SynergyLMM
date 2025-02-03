#' @import nlme
#' @importFrom rlang .data
NULL

#' @title Linear Mixed Effect Model for Tumor Growth
#' 
#' @description
#' `lmmModel()` fits a linear mixed effect model from a tumor growth dataset. The input data frame must be in long format and include at least the following columns: column with the sample ids,
#' column with the time at which each measurement has been done, a column indicating the treatment group, and a column with the tumor measurement (e.g. tumor volume).
#' @param data A data frame with the tumor growth data, in long format.
#' It must contain at least the following columns: mice IDs, time of follow-up (numeric number), treatment and tumor volume (numeric number).
#' @param sample_id String indicating the name of the column with the mice IDs.
#' @param time String indicating the name of the column with the time of follow-up.
#' @param treatment String indicating the name of the column with the treatment corresponding to each sample.
#' @param tumor_vol String indicating the name of the column with the tumor volume (or any other measurement representing the tumor growth).
#' @param trt_control String indicating the name assigned to the 'Control' group.
#' @param drug_a String indicating the name assigned to the 'Drug A' group.
#' @param drug_b String indicating the name assigned to the 'Drug B' group.
#' @param drug_c String indicating the name assigned to the 'Drug C' group (if present).
#' @param combination String indicating the name assigned to the Combination ('Drug A' + 'Drug B', or 'Drug A' + 'Drug B' + 'Drug C') group.
#' @param time_start Numeric value indicating the time point at which the treatment started. If not
#' specified, the minimum value in the `time` column is used as the starting time point.
#' @param time_end Numeric value indicating the last time point to be included in the analysis. If not
#' specified, the maximum value in the `time` column is used as the final time point.
#' @param min_observations Minimum number of observation for each sample to be included in the analysis. 
#' @param show_plot Logical indicating if a plot for the log of the relative tumor volume (RTV) vs Time for each sample, 
#' and the model calculated marginal slope for each treatment, should be produced.
#' @param tum_vol_0 Select the behavior of the function regarding measurements in which the tumor measurement is 0, and therefore the logarithmic 
#' transformation is not possible. Possible options are 'ignore' (default), to ignore these measurements, or 'transform', to add 1 unit to all
#' measurements before the log transformation.
#' @param ... Additional arguments to be passed to \code{nlme::\link[nlme:lme]{lme}}.
#' 
#' @details
#' 
#' `lmmModel()` relies in the assumption that tumor growth follows an exponential kinetics. Any departure from this assumption can be tested using the diagnostics functions [`ranefDiagnostics()`],
#'  [`residDiagnostics()`], and [`ObsvsPred()`].
#' 
#' The model formula for the fitted model is:
#' \deqn{\log RTV_{i}(t) = \beta_{Trt_i} \cdot t + b_i \cdot t + \varepsilon_{i} (t).}
#' 
#' The term \eqn{\log RTV_{i}(t)} denotes the value of the logarithm of the relative tumor volume measured for subject \eqn{i} at time \eqn{t}. In the fixed-effect part of the model, 
#' \eqn{\beta_{Trt_i}} represents the coefficients corresponding to the specific growth rate for each treatment \eqn{Trt_i}, where \eqn{Trt_i \in \{control, A, B, combination\}} in the case of
#' two-drugs combination experiments, or \eqn{Trt_i \in \{control, A, B, C, combination\}} in the case of three-drugs combination experiments. That means that \eqn{\beta_{control}} corresponds
#' to the control group, \eqn{\beta_{A}} to Drug A, \eqn{\beta_{B}} to Drug B, \eqn{\beta_{C}} to Drug C (if present), and \eqn{\beta_{combination}} to the combination of Drug A + Drug B (or
#' Drug A + Drug B + Drug C, in the case of a three-drugs experiment). 
#' 
#' \eqn{t} is the time of measurement. 
#' 
#' In the random-effects part of the model, \eqn{b_i} is the random slope associated with time representing the subject-specific growth-rate. 
#' Finally, \eqn{\varepsilon_{i}(t)} is the residual random error.
#' 
#' The implementation of the linear mixed model in `lmmModel()` is done using \code{nlme::\link[nlme:lme]{lme}}, which also allows for the 
#' specification of within-group correlations structures and/or unequal variances. These, and additional parameters,
#' can be passed to the \code{nlme::\link[nlme:lme]{lme}} function through the `...` argument for fitting the model (see examples below).
#' 
#' 
#' @return An object of class "lme" (see \code{nlme::\link[nlme:lme]{lme}} for details) representing the linear mixed-effects model fit. If `show_plot = TRUE`, the plot
#' of the tumor growth data obtained with [plot_lmmModel()] is also shown. 
#' @references
#' - Pinheiro JC, Bates DM (2000). _Mixed-Effects Models in S and S-PLUS_. Springer, New York. doi:10.1007/b98882 \doi{doi:10.1007/b98882}.
#' - Pinheiro J, Bates D, R Core Team (2024). _nlme: Linear and Nonlinear Mixed Effects Models_. R package version 3.1-166, <https://CRAN.R-project.org/package=nlme>.
#' - Andrzej Galecki & Tomasz Burzykowski (2013) _Linear Mixed-Effects Models Using R: A Step-by-Step Approach_ First Edition. Springer, New York. ISBN 978-1-4614-3899-1
#' 
#' @examples
#' data("grwth_data")
#' # Most simple model
#' lmmModel(
#'  data = grwth_data,
#'  sample_id = "subject",
#'  time = "Time",
#'  treatment = "Treatment",
#'  tumor_vol = "TumorVolume",
#'  trt_control = "Control",
#'  drug_a = "DrugA",
#'  drug_b = "DrugB",
#'  combination = "Combination"
#'  )

#'# Changing the last time point of follow-up
#'lmmModel(
#'  data = grwth_data,
#'  sample_id = "subject",
#'  time = "Time",
#'  treatment = "Treatment",
#'  tumor_vol = "TumorVolume",
#'  trt_control = "Control",
#'  drug_a = "DrugA",
#'  drug_b = "DrugB",
#'  combination = "Combination",
#'  time_end = 21
#'  )

#'# Adding additional parameters for model fitting
#'lmmModel(
#'  data = grwth_data,
#'  sample_id = "subject",
#'  time = "Time",
#'  treatment = "Treatment",
#'  tumor_vol = "TumorVolume",
#'  trt_control = "Control",
#'  drug_a = "DrugA",
#'  drug_b = "DrugB",
#'  combination = "Combination",
#'  # Adding variance function to represent a different variance per subject
#'  weights = nlme::varIdent(form = ~1|SampleID),
#'  # Specifiying control values for lme Fit (useful when convergence problems appear)
#'  control = nlme::lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 100, msMaxEval = 1000)
#'  )
#' 
#' # Fit a model specifying a different variance per Time
#'lmmModel(
#'  data = grwth_data,
#'  sample_id = "subject",
#'  time = "Time",
#'  treatment = "Treatment",
#'  tumor_vol = "TumorVolume",
#'  trt_control = "Control",
#'  drug_a = "DrugA",
#'  drug_b = "DrugB",
#'  combination = "Combination",
#'  # Adding variance function to represent a different variance per Time
#'  weights = nlme::varIdent(form = ~1|Time)
#'  )
#' 
#' 
#' 
#' @export
lmmModel <- function(data,
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
                     ...) {
  
  
  # Check if required columns are present
  tryCatch({
    required_columns <- c(sample_id, time, treatment, tumor_vol)
    missing_columns <- setdiff(required_columns, names(data))
    if (length(missing_columns) > 0) {
      stop(
        "The following required columns are missing from the data: ",
        paste(missing_columns, collapse = ", ")
      )
    }
  }, error = function(e) {
    stop("Error in column checking: ", e$message)
  })
  
  # Check if there are exactly 4 or 5 treatments in the treatment column
  
  tryCatch({
    if (is.na(drug_c)) {
      expected_treatments <- c(trt_control, drug_a, drug_b, combination)
      actual_treatments <- unique(data[[treatment]])
    } else {
      expected_treatments <- c(trt_control, drug_a, drug_b, drug_c, combination)
      actual_treatments <- unique(data[[treatment]])
    }
    
    missing_treatments <- setdiff(expected_treatments, actual_treatments)
    
    if (length(missing_treatments) > 0) {
      stop(
        "The treatment column is missing expected treatments: ",
        paste(missing_treatments, collapse = ", ")
      )
    }
    
    unrecognized_treatments <- setdiff(actual_treatments, expected_treatments)
    
    if (length(unrecognized_treatments) > 0) {
      stop(
        "The treatment column contains unrecognized treatments: ",
        paste(unrecognized_treatments, collapse = ", ")
      )
    }
    
  }, error = function(e) {
    stop("Error in treatment checking: ", e$message)
  })
  
  tryCatch({
    if (!is.numeric(min_observations) || length(min_observations) != 1 || min_observations <= 0) {
      stop("The `min_observations` parameter must be a positive numeric value.")
    }
  }, error = function(e){
    stop(e$message)
  })
  
  col.names <- c(sample_id, time, treatment, tumor_vol)
  TV.df <- data %>% dplyr::select(dplyr::all_of(col.names))
  colnames(TV.df) <- c("SampleID", "Time", "Treatment", "TV")
  
  if (!is.na(drug_c)){
    TV.df$Treatment <- factor(TV.df$Treatment,
                              levels = c(trt_control, drug_a, drug_b, drug_c, combination))
  } else {
    TV.df$Treatment <- factor(TV.df$Treatment,
                              levels = c(trt_control, drug_a, drug_b, combination))
  }
  
  
  TV.df$Time <- as.numeric(TV.df$Time)
  
  TV.df$TV <- as.numeric(TV.df$TV)
  
  # First, we will remove those rows for which we don't have the total volume, and
  # we will use only the data after the treatment start
  
  # If time of treatment start is not defined, use the minimum value:
  
  if(is.null(time_start)) {
    time_start <- min(TV.df$Time)
  }
  
  TV.df <- TV.df %>% dplyr::filter(.data$Time >= time_start)
  
  # Remove those samples for which TV0 == 0
  # (and therefore, no RTV can be calculated)
  
  samples0 <- TV.df %>% dplyr::filter(.data$Time == time_start & .data$TV == 0) %>% dplyr::select(.data$SampleID)
  
  if (length(samples0$SampleID) > 0) {
    warning(paste(paste(samples0$SampleID, collapse = ","), 
                  " subjects have measurements with value 0 at the initial time point: ",time_start,". ",
                  "These subjects will be removed from the analysis.", sep = ""))
  }
  
  samples <- TV.df %>% dplyr::filter(.data$Time == time_start & .data$TV != 0) %>% dplyr::select(.data$SampleID)
  
  TV.df <- TV.df %>% dplyr::filter(.data$SampleID %in% samples$SampleID)
  
  
  # Handling cases in which the tumor volume is 0
  
  if (sum(TV.df$TV == 0, na.rm = TRUE) > 0) {
    warning("Some tumor measurements are 0.")
    if (tum_vol_0 == "ignore") {
      warning("Tumor measurements that are 0 will be ignored to avoid obtaining Inf values when taking log.",
              "\nSet 'tum_vol_0' argument to 'transform' to change this behavior.")
      TV.df$TV[TV.df$TV == 0] <- NA
    } else if (tum_vol_0 == "transform") {
      warning("All measurements will be transformed adding 1 unit to avoid obtaining Inf values when taking log.",
              "\nSet 'tum_vol_0' argument to 'ignore' to change this behavior.")
      TV.df$TV <- TV.df$TV + 1
    }
  }
  
  # Remove mesurements in which tumor measurement does not have a value assigned
  
  TV.df <- TV.df %>% dplyr::filter(!is.na(.data$TV))
  
  # Filter data until the maximum day specified in time_end
  if(!is.null(time_end)){
    TV.df <- TV.df %>% dplyr::filter(.data$Time <= time_end)
  }
  
  # Remove samples with less than the minimum of observations specified
  
  samples <- TV.df %>% dplyr::count(.data$SampleID, .by = .data$SampleID) %>%
    dplyr::filter(.data$n >= min_observations) %>% dplyr::select(.data$SampleID)
  
  TV.df <- TV.df %>% dplyr::filter(.data$SampleID %in% samples$SampleID)
  
  # Calculate the relative tumor volume
  TV.df <- getRTV(TV.df, time_start)
  
  # Convert SampleID to factor
  
  TV.df$SampleID <- as.factor(TV.df$SampleID)
  
  # Data frame for visualization
  TV.plot <- TV.df
  
  TV.df <- TV.df %>% dplyr::filter(.data$Time != time_start)
  
  TV.df$Time <- TV.df$Time - time_start
  TV.plot$Time <- TV.plot$Time - time_start
  
  TV.df <- as.data.frame(TV.df)
  
  args <- list(...)
  
  model <- do.call(nlme::lme, c(
    list(
      fixed = logRTV ~  0 + Time:Treatment,
      random = ~ 0 + Time | SampleID,
      data = TV.df
    ),
    args
  ))
  
  model$dt1 <- TV.plot
  
  if (show_plot) {
    plot(
      plot_lmmModel(
        model = model,
        trt_control = trt_control,
        drug_a = drug_a,
        drug_b = drug_b,
        drug_c = drug_c,
        combination = combination
      )
    )
  }
  attr(model, "SynergyLMM") <- "lmmModel"
  return(model)
}
