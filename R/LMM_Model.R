#' @import nlme
NULL

#' @title Helper function to calculate the relative tumor volume from an imput data frame of tumor growth
#' @param data Data frame with the tumor growth data.
#' @param day_start Numeric value indicating the day at which the treatment started.
#' @returns The function returns the original data frame of tumor growth data, with 3 additional columns, corresponding to:
#' - RTV: Relative tumor volume to the tumor volume at `start_day`.
#' - logRTV: Logarithm of RTV column.
#' - TV0: Tumor volume at `start_day`.
#' @export
.getRTV <- function(data, day_start){
  
  TV.df <- data
  
  # df with the initial tumor volume.
  
  TV0 <- as.data.frame(TV.df %>% 
                         dplyr::filter(.data$Day == day_start) %>% 
                         dplyr::select(.data$Mouse, .data$TV))
  
  # Create the vectors for the relative tumor volumes
  
  samples <- unique(TV.df$Mouse)
  
  RTV.df <- data.frame(Mouse = character(0), RTV = numeric(0))
  
  # Relative Tumor Volume
  
  for (i in samples) {
    if (i %in% TV0$Mouse) {
      rtv <- TV.df %>% dplyr::filter(.data$Mouse == i) %>% dplyr::select(.data$Mouse, .data$TV)
      rtv$RTV <- rtv$TV / TV0[TV0$Mouse == i, "TV"]
      RTV.df <- rbind(RTV.df, rtv[, c("Mouse", "RTV")])
    } else {
      rtv <- TV.df %>% dplyr::filter(.data$Mouse == i) %>% dplyr::select(.data$Mouse, .data$TV)
      rtv$RTV <- NA
      RTV.df <- rbind(RTV.df, rtv[, c("Mouse", "RTV")])
    }
  }
  
  TV.df$RTV <- RTV.df$RTV
  
  TV.df$logRTV <- log(TV.df$RTV)
  
  colnames(TV0) <- c("Mouse","TV0")
  
  TV.df <- dplyr::left_join(TV.df, TV0, by = "Mouse")
  return(TV.df)
}

#' @title Generate Linear Mixed Model for synergy calculation
#' 
#' @param data A data frame with the tumor growth data, in long format.
#' It must contain at least contain the following columns: mice IDs, days of follow-up (numeric number), treatment and tumor volume (numeric number).
#' @param mouse_id String indicating the name of the column with the mice IDs.
#' @param day String indicating the name of the column with the days of follow-up.
#' @param treatment String indicating the name of the column with the treatment corresponding to each mouse.
#' @param tumor_vol String indicating the name of the column with the tumor volume (or any other measurement representing the tumor growth).
#' @param trt_control String indicating the name assigned to the 'Control' group.
#' @param drug_a String indicating the name assigned to the 'Drug A' group.
#' @param drug_b String indicating the name assigned to the 'Drug B' group.
#' @param drug_ab String indicating the name assigned to the Combination ('Drug A' + 'Drug B') group.
#' @param day_start Numeric value indicating the day at which the treatment started.
#' @param min_observations Minimum number of observation for each mouse to be included in the analysis.
#' @param show_plot Logical indicating if a plot for the log of the relative tumor volume (RTV) vs Day for each mouse, 
#' and the model calculated marginal slope for each treatment, should be produced.
#' @param ... Additional arguments to be passed to \code{nlme::\link[nlme:lme]{lme}}.
#' @return An object of class "lme" (see \code{nlme::\link[nlme:lme]{lme}} for details) representing the linear mixed-effects model fit. 
#' @export
lmm_model <- function(data,
                      mouse_id = "Mouse",
                      day = "Day",
                      treatment = "Treatment",
                      tumor_vol = "TV",
                      trt_control = "Control",
                      drug_a = "Drug_A",
                      drug_b = "Drug_B",
                      drug_ab = "Drug_AB",
                      day_start = 0,
                      min_observations = 1,
                      show_plot = TRUE,
                      ...) {
  
  col.names <- c(mouse_id, day, treatment, tumor_vol)
  TV.df <- data %>% dplyr::select(dplyr::all_of(col.names))
  colnames(TV.df) <- c("Mouse", "Day", "Treatment", "TV")
  
  TV.df$Treatment <- factor(TV.df$Treatment,
                            levels = c(trt_control, drug_a, drug_b, drug_ab))
  
  TV.df$Day <- as.numeric(TV.df$Day)
  
  TV.df$TV <- as.numeric(TV.df$TV)
  
  # First, we will remove those rows for which we don't have the total volume, and
  # we will use only the data after the treatment start
  
  TV.df <- TV.df %>% dplyr::filter(.data$Day >= day_start & !is.na(.data$TV))
  
  # Calculate the relative tumor volume
  TV.df <- .getRTV(TV.df, day_start)
  
  # Data frame for visualization
  TV.plot <- TV.df
  
  TV.df <- TV.df %>% dplyr::filter(.data$Day != day_start)
  
  TV.df$Day <- TV.df$Day - day_start
  TV.plot$Day <- TV.plot$Day - day_start
  
  TV.df <- as.data.frame(TV.df)
  
  args <- list(...)
  
  model <- do.call(nlme::lme, c(
    list(
      fixed = logRTV ~  0 + Day:Treatment,
      random = ~ 0 + Day | Mouse,
      data = TV.df
    ),
    args
  ))
  
  model$dt1 <- TV.plot
  
  if (show_plot) {
    print(
      Plot_LMM_Model(
        model = model,
        trt_control = trt_control,
        drug_a = drug_a,
        drug_b = drug_b,
        drug_ab = drug_ab
      )
    )
  }
  
  return(model)
}

