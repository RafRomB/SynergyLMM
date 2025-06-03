#' @importFrom ggplot2 aes facet_wrap geom_hline geom_line geom_point geom_segment ggplot scale_color_manual scale_x_continuous xlab ylab
#' @importFrom stats na.omit
#' @importFrom cowplot plot_grid
#' @importFrom rlang .data
NULL

#' @title Plotting of tumor growth data from a fitted model
#' @description
#'  Vizualization of tumor growth data and linear mixed model fitted regression line for the fixed effects. This functions returns a [ggplot2](https://ggplot2.tidyverse.org/) plot, allowing for
#'  further personalization.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param trt_control String indicating the name assigned to the 'Control' group.
#' @param drug_a String indicating the name assigned to the 'Drug A' group.
#' @param drug_b String indicating the name assigned to the 'Drug B' group.
#' @param drug_c String indicating the name assigned to the 'Drug C' group (if present).
#' @param combination String indicating the name assigned to the Combination ('Drug A' + 'Drug B', or 'Drug A' + 'Drug B' + 'Drug C') group.
#' @returns A ggplot2 plot (see [ggplot2::ggplot()] for more details) showing the tumor growth data represented as log(relative tumor volume) versus time since treatment initiation. 
#' The regression lines corresponding to the fixed effects for each treatment group are also plotted.
#' @examples
#' data(grwth_data)
#' # Fit the model
#' lmm <- lmmModel(
#'   data = grwth_data,
#'   sample_id = "subject",
#'   time = "Time",
#'   treatment = "Treatment",
#'   tumor_vol = "TumorVolume",
#'   trt_control = "Control",
#'   drug_a = "DrugA",
#'   drug_b = "DrugB",
#'   combination = "Combination",
#'   show_plot = FALSE
#'   )
#' # Default plot
#' plot_lmmModel(lmm,
#' trt_control = "Control",
#' drug_a = "DrugA",
#' drug_b = "DrugB",
#' combination = "Combination"
#' )
#' # Adding ggplot2 elements
#' plot_lmmModel(lmm,
#' trt_control = "Control",
#' drug_a = "DrugA",
#' drug_b = "DrugB",
#' combination = "Combination"
#' ) + ggplot2::labs(title = "Example Plot") + ggplot2::theme(legend.position = "top")
#' 
#' @export
plot_lmmModel <- function(model,
                          trt_control = "Control",
                          drug_a = "Drug_A",
                          drug_b = "Drug_B",
                          drug_c = NA,
                          combination = "Combination") {
  
  if (sum(!na.omit(c(trt_control, drug_a, drug_b, drug_c, combination)) %in% unique(model$dt1$Treatment)) > 0) {
    stop("Treatment group names provided do not coincide with treatment group names in the model. Please, provide the correct treatment group names ",
         "in the arguments.")
  }
  
  hline <- data.frame(yintercept = 0)
  trt_col <- c("#3c3c3b", "#d50c52", "#00a49c", "#ff7f55","#601580")
  
  if (is.na(drug_c)) {
    trt_col <- trt_col[c(1:3,5)]
  }
  
  plot_df <- model$data %>%
    dplyr::mutate(fixed_effect = predict(model, level = 0))  # level = 0 = fixed only
  
  # Add manually data at time zero
  initial_points <- model$data %>%
    dplyr::distinct(.data$SampleID, .data$Treatment) %>%
    dplyr::mutate(Time = 0, logRTV = 0)
  initial_points <- initial_points %>%
    dplyr::mutate(fixed_effect = 0)
  
  plot_df <- dplyr::bind_rows(plot_df, initial_points)
  
  
  p <- plot_df %>% 
    ggplot(aes(.data$Time, .data$logRTV, color = .data$Treatment)) +
    geom_line(aes(group = .data$SampleID), alpha = 0.33) + geom_point(aes(group = .data$SampleID)) +
    geom_line(aes(y = .data$fixed_effect, group = .data$SampleID), lwd = 1, alpha = 0.8) +
    ylab("log (RTV)") + 
    xlab("Time since start of treatment") + 
    scale_x_continuous(breaks = unique(model$dt1$Time)) + 
    cowplot::theme_cowplot() + facet_wrap(~Treatment) +
    geom_hline(data = hline, aes(yintercept = .data$yintercept), linetype = "dashed") +
    scale_color_manual(values = trt_col)
  return(p)
}

