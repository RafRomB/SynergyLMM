#' @importFrom ggplot2 aes annotate arrow geom_hline geom_point geom_segment ggplot labs scale_fill_gradient2 scale_x_continuous unit xlab ylab
#' @importFrom rlang .data
#' @import graphics
NULL

#' @title Plotting synergy results
#' @description Visualization of synergy results obtained by [lmmSynergy()].  This functions returns a [ggplot2] plot, allowing for
#'  further personalization.
#' @param syn_data Object obtained by [lmmSynergy()] with the results of synergy calculation using linear mixed models.
#' @details
#' `plot_lmmSynergy` produces a [ggplot2] plot with the results of the synergy calculation. Each dot represents the estimated combination index
#' or synergy score, and the gray lines represent the 95% confidence intervals, for each day. Each dot is colored based on the \eqn{- \log_{10} (p-value)}, with
#' purple colors indicating a \eqn{-\log_{10} (p-value) < 1.3; (p-value > 0.05)}, and green colors indicating a \eqn{-\log_{10} (p-value) > 1.3; (p-value < 0.05)}.
#' @returns A list with ggplot2 plots (see [ggplot2::ggplot()] for more details) with the combination index (CI) and synergy score (SS)
#' estimates, confidence intervals and p-values for the synergy calculation using linear mixed models.
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
#'   combination = "Combination"
#'   )
#' # Obtain synergy results
#' lmmSyn <- lmmSynergy(lmm)
#' # Plot synergy results
#' plot_lmmSynergy(lmmSyn)
#' # Accessing to the combination index plot
#' plot_lmmSynergy(lmmSyn)$CI
#' # Accessing to only synergy score plot
#' plot_lmmSynergy(lmmSyn)$SS
#' # Accessing to the grid of both plots side by side
#' plot_lmmSynergy(lmmSyn)$CI_SS
#' # Adding ggplot2 elements
#' plot_lmmSynergy(lmmSyn)$CI + 
#' ggplot2::labs(title = "Synergy Calculation for Bliss") + 
#' ggplot2::theme(legend.position = "top")  
#' 
#' @export
plot_lmmSynergy <- function(syn_data){
  syn_data <- syn_data$Synergy
  
  syn_data$Metric[syn_data$Metric == "CI"] <- "Combination Index"
  syn_data$Metric[syn_data$Metric == "SS"] <- "Synergy Score"
  
  Model <- unique(syn_data$Model)
  
  CI <- syn_data %>% dplyr::filter(.data$Metric == "Combination Index") %>% ggplot(aes(x = .data$Time, y = .data$Estimate)) +
    geom_segment(aes(x= .data$Time, y = .data$lwr, yend = .data$upr), color = "gray60", lwd = 1, 
                 arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) + cowplot::theme_cowplot() +
    geom_point(aes(fill  = -log10(.data$pval)), size = 5, shape = 23, color = "gray60") +
    scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3, na.value = "firebrick3") +
    ylab("Combination Index") + xlab("Time since start of treatment") + 
    scale_x_continuous(breaks = unique(syn_data$Time)) + 
    geom_hline(yintercept = 1, lty = "dashed") + #facet_wrap(~Metric) + theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold"))
    labs(title = paste(Model, "Combination Index", sep = " ")) + 
    annotate(geom = "text", x = (min(syn_data$Time)-(syn_data$Time[2]-syn_data$Time[1])), 
             y = 0.95, angle = 90, hjust = 1, label = "Synergy", fontface = "bold", color = "#1f78b4") +
    annotate(geom = "text", x =(min(syn_data$Time)-(syn_data$Time[2]-syn_data$Time[1])), 
             y = 1.05, angle = 90, hjust = 0, label = "Antagonism", fontface = "bold", color = "#c21d2f")
  
  SS <- syn_data %>% dplyr::filter(.data$Metric == "Synergy Score") %>% ggplot(aes(x = .data$Time, y = .data$Estimate)) +
    geom_segment(aes(x= .data$Time, y = .data$lwr, yend = .data$upr), color = "gray60", lwd = 1,
                 arrow = arrow(angle = 90, length = unit(0.01, "npc"),ends = "both")) + cowplot::theme_cowplot() +
    geom_point(aes(fill  = -log10(.data$pval)), size = 5, shape = 23, color = "gray60") +
    scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3, na.value = "firebrick3") +
    ylab("Synergy Score") + xlab("Time since start of treatment") +
    scale_x_continuous(breaks = unique(syn_data$Time)) + 
    geom_hline(yintercept = 0, lty = "dashed") + #facet_wrap(~Metric) + theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold"))
    labs(title = paste(Model, "Synergy Score", sep = " ")) + 
    annotate(geom = "text", x = (min(syn_data$Time)-(syn_data$Time[2]-syn_data$Time[1])), 
             y = 0.33, angle = 90, hjust = 0, label = "Synergy", fontface = "bold", color = "#1f78b4") +
    annotate(geom = "text", x = (min(syn_data$Time)-(syn_data$Time[2]-syn_data$Time[1])), 
             y = -0.33, angle = 90, hjust = 1, label = "Antagonism", fontface = "bold", color = "#c21d2f")
  
  if (sum(syn_data$pval == 0) > 1) {
    ndec <- nchar(strsplit(as.character(min(syn_data$pval[syn_data$pval!=0])), "\\.")[[1]][2])
    apx_p <- paste("p<",ndec/ndec*10^-(ndec), sep = "")
    CI <- CI +
      annotate(
        geom = "point",
        x = Inf, y = min(syn_data$lwr[syn_data$Metric == "Combination Index"]),  # Position the diamond outside the plot (adjust as needed)
        shape = 23, size = 4, fill = "firebrick3", color = "gray60"
      ) +
      annotate(
        geom = "text",
        x = Inf, y = min(syn_data$lwr[syn_data$Metric == "Combination Index"]),  # Place the text slightly below the diamond
        label = apx_p, hjust = -0.25, vjust = 0.25, 
      ) +
      coord_cartesian(clip = "off")
    SS <- SS +
      annotate(
        geom = "point",
        x = Inf, y = min(syn_data$lwr[syn_data$Metric == "Synergy Score"]),  # Position the diamond outside the plot (adjust as needed)
        shape = 23, size = 4, fill = "firebrick3", color = "gray60"
      ) +
      annotate(
        geom = "text",
        x = Inf, y = min(syn_data$lwr[syn_data$Metric == "Synergy Score"]),  # Place the text slightly below the diamond
        label = apx_p, hjust = -0.25, vjust = 0.25, 
      ) +
      coord_cartesian(clip = "off")
  }
  
  CI_SS <- cowplot::plot_grid(CI, SS)
  return(list(CI = CI, SS = SS, CI_SS = CI_SS))
}
