#' @import ggplot2
#' @import graphics
#' @import stats
#' @import lattice
#' @importFrom cowplot plot_grid
NULL

#' @title Plotting of tumor growth data from a fitted model
#' @description
#'  Vizualization of tumor growth data and linear mixed model fitted regression line for the fixed effects. This functions returns a [ggplot2] plot, allowing for
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
#' #' data(grwth_data)
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
  
  if (!is.na(drug_c)) {
    segment_data <- data.frame(x = rep(0,5), 
                               xend = model$dt1 %>% dplyr::group_by(.data$Treatment) %>% dplyr::summarise(Max = max(.data$Time)) %>% dplyr::select(.data$Max),
                               y = rep(0, 5), 
                               yend = nlme::fixef(model))
  } else {
    segment_data <- data.frame(x = rep(0,4), 
                               xend = model$dt1 %>% dplyr::group_by(.data$Treatment) %>% dplyr::summarise(Max = max(.data$Time)) %>% dplyr::select(.data$Max),
                               y = rep(0, 4), 
                               yend = nlme::fixef(model))
    }
  
  
  
  segment_data$yend <- segment_data$Max*segment_data$yend
  colnames(segment_data) <- c("x", "xend", "y", "yend")
  if (!is.na(drug_c)){
    segment_data$Treatment <- factor(x = c(trt_control, drug_a, drug_b, drug_c, combination), levels = c(trt_control, drug_a, drug_b, drug_c, combination))
  } else {
    segment_data$Treatment <- factor(x = c(trt_control, drug_a, drug_b, combination), levels = c(trt_control, drug_a, drug_b, combination))
  }
  hline <- data.frame(yintercept = 0)
  trt_col <- c("#3c3c3b", "#d50c52", "#00a49c", "#ff7f55","#601580")
  
  if (is.na(drug_c)) {
    trt_col <- trt_col[c(1:3,5)]
  }
  
  p <- model$dt1 %>% 
          ggplot(aes(.data$Time, .data$logRTV, color = .data$Treatment)) +
          geom_line(aes(group = .data$SampleID), alpha = 0.33) + geom_point(aes(group = .data$SampleID)) +
          ylab("Log (RTV)") + 
          xlab("Time since start of treatment") + 
          scale_x_continuous(breaks = unique(model$dt1$Time)) + 
          cowplot::theme_cowplot() + facet_wrap(~Treatment) +
          geom_segment(data = segment_data, 
                       aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend), 
                       lwd = 1.25, alpha = 0.75) + 
          geom_hline(data = hline, aes(yintercept = .data$yintercept), linetype = "dashed") +
          scale_color_manual(values = trt_col)
  return(p)
}



#' @title Plots for random effects diagnostics
#' @description
#' Visualization of random effects diagnostics for a fitted linear mixed model of tumor growth data.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @returns A list with different plots for evaluating the normality and homoscedasticity of the random effects, including:
#' - A normal Q-Q plot of the random effects of the model.
#' - A normal Q-Q plot of the residuals by sample.
#' - Boxplots of the conditional residuals by sample.
#' - Dotplots of the pearson residuals vs fitted values by sample.
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
#' # Generate plots 
#' plot_ranefDiagnostics(lmm)
#' # Access to specific plots
#' plot_ranefDiagnostics(lmm)[[1]]
#' plot_ranefDiagnostics(lmm)[[2]]
#' @export
plot_ranefDiagnostics <- function(model){
  # Residuals to draw QQ line
  
  # Individual Plots
  p1 <- ggplot(nlme::ranef(model), aes(sample = nlme::ranef(model)$Time)) + stat_qq(col = "gray20") + stat_qq_line() +
    labs(title = "Normal Q-Q Plot of Random Effects") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + cowplot::theme_cowplot()+
    theme(plot.title = element_text(size = 10, hjust = 0.5), axis.title = element_text(size = 12))
  p2 <- qqnorm(model, ~resid(., type = "normalized")|SampleID, pch=20, cex = 0.5, col = "gray20",
               main = list("Normal Q-Q Plot of Normalized Residuals by Sample", cex = 0.8), 
               par.strip.text=list(col="black", cex=0.8), xlab = "Normalized Residuals", abline = c(0,1))
  p3 <- plot(model, SampleID ~ resid(., type = "response"), abline = 0, main = list("Raw Residuals by Subject", cex = 0.8),
             xlab = "Residuals")
  p4 <- plot(model, residuals(., type = "pearson") ~ fitted(.)|SampleID, id = 0.05, adj = -0.03, pch = 20, col = "slateblue4", cex=0.75,
             main = list("Pearson Residuals vs Fitted Values by Sample", cex =0.8),par.strip.text=list(col="black", cex=0.8), idLabels = ~Time,
             abline = 0, ylab = "Pearson residuals")
  # Arranged plots
  p5 <- cowplot::plot_grid(p1,p2,p3,p4, ncol = 2)
  return(list(p1,p2,p3,p4,p5))
}

#' @title Plots for residuals diagnostics
#' @description
#' Visualization of residuals diagnostics for a fitted linear mixed model of tumor growth data.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @returns A list with different plots for evaluating the normality and homoscedasticity of the random effects, including:
#' - A normal Q-Q plot of the normalized residuals of the model.
#' - A normal Q-Q plot of the normalized residuals of the model by Time.
#' - A normal Q-Q plot of the normalized residuals of the model by Treatment.
#' - A dotplot of pearson residuals vs fitted values.
#' - A dotplot of the pearson residuals by Time and Treatment.
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
#' # Generate plots 
#' plot_residDiagnostics(lmm)
#' # Access to specific plots
#' plot_residDiagnostics(lmm)[[1]]
#' plot_residDiagnostics(lmm)[[2]]
#' @export
plot_residDiagnostics <- function(model){
  # Individual Plots
  p1 <- qqnorm(model, ~resid(., type = "normalized"),pch = 20, main = "Q-Q Plot of Normalized Residuals",
               xlab = "Normalized residuals", abline = c(0,1))
  p2 <- qqnorm(model, ~resid(., type = "normalized")|Time,pch = 20, main = "Q-Q Plot of Normalized Residuals by Time",
               par.strip.text=list(col="black", cex=1),
               xlab = "Normalized residuals", abline = c(0,1))
  p3 <- qqnorm(model, ~resid(., type = "normalized")|Treatment, pch = 20, main = "Q-Q Plot of Normalized Residuals by Treatment",
               par.strip.text=list(col="black", cex=1),
               xlab = "Normalized residuals", abline = c(0,1))
  p4 <- plot(model,main = "Pearson Residuals vs Fitted Values", pch = 20, ylab = "Pearson residuals")
  p5 <- plot(model, resid(., type = "pearson") ~Time|Treatment, id = 0.05, pch=20,
             adj = -0.03, cex = 1, main = "Pearson Residuals per Time and Treatment", 
             par.strip.text=list(col="black", cex=1), ylab = "Pearson residuals", abline = 0)
  # Arranged plots
  p6 <- cowplot::plot_grid(p1,p2,p3,p4,p5, nrow = 3, ncol = 2, align = "hv")
  return(list(p1,p2,p3,p4,p5,p6))
}

#' @title Plots of Observed vs Predicted Values
#' @description
#' Visualization of observed vs predicted values by a fitted linear mixed model of tumor growth data.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param nrow Number of rows of the layout to organize the observed vs predicted plots.
#' @param ncol Number of columns of the layout to organize the observed vs predicted plots.
#' @returns A layout (arranged in `nrow` rows and `ncol` columns) of the observed and predicted values of \eqn{log}(relative tumor volume) vs Time for each SampleID (i.e. subject), 
#' with the actual measurements, the regression line for each SampleID, and the marginal, treatment-specific, 
#' regression line for each treatment group.
#' @examples
#' #' data(grwth_data)
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
#' # Obtain the plots
#' plot_ObsvsPred(lmm, nrow = 4, ncol = 8)    
#' 
#' @export
plot_ObsvsPred <- function(model, nrow = 4, ncol = 5){
  TV.df <- model$data
  aug.Pred <- nlme::augPred(model, primary = ~Time, level = 0:1, length.out = 2, minimum = 0)
  plot(aug.Pred, layout = c(ncol, nrow, 1), lty = c(1,2),
       key = list(lines = list(lty = c(1,2), col = c("slateblue", "orange")),
                  text = list(c("Marginal", "Subject-specific")),
                  columns = 2,
                  space="top"), 
       pch = 20, lwd = 1.5, main = "Observed and Predicted Values by Time",
       par.strip.text=list(col="black", cex=1),
       xlab = list("Time", cex = 1.2), ylab = list("log RTV", cex = 1.2))
}

#' @title Plotting synergy results
#' @description Visualization of synergy results obtained by [lmmSynergy()].  This functions returns a [ggplot2] plot, allowing for
#'  further personalization.
#' @param syn_data Object obtained by [lmmSynergy()] with the results of synergy calculation using linear mixed models.
#' @details
#' `plot_lmmSynergy` produces a [ggplot2] plot with the results of the synergy calculation. Each dot represents the estimated combination index
#' or synergy score, and the gray lines represent the 95% confidence intervals, for each day. Each dot is colored based on the \eqn{- \log_{10} (p-value)}, with
#' purple colors indicating a \eqn{-\log_{10} (p-value) < 1.3; (p-value > 0.05)}, and green colors indicating a \eqn{-\log_{10} (p-value) > 1.3; (p-value < 0.05)}.
#' @returns A ggplot2 plot (see [ggplot2::ggplot()] for more details) with the combination index (CI) and synergy score (SS)
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
#' # Adding ggplot2 elements
#' plot_lmmSynergy(lmmSyn) + 
#' ggplot2::labs(title = "Synergy Calculation for Bliss") + 
#' ggplot2::theme(legend.position = "top")  
#' 
#' @export
plot_lmmSynergy <- function(syn_data){
  syn_data <- syn_data$Synergy
  
  syn_data$Metric[syn_data$Metric == "CI"] <- "Combination Index"
  syn_data$Metric[syn_data$Metric == "SS"] <- "Synergy Score"
  
  
  CI <- syn_data %>% dplyr::filter(.data$Metric == "Combination Index") %>% ggplot(aes(x = .data$Time, y = .data$Estimate)) +
    geom_segment(aes(x= .data$Time, y = .data$lwr, yend = .data$upr), color = "gray60", lwd = 1) + cowplot::theme_cowplot() +
    geom_point(aes(fill  = -log10(.data$pval)), size = 5, shape = 23, color = "gray60") +
    scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3) +
    ylab("Combination Index") + scale_x_continuous(breaks = unique(syn_data$Time)) + 
    geom_hline(yintercept = 1, lty = "dashed") + facet_wrap(~Metric) + theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold"))
  
  SS <- syn_data %>% dplyr::filter(.data$Metric == "Synergy Score") %>% ggplot(aes(x = .data$Time, y = .data$Estimate)) +
    geom_segment(aes(x= .data$Time, y = .data$lwr, yend = .data$upr), color = "gray60", lwd = 1) + cowplot::theme_cowplot() +
    geom_point(aes(fill  = -log10(.data$pval)), size = 5, shape = 23, color = "gray60") +
    scale_fill_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3) +
    ylab("Synergy Score") + scale_x_continuous(breaks = unique(syn_data$Time)) + 
    geom_hline(yintercept = 0, lty = "dashed") + facet_wrap(~Metric) + theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold"))
  
  
  cowplot::plot_grid(CI, SS)
}

#' @title Helper function to plot exemplary data for power calculation
#' @description
#' `plot_exmpDt` plots the regression lines of any exemplary data produced for a priori power
#' calculation.
#' @param exmpDt Data frame with exemplary data obtained with [APrioriPwr()].
#' @param grwrControl Value for the label of the coefficient for Control treatment group tumor growth rate.
#' @param grwrA Value for the label of the coefficient for Drug A treatment group tumor growth rate.
#' @param grwrB Value for the label of the coefficient for Drug B treatment group tumor growth rate.
#' @param grwrComb Value for the label of the coefficient for Combination (Drug A + Drug B) treatment group tumor growth rate.
#' @param sd_ranef Value for the label of the random effects standard deviation of the model.
#' @param sgma Value for the label of the residuals standard deviation of the model.
#' @returns A ggplot2 plot (see [ggplot2::ggplot()] for more details) showing the regression lines corresponding 
#' to the fixed effects for each treatment of the exemplary data for power calculations.
#' @keywords internal
#' @noRd
.plot_exmpDt <- function(exmpDt, grwrControl = NULL, grwrA = NULL, grwrB=NULL, grwrComb=NULL, sd_ranef=NULL, sgma=NULL){
  # Ploting exemplary data
  
  selDt <- with(exmpDt,{
    lvls <- levels(Treatment)
    i <- match(lvls, Treatment)
    subj <- subject[i]
    subset(exmpDt, subject %in% subj)
  })
  
  selDt %>% ggplot(aes(x = .data$Time, y = .data$mA)) + geom_line(aes(colour = .data$Treatment), lwd = 2) + 
    labs(title = "Exemplary Data") + ylab("logRTV") + xlab("Time since treatment start") + cowplot::theme_cowplot() +
    annotate(geom = "text", x = max(selDt$Time)+3, y = max(selDt$mA), label = paste("GR Control=",grwrControl), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Time)+3, y = max(selDt$mA)*0.95, label = paste("GR Drug A=",grwrA), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Time)+3, y = max(selDt$mA)*0.9, label = paste("GR Drug B=",grwrB), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Time)+3, y = max(selDt$mA)*0.85, label = paste("GR Combination=", grwrComb), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Time)+3, y = max(selDt$mA)*0.8, label = paste("SD=",sd_ranef), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Time)+3, y = max(selDt$mA)*0.75, label = paste("Sigma=",sgma), hjust = -0.05, size = 4) +
    coord_cartesian(xlim = c(0, max(selDt$Time)+5), clip = "off") +
    scale_color_manual(values = c("#3c3c3b", "#d50c52", "#00a49c", "#601580"))
}


