#' @import ggplot2
#' @import graphics
#' @import stats
#' @import lattice
#' @importFrom cowplot plot_grid
NULL


#' @title Vizualization of tumor growth data and linear mixed model fitted regression line.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param trt_control String indicating the name assigned to the 'Control' group.
#' @param drug_a String indicating the name assigned to the 'Drug A' group.
#' @param drug_b String indicating the name assigned to the 'Drug B' group.
#' @param drug_ab String indicating the name assigned to the Combination ('Drug A' + 'Drug B') group.
#' @returns A ggplot2 plot (see [ggplot2::ggplot()] for more details) showing the tumor growth data represented as log(relative tumor volume) versus time since treatment initiation. 
#' The regression lines corresponding to the fixed effects for each treatment group are also plotted.
#' @export
plot_lmmModel <- function(model, trt_control = "Control", drug_a = "Drug_A", drug_b = "Drug_B", drug_ab= "Drug_AB"){
  segment_data <- data.frame(x = rep(0,4), 
                             xend = model$dt1 %>% dplyr::group_by(.data$Treatment) %>% dplyr::summarise(Max = max(.data$Day)) %>% dplyr::select(.data$Max),
                             y = rep(0, 4), 
                             yend = nlme::fixef(model))
  
  segment_data$yend <- segment_data$Max*segment_data$yend
  colnames(segment_data) <- c("x", "xend", "y", "yend")
  segment_data$Treatment <- factor(x = c(trt_control, drug_a, drug_b, drug_ab), levels = c(trt_control, drug_a, drug_b, drug_ab))
  
  p <- model$dt1 %>% 
          ggplot(aes(.data$Day, .data$logRTV, color = .data$Treatment)) +
          geom_line(aes(group = .data$Mouse), alpha = 0.33) + geom_point(aes(group = .data$Mouse)) +
          ylab("Log (RTV)") + 
          xlab("Days since start of treatment") + 
          scale_x_continuous(breaks = unique(model$dt1$Day)) + 
          cowplot::theme_cowplot() + facet_wrap(~Treatment) +
          geom_segment(data = segment_data, 
                       aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend), 
                       lwd = 1.25, alpha = 0.75)
  return(p)
}


#' @title Visualization of random effects diagnostics for a fitted linear mixed model of tumor growth data.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @returns A list with different plots for evaluating the normality and homoscedasticity of the random effects, including:
#' - A normal Q-Q plot of the random effects of the model.
#' - A normal Q-Q plot of the residuals by mouse.
#' - Boxplots of the conditional residuals by mouse.
#' - Dotplots of the pearson residuals vs fitted values by mouse.
#' @export
plot_ranefDiagnostics <- function(model){
  # Individual Plots
  p1 <- ggplot(nlme::ranef(model), aes(sample = nlme::ranef(model)$Day)) + stat_qq(col = "gray20") + stat_qq_line() +
    labs(title = "Normal Q-Q Plot of Random Effects") + xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + cowplot::theme_cowplot()+
    theme(plot.title = element_text(size = 10, hjust = 0.5), axis.title = element_text(size = 12))
  p2 <- qqnorm(model, ~resid(.)|Mouse, pch=20, cex = 0.5, col = "gray20",
               main = list("Normal Q-Q Plot by Mouse", cex = 0.8), par.strip.text=list(col="black", cex=.66))
  p3 <- plot(model, Mouse ~ resid(., type = "response"), abline = 0, main = list("Conditional Residuals by Mouse", cex = 0.8))
  p4 <- plot(model, residuals(., type = "pearson") ~ fitted(.)|Mouse, id = 0.05, adj = -0.03, pch = 20, col = "slateblue4", cex=0.5,
             main = list("Pearson Residuals vs Fitted Values by Mouse", cex =0.8),par.strip.text=list(col="black", cex=.66))
  # Arranged plots
  p5 <- cowplot::plot_grid(p1,p2,p3,p4, ncol = 2)
  return(list(p1,p2,p3,p4,p5))
}

#' @title Visualization of residuals diagnostics for a fitted linear mixed model of tumor growth data.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @returns A list with different plots for evaluating the normality and homoscedasticity of the random effects, including:
#' - A normal Q-Q plot of the normalized residuals of the model.
#' - A normal Q-Q plot of the normalized residuals of the model by Day.
#' - A normal Q-Q plot of the normalized residuals of the model by Treatment.
#' - A dotplot of pearson residuals vs fitted values.
#' - A dotplot of the pearson residuals by Day and Treatment.
#' @export
plot_residDiagnostics <- function(model){
  # Individual Plots
  p1 <- qqnorm(model, ~resid(., type = "normalized"),pch = 20, main = "Q-Q Plot of Normalized Residuals")
  p2 <- qqnorm(model, ~resid(., type = "normalized")|Day,pch = 20, main = "Q-Q Plot of Normalized Residuals by Day",
               par.strip.text=list(col="black", cex=.5))
  p3 <- qqnorm(model, ~resid(., type = "normalized")|Treatment, pch = 20, main = "Q-Q Plot of Normalized Residuals by Treatment",
               par.strip.text=list(col="black", cex=.5))
  p4 <- plot(model,main = "Pearson Residuals vs Fitted Values", pch = 20)
  p5 <- plot(model, resid(., type = "pearson") ~Day|Treatment, id = 0.05, pch=20,
             adj = -0.03, cex = 0.6, main = "Pearson Residuals per Time and Treatment", 
             par.strip.text=list(col="black", cex=.5))
  # Arranged plots
  p6 <- cowplot::plot_grid(p1,p2,p3,p4,p5, nrow = 3, ncol = 2, align = "hv")
  return(list(p1,p2,p3,p4,p5,p6))
}

#' @title Visualization of observed vs predicted values by a fitted linear mixed model of tumor growth data.
#' @param model An object of class "lme" representing the linear mixed-effects model fitted by [`lmmModel()`].
#' @param nrow Number of rows of the layout to organize the observed vs predicted plots.
#' @param ncol Number of columns of the layout to organize the observed vs predicted plots.
plot_ObsvsPred <- function(model, nrow, ncol){
  TV.df <- model$data
  aug.Pred <- nlme::augPred(model, primary = ~Day, level = 0:1, length.out = 2)
  plot(aug.Pred, layout = c(ncol, nrow, 1), lty = c(1,2),
       key = list(lines = list(lty = c(1,2), col = c("slateblue", "orange")),
                  text = list(c("Marginal", "Mouse-specific")),
                  columns = 2,
                  space="top"), 
       pch = 20, lwd = 1.5, main = "Observed and Predicted Values by Day",
       par.strip.text=list(col="black", cex=.5))
}

#' @title Visualization of synergy results obtained by [lmmSynergy()].
#' @param syn_data Object obtained by [lmmSynergy()] with the results of synergy calculation using linear mixed models.
#' @returns A ggplot2 plot (see [ggplot2::ggplot()] for more details) with the combination index (CI) and synergy score (SS)
#' estimates, confidence intervals and p-values for the synergy calculation using linear mixed models.
#' @export
plot_lmmSynergy <- function(syn_data){
  syn_data <- syn_data$Synergy
  hline <- data.frame(Metric = c("CI", "SS"), yintercept = c(1,0))
  syn_data %>% ggplot(aes(x = .data$Day, y = .data$Estimate)) +
    geom_segment(aes(x= .data$Day, y = .data$lwr, yend = .data$upr), color = "gray70", lwd = 1) + cowplot::theme_cowplot() +
    geom_point(aes(colour  = -log10(.data$pval)), size = 5, shape = 18) +
    scale_color_gradient2(name = "-log10\np-value", low = "darkorchid4",mid = "gray90", high = "darkcyan",midpoint = 1.3) +
    ylab("Value") + scale_x_continuous(breaks = unique(syn_data$Day)) +
    facet_wrap(~Metric, scales = "free") + theme(strip.background = element_rect(fill = "cyan4"), strip.text = element_text(color = "white", face = "bold")) + 
    labs(title = paste("Combination Index and Synergy Score for", unique(syn_data$Model), "Synergy")) +
    geom_hline(data = hline, aes(yintercept = .data$yintercept), linetype = "dashed")
}

#' @title Plot of exemplary data for power calculation
#' @param exmpDt Data frame with exemplary data obtained with [APrioriPwr()].
#' @param grwrControl Value for the label of the coefficient for Control treatment group tumor growth rate.
#' @param grwrA Value for the label of the coefficient for Drug A treatment group tumor growth rate.
#' @param grwrB Value for the label of the coefficient for Drug B treatment group tumor growth rate.
#' @param grwrComb Value for the label of the coefficient for Combination (Drug A + Drug B) treatment group tumor growth rate.
#' @param sd_ranef Value for the label of the random effects standard deviation of the model.
#' @param sgma Value for the label of the residuals standard deviation of the model.
#' @returns A ggplot2 plot (see [ggplot2::ggplot()] for more details) showing the regression lines corresponding 
#' to the fixed effects for each treatment of the exemplary data for power calculations.
#' @export
plot_exmpDt <- function(exmpDt, grwrControl = NULL, grwrA = NULL, grwrB=NULL, grwrComb=NULL, sd_ranef=NULL, sgma=NULL){
  # Ploting exemplary data
  
  selDt <- with(exmpDt,{
    lvls <- levels(Treatment)
    i <- match(lvls, Treatment)
    subj <- subject[i]
    subset(exmpDt, subject %in% subj)
  })
  
  selDt %>% ggplot(aes(x = .data$Day, y = .data$mA)) + geom_line(aes(colour = .data$Treatment), lwd = 2) + 
    labs(title = "Exemplary Data") + ylab("logRTV") + xlab("Days") + cowplot::theme_cowplot() +
    annotate(geom = "text", x = max(selDt$Day)+3, y = max(selDt$mA), label = paste("GR Control=",grwrControl), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Day)+3, y = max(selDt$mA)*0.95, label = paste("GR Drug A=",grwrA), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Day)+3, y = max(selDt$mA)*0.9, label = paste("GR Drug B=",grwrB), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Day)+3, y = max(selDt$mA)*0.85, label = paste("GR Combination=", grwrComb), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Day)+3, y = max(selDt$mA)*0.8, label = paste("SD=",sd_ranef), hjust = -0.05, size = 4) +
    annotate(geom = "text", x = max(selDt$Day)+3, y = max(selDt$mA)*0.75, label = paste("Sigma=",sgma), hjust = -0.05, size = 4) +
    coord_cartesian(xlim = c(0, max(selDt$Day)+5), clip = "off")
}


