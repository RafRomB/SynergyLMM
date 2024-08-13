
#' @importFrom ggplot2 aes geom_line geom_point ylab xlab scale_x_continuous facet_wrap geom_segment
#' @export
LMM_Model <- function(data, Mouse, Day, Treatment, TV, C, A, B, AB, day_start = 0, min_obs = 1, plot = TRUE, ...){
 col.names <- c(Mouse, Day, Treatment, TV)
 TV.df <- data %>% dplyr::select(dplyr::all_of(col.names))
 colnames(TV.df) <- c("Mouse", "Day", "Treatment", "TV")
 
 TV.df$Treatment <- factor(TV.df$Treatment, levels = c(C, A, B, AB))
 
 TV.df$Day <- as.numeric(TV.df$Day)
 
 TV.df$TV <- as.numeric(TV.df$TV) 
 
 # First, we will remove those rows for which we don't have the total volume, and
 # we will use only the data after the treatment start
 
 TV.df <- TV.df %>% dplyr::filter(Day >= day_start & !is.na(TV))
 
 # Remove samples with less than the minimum of observations specified
 
 samples <- TV.df %>% dplyr::count(Mouse, .by = Mouse) %>% 
   dplyr::filter(n > min_obs) %>% dplyr::select(Mouse)
 
 TV.df <- TV.df %>% dplyr::filter(Mouse %in% samples$Mouse)
 
 TV.df$Mouse <- as.factor(TV.df$Mouse)
 
 # df with the initial tumor volume.
 
 TV0 <- as.data.frame(TV.df %>% dplyr::filter(Day == day_start) %>% dplyr::select(Mouse, TV))
 
 # Create the vectors for the relative tumor volumes
 
 samples <- unique(TV.df$Mouse)
 
 RTV.df <- data.frame(Mouse = character(0), RTV = numeric(0))
 
 # Relative Tumor Volume
 
 for (i in samples) {
   if (i %in% TV0$Mouse) {
     rtv <- TV.df %>% dplyr::filter(Mouse == i) %>% dplyr::select(Mouse, TV) 
     rtv$RTV <- rtv$TV / TV0[TV0$Mouse == i, "TV"]
     RTV.df <- rbind(RTV.df, rtv[,c("Mouse","RTV")])  
   } else {
     rtv <- TV.df %>% dplyr::filter(Mouse == i) %>% dplyr::select(Mouse, TV) 
     rtv$RTV <- NA
     RTV.df <- rbind(RTV.df, rtv[,c("Mouse","RTV")])
   }
 }
 
 TV.df$RTV <- RTV.df$RTV
 
 TV.df$logRTV <- log(TV.df$RTV)
 
 colnames(TV0) <- c("Mouse","TV0")
 
 TV.df <- dplyr::left_join(TV.df, TV0, by = "Mouse")
 
 TV.plot <- TV.df
 
 TV.df <- TV.df %>% dplyr::filter(Day != day_start)
 
 TV.df$Day <- TV.df$Day - day_start
 TV.plot$Day <- TV.plot$Day - day_start
 
 TV.df <- as.data.frame(TV.df)
 
 args <- substitute(list(...))

 model <- nlme::lme(fixed = logRTV ~  0 + Day:Treatment, random = ~0 + Day|Mouse, data = TV.df, ...)
 
 if("weights" %in% names(args)){
   model$call$weights <- args$weights
 }
 if("control" %in% names(args)){
   model$call$control <- args$control
 }
 
 if(plot){
   
   segment_data <- data.frame(x = rep(0,4), 
                              xend = TV.plot %>% dplyr::group_by(Treatment) %>% dplyr::summarise(Max = max(Day)) %>% dplyr::select(Max),
                              y = rep(0, 4), 
                              yend = nlme::fixef(model))
   
   segment_data$yend <- segment_data$Max*segment_data$yend
   colnames(segment_data) <- c("x", "xend", "y", "yend")
   segment_data$Treatment <- factor(x = c(C, A, B, AB), levels = c(C, A, B, AB))
   
   print(TV.plot %>% 
     ggplot2::ggplot(aes(Day, logRTV, color = Treatment)) +
     geom_line(aes(group = Mouse), alpha = 0.33)+ geom_point(aes(group = Mouse)) +
     ylab("Log (RTV)") + 
     xlab("Days since start of treatment") + 
     scale_x_continuous(breaks = unique(TV.plot$Day)) + 
     cowplot::theme_cowplot() + facet_wrap(~Treatment) +
     geom_segment(data = segment_data, 
                  aes(x = x, xend = xend, y = y, yend = yend), 
                  lwd = 1.25, alpha = 0.75))
 }
 
 return(model)
}

