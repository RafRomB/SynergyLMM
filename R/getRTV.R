#' @importFrom rlang .data
NULL

#' @title Helper function to calculate the relative tumor volume from an imput data frame of tumor growth
#' 
#' @description
#' `getRTV` is a helper function used inside [`lmmModel()`] to obtain a dataframe with a column _RTV_ corresponding
#' to the relative tumor volume to time `time_start`, and a column _logRTV_ with the logarithm of _RTV_. 
#' 
#' @param data Data frame with the tumor growth data. The input data frame columns have to have the following names:
#' - `SampleID`: Column with the identifiers for each sample.
#' - `Time`: Column with the time for each measurement.
#' - `TV`: Column with the tumor volume measurement.
#' @param time_start Numeric value indicating the time at which the treatment started.
#' @returns The function returns the original data frame of tumor growth data, with 3 additional columns, corresponding to:
#' - RTV: Relative tumor volume to the tumor volume at `start_time`.
#' - logRTV: Logarithm of RTV column.
#' - TV0: Tumor volume at `start_time`.
#' @examples
#' # Load example dataset
#' data("grwth_data")
#' # Change column names
#' colnames(grwth_data) <- c("SampleID", "Time", "Treatment", "TV")
#' # Calculate relative tumor volume
#' getRTV(data = grwth_data, time_start = 0)

#' 
#' @export
getRTV <- function(data, time_start) {
  TV.df <- data
  
  # df with the initial tumor volume.
  
  TV0 <- as.data.frame(
    TV.df %>%
      dplyr::filter(.data$Time == time_start) %>%
      dplyr::select(.data$SampleID, .data$TV)
  )
  
  # Create the vectors for the relative tumor volumes
  
  samples <- unique(TV.df$SampleID)
  
  RTV.df <- data.frame(SampleID = character(0), RTV = numeric(0))
  
  # Relative Tumor Volume
  
  for (i in samples) {
    if (i %in% TV0$SampleID) {
      rtv <- TV.df %>% dplyr::filter(.data$SampleID == i) %>% 
        dplyr::select(.data$SampleID, .data$TV)
      rtv$RTV <- rtv$TV / TV0[TV0$SampleID == i, "TV"]
      RTV.df <- rbind(RTV.df, rtv[, c("SampleID", "RTV")])
    } else {
      rtv <- TV.df %>% dplyr::filter(.data$SampleID == i) %>% 
        dplyr::select(.data$SampleID, .data$TV)
      rtv$RTV <- NA
      RTV.df <- rbind(RTV.df, rtv[, c("SampleID", "RTV")])
    }
  }
  
  TV.df$RTV <- RTV.df$RTV
  
  TV.df$logRTV <- log(TV.df$RTV)
  
  colnames(TV0) <- c("SampleID", "TV0")
  
  TV.df <- dplyr::left_join(TV.df, TV0, by = "SampleID")
  
  return(TV.df)
}
