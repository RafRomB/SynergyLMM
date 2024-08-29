# Tests for .simulateTumorGrowth ----

test_that(".simulateTumorGrowth returns a data frame with correct columns", {
  result <- .simulateTumorGrowth(
    npg = 5,
    timepoints = c(0, 3, 5, 10),
    initial_volume = 100,
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.04,
    sd = 0.1
  )
  
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 4)  # subject, Treatment, Day, TumorVolume
  expect_true(all(c("subject", "Treatment", "Day", "TumorVolume") %in% colnames(result)))
})

test_that(".simulateTumorGrowth returns correct number of rows for given npg and timepoints", {
  npg <- 5
  timepoints <- c(0, 3, 5, 10)
  result <- .simulateTumorGrowth(
    npg = npg,
    timepoints = timepoints,
    initial_volume = 100,
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.04,
    sd = 0.1
  )
  
  # 4 groups, `npg` subjects per group, `length(timepoints)` measurements per subject
  expect_equal(nrow(result), 4 * npg * length(timepoints))
})

test_that(".simulateTumorGrowth handles different growth rates correctly", {
  result <- .simulateTumorGrowth(
    npg = 5,
    timepoints = c(0, 3, 5, 10),
    initial_volume = 100,
    grwrControl = 0.2,
    grwrA = 0.15,
    grwrB = 0.10,
    grwrComb = 0.05,
    sd = 0.1
  )
  
  # Check if the mean tumor volume is higher for groups with higher growth rates at the last time point
  last_time <- max(result$Day)
  last_time <- dplyr::filter(result, Day == last_time)
  means <- dplyr::summarize(last_time, .by = Treatment, mean_volume = mean(TumorVolume))
  
  expect_true(means$mean_volume[means$Treatment == "DrugA"] < means$mean_volume[means$Treatment == "Control"])
  expect_true(means$mean_volume[means$Treatment == "DrugB"] < means$mean_volume[means$Treatment == "DrugA"])
  expect_true(means$mean_volume[means$Treatment == "Combination"] < means$mean_volume[means$Treatment == "DrugB"])
})

test_that(".simulateTumorGrowth handles different standard deviations correctly", {
  result_sd_low <- .simulateTumorGrowth(
    npg = 5,
    timepoints = c(0, 3, 5, 10),
    initial_volume = 100,
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.04,
    sd = 0.01
  )
  
  result_sd_high <- .simulateTumorGrowth(
    npg = 5,
    timepoints = c(0, 3, 5, 10),
    initial_volume = 100,
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.04,
    sd = 0.5
  )
  
  # The variance of tumor volume should be higher with a higher sd
  var_low <- var(result_sd_low$TumorVolume)
  var_high <- var(result_sd_high$TumorVolume)
  
  expect_true(var_high > var_low)
})

test_that(".simulateTumorGrowth handles single timepoint correctly", {
  result <- .simulateTumorGrowth(
    npg = 5,
    timepoints = c(0),
    initial_volume = 100,
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.04,
    sd = 0.1
  )
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4 * 5)  # 4 groups, 5 subjects per group, 1 timepoint
})
