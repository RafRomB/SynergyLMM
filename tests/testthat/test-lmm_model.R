# Dummy dataset for testing
set.seed(123)
test_data <- data.frame(
  Mouse = rep(1:5, each = 5),
  Day = rep(0:4, times = 5),
  TV = c(100, 120, 140, 160, 180,   # Mouse 1
         80, 100, 120, 150, 190,    # Mouse 2
         90, 95, 100, 110, 130,     # Mouse 3
         110, 115, 120, 125, 130,   # Mouse 4
         70, 75, 85, 90, 100)       # Mouse 5
)

# Test if the function correctly calculates RTV and logRTV
test_that(".getRTV correctly calculates RTV and logRTV", {
  result <- .getRTV(test_data, day_start = 0)
  
  # Expected RTV for each Mouse at each Day
  expected_RTV <- c(1, 1.2, 1.4, 1.6, 1.8,   # Mouse 1
                    1, 1.25, 1.5, 1.875, 2.375,  # Mouse 2
                    1, 1.055556, 1.111111, 1.222222, 1.444444, # Mouse 3
                    1, 1.045455, 1.090909, 1.136364, 1.181818, # Mouse 4
                    1, 1.071429, 1.214286, 1.285714, 1.428571) # Mouse 5
  
  # Check if RTV is correctly calculated
  expect_equal(result$RTV, expected_RTV, tolerance = 1e-5)
  
  # Check if logRTV is correctly calculated
  expect_equal(result$logRTV, log(expected_RTV), tolerance = 1e-5)
})

# Test if the function correctly adds the TV0 column
test_that(".getRTV correctly adds TV0 column", {
  result <- .getRTV(test_data, day_start = 0)
  
  # Expected TV0 for each Mouse
  expected_TV0 <- rep(c(100, 80, 90, 110, 70), each = 5)
  
  # Check if TV0 is correctly added
  expect_equal(result$TV0, expected_TV0)
})

# Test if the function handles a case where some mice don't have data at day_start
test_that(".getRTV handles cases with missing TV at day_start", {
  missing_data <- test_data %>% dplyr::filter(!(Mouse == 2 & Day == 0))
  result <- .getRTV(missing_data, day_start = 0)
  
  # Mouse 2 should have NA for RTV and logRTV because there is no Day 0 record
  expect_true(all(is.na(result$RTV[result$Mouse == 2])))
  expect_true(all(is.na(result$logRTV[result$Mouse == 2])))
  
  # Other mice should have calculated RTV and logRTV
  expect_false(any(is.na(result$RTV[result$Mouse != 2])))
  expect_false(any(is.na(result$logRTV[result$Mouse != 2])))
})

# Test if the function handles an empty dataset
test_that(".getRTV handles empty dataset", {
  empty_data <- data.frame(Mouse = integer(0), Day = integer(0), TV = numeric(0))
  
  result <- .getRTV(empty_data, day_start = 0)
  
  # The result should be an empty data frame with the expected columns
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
  expect_equal(colnames(result), c("Mouse", "Day", "TV", "RTV", "logRTV", "TV0"))
})

# Test if the function handles a dataset with a single mouse
test_that(".getRTV handles a dataset with a single mouse", {
  single_mouse_data <- test_data[test_data$Mouse == 1, ]
  
  result <- .getRTV(single_mouse_data, day_start = 0)
  
  # RTV should be correctly calculated for the single mouse
  expected_RTV <- c(1, 1.2, 1.4, 1.6, 1.8)
  expect_equal(result$RTV, expected_RTV, tolerance = 1e-5)
  
  # logRTV should be correctly calculated for the single mouse
  expect_equal(result$logRTV, log(expected_RTV), tolerance = 1e-5)
  
  # TV0 should be the same for all rows
  expect_equal(result$TV0, rep(100, 5))
})