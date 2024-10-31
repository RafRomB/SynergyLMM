# Dummy dataset for testing .getRTV ----
set.seed(123)
test_data <- data.frame(
  SampleID = rep(1:5, each = 5),
  Time = rep(0:4, times = 5),
  TV = c(100, 120, 140, 160, 180,   # Mouse 1
         80, 100, 120, 150, 190,    # Mouse 2
         90, 95, 100, 110, 130,     # Mouse 3
         110, 115, 120, 125, 130,   # Mouse 4
         70, 75, 85, 90, 100)       # Mouse 5
)

# Tests for .getRTV function ----
# Test if .getRTV function correctly calculates RTV and logRTV
test_that(".getRTV correctly calculates RTV and logRTV", {
  result <- .getRTV(test_data, time_start = 0)
  
  # Expected RTV for each Mouse at each Time
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

# Test if .getRTV function correctly adds the TV0 column
test_that(".getRTV correctly adds TV0 column", {
  result <- .getRTV(test_data, time_start = 0)
  
  # Expected TV0 for each Mouse
  expected_TV0 <- rep(c(100, 80, 90, 110, 70), each = 5)
  
  # Check if TV0 is correctly added
  expect_equal(result$TV0, expected_TV0)
})

# Test if .getRTV function handles a case where some mice don't have data at time_start
test_that(".getRTV handles cases with missing TV at time_start", {
  missing_data <- test_data %>% dplyr::filter(!(SampleID == 2 & Time == 0))
  result <- .getRTV(missing_data, time_start = 0)
  
  # Mouse 2 should have NA for RTV and logRTV because there is no Time 0 record
  expect_true(all(is.na(result$RTV[result$SampleID == 2])))
  expect_true(all(is.na(result$logRTV[result$SampleID == 2])))
  
  # Other mice should have calculated RTV and logRTV
  expect_false(any(is.na(result$RTV[result$SampleID != 2])))
  expect_false(any(is.na(result$logRTV[result$SampleID != 2])))
})

# Test if .getRTV function handles an empty dataset
test_that(".getRTV handles empty dataset", {
  empty_data <- data.frame(SampleID = integer(0), Time = integer(0), TV = numeric(0))
  
  result <- .getRTV(empty_data, time_start = 0)
  
  # The result should be an empty data frame with the expected columns
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
  expect_equal(colnames(result), c("SampleID", "Time", "TV", "RTV", "logRTV", "TV0"))
})

# Test if .getRTV function handles a dataset with a single mouse
test_that(".getRTV handles a dataset with a single mouse", {
  single_mouse_data <- test_data[test_data$SampleID == 1, ]
  
  result <- .getRTV(single_mouse_data, time_start = 0)
  
  # RTV should be correctly calculated for the single mouse
  expected_RTV <- c(1, 1.2, 1.4, 1.6, 1.8)
  expect_equal(result$RTV, expected_RTV, tolerance = 1e-5)
  
  # logRTV should be correctly calculated for the single mouse
  expect_equal(result$logRTV, log(expected_RTV), tolerance = 1e-5)
  
  # TV0 should be the same for all rows
  expect_equal(result$TV0, rep(100, 5))
})

# Dummy dataset for testing lmmModel ----
set.seed(123)
test_data <- data.frame(
  SampleID = rep(1:10, each = 10),
  Time = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

# Tests for lmmModel function ----

# Test if lmmModel function runs without errors on valid input
test_that("lmmModel runs without error on valid input", {
  result <- lmmModel(
    data = test_data,
    sample_id = "SampleID",
    time = "Time",
    treatment = "Treatment",
    tumor_vol = "TV",
    trt_control = "Control",
    drug_a = "Drug_A",
    drug_b = "Drug_B",
    drug_ab = "Drug_AB",
    time_start = 0,
    time_end = NULL,
    min_observations = 1,
    show_plot = FALSE
  )
  
  expect_s3_class(result, "lme")
})


# Test if lmmModel function returns an error when required columns are missing
test_that("lmmModel throws an error when required columns are missing", {
  missing_data <- test_data[, -1] # Removing the 'SampleID' column
  
  expect_error(
    lmmModel(
      data = missing_data,
      sample_id = "SampleID",
      time = "Time",
      treatment = "Treatment",
      tumor_vol = "TV",
      trt_control = "Control",
      drug_a = "Drug_A",
      drug_b = "Drug_B",
      drug_ab = "Drug_AB",
      time_start = 0,
      min_observations = 1,
      show_plot = FALSE
    ),
    "The following required columns are missing from the data: SampleID"
  )
})

# Test if lmmModel function returns an error when treatment column contains unrecognized treatments
test_that("lmmModel throws an error when treatment column contains unrecognized treatments", {
  trt_data <- test_data
  trt_data$Treatment[4] <- "unk_trt"
  
  expect_error(
    lmmModel(
      data = trt_data,
      sample_id = "SampleID",
      time = "Time",
      treatment = "Treatment",
      tumor_vol = "TV",
      trt_control = "Control",
      drug_a = "Drug_A",
      drug_b = "Drug_B",
      drug_ab = "Drug_AB",
      time_start = 0,
      min_observations = 1,
      show_plot = FALSE
    ),
    "The treatment column contains unrecognized treatments: unk_trt"
  )
})

# Test if lmmModel function returns an error when an expected treatment is missing
test_that("lmmModel throws an error when an expected treatment is missing", {

  expect_error(
    lmmModel(
      data = test_data,
      sample_id = "SampleID",
      time = "Time",
      treatment = "Treatment",
      tumor_vol = "TV",
      trt_control = "Control",
      drug_a = "Drug_A",
      drug_b = "Drug_B",
      drug_ab = "Drug_X",
      time_start = 0,
      min_observations = 1,
      show_plot = FALSE
    ),
    "The treatment column is missing expected treatments: Drug_X"
  )
})

test_that("lmmModel throws an error when 'min_observations' is a negative value", {
  expect_error(
    lmmModel(
    data = test_data,
    sample_id = "SampleID",
    time = "Time",
    treatment = "Treatment",
    tumor_vol = "TV",
    trt_control = "Control",
    drug_a = "Drug_A",
    drug_b = "Drug_B",
    drug_ab = "Drug_AB",
    time_start = 0,
    min_observations = -1,
    show_plot = FALSE
  ),
  "The `min_observations` parameter must be a positive numeric value."
  )
})

test_that("lmmModel correctly filters data based on time_start", {
  result <- lmmModel(
    data = test_data,
    sample_id = "SampleID",
    time = "Time",
    treatment = "Treatment",
    tumor_vol = "TV",
    trt_control = "Control",
    drug_a = "Drug_A",
    drug_b = "Drug_B",
    drug_ab = "Drug_AB",
    time_start = 5,  # Change time_start to 5
    min_observations = 1,
    show_plot = FALSE
  )
  
  filtered_data <- result$dt1
  expect_true(all(filtered_data$Time >= 0))  # Since time_start was subtracted
})

test_that("lmmModel correctly filters data based on time_end", {
  result <- lmmModel(
    data = test_data,
    sample_id = "SampleID",
    time = "Time",
    treatment = "Treatment",
    tumor_vol = "TV",
    trt_control = "Control",
    drug_a = "Drug_A",
    drug_b = "Drug_B",
    drug_ab = "Drug_AB",
    time_start = 0,
    time_end = 5, # Change time_end to 5
    min_observations = 1,
    show_plot = FALSE
  )
  
  filtered_data <- result$dt1
  expect_true(all(filtered_data$Time <= 5))  # No time beyond 5 should be present
})


test_that("lmmModel excludes samples with TV0 == 0", {
  # Create a small example dataset
  set.seed(123)
  test_data <- data.frame(
    SampleID = rep(1:5, each = 5),
    Time = rep(0:4, times = 5),
    Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), 
                    each = 5, length.out = 25),
    TV = c(0, 120, 140, 160, 180,   # Mouse 1
           80, 100, 120, 150, 190,    # Mouse 2
           90, 95, 100, 110, 130,     # Mouse 3
           0, 115, 120, 125, 130,   # Mouse 4
           70, 75, 85, 90, 100)       # Mouse 5
  )
  
  # Run the lmmModel function with this dataset
  model <- lmmModel(
    data = test_data,
    sample_id = "SampleID",
    time = "Time",
    treatment = "Treatment",
    tumor_vol = "TV",
    trt_control = "Control",
    drug_a = "Drug_A",
    drug_b = "Drug_B",
    drug_ab = "Drug_AB",
    time_start = 0,
    min_observations = 1,
    show_plot = FALSE
  )
  
  # Extract the data used in the model
  model_data <- model$dt1
  
  # Check that Mouse1 and Mouse3 are not in the final data because their TV0 is 0
  expect_false(1 %in% model_data$SampleID)
  expect_false(4 %in% model_data$SampleID)
  
  # Check that only Mouse2 remains
  expect_true(2 %in% model_data$SampleID)
  expect_equal(levels(model_data$SampleID), as.character(c(2,3,5)))
  
  expect_warning(lmmModel(
    data = test_data,
    sample_id = "SampleID",
    time = "Time",
    treatment = "Treatment",
    tumor_vol = "TV",
    trt_control = "Control",
    drug_a = "Drug_A",
    drug_b = "Drug_B",
    drug_ab = "Drug_AB",
    time_start = 0,
    min_observations = 1,
    show_plot = FALSE
  ), "1,4 subjects have measurements with value 0 at time 0. These subjects will be removed from the analysis.")
  
})

test_that("lmmModel handle samples with TV == 0", {
  # Create a small example dataset
  set.seed(123)
  test_data <- data.frame(
    SampleID = rep(1:5, each = 5),
    Time = rep(0:4, times = 5),
    Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), 
                    each = 5, length.out = 25),
    TV = c(65, 120, 140, 160, 180,   # Mouse 1
           80, 100, 120, 0, 190,    # Mouse 2
           90, 95, 100, 110, 130,     # Mouse 3
           50, 115, 120, 125, 130,   # Mouse 4
           70, 75, 85, 90, 0)       # Mouse 5
  )
  
  # Run the lmmModel function with this dataset
  model <- lmmModel(
    data = test_data,
    sample_id = "SampleID",
    time = "Time",
    treatment = "Treatment",
    tumor_vol = "TV",
    trt_control = "Control",
    drug_a = "Drug_A",
    drug_b = "Drug_B",
    drug_ab = "Drug_AB",
    time_start = 0,
    min_observations = 1,
    show_plot = FALSE
  )
  
  # Extract the data used in the model
  model_data <- model$dt1
  
  # Check that there are no values with TV 0 in the final data
  expect_equal(sum(model_data[model_data$SampleID == 2,"TV"] == 0), 0)
  expect_equal(sum(model_data[model_data$SampleID == 5,"TV"] == 0), 0)
})

# Test if lmmModel function respects the min_observations parameter
test_that("lmmModel respects the min_observations parameter", {
  min_obs_data <- test_data
  min_obs_data <- min_obs_data[-10, ]  # Remove last time measurement from Mouse 1
  
  result <- lmmModel(
    data = min_obs_data,
    sample_id = "SampleID",
    time = "Time",
    treatment = "Treatment",
    tumor_vol = "TV",
    trt_control = "Control",
    drug_a = "Drug_A",
    drug_b = "Drug_B",
    drug_ab = "Drug_AB",
    time_start = 0,
    min_observations = 10,  # Require at least 10 observations per mouse
    show_plot = FALSE
  )
  
  # Expect that Mouse 1 is not in the result because it was removed
  expect_false(any(result$data$SampleID == 1))
})


# Test if lmmModel function produces a plot when show_plot is TRUE
test_that("lmmModel produces plot when show_plot is TRUE", {
  expect_output(
    lmmModel(
      data = test_data,
      sample_id = "SampleID",
      time = "Time",
      treatment = "Treatment",
      tumor_vol = "TV",
      trt_control = "Control",
      drug_a = "Drug_A",
      drug_b = "Drug_B",
      drug_ab = "Drug_AB",
      time_start = 0,
      min_observations = 1,
      show_plot = TRUE
    ),
    NA  # Check that no errors occur when plotting (manual visual check might be needed)
  )
})

# Test if lmmModel function passes additional arguments correctly to nlme::lme
test_that("lmmModel passes additional arguments correctly to nlme::lme", {
  result <- lmmModel(
    data = test_data,
    sample_id = "SampleID",
    time = "Time",
    treatment = "Treatment",
    tumor_vol = "TV",
    trt_control = "Control",
    drug_a = "Drug_A",
    drug_b = "Drug_B",
    drug_ab = "Drug_AB",
    time_start = 0,
    min_observations = 1,
    show_plot = FALSE,
    control = nlme::lmeControl(opt = "optim")  # Passing additional argument
  )
  
  expect_equal(result$call$control$opt, "optim")  # Check if the control argument was passed correctly
})

# Test for lmmModel_estimates function ----

set.seed(123)
test_data <- data.frame(
  SampleID = rep(1:10, each = 10),
  Time = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)


test_that("lmmModel_estimates returns a data frame with correct structure", {
  model <- lmmModel(test_data)
  result <- lmmModel_estimates(model)
  
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 6)  # control, drug_a, drug_b, combination, sd_ranef, sd_resid
  expect_equal(rownames(result), "estimate")
  expect_equal(colnames(result), c("control", "drug_a", "drug_b", "combination", "sd_ranef", "sd_resid"))
})

test_that("lmmModel_estimates returns correct values for coefficients and standard deviations", {
  model <- lmmModel(test_data)
  result <- lmmModel_estimates(model)
  
  # Check that the coefficients match the model's fixed effects
  expect_equal(result$control, model$coefficients$fixed[[1]])
  expect_equal(result$drug_a, model$coefficients$fixed[[2]])
  expect_equal(result$drug_b, model$coefficients$fixed[[3]])
  expect_equal(result$combination, model$coefficients$fixed[[4]])
  
  # Check that the standard deviations match the model's random effects and residuals
  expect_equal(result$sd_ranef, sqrt(model$modelStruct$reStruct[[1]][1]))
  expect_equal(result$sd_resid, model$sigma)
})
