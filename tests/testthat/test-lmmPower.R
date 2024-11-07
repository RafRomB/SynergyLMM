# Tests for PostHocPwr ----

# Example data and model for testing
set.seed(123)
test_data <- data.frame(
  Mouse = rep(1:10, each = 10),
  Time = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

model <- lmmModel(
  data = test_data,
  sample_id = "Mouse",
  time = "Time",
  treatment = "Treatment",
  tumor_vol = "TV",
  trt_control = "Control",
  drug_a = "Drug_A",
  drug_b = "Drug_B",
  combination = "Drug_AB",
  time_start = 0,
  min_observations = 1,
  show_plot = FALSE
)

test_that("Test PostHocPwr with valid input and default parameters (Bliss method)", {
  # Call the function with default method ("Bliss")
  result <- PostHocPwr(model, nsim = 10) # Use small nsim for quicker testing
  
  # Check that the result is numeric
  expect_type(result, "double")
  
  # Check that the result is between 0 and 1
  expect_true(result >= 0 && result <= 1)
})

test_that("Test PostHocPwr with HSA method", {
  # Call the function with method = "HSA"
  result <- PostHocPwr(model, nsim = 10, method = "HSA") # Use small nsim for quicker testing
  
  # Check that the result is numeric
  expect_type(result, "double")
  
  # Check that the result is between 0 and 1
  expect_true(result >= 0 && result <= 1)
})

test_that("Test PostHocPwr with HSA method independently of the order of treatment input", {
  
  model <- lmmModel(
    data = test_data,
    sample_id = "Mouse",
    time = "Time",
    treatment = "Treatment",
    tumor_vol = "TV",
    trt_control = "Control",
    drug_a = "Drug_B", # Change order of treatments
    drug_b = "Drug_A", # Change order of treatments
    combination = "Drug_AB",
    time_start = 0,
    min_observations = 1,
    show_plot = FALSE
  )
  
  # Call the function with method = "HSA"
  result <- PostHocPwr(model, nsim = 10, method = "HSA") # Use small nsim for quicker testing
  
  # Check that the result is numeric
  expect_type(result, "double")
  
  # Check that the result is between 0 and 1
  expect_true(result >= 0 && result <= 1)
})


test_that("Test PostHocPwr with invalid method input", {
  # Expect an error when an invalid method is provided
  expect_error(PostHocPwr(model, method = "InvalidMethod"),
               "Invalid 'method' provided. Choose from 'Bliss' or 'HSA'.")
})

test_that("Test PostHocPwr with different values of nsim", {
  # Test with nsim = 1 (edge case)
  result <- PostHocPwr(model, nsim = 1)
  expect_type(result, "double")
  expect_true(result >= 0 && result <= 1)
  
  # Test with nsim = 100 (larger number of simulations)
  result <- PostHocPwr(model, nsim = 100)
  expect_type(result, "double")
  expect_true(result >= 0 && result <= 1)
})

test_that("Test PostHocPwr with different pvalue thresholds", {
  # Test with pvalue = 0.01
  result <- PostHocPwr(model, pvalue = 0.01, nsim = 10)
  expect_type(result, "double")
  expect_true(result >= 0 && result <= 1)
  
  # Test with pvalue = 0.10
  result2 <- PostHocPwr(model, pvalue = 0.50, nsim = 10)
  expect_type(result2, "double")
  expect_true(result2 >= 0 && result <= 1)
  
  # result with lower p-value threshold should be smaller
  expect_true(result < result2)
  
})

test_that("Test PostHocPwr setting a value for time", {
  # Test with time = 5
  result <- PostHocPwr(model, nsim = 100, time = 5)
  expect_type(result, "double")
  expect_true(result >= 0 && result <= 1)
  
  # Test with time = 9
  result <- PostHocPwr(model, nsim = 100, time = 9)
  expect_type(result, "double")
  expect_true(result >= 0 && result <= 1)
  
  
})


test_that("Test PostHocPwr passess arguments to simulateY function", {
  
  # Passing the same seed to initiate the random number generator
  result.seed <- PostHocPwr(model, pvalue = 0.5, nsim = 50, seed = 123)
  result2.seed <- PostHocPwr(model, pvalue = 0.5, nsim = 50, seed = 123)
  expect_true(result.seed == result2.seed)
  
})

# Example data and model for testing
set.seed(123)
test_data <- data.frame(
  Mouse = rep(1:10, each = 10),
  Time = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_Z","Drug_ABZ"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

model <- lmmModel(
  data = test_data,
  sample_id = "Mouse",
  time = "Time",
  treatment = "Treatment",
  tumor_vol = "TV",
  trt_control = "Control",
  drug_a = "Drug_A",
  drug_b = "Drug_B",
  drug_c = "Drug_Z",
  combination = "Drug_ABZ",
  time_start = 0,
  min_observations = 1,
  show_plot = FALSE
)

test_that("Test PostHocPwr with 3 drugs (Bliss method)", {
  # Call the function with default method ("Bliss")
  result <- PostHocPwr(model, nsim = 10) # Use small nsim for quicker testing
  
  # Check that the result is numeric
  expect_type(result, "double")
  
  # Check that the result is between 0 and 1
  expect_true(result >= 0 && result <= 1)
})

test_that("Test PostHocPwr with HSA method", {
  # Call the function with method = "HSA"
  result <- PostHocPwr(model, nsim = 10, method = "HSA") # Use small nsim for quicker testing
  
  # Check that the result is numeric
  expect_type(result, "double")
  
  # Check that the result is between 0 and 1
  expect_true(result >= 0 && result <= 1)
})

# Tests for APrioriPwr ----

test_that("APrioriPwr throws an error when neither sd_eval and sgma_eval nor grwrComb_eval are provided", {
  expect_error(
    APrioriPwr(grwrControl = 0.1, grwrA = 0.1, grwrB = 0.1, grwrComb = 0.1, sd_ranef = 0.5, sgma = 0.5),
    "One of the following, 'sd_eval' and 'sgma_eval', or 'grwrComb_eval', arguments must be specified"
  )
})

test_that("APrioriPwr throws an error when only one of sd_eval or sgma_eval is provided", {
  expect_error(
    APrioriPwr(grwrControl = 0.1, grwrA = 0.1, grwrB = 0.1, grwrComb = 0.1, sd_ranef = 0.5, sgma = 0.5, sd_eval = c(0.2, 0.3)),
    "Both, 'sd_eval' and 'sgma_eval', must be specified"
  )
  expect_error(
    APrioriPwr(grwrControl = 0.1, grwrA = 0.1, grwrB = 0.1, grwrComb = 0.1, sd_ranef = 0.5, sgma = 0.5, sgma_eval = c(0.2, 0.3)),
    "Both, 'sd_eval' and 'sgma_eval', must be specified"
  )
})

test_that("APrioriPwr handles incorrect method input gracefully", {
  expect_error(
    APrioriPwr(
      grwrControl = 0.1,
      grwrA = 0.1,
      grwrB = 0.1,
      grwrComb = 0.1,
      sd_ranef = 0.5,
      sgma = 0.5,
      sd_eval = c(0.2, 0.3),
      sgma_eval = c(0.2, 0.3),
      method = "RA"
    ),
    "Invalid 'method' provided. Choose from 'Bliss' or 'HSA'."
  )
})


test_that("APrioriPwr provides correct output structure when sd_eval and sgma_eval are provided", {
  result <- APrioriPwr(
    grwrControl = 0.1, grwrA = 0.1, grwrB = 0.1, grwrComb = 0.1,
    sd_ranef = 0.5, sgma = 0.5, sd_eval = c(0.2, 0.3), sgma_eval = c(0.2, 0.3)
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("subject", "Treatment", "Time", "mA", "m0") %in% colnames(result)))
  expect_equal(nrow(result), 4 * 5 * 4)  # 4 treatments, 5 subjects per group, 4 times
})

test_that("APrioriPwr provides correct output structure when grwrComb_eval is provided", {
  result <- APrioriPwr(
    grwrControl = 0.1, grwrA = 0.1, grwrB = 0.1, grwrComb = 0.1,
    sd_ranef = 0.5, sgma = 0.5, grwrComb_eval = c(0.1, 0.2, 0.3)
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("subject", "Treatment", "Time", "mA", "m0") %in% colnames(result)))
})

test_that("APrioriPwr plots are generated correctly for different methods", {
  
  result_bliss <- APrioriPwr(
    grwrControl = 0.1,
    grwrA = 0.1,
    grwrB = 0.1,
    grwrComb = 0.1,
    sd_ranef = 0.5,
    sgma = 0.5,
    sd_eval = c(0.2, 0.3),
    sgma_eval = c(0.2, 0.3),
    method = "Bliss"
  )
  
  expect_s3_class(result_bliss, "data.frame")
  expect_true(all(c("subject", "Treatment", "Time", "mA", "m0") %in% colnames(result_bliss)))
  
  result_hsa <- APrioriPwr(
    grwrControl = 0.1,
    grwrA = 0.1,
    grwrB = 0.1,
    grwrComb = 0.1,
    sd_ranef = 0.5,
    sgma = 0.5,
    sd_eval = c(0.2, 0.3),
    sgma_eval = c(0.2, 0.3),
    method = "HSA"
  )
  
  expect_s3_class(result_hsa, "data.frame")
  expect_true(all(c("subject", "Treatment", "Time", "mA", "m0") %in% colnames(result_hsa)))
  
})

test_that("APrioriPwr handles correctly HSA method independently of the order of drug definition in the input", {
  
  expc_output <- "1    56 3.266493 3.266493 0.4272459"
  
  expect_output(APrioriPwr(
    grwrControl = 0.5, grwrA = 0.4, grwrB = 0.3, grwrComb = 0.01,
    sd_ranef = 0.5, sgma = 0.5, grwrComb_eval = c(0.1, 0.2, 0.3),
    method = "HSA"
  ), expc_output)
  
  expect_output(APrioriPwr(
    grwrControl = 0.5, grwrA = 0.3, grwrB = 0.4, grwrComb = 0.01, # Change order of drugs growth rates
    sd_ranef = 0.5, sgma = 0.5, grwrComb_eval = c(0.1, 0.2, 0.3),
    method = "HSA"
  ), expc_output)
})


test_that("APrioriPwr works with edge cases with single values in evaluation vectors", {
  result <- APrioriPwr(
    grwrControl = 0.1, grwrA = 0.1, grwrB = 0.1, grwrComb = 0.1,
    sd_ranef = 0.5, sgma = 0.5, sd_eval = c(0.3), sgma_eval = c(0.3)
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("subject", "Treatment", "Time", "mA", "m0") %in% colnames(result)))
})

test_that("APrioriPwr behaves correctly when all evaluation parameters are provided", {
  result <- APrioriPwr(
    grwrControl = 0.1, grwrA = 0.1, grwrB = 0.1, grwrComb = 0.1,
    sd_ranef = 0.5, sgma = 0.5, sd_eval = c(0.2, 0.3), sgma_eval = c(0.2, 0.3), grwrComb_eval = c(0.1, 0.2)
  )
  expect_s3_class(result, "data.frame")
  expect_true(all(c("subject", "Treatment", "Time", "mA", "m0") %in% colnames(result)))
})

# Tests for PwrSampleSize ----

test_that("PwrSampleSize handles different sample sizes correctly and returns a data frame", {
  result <- PwrSampleSize(
    npg = c(5, 8, 10),
    time = c(0, 3, 5, 10),
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.03,
    sd_ranef = 0.01,
    sgma = 0.1,
    method = "Bliss"
  )
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)  # N and Power
  expect_equal(colnames(result), c("N", "Power"))
  expect_equal(nrow(result), length(c(5, 8, 10)))
})

test_that("PwrSampleSize handles a single value for npg correctly", {
  result <- PwrSampleSize(
    npg = c(5),
    time = c(0, 3, 5, 10),
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.03,
    sd_ranef = 0.01,
    sgma = 0.1,
    method = "Bliss"
  )
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 1)
  expect_equal(result$N, 5)
})

test_that("PwrSampleSize handles the 'HSA' method correctly", {
  result <- PwrSampleSize(
    npg = c(5, 8),
    time = c(0, 3, 5, 10),
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.03,
    sd_ranef = 0.01,
    sgma = 0.1,
    method = "HSA"
  )
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), length(c(5, 8)))
})

test_that("PwrSampleSize generates the correct plots without errors for 'Bliss' method", {
  expect_silent(
    PwrSampleSize(
      npg = c(5, 8, 10),
      time = c(0, 3, 5, 10),
      grwrControl = 0.08,
      grwrA = 0.07,
      grwrB = 0.06,
      grwrComb = 0.03,
      sd_ranef = 0.01,
      sgma = 0.1,
      method = "Bliss"
    )
  )
})

test_that("PwrSampleSize generates the correct plots without errors for 'HSA' method", {
  expect_silent(
    PwrSampleSize(
      npg = c(5, 8, 10),
      time = c(0, 3, 5, 10),
      grwrControl = 0.08,
      grwrA = 0.07,
      grwrB = 0.06,
      grwrComb = 0.03,
      sd_ranef = 0.01,
      sgma = 0.1,
      method = "HSA"
    )
  )
})

test_that("PwrSampleSize handles different numbers of times correctly", {
  result <- PwrSampleSize(
    npg = c(5, 8, 10),
    time = c(0, 5, 10),
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.03,
    sd_ranef = 0.01,
    sgma = 0.1,
    method = "Bliss"
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), length(c(5, 8, 10)))
})

test_that("PwrSampleSize handles incorrect method input gracefully", {
  expect_error(
    PwrSampleSize(
      npg = c(5, 8, 10),
      time = c(0, 3, 5, 10),
      grwrControl = 0.08,
      grwrA = 0.07,
      grwrB = 0.06,
      grwrComb = 0.03,
      sd_ranef = 0.01,
      sgma = 0.1,
      method = "RA"
    ),
    "Invalid 'method' provided. Choose from 'Bliss' or 'HSA'."
  )
})

# Tests for PwrTime ----

test_that("PwrTime returns a data frame with correct columns for 'max' type", {
  result <- PwrTime(
    npg = 5,
    time = list(seq(0, 9, 3), seq(0, 21, 3), seq(0, 30, 3)),
    type = "max",
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.03,
    sd_ranef = 0.01,
    sgma = 0.1,
    method = "Bliss"
  )
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)  # Time and Power
  expect_equal(colnames(result), c("Time", "Power"))
  expect_equal(nrow(result), 3)  # 3 different max follow-up times
})

test_that("PwrTime returns a data frame with correct columns for 'freq' type", {
  result <- PwrTime(
    npg = 5,
    time = list(seq(0, 30, 10), seq(0, 30, 5), seq(0, 30, 3)),
    type = "freq",
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.03,
    sd_ranef = 0.01,
    sgma = 0.1,
    method = "Bliss"
  )
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)  # Time and Power
  expect_equal(colnames(result), c("Time", "Power"))
  expect_equal(nrow(result), 3)  # 3 different frequencies
})

test_that("Warning is thrown when 'type' is 'max' and times have the same maximum", {
  expect_warning(
    PwrTime(
      npg = 5,
      time = list(seq(0, 21, 3), seq(0, 21, 5), seq(0, 21, 7)),
      type = "max",
      grwrControl = 0.08,
      grwrA = 0.07,
      grwrB = 0.06,
      grwrComb = 0.03,
      sd_ranef = 0.01,
      sgma = 0.1,
      method = "Bliss"
    ),
    "Your list 'time' has several vectors with the same maximum time of follow-up."
  )
})

test_that("PwrTime handles the 'HSA' method correctly", {
  result <- PwrTime(
    npg = 5,
    time = list(seq(0, 9, 3), seq(0, 21, 3)),
    type = "max",
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.03,
    sd_ranef = 0.01,
    sgma = 0.1,
    method = "HSA"
  )
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 2)
})

test_that("PwrTime handles incorrect method input gracefully", {
  expect_error(
    PwrTime(
      npg = 5,
      time = list(seq(0, 30, 10), seq(0, 30, 5), seq(0, 30, 3)),
      type = "freq",
      grwrControl = 0.08,
      grwrA = 0.07,
      grwrB = 0.06,
      grwrComb = 0.03,
      sd_ranef = 0.01,
      sgma = 0.1,
      method = "RA"
    ),
    "Invalid 'method' provided. Choose from 'Bliss' or 'HSA'."
  )
})

test_that("PwrTime generates the correct plots without errors for 'Bliss' method", {
  expect_silent(PwrTime(
    npg = 5,
    time = list(seq(0, 9, 3), seq(0, 21, 3), seq(0, 30, 3)),
    type = "max",
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.03,
    sd_ranef = 0.01,
    sgma = 0.1,
    method = "Bliss"
  ))
})

test_that("PwrTime generates the correct plots without errors for 'HSA' method", {
  expect_silent(PwrTime(
    npg = 5,
    time = list(seq(0, 9, 3), seq(0, 21, 3), seq(0, 30, 3)),
    type = "freq",
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.03,
    sd_ranef = 0.01,
    sgma = 0.1,
    method = "HSA"
  ))
})

test_that("PwrTime handles incorrect 'type' input gracefully", {
  expect_error(
    PwrTime(
      npg = 5,
      time = list(seq(0, 9, 3), seq(0, 21, 3), seq(0, 30, 3)),
      type = "invalidType",
      grwrControl = 0.08,
      grwrA = 0.07,
      grwrB = 0.06,
      grwrComb = 0.03,
      sd_ranef = 0.01,
      sgma = 0.1,
      method = "Bliss"
    ),
    "invalidType: Invalid 'type' provided. Choose from 'max' or 'freq'."
  )
})

