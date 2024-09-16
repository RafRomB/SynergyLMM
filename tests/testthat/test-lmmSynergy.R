# Tests for lmmSynergy ----

# Example data and model for testing
set.seed(123)
test_data <- data.frame(
  Mouse = rep(1:10, each = 10),
  Day = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

model <- lmmModel(
  data = test_data,
  sample_id = "Mouse",
  time = "Day",
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

test_that("Test lmmSynergy with valid input and default parameters (Bliss method)", {
  # Call the function with default method ("Bliss")
  result <- lmmSynergy(model, show_plot = FALSE)
  
  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 2)
  expect_named(result, c("Contrasts", "Synergy"))
  
  # Check that "Contrasts" is a list and "Synergy" is a data frame
  expect_type(result$Contrasts, "list")
  expect_s3_class(result$Synergy, "data.frame")
})

test_that("Test lmmSynergy with HSA method", {
  # Call the function with method = "HSA"
  result <- lmmSynergy(model, method = "HSA", show_plot = FALSE)
  
  # Check that the contrast used is one of the HSA contrasts
  expect_true(any(grepl("b4 = b2", result$Contrasts[[1]]$term)) || 
                any(grepl("b4 = b3", result$Contrasts[[1]]$term)))
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
})

test_that("Test lmmSynergy with RA method and normality tests", {
  
  # Check that the Shapiro - Wilk normality test output is printed
  expect_output(lmmSynergy(model, method = "RA", norm_test = "shapiroTest", show_plot = FALSE), "Shapiro - Wilk Normality Test")
  
  # Check that the D'Agostino normality test output is printed
  expect_output(lmmSynergy(model, method = "RA", norm_test = "dagoTest", show_plot = FALSE), "D'Agostino Normality Test")
  
  # Check that the Anderson-Darling normality test output is printed
  expect_output(lmmSynergy(model, method = "RA", norm_test = "adTest", show_plot = FALSE), "Anderson - Darling Normality Test")
  
  # Call the function with method = "RA"
  result <- lmmSynergy(model, method = "RA", show_plot = FALSE)
  
  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 2)
  expect_named(result, c("Contrasts", "Synergy"))
  
  # Check that "Contrasts" is a list and "Synergy" is a data frame
  expect_type(result$Contrasts, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
})

test_that("Test lmmSynergy with robust standard errors (robustSE = TRUE)", {
  # Call the function with robustSE = TRUE and type = "CR1"
  result <- lmmSynergy(model, robustSE = TRUE, type = "CR1", show_plot = FALSE)
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Call the function with robustSE = TRUE and type = "CR2"
  result <- lmmSynergy(model, robustSE = TRUE, type = "CR2", show_plot = FALSE)
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
})

test_that("Test lmmSynergy with different values of min_time", {
  # Call the function with min_time = 5
  result <- lmmSynergy(model, min_time = 5)
  
  # Check that only times >= 5 are included
  expect_true(all(result$Synergy$Day >= 5))
})

test_that("Test lmmSynergy plotting functionality with show_plot = TRUE", {
  # Check that no error is thrown and a plot is generated
  expect_silent(lmmSynergy(model, show_plot = TRUE))
})

test_that("Test lmmSynergy with incorrect method input", {
  # Expect an error when an invalid method is provided
  expect_error(lmmSynergy(model, method = "InvalidMethod"),
               "Invalid 'method' provided. Choose from 'Bliss', 'HSA', or 'RA'.")
})

