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
  combination = "Drug_AB",
  time_start = 0,
  min_observations = 1,
  show_plot = FALSE
)

test_that("Test lmmSynergy with valid input and default parameters (Bliss method)", {
  # Call the function with default method ("Bliss")
  result <- lmmSynergy(model, show_plot = FALSE)
  
  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_named(result, c("Contrasts", "Synergy", "Estimates","nsim"))
  
  # Check that "Contrasts" is a list and "Synergy" is a data frame
  expect_type(result$Contrasts, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj", "Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with HSA method", {
  # Call the function with method = "HSA"
  result <- lmmSynergy(model, method = "HSA", show_plot = FALSE)
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj","Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with RA method", {
  
  # Call the function with method = "RA"
  result <- lmmSynergy(model, method = "RA", show_plot = FALSE)
  
  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_named(result, c("Contrasts", "Synergy", "Estimates", "nsim"))
  
  # Check that "Contrasts" is a NULL and "Synergy" is a data frame
  expect_null(result$Contrasts)
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj","Time") %in% colnames(synergy)))
  
})

test_that("Test lmmSynergy with robust sandwich estimators (robust = TRUE)", {
  # Call the function with robust = TRUE and type = "CR1"
  result <- lmmSynergy(model, robust = TRUE, type = "CR1", show_plot = FALSE)
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Call the function with robust = TRUE and type = "CR2"
  result <- lmmSynergy(model, robust = TRUE, type = "CR2", show_plot = FALSE)
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj" ,"Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy method = 'RA' robust = TRUE works correctly", {
  result <- lmmSynergy(model, method = "RA", min_time = 0, ra_nsim = 1000, robust = TRUE, type = "CR2", show_plot = FALSE)
  
  expect_type(result, "list")
  expect_named(result, c("Contrasts", "Synergy", "Estimates", "nsim"))
  
  # Check that "Contrasts" is a NULL and "Synergy" is a data frame
  expect_null(result$Contrasts)
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj","Time") %in% colnames(synergy)))
  
  # Ensure that robust used clubSandwich to calculate the variance-covariance matrix
  expect_true(all(synergy$Model == "RA"))
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

# Example data and model for testing
set.seed(123)
test_data <- data.frame(
  Mouse = rep(1:10, each = 10),
  Day = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_Z","Drug_ABZ"), each = 10, length.out = 100),
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
  drug_c = "Drug_Z",
  combination = "Drug_ABZ",
  time_start = 0,
  min_observations = 1,
  show_plot = FALSE
)

test_that("Test lmmSynergy with 3 drugs (Bliss method)", {
  # Call the function with default method ("Bliss")
  result <- lmmSynergy(model, show_plot = FALSE)
  
  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_named(result, c("Contrasts", "Synergy", "Estimates", "nsim"))
  
  # Check that "Contrasts" is a list and "Synergy" is a data frame
  expect_type(result$Contrasts, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj","Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with 3 drugs with HSA method", {
  # Call the function with method = "HSA"
  result <- lmmSynergy(model, method = "HSA", show_plot = FALSE)
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj", "Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with RA method", {
  
  # Call the function with method = "RA"
  result <- lmmSynergy(model, method = "RA", show_plot = FALSE)
  
  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_named(result, c("Contrasts", "Synergy", "Estimates","nsim"))
  
  # Check that "Contrasts" is a NULL and "Synergy" is a data frame
  expect_null(result$Contrasts)
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj", "Time") %in% colnames(synergy)))
  
})

test_that("Test lmmSynergy with RA method and 'robust' = TRUE works correctly", {
  
  # Call the function with method = "RA"
  result <- lmmSynergy(model, method = "RA", robust = TRUE, show_plot = FALSE)
  
  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_named(result, c("Contrasts", "Synergy", "Estimates","nsim"))
  
  # Check that "Contrasts" is a NULL and "Synergy" is a data frame
  expect_null(result$Contrasts)
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj", "Time") %in% colnames(synergy)))
  
})

test_that("Test lmmSynergy prints warning message if p-value = 0", {
  # Check that the warning message is printed
  set.seed(245)
  expect_warning(lmmSynergy(model, method = "RA", robust = TRUE, nsim = 10),
                 "p-values below p<1e-01 are approximated to 0. If you used method = 'RA' consider increasing 'nsim' value for more precise p-values.")
})


# Tests for lmmSynergy and Gompertz model----

# Example data and model for testing
set.seed(123)
test_data <- data.frame(
  Mouse = rep(1:10, each = 10),
  Day = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rbeta(10, 3, 1)
)

model <- lmmModel(
  data = test_data,
  grwth_model = "gompertz",
  sample_id = "Mouse",
  time = "Day",
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

test_that("Test lmmSynergy with valid input and default parameters (Bliss method)", {
  # Call the function with default method ("Bliss")
  result <- lmmSynergy(model, show_plot = FALSE)
  
  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 3)
  expect_named(result, c("Synergy", "Estimates","nsim"))
  
  # Check that "Synergy" is a data frame
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj", "Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with HSA method", {
  # Call the function with method = "HSA"
  expect_warning(
    result <- lmmSynergy(model, method = "HSA", show_plot = FALSE))
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj","Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with RA method", {
  
  # Call the function with method = "RA"
  result <- lmmSynergy(model, method = "RA", show_plot = FALSE)
  
  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 3)
  expect_named(result, c("Synergy", "Estimates", "nsim"))
  
  # Check that "Synergy" is a data frame
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj","Time") %in% colnames(synergy)))
  
})

test_that("Test lmmSynergy with robust sandwich estimators (robust = TRUE) returns a warning message", {
  # Call the function with robust = TRUE and type = "CR1"
  expect_warning(
    lmmSynergy(model, robust = TRUE, type = "CR1", show_plot = FALSE),
    "Sandwich-based robust estimators are only available for exponential growth models.")
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

# Example data and model for testing
set.seed(123)
test_data <- data.frame(
  Mouse = rep(1:10, each = 10),
  Day = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_Z","Drug_ABZ"), each = 10, length.out = 100),
  TV = rbeta(10,3,1)
)

model <- lmmModel(
  data = test_data,
  grwth_model = "gompertz",
  sample_id = "Mouse",
  time = "Day",
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

test_that("Test lmmSynergy with 3 drugs (Bliss method)", {
  # Call the function with default method ("Bliss")
  result <- lmmSynergy(model, show_plot = FALSE)
  
  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 3)
  expect_named(result, c("Synergy", "Estimates", "nsim"))
  
  # Check that "Synergy" is a data frame
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj","Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with 3 drugs with HSA method", {
  # Call the function with method = "HSA"
  result <- lmmSynergy(model, method = "HSA", show_plot = FALSE)
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj", "Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with RA method", {
  
  # Call the function with method = "RA"
  result <- lmmSynergy(model, method = "RA", show_plot = FALSE)
  
  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 3)
  expect_named(result, c("Synergy", "Estimates","nsim"))
  
  # Check that "Synergy" is a data frame
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj", "Time") %in% colnames(synergy)))
  
})


test_that("Test lmmSynergy prints warning message if p-value = 0", {
  # Check that the warning message is printed
  set.seed(245)
  expect_warning(lmmSynergy(model, method = "RA", nsim = 10),
                 "p-values below p<1e-01 are approximated to 0. If you used a Gompertz model, consider increasing 'nsim' value for more precise p-values.")
})


