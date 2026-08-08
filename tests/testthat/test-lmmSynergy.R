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

test_that("Test lmmSynergy with valid input and default parameters (exponential, Bliss method)", {
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
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with p-value adjustment (exponential)", {
  # Call the function with default method ("Bliss")
  result <- lmmSynergy(model, show_plot = FALSE, padj = "BH")
  
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


test_that("Test lmmSynergy with HSA method (exponential)", {
  # Call the function with method = "HSA"
  result <- lmmSynergy(model, method = "HSA", show_plot = FALSE)
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval","Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with RA method (exponential)", {
  
  # Call the function with method = "RA"
  result <- lmmSynergy(model, method = "RA", nsim = 100, show_plot = FALSE)

  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_named(result, c("Contrasts", "Synergy", "Estimates", "nsim"))

  # Check that "Contrasts" is a NULL and "Synergy" is a data frame
  expect_null(result$Contrasts)
  expect_s3_class(result$Synergy, "data.frame")

  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval","Time") %in% colnames(synergy)))

})

test_that("Test lmmSynergy with robust sandwich estimators (exponential, robust = TRUE)", {
  # Call the function with robust = TRUE and type = "CR1"
  result <- lmmSynergy(model, robust = TRUE, type = "CR1", show_plot = FALSE)
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Call the function with robust = TRUE and type = "CR2"
  result <- lmmSynergy(model, robust = TRUE, type = "CR2", show_plot = FALSE, padj = "BH")
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj" ,"Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy method = 'RA' robust = TRUE works correctly (exponential)", {
  result <- lmmSynergy(model, method = "RA", min_time = 0, nsim = 100, robust = TRUE, type = "CR2", padj = "BH", show_plot = FALSE)
  
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

test_that("Test lmmSynergy with different values of min_time (exponential)", {
  # Call the function with min_time = 5
  result <- lmmSynergy(model, min_time = 5, show_plot = FALSE)

  # Check that only times >= 5 are included
  expect_true(all(result$Synergy$Time >= 5))
})

test_that("Test lmmSynergy plotting functionality with show_plot = TRUE (exponential)", {
  # Check that no error is thrown and a plot is generated. Messages and warnings
  # are suppressed on purpose: whether the '+Inf' capping message or the
  # 'p-values approximated to 0' warning is emitted depends on the simulated
  # draws, which are not reproducible across BLAS implementations.
  expect_no_error(suppressWarnings(suppressMessages(
    lmmSynergy(model, show_plot = TRUE)
  )))
})

test_that("Test lmmSynergy with incorrect method input (exponential)", {
  # Expect an error when an invalid method is provided
  expect_error(lmmSynergy(model, method = "InvalidMethod"),
               "Invalid 'method' provided. Choose from 'Bliss', 'HSA', or 'RA'.",
               fixed = TRUE)
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

test_that("Test lmmSynergy with 3 drugs (exponential, Bliss method)", {
  # Call the function with default method ("Bliss")
  result <- lmmSynergy(model, padj = "BH", show_plot = FALSE)
  
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

test_that("Test lmmSynergy with 3 drugs with HSA method (exponential)", {
  # Call the function with method = "HSA"
  result <- lmmSynergy(model, method = "HSA", padj = "BH", show_plot = FALSE)
  
  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "padj", "Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with RA method (exponential, 3 drugs)", {

  # Call the function with method = "RA"
  result <- lmmSynergy(model, method = "RA", nsim = 100, show_plot = FALSE)

  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_named(result, c("Contrasts", "Synergy", "Estimates","nsim"))

  # Check that "Contrasts" is a NULL and "Synergy" is a data frame
  expect_null(result$Contrasts)
  expect_s3_class(result$Synergy, "data.frame")

  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Time") %in% colnames(synergy)))

})

test_that("Test lmmSynergy with RA method and 'robust' = TRUE works correctly (exponential, 3 drugs)", {

  # Call the function with method = "RA"
  result <- lmmSynergy(model, method = "RA", nsim = 100, robust = TRUE, show_plot = FALSE)
  
  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_named(result, c("Contrasts", "Synergy", "Estimates","nsim"))
  
  # Check that "Contrasts" is a NULL and "Synergy" is a data frame
  expect_null(result$Contrasts)
  expect_s3_class(result$Synergy, "data.frame")
  
  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Time") %in% colnames(synergy)))
  
})

test_that("Test lmmSynergy warns about p-values = 0 (exponential, 3 drugs, robust)", {
  # Whether any simulated p-value is exactly 0 depends on the Monte Carlo draws,
  # which are not reproducible across BLAS/LAPACK implementations. Assert the
  # link between a zero p-value and the warning instead of the warning itself.
  # The exact wording of the message is tested in test-utils.R.
  w <- capture_warnings(
    result <- lmmSynergy(model, method = "RA", robust = TRUE, nsim = 10, show_plot = FALSE)
  )
  expect_equal(any(grepl("are approximated to 0", w, fixed = TRUE)),
               any(result$Synergy$pval == 0))
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

test_that("Test lmmSynergy with valid input and default parameters (Gompertz, Bliss method)", {
  # Call the function with default method ("Bliss")
  result <- lmmSynergy(model, nsim = 10, padj = "BH", show_plot = FALSE)

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

test_that("Test lmmSynergy with HSA method (Gompertz)", {
  # Call the function with method = "HSA". All Gompertz p-values are simulated,
  # so a 'p-values approximated to 0' warning may or may not be emitted
  # depending on the draws; it must not be asserted (see test-utils.R).
  w <- capture_warnings(
    result <- lmmSynergy(model, method = "HSA", nsim = 10, show_plot = FALSE))
  expect_equal(any(grepl("are approximated to 0", w, fixed = TRUE)),
               any(result$Synergy$pval == 0))

  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")

  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval","Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with RA method (Gompertz)", {

  # Call the function with method = "RA"
  result <- lmmSynergy(model, method = "RA", nsim = 10, show_plot = FALSE)

  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 3)
  expect_named(result, c("Synergy", "Estimates", "nsim"))

  # Check that "Synergy" is a data frame
  expect_s3_class(result$Synergy, "data.frame")

  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval","Time") %in% colnames(synergy)))

})

test_that("Test lmmSynergy with robust sandwich estimators (Gompertz) returns a warning message", {
  # Call the function with robust = TRUE and type = "CR1"
  expect_warning(
    lmmSynergy(model, robust = TRUE, nsim = 10, type = "CR1", show_plot = FALSE),
    "Sandwich-based robust estimators are only available for exponential growth models",
    fixed = TRUE
    )
  })


test_that("Test lmmSynergy with different values of min_time (Gompertz)", {
  # Call the function with min_time = 5
  result <- suppressWarnings(lmmSynergy(model, min_time = 5, nsim = 10, show_plot = FALSE))

  # Check that only times >= 5 are included
  expect_true(all(result$Synergy$Time >= 5))
})

test_that("Test lmmSynergy plotting functionality with show_plot = TRUE (Gompertz)", {
  # Check that no error is thrown and a plot is generated. See the note on the
  # equivalent exponential test: emitted messages/warnings depend on the draws.
  expect_no_error(suppressWarnings(suppressMessages(
    lmmSynergy(model, nsim = 10, show_plot = TRUE)
  )))
})

test_that("Test lmmSynergy with incorrect method input (Gompertz)", {
  # Expect an error when an invalid method is provided
  expect_error(lmmSynergy(model, method = "InvalidMethod"),
               "Invalid 'method' provided. Choose from 'Bliss', 'HSA', or 'RA'.",
               fixed = TRUE)
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

test_that("Test lmmSynergy with 3 drugs (Gompertz, Bliss method)", {
  # Call the function with default method ("Bliss")
  result <- lmmSynergy(model, nsim = 10, padj = "BH", show_plot = FALSE)

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

test_that("Test lmmSynergy with 3 drugs with HSA method (Gompertz)", {
  # Call the function with method = "HSA"
  result <- lmmSynergy(model, method = "HSA", nsim = 10, show_plot = FALSE)

  # Check that the result is structured as expected
  expect_type(result, "list")
  expect_s3_class(result$Synergy, "data.frame")

  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Time") %in% colnames(synergy)))
})

test_that("Test lmmSynergy with RA method (Gompertz, 3 drugs)", {

  # Call the function with method = "RA"
  result <- lmmSynergy(model, method = "RA", nsim = 10, show_plot = FALSE)

  # Check that the result is a list with two elements
  expect_type(result, "list")
  expect_equal(length(result), 3)
  expect_named(result, c("Synergy", "Estimates","nsim"))

  # Check that "Synergy" is a data frame
  expect_s3_class(result$Synergy, "data.frame")

  # Check the structure of 'Synergy' dataframe
  synergy <- result$Synergy
  expect_true(all(c("Model", "Metric", "Estimate", "lwr", "upr", "pval", "Time") %in% colnames(synergy)))

})


test_that("Test lmmSynergy warns about p-values = 0 (Gompertz, 3 drugs)", {
  # Whether any simulated p-value is exactly 0 depends on the Monte Carlo draws,
  # which are not reproducible across BLAS/LAPACK implementations. Assert the
  # link between a zero p-value and the warning instead of the warning itself.
  # The exact wording of the message is tested in test-utils.R.
  w <- capture_warnings(
    result <- lmmSynergy(model, method = "RA", nsim = 10, show_plot = FALSE)
  )
  expect_equal(any(grepl("are approximated to 0", w, fixed = TRUE)),
               any(result$Synergy$pval == 0))
})


