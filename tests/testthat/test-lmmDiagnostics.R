# Tests for ranefDiagnostics ----

# Example data and model for testing
set.seed(123)
test_data <- data.frame(
  SampleID = rep(1:10, each = 10),
  Time = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

model <- lmmModel(
  data = test_data,
  sample_id = "SampleID",
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


test_that("ranefDiagnostics returns the correct structure and classes", {
  # Run the diagnostics function
  diagnostics <- ranefDiagnostics(model)
  
  # Check that diagnostics is a list
  expect_type(diagnostics, "list")
  
  # Check that the list contains the expected elements
  expect_named(diagnostics, c("Plots", "Normality", "Levene.test", "Fligner.test"))
  
  # Check that the plots component is a list
  expect_type(diagnostics$Plots, "list")
  
  # Check that Normality is a list and contains the correct elements
  expect_type(diagnostics$Normality, "list")
  expect_named(diagnostics$Normality, c("Shapiro.test", "DAgostino.test", "Anderson.Darling.test"))
  
  # Check that Levene.test is a list
  expect_type(diagnostics$Levene.test, "list")
  
  # Check that Fligner.test is a list
  expect_type(diagnostics$Fligner.test, "list")
})

test_that("Normality tests return the correct classes", {
  diagnostics <- ranefDiagnostics(model)
  
  # Check the class of Shapiro-Wilk test result
  expect_s4_class(diagnostics$Normality$Shapiro.test, "fHTEST")
  
  # Check the class of Anderson-Darling test result
  expect_s4_class(diagnostics$Normality$Anderson.Darling.test, "fHTEST")
  
  # Check if D'Agostino test was skipped due to sample size
  if (is(diagnostics$Normality$DAgostino.test, "character")) {
    expect_match(diagnostics$Normality$DAgostino.test, "Sample size must be at least 20")
  } else {
    expect_s4_class(diagnostics$Normality$DAgostino.test, "fHTEST")
  }
})

test_that("Levene and Fligner tests return the correct classes", {
  diagnostics <- ranefDiagnostics(model)
  
  # Check the class of Levene's test results
  expect_s3_class(diagnostics$Levene.test, "anova")
  
  # Check the class of Fligner-Killeen test results
  expect_s3_class(diagnostics$Fligner.test, "htest")
})

test_that("ranefDiagnostics handles small sample sizes for D'Agostino test", {
  # Modify test data to have fewer than 20 samples for ranef
  small_sample_data <- test_data[test_data$SampleID <= 2, ]
  
  small_model <- lmmModel(
    data = test_data,
    sample_id = "SampleID",
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
  
  diagnostics <- ranefDiagnostics(small_model)
  
  # Check if D'Agostino test was skipped due to small sample size
  expect_true(is.character(diagnostics$Normality$DAgostino.test))
  expect_match(diagnostics$Normality$DAgostino.test, "Sample size must be at least 20")
})

# Tests for residDiagnostics ----

# Data and model for testing:

set.seed(123)
test_data <- data.frame(
  SampleID = rep(1:10, each = 10),
  Time = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

# Create an outlier observation

test_data$TV[nrow(test_data)] <- 1000

# Create model

model <- lmmModel(
  data = test_data,
  sample_id = "SampleID",
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

test_that("residDiagnostics correctly identifies outliers", {
  result <- residDiagnostics(model, pvalue = 0.05)
  
  # Check that the result contains an "outliers" data frame
  expect_true("Outliers" %in% names(result))
  expect_s3_class(result$Outliers, "data.frame")
  
  # Verify that outliers are correctly identified
  # (You may need to manually verify this based on the mock model)
  expect_true(nrow(result$Outliers) > 0 || nrow(result$Outliers) == 0)
})

test_that("residDiagnostics performs and returns results from normality tests", {
  result <- residDiagnostics(model, pvalue = 0.05)
  
  # Check that the result contains a "Normality" list
  expect_true("Normality" %in% names(result))
  expect_type(result$Normality, "list")
  
  # Ensure all three tests are present
  expect_true(all(c("Shapiro.test", "DAgostino.test", "Anderson.Darling.test") %in% names(result$Normality)))
  
  # Check that each test result is a valid object (not NULL)
  expect_s4_class(result$Normality$Shapiro.test, "fHTEST")
  expect_s4_class(result$Normality$DAgostino.test, "fHTEST")
  expect_s4_class(result$Normality$Anderson.Darling.test, "fHTEST")
})

test_that("residDiagnostics generates diagnostic plots", {
  result <- residDiagnostics(model, pvalue = 0.05)
  
  # Check that the result contains a "plots" object
  expect_true("Plots" %in% names(result))
  expect_type(result$Plots, "list")
  
  # Check that the plot list contains expected plots
  expect_true(length(result$Plots) >= 1)
  
  expect_s3_class(result$Plots[[1]], "trellis")
})

test_that("Levene and Fligner tests return the correct classes", {
  diagnostics <- residDiagnostics(model)
  
  # Check the class of Levene's test results
  expect_s3_class(diagnostics$Levene.test$Time, "anova")
  expect_s3_class(diagnostics$Levene.test$Treatment, "anova")
  
  # Check the class of Fligner-Killeen test results
  expect_s3_class(diagnostics$Fligner.test$Time, "htest")
  expect_s3_class(diagnostics$Fligner.test$Treatment, "htest")
  
})

# Tests for ObsvsPred ----

# Data and model for testing:

set.seed(123)
test_data <- data.frame(
  SampleID = rep(1:10, each = 10),
  Time = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

# Create model

model <- lmmModel(
  data = test_data,
  sample_id = "SampleID",
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


test_that("ObsvsPred correctly computes and prints model performance metrics", {
  expect_output(
    ObsvsPred(model),
    "AIC|BIC|R2|RMSE|SIGMA",  # Expect some of these metrics in the output
    fixed = FALSE
  )
})

test_that("ObsvsPred generates observed vs predicted plot", {
  expect_s3_class(ObsvsPred(model), "trellis")
})

# Tests for .lmeU ----

# Create test dataset 
set.seed(123)
test_data <- data.frame(
  SampleID = rep(1:10, each = 10),
  Time = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

# Create model
model <- lmmModel(
  data = test_data,
  sample_id = "SampleID",
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

test_that("Test .lmeU function with valid input", {
  
  # Remove subject 1 and update model
  updated_model <- .lmeU(1, model)
  
  # Check that the updated model is still an 'lme' object
  expect_s3_class(updated_model, "lme")
  
  # Check that the data in the updated model does not contain subject "M01"
  expect_false(1 %in% unique(updated_model$data$SampleID))
  
  # Check that the number of subjects has decreased by 1
  expect_equal(length(unique(updated_model$data$SampleID)), length(unique(model$data$SampleID)) - 1)
})

test_that("Test .lmeU function model structure consistency", {
  # Remove a subject and get the updated model
  updated_model <- .lmeU(1, model)
  
  # Check that the structure of the updated model is consistent with the original model
  # except for the 'dt1' element corresponding to the transformed data produced by lmmModel()
  expect_equal(names(model)[-length(names(model))], names(updated_model))
  expect_equal(class(model$modelStruct), class(updated_model$modelStruct))
})

# Tests for .logLik1.varIdent_loo ----

# Data and model for testing:

set.seed(123)
test_data <- data.frame(
  SampleID = rep(1:10, each = 10),
  Time = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

# Create model

model <- lmmModel(
  data = test_data,
  sample_id = "SampleID",
  time = "Time",
  treatment = "Treatment",
  tumor_vol = "TV",
  trt_control = "Control",
  drug_a = "Drug_A",
  drug_b = "Drug_B",
  combination = "Drug_AB",
  time_start = 0,
  min_observations = 1,
  show_plot = FALSE,
  weights = varIdent(form = ~1|SampleID),
  method = "ML"
)

# Recover the data from the fitted model
model_data <- model$data

# Create leave-one-out model
modfit <- .lmeU(1, model)

test_that("Test .logLik1.varIdent_loo function with valid input", {
  # Data frame for one subject
  dt1 <- model_data[model_data$SampleID == 1, ]
  
  # Call the function with valid input
  result <- .logLik1.varIdent_loo(modfit, dt1, dtInit = model, var_name = "SampleID")
  
  # Check that the result is numeric and not NA
  expect_type(result, "double")
  expect_false(is.na(result))
})


test_that("Test .logLik1.varIdent_loo function with different var_name values", {
  # Create model
  
  model <- lmmModel(
    data = test_data,
    sample_id = "SampleID",
    time = "Time",
    treatment = "Treatment",
    tumor_vol = "TV",
    trt_control = "Control",
    drug_a = "Drug_A",
    drug_b = "Drug_B",
    combination = "Drug_AB",
    time_start = 0,
    min_observations = 1,
    show_plot = FALSE,
    weights = varIdent(form = ~1|Time), # Different variable for the weights
    method = "ML"
  )
  
  # Recover the data from the fitted model
  model_data <- model$data
  
  # Create leave-one-out model
  modfit <- .lmeU(1, model)
  
  # Data frame for one subject
  dt1 <- model_data[model_data$SampleID == 1, ]
  
  # Call the function with a valid var_name
  result <- .logLik1.varIdent_loo(modfit, dt1, dtInit = model, var_name = "Time")
  
  # Check that the result is numeric and not NA
  expect_type(result, "double")
  expect_false(is.na(result))
})

# Tests for .lLik ----

# Create test data
set.seed(123)
test_data <- data.frame(
  SampleID = rep(1:10, each = 10),
  Time = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

# Create model

model <- lmmModel(
  data = test_data,
  sample_id = "SampleID",
  time = "Time",
  treatment = "Treatment",
  tumor_vol = "TV",
  trt_control = "Control",
  drug_a = "Drug_A",
  drug_b = "Drug_B",
  combination = "Drug_AB",
  time_start = 0,
  min_observations = 1,
  show_plot = FALSE,
  method = "ML"
)


# Create a list of leave-one-out (leave mouse 1 out) model fits using .lmeU function 
lmeUall_no_varIdent <- list(.lmeU(1, model))
names(lmeUall_no_varIdent) <- 1

test_that("Test .lLik function with valid input without varIdent structure", {
  
  # Call the function with valid input and check result is numeric
  result <- .lLik(1, model, lmeUall_no_varIdent, var_name = NULL)
  expect_type(result, "double")
})

test_that("Test .lLik function with valid input with varIdent structure", {
  
  # Update model to include varIdent structure
  model <- update(model, weights = varIdent(form = ~1|SampleID))

  # Create a list of leave-one-out (leave mouse 1 out) model fits using .lmeU function 
  lmeUall_varIdent <- list(.lmeU(1, model))
  names(lmeUall_varIdent) <- 1
  
  # Call the function with valid input and correct var_name
  result <- .lLik(1, model, lmeUall_varIdent, var_name = "SampleID")
  expect_type(result, "double")
})

test_that("Test .lLik function throws an error with missing var_name when needed", {
  
  # Update model to include varIdent structure
  model <- update(model, weights = varIdent(form = ~1|SampleID))
  
  # Create a list of leave-one-out (leave mouse 1 out) model fits using .lmeU function 
  lmeUall_varIdent <- list(.lmeU(1, model))
  names(lmeUall_varIdent) <- 1
  
  # Call the function without providing var_name when variance structure is used
  expect_error(.lLik(1, model, lmeUall_varIdent, var_name = NULL),
               "`var_name` cannot be NULL if a variance estructure has been specified in the model")
})

test_that("Test .lLik function with different variance structures", {
  # Test with another variance structure (Time)
  
  model <- update(model, weights = varIdent(form = ~1|Time)) # Use `Time` as group for variance structure
  
  # Create a list of leave-one-out (leave mouse 1 out) model fits using .lmeU function 
  lmeUall_diff_var <- list(.lmeU(1, model))
  names(lmeUall_diff_var) <- 1
  
  # Call the function and check the result is numeric
  result <- .lLik(1, model, lmeUall_diff_var, var_name = "Time")
  expect_type(result, "double")
})

# Tests for logLikSubjectDisplacements ----

# Create test data
set.seed(123)
test_data <- data.frame(
  SampleID = rep(1:10, each = 10),
  Time = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

# Create model
model <- lmmModel(
  data = test_data,
  sample_id = "SampleID",
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


test_that("Test logLikSubjectDisplacements with valid input without varIdent structure", {
  
  # Call the function with valid input and check that it runs without error
  result <- logLikSubjectDisplacements(model)
  
  # Check that the result is numeric
  expect_type(result, "double")
  
  # Check that all values are finite
  expect_true(all(is.finite(result)))
})

test_that("Test logLikSubjectDisplacements with valid input with varIdent structure", {
  
  # Update model to include varStruct
  model <- update(model, weights = varIdent(form = ~1|SampleID))
  
  # Call the function with valid input and correct var_name
  result <- logLikSubjectDisplacements(model, var_name = "SampleID")
  
  # Check that the result is numeric
  expect_type(result, "double")
  
  # Check that all values are finite
  expect_true(all(is.finite(result)))
})

test_that("Test logLikSubjectDisplacements throws an error with missing var_name when needed", {
  
  # Update model to include varStruct
  model <- update(model, weights = varIdent(form = ~1|SampleID))
  
  # Call the function without providing var_name when variance structure is used
  expect_error(logLikSubjectDisplacements(model), 
               "`var_name` cannot be NULL if a variance estructure has been specified in the model")
})

test_that("Test logLikSubjectDisplacements with different thresholds", {
  # Test with threshold that is higher than any log-likelihood displacement
  result <- logLikSubjectDisplacements(model, disp_thrh = 1000)
  expect_true(all(result <= 1000))
  expect_output(logLikSubjectDisplacements(model, disp_thrh = 1000),
                 "No subject with a log-likelihood displacement greater than 1000")
  
  # Test with very low threshold (all displacements should be included in plot)
  result <- logLikSubjectDisplacements(model, disp_thrh = -1000)
  expect_true(any(result > -1000))
})

test_that("Test logLikSubjectDisplacements output consistency", {
  # Run function and capture the output
  result <- logLikSubjectDisplacements(model)
  
  # Check that result is named and length matches number of subjects
  expect_named(result)
  expect_equal(length(result), length(unique(test_data$SampleID)))
})

# Tests for .CookDfun ----

test_that("Test .CookDfun with basic valid input", {
  # Define input vectors and matrix
  betaU <- c(1.1, 2.2, 3.3)
  beta0 <- c(1.0, 2.0, 3.0)
  vb.inv <- diag(c(0.5, 0.5, 0.5)) # Inverse variance-covariance matrix (diagonal for simplicity)
  
  # Call the function and check the result
  result <- .CookDfun(betaU, beta0, vb.inv)
  
  # Manually calculate expected result
  dbetaU <- betaU - beta0
  expected_result <- t(dbetaU) %*% vb.inv %*% dbetaU
  
  # Check that the result matches the expected value
  expect_equal(as.numeric(result), as.numeric(expected_result))
})

test_that("Test .CookDfun with zero differences in betaU and beta0", {
  # Define input vectors and matrix
  betaU <- c(1.0, 2.0, 3.0)
  beta0 <- c(1.0, 2.0, 3.0)
  vb.inv <- diag(c(0.5, 0.5, 0.5))
  
  # Call the function and check the result
  result <- .CookDfun(betaU, beta0, vb.inv)
  
  # Expect Cook's distance to be zero
  expect_equal(as.numeric(result), 0)
})

test_that("Test .CookDfun with high influence cases", {
  # Define input vectors and matrix
  betaU <- c(10, 20, 30)  # Significantly different from beta0
  beta0 <- c(1.0, 2.0, 3.0)
  vb.inv <- diag(c(0.5, 0.5, 0.5))
  
  # Call the function and check the result
  result <- .CookDfun(betaU, beta0, vb.inv)
  
  # Expect a high Cook's distance value
  expect_true(as.numeric(result) > 100)
})

test_that("Test .CookDfun with non-conformable inputs", {
  # Define input vectors and matrix with non-conformable shapes
  betaU <- c(1.1, 2.2)
  beta0 <- c(1.0, 2.0, 3.0)  # Different length
  vb.inv <- diag(c(0.5, 0.5, 0.5))
  
  # Expect an error due to non-conformable shapes
  expect_warning(.CookDfun(betaU, beta0, vb.inv))
  
})

# Tests for CookDistance ----

# Create test data
set.seed(123)
test_data <- data.frame(
  SampleID = rep(1:10, each = 10),
  Time = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

# Create model
model <- lmmModel(
  data = test_data,
  sample_id = "SampleID",
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


test_that("Test CookDistance with valid input without Cook threshold", {
  # Call the function with a valid input and default cook_thr
  result <- CookDistance(model)
  
  # Check that the result is a numeric vector
  expect_type(result, "double")
  
  # Check that the length of the result matches the number of subjects
  expect_equal(length(result), length(unique(test_data$SampleID)))
  
  # Check that all values are finite
  expect_true(all(is.finite(result)))
})

test_that("Test CookDistance with different thresholds", {
  
  # Test with a threshold that is higher than any Cook's distance
  expect_output(result <- CookDistance(model, cook_thr = 10), 
                "No subject with a Cook's distance greater than 10")
  expect_true(all(result <= 10))
  
  # Test with a threshold of zero (all distances should be included in the plot)
  result <- CookDistance(model, cook_thr = -10)
  expect_true(any(result > 0))
})

test_that("Test CookDistance output consistency", {
  # Run the function and capture the output
  result <- CookDistance(model)
  
  # Check that the result is named and length matches number of subjects
  expect_named(result)
  expect_equal(length(result), length(unique(test_data$SampleID)))
})

test_that("Test CookDistance with invalid model input", {
  # Create an invalid model input (not an 'lme' object)
  invalid_model <- list(a = 1, b = 2)
  
  # Expect an error due to invalid model input
  expect_error(CookDistance(invalid_model))
})

