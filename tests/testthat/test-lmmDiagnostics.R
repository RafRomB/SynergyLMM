
# Example data and model for testing ----
set.seed(123)
test_data <- data.frame(
  Mouse = rep(1:10, each = 10),
  Day = rep(0:9, times = 10),
  Treatment = rep(c("Control", "Drug_A", "Drug_B", "Drug_AB"), each = 10, length.out = 100),
  TV = rnorm(100, mean = 100, sd = 20)
)

model <- lmmModel(
  data = test_data,
  mouse_id = "Mouse",
  day = "Day",
  treatment = "Treatment",
  tumor_vol = "TV",
  trt_control = "Control",
  drug_a = "Drug_A",
  drug_b = "Drug_B",
  drug_ab = "Drug_AB",
  day_start = 0,
  min_observations = 1,
  show_plot = FALSE
)

# Tests for ranefDiagnostics function ----

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
  expect_s3_class(diagnostics$Levene.test$conditional_resid, "anova")
  expect_s3_class(diagnostics$Levene.test$marginal_resid, "anova")
  expect_s3_class(diagnostics$Levene.test$normalized_resid, "anova")
  
  # Check the class of Fligner-Killeen test results
  expect_s3_class(diagnostics$Fligner.test$conditional_resid, "htest")
  expect_s3_class(diagnostics$Fligner.test$marginal_resid, "htest")
  expect_s3_class(diagnostics$Fligner.test$normalized_resid, "htest")
})

test_that("ranefDiagnostics handles small sample sizes for D'Agostino test", {
  # Modify test data to have fewer than 20 samples for ranef
  small_sample_data <- test_data[test_data$Mouse <= 2, ]
  
  small_model <- lmmModel(
    data = test_data,
    mouse_id = "Mouse",
    day = "Day",
    treatment = "Treatment",
    tumor_vol = "TV",
    trt_control = "Control",
    drug_a = "Drug_A",
    drug_b = "Drug_B",
    drug_ab = "Drug_AB",
    day_start = 0,
    min_observations = 1,
    show_plot = FALSE
  )
  
  diagnostics <- ranefDiagnostics(small_model)
  
  # Check if D'Agostino test was skipped due to small sample size
  expect_true(is.character(diagnostics$Normality$DAgostino.test))
  expect_match(diagnostics$Normality$DAgostino.test, "Sample size must be at least 20")
})
