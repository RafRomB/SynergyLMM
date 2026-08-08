# Tests for simulateTumorGrowth ----

# simulateTumorGrowth() adds multiplicative rnorm() noise, so every test that
# compares simulated groups seeds the generator explicitly. Without this the
# results depend on how much randomness the previously run test files consumed.
set.seed(123)

test_that("simulateTumorGrowth returns a data frame with correct columns", {
  result <- simulateTumorGrowth(
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
  expect_equal(ncol(result), 4)  # subject, Treatment, Time, TumorVolume
  expect_true(all(c("subject", "Treatment", "Time", "TumorVolume") %in% colnames(result)))
})

test_that("simulateTumorGrowth returns correct number of rows for given npg and timepoints", {
  npg <- 5
  timepoints <- c(0, 3, 5, 10)
  result <- simulateTumorGrowth(
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

test_that("simulateTumorGrowth handles different growth rates correctly", {
  set.seed(123)
  result <- simulateTumorGrowth(
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
  last_time <- max(result$Time)
  last_time <- dplyr::filter(result, Time == last_time)
  means <- dplyr::summarize(last_time, .by = Treatment, mean_volume = mean(TumorVolume))
  
  expect_true(means$mean_volume[means$Treatment == "DrugA"] < means$mean_volume[means$Treatment == "Control"])
  expect_true(means$mean_volume[means$Treatment == "DrugB"] < means$mean_volume[means$Treatment == "DrugA"])
  expect_true(means$mean_volume[means$Treatment == "Combination"] < means$mean_volume[means$Treatment == "DrugB"])
})

test_that("simulateTumorGrowth handles different standard deviations correctly", {
  set.seed(123)
  result_sd_low <- simulateTumorGrowth(
    npg = 5,
    timepoints = c(0, 3, 5, 10),
    initial_volume = 100,
    grwrControl = 0.08,
    grwrA = 0.07,
    grwrB = 0.06,
    grwrComb = 0.04,
    sd = 0.01
  )
  
  set.seed(456)
  result_sd_high <- simulateTumorGrowth(
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

test_that("simulateTumorGrowth handles single timepoint correctly", {
  set.seed(123)
  result <- simulateTumorGrowth(
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


# Tests for warn_zero_pval ----

# These tests cover the exact wording of the warning emitted by lmmSynergy()
# when a simulated p-value is exactly 0. They are deterministic by design: the
# integration tests in test-lmmSynergy.R cannot assert the warning, because
# whether a Monte Carlo p-value hits exactly 0 depends on the BLAS/LAPACK
# implementation used for the eigendecomposition in MASS::mvrnorm().

test_that("warn_zero_pval warns only when a p-value is exactly 0", {
  expect_warning(warn_zero_pval(c(0.5, 0), nsim = 10),
                 "are approximated to 0", fixed = TRUE)
  expect_no_warning(warn_zero_pval(c(0.5, 0.1), nsim = 10))
  # NA p-values must not trigger the warning nor an error
  expect_no_warning(warn_zero_pval(c(0.5, NA), nsim = 10))
  expect_warning(warn_zero_pval(c(NA, 0), nsim = 10),
                 "are approximated to 0", fixed = TRUE)
})

test_that("warn_zero_pval returns whether it warned", {
  expect_true(suppressWarnings(warn_zero_pval(c(0.5, 0), nsim = 10)))
  expect_false(warn_zero_pval(c(0.5, 0.1), nsim = 10))
})

test_that("warn_zero_pval reports the resolution allowed by nsim", {
  expect_warning(
    warn_zero_pval(c(0.5, 0), nsim = 10),
    "p-values below p<1e-01 are approximated to 0. If you used method = 'RA' consider increasing 'nsim' value for more precise p-values.",
    fixed = TRUE
  )
  expect_warning(
    warn_zero_pval(c(0.5, 0), nsim = 10, gompertz = TRUE),
    "p-values below p<1e-01 are approximated to 0. If you used a Gompertz model, consider increasing 'nsim' value for more precise p-values.",
    fixed = TRUE
  )
  expect_warning(warn_zero_pval(c(0.5, 0), nsim = 1000), "p<1e-03", fixed = TRUE)
  expect_warning(warn_zero_pval(c(0.5, 0), nsim = 10000), "p<1e-04", fixed = TRUE)
})

test_that("approx_pval_label describes the smallest resolvable p-value", {
  expect_equal(approx_pval_label(10), "p<1e-01")
  expect_equal(approx_pval_label(1000), "p<1e-03")
  expect_equal(approx_pval_label(10000), "p<1e-04")
  # Non powers of ten are reported correctly rather than rounded to a power
  expect_equal(approx_pval_label(500), "p<2e-03")
  # Insensitive to the user's 'scipen' option
  op <- options(scipen = 1000)
  on.exit(options(op), add = TRUE)
  expect_equal(approx_pval_label(1000), "p<1e-03")
})
