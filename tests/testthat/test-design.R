small_design <- function(k1 = 14) {
  parallel_crt(icc = 0.05, cac = 0.8, delta = 0.25, k1 = k1, m1 = 30,
               k2 = 0:1, m2 = c(20, 30), rho = 30)
}

planned_design <- function() {
  suppressWarnings(adaptive_analysis(small_design(), target_power = 0.7, tol = 0.02,
                                     method = "lambda", verbose = FALSE))
}

test_that("interim_sensitivity works with its default z1 range", {
  skip_on_cran()
  sens <- interim_sensitivity(planned_design(),
                              scenarios = list(low = list(icc = 0.03)), n_z1 = 40)
  expect_s3_class(sens, "interim_sensitivity")
  expect_true(all(is.finite(sens$z1_grid)))
})

test_that("summary of explored results shows the stage 1 parameters", {
  skip_on_cran()
  explored <- suppressWarnings(
    adaptive_analysis(small_design(k1 = 13:14), target_power = 0.7, tol = 0.02,
                      method = "lambda", verbose = FALSE)
  )
  out <- NULL
  utils::capture.output(out <- summary(explored))
  expect_s3_class(out, "data.frame")
  printed <- utils::capture.output(summary(explored))
  expect_true(any(grepl("k1", printed)))
  expect_true(any(grepl("m1", printed)))
})

test_that("compare_to_fixed reports a valid comparison", {
  skip_on_cran()
  comparison <- suppressWarnings(compare_to_fixed(planned_design(), verbose = FALSE))
  expect_true(comparison$valid)
  expect_false(any(grepl("WARNING", utils::capture.output(print(comparison)))))
})

test_that("a stored parallel model gives the same information as a fresh one", {
  skip_on_cran()
  fixed <- c(list(icc = 0.05, cac = 0.8, delta = 0.25, family = "gaussian",
                  baseline = 0, k1 = 14, m1 = 30, cache = new.env()))
  first <- parallel_model_builder(list(k2 = 1, m2 = 20), fixed)
  changed <- parallel_model_builder(list(k2 = 1, m2 = 30), fixed)
  fresh <- fixed
  fresh$cache <- new.env()
  expect_equal(changed$I_eff, parallel_model_builder(list(k2 = 1, m2 = 30), fresh)$I_eff)
  expect_false(isTRUE(all.equal(first$I_eff, changed$I_eff)))
})

test_that("the probability reshaping keeps every row and value", {
  df <- data.frame(k1 = 1:2, prob_efficacy = c(0.1, 0.2), prob_futility = c(0.3, 0.4),
                   prob_continue = c(0.6, 0.4))
  long <- acrt:::.stack_probabilities(df)
  expect_equal(nrow(long), 6)
  expect_equal(long$probability[long$outcome == "futility" & long$k1 == 2], 0.4)
})
