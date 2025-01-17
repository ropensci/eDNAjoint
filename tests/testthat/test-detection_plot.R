test_that("detection_plot input checks work", {
  testthat::skip_on_cran()
  #' @srrstats {G5.2,G5.2b} Tests the assure function input checks are behaving
  #'   as expected.
  # run joint model to do tests with
  model1 <- suppressWarnings({
    joint_model(data = goby_data, cov = c("Filter_time", "Salinity"),
                n_chain = 1, n_warmup = 25, n_iter = 100, multicore = FALSE)
  })

  model2 <- suppressWarnings({
    joint_model(data = green_crab_data, family = "negbin",
                n_chain = 1, n_warmup = 25, n_iter = 100, multicore = FALSE)
  })

  #1. make sure model fit is of class stanfit
  expect_error(detection_plot(as.matrix(model1$model), mu_min = 0.1,
                              mu_max = 1, cov_val = c(0, 0)),
               "model_fit must be of class 'stanfit'.")

  #2. make sure mu_min is a numeric value
  expect_error(detection_plot(model1$model, mu_min = "0.1", mu_max = 1,
                              cov_val = c(0, 0)),
               "mu_min must be a numeric value greater than 0")

  #3. make sure mu_min is a numeric value
  expect_error(detection_plot(model1$model, mu_min = 0, mu_max = 1,
                              cov_val = c(0, 0)),
               "mu_min must be a numeric value greater than 0")

  #4. make sure mu_max is a numeric value
  expect_error(detection_plot(model1$model, mu_min = 0.1, mu_max = "1",
                              cov_val = c(0, 0)),
               "mu_max must be a numeric value greater than mu_min")

  #5. make sure mu_max is a numeric value
  expect_error(detection_plot(model1$model, mu_min = 0.1, mu_max = 0.1,
                              cov_val = c(0, 0)),
               "mu_max must be a numeric value greater than mu_min")

  #6. make sure mu_max is a numeric value
  expect_error(detection_plot(model1$model, mu_min = 0.1, mu_max = 1,
                              cov_val = c(0, 0), probability = 1.1),
               "probability must be a numeric value between 0 and 1")

  #7. cov_val is numeric, if provided
  expect_error(detection_plot(model1$model, mu_min = 0.1, mu_max = 1,
                              cov_val = c(0, "0")),
               "cov_val must be a numeric vector")

  #8. Only include input cov_val if covariates are included in model
  expect_error(detection_plot(model2$model, mu_min = 0.1, mu_max = 1,
                              cov_val = c(0, 0)),
               paste0("cov_val must be NULL if the model does not contain ",
                      "site-level covariates."))

  #9. Input cov_val is the same length as the number of estimated covariates.
  expect_error(detection_plot(model1$model, mu_min = 0.1, mu_max = 1,
                              cov_val = c(0, 0, 0)),
               paste0("cov_val must be of the same length as the number of ",
                      "non-intercept site-level coefficients in the model."))

  #10. If covariates are in model, cov_val must be provided
  expect_error(detection_plot(model1$model, mu_min = 0.1, mu_max = 1),
               paste0("cov_val must be provided if the model contains ",
                      "site-level covariates."))

  #11. pcr_n must be an integer
  expect_error(detection_plot(model1$model, mu_min = 0.1, mu_max = 1,
                              cov_val = c(0, 0), pcr_n = 6.8),
               "pcr_n should be an integer.")
})
