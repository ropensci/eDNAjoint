test_that("detection_calculate input checks work", {
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
  expect_error(detection_calculate(as.matrix(model1$model), mu = c(0.1, 0.5),
                                   cov_val = c(0, 0)),
               "model_fit must be of class 'stanfit'.")

  #2. make sure ci is valid
  expect_error(detection_calculate(model1$model, mu = c("0", 0.5),
                                   cov_val = c(0, 0)),
               "mu must be a numeric vector of positive values")

  #3. make sure mu is a numeric vector of positive values
  expect_error(detection_calculate(model1$model, mu = c(0, 0.5),
                                   cov_val = c(0, 0)),
               "mu must be a numeric vector of positive values")

  #4. make sure probability is a numeric value between 0 and 1
  expect_error(detection_calculate(model1$model, mu = c(0.1, 0.5),
                                   cov_val = c(0, 0), probability = 1.05),
               "probability must be a numeric value between 0 and 1")

  #5. cov_val is numeric, if provided
  expect_error(detection_calculate(model1$model, mu = c(0.1, 0.5),
                                   cov_val = c("0", 0)),
               "cov_val must be a numeric vector")

  #6. Only include input cov_val if covariates are included in model
  expect_error(detection_calculate(model2$model, mu = c(0.1, 0.5),
                                   cov_val = c(0, 0)),
               paste0("cov_val must be NULL if the model does not ",
                      "contain site-level covariates."))

  #7. Input cov_val is the same length as the number of estimated covariates.
  expect_error(detection_calculate(model1$model, mu = c(0.1, 0.5),
                                   cov_val = c(0, 0, 0)),
               paste0("cov_val must be of the same length as the number of ",
                      "non-intercept site-level coefficients in the model."))

  #8. If covariates are in model, cov_val must be provided
  expect_error(detection_calculate(model1$model, mu = c(0.1, 0.5)),
               paste0("cov_val must be provided if the model contains ",
                      "site-level covariates."))

  #9. pcr_n must be an integer
  expect_error(detection_calculate(model1$model, mu = c(0.1, 0.5),
                                   cov_val = c(0, 0), pcr_n = 6.8),
               "pcr_n should be an integer.")
})
