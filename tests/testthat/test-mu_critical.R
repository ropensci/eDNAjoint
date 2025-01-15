test_that("mu_critical input checks work", {
  testthat::skip_on_cran()
  #' @srrstats {G5.2,G5.2b} Tests the assure function input checks are
  #'   behaving as expected.

  # run joint model to do tests with
  model1 <- suppressWarnings({
    joint_model(data = goby_data, cov = c("Filter_time", "Salinity"),
                multicore = FALSE, n_chain = 1,
                n_warmup = 25, n_iter = 100)
  })

  model2 <- suppressWarnings({
    traditional_model(data = goby_data, multicore = FALSE,
                      n_chain = 1, n_warmup = 25, n_iter = 100)
  })

  model3 <- suppressWarnings({
    joint_model(data = green_crab_data, family = "negbin", multicore = FALSE,
                n_chain = 1, n_warmup = 25, n_iter = 100)
  })

  #1. make sure model fit is of class stanfit
  expect_error(mu_critical(as.matrix(model1$model), cov_val = c(0, 0)),
               "model_fit must be of class 'stanfit'.")

  #2. make sure ci is valid
  expect_error(mu_critical(model1$model, ci = 1, cov_val = c(0, 0)),
               "ci must be a numeric value >0 and <1.")

  #3. make sure model fit contains p10 parameter
  expect_error(mu_critical(model2$model),
               "model_fit must contain 'p10' parameter.")

  #4. if modelfit contains alpha, cov_val must be provided
  expect_error(mu_critical(model1$model),
               paste0("If model_fit contains site-level covariates, values ",
                      "must be provided for cov_val"))

  #5. cov_val is numeric, if provided
  expect_error(mu_critical(model1$model, cov_val = c(TRUE, TRUE)),
               "cov_val must be a numeric vector")

  #6. Only include input cov_val if covariates are included in model
  expect_error(mu_critical(model3$model, cov_val = c(0, 0)),
               paste0("cov_val must be NULL if the model does not ",
                      "contain site-level covariates."))

  #7. Only include input cov_val if covariates are included in model
  expect_error(mu_critical(model1$model, cov_val = c(0, 0, 0)),
               paste0("cov_val must be of the same length as the number of ",
                      "non-intercept site-level coefficients in the model."))

})

test_that("mu_critical output check", {
  testthat::skip_on_cran()

  # fit two models
  fit_cov <- suppressWarnings({
    joint_model(data = goby_data, cov = c("Filter_time", "Salinity"),
                family = "poisson", p10_priors = c(1, 20),
                q = FALSE, multicore = FALSE, n_chain = 1,
                n_warmup = 25, n_iter = 100)
  })

  fit_q <- suppressWarnings({
    joint_model(data = green_crab_data, cov = NULL,
                family = "negbin", p10_priors = c(1, 20),
                q = TRUE, multicore = FALSE, n_chain = 1,
                n_warmup = 25, n_iter = 100)
  })

  # check lengths of mu_critical output
  expect_equal(length(mu_critical(fit_cov$model, cov_val = c(0, 0), ci = 0.9)),
               3)
  expect_equal(dim(mu_critical(fit_q$model, ci = 0.9)),
               c(3, 2))

  # check names of mu_critical output
  expect_equal(names(mu_critical(fit_cov$model, cov_val = c(0, 0), ci = 0.9)),
               c("median", "lower_ci", "upper_ci"))
  expect_equal(rownames(mu_critical(fit_q$model, ci = 0.9)),
               c("median", "lower_ci", "upper_ci"))
  expect_equal(names(mu_critical(fit_q$model, ci = 0.9)),
               c("gear_1", "gear_2"))

  # check values of mu_critical output
  expect_equal(is.numeric(mu_critical(fit_cov$model,
                                      cov_val = c(0, 0), ci = 0.9)$median),
               TRUE)
  expect_equal(is.numeric(mu_critical(fit_cov$model,
                                      cov_val = c(0, 0),
                                      ci = 0.9)$lower_ci$CI_low),
               TRUE)
  expect_equal(is.numeric(mu_critical(fit_cov$model,
                                      cov_val = c(0, 0),
                                      ci = 0.9)$upper_ci$CI_high),
               TRUE)
  expect_equal(is.numeric(mu_critical(fit_q$model, ci = 0.9)$gear_1),
               TRUE)
  expect_equal(is.numeric(mu_critical(fit_q$model, ci = 0.9)$gear_2),
               TRUE)

})
