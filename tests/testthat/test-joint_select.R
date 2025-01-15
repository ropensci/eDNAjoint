test_that("joint_select input checks work", {
  testthat::skip_on_cran()
  #' @srrstats {G5.2,G5.2b} Tests the assure function input checks are behaving
  #'   as expected.

  # run joint model and traditional model to do tests with

  model1 <- suppressWarnings({
    traditional_model(data = green_crab_data, family = "negbin",
                      multicore = FALSE, n_chain = 1,
                      n_warmup = 25, n_iter = 100)
  })

  model2 <- suppressWarnings({
    joint_model(data = green_crab_data, family = "negbin",
                multicore = FALSE, n_chain = 1,
                n_warmup = 25, n_iter = 100)
  })


  # 1. Check that model inputs are a list
  expect_error(joint_select(c(as.matrix(model1$model),
                              as.matrix(model2$model))),
               "model_fits must be a list.")

  # 2. Check that all data objects are of class stanfit
  expect_error(joint_select(list(as.matrix(model1$model),
                                 as.matrix(model2$model))),
               "Model fits in model_fits input must be of class 'stanfit'.")

  # 3. Check that all models are of the same type
  expect_error(joint_select(list(model1$model, model2$model)),
               paste0("All model_fits must be fit with either ",
                      "joint_model\\(\\) or all with traditional_model\\(\\)."))
})

test_that("joint_select output check", {
  testthat::skip_on_cran()
  # fit models
  fit_q1 <- joint_model(data = green_crab_data, family = "negbin",
                        p10_priors = c(1, 20), q = TRUE, multicore = FALSE,
                        n_chain = 1, n_iter = 1000)
  fit_q2 <- joint_model(data = green_crab_data, family = "negbin",
                        p10_priors = c(1, 50), q = TRUE, multicore = FALSE,
                        n_chain = 1, n_iter = 1000)

  # run select
  select <- joint_select(model_fits = list(fit_q1$model, fit_q2$model))

  # check dimensions
  expect_equal(dim(select)[1], 2)

  # check numeric
  expect_equal(is.numeric(select[, 1]), TRUE)
  expect_equal(is.numeric(select[, 2]), TRUE)

})
