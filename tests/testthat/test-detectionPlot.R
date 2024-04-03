test_that("detectionPlot input checks work", {
  #' @srrstats {G5.2,G5.2b} Tests the assure function input checks are behaving
  #'   as expected.
  # run joint model to do tests with
  model1 <- jointModel(data=gobyData, cov=c('Filter_time','Salinity'),
                       multicore=FALSE)

  model2 <- jointModel(data=greencrabData,family='negbin',multicore=FALSE)

  #1. make sure model fit is of class stanfit
  expect_error(detectionPlot(as.matrix(model1$model), mu.min = 0.1,
                             mu.max = 1, cov.val = c(0,0)),
               "modelfit must be of class 'stanfit'.")

  #2. make sure mu.min is a numeric value
  expect_error(detectionPlot(model1$model, mu.min = '0.1', mu.max = 1,
                             cov.val = c(0,0)),
               "mu.min must be a numeric value greater than 0")

  #3. make sure mu.min is a numeric value
  expect_error(detectionPlot(model1$model, mu.min = 0, mu.max = 1,
                             cov.val = c(0,0)),
               "mu.min must be a numeric value greater than 0")

  #4. make sure mu.max is a numeric value
  expect_error(detectionPlot(model1$model, mu.min = 0.1, mu.max = '1',
                             cov.val = c(0,0)),
               "mu.max must be a numeric value greater than mu.min")

  #5. make sure mu.max is a numeric value
  expect_error(detectionPlot(model1$model, mu.min = 0.1, mu.max = 0.1,
                             cov.val = c(0,0)),
               "mu.max must be a numeric value greater than mu.min")

  #6. make sure mu.max is a numeric value
  expect_error(detectionPlot(model1$model, mu.min = 0.1, mu.max = 1,
                             cov.val = c(0,0), probability = 1.1),
               "probability must be a numeric value between 0 and 1")

  #7. cov.val is numeric, if provided
  expect_error(detectionPlot(model1$model, mu.min = 0.1, mu.max = 1,
                             cov.val = c(0,'0')),
               "cov.val must be a numeric vector")

  #8. Only include input cov.val if covariates are included in model
  expect_error(detectionPlot(model2$model, mu.min = 0.1, mu.max = 1,
                             cov.val = c(0,0)),
               paste0("cov.val must be 'None' if the model does not contain ",
                      "site-level covariates."))

  #9. Input cov.val is the same length as the number of estimated covariates.
  expect_error(detectionPlot(model1$model, mu.min = 0.1, mu.max = 1,
                             cov.val = c(0,0,0)),
               paste0("cov.val must be of the same length as the number of ",
                      "non-intercept site-level coefficients in the model."))

  #10. If covariates are in model, cov.val must be provided
  expect_error(detectionPlot(model1$model, mu.min = 0.1, mu.max = 1),
               paste0("cov.val must be provided if the model contains ",
                      "site-level covariates."))

  #11. qPCR.N must be an integer
  expect_error(detectionPlot(model1$model, mu.min = 0.1, mu.max = 1,
                             cov.val = c(0,0), qPCR.N = 6.8),
               "qPCR.N should be an integer.")
})


