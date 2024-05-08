test_that("muCritical input checks work", {
  #' @srrstats {G5.2,G5.2b} Tests the assure function input checks are
  #'   behaving as expected.

  # run joint model to do tests with
  model1 <- suppressWarnings({jointModel(data=gobyData,
                                         cov=c('Filter_time','Salinity'),
                                         multicore=FALSE, n.chain=1,
                                         n.iter.burn = 25,
                                         n.iter.sample = 75)})

  model2 <- suppressWarnings({traditionalModel(data=gobyData,multicore=FALSE,
                                               n.chain=1,n.iter.burn = 25,
                                               n.iter.sample = 75)})

  model3 <- suppressWarnings({jointModel(data=greencrabData,family='negbin',
                                         multicore=FALSE,
                                         n.chain=1,n.iter.burn = 25,
                                         n.iter.sample = 75)})

  #1. make sure model fit is of class stanfit
  expect_error(muCritical(as.matrix(model1$model), cov.val = c(0,0)),
               "modelfit must be of class 'stanfit'.")

  #2. make sure ci is valid
  expect_error(muCritical(model1$model, ci = 1, cov.val = c(0,0)),
               "ci must be a numeric value >0 and <1.")

  #3. make sure model fit contains p10 parameter
  expect_error(muCritical(model2$model),
               "modelfit must contain 'p10' parameter.")

  #4. if modelfit contains alpha, cov.val must be provided
  expect_error(muCritical(model1$model),
               paste0("If modelfit contains site-level covariates, values ",
                      "must be provided for cov.val"))

  #5. cov.val is numeric, if provided
  expect_error(muCritical(model1$model, cov.val=c(TRUE,TRUE)),
               "cov.val must be a numeric vector")

  #6. Only include input cov.val if covariates are included in model
  expect_error(muCritical(model3$model, cov.val=c(0,0)),
               paste0("cov.val must be NULL if the model does not ",
                      "contain site-level covariates."))

  #7. Only include input cov.val if covariates are included in model
  expect_error(muCritical(model1$model, cov.val=c(0,0,0)),
               paste0("cov.val must be of the same length as the number of ",
                      "non-intercept site-level coefficients in the model."))

})

test_that("muCritical output check", {

  # fit two models
  fit.cov = suppressWarnings({jointModel(data=gobyData,
                                         cov=c('Filter_time','Salinity'),
                                         family="poisson", p10priors=c(1,20),
                                         q=FALSE, multicore=FALSE, n.chain=1,
                                         n.iter.burn = 25,
                                         n.iter.sample = 75)})
  fit.q = suppressWarnings({jointModel(data=greencrabData, cov=NULL,
                                       family="negbin",p10priors=c(1,20),
                                       q=TRUE, multicore=FALSE, n.chain=1,
                                       n.iter.burn = 25,
                                       n.iter.sample = 75)})

  # check lengths of muCritical output
  expect_equal(length(muCritical(fit.cov$model, cov.val=c(0,0), ci = 0.9)),
               3)
  expect_equal(dim(muCritical(fit.q$model, ci = 0.9)),
               c(3,2))

  # check names of muCritical output
  expect_equal(names(muCritical(fit.cov$model, cov.val=c(0,0), ci = 0.9)),
               c('median','lower_ci','upper_ci'))
  expect_equal(rownames(muCritical(fit.q$model, ci = 0.9)),
               c('median','lower_ci','upper_ci'))
  expect_equal(names(muCritical(fit.q$model, ci = 0.9)),
               c('gear_1','gear_2'))

  # check values of muCritical output
  expect_equal(is.numeric(muCritical(fit.cov$model,
                                     cov.val=c(0,0), ci = 0.9)$median),
               TRUE)
  expect_equal(is.numeric(muCritical(fit.cov$model,
                                     cov.val=c(0,0),
                                     ci = 0.9)$lower_ci$CI_low),
               TRUE)
  expect_equal(is.numeric(muCritical(fit.cov$model,
                                     cov.val=c(0,0),
                                     ci = 0.9)$upper_ci$CI_high),
               TRUE)
  expect_equal(is.numeric(muCritical(fit.q$model, ci = 0.9)$gear_1),
               TRUE)
  expect_equal(is.numeric(muCritical(fit.q$model, ci = 0.9)$gear_2),
               TRUE)

})
