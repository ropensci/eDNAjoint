test_that("jointSelect input checks work", {
  #' @srrstats {G5.2,G5.2b} Tests the assure function input checks are behaving
  #'   as expected.

  # run joint model and traditional model to do tests with
  data <- data("greencrabData")

  model1 <- suppressWarnings({traditionalModel(data=greencrabData,
                                               family='negbin',
                             multicore=FALSE, n.chain=1,n.iter.burn = 25,
                             n.iter.sample = 75)})
  model2 <- suppressWarnings({jointModel(data=greencrabData, family='negbin',
                                         multicore=FALSE,n.chain=1,
                                         n.iter.burn = 25,
                                         n.iter.sample = 75,)})


  # 1. Check that model inputs are a list
  expect_error(jointSelect(c(as.matrix(model1$model),as.matrix(model2$model))),
               "modelfits must be a list.")

  # 2. Check that all data objects are of class stanfit
  expect_error(jointSelect(list(as.matrix(model1$model),
                                as.matrix(model2$model))),
               "Model fits in modelfits input must be of class 'stanfit'.")

  # 3. Check that all models are of the same type
  expect_error(jointSelect(list(model1$model,model2$model)),
               paste0("All modelfits must be fit with either jointModel\\(\\) ",
                      "or all with traditionalModel\\(\\)."))
})

test_that("jointSelect output check", {
  # fit models
  fit.q1 = jointModel(data=greencrabData, family="negbin",
                     p10priors=c(1,20), q=TRUE, multicore=FALSE,
                     n.chain=1, n.iter.sample = 1000)
  fit.q2 = jointModel(data=greencrabData, family="negbin",
                      p10priors=c(1,50), q=TRUE, multicore=FALSE,
                      n.chain=1, n.iter.sample = 1000)

  # run select
  select <- jointSelect(modelfits=list(fit.q1$model, fit.q2$model))

  # check dimensions
  expect_equal(dim(select)[1],2)

  # check numeric
  expect_equal(is.numeric(select[,1]),TRUE)
  expect_equal(is.numeric(select[,2]),TRUE)

})
