test_that("jointSelect input checks work", {
  #' @srrstats {G5.2,G5.2b} Tests the assure function input checks are behaving
  #' as expected.

  # run joint model and traditional model to do tests with
  data <- data("greencrabData")

  model1 <- traditionalModel(data=greencrabData, family='negbin',
                             multicore=FALSE0)
  model2 <- jointModel(data=greencrabData, family='negbin',multicore=FALSE)


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
