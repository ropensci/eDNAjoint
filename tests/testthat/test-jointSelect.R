test_that("jointSelect input checks work", {

  # run joint model and traditional model to do tests with
  data <- data("greencrabData")

  model1 <- traditionalModel(data=greencrabData, family='negbin',multicore=FALSE,
                            n.chain=2,n.iter.burn=0,n.iter.sample=200)
  model2 <- jointModel(data=greencrabData, family='negbin',multicore=FALSE,
                      n.chain=2,n.iter.burn=0,n.iter.sample=200)


  # 1. Check that model inputs are a list
  expect_error(jointSelect(c(as.matrix(model1),as.matrix(model2))),
               "modelfits must be a list.")

  # 2. Check that all data objects are of class stanfit
  expect_error(jointSelect(list(as.matrix(model1),as.matrix(model2))),
               "Model fits in modelfits input must be of class 'stanfit'.")

  # 3. Check that all models are of the same type
  expect_error(jointSelect(list(model1,model2)),
               "All modelfits must be fit with either jointModel\\(\\) or all with traditionalModel\\(\\).")
})

