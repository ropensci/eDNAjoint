test_that("jointSummarize input checks work", {
  #' @srrstats {G5.2,G5.2b} Tests the assure function input checks are
  #' behaving as expected.

  # run traditional model to do tests with
  data <- data("greencrabData")

  out <- traditionalModel(data=greencrabData, family='negbin',multicore=FALSE,
                          n.chain=2,n.iter.burn=50,n.iter.sample=200)

  #1. make sure model fit is of class stanfit
  data <- data.frame(y=c(1,2,3),x=c(1,2,3))
  lm_out <- lm(y ~ x, data=data)
  expect_error(jointSummarize(lm_out$model),
               "modelfit must be of class 'stanfit'.")

  #2. make sure probs is a numeric vector
  expect_error(jointSummarize(out$model,probs=c('95%')),
               "probs must be a numeric vector.")

  #3. make sure all values of probs are between 0 and 1
  expect_error(jointSummarize(out$model,probs=c(5,95)),
               "probs must be between 0 and 1.")

  #4. make sure par is a character vector
  expect_error(jointSummarize(out$model,par=c(1,2,3)),
               "par must be a character vector.")

  #5. make sure model fit contains all par input
  expect_error(jointSummarize(out$model,par='alpha'),
               "modelfit must contain all selected parameters: alpha")
})


