test_that("summarize input checks work", {

  #create model output
  test <- list(qPCR.N=greencrabData$qPCR.N[1:10,],
               qPCR.K=greencrabData$qPCR.K[1:10,],
               count=greencrabData$count[1:10,],
               count.type=greencrabData$count.type[1:10,])
  out <- jointModel(data=test,q=TRUE,family='negbin',n.chain=1,n.iter.burn=500,
                    n.iter.sample=500)

  #1. make sure model fit is of class stanfit
  data <- data.frame(y=c(1,2,3),x=c(1,2,3))
  lm_out <- lm(y ~ x, data=data)
  expect_error(summarize(lm_out),
               "modelfit must be of class 'stanfit'.")

  #2. make sure probs is a numeric vector
  expect_error(summarize(out,probs=c('95%')),
               "probs must be a numeric vector.")

  #3. make sure all values of probs are between 0 and 1
  expect_error(summarize(out,probs=c(5,95)),
               "probs must be between 0 and 1.")

  #4. make sure par is a character vector
  expect_error(summarize(out,par=c(1,2,3)),
               "par must be a character vector.")

  #5. make sure model fit contains all par input
  expect_error(summarize(out,par='alpha'),
               "modelfit must contain all selected parameters: alpha")
})
