test_that("jointModel input checks work", {
  #1. input tags are valid, q = FALSE, cov = FALSE
  expect_error(jointModel(data=list(qPCR.n=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.k=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)))),
               "Data should include 'qPCR.N', 'qPCR.K', and 'count'.")
  #2. input tags are valid, q = FALSE, cov = TRUE
  site.cov=rbind(c(1,0),c(0.4,-0.4))
  colnames(site.cov)=c('var_a','var_b')
  expect_error(jointModel(data=list(qPCR.n=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.k=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b')),
               "Data should include 'qPCR.N', 'qPCR.K', 'count', and 'site.cov'.")
  #3. input tags are valid, q = TRUE, cov = FALSE
  expect_error(jointModel(data=list(qPCR.n=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.k=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,2,1),c(1,2,NA))),
                          q=TRUE),
               "Data should include 'qPCR.N', 'qPCR.K', 'count', and 'count.type'.")
  #4. input tags are valid, q = TRUE, cov = TRUE
  site.cov=rbind(c(1,0),c(0.4,-0.4))
  colnames(site.cov)=c('var_a','var_b')
  expect_error(jointModel(data=list(qPCR.n=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.k=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,2,1),c(1,2,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b'),q=TRUE),
               "Data should include 'qPCR.N', 'qPCR.K', 'count', 'count.type', and 'site.cov'.")
})
