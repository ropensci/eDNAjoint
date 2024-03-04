test_that("jointModel input checks work", {
  #1. input tags are valid, q = FALSE, cov = FALSE
  expect_error(jointModel(data=list(qPCR.n=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.k=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)))),
               "Data should include 'qPCR.N', 'qPCR.K', and 'count'.")
  #2. input tags are valid, q = FALSE, cov = TRUE
  site.cov=cbind(c(1,0),c(0.4,-0.4))
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
  site.cov=cbind(c(1,0),c(0.4,-0.4))
  colnames(site.cov)=c('var_a','var_b')
  expect_error(jointModel(data=list(qPCR.n=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.k=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,2,1),c(1,2,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b'),q=TRUE),
               "Data should include 'qPCR.N', 'qPCR.K', 'count', 'count.type', and 'site.cov'.")

  #5. make sure dimensions of qPCR.N and qPCR.K are equal
  #' @srrstats {BS2.1a} Test to ensure pre-processing routines to ensure all input data is dimensionally commensurate
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1,1),c(1,1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)))),
               "Dimensions of qPCR.N and qPCR.K do not match.")

  #6. make sure dimensions of count and count.type are equal, if count.type is present
  #' @srrstats {BS2.1a} Test to ensure pre-processing routines to ensure all input data is dimensionally commensurate
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,2),c(1,2))),
                          q=TRUE),
               "Dimensions of count and count.type do not match.")

  #7. make sure number of rows in count = number of rows in qPCR.N and qPCR.K
  #' @srrstats {BS2.1a} Test to ensure pre-processing routines to ensure all input data is dimensionally commensurate
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA),c(4,1,1)))),
               "Number of sites \\(rows\\) in qPCR data and traditional survey count data do not match.")

  #8. make sure all data is numeric -- if q == TRUE
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c('NA',2,2),c(1,2,2))),
                          q=TRUE),
               "Data should be numeric.")

  #9. make sure all data is numeric -- if q == FALSE
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,'NA')),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)))),
               "Data should be numeric.")

  #10. make sure locations of NAs in count data match locations of NAs in count.type data
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(NA,2,2),c(1,2,2))),
                          q=TRUE),
               "Empty data cells \\(NA\\) in count data should match empty data cells \\(NA\\) in count.type data.")

  #11. make sure locations of NAs in qPCR.N data match locations of NAs in qPCR.K data
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,NA,1)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)))),
               "Empty data cells \\(NA\\) in qPCR.N data should match empty data cells \\(NA\\) in qPCR.K data.")

  #12. make sure family is either 'poisson', 'negbin', or 'gamma'
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA))),
                          family='normal'),
               "Invalid family. Options include 'poisson', 'negbin', and 'gamma'.")

  #13. p10 priors is a vector of two numeric values
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA))),
                          p10priors=c(1,1,2)),
               "p10priors should be a vector of two positive numeric values. ex. c\\(1,20\\)")

  #14. p10 priors is a vector of two numeric values
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA))),
                          p10priors='1,20'),
               "p10priors should be a vector of two positive numeric values. ex. c\\(1,20\\)")

  #15. p10 priors is a vector of two numeric values
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA))),
                          p10priors=c(0,20)),
               "p10priors should be a vector of two positive numeric values. ex. c\\(1,20\\)")

  #16. phi priors is a vector of two numeric values
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA))),
                          phipriors=c(0,1)),
               "phipriors should be a vector of two positive numeric values. ex. c\\(0.25,0.25\\)")

  #17. the smallest count.type is 1
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(0,1,2),c(1,2,NA))),
                          q=TRUE),
               "The first gear type should be referenced as 1 in count.type. Subsequent gear types should be referenced 2, 3, 4, etc.")

  #18. count are integers
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4.1,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,1,2),c(1,2,NA))),
                          q=TRUE, family = 'negbin'),
               "All values in count should be non-negative integers. Use family = 'gamma' if count is continuous.")

  #19. qPCR.N are integers
  expect_error(jointModel(data=list(qPCR.N=rbind(c(0.99,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,1,2),c(1,2,NA))),
                          q=TRUE),
               "All values in qPCR.N should be non-negative integers.")

  #20. qPCR.K are integers
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3.1,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,1,2),c(1,2,NA))),
                          q=TRUE),
               "All values in qPCR.K should be non-negative integers.")

  #21. count.type are integers
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1.1,1,2),c(1,2,NA))),
                          q=TRUE),
               "All values in count.type should be integers.")

  #22. site.cov is numeric, if present
  site.cov=cbind(c('high','low'),c(0.4,-0.4))
  colnames(site.cov)=c('var_a','var_b')
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,2,1),c(1,2,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b'),q=TRUE),
               "site.cov should be numeric.")

  #23. cov values match column names in site.cov
  site.cov=cbind(c(0,1),c(0.4,-0.4))
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,2,1),c(1,2,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b'),q=TRUE),
               "cov values should be listed in the column names of site.cov in the data.")

  #24. site.cov has same number of rows as qPCR.N and count, if present
  site.cov=cbind(c(0,1,1),c(0.4,-0.4,1))
  colnames(site.cov)=c('var_a','var_b')
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,2,1),c(1,2,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b'),q=TRUE),
               "The number of rows in site.cov matrix should match the number of rows in all other matrices.")

})
