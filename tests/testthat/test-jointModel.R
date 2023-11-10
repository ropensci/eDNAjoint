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
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1,1),c(1,1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)))),
               "Dimensions of qPCR.N and qPCR.K do not match.")

  #6. make sure dimensions of count and count.type are equal, if count.type is present
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,2),c(1,2))),
                          q=TRUE),
               "Dimensions of count and count.type do not match.")

  #7. make sure number of rows in count = number of rows in qPCR.N and qPCR.K
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
               "Data should be numeric \\(i.e. contain integers or NA\\).")

  #9. make sure all data is numeric -- if q == FALSE
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,'NA')),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)))),
               "Data should be numeric \\(i.e. contain integers or NA\\).")

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

  #12. make sure family is either 'poisson' or 'negbin'
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA))),
                          family='normal'),
               "Invalid family. Options include 'poisson' and 'negbin'.")

  #13. p10 priors is a vector of two integers
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA))),
                          p10priors=c(1,1,2)),
               "p10priors should be a vector of two integers. ex. c\\(1,20\\)")

  #14. p10 priors is a vector of two integers
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA))),
                          p10priors='1,20'),
               "p10priors should be a vector of two integers. ex. c\\(1,20\\)")

  #15. the smallest count.type is 1
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(0,1,2),c(1,2,NA))),
                          q=TRUE),
               "The first gear type should be referenced as 1 in count.type. Subsequent gear types should be referenced 2, 3, 4, etc.")

  #16. count are integers
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4.1,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,1,2),c(1,2,NA))),
                          q=TRUE),
               "All values in count should be integers.")

  #17. qPCR.N are integers
  expect_error(jointModel(data=list(qPCR.N=rbind(c(0.99,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,1,2),c(1,2,NA))),
                          q=TRUE),
               "All values in qPCR.N should be integers.")

  #18. qPCR.K are integers
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3.1,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,1,2),c(1,2,NA))),
                          q=TRUE),
               "All values in qPCR.K should be integers.")

  #19. count.type are integers
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1.1,1,2),c(1,2,NA))),
                          q=TRUE),
               "All values in count.type should be integers.")

  #20. site.cov is numeric, if present
  site.cov=cbind(c('high','low'),c(0.4,-0.4))
  colnames(site.cov)=c('var_a','var_b')
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,2,1),c(1,2,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b'),q=TRUE),
               "site.cov should be numeric.")

  #21. cov values match column names in site.cov
  site.cov=cbind(c(0,1),c(0.4,-0.4))
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,2,1),c(1,2,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b'),q=TRUE),
               "cov values should be listed in the column names of site.cov in the data.")

  #22. site.cov has same number of rows as qPCR.N and count, if present
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

test_that("joint_catchability_pois works", {

  test <- list(qPCR.N=greencrabData$qPCR.N[1:10,],
               qPCR.K=greencrabData$qPCR.K[1:10,],
               count=greencrabData$count[1:10,],
               count.type=greencrabData$count.type[1:10,])

  out <- jointModel(data=test,q=TRUE,n.chain=1,n.iter.burn=500,
                    n.iter.sample=500)
  expect_lt(sum(colMeans(rstan::extract(out,par='log_lik')$log_lik)), 0)

})

test_that("joint_catchability_negbin works", {

  test <- list(qPCR.N=greencrabData$qPCR.N[1:10,],
               qPCR.K=greencrabData$qPCR.K[1:10,],
               count=greencrabData$count[1:10,],
               count.type=greencrabData$count.type[1:10,])

  out <- jointModel(data=test,q=TRUE,family='negbin',n.chain=1,n.iter.burn=500,
                    n.iter.sample=500)
  expect_lt(sum(colMeans(rstan::extract(out,par='log_lik')$log_lik)), 0)

})

test_that("joint_pois works", {

  test <- list(qPCR.N=greencrabData$qPCR.N[1:10,],
               qPCR.K=greencrabData$qPCR.K[1:10,],
               count=greencrabData$count[1:10,],
               count.type=greencrabData$count.type[1:10,])

  out <- jointModel(data=test,q=FALSE,n.chain=1,n.iter.burn=500,
                    n.iter.sample=500)
  expect_lt(sum(colMeans(rstan::extract(out,par='log_lik')$log_lik)), 0)

})

test_that("joint_negbin works", {

  test <- list(qPCR.N=greencrabData$qPCR.N[1:10,],
               qPCR.K=greencrabData$qPCR.K[1:10,],
               count=greencrabData$count[1:10,],
               count.type=greencrabData$count.type[1:10,])

  out <- jointModel(data=test,q=FALSE,family='negbin',n.chain=1,n.iter.burn=500,
                    n.iter.sample=500)
  expect_lt(sum(colMeans(rstan::extract(out,par='log_lik')$log_lik)), 0)

})

test_that("joint_cov_negbin works", {

  test <- list(qPCR.N=gobyData$qPCR.N[1:10,],
               qPCR.K=gobyData$qPCR.K[1:10,],
               count=gobyData$count[1:10,],
               site.cov=gobyData$site.cov[1:10,])

  out <- jointModel(data=test,q=FALSE,family='negbin',
                    cov=c('Filter_time','Salinity'),
                    n.chain=1,n.iter.burn=500,
                    n.iter.sample=500)
  expect_lt(sum(colMeans(rstan::extract(out,par='log_lik')$log_lik)), 0)

})

test_that("joint_cov_pois works", {

  test <- list(qPCR.N=gobyData$qPCR.N[1:10,],
               qPCR.K=gobyData$qPCR.K[1:10,],
               count=gobyData$count[1:10,],
               site.cov=gobyData$site.cov[1:10,])

  out <- jointModel(data=test,q=FALSE,
                    cov=c('Filter_time','Salinity'),
                    n.chain=1,n.iter.burn=500,
                    n.iter.sample=500)
  expect_lt(sum(colMeans(rstan::extract(out,par='log_lik')$log_lik)), 0)

})

test_that("joint_cov_catchability_negbin works", {

  count <- gobyData$count[1:10,]
  count.type <- matrix(NA,10,22)
  for(i in 1:nrow(count.type)){
      n <- length(count[i,][!is.na(count[i,])])
      n_1 <- round(n/2)
      n_2 <- n - n_1
      count.type[i,1:n] <- c(rep(1,n_1),rep(2,n_2))
  }

  test <- list(qPCR.N=gobyData$qPCR.N[1:10,],
               qPCR.K=gobyData$qPCR.K[1:10,],
               count=count,
               count.type=count.type,
               site.cov=gobyData$site.cov[1:10,])

  out <- jointModel(data=test,q=TRUE,family='negbin',
                    cov=c('Filter_time','Salinity'),
                    n.chain=1,n.iter.burn=500,
                    n.iter.sample=500)
  expect_lt(sum(colMeans(rstan::extract(out,par='log_lik')$log_lik)), 0)

})

test_that("joint_cov_catchability_pois works", {

  count <- gobyData$count[1:10,]
  count.type <- matrix(NA,10,22)
  for(i in 1:nrow(count.type)){
    n <- length(count[i,][!is.na(count[i,])])
    n_1 <- round(n/2)
    n_2 <- n - n_1
    count.type[i,1:n] <- c(rep(1,n_1),rep(2,n_2))
  }

  test <- list(qPCR.N=gobyData$qPCR.N[1:10,],
               qPCR.K=gobyData$qPCR.K[1:10,],
               count=count,
               count.type=count.type,
               site.cov=gobyData$site.cov[1:10,])

  out <- jointModel(data=test,q=TRUE,
                    cov=c('Filter_time','Salinity'),
                    n.chain=1,n.iter.burn=500,
                    n.iter.sample=500)
  expect_lt(sum(colMeans(rstan::extract(out,par='log_lik')$log_lik)), 0)

})
