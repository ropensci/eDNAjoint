test_that("jointModel input checks work", {
  #' @srrstats {G5.2,G5.2b,BS2.15} Tests the assure function input checks are behaving as expected.
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

  #25. make sure count.type is not zero-length
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=matrix(NA,ncol=3,nrow=0)),
                          q=TRUE),
               "count.type contains zero-length data.")

  #26. make sure no column is entirely NA in count.type
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(4,1,NA),c(1,1,NA))),
                          q=TRUE),
               "count.type contains a column with all NA.")

  #27. make sure no column is entirely NA in qPCR.N
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,NA),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)))),
               "qPCR.N contains a column with all NA.")

  #28. make sure no column is entirely NA in qPCR.K
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,NA),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)))),
               "qPCR.K contains a column with all NA.")

  #29. make sure no column is entirely NA in count
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,NA),c(1,1,NA)))),
               "count contains a column with all NA.")

  #30. make sure no data are undefined
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,Inf),c(1,1,NA)))),
               "count contains undefined values \\(i.e., Inf or -Inf\\)")

  #31. make sure no data are undefined
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,Inf),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)))),
               "qPCR.K contains undefined values \\(i.e., Inf or -Inf\\)")

  #32. make sure no data are undefined
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,Inf),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)))),
               "qPCR.N contains undefined values \\(i.e., Inf or -Inf\\)")

  #33. make sure site.cov is not zero-length
  site.cov=matrix(NA,ncol=2,nrow=0)
  colnames(site.cov)=c('var_a','var_b')
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b')),
               "site.cov contains zero-length data.")

  #34. make sure no column is entirely NA in site.cov
  site.cov=rbind(c(4,1,NA),c(1,1,NA))
  colnames(site.cov)=c('var_a','var_b','var_c')
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b')),
               "site.cov contains a column with all NA.")

  #35. make sure no data are undefined
  site.cov=rbind(c(4,1,Inf),c(1,1,NA))
  colnames(site.cov)=c('var_a','var_b','var_c')
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b')),
               "site.cov contains undefined values \\(i.e., Inf or -Inf\\)")

  #36. length of initial values is equal to the number of chains
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu = stats::runif(3, 0.01, 5),
      p10 = stats::runif(1,log(0.0001),log(0.08)),
      alpha = rep(0.1,3)
    )
  }
  site.cov=rbind(c(4,1),c(1,1))
  colnames(site.cov)=c('var_a','var_b')
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b'), initial_values=inits,
                          n.chain=5),
               "The length of the list of initial values should equal the number of chains \\(n.chain, default is 4\\).")

  #37. initial values check: if mu is numeric
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu = stats::runif(3, -1, 0),
      p10 = stats::runif(1,log(0.0001),log(0.08)),
      alpha = rep(0.1,3)
    )
  }
  site.cov=rbind(c(4,1),c(1,1))
  colnames(site.cov)=c('var_a','var_b')
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b'), initial_values=inits),
               "Initial values for 'mu' should be numeric values > 0.")

  #38. initial values check: mu length
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu = stats::runif(4, 0.1, 1),
      p10 = stats::runif(1,log(0.0001),log(0.08)),
      alpha = rep(0.1,3)
    )
  }
  site.cov=rbind(c(4,1),c(1,1))
  colnames(site.cov)=c('var_a','var_b')
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b'), initial_values=inits),
               "The length of initial values for 'mu' should equal the number of sites.")


  #39. initial values check: p10 length
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu = stats::runif(2,0,1),
      p10 = stats::runif(2,log(0.0001),log(0.08)),
      alpha = rep(0.1,3)
    )
  }
  site.cov=rbind(c(4,1),c(1,1))
  colnames(site.cov)=c('var_a','var_b')
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b'), initial_values=inits),
               "The length of initial values for 'p10' should equal 1.")

  #40. initial values check: beta length
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu = stats::runif(2,0,1),
      p10 = stats::runif(1,log(0.0001),log(0.08)),
      beta = c(1,0)
    )
  }
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA))),
                          initial_values=inits),
               "The length of initial values for 'beta' should equal 1.")

  #41. initial values check: alpha length
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu = stats::runif(2,0,1),
      p10 = stats::runif(1,log(0.0001),log(0.08)),
      alpha = rep(0.1,2)
    )
  }
  site.cov=rbind(c(4,1),c(1,1))
  colnames(site.cov)=c('var_a','var_b')
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    site.cov=site.cov),
                          cov=c('var_a','var_b'), initial_values=inits),
               "The length of initial values for 'alpha' should equal\\: \\# covariates \\+ 1 \\(i.e., including intercept\\).")

  #42. initial values check: beta length
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      q = c(0.1,0.1)
    )
  }
  expect_error(jointModel(data=list(qPCR.N=rbind(c(1,1,1),c(1,1,NA)),
                                    qPCR.K=rbind(c(3,3,3),c(3,3,NA)),
                                    count=rbind(c(4,1,1),c(1,1,NA)),
                                    count.type=rbind(c(1,2,1),c(1,1,NA))),
                          initial_values=inits),
               "The length of initial values for 'q' should equal: \\# unique gear types \\- 1 \\(i.e., q for reference type = 1\\).")





})


# correctness and parameter recovery tests
#' @srrstats {G5.4, G5.6} Correctness/parameter recovery tests to test that the implementation produces expected results given data with known properties
test_that("jointModel parameter recovery tests work",{

  ################################
  # model run 1: smaller dataset #
  ################################

  # simulate data: seed 123
  #' @srrstats {G5.5, G5.6} Running correctness test with a fixed random seed
  set.seed(123)
  # constants
  nsite <- 20
  nobs_count <- 200
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    count[i,] <- rpois(nobs_count,mu[i])
  }
  # site-level covariates
  mat_site <- matrix(NA,nrow=nsite,ncol=length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] = mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] = min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    site.cov = mat_site
  )
  # run model
  fit <- jointModel(data=data, cov=c('var_a','var_b'))
  # summary
  summary <- as.data.frame(rstan::summary(fit$model, pars = c('mu','alpha','log_p10'), probs = c(0.025, 0.975))$summary)

  # set up empty vector to check if true values are in 95% interval of posterior estimates
  check <- rep(NA,length(mu)+length(alpha)+length(log_p10))
  # check mu
  for(i in 1:nsite){
    par <- paste0('mu[',i,']')
    if(mu[i] > summary[par,'2.5%'] && mu[i] < summary[par,'97.5%']){
      check[i] <- TRUE
    } else {
      check[i] <- FALSE
    }
  }
  # check alpha
  for(i in 1:length(alpha)){
    par <- paste0('alpha[',i,']')
    if(alpha[i] > summary[par,'2.5%'] && alpha[i] < summary[par,'97.5%']){
      check[nsite+i] <- TRUE
    } else {
      check[nsite+i] <- FALSE
    }
  }
  # check p10
  if(log_p10 > summary['log_p10','2.5%'] && log_p10 < summary['log_p10','97.5%']){
    check[nsite+length(alpha)+1] <- TRUE
  } else {
    check[nsite+length(alpha)+1] <- FALSE
  }

  #' @srrstats {G3.0, G5.6a} Instead of comparing floating point values for equality, here the model is tested to determine if the true parameter values are within the 95% quantiles of the posterior
  # all should be equal to true
  expect_equal(check,rep(TRUE,nsite+length(alpha)+length(log_p10)))

  ##########################################
  # model run 2: smaller dataset, new seed #
  ##########################################

  #' @srrstats {G5.6b,G5.9b} Run test for multiple seeds with same data
  # simulate data: seed 222
  set.seed(222)
  # constants
  nsite <- 20
  nobs_count <- 200
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    count[i,] <- rpois(nobs_count,mu[i])
  }
  # site-level covariates
  mat_site <- matrix(NA,nrow=nsite,ncol=length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] = mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] = min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    site.cov = mat_site
  )
  # run model
  fit <- jointModel(data=data, cov=c('var_a','var_b'))
  # summary
  summary <- as.data.frame(rstan::summary(fit$model, pars = c('mu','alpha','log_p10'), probs = c(0.025, 0.975))$summary)

  # set up empty vector to check if true values are in 95% interval of posterior estimates
  check <- rep(NA,length(mu)+length(alpha)+length(log_p10))
  # check mu
  for(i in 1:nsite){
    par <- paste0('mu[',i,']')
    if(mu[i] > summary[par,'2.5%'] && mu[i] < summary[par,'97.5%']){
      check[i] <- TRUE
    } else {
      check[i] <- FALSE
    }
  }
  # check alpha
  for(i in 1:length(alpha)){
    par <- paste0('alpha[',i,']')
    if(alpha[i] > summary[par,'2.5%'] && alpha[i] < summary[par,'97.5%']){
      check[nsite+i] <- TRUE
    } else {
      check[nsite+i] <- FALSE
    }
  }
  # check p10
  if(log_p10 > summary['log_p10','2.5%'] && log_p10 < summary['log_p10','97.5%']){
    check[nsite+length(alpha)+1] <- TRUE
  } else {
    check[nsite+length(alpha)+1] <- FALSE
  }

  #' @srrstats {G3.0, G5.6a} Instead of comparing floating point values for equality, here the model is tested to determine if the true parameter values are within the 95% quantiles of the posterior
  # all should be equal to true
  expect_equal(check,rep(TRUE,nsite+length(alpha)+length(log_p10)))

  ################################
  # model run 3: larger dataset #
  ################################

  # simulate data: seed 123
  set.seed(123)
  # constants
  nsite <- 40 # increased dataset size
  nobs_count <- 400 # increased dataset size
  nobs_pcr <- 20 # increased dataset size
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    count[i,] <- rpois(nobs_count,mu[i])
  }
  # site-level covariates
  mat_site <- matrix(NA,nrow=nsite,ncol=length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] = mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] = min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    site.cov = mat_site
  )
  # run model
  fit_large <- jointModel(data=data, cov=c('var_a','var_b'))
  # summary
  summary_large <- as.data.frame(rstan::summary(fit_large$model, pars = c('mu','alpha','log_p10'), probs = c(0.025, 0.975))$summary)

  # set up empty vector to check if true values are in 95% interval of posterior estimates
  check <- rep(NA,length(mu)+length(alpha)+length(log_p10))
  # check mu
  for(i in 1:nsite){
    par <- paste0('mu[',i,']')
    if(mu[i] > summary_large[par,'2.5%'] && mu[i] < summary_large[par,'97.5%']){
      check[i] <- TRUE
    } else {
      check[i] <- FALSE
    }
  }
  # check alpha
  for(i in 1:length(alpha)){
    par <- paste0('alpha[',i,']')
    if(alpha[i] > summary[par,'2.5%'] && alpha[i] < summary[par,'97.5%']){
      check[nsite+i] <- TRUE
    } else {
      check[nsite+i] <- FALSE
    }
  }
  # check p10
  if(log_p10 > summary['log_p10','2.5%'] && log_p10 < summary['log_p10','97.5%']){
    check[nsite+length(alpha)+1] <- TRUE
  } else {
    check[nsite+length(alpha)+1] <- FALSE
  }

  #' @srrstats {G3.0, G5.6a} Instead of comparing floating point values for equality, here the model is tested to determine if the true parameter values are within the 95% quantiles of the posterior
  # all should be equal to true
  expect_equal(check,rep(TRUE,nsite+length(alpha)+length(log_p10)))

})
