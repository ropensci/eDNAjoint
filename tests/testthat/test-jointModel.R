test_that("jointModel input checks work", {
  testthat::skip_on_cran()
  #' @srrstats {G5.2,G5.2b,BS2.15} Tests the assure function input checks are
  #'   behaving as expected.
  #1. input tags are valid, q = FALSE, cov = FALSE
  expect_error(jointModel(data = list(qPCR.k = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.n = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          multicore = FALSE),
               paste("Data should include 'qPCR.N', 'qPCR.K', and 'count'.",
                      'See the eDNAjoint guide for data formatting help: ',
                      paste0('https://ednajoint.netlify.app',
                             '/usecase1.html#prepare-the-data'),
                     sep='\n'))
  #2. input tags are valid, q = FALSE, cov = TRUE
  site.cov <- cbind(c(1,0),c(0.4,-0.4))
  colnames(site.cov) <- c('var_a','var_b')
  expect_error(jointModel(data = list(qPCR.k = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.n = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'),
                          multicore = FALSE),
               paste(paste0("Data should include 'qPCR.N', 'qPCR.K', ",
                            "'count', and 'site.cov'."),
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase2.html#prepare-the-data'),
                     sep='\n'))
  #3. input tags are valid, q = TRUE, cov = FALSE
  expect_error(jointModel(data = list(qPCR.k = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.n = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,2,NA))),
                          q = TRUE,
                          multicore = FALSE),
               paste(paste0("Data should include 'qPCR.N', 'qPCR.K', ",
                            "'count', and 'count.type'."),
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase3.html#prepare-the-data'),
                     sep='\n'))
  #4. input tags are valid, q = TRUE, cov = TRUE
  site.cov <- cbind(c(1,0),c(0.4,-0.4))
  colnames(site.cov) <- c('var_a','var_b')
  expect_error(jointModel(data = list(qPCR.k = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.n = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,2,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'),q = TRUE,
                          multicore = FALSE),
               paste(paste0("Data should include 'qPCR.N', 'qPCR.K', 'count', ",
                            "'count.type', and 'site.cov'."),
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase3.html#prepare-the-data'),
                     sep='\n'))

  #5. make sure dimensions of qPCR.N and qPCR.K are equal
  #' @srrstats {BS2.1a} Test to ensure pre-processing routines to ensure all
  #'   input data is dimensionally commensurate
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1,1),c(1,1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          multicore = FALSE),
               paste("Dimensions of qPCR.N and qPCR.K do not match.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n'))

  #6. make sure dimensions of count and count.type are equal, if count.type is
  # present
  #' @srrstats {BS2.1a} Test to ensure pre-processing routines to ensure all
  #'   input data is dimensionally commensurate
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2),c(1,2))),
                          q = TRUE,
                          multicore = FALSE),
               paste("Dimensions of count and count.type do not match.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase3.html#prepare-the-data'),
                     sep='\n'))

  #7. make sure number of rows in count = number of rows in qPCR.N and qPCR.K
  #' @srrstats {BS2.1a} Test to ensure pre-processing routines to ensure all
  #'   input data is dimensionally commensurate
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA),
                                                    c(4,1,1))),
                          multicore = FALSE),
               paste(paste0("Number of sites \\(rows\\) in qPCR data and ",
                            "traditional survey count data do not match."),
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n'))

  #8. make sure all data is numeric -- if q == TRUE
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c('NA',2,2),c(1,2,2))),
                          q = TRUE,
                          multicore = FALSE),
               paste("Data should be numeric.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase3.html#prepare-the-data'),
                     sep='\n'))

  #9. make sure all data is numeric -- if q == FALSE
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,'NA')),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          multicore = FALSE),
               paste("Data should be numeric.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n'))

  #10. make sure locations of NAs in count data match locations of NAs in
  # count.type data
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(NA,2,2),c(1,2,2))),
                          q = TRUE,
                          multicore = FALSE),
               paste(paste0("Empty data cells \\(NA\\) in count data should ",
                            "match ",
                            "empty data cells \\(NA\\) in count.type data."),
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase3.html#prepare-the-data'),
                     sep='\n'))

  #11. make sure locations of NAs in qPCR.N data match locations of NAs in
  # qPCR.K data
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,NA,1)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          multicore = FALSE),
               paste(paste0("Empty data cells \\(NA\\) in qPCR.N data should ",
                            "match ",
                            "empty data cells \\(NA\\) in qPCR.K data."),
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n'))

  #12. make sure family is either 'poisson', 'negbin', or 'gamma'
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          family = 'normal',
                          multicore = FALSE),
               paste0("Invalid family. Options include 'poisson', ",
                      "'negbin', and 'gamma'."))

  #13. p10 priors is a vector of two numeric values
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          p10priors = c(1,1,2),
                          multicore = FALSE),
               paste0("p10priors should be a vector of two positive numeric ",
                      "values. ex. c\\(1,20\\)"))

  #14. p10 priors is a vector of two numeric values
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          multicore = FALSE,
                          p10priors = '1,20'),
               paste0("p10priors should be a vector of two positive numeric ",
                      "values. ex. c\\(1,20\\)"))

  #15. p10 priors is a vector of two numeric values
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          p10priors = c(0,20),
                          multicore = FALSE),
               paste0("p10priors should be a vector of two positive numeric ",
                      "values. ex. c\\(1,20\\)"))

  #16. phi priors is a vector of two numeric values
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          phipriors = c(0,1),family = 'negbin',
                          multicore = FALSE),
               paste0("phipriors should be a vector of two positive numeric ",
                      "values. ex. c\\(0.25,0.25\\)"))

  #17. the smallest count.type is 1
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(0,1,2),c(1,2,NA))),
                          q = TRUE,
                          multicore = FALSE),
               paste(paste0("The first gear type should be referenced as 1 in ",
                            "count.type. Subsequent gear types should be ",
                            "referenced 2, 3, 4, etc."),
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase3.html#prepare-the-data'),
                     sep='\n'))

  #18. count are integers
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4.1,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,1,2),c(1,2,NA))),
                          q = TRUE, family  =  'negbin',
                          multicore = FALSE),
               paste0("All values in count should be non-negative integers. ",
                      "Use family = 'gamma' if count is continuous."))

  #19. qPCR.N are integers
  expect_error(jointModel(data = list(qPCR.K = rbind(c(0.99,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,1,2),c(1,2,NA))),
                          q = TRUE,
                          multicore = FALSE),
               paste("All values in qPCR.K should be non-negative integers.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n'))

  #20. qPCR.K are integers
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3.1,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,1,2),c(1,2,NA))),
                          q = TRUE,
                          multicore = FALSE),
               paste("All values in qPCR.N should be non-negative integers.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n'))

  #21. count.type are integers
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1.1,1,2),c(1,2,NA))),
                          q = TRUE,
                          multicore = FALSE),
               paste("All values in count.type should be integers.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase3.html#prepare-the-data'),
                     sep='\n'))

  #22. site.cov is numeric, if present
  site.cov <- cbind(c('high','low'),c(0.4,-0.4))
  colnames(site.cov) <- c('var_a','var_b')
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,2,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'),q = TRUE,
                          multicore = FALSE),
               paste("site.cov should be numeric.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app','/usecase2.html'),
                     sep='\n'))

  #23. cov values match column names in site.cov
  site.cov <- cbind(c(0,1),c(0.4,-0.4))
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,2,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'),q = TRUE,
                          multicore = FALSE),
               paste(paste0("cov values should be listed in the column names of ",
                            "site.cov in the data."),
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app','/usecase2.html'),
                     sep='\n'))

  #24. site.cov has same number of rows as qPCR.N and count, if present
  site.cov <- cbind(c(0,1,1),c(0.4,-0.4,1))
  colnames(site.cov) <- c('var_a','var_b')
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,2,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'),q = TRUE,
                          multicore = FALSE),
               paste(paste0("The number of rows in site.cov matrix should match ",
                            "the number of rows in all other matrices."),
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app','/usecase2.html'),
                     sep='\n'))

  #25. make sure count.type is not zero-length
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = matrix(NA,ncol = 3,
                                                          nrow = 0)),
                          q = TRUE,
                          multicore = FALSE),
               paste("count.type contains zero-length data.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase3.html#prepare-the-data'),
                     sep='\n'))

  #26. make sure no column is entirely NA in count.type
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(4,1,NA),c(1,1,NA))),
                          q = TRUE,
                          multicore = FALSE),
               paste("count.type contains a column with all NA.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase3.html#prepare-the-data'),
                     sep='\n'))

  #27. make sure no column is entirely NA in qPCR.K
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,NA),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          multicore = FALSE),
               paste("qPCR.K contains a column with all NA.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n'))

  #28. make sure no column is entirely NA in qPCR.N
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,NA),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          multicore = FALSE),
               paste("qPCR.N contains a column with all NA.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n'))

  #29. make sure no column is entirely NA in count
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,NA),c(1,1,NA))),
                          multicore = FALSE),
               paste("count contains a column with all NA.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n'))

  #30. make sure no data are undefined
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,Inf),c(1,1,NA))),
                          multicore = FALSE),
               paste("count contains undefined values \\(i.e., Inf or -Inf\\)",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n'))

  #31. make sure no data are undefined
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,Inf),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          multicore = FALSE),
               paste("qPCR.N contains undefined values \\(i.e., Inf or -Inf\\)",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n'))

  #32. make sure no data are undefined
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,Inf),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          multicore = FALSE),
               paste("qPCR.K contains undefined values \\(i.e., Inf or -Inf\\)",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n'))

  #33. make sure site.cov is not zero-length
  site.cov <- matrix(NA,ncol = 2,nrow = 0)
  colnames(site.cov) <- c('var_a','var_b')
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'),
                          multicore = FALSE),
               paste("site.cov contains zero-length data.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app','/usecase2.html'),
                     sep='\n'))

  #34. make sure no column is entirely NA in site.cov
  site.cov <- rbind(c(4,1,NA),c(1,1,NA))
  colnames(site.cov) <- c('var_a','var_b','var_c')
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'),
                          multicore = FALSE),
               paste("site.cov contains a column with all NA.",
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app','/usecase2.html'),
                     sep='\n'))

  #35. make sure no data are undefined
  site.cov <- rbind(c(4,1,Inf),c(1,1,NA))
  colnames(site.cov) <- c('var_a','var_b','var_c')
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'),
                          multicore = FALSE),
               paste(paste0("site.cov contains undefined values \\(i.e., ",
                            "Inf or -Inf\\)"),
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app','/usecase2.html'),
                     sep='\n'))

  #36. length of initial values is equal to the number of chains
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu <- stats::runif(3, 0.01, 5),
      p10 <- exp(stats::runif(1,log(0.0001),log(0.08))),
      alpha <- rep(0.1,3)
    )
  }
  site.cov <- rbind(c(4,1),c(1,1))
  colnames(site.cov) <- c('var_a','var_b')
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'), initial_values = inits,
                          n.chain = 5,
                          multicore = FALSE),
               paste(paste0("The length of the list of initial values ",
                            "should equal the number of chains \\(n.chain, ",
                            "default is 4\\)."),
                     paste0('See the eDNAjoint guide for help formatting ',
                            'initial values: '),
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#initialvalues'),
                     sep='\n'))

  #37. initial values check: if mu is numeric
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu <- stats::runif(3, -1, 0),
      p10 <- exp(stats::runif(1,log(0.0001),log(0.08))),
      alpha <- rep(0.1,3)
    )
    names(inits[[i]]) <- c('mu','p10','alpha')
  }
  site.cov <- rbind(c(4,1),c(1,1))
  colnames(site.cov) <- c('var_a','var_b')
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'), initial_values = inits,
                          multicore = FALSE),
               paste("Initial values for 'mu' should be numeric values > 0.",
                     paste0('See the eDNAjoint guide for help formatting ',
                            'initial values: '),
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#initialvalues'),
                     sep='\n'))

  #38. initial values check: mu length
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu <- stats::runif(4, 0.1, 1),
      p10 <- exp(stats::runif(1,log(0.0001),log(0.08))),
      alpha <- rep(0.1,3)
    )
    names(inits[[i]]) <- c('mu','p10','alpha')
  }
  site.cov <- rbind(c(4,1),c(1,1))
  colnames(site.cov) <- c('var_a','var_b')
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'), initial_values = inits,
                          multicore = FALSE),
               paste(paste0("The length of initial values for 'mu' should ",
                            "equal the number of sites."),
                      paste0('See the eDNAjoint guide for help formatting ',
                             'initial values: '),
                      paste0('https://ednajoint.netlify.app',
                             '/usecase1.html#initialvalues'),
                      sep='\n'))


  #39. initial values check: p10 numeric
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu <- stats::runif(2,0,1),
      p10 <- "-1",
      alpha <- rep(0.1,3)
    )
    names(inits[[i]]) <- c('mu','p10','alpha')
  }
  site.cov <- rbind(c(4,1),c(1,1))
  colnames(site.cov) <- c('var_a','var_b')
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'), initial_values = inits,
                          multicore = FALSE),
               paste("Initial values for 'p10' should be numeric.",
                     paste0('See the eDNAjoint guide for help formatting ',
                            'initial values: '),
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#initialvalues'),
                     sep='\n'))

  #40. initial values check: p10 length
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu <- stats::runif(2,0,1),
      p10 <- exp(stats::runif(2,log(0.0001),log(0.08))),
      alpha <- rep(0.1,3)
    )
    names(inits[[i]]) <- c('mu','p10','alpha')
  }
  site.cov <- rbind(c(4,1),c(1,1))
  colnames(site.cov) <- c('var_a','var_b')
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'), initial_values = inits,
                          multicore = FALSE),
               paste("The length of initial values for 'p10' should equal 1.",
                     paste0('See the eDNAjoint guide for help formatting ',
                            'initial values: '),
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#initialvalues'),
                     sep='\n'))

  #41. initial values check: alpha numeric
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu <- stats::runif(2,0,1),
      p10 <- exp(stats::runif(1,log(0.0001),log(0.08))),
      alpha <- "1"
    )
    names(inits[[i]]) <- c('mu','p10','alpha')
  }
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          initial_values = inits,
                          multicore = FALSE),
               paste("Initial values for 'alpha' should be numeric.",
                     paste0('See the eDNAjoint guide for help formatting ',
                            'initial values: '),
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#initialvalues'),
                     sep='\n'))

  #42. initial values check: alpha length
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu <- stats::runif(2,0,1),
      p10 <- exp(stats::runif(1,log(0.0001),log(0.08))),
      alpha <- c(1,0)
    )
    names(inits[[i]]) <- c('mu','p10','alpha')
  }
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA))),
                          initial_values = inits,
                          multicore = FALSE),
               paste("The length of initial values for 'alpha' should equal 1.",
                     paste0('See the eDNAjoint guide for help formatting ',
                            'initial values: '),
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#initialvalues'),
                     sep='\n'))

  #43. initial values check: alpha numeric
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu <- stats::runif(2,0,1),
      p10 <- exp(stats::runif(1,log(0.0001),log(0.08))),
      alpha <- c("1","2")
    )
    names(inits[[i]]) <- c('mu','p10','alpha')
  }
  site.cov <- rbind(c(4,1),c(1,1))
  colnames(site.cov) <- c('var_a','var_b')
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'), initial_values = inits,
                          multicore = FALSE),
               paste("Initial values for 'alpha' should be numeric.",
                     paste0('See the eDNAjoint guide for help formatting ',
                            'initial values: '),
                     paste0('https://ednajoint.netlify.app',
                            '/usecase2.html#initialvalues'),
                     sep='\n'))

  #44. initial values check: alpha length
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      mu <- stats::runif(2,0,1),
      p10 <- exp(stats::runif(1,log(0.0001),log(0.08))),
      alpha <- rep(0.1,2)
    )
    names(inits[[i]]) <- c('mu','p10','alpha')
  }
  site.cov <- rbind(c(4,1),c(1,1))
  colnames(site.cov) <- c('var_a','var_b')
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      site.cov = site.cov),
                          cov = c('var_a','var_b'), initial_values = inits,
                          multicore = FALSE),
               paste(paste0("The length of initial values for 'alpha' should ",
                            "equal\\: \\# covariates \\+ 1 \\(i.e., ",
                            "including intercept\\)."),
                     paste0('See the eDNAjoint guide for help formatting ',
                            'initial values: '),
                     paste0('https://ednajoint.netlify.app',
                            '/usecase2.html#initialvalues'),
                     sep='\n'))

  #45. initial values check: q numeric
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      q <- "0.1"
    )
    names(inits[[i]]) <- c('q')
  }
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          initial_values = inits,
                          multicore = FALSE),
               paste("Initial values for 'q' should be numeric.",
                     paste0('See the eDNAjoint guide for help formatting ',
                            'initial values: '),
                     paste0('https://ednajoint.netlify.app',
                            '/usecase2.html#initialvalues'),
                     sep='\n'))

  #46. initial values check: q length
  n.chain <- 4
  inits <- list()
  for(i in 1:n.chain){
    inits[[i]] <- list(
      q <- c(0.1,0.1)
    )
    names(inits[[i]]) <- c('q')
  }
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          initial_values = inits,
                          multicore = FALSE),
               paste(paste0("The length of initial values for 'q' should ",
                            "equal: \\# unique gear types \\- 1 \\(i.e., q ",
                            "for reference type = 1\\)."),
                     paste0('See the eDNAjoint guide for help formatting ',
                            'initial values: '),
                     paste0('https://ednajoint.netlify.app',
                            '/usecase2.html#initialvalues'),
                     sep='\n'))

  #47. check length and range of n.chain
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          n.chain  =  c(1,1),multicore = FALSE),
               paste0("n.chain should be an integer > 0 and of length 1."))
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          n.chain = 0,multicore = FALSE),
               paste0("n.chain should be an integer > 0 and of length 1."))

  #48. check length and range of n.iter.sample
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          n.iter.sample = c(1,1),multicore = FALSE),
               paste0("n.iter.sample should be an integer > 0 and of length 1.")
               )
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          n.iter.sample = 0,multicore = FALSE),
               paste0("n.iter.sample should be an integer > 0 and of length 1.")
               )

  #49. check length and range of n.iter.burn
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          n.iter.burn = c(1,1),multicore = FALSE),
               paste0("n.iter.burn should be an integer > 0 and of length 1.")
  )
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          n.iter.burn = 0,multicore = FALSE),
               paste0("n.iter.burn should be an integer > 0 and of length 1.")
  )

  #50. check length and range of thin
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          thin = c(1,1),multicore = FALSE),
               paste0("thin should be an integer > 0 and of length 1.")
  )
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          thin = 0,multicore = FALSE),
               paste0("thin should be an integer > 0 and of length 1.")
  )

  #51. check length and range of adapt_delta
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          adapt_delta = c(0.9,0.9),multicore = FALSE),
               paste0("adapt_delta should be a numeric value > 0 and < 1 and ",
                      "of length 1.")
  )
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          adapt_delta = 1.2,multicore = FALSE),
               paste0("adapt_delta should be a numeric value > 0 and < 1 and ",
                      "of length 1.")
  )

  #52. check length of seed
  expect_error(jointModel(data = list(qPCR.K = rbind(c(1,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                                      count.type = rbind(c(1,2,1),c(1,1,NA))),
                          seed = c(1,2),multicore = FALSE),
               paste0("seed should be an integer of length 1.")
  )

  #53. check K <= N
  expect_error(jointModel(data = list(qPCR.K = rbind(c(4,1,1),c(1,1,NA)),
                                      qPCR.N = rbind(c(3,3,3),c(3,3,NA)),
                                      count = rbind(c(4,1,1),c(1,1,NA)),
                          multicore = FALSE)),
               paste(paste0("N should be >= K in qPCR data. N is the number ",
                            "of qPCR replicates per sample, and K is the ",
                            "number of positive detections among replicates."),
                     'See the eDNAjoint guide for data formatting help: ',
                     paste0('https://ednajoint.netlify.app',
                            '/usecase1.html#prepare-the-data'),
                     sep='\n')
  )


})


# correctness and parameter recovery tests
#' @srrstats {G5.4, G5.6} Correctness/parameter recovery tests to test that
#'   the implementation produces expected results given data with known
#'   properties
#' @srrstats {PD4.0} These tests do not test for numeric equality of outputs,
#'   but rather test for the recovery of parameter values. In these tests, I
#'   simulate data with known parameter values, run the main statistical
#'   function (jointModel()) with the data, and then test if the known parameter
#'   values are within the 95% credibility interval of the parameters'
#'   posteriors in the function output.
test_that("jointModel parameter recovery tests work",{
  testthat::skip_on_cran()

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
  mu <- rlnorm(nsite,meanlog = log(1),sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  # count
  count <- matrix(NA,nrow = nsite,ncol = nobs_count)
  for(i in 1:nsite){
    count[i,] <- rpois(nobs_count,mu[i])
  }
  # site-level covariates
  mat_site <- matrix(NA,nrow = nsite,ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability of
  # eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA,nrow = nsite,ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow = nsite,ncol = nobs_pcr)
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
  fit1 <- suppressWarnings({jointModel(data = data, cov = c('var_a','var_b'),
                    multicore = FALSE, seed = 10)})
  # summary
  summary1 <- as.data.frame(rstan::summary(fit1$model,
                                           pars = c('mu','alpha','log_p10'),
                                           probs = c(0.025, 0.975))$summary)

  # set up empty vector to check if true values are in 95% interval of
  # posterior estimates
  check <- rep(NA,length(alpha)+length(log_p10))
  # check alpha
  for(i in seq_along(alpha)){
    par <- paste0('alpha[',i,']')
    if(alpha[i] > summary1[par,'2.5%'] && alpha[i] < summary1[par,'97.5%']){
      check[i] <- TRUE
    } else {
      check[i] <- FALSE
    }
  }
  # check p10
  if(log_p10 > summary1['log_p10','2.5%'] &&
     log_p10 < summary1['log_p10','97.5%']){
    check[length(alpha)+1] <- TRUE
  } else {
    check[length(alpha)+1] <- FALSE
  }

  #' @srrstats {G3.0, G5.6a} Instead of comparing floating point values for
  #'   equality, here the model is tested to determine if the true parameter
  #'   values are within the 95% quantiles of the posterior
  # all should be equal to true
  expect_equal(check,rep(TRUE,length(alpha)+length(log_p10)))

  # test that output values are on the same scale as the data
  mu_estimates <- rep(NA,nsite)
  for(i in 1:nsite){
    mu_estimates[i] <- summary1[paste0('mu[',i,',1]'),'mean']
  }

  # get mean of input count data at each site
  data_mean <- rep(NA,nsite)
  for(i in 1:nsite){
    data_mean[i] <- mean(data$count[i,])
  }
  #' @srrstats {BS7.4, BS7.4a} Check to ensure that mean posterior estimates
  #'   are on the same scale as the mean of the input data (here checking
  #'   estimates of mu, i.e., expected catch rate at each site)
  # check that estimates are on same scale as data
  expect_equal(round(mu_estimates,0),round(data_mean,0))

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
  mu <- rlnorm(nsite,meanlog = log(1),sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  # count
  count <- matrix(NA,nrow = nsite,ncol = nobs_count)
  for(i in 1:nsite){
    count[i,] <- rpois(nobs_count,mu[i])
  }
  # site-level covariates
  mat_site <- matrix(NA,nrow = nsite,ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability of
  # eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA,nrow = nsite,ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow = nsite,ncol = nobs_pcr)
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
  fit2 <- suppressWarnings({jointModel(data = data, cov = c('var_a','var_b'),
                    multicore = FALSE, seed = 2)})
  # summary
  summary2 <- as.data.frame(rstan::summary(fit2$model,
                                           pars = c('mu','alpha','log_p10'),
                                           probs = c(0.025, 0.975))$summary)

  # set up empty vector to check if true values are in 95% interval of
  # posterior estimates
  check <- rep(NA,length(alpha)+length(log_p10))
  # check alpha
  for(i in seq_along(alpha)){
    par <- paste0('alpha[',i,']')
    if(alpha[i] > summary2[par,'2.5%'] && alpha[i] < summary2[par,'97.5%']){
      check[i] <- TRUE
    } else {
      check[i] <- FALSE
    }
  }
  # check p10
  if(log_p10 > summary2['log_p10','2.5%'] &&
     log_p10 < summary2['log_p10','97.5%']){
    check[length(alpha)+1] <- TRUE
  } else {
    check[length(alpha)+1] <- FALSE
  }

  #' @srrstats {G3.0, G5.6a} Instead of comparing floating point values for
  #'   equality, here the model is tested to determine if the true parameter
  #'   values are within the 95% quantiles of the posterior
  # all should be equal to true
  expect_equal(check,rep(TRUE,length(alpha)+length(log_p10)))

  ################################
  # model run 3: larger dataset #
  ################################

  # simulate data: seed 123
  set.seed(123)
  # constants
  nsite <- 20
  nobs_count <- 200 # increased dataset size (doubled)
  nobs_pcr <- 16 # increased dataset size (doubled)
  # params
  mu <- rlnorm(nsite,meanlog = log(1),sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  # count
  count <- matrix(NA,nrow = nsite,ncol = nobs_count)
  for(i in 1:nsite){
    count[i,] <- rpois(nobs_count,mu[i])
  }
  # site-level covariates
  mat_site <- matrix(NA,nrow = nsite,ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA,nrow = nsite,ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow = nsite,ncol = nobs_pcr)
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
  fit_large <- suppressWarnings({jointModel(data = data,
                                            cov = c('var_a','var_b'),
                                            multicore = FALSE, seed = 10)})
  # summary
  summary_large <- as.data.frame(rstan::summary(fit_large$model,
                                                pars = 'alpha',
                                                probs = c(0.025,
                                                          0.975))$summary)

  # compare standard errors of estimates of alpha coefficients with
  # implementation with a smaller dataset
  se_large <- c(summary_large['alpha[2]','se_mean'],
                summary_large['alpha[3]','se_mean']
  )

  # get values for smaller dataset
  se_small <- c(summary1['alpha[2]','se_mean'],
                summary1['alpha[3]','se_mean']
  )

  # set up empty vector to check if standard errors are smaller with
  # larger dataset
  check <- rep(NA,length(alpha)-1)
  for(i in seq_along(check)){
    check[i] <- se_large[i] < se_small[i]
  }

  #' @srrstats {G5.7} Check to see that implementation performs as expected
  #'   properties of data change (i.e., standard error of posteriors is smaller
  #'   if there are more data observations)
  # all should be equal to true
  expect_equal(check,rep(TRUE,length(alpha)-1))

})


test_that("jointModel probability distribution tests work",{
  testthat::skip_on_cran()

  #' @srrstats {PD4.2} This package fits models with fixed distributions to
  #'   data. Users can choose the distribution used to represent the
  #'   data generating process for traditional survey data through an input
  #'   argument to the function (family). This test assesses whether estimates
  #'   of mu, or the site-level expected catch rate, is different when a
  #'   negative binomial or poisson distribution is used to represent the data
  #'   generating process.

  #######################################################
  # generate data with a negative binomial distribution #
  set.seed(222)
  # constants
  nsite <- 20
  nobs_count <- 200
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1),sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5
  phi <- 0.1
  # count
  count <- matrix(NA,nrow = nsite,ncol = nobs_count)
  for(i in 1:nsite){
    count[i,] <- rnbinom(nobs_count,mu = mu[i],size = phi)
  }
  # p11 (probability of true positive eDNA detection) and p (probability of
  # eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA,nrow = nsite,ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow = nsite,ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count
  )

  ##############################################################
  # fit data with a negative binomial and poisson distribution #

  # run model -- negative binomial distribution
  negbin.fit <- suppressWarnings({jointModel(data = data, multicore = FALSE,
                                             family = 'negbin',seed = 2)})
  # run model -- poisson distribution
  pois.fit <- suppressWarnings({jointModel(data = data, multicore = FALSE,
                                           family = 'poisson',seed = 2)})

  # summarize outputs
  negbin_summary <- as.data.frame(rstan::summary(negbin.fit$model,
                                                 pars = 'mu',
                                                 probs = c(0.025,
                                                         0.975))$summary)
  pois_summary <- as.data.frame(rstan::summary(pois.fit$model,
                                                 pars = 'mu',
                                                 probs = c(0.025,
                                                         0.975))$summary)

  # summarize differences in mean estimates of mu
  summary_table <- table(negbin_summary$mean > pois_summary$mean)

  # calculate percentage of mean estimates of mu that are greater when
  # a negative binomial distribution are used rather than a poisson distribution
  percent <- summary_table[[2]]/(summary_table[[2]]+summary_table[[1]])


  # test that estimates of mu when a negative binomial distribution is used is
  # greater than estimates of mu when a poisson distribution is used for
  # >80% of sites
  expect_gt(percent,0.8)


})

test_that("semi-paired model works", {
  testthat::skip_on_cran()

  ## 1.
  # model includes 'p10','q','phi','beta'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5
  q <- 2
  phi <- 1.2
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rnbinom(n = 1, mu = mu[i], size = phi)
      } else {
        count[i,j] <- rnbinom(n = 1, mu = mu[i]*q, size = phi)
      }
    }
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2,4),] <- NA
  count_type[c(2,4),] <- NA

  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    count.type = count_type
  )

  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = beta,
    phi = phi,
    q = q
  )
  names(inits[[1]]) <- c('mu','p10','alpha','phi','q')

  # run model
  fit <- suppressWarnings({jointModel(data = data, q = TRUE, family = 'negbin',
                                      initial_values = inits, n.iter.burn = 25,
                                      n.iter.sample = 75,
                                      n.chain = 1, multicore = FALSE, seed = 10
  )})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','q[1]','phi','alpha[1]') %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model,par='mu')$summary[,1])
  expect_equal(length(mu_estimates),nsite*2)
  expect_true(is.numeric(mu_estimates))


  ## 2.
  # model includes 'p10','beta','q'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5
  q <- 2
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rpois(1, mu[i])
      } else {
        count[i,j] <- rpois(1, mu[i]*q)
      }
    }
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i], nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2,4),] <- NA
  count_type[c(2,4),] <- NA

  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    count.type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = beta,
    q = 2
  )
  names(inits[[1]]) <- c('mu','p10','alpha','q')
  # run model
  fit <- suppressWarnings({jointModel(data = data, q = TRUE,
                                      n.chain = 1, multicore = FALSE, seed = 10,
                                      initial_values = inits, n.iter.burn = 25,
                                      n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','q[1]','alpha[1]') %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model,par='mu')$summary[,1])
  expect_equal(length(mu_estimates),nsite*2)
  expect_true(is.numeric(mu_estimates))

  ## 3.
  # model includes 'p10','beta','q', 'alpha_gamma','beta_gamma'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5
  q <- 2
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rgamma(1,shape = alpha_gamma[i],rate = beta_gamma)
      } else {
        count[i,j] <- rgamma(1,shape = alpha_gamma[i]*q,rate = beta_gamma)
      }
    }
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2,4),] <- NA
  count_type[c(2,4),] <- NA

  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    count.type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1,length(mu)),
    p10 = exp(log_p10),
    alpha = beta,
    q = q
  )
  names(inits[[1]]) <- c('alpha_gamma','beta_gamma','p10','alpha','q')
  # run model
  fit <- suppressWarnings({jointModel(data = data, q = TRUE, family = 'gamma',
                                      n.chain = 1, multicore = FALSE, seed = 10,
                                      initial_values = inits, n.iter.burn = 25,
                                      n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','q[1]','alpha[1]') %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model,par='mu')$summary[,1])
  expect_equal(length(mu_estimates),nsite*2)
  expect_true(is.numeric(mu_estimates))

  ## 4.
  # model includes 'p10','phi','beta'

  # constants
  nsite <- 50
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5
  phi <- 1.2
  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    count[i,] <- rnbinom(n = nobs_count, mu = mu[i], size = phi)
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count
  )

  # add NA to sites 2 and 4
  count[c(2,4),] <- NA

  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = beta,
    phi = phi
  )
  names(inits[[1]]) <- c('mu','p10','alpha','phi')
  # run model
  fit <- suppressWarnings({jointModel(data = data, family = 'negbin',
                                      initial_values = inits,
                                      n.chain = 1, multicore = FALSE, seed = 10,
                                      n.iter.burn = 25,
                                      n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','phi','alpha[1]') %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model,par='mu')$summary[,1])
  expect_equal(length(mu_estimates),nsite)
  expect_true(is.numeric(mu_estimates))


  ## 5.
  # model includes 'p10','beta'

  # constants
  nsite <- 50
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    count[i,] <- rpois(nobs_count, mu[i])
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2,4),] <- NA

  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = beta
  )
  names(inits[[1]]) <- c('mu','p10','alpha')
  # run model
  fit <- suppressWarnings({jointModel(data = data, initial_values = inits,
                                      n.chain = 1, multicore = FALSE, seed = 10,
                                      n.iter.burn = 25,
                                      n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','alpha[1]') %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model,par='mu')$summary[,1])
  expect_equal(length(mu_estimates),nsite)
  expect_true(is.numeric(mu_estimates))


  ## 6.
  # model includes 'p10','beta','alpha_gamma','alpha_beta'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    count[i,] <- rgamma(nobs_count,shape = alpha_gamma[i],rate = beta_gamma)
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2,4),] <- NA

  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1,length(mu)),
    p10 = exp(log_p10),
    alpha = beta
  )
  names(inits[[1]]) <- c('alpha_gamma','beta_gamma','p10','alpha')
  # run model
  fit <- suppressWarnings({jointModel(data = data, initial_values = inits,
                                      family = 'gamma',
                                      n.chain = 1, multicore = FALSE, seed = 10,
                                      n.iter.burn = 25,
                                      n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','alpha[1]') %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model,par='mu')$summary[,1])
  expect_equal(length(mu_estimates),nsite)
  expect_true(is.numeric(mu_estimates))


  ## 7.
  # model includes 'p10','q','phi','alpha'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  q <- 2
  phi <- 10
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rnbinom(n = 1, mu = mu[i], size = phi)
      } else {
        count[i,j] <- rnbinom(n = 1, mu = mu[i]*q, size = phi)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2,4),] <- NA
  count_type[c(2,4),] <- NA

  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    count.type = count_type,
    site.cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha,
    q = q,
    phi = phi
  )
  names(inits[[1]]) <- c('mu','p10','alpha','q','phi'
  )
  # run model
  fit <- suppressWarnings({jointModel(data = data, family = 'negbin', q = TRUE,
                                      cov = c('var_a','var_b'),
                                      n.chain = 1, multicore = FALSE, seed = 10,
                                      initial_values = inits, adapt_delta = 0.99,
                                      n.iter.burn = 25,
                                      n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','q[1]','phi','alpha[1]',
                    'alpha[2]','alpha[3]') %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model,par='mu')$summary[,1])
  expect_equal(length(mu_estimates),nsite*2)
  expect_true(is.numeric(mu_estimates))


  ## 8.
  # model includes 'p10','alpha','q'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  q <- 2
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rpois(n = 1, mu[i])
      } else {
        count[i,j] <- rpois(n = 1, mu[i]*q)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # add NA to sites 2 and 4
  count[c(2,4),] <- NA
  count_type[c(2,4),] <- NA

  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    count.type = count_type,
    site.cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha
  )
  names(inits[[1]]) <- c('mu','p10','alpha'
  )
  # run model
  fit <- suppressWarnings({jointModel(data = data, q = TRUE,
                                      cov = c('var_a','var_b'),
                                      n.chain = 1, multicore = FALSE, seed = 10,
                                      initial_values = inits, n.iter.burn = 25,
                                      n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','q[1]','alpha[1]',
                    'alpha[2]','alpha[3]') %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model,par='mu')$summary[,1])
  expect_equal(length(mu_estimates),nsite*2)
  expect_true(is.numeric(mu_estimates))


  ## 9.
  # model includes 'p10','alpha','q', 'alpha_gamma', 'alpha_beta'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  q <- 2
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rgamma(1,shape = alpha_gamma[i],rate = beta_gamma)
      } else {
        count[i,j] <- rgamma(1,shape = alpha_gamma[i]*q,rate = beta_gamma)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2,4),] <- NA
  count_type[c(2,4),] <- NA

  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    count.type = count_type,
    site.cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1,length(mu)),
    p10 = exp(log_p10),
    alpha = alpha,
    q = q
  )
  names(inits[[1]]) <- c('alpha_gamma','beta_gamma','p10','alpha','q'
  )
  # run model
  fit <- suppressWarnings({jointModel(data = data, q = TRUE, family = 'gamma',
                                      cov = c('var_a','var_b'),
                                      n.chain = 1, multicore = FALSE, seed = 10,
                                      initial_values = inits, n.iter.burn = 25,
                                      n.iter.sample = 75
  )})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','q[1]','alpha[1]') %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model,par='mu')$summary[,1])
  expect_equal(length(mu_estimates),nsite*2)
  expect_true(is.numeric(mu_estimates))


  ## 10.
  # model includes 'p10','phi','alpha'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  phi <- 10

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    count[i,] <- rnbinom(n = nobs_count, mu = mu[i], size = phi)
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2,4),] <- NA

  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    site.cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha,
    phi = phi
  )
  names(inits[[1]]) <- c('mu','p10','alpha',
                         'phi')
  # run model
  fit <- suppressWarnings({jointModel(data = data, family = 'negbin',
                                      cov = c('var_a','var_b'),n.iter.burn = 25,
                                      n.iter.sample = 75,
                                      n.chain = 1, multicore = FALSE, seed = 10,
                                      initial_values = inits)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','phi','alpha[1]') %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model,par='mu')$summary[,1])
  expect_equal(length(mu_estimates),nsite)
  expect_true(is.numeric(mu_estimates))


  ## 11.
  # model includes 'p10','alpha'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    count[i,] <- rpois(nobs_count, mu[i])
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2,4),] <- NA

  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    site.cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha
  )
  names(inits[[1]]) <- c('mu','p10','alpha'
  )
  # run model
  fit <- suppressWarnings({jointModel(data = data,
                                      cov = c('var_a','var_b'),
                                      n.chain = 1, multicore = FALSE, seed = 10,
                                      initial_values = inits, n.iter.burn = 25,
                                      n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','alpha[1]',
                    'alpha[2]','alpha[3]') %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model,par='mu')$summary[,1])
  expect_equal(length(mu_estimates),nsite)
  expect_true(is.numeric(mu_estimates))


  ## 12.
  # model includes 'p10','alpha','alpha_gamma','beta_gamma'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    count[i,] <- rgamma(nobs_count,shape = alpha_gamma[i],rate = beta_gamma)
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2,4),] <- NA

  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    site.cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1,length(mu)),
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha
  )
  names(inits[[1]]) <- c('alpha_gamma','beta_gamma','mu','p10','alpha'
  )
  # run model
  fit <- suppressWarnings({jointModel(data = data,family = 'gamma',
                                      cov = c('var_a','var_b'),
                                      n.chain = 1, multicore = FALSE, seed = 10,
                                      initial_values = inits, n.iter.burn = 25,
                                      n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','alpha[1]',
                    'alpha[2]','alpha[3]') %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model,par='mu')$summary[,1])
  expect_equal(length(mu_estimates),nsite)
  expect_true(is.numeric(mu_estimates))


})

