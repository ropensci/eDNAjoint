test_that("traditional_model input checks work", {
  testthat::skip_on_cran()
  #' @srrstats {G5.2,G5.2b,BS2.15} Tests the assure function input checks are
  #'   behaving as expected.
  #1. input tags are valid, q = TRUE
  expect_error(traditional_model(data = list(Count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2, 1),
                                                                c(1, 2, NA))),
                                 q = TRUE,
                                 multicore = FALSE),
               "Data should include 'count' and 'count_type'.")

  #2. input tags are valid, q = FALSE
  expect_error(traditional_model(data = list(Count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA))),
                                 multicore = FALSE),
               "Data should include 'count'.")

  #3. make sure dimensions of count and count_type are equal, if count_type is
  # present
  #' @srrstats {BS2.1a} Test to ensure pre-processing routines to ensure all
  #'   input data is dimensionally commensurate
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2),
                                                                c(1, 2))),
                                 q = TRUE,
                                 multicore = FALSE),
               "Dimensions of count and count_type do not match.")

  #4. make sure all data is numeric -- if q == TRUE
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c("NA", 2, 2),
                                                                c(1, 2, 2))),
                                 q = TRUE,
                                 multicore = FALSE),
               "Data should be numeric.")

  #5. make sure all data is numeric -- if q == FALSE
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, "NA"))),
                                 multicore = FALSE),
               "Data should be numeric.")

  #6. make sure locations of NAs in count data match locations of NAs in
  # count_type data
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(NA, 2, 2),
                                                                c(1, 2, 2))),
                                 q = TRUE,
                                 multicore = FALSE),
               paste0("Empty data cells \\(NA\\) in count data should match ",
                      "empty data cells \\(NA\\) in count_type data."))

  #7. make sure family is either 'poisson', 'negbin', or 'gamma'
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA))),
                                 family = "normal",
                                 multicore = FALSE),
               paste0("Invalid family. Options include 'poisson', ",
                      "'negbin', or 'gamma'."))

  #8. the smallest count_type is 1
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(0, 1, 2),
                                                                c(1, 2, NA))),
                                 q = TRUE,
                                 multicore = FALSE),
               paste0("The first gear type should be referenced as 1 in ",
                      "count_type. Subsequent gear types should be ",
                      "referenced 2, 3, 4, etc."))

  #9. count are integers
  expect_error(traditional_model(data = list(count = rbind(c(4.1, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 1, 2),
                                                                c(1, 2, NA))),
                                 q = TRUE, family = "negbin",
                                 multicore  = FALSE),
               paste0("All values in count should be non-negative integers. ",
                      "Use family = 'gamma' if count is continuous."))


  #10. count_type are integers
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1.1, 1, 2),
                                                                c(1, 2, NA))),
                                 q = TRUE,
                                 multicore = FALSE),
               "All values in count_type should be integers.")

  #11. phi priors is a vector of two numeric values
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA))),
                                 phi_priors = c(0, 1), family = "negbin",
                                 multicore = FALSE),
               paste0("phi_priors should be a vector of two positive ",
                      "numeric values. ex. c\\(0.25,0.25\\)"))

  #12. make sure no column is entirely NA in count
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, NA),
                                                           c(1, 1, NA))),
                                 multicore = FALSE),
               "count contains a column with all NA.")

  #13. make sure no data are undefined
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, Inf),
                                                           c(1, 1, NA))),
                                 multicore = FALSE),
               "count contains undefined values \\(i.e., Inf or -Inf\\)")

  #14. length of initial values is equal to the number of chains
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      mu <- stats::runif(2, 0, 1)
    )
    names(inits[[i]]) <- c("mu")
  }
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA))),
                                 n_chain = 5, initial_values = inits,
                                 multicore = FALSE),
               paste0("The length of the list of initial values should equal ",
                      "the number of chains \\(n_chain, default is 4\\)."))

  #15. initial values check mu length
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      mu <- stats::runif(3, 0, 1)
    )
    names(inits[[i]]) <- c("mu")
  }
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA))),
                                 initial_values = inits,
                                 multicore = FALSE),
               paste0("The length of initial values for 'mu' should ",
                      "equal the number of sites."))

  #16. initial values check mu is positive numeric
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      mu <- stats::runif(3, -1, 0)
    )
    names(inits[[i]]) <- c("mu")
  }
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA))),
                                 initial_values = inits,
                                 multicore = FALSE),
               "Initial values for 'mu' should be numeric values > 0.")

  #17. initial values check q
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      q <- c(0.1, 0.1)
    )
    names(inits[[i]]) <- c("q")
  }
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 1, 2),
                                                                c(1, 2, NA))),
                                 initial_values = inits,
                                 multicore = FALSE),
               paste0("The length of initial values for 'q' should equal: ",
                      "\\# unique gear types \\- 1 \\(i.e., q for reference ",
                      "type = 1\\)."))

  #18. check length and range of n_chain
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2, 1),
                                                                c(1, 1, NA))),
                                 n_chain = c(1, 1), multicore = FALSE),
               paste0("n_chain should be an integer > 0 and of length 1."))
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2, 1),
                                                                c(1, 1, NA))),
                                 n_chain = 0, multicore = FALSE),
               paste0("n_chain should be an integer > 0 and of length 1."))

  #19. check length and range of n_iter
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2, 1),
                                                                c(1, 1, NA))),
                                 n_iter = c(1, 1), multicore = FALSE),
               paste0("n_iter should be an integer > 0 and of length 1."))

  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2, 1),
                                                                c(1, 1, NA))),
                                 n_iter = 0, multicore = FALSE),
               paste0("n_iter should be an integer > 0 and of length 1."))

  #20. check length and range of n_warmup
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2, 1),
                                                                c(1, 1, NA))),
                                 n_warmup = c(1, 1), multicore = FALSE),
               "n_warmup should be an integer > 0 and of length 1.")

  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2, 1),
                                                                c(1, 1, NA))),
                                 n_warmup = 0, multicore = FALSE),
               "n_warmup should be an integer > 0 and of length 1.")

  #21. check length and range of thin
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2, 1),
                                                                c(1, 1, NA))),
                                 thin = c(1, 1), multicore = FALSE),
               "thin should be an integer > 0 and of length 1.")

  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2, 1),
                                                                c(1, 1, NA))),
                                 thin = 0, multicore = FALSE),
               "thin should be an integer > 0 and of length 1.")

  #22. check length and range of adapt_delta
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2, 1),
                                                                c(1, 1, NA))),
                                 adapt_delta = c(0.9, 0.9), multicore = FALSE),
               paste0("adapt_delta should be a numeric value > 0 and < 1 and ",
                      "of length 1."))

  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2, 1),
                                                                c(1, 1, NA))),
                                 adapt_delta = 1.2, multicore = FALSE),
               paste0("adapt_delta should be a numeric value > 0 and < 1 and ",
                      "of length 1."))

  #23. check length of seed
  expect_error(traditional_model(data = list(count = rbind(c(4, 1, 1),
                                                           c(1, 1, NA)),
                                             count_type = rbind(c(1, 2, 1),
                                                                c(1, 1, NA))),
                                 seed = c(1, 2), multicore = FALSE),
               "seed should be an integer of length 1.")


})
