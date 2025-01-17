test_that("joint_model input checks work - pt 1", {
  testthat::skip_on_cran()
  #' @srrstats {G5.2,G5.2b,BS2.15} Tests the assure function input checks are
  #'   behaving as expected.
  #1. input tags are valid, q = FALSE, cov = FALSE
  expect_error(joint_model(data = list(pcr.k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr.n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           multicore = FALSE),
               paste("Data should include 'pcr_n', 'pcr_k', and 'count'.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))
  #2. input tags are valid, q = FALSE, cov = TRUE
  site_cov <- cbind(c(1, 0), c(0.4, -0.4))
  colnames(site_cov) <- c("var_a", "var_b")
  expect_error(joint_model(data = list(pcr.k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr.n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       site.cov = site_cov),
                           cov = c("var_a", "var_b"),
                           multicore = FALSE),
               paste(paste0("Data should include 'pcr_n', 'pcr_k', ",
                            "'count', and 'site_cov'."),
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase2.html#prepare-the-data"),
                     sep = "\n"))
  #3. input tags are valid, q = TRUE, cov = FALSE
  expect_error(joint_model(data = list(pcr.k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr.n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count.type = rbind(c(1, 2, 1),
                                                          c(1, 2, NA))),
                           q = TRUE,
                           multicore = FALSE),
               paste(paste0("Data should include 'pcr_n', 'pcr_k', ",
                            "'count', and 'count_type'."),
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase3.html#prepare-the-data"),
                     sep = "\n"))
  #4. input tags are valid, q = TRUE, cov = TRUE
  site_cov <- cbind(c(1, 0), c(0.4, -0.4))
  colnames(site_cov) <- c("var_a", "var_b")
  expect_error(joint_model(data = list(pcr.k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr.n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count.type = rbind(c(1, 2, 1),
                                                          c(1, 2, NA)),
                                       site.cov = site_cov),
                           cov = c("var_a", "var_b"), q = TRUE,
                           multicore = FALSE),
               paste(paste0("Data should include 'pcr_n', 'pcr_k', 'count', ",
                            "'count_type', and 'site_cov'."),
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase3.html#prepare-the-data"),
                     sep = "\n"))

  #5. make sure dimensions of pcr_n and pcr_k are equal
  #' @srrstats {BS2.1a} Test to ensure pre-processing routines to ensure all
  #'   input data is dimensionally commensurate
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1, 1),
                                                     c(1, 1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           multicore = FALSE),
               paste("Dimensions of pcr_n and pcr_k do not match.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))

  #6. make sure dimensions of count and count_type are equal, if count_type is
  # present
  #' @srrstats {BS2.1a} Test to ensure pre-processing routines to ensure all
  #'   input data is dimensionally commensurate
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2), c(1, 2))),
                           q = TRUE,
                           multicore = FALSE),
               paste("Dimensions of count and count_type do not match.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase3.html#prepare-the-data"),
                     sep = "\n"))

  #7. make sure number of rows in count = number of rows in pcr_n and pcr_k
  #' @srrstats {BS2.1a} Test to ensure pre-processing routines to ensure all
  #'   input data is dimensionally commensurate
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA),
                                                     c(4, 1, 1))),
                           multicore = FALSE),
               paste(paste0("Number of sites \\(rows\\) in pcr data and ",
                            "traditional survey count data do not match."),
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))

  #8. make sure all data is numeric -- if q == TRUE
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c("NA", 2, 2),
                                                          c(1, 2, 2))),
                           q = TRUE,
                           multicore = FALSE),
               paste("Data should be numeric.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase3.html#prepare-the-data"),
                     sep = "\n"))

  #9. make sure all data is numeric -- if q == FALSE
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, "NA")),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           multicore = FALSE),
               paste("Data should be numeric.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))

  #10. make sure locations of NAs in count data match locations of NAs in
  # count_type data
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(NA, 2, 2),
                                                          c(1, 2, 2))),
                           q = TRUE,
                           multicore = FALSE),
               paste(paste0("Empty data cells \\(NA\\) in count data should ",
                            "match ",
                            "empty data cells \\(NA\\) in count_type data."),
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase3.html#prepare-the-data"),
                     sep = "\n"))

  #11. make sure locations of NAs in pcr_n data match locations of NAs in
  # pcr_k data
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, NA, 1)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           multicore = FALSE),
               paste(paste0("Empty data cells \\(NA\\) in pcr_n data should ",
                            "match ",
                            "empty data cells \\(NA\\) in pcr_k data."),
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))

  #12. make sure family is either "poisson", "negbin", or "gamma"
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           family = "normal",
                           multicore = FALSE),
               paste0("Invalid family. Options include 'poisson', ",
                      "'negbin', and 'gamma'."))

  #13. p10 priors is a vector of two numeric values
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           p10_priors = c(1, 1, 2),
                           multicore = FALSE),
               paste0("p10_priors should be a vector of two positive numeric ",
                      "values. ex. c\\(1,20\\)"))

  #14. p10 priors is a vector of two numeric values
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           multicore = FALSE,
                           p10_priors = "1,20"),
               paste0("p10_priors should be a vector of two positive numeric ",
                      "values. ex. c\\(1,20\\)"))

  #15. p10 priors is a vector of two numeric values
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           p10_priors = c(0, 20),
                           multicore = FALSE),
               paste0("p10_priors should be a vector of two positive numeric ",
                      "values. ex. c\\(1,20\\)"))

  #16. phi priors is a vector of two numeric values
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           phi_priors = c(0, 1), family = "negbin",
                           multicore = FALSE),
               paste0("phi_priors should be a vector of two positive numeric ",
                      "values. ex. c\\(0.25,0.25\\)"))

  #17. the smallest count_type is 1
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(0, 1, 2),
                                                          c(1, 2, NA))),
                           q = TRUE,
                           multicore = FALSE),
               paste(paste0("The first gear type should be referenced as 1 in ",
                            "count_type. Subsequent gear types should be ",
                            "referenced 2, 3, 4, etc."),
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase3.html#prepare-the-data"),
                     sep = "\n"))

  #18. count are integers
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4.1, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 1, 2),
                                                          c(1, 2, NA))),
                           q = TRUE, family  =  "negbin",
                           multicore = FALSE),
               paste0("All values in count should be non-negative integers. ",
                      "Use family = 'gamma' if count is continuous."))

  #19. pcr_n are integers
  expect_error(joint_model(data = list(pcr_k = rbind(c(0.99, 1, 1),
                                                     c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 1, 2),
                                                          c(1, 2, NA))),
                           q = TRUE,
                           multicore = FALSE),
               paste("All values in pcr_k should be non-negative integers.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))

  #20. pcr_k are integers
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3.1, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 1, 2),
                                                          c(1, 2, NA))),
                           q = TRUE,
                           multicore = FALSE),
               paste("All values in pcr_n should be non-negative integers.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))

  #21. count_type are integers
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1.1, 1, 2),
                                                          c(1, 2, NA))),
                           q = TRUE,
                           multicore = FALSE),
               paste("All values in count_type should be integers.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase3.html#prepare-the-data"),
                     sep = "\n"))

  #22. site_cov is numeric, if present
  site_cov <- cbind(c("high", "low"), c(0.4, -0.4))
  colnames(site_cov) <- c("var_a", "var_b")
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 2, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"), q = TRUE,
                           multicore = FALSE),
               paste("site_cov should be numeric.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app", "/usecase2.html"),
                     sep = "\n"))

  #23. cov values match column names in site_cov
  site_cov <- cbind(c(0, 1), c(0.4, -0.4))
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 2, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"), q = TRUE,
                           multicore = FALSE),
               paste(paste0("cov values should be listed in the column names ",
                            "of site_cov in the data."),
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app", "/usecase2.html"),
                     sep = "\n"))

  #24. site_cov has same number of rows as pcr_n and count, if present
  site_cov <- cbind(c(0, 1, 1), c(0.4, -0.4, 1))
  colnames(site_cov) <- c("var_a", "var_b")
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 2, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"), q = TRUE,
                           multicore = FALSE),
               paste(paste0("The number of rows in site_cov matrix should ",
                            "match the number of rows in all other matrices."),
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app", "/usecase2.html"),
                     sep = "\n"))

  #25. make sure count_type is not zero-length
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = matrix(NA, ncol = 3,
                                                           nrow = 0)),
                           q = TRUE,
                           multicore = FALSE),
               paste("count_type contains zero-length data.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase3.html#prepare-the-data"),
                     sep = "\n"))

  #26. make sure no column is entirely NA in count_type
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(4, 1, NA),
                                                          c(1, 1, NA))),
                           q = TRUE,
                           multicore = FALSE),
               paste("count_type contains a column with all NA.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase3.html#prepare-the-data"),
                     sep = "\n"))

  #27. make sure no column is entirely NA in pcr_k
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, NA), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           multicore = FALSE),
               paste("pcr_k contains a column with all NA.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))

  #28. make sure no column is entirely NA in pcr_n
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, NA), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           multicore = FALSE),
               paste("pcr_n contains a column with all NA.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))

  #29. make sure no column is entirely NA in count
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, NA), c(1, 1, NA))),
                           multicore = FALSE),
               paste("count contains a column with all NA.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))

  #30. make sure no data are undefined
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, Inf),
                                                     c(1, 1, NA))),
                           multicore = FALSE),
               paste("count contains undefined values \\(i.e., Inf or -Inf\\)",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))

  #31. make sure no data are undefined
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, Inf), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           multicore = FALSE),
               paste("pcr_n contains undefined values \\(i.e., Inf or -Inf\\)",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))

  #32. make sure no data are undefined
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, Inf), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           multicore = FALSE),
               paste("pcr_k contains undefined values \\(i.e., Inf or -Inf\\)",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))

  #33. make sure site_cov is not zero-length
  site_cov <- matrix(NA, ncol = 2, nrow = 0)
  colnames(site_cov) <- c("var_a", "var_b")
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"),
                           multicore = FALSE),
               paste("site_cov contains zero-length data.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app", "/usecase2.html"),
                     sep = "\n"))

  #34. make sure no column is entirely NA in site_cov
  site_cov <- rbind(c(4, 1, NA), c(1, 1, NA))
  colnames(site_cov) <- c("var_a", "var_b", "var_c")
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"),
                           multicore = FALSE),
               paste("site_cov contains a column with all NA.",
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app", "/usecase2.html"),
                     sep = "\n"))

  #35. make sure no data are undefined
  site_cov <- rbind(c(4, 1, Inf), c(1, 1, NA))
  colnames(site_cov) <- c("var_a", "var_b", "var_c")
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"),
                           multicore = FALSE),
               paste(paste0("site_cov contains undefined values \\(i.e., ",
                            "Inf or -Inf\\)"),
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app", "/usecase2.html"),
                     sep = "\n"))

  #36. length of initial values is equal to the number of chains
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      mu <- stats::runif(3, 0.01, 5),
      p10 <- exp(stats::runif(1, log(0.0001), log(0.08))),
      alpha <- rep(0.1, 3)
    )
  }
  site_cov <- rbind(c(4, 1), c(1, 1))
  colnames(site_cov) <- c("var_a", "var_b")
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"), initial_values = inits,
                           n_chain = 5,
                           multicore = FALSE),
               paste(paste0("The length of the list of initial values ",
                            "should equal the number of chains \\(n_chain, ",
                            "default is 4\\)."),
                     paste0("See the eDNAjoint guide for help formatting ",
                            "initial values: "),
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#initialvalues"),
                     sep = "\n"))

  #37. initial values check: if mu is numeric
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      mu <- stats::runif(3, -1, 0),
      p10 <- exp(stats::runif(1, log(0.0001), log(0.08))),
      alpha <- rep(0.1, 3)
    )
    names(inits[[i]]) <- c("mu", "p10", "alpha")
  }
  site_cov <- rbind(c(4, 1), c(1, 1))
  colnames(site_cov) <- c("var_a", "var_b")
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"), initial_values = inits,
                           multicore = FALSE),
               paste("Initial values for 'mu' should be numeric values > 0.",
                     paste0("See the eDNAjoint guide for help formatting ",
                            "initial values: "),
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#initialvalues"),
                     sep = "\n"))

  #38. initial values check: mu length
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      mu <- stats::runif(4, 0.1, 1),
      p10 <- exp(stats::runif(1, log(0.0001), log(0.08))),
      alpha <- rep(0.1, 3)
    )
    names(inits[[i]]) <- c("mu", "p10", "alpha")
  }
  site_cov <- rbind(c(4, 1), c(1, 1))
  colnames(site_cov) <- c("var_a", "var_b")
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"), initial_values = inits,
                           multicore = FALSE),
               paste(paste0("The length of initial values for 'mu' should ",
                            "equal the number of sites."),
                     paste0("See the eDNAjoint guide for help formatting ",
                            "initial values: "),
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#initialvalues"),
                     sep = "\n"))


  #39. initial values check: p10 numeric
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      mu <- stats::runif(2, 0, 1),
      p10 <- "-1",
      alpha <- rep(0.1, 3)
    )
    names(inits[[i]]) <- c("mu", "p10", "alpha")
  }
  site_cov <- rbind(c(4, 1), c(1, 1))
  colnames(site_cov) <- c("var_a", "var_b")
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"), initial_values = inits,
                           multicore = FALSE),
               paste("Initial values for 'p10' should be numeric.",
                     paste0("See the eDNAjoint guide for help formatting ",
                            "initial values: "),
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#initialvalues"),
                     sep = "\n"))

})

test_that("joint_model input checks work - pt 2", {
  testthat::skip_on_cran()

  #40. initial values check: p10 length
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      mu <- stats::runif(2, 0, 1),
      p10 <- exp(stats::runif(2, log(0.0001), log(0.08))),
      alpha <- rep(0.1, 3)
    )
    names(inits[[i]]) <- c("mu", "p10", "alpha")
  }
  site_cov <- rbind(c(4, 1), c(1, 1))
  colnames(site_cov) <- c("var_a", "var_b")
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"), initial_values = inits,
                           multicore = FALSE),
               paste("The length of initial values for 'p10' should equal 1.",
                     paste0("See the eDNAjoint guide for help formatting ",
                            "initial values: "),
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#initialvalues"),
                     sep = "\n"))

  #41. initial values check: alpha numeric
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      mu <- stats::runif(2, 0, 1),
      p10 <- exp(stats::runif(1, log(0.0001), log(0.08))),
      alpha <- "1"
    )
    names(inits[[i]]) <- c("mu", "p10", "alpha")
  }
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           initial_values = inits,
                           multicore = FALSE),
               paste("Initial values for 'alpha' should be numeric.",
                     paste0("See the eDNAjoint guide for help formatting ",
                            "initial values: "),
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#initialvalues"),
                     sep = "\n"))

  #42. initial values check: alpha length
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      mu <- stats::runif(2, 0, 1),
      p10 <- exp(stats::runif(1, log(0.0001), log(0.08))),
      alpha <- c(1, 0)
    )
    names(inits[[i]]) <- c("mu", "p10", "alpha")
  }
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           initial_values = inits,
                           multicore = FALSE),
               paste("The length of initial values for 'alpha' should equal 1.",
                     paste0("See the eDNAjoint guide for help formatting ",
                            "initial values: "),
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#initialvalues"),
                     sep = "\n"))

  #43. initial values check: alpha numeric
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      mu <- stats::runif(2, 0, 1),
      p10 <- exp(stats::runif(1, log(0.0001), log(0.08))),
      alpha <- c("1", "2")
    )
    names(inits[[i]]) <- c("mu", "p10", "alpha")
  }
  site_cov <- rbind(c(4, 1), c(1, 1))
  colnames(site_cov) <- c("var_a", "var_b")
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"), initial_values = inits,
                           multicore = FALSE),
               paste("Initial values for 'alpha' should be numeric.",
                     paste0("See the eDNAjoint guide for help formatting ",
                            "initial values: "),
                     paste0("https://ednajoint.netlify.app",
                            "/usecase2.html#initialvalues"),
                     sep = "\n"))

  #44. initial values check: alpha length
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      mu <- stats::runif(2, 0, 1),
      p10 <- exp(stats::runif(1, log(0.0001), log(0.08))),
      alpha <- rep(0.1, 2)
    )
    names(inits[[i]]) <- c("mu", "p10", "alpha")
  }
  site_cov <- rbind(c(4, 1), c(1, 1))
  colnames(site_cov) <- c("var_a", "var_b")
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       site_cov = site_cov),
                           cov = c("var_a", "var_b"), initial_values = inits,
                           multicore = FALSE),
               paste(paste0("The length of initial values for 'alpha' should ",
                            "equal\\: \\# covariates \\+ 1 \\(i.e., ",
                            "including intercept\\)."),
                     paste0("See the eDNAjoint guide for help formatting ",
                            "initial values: "),
                     paste0("https://ednajoint.netlify.app",
                            "/usecase2.html#initialvalues"),
                     sep = "\n"))

  #45. initial values check: q numeric
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      q <- "0.1"
    )
    names(inits[[i]]) <- c("q")
  }
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           initial_values = inits,
                           multicore = FALSE),
               paste("Initial values for 'q' should be numeric.",
                     paste0("See the eDNAjoint guide for help formatting ",
                            "initial values: "),
                     paste0("https://ednajoint.netlify.app",
                            "/usecase2.html#initialvalues"),
                     sep = "\n"))

  #46. initial values check: q length
  n_chain <- 4
  inits <- list()
  for (i in 1:n_chain) {
    inits[[i]] <- list(
      q <- c(0.1, 0.1)
    )
    names(inits[[i]]) <- c("q")
  }
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           initial_values = inits,
                           multicore = FALSE),
               paste(paste0("The length of initial values for 'q' should ",
                            "equal: \\# unique gear types \\- 1 \\(i.e., q ",
                            "for reference type = 1\\)."),
                     paste0("See the eDNAjoint guide for help formatting ",
                            "initial values: "),
                     paste0("https://ednajoint.netlify.app",
                            "/usecase2.html#initialvalues"),
                     sep = "\n"))

  #47. check length and range of n_chain
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           n_chain = c(1, 1), multicore = FALSE),
               paste0("n_chain should be an integer > 0 and of length 1."))
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           n_chain = 0, multicore = FALSE),
               paste0("n_chain should be an integer > 0 and of length 1."))

  #48. check length and range of n_iter
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           n_iter = c(1, 1), multicore = FALSE),
               paste0("n_iter should be an integer > 0 and of length 1."))

  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           n_iter = 0, multicore = FALSE),
               paste0("n_iter should be an integer > 0 and of length 1."))

  #49. check length and range of n_warmup
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           n_warmup = c(1, 1), multicore = FALSE),
               "n_warmup should be an integer > 0 and of length 1.")

  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           n_warmup = 0, multicore = FALSE),
               "n_warmup should be an integer > 0 and of length 1.")

  #50. check length and range of thin
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           thin = c(1, 1), multicore = FALSE),
               "thin should be an integer > 0 and of length 1.")

  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           thin = 0, multicore = FALSE),
               "thin should be an integer > 0 and of length 1.")

  #51. check length and range of adapt_delta
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           adapt_delta = c(0.9, 0.9), multicore = FALSE),
               paste0("adapt_delta should be a numeric value > 0 and < 1 and ",
                      "of length 1."))

  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           adapt_delta = 1.2, multicore = FALSE),
               paste0("adapt_delta should be a numeric value > 0 and < 1 and ",
                      "of length 1."))

  #52. check length of seed
  expect_error(joint_model(data = list(pcr_k = rbind(c(1, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       count_type = rbind(c(1, 2, 1),
                                                          c(1, 1, NA))),
                           seed = c(1, 2), multicore = FALSE),
               "seed should be an integer of length 1.")

  #53. check K <= N
  expect_error(joint_model(data = list(pcr_k = rbind(c(4, 1, 1), c(1, 1, NA)),
                                       pcr_n = rbind(c(3, 3, 3), c(3, 3, NA)),
                                       count = rbind(c(4, 1, 1), c(1, 1, NA))),
                           multicore = FALSE),
               paste(paste0("N should be >= K in pcr data. N is the number ",
                            "of pcr replicates per sample, and K is the ",
                            "number of positive detections among replicates."),
                     "See the eDNAjoint guide for data formatting help: ",
                     paste0("https://ednajoint.netlify.app",
                            "/usecase1.html#prepare-the-data"),
                     sep = "\n"))


})


# correctness and parameter recovery tests
#' @srrstats {G5.4, G5.6} Correctness/parameter recovery tests to test that
#'   the implementation produces expected results given data with known
#'   properties
#' @srrstats {PD4.0} These tests do not test for numeric equality of outputs,
#'   but rather test for the recovery of parameter values. In these tests, I
#'   simulate data with known parameter values, run the main statistical
#'   function (joint_model()) with the data, and then test if the known
#'   parameter values are within the 95% credibility interval of the parameters'
#'   posteriors in the function output.
test_that("joint_model parameter recovery tests work - pt 1", {
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
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability of
  # eDNA detection)
  p11 <- mu / (mu + exp(mat_site %*% alpha))
  p <- pmin(p11 + exp(log_p10), 1)

  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    count[i, ] <- rpois(nobs_count, mu[i])
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    site_cov = mat_site
  )
  # run model
  fit1 <- suppressWarnings({
    joint_model(data = data, cov = c("var_a", "var_b"),
                multicore = FALSE, seed = 10)
  })
  # summary
  summary1 <- as.data.frame(rstan::summary(fit1$model,
                                           pars = c("mu", "alpha", "log_p10"),
                                           probs = c(0.025, 0.975))$summary)

  # set up empty vector to check if true values are in 95% interval of
  # posterior estimates
  par <- rep(NA, length(alpha))
  # check alpha
  for (i in seq_along(alpha)) {
    par[i] <- paste0("alpha[", i, "]")
  }
  check <- c(alpha > summary1[par, "2.5%"], alpha < summary1[par, "97.5%"],
             log_p10 > summary1["log_p10", "2.5%"],
             log_p10 < summary1["log_p10", "97.5%"])

  #' @srrstats {G3.0, G5.6a} Instead of comparing floating point values for
  #'   equality, here the model is tested to determine if the true parameter
  #'   values are within the 95% quantiles of the posterior
  # all should be equal to true
  expect_equal(check, rep(TRUE, length(check)))

  # test that output values are on the same scale as the data
  mu_estimates <- rep(NA, nsite)
  for (i in 1:nsite) {
    mu_estimates[i] <- summary1[paste0("mu[", i, ",1]"), "mean"]
  }

  # get mean of input count data at each site
  data_mean <- rowMeans(data$count)
  #' @srrstats {BS7.4, BS7.4a} Check to ensure that mean posterior estimates
  #'   are on the same scale as the mean of the input data (here checking
  #'   estimates of mu, i.e., expected catch rate at each site)
  # check that estimates are on same scale as data
  expect_equal(round(mu_estimates, 0), round(data_mean, 0))

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
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(mat_site %*% alpha))
  p <- pmin(p11 + exp(log_p10), 1)

  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
    count[i, ] <- rpois(nobs_count, mu[i])
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    site_cov = mat_site
  )
  # run model
  fit_large <- suppressWarnings({
    joint_model(data = data, cov = c("var_a", "var_b"),
                multicore = FALSE, seed = 10)
  })
  # summary
  summary_large <- as.data.frame(rstan::summary(fit_large$model,
                                                pars = "alpha",
                                                probs = c(0.025,
                                                          0.975))$summary)

  # compare standard errors of estimates of alpha coefficients with
  # implementation with a smaller dataset
  se_large <- c(
    summary_large["alpha[2]", "se_mean"], summary_large["alpha[3]", "se_mean"]
  )

  # get values for smaller dataset
  se_small <- c(
    summary1["alpha[2]", "se_mean"], summary1["alpha[3]", "se_mean"]
  )

  # set up empty vector to check if standard errors are smaller with
  # larger dataset
  check <- se_large < se_small

  #' @srrstats {G5.7} Check to see that implementation performs as expected
  #'   properties of data change (i.e., standard error of posteriors is smaller
  #'   if there are more data observations)
  # all should be equal to true
  expect_equal(check, rep(TRUE, length(check)))

})

test_that("joint_model parameter recovery tests work - pt 2", {
  testthat::skip_on_cran()

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
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    count[i, ] <- rpois(nobs_count, mu[i])
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability of
  # eDNA detection)
  p11 <- mu / (mu + exp(mat_site %*% alpha))
  p <- pmin(p11 + exp(log_p10), 1)

  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    site_cov = mat_site
  )
  # run model
  fit2 <- suppressWarnings({
    joint_model(data = data, cov = c("var_a", "var_b"),
                multicore = FALSE, seed = 2)
  })
  # summary
  summary2 <- as.data.frame(rstan::summary(fit2$model,
                                           pars = c("mu", "alpha", "log_p10"),
                                           probs = c(0.025, 0.975))$summary)

  # set up empty vector to check if true values are in 95% interval of
  # posterior estimates
  par <- rep(NA, length(alpha))
  # check alpha
  for (i in seq_along(alpha)) {
    par[i] <- paste0("alpha[", i, "]")
  }
  check <- c(alpha > summary2[par, "2.5%"], alpha < summary2[par, "97.5%"],
             log_p10 > summary2["log_p10", "2.5%"],
             log_p10 < summary2["log_p10", "97.5%"])

  #' @srrstats {G3.0, G5.6a} Instead of comparing floating point values for
  #'   equality, here the model is tested to determine if the true parameter
  #'   values are within the 95% quantiles of the posterior
  # all should be equal to true
  expect_equal(check, rep(TRUE, length(check)))

})


test_that("joint_model probability distribution tests work", {
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
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5
  phi <- 0.1
  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    count[i, ] <- rnbinom(nobs_count, mu = mu[i], size = phi)
  }
  # p11 (probability of true positive eDNA detection) and p (probability of
  # eDNA detection)
  p11 <- mu / (mu + exp(beta))
  p <- pmin(p11 + exp(log_p10), 1)

  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count
  )

  ##############################################################
  # fit data with a negative binomial and poisson distribution #

  # run model -- negative binomial distribution
  negbin_fit <- suppressWarnings({
    joint_model(data = data, multicore = FALSE, family = "negbin", seed = 2)
  })
  # run model -- poisson distribution
  pois_fit <- suppressWarnings({
    joint_model(data = data, multicore = FALSE, family = "poisson", seed = 2)
  })

  # summarize outputs
  negbin_summary <- as.data.frame(rstan::summary(negbin_fit$model,
                                                 pars = "mu",
                                                 probs = c(0.025,
                                                           0.975))$summary)
  pois_summary <- as.data.frame(rstan::summary(pois_fit$model,
                                               pars = "mu",
                                               probs = c(0.025,
                                                         0.975))$summary)

  # summarize differences in mean estimates of mu
  summary_table <- table(negbin_summary$mean > pois_summary$mean)

  # calculate percentage of mean estimates of mu that are greater when
  # a negative binomial distribution are used rather than a poisson distribution
  percent <- summary_table[[2]] / (summary_table[[2]] + summary_table[[1]])


  # test that estimates of mu when a negative binomial distribution is used is
  # greater than estimates of mu when a poisson distribution is used for
  # >80% of sites
  expect_gt(percent, 0.8)


})

test_that("semi-paired model works - p10, q, phi, beta", {
  testthat::skip_on_cran()

  ## 1.
  # model includes "p10", "q", "phi", "beta"

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
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rnbinom(n = 1, mu = mu[i], size = phi)
      } else {
        count[i, j] <- rnbinom(n = 1, mu = mu[i] * q, size = phi)
      }
    }
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(beta))
  p <- pmin(p11 + exp(log_p10), 1)
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2, 4), ] <- NA
  count_type[c(2, 4), ] <- NA

  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    count_type = count_type
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
  names(inits[[1]]) <- c("mu", "p10", "alpha", "phi", "q")

  # run model
  fit <- suppressWarnings({
    joint_model(data = data, q = TRUE, family = "negbin",
                initial_values = inits, n_warmup = 25,
                n_iter = 100, n_chain = 1, multicore = FALSE, seed = 10)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "q[1]", "phi", "alpha[1]") %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model, par = "mu")$summary[, 1])
  expect_equal(length(mu_estimates), nsite * 2)
  expect_true(is.numeric(mu_estimates))

})

test_that("semi-paired model works - p10, beta, q", {
  testthat::skip_on_cran()


  ## 2.
  # model includes "p10", "beta", "q"

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
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rpois(1, mu[i])
      } else {
        count[i, j] <- rpois(1, mu[i] * q)
      }
    }
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(beta))
  p <- pmin(p11 + exp(log_p10), 1)
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2, 4), ] <- NA
  count_type[c(2, 4), ] <- NA

  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    count_type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = beta,
    q = 2
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha", "q")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, q = TRUE, n_chain = 1, multicore = FALSE,
                seed = 10, initial_values = inits, n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "q[1]", "alpha[1]") %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model, par = "mu")$summary[, 1])
  expect_equal(length(mu_estimates), nsite * 2)
  expect_true(is.numeric(mu_estimates))

})

test_that("semi-paired model works - gamma, p10, beta, q", {
  testthat::skip_on_cran()

  ## 3.
  # model includes "p10", "beta", "q", "alpha_gamma", "beta_gamma"

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
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rgamma(1, shape = alpha_gamma[i], rate = beta_gamma)
      } else {
        count[i, j] <- rgamma(1, shape = alpha_gamma[i] * q, rate = beta_gamma)
      }
    }
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(beta))
  p <- pmin(p11 + exp(log_p10), 1)
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2, 4), ] <- NA
  count_type[c(2, 4), ] <- NA

  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    count_type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1, length(mu)),
    p10 = exp(log_p10),
    alpha = beta,
    q = q
  )
  names(inits[[1]]) <- c("alpha_gamma", "beta_gamma", "p10", "alpha", "q")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, q = TRUE, family = "gamma",
                n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits, n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "q[1]", "alpha[1]") %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model, par = "mu")$summary[, 1])
  expect_equal(length(mu_estimates), nsite * 2)
  expect_true(is.numeric(mu_estimates))

})

test_that("semi-paired model works - p10, phi, beta", {
  testthat::skip_on_cran()

  ## 4.
  # model includes "p10", "phi", "beta"

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
  for (i in 1:nsite) {
    count[i, ] <- rnbinom(n = nobs_count, mu = mu[i], size = phi)
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(beta))
  p <- pmin(p11 + exp(log_p10), 1)
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2, 4), ] <- NA

  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count
  )

  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = beta,
    phi = phi
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha", "phi")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, family = "negbin", initial_values = inits,
                n_chain = 1, multicore = FALSE, seed = 10,
                n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "phi", "alpha[1]") %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model, par = "mu")$summary[, 1])
  expect_equal(length(mu_estimates), nsite)
  expect_true(is.numeric(mu_estimates))

})

test_that("semi-paired model works - p10, beta", {
  testthat::skip_on_cran()


  ## 5.
  # model includes "p10", "beta"

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
  for (i in 1:nsite) {
    count[i, ] <- rpois(nobs_count, mu[i])
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(beta))
  p <- pmin(p11 + exp(log_p10), 1)
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2, 4), ] <- NA

  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = beta
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, initial_values = inits,
                n_chain = 1, multicore = FALSE, seed = 10,
                n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "alpha[1]") %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model, par = "mu")$summary[, 1])
  expect_equal(length(mu_estimates), nsite)
  expect_true(is.numeric(mu_estimates))

})

test_that("semi-paired model works - gamma, p10, beta", {
  testthat::skip_on_cran()


  ## 6.
  # model includes "p10", "beta", "alpha_gamma", "alpha_beta"

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
  for (i in 1:nsite) {
    count[i, ] <- rgamma(nobs_count, shape = alpha_gamma[i], rate = beta_gamma)
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(beta))
  p <- pmin(p11 + exp(log_p10), 1)
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2, 4), ] <- NA

  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1, length(mu)),
    p10 = exp(log_p10),
    alpha = beta
  )
  names(inits[[1]]) <- c("alpha_gamma", "beta_gamma", "p10", "alpha")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, initial_values = inits,
                family = "gamma",
                n_chain = 1, multicore = FALSE, seed = 10,
                n_warmup = 25,
                n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "alpha[1]") %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model, par = "mu")$summary[, 1])
  expect_equal(length(mu_estimates), nsite)
  expect_true(is.numeric(mu_estimates))

})

test_that("semi-paired model works - p10, q, phi, alpha", {
  testthat::skip_on_cran()


  ## 7.
  # model includes "p10", "q", "phi", "alpha"

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
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rnbinom(n = 1, mu = mu[i], size = phi)
      } else {
        count[i, j] <- rnbinom(n = 1, mu = mu[i] * q, size = phi)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(mat_site %*% alpha))
  p <- pmin(p11 + exp(log_p10), 1)
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2, 4), ] <- NA
  count_type[c(2, 4), ] <- NA

  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    count_type = count_type,
    site_cov = mat_site
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
  names(inits[[1]]) <- c("mu", "p10", "alpha", "q", "phi"
  )
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, family = "negbin", q = TRUE,
                cov = c("var_a", "var_b"),
                n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits,
                adapt_delta = 0.99,
                n_warmup = 25,
                n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "q[1]", "phi", "alpha[1]",
                    "alpha[2]", "alpha[3]") %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model, par = "mu")$summary[, 1])
  expect_equal(length(mu_estimates), nsite * 2)
  expect_true(is.numeric(mu_estimates))

})

test_that("semi-paired model works - p10, alpha, q", {
  testthat::skip_on_cran()


  ## 8.
  # model includes "p10", "alpha", "q"

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
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rpois(n = 1, mu[i])
      } else {
        count[i, j] <- rpois(n = 1, mu[i] * q)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(mat_site %*% alpha))
  p <- pmin(p11 + exp(log_p10), 1)
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # add NA to sites 2 and 4
  count[c(2, 4), ] <- NA
  count_type[c(2, 4), ] <- NA

  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    count_type = count_type,
    site_cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha"
  )
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, q = TRUE,
                cov = c("var_a", "var_b"),
                n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits, n_warmup = 25,
                n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "q[1]", "alpha[1]",
                    "alpha[2]", "alpha[3]") %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model, par = "mu")$summary[, 1])
  expect_equal(length(mu_estimates), nsite * 2)
  expect_true(is.numeric(mu_estimates))

})

test_that("semi-paired model works - gamma, p10, alpha, q", {
  testthat::skip_on_cran()


  ## 9.
  # model includes "p10", "alpha", "q", "alpha_gamma", "alpha_beta"

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
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rgamma(1, shape = alpha_gamma[i], rate = beta_gamma)
      } else {
        count[i, j] <- rgamma(1, shape = alpha_gamma[i] * q, rate = beta_gamma)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(mat_site %*% alpha))
  p <- pmin(p11 + exp(log_p10), 1)
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2, 4), ] <- NA
  count_type[c(2, 4), ] <- NA

  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    count_type = count_type,
    site_cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1, length(mu)),
    p10 = exp(log_p10),
    alpha = alpha,
    q = q
  )
  names(inits[[1]]) <- c("alpha_gamma", "beta_gamma", "p10", "alpha", "q")

  # run model
  fit <- suppressWarnings({
    joint_model(data = data, q = TRUE, family = "gamma",
                cov = c("var_a", "var_b"),
                n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits, n_warmup = 25,
                n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "q[1]", "alpha[1]") %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model, par = "mu")$summary[, 1])
  expect_equal(length(mu_estimates), nsite * 2)
  expect_true(is.numeric(mu_estimates))

})

test_that("semi-paired model works - p10, phi, alpha", {
  testthat::skip_on_cran()


  ## 10.
  # model includes "p10", "phi", "alpha"

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
  for (i in 1:nsite) {
    count[i, ] <- rnbinom(n = nobs_count, mu = mu[i], size = phi)
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(mat_site %*% alpha))
  p <- pmin(p11 + exp(log_p10), 1)
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2, 4), ] <- NA

  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    site_cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha,
    phi = phi
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha",
                         "phi")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, family = "negbin",
                cov = c("var_a", "var_b"), n_warmup = 25, n_iter = 100,
                n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "phi", "alpha[1]") %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model, par = "mu")$summary[, 1])
  expect_equal(length(mu_estimates), nsite)
  expect_true(is.numeric(mu_estimates))

})

test_that("semi-paired model works - p10, alpha", {
  testthat::skip_on_cran()


  ## 11.
  # model includes "p10", "alpha"

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
  for (i in 1:nsite) {
    count[i, ] <- rpois(nobs_count, mu[i])
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(mat_site %*% alpha))
  p <- pmin(p11 + exp(log_p10), 1)
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2, 4), ] <- NA

  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    site_cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha"
  )
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, cov = c("var_a", "var_b"),
                n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits, n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "alpha[1]",
                    "alpha[2]", "alpha[3]") %in% output_params))

  # check length of of mu
  mu_estimates <- as.vector(rstan::summary(fit$model, par = "mu")$summary[, 1])
  expect_equal(length(mu_estimates), nsite)
  expect_true(is.numeric(mu_estimates))

})

test_that("semi-paired model works - gamma, p10, alpha", {
  testthat::skip_on_cran()


  ## 12.
  # model includes "p10", "alpha", "alpha_gamma", "beta_gamma"

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
  for (i in 1:nsite) {
    count[i, ] <- rgamma(nobs_count, shape = alpha_gamma[i], rate = beta_gamma)
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- mu / (mu + exp(mat_site %*% alpha))
  p <- pmin(p11 + exp(log_p10), 1)
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(3, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }

  # add NA to sites 2 and 4
  count[c(2, 4), ] <- NA

  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    site_cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1, length(mu)),
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha
  )
  names(inits[[1]]) <- c("alpha_gamma", "beta_gamma", "mu", "p10", "alpha"
  )
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, family = "gamma",
                cov = c("var_a", "var_b"),
                n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits, n_warmup = 25,
                n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "alpha[1]",
                    "alpha[2]", "alpha[3]") %in% output_params))

  # check length of mu
  mu_estimates <- as.vector(rstan::summary(fit$model, par = "mu")$summary[, 1])
  expect_equal(length(mu_estimates), nsite)
  expect_true(is.numeric(mu_estimates))


})
