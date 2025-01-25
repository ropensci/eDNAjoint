## helper functions for joint_model.R ##

####################################
# helper functions: get model type #
####################################

#' @noRd
is_catch_type <- function(q) {
  out <- ifelse(q == TRUE, TRUE, FALSE)
  return(out)
}
#' @noRd
get_family_index <- function(family) {
  if (family == "poisson") {
    index <- 1
  } else if (family == "negbin") {
    index <- 2
  } else if (family == "gamma") {
    index <- 3
  }
  return(index)
}

####################################
# helper functions: get stan model #
####################################

#' @noRd
get_stan_model <- function(family) {

  index <- ifelse(family %in% c("poisson", "negbin"), 1, 2)

  return(index)
}


####################################
# helper functions: initial values #
####################################
#' @srrstats {BS2.7, BS2.11} Option for user to provide initial values for each
#'   chain

#' @noRd
get_inits <- function(n_chain, pcr_all, initial_values, cov, l_match_trad,
                      l_match_dna, data, q_names) {
  if (!is.null(q_names)) {
    inits <- init_joint_cov_catchability(n_chain, pcr_all, q_names, cov,
                                         initial_values, l_match_trad,
                                         l_match_dna, data)
  } else {
    inits <- init_joint_cov(n_chain, pcr_all, cov, initial_values,
                            l_match_trad, l_match_dna, data)
  }

  return(inits)
}

#' @noRd
init_joint_cov <- function(n_chain, pcr_all, cov, initial_values,
                           l_match_trad, l_match_dna, data) {

  # get mu means
  mu_means_trad <- as.vector(stats::na.omit(rowMeans(data$count,
                                                     na.rm = TRUE) + 0.01))
  mu_means_all <- rep(NA, dim(l_match_dna)[1] + dim(l_match_trad)[1])
  mu_means_all[l_match_trad$L_ind] <- mu_means_trad
  if (dim(l_match_dna)[1] > 0) {
    mu_means_all[l_match_dna$L_ind] <- rep(mean(mu_means_trad),
                                           dim(l_match_dna)[1])
  }

  # helper function
  # joint model, catchability coefficient, site covariates
  inits_list <- list()
  if (all(!is.null(initial_values))) {
    for (i in 1:n_chain) {
      inits_list[[i]] <- list(
        mu_trad = if ("mu" %in% names(initial_values[[i]])) {
          initial_values[[i]]$mu[l_match_trad$L_ind]
        } else {
          mu_means_trad
        },
        mu = if ("mu" %in% names(initial_values[[i]])) {
          initial_values[[i]]$mu
        } else {
          mu_means_all
        },
        log_p10 = if ("p10" %in% names(initial_values[[i]])) {
          log(initial_values[[i]]$p10)
        } else {
          stats::runif(1, log(0.0001), log(0.01))
        },
        alpha = if ("alpha" %in% names(initial_values[[i]])) {
          as.array(initial_values[[i]]$alpha)
        } else {
          as.array(c(3.5, rep(0, length(cov))))
        },
        p_dna = rep(0.4, dim(l_match_dna)[1]),
        p11_dna = rep(0.4, dim(l_match_dna)[1]) - 0.01
      )
    }
  } else {
    for (i in 1:n_chain) {
      inits_list[[i]] <- list(
        mu_trad = mu_means_trad,
        mu = mu_means_all,
        log_p10 = stats::runif(1, log(0.0001), log(0.01)),
        alpha = as.array(c(3.5, rep(0, length(cov)))),
        p_dna = rep(0.4, dim(l_match_dna)[1]),
        p11_dna = rep(0.4, dim(l_match_dna)[1]) - 0.01
      )
    }
  }

  return(inits_list)
}

#' @noRd
init_joint_cov_catchability <- function(n_chain, pcr_all, q_names, cov,
                                        initial_values, l_match_trad,
                                        l_match_dna, data) {

  # get mu means
  mu_means_trad <- as.vector(stats::na.omit(rowMeans(data$count,
                                                     na.rm = TRUE) + 0.01))
  mu_means_all <- rep(NA, dim(l_match_dna)[1] + dim(l_match_trad)[1])
  mu_means_all[l_match_trad$L_ind] <- mu_means_trad
  if (dim(l_match_dna)[1] > 0) {
    mu_means_all[l_match_dna$L_ind] <- rep(mean(mu_means_trad),
                                           dim(l_match_dna)[1])
  }

  # helper function
  # joint model, catchability coefficient, site covariates
  inits_list <- list()
  if (all(!is.null(initial_values))) {
    for (i in 1:n_chain) {
      inits_list[[i]] <- list(
        mu_trad = if ("mu" %in% names(initial_values[[i]])) {
          initial_values[[i]]$mu[l_match_trad$L_ind]
        } else {
          mu_means_trad
        },
        mu = if ("mu" %in% names(initial_values[[i]])) {
          initial_values[[i]]$mu
        } else {
          mu_means_all
        },
        log_p10 = if ("p10" %in% names(initial_values[[i]])) {
          log(initial_values[[i]]$p10)
        } else {
          stats::runif(1, log(0.0001), log(0.01))
        },
        alpha = if ("alpha" %in% names(initial_values[[i]])) {
          as.array(initial_values[[i]]$alpha)
        } else {
          as.array(c(3.5, rep(0, length(cov))))
        },
        q_trans = if ("q" %in% names(initial_values[[i]])) {
          as.data.frame(initial_values[[i]]$q)
        } else {
          as.data.frame(stats::runif(length(q_names), 0.01, 1))
        },
        p_dna = rep(0.4, dim(l_match_dna)[1]),
        p11_dna = rep(0.4, dim(l_match_dna)[1]) - 0.01
      )
    }
  } else {
    for (i in 1:n_chain) {
      inits_list[[i]] <- list(
        mu_trad = mu_means_trad,
        mu = mu_means_all,
        log_p10 = stats::runif(1, log(0.0001), log(0.01)),
        alpha = as.array(c(3.5, rep(0, length(cov)))),
        q_trans = as.data.frame(stats::runif(length(q_names), 0.01, 1)),
        p_dna = rep(0.4, dim(l_match_dna)[1]),
        p11_dna = rep(0.4, dim(l_match_dna)[1]) - 0.01
      )
    }
  }

  return(inits_list)

}



################
# input checks #
################
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#'   messages

#' @noRd
# input checks if catchabilty coefficients are used
catchability_checks_1 <- function(data, cov) {

  #' @srrstats {G2.13} Pre-processing routines to check for missing data
  ## All tags in data are valid (i.e., include pcr_n, pcr_k, count,
  ## count_type, and site_cov)
  # no covariates
  if (all(is.null(cov)) && !all(c("pcr_n", "pcr_k",
                                  "count", "count_type") %in% names(data))) {
    err_msg1 <- paste0("Data should include 'pcr_n', 'pcr_k', ",
                       "'count', and 'count_type'.")
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase3.html#prepare-the-data")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }
  # q=TRUE and cov != "None"
  if (all(!is.null(cov)) && !all(c("pcr_n", "pcr_k", "count",
                                   "count_type",
                                   "site_cov") %in% names(data))) {
    err_msg1 <- paste0("Data should include 'pcr_n', 'pcr_k', ",
                       "'count', 'count_type', and 'site_cov'.")
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase3.html#prepare-the-data")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }

  ## make sure count_type is not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for zero-length
  #'   data
  if (dim(data$count_type)[1] == 0) {
    err_msg1 <- "count_type contains zero-length data."
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase3.html#prepare-the-data")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }
  ## make sure no column is entirely NA in count_type
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with
  #' all NA
  if (any(apply(data$count_type, 2, function(col) all(is.na(col))))) {
    err_msg1 <- "count_type contains a column with all NA."
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase3.html#prepare-the-data")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }

  ## make sure dimensions of count and count_type are equal, if
  ## count_type is present
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #'   is dimensionally commensurate
  if (dim(data$count)[1] != dim(data$count_type)[1] ||
        dim(data$count)[2] != dim(data$count_type)[2]) {
    err_msg1 <- "Dimensions of count and count_type do not match."
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase3.html#prepare-the-data")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }
}

catchability_checks_2 <- function(data, cov) {

  ## make sure all data is numeric -- if q == TRUE
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if (is.numeric(data$pcr_k) == FALSE ||
        is.numeric(data$pcr_n) == FALSE ||
        is.numeric(data$count) == FALSE ||
        is.numeric(data$count_type) == FALSE) {
    err_msg1 <- "Data should be numeric."
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase3.html#prepare-the-data")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }
  ## make sure locations of NAs in count data match locations of NAs in
  ## count_type data
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input
  #'   data is dimensionally commensurate
  if (
    any((which(is.na(data$count)) == which(is.na(data$count_type))) == FALSE)
  ) {
    err_msg1 <- paste0("Empty data cells (NA) in count data should match ",
                       "empty data cells (NA) in count_type data.")
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase3.html#prepare-the-data")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }
  ## the smallest count_type is 1
  if (min(data$count_type, na.rm =  TRUE) != 1) {
    err_msg1 <- paste0("The first gear type should be referenced as 1 in ",
                       "count_type. Subsequent gear types should be ",
                       "referenced 2, 3, 4, etc.")
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase3.html#prepare-the-data")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }

  ## count_type are integers
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if (!all(data$count_type %% 1 %in% c(0, NA))) {
    err_msg1 <- "All values in count_type should be integers."
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase3.html#prepare-the-data")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }
}

#' @noRd
# input checks if no catchabilty coefficients are used
no_catchability_checks <- function(data, cov) {

  #' @srrstats {G2.13} Pre-processing routines to check for missing data
  ## All tags in data are valid (i.e., include pcr_n, pcr_k, count,
  ## and site_cov)
  # no covariates
  if (all(is.null(cov)) &&
        !all(c("pcr_n", "pcr_k", "count") %in% names(data))) {
    err_msg1 <- "Data should include 'pcr_n', 'pcr_k', and 'count'."
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase1.html#prepare-the-data")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }
  # with covariates
  if (all(!is.null(cov)) &&
        !all(c("pcr_n", "pcr_k", "count", "site_cov") %in% names(data))) {
    err_msg1 <- "Data should include 'pcr_n', 'pcr_k', 'count', and 'site_cov'."
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase2.html#prepare-the-data")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }

  ## make sure all data is numeric -- if q == FALSE
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if (is.numeric(data$pcr_k) == FALSE ||
        is.numeric(data$pcr_n) == FALSE ||
        is.numeric(data$count) == FALSE) {
    err_msg1 <- "Data should be numeric."
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase1.html#prepare-the-data")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }
}

#' @noRd
# Helper function to validate conditions
validate_condition <- function(condition, message) {
  if (condition) stop(message)
}

#' @noRd
# input checks for all variations
all_checks_1 <- function(data, cov, family, p10_priors, phi_priors, n_chain,
                         n_warmup, n_iter, thin, adapt_delta, seed) {


  ## make sure count, pcr_n, and pcr_k are not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for
  #'   zero-length data
  err_msg1 <- "Input data contains zero-length data."
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(dim(data$pcr_n)[1] == 0 || dim(data$pcr_k)[1] == 0 ||
                       dim(data$count)[1] == 0,
                     err_msg)

  ## make sure no column is entirely NA in pcr_n
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column
  #'   with all NA
  err_msg1 <- "pcr_n contains a column with all NA."
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(any(apply(data$pcr_n, 2, function(col) all(is.na(col)))),
                     err_msg)

  ## make sure no column is entirely NA in pcr_k
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with
  #'   all NA
  err_msg1 <- "pcr_k contains a column with all NA."
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(any(apply(data$pcr_k, 2, function(col) all(is.na(col)))),
                     err_msg)

  ## make sure no column is entirely NA in count
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with
  #'   all NA
  err_msg1 <- "count contains a column with all NA."
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(any(apply(data$count, 2, function(col) all(is.na(col)))),
                     err_msg)

  ## make sure dimensions of pcr_n and pcr_k are equal
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #'   is dimensionally commensurate
  err_msg1 <- "Dimensions of pcr_n and pcr_k do not match."
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(dim(data$pcr_n)[1] != dim(data$pcr_k)[1] ||
                       dim(data$pcr_n)[2] != dim(data$pcr_k)[2],
                     err_msg)

  ## make sure number of rows in count = number of rows in pcr_n and pcr_k
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #'   is dimensionally commensurate
  err_msg1 <- paste0("Number of sites (rows) in pcr data and traditional ",
                     "survey count data do not match.")
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(dim(data$pcr_n)[1] != dim(data$count)[1],
                     err_msg)

  ## make sure locations of NAs in pcr_n data match locations of NAs in
  ## pcr_k data
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #'   is dimensionally commensurate
  err_msg1 <- paste0("Empty data cells (NA) in pcr_n data should match ",
                     "empty data cells (NA) in pcr_k data.")
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(
    any((which(is.na(data$pcr_n)) == which(is.na(data$pcr_k))) == FALSE),
    err_msg
  )

  ## make sure family is either 'poisson', 'negbin', or 'gamma'
  #' @srrstats {G2.3,G2.3a,G2.3b} Permit only expected univariate
  #'   (case-insensitive) parameter values
  validate_condition(!c(tolower(family) %in% c("poisson", "negbin", "gamma")),
                     paste0("Invalid family. Options include 'poisson', ",
                            "'negbin', and 'gamma'."))

  ## p10_priors is a vector of two integers
  #' @srrstats {G2.0,BS2.2,BS2.3,BS2.4,BS2.5,BS2.6} Checks of vector length and
  #'   appropriateness of distributional parameters (i.e., vector of length 2,
  #'   numeric values > 0), implemented prior to analytic routines
  validate_condition(!is.numeric(p10_priors) || length(p10_priors) != 2 ||
                       any(p10_priors <= 0),
                     paste0("p10_priors should be a vector of two positive ",
                            "numeric values. ex. c(1,20)"))

  ## phi_priors is a vector of two numeric values
  #' @srrstats {G2.0,BS2.2,BS2.3,BS2.4,BS2.5,BS2.6} Checks of vector length
  #'   and appropriateness of distributional parameters (i.e., vector of length
  #'   2, numeric values > 0), implemented prior to analytic routines
  if (family == "negbin") {
    validate_condition(!is.numeric(phi_priors) || length(phi_priors) != 2 ||
                         any(phi_priors <= 0),
                       paste0("phi_priors should be a vector of two positive ",
                              "numeric values. ex. c(0.25,0.25)"))
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  err_msg1 <- "count contains undefined values (i.e., Inf or -Inf)"
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(any(data$count == Inf, na.rm =  TRUE) ||
                       any(data$count == -Inf, na.rm =  TRUE),
                     err_msg)
}

all_checks_2 <- function(data, cov, family, p10_priors, phi_priors, n_chain,
                         n_warmup, n_iter, thin, adapt_delta, seed) {

  ## count are integers, if family is poisson or negbin
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #'   for distributional parameters (i.e., count data must be an integer if a
  #'   poisson or negative binomial distribution is used), implemented prior to
  #'   analytic routines
  if (tolower(family) %in% c("poisson", "negbin")) {
    validate_condition(!all(data$count %% 1 %in% c(0, NA)) ||
                         any(data$count < 0, na.rm =  TRUE),
                       paste0("All values in count should be non-negative ",
                              "integers. ",
                              "Use family = 'gamma' if count is continuous."))
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  err_msg1 <- "pcr_n contains undefined values (i.e., Inf or -Inf)"
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(any(data$pcr_n == Inf, na.rm = TRUE) ||
                       any(data$pcr_n == -Inf, na.rm =  TRUE),
                     err_msg)

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  err_msg1 <- "pcr_k contains undefined values (i.e., Inf or -Inf)"
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(any(data$pcr_k == Inf, na.rm =  TRUE) ||
                       any(data$pcr_k == -Inf, na.rm =  TRUE),
                     err_msg)

  ## pcr_n are integers
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #'   for distributional parameters (i.e., pcr data are non-negative
  #'   integers), implemented prior to analytic routines
  err_msg1 <- "All values in pcr_n should be non-negative integers."
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(!all(data$pcr_n %% 1 %in% c(0, NA)) || any(data$pcr_n < 0,
                                                                na.rm =  TRUE),
                     err_msg)

  ## pcr_k are integers
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #'   for distributional parameters (i.e., pcr data are non-negative
  #'   integers), implemented prior to analytic routines
  err_msg1 <- "All values in pcr_k should be non-negative integers."
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(!all(data$pcr_k %% 1 %in% c(0, NA)) || any(data$pcr_k < 0,
                                                                na.rm =  TRUE),
                     err_msg)

  ## check length and range of n_chain
  validate_condition(any(length(as.integer(n_chain)) > 1 || n_chain < 1),
                     "n_chain should be an integer > 0 and of length 1.")

  ## check length and range of n_iter
  validate_condition(any(length(as.integer(n_iter)) > 1 || n_iter < 1),
                     "n_iter should be an integer > 0 and of length 1.")

  ## check length and range of n_warmup
  validate_condition(any(length(as.integer(n_warmup)) > 1 || n_warmup < 1),
                     "n_warmup should be an integer > 0 and of length 1.")

  ## check length and range of thin
  validate_condition(any(length(as.integer(thin)) > 1 || thin < 1),
                     "thin should be an integer > 0 and of length 1.")

  ## check length and range of adapt_delta
  validate_condition(any(length(adapt_delta) > 1 || adapt_delta < 0 ||
                           adapt_delta > 1),
                     paste0("adapt_delta should be a numeric value > 0 and ",
                            "< 1 and of length 1."))

  ## check length of seed
  if (!is.null(seed)) {
    validate_condition(length(as.integer(seed)) > 1,
                       "seed should be an integer of length 1.")
  }

  ## check that N >= K
  err_msg1 <- paste0("N should be >= K in pcr data. N is the number of pcr ",
                     "replicates per sample, and K is the number of ",
                     "positive detections among replicates.")
  err_msg2 <- "See the eDNAjoint guide for data formatting help: "
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#prepare-the-data")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(any(data$pcr_k > data$pcr_n, na.rm =  TRUE),
                     err_msg)
}

#' @noRd
# input checks if site-level covariates are used
covariate_checks <- function(data, cov) {

  ## make sure site_cov is not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for zero-length
  #'   data
  if (dim(data$site_cov)[1] == 0) {
    err_msg1 <- "site_cov contains zero-length data."
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase2.html")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }

  ## check for NA in site_cov
  if (any(is.na(data$site_cov))) {
    err_msg1 <- "site_cov should not contain missing values (i.e., NA)."
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase2.html")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }

  ## site_cov is numeric, if present
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if (!is.numeric(data$site_cov)) {
    err_msg1 <- "site_cov should be numeric."
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase2.html")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  if (any(data$site_cov == Inf) || any(data$site_cov == -Inf)) {
    err_msg1 <- "site_cov contains undefined values (i.e., Inf or -Inf)"
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase2.html")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }

  ## cov values match column names in site_cov
  if (!all(cov %in% colnames(data$site_cov))) {
    err_msg1 <- paste0("cov values should be listed in the column names of ",
                       "site_cov in the data.")
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase2.html")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }

  ## site_cov has same number of rows as pcr_n and count, if present
  #' @srrstats {BS2.1} Pre-processing routines to ensure all input data is
  #'   dimensionally commensurate
  if (dim(data$pcr_n)[1] != dim(data$site_cov)[1]) {
    err_msg1 <- paste0("The number of rows in site_cov matrix should match ",
                       "the number of rows in all other matrices.")
    err_msg2 <- "See the eDNAjoint guide for data formatting help: "
    err_msg3 <- paste0("https://ednajoint.netlify.app",
                       "/usecase2.html")
    err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
    stop(err_msg)
  }

  ## add warning if number of covariates is greater than the number of sites
  #' @srrstats {G5.8d} Pre-processing routines to check if data is outside
  #'   scope of algorithm (i.e., # site-level covariates is greater than the
  #'   number of sites)
  if (length(cov) > dim(data$site_cov)[2]) {
    warn_msg <- paste0("The number of site-level covariates exceeds ",
                       "the number of sites (i.e., n < p).")
    warning(warn_msg)
  }

  ## add warning if number of site-covariate data has perfect collinearity
  #' @srrstats {BS3.1} Pre-processing routines to check if site covariate
  #'   data has perfect collinearity
  rank_mat <- qr(data$site_cov)$rank
  if (rank_mat < ncol(data$site_cov)) {
    warn_msg <- "Data in site_cov exhibits perfect collinearity."
    warning(warn_msg)
  }
}

#' @noRd
# checks if initial values are provided
initial_values_checks_1 <- function(initial_values, data, cov, n_chain) {

  ## length of initial values is equal to the number of chains
  err_msg1 <- paste0("The length of the list of initial values should equal ",
                     "the number of chains (n_chain, default is 4).")
  err_msg2 <- paste0("See the eDNAjoint guide for help formatting ",
                     "initial values: ")
  err_msg3 <- paste0("https://ednajoint.netlify.app",
                     "/usecase1.html#initialvalues")
  err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
  validate_condition(!is.null(initial_values) &&
                       length(initial_values) != n_chain,
                     err_msg)
}
initial_values_checks_2 <- function(initial_values, data, cov, n_chain) {

  for (i in 1:n_chain) {

    ## check mu input
    if ("mu" %in% names(initial_values[[i]])) {
      ## if mu is numeric
      err_msg1 <- "Initial values for 'mu' should be numeric values > 0."
      err_msg2 <- paste0("See the eDNAjoint guide for help formatting ",
                         "initial values: ")
      err_msg3 <- paste0("https://ednajoint.netlify.app",
                         "/usecase1.html#initialvalues")
      err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
      validate_condition(any(!is.numeric(initial_values[[i]]$mu)) ||
                           any(initial_values[[i]]$mu < 0),
                         err_msg)

      ## check mu length
      err_msg1 <- paste0("The length of initial values for 'mu' should ",
                         "equal the number of sites.")
      err_msg2 <- paste0("See the eDNAjoint guide for help formatting ",
                         "initial values: ")
      err_msg3 <- paste0("https://ednajoint.netlify.app",
                         "/usecase1.html#initialvalues")
      err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
      validate_condition(length(initial_values[[i]]$mu) != dim(data$count)[1],
                         err_msg)
    }

    ## check p10 input
    if ("p10" %in% names(initial_values[[i]])) {
      ## if p10 is numeric
      err_msg1 <- "Initial values for 'p10' should be numeric."
      err_msg2 <- paste0("See the eDNAjoint guide for help formatting ",
                         "initial values: ")
      err_msg3 <- paste0("https://ednajoint.netlify.app",
                         "/usecase1.html#initialvalues")
      err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
      validate_condition(!is.numeric(initial_values[[i]]$p10),
                         err_msg)

      ## check p10 length
      err_msg1 <- "The length of initial values for 'p10' should equal 1."
      err_msg2 <- paste0("See the eDNAjoint guide for help formatting ",
                         "initial values: ")
      err_msg3 <- paste0("https://ednajoint.netlify.app",
                         "/usecase1.html#initialvalues")
      err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
      validate_condition(length(initial_values[[i]]$p10) != 1,
                         err_msg)
    }

    ## check alpha input -- no covariates
    if ("alpha" %in% names(initial_values[[i]]) && is.null(cov)) {
      ## if alpha is numeric
      err_msg1 <- "Initial values for 'alpha' should be numeric."
      err_msg2 <- paste0("See the eDNAjoint guide for help formatting ",
                         "initial values: ")
      err_msg3 <- paste0("https://ednajoint.netlify.app",
                         "/usecase1.html#initialvalues")
      err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
      validate_condition(!is.numeric(initial_values[[i]]$alpha),
                         err_msg)

      ## check alpha length
      err_msg1 <- "The length of initial values for 'alpha' should equal 1."
      err_msg2 <- paste0("See the eDNAjoint guide for help formatting ",
                         "initial values: ")
      err_msg3 <- paste0("https://ednajoint.netlify.app",
                         "/usecase1.html#initialvalues")
      err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
      validate_condition(length(initial_values[[i]]$alpha) != 1,
                         err_msg)
    }

    ## check alpha input -- covariates
    if ("alpha" %in% names(initial_values[[i]]) && !is.null(cov)) {
      ## if alpha is numeric
      err_msg1 <- "Initial values for 'alpha' should be numeric."
      err_msg2 <- paste0("See the eDNAjoint guide for help formatting ",
                         "initial values: ")
      err_msg3 <- paste0("https://ednajoint.netlify.app",
                         "/usecase2.html#initialvalues")
      err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
      validate_condition(any(!is.numeric(initial_values[[i]]$alpha)),
                         err_msg)
      ## check alpha length
      err_msg1 <- paste0("The length of initial values for 'alpha' should ",
                         "equal: # covariates + 1 (i.e., including ",
                         "intercept).")
      err_msg2 <- paste0("See the eDNAjoint guide for help formatting ",
                         "initial values: ")
      err_msg3 <- paste0("https://ednajoint.netlify.app",
                         "/usecase2.html#initialvalues")
      err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
      validate_condition(length(initial_values[[i]]$alpha) != (length(cov) + 1),
                         err_msg)
    }

    ## check q input
    if ("q" %in% names(initial_values[[i]])) {
      ## if q is numeric
      err_msg1 <- "Initial values for 'q' should be numeric."
      err_msg2 <- paste0("See the eDNAjoint guide for help formatting ",
                         "initial values: ")
      err_msg3 <- paste0("https://ednajoint.netlify.app",
                         "/usecase2.html#initialvalues")
      err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
      validate_condition(any(!is.numeric(initial_values[[i]]$q)) ||
                           any(initial_values[[i]]$q < 0),
                         err_msg)

      ## check q length
      err_msg1 <- paste0("The length of initial values for 'q' should equal:",
                         " # unique gear types - 1 (i.e., q for reference ",
                         "type = 1).")
      err_msg2 <- paste0("See the eDNAjoint guide for help formatting ",
                         "initial values: ")
      err_msg3 <- paste0("https://ednajoint.netlify.app",
                         "/usecase2.html#initialvalues")
      err_msg <- paste(err_msg1, err_msg2, err_msg3, sep = "\n")
      validate_condition(length(initial_values[[i]]$q) !=
                           (length(table(data$count_type)) - 1),
                         err_msg)
    }
  }
}
