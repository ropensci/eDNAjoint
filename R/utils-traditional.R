## helper functions for traditionalModel.R ##

###########
#helper functions: initial values
###########
#' @noRd
#' @srrstats {BS2.7,BS2.11} Option for user to provide initial values for
#'   each chain
init_trad_catchability <- function(n_chain, count_all, q_names,
                                   initial_values) {
  #helper function
  #traditional model, catchability coefficient
  init_list <- list()
  if (all(!is.null(initial_values))) {
    for (i in 1:n_chain) {
      init_list[[i]] <- list(
        mu = if ("mu" %in% names(initial_values[[i]])) {
          initial_values[[i]]$mu
        } else {
          stats::runif(length(unique(count_all$L_ind)), 0.01, 5)
        },
        q = if ("q" %in% names(initial_values[[i]])) {
          as.data.frame(initial_values[[i]]$q)
        } else {
          as.data.frame(stats::runif(length(q_names), 0.01, 1))
        }
      )
    }
  } else {
    for (i in 1:n_chain) {
      init_list[[i]] <- list(
        mu = stats::runif(length(unique(count_all$L_ind)), 0.01, 5),
        q = as.data.frame(stats::runif(length(q_names), 0.01, 1))
      )
    }
  }

  return(init_list)
}

#' @noRd
init_trad <- function(n_chain, count_all, initial_values) {
  #helper function
  #traditional model
  init_list <- list()
  if (all(!is.null(initial_values))) {
    for (i in 1:n_chain) {
      init_list[[i]] <- list(
        mu = if ("mu" %in% names(initial_values[[i]])) {
          initial_values[[i]]$mu
        } else {
          stats::runif(length(unique(count_all$L_ind)), 0.01, 5)
        }
      )
    }
  } else {
    for (i in 1:n_chain) {
      init_list[[i]] <- list(
        mu = stats::runif(length(unique(count_all$L_ind)), 0.01, 5)
      )
    }
  }

  return(init_list)

}

#' @noRd
# Helper function to validate conditions
validate_condition <- function(condition, message) {
  if (condition) stop(message)
}

# function for input checks
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#'   messages
#' @noRd
trad_model_input_checks_1 <- function(data, family, q, phi_priors, n_chain,
                                      n_warmup, n_iter,
                                      thin, adapt_delta, seed) {

  ## make sure all data tags are valid -- if q == TRUE
  #' @srrstats {G2.13} Pre-processing routines to check for missing data
  validate_condition(q == TRUE &&
                       !all(c("count_type", "count") %in% names(data)),
                     "Data should include 'count' and 'count_type'.")

  ## make sure all data tags are valid -- if q == FALSE
  #' @srrstats {G2.13} Pre-processing routines to check for missing data
  validate_condition(q == FALSE && !all(c("count") %in% names(data)),
                     "Data should include 'count'.")

  ## make sure count is not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for
  #'   zero-length data
  validate_condition(dim(data$count)[1] == 0,
                     "count contains zero-length data.")

  ## make sure no column is entirely NA in count
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column
  #'   with all NA
  validate_condition(any(apply(data$count, 2, function(col) all(is.na(col)))),
                     "count contains a column with all NA.")

  ## make sure dimensions of count and count_type are equal, if count_type is
  ## present
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #'   is dimensionally commensurate
  if (q == TRUE) {
    validate_condition(dim(data$count)[1] != dim(data$count_type)[1] ||
                         dim(data$count)[2] != dim(data$count_type)[2],
                       "Dimensions of count and count_type do not match.")
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  validate_condition(any(data$count == Inf, na.rm = TRUE) ||
                       any(data$count == -Inf, na.rm = TRUE),
                     "count contains undefined values (i.e., Inf or -Inf)")

  if (q == TRUE) {
    ## make sure all data is numeric -- if q == TRUE
    #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
    #'   for distributional parameters (i.e., count data must numeric),
    #'   implemented prior to analytic routines
    #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
    #'   unsupported type
    validate_condition(is.numeric(data$count) == FALSE ||
                         is.numeric(data$count_type) == FALSE,
                       "Data should be numeric.")

    ## make sure locations of NAs in count data match locations of NAs in
    ## count_type data
    #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input
    #'   data is dimensionally commensurate
    validate_condition(
      any((which(is.na(data$count_type)) == which(is.na(data$count))) == FALSE),
      paste0("Empty data cells (NA) in count data should match ",
             "empty data cells (NA) in count_type data.")
    )

    ## make sure count_type is not zero-length
    #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for
    #'   zero-length data
    validate_condition(dim(data$count_type)[1] == 0,
                       "count_type contains zero-length data.")

    ## make sure no column is entirely NA in count_type
    #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column
    #'   with all NA
    validate_condition(any(apply(data$count_type, 2,
                                 function(col) all(is.na(col)))),
                       "count_type contains a column with all NA.")
  }

  ## make sure all data is numeric -- if q == FALSE
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #'   for distributional parameters (i.e., count data must positive and
  #'   numeric), implemented prior to analytic routines
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if (q == FALSE) {
    validate_condition(is.numeric(data$count) == FALSE ||
                         any(data$count < 0, na.rm = TRUE),
                       "Data should be numeric.")
  }

  ## make sure family is either 'poisson', 'negbin', or 'gamma'
  #' @srrstats {G2.3,G2.3a,G2.3b} Permit only expected univariate
  #'   (case-insensitive) parameter values
  validate_condition(!c(tolower(family) %in% c("poisson", "negbin", "gamma")),
                     paste0("Invalid family. Options include 'poisson', ",
                            "'negbin', or 'gamma'."))

  ## the smallest count_type is 1
  validate_condition(q == TRUE && min(data$count_type, na.rm = TRUE) != 1,
                     paste0("The first gear type should be referenced as 1 in ",
                            "count_type. Subsequent gear types should be ",
                            "referenced 2, 3, 4, etc."))

  ## count are integers, if family is poisson or negbin
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #'   for distributional parameters (i.e., count data must be non-negative
  #'   integers if a poisson or negative binomial distribution is used),
  #'   implemented prior to analytic routines
  if (tolower(family) %in% c("poisson", "negbin")) {
    validate_condition(!all(data$count %% 1 %in% c(0, NA)) ||
                         any(data$count < 0, na.rm = TRUE),
                       paste0("All values in count should be non-negative ",
                              "integers. Use family = 'gamma' if count is ",
                              "continuous."))
  }

}

trad_model_input_checks_2 <- function(data, family, q, phi_priors, n_chain,
                                      n_warmup, n_iter,
                                      thin, adapt_delta, seed) {

  ## count_type are integers
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  validate_condition(q == TRUE && !all(data$count_type %% 1 %in% c(0, NA)),
                     "All values in count_type should be integers.")

  ## phi_priors is a vector of two numeric values
  #' @srrstats {G2.0,BS2.2,BS2.3,BS2.4,BS2.5} Checks of vector length and
  #'   appropriateness of distributional parameters (i.e., vector of length 2,
  #'   numeric values > 0), implemented prior to analytic routines
  if (family == "negbin") {
    validate_condition(!is.numeric(phi_priors) || length(phi_priors) != 2 ||
                         any(phi_priors <= 0),
                       paste0("phi_priors should be a vector of two positive ",
                              "numeric values. ex. c(0.25,0.25)"))
  }

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
}

# checks if initial values are provided
#' @noRd
initial_values_checks_trad <- function(initial_values, data, n_chain) {

  ## length of initial values is equal to the number of chains
  validate_condition(length(initial_values) != n_chain,
                     paste0("The length of the list of initial values should ",
                            "equal the number of chains (n_chain, default ",
                            "is 4)."))

  for (i in 1:n_chain) {

    ## check mu input
    if ("mu" %in% names(initial_values[[i]])) {
      ## if mu is numeric
      validate_condition(any(!is.numeric(initial_values[[i]]$mu)) ||
                           any(initial_values[[i]]$mu < 0),
                         paste0("Initial values for 'mu' should be numeric ",
                                "values > 0."))
      ## check mu length
      validate_condition(length(initial_values[[i]]$mu) != dim(data$count)[1],
                         paste0("The length of initial values for 'mu' should ",
                                "equal the number of sites."))
    }

    ## check q input
    if ("q" %in% names(initial_values[[i]])) {
      ## if q is numeric
      validate_condition(any(!is.numeric(initial_values[[i]]$q)) ||
                           any(initial_values[[i]]$q < 0),
                         "Initial values for 'q' should be numeric.")
      ## check q length
      validate_condition(length(initial_values[[i]]$q) !=
                           (length(table(data$count_type)) - 1),
                         paste0("The length of initial values for 'q' should ",
                                "equal: # unique gear types - 1 (i.e., q for ",
                                "reference type = 1)."))
    }
  }
}
