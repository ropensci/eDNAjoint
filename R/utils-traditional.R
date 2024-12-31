## helper functions for traditionalModel.R ##

###########
#helper functions: initial values
###########
#' @noRd
#' @srrstats {BS2.7,BS2.11} Option for user to provide initial values for
#'   each chain
init_trad_catchability <- function(n.chain, count_all, q_names, initial_values){
  #helper function
  #traditional model, catchability coefficient
  A <- list()
  if(all(!is.null(initial_values))){
    for(i in 1:n.chain){
      A[[i]] <- list(
        if('mu' %in% names(initial_values[[i]])){
          mu <- initial_values[[i]]$mu
        } else {
          mu <- stats::runif(length(unique(count_all$L_ind)), 0.01, 5)
        },

        if('q' %in% names(initial_values[[i]])){
          q <- as.data.frame(initial_values[[i]]$q)
        } else {
          q <- as.data.frame(stats::runif(length(q_names),0.01,1))
        }
      )
      names(A[[i]]) <- c('mu','q')
    }
  } else {
    for(i in 1:n.chain){
      A[[i]] <- list(
        mu <- stats::runif(length(unique(count_all$L_ind)), 0.01, 5),
        q <- as.data.frame(stats::runif(length(q_names),0.01,1))
      )
      names(A[[i]]) <- c('mu','q')
    }
  }

  return(A)
}

#' @noRd
init_trad <- function(n.chain, count_all, initial_values){
  #helper function
  #traditional model
  A <- list()
  if(all(!is.null(initial_values))){
    for(i in 1:n.chain){
      A[[i]] <- list(
        if('mu' %in% names(initial_values[[i]])){
          mu <- initial_values[[i]]$mu
        } else {
          mu <- stats::runif(length(unique(count_all$L_ind)), 0.01, 5)
        }
      )
      names(A[[i]]) <- 'mu'
    }
  } else {
    for(i in 1:n.chain){
      A[[i]] <- list(
        mu <- stats::runif(length(unique(count_all$L_ind)), 0.01, 5)
      )
      names(A[[i]]) <- 'mu'
    }
  }

  return(A)
}

# function for input checks
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#'   messages
#' @noRd
traditionalModel_input_checks <- function(data, family, q, phipriors, n.chain,
                                          n.warmup, n.iter,
                                          thin, adapt_delta, seed){

  ## make sure all data tags are valid -- if q == TRUE
  #' @srrstats {G2.13} Pre-processing routines to check for missing data
  if (q == TRUE && !all(c('count.type','count') %in% names(data))){
    errMsg <- "Data should include 'count' and 'count.type'."
    stop(errMsg)
  }

  ## make sure all data tags are valid -- if q == FALSE
  #' @srrstats {G2.13} Pre-processing routines to check for missing data
  if (q == FALSE && !all(c('count') %in% names(data))){
    errMsg <- "Data should include 'count'."
    stop(errMsg)
  }

  ## make sure count is not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for
  #'   zero-length data
  if (dim(data$count)[1] == 0) {
    errMsg <- "count contains zero-length data."
    stop(errMsg)
  }
  ## make sure no column is entirely NA in count
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column
  #'   with all NA
  if (any(apply(data$count, 2, function(col) all(is.na(col))))) {
    errMsg <- "count contains a column with all NA."
    stop(errMsg)
  }

  ## make sure dimensions of count and count.type are equal, if count.type is
  ## present
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #'   is dimensionally commensurate
  if (q == TRUE){
    if(dim(data$count)[1] != dim(data$count.type)[1]|
       dim(data$count)[2] != dim(data$count.type)[2]) {
      errMsg <- "Dimensions of count and count.type do not match."
      stop(errMsg)
    }
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  if(any(data$count == Inf,na.rm=TRUE) | any(data$count == -Inf,na.rm=TRUE)){
    errMsg <- "count contains undefined values (i.e., Inf or -Inf)"
    stop(errMsg)
  }

  ## make sure all data is numeric -- if q == TRUE
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #'   for distributional parameters (i.e., count data must numeric),
  #'   implemented prior to analytic routines
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if (q == TRUE) {
    if(is.numeric(data$count) == FALSE |
       is.numeric(data$count.type) == FALSE) {
      errMsg <- "Data should be numeric."
      stop(errMsg)
    }
  }

  ## make sure all data is numeric -- if q == FALSE
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #'   for distributional parameters (i.e., count data must positive and
  #'   numeric), implemented prior to analytic routines
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if (q == FALSE) {
    if(is.numeric(data$count) == FALSE | any(data$count < 0, na.rm=TRUE)) {
      errMsg <- "Data should be numeric."
      stop(errMsg)
    }
  }

  if(q == TRUE){
    ## make sure locations of NAs in count data match locations of NAs in
    ## count.type data
    #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input
    #'   data is dimensionally commensurate
    if(any((which(is.na(data$count.type)) == which(is.na(data$count))) == FALSE)
    ){
      errMsg <- paste0("Empty data cells (NA) in count data should match ",
                       "empty data cells (NA) in count.type data.")
      stop(errMsg)
    }
    ## make sure count.type is not zero-length
    #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for
    #'   zero-length data
    if (dim(data$count.type)[1] == 0) {
      errMsg <- "count.type contains zero-length data."
      stop(errMsg)
    }
    ## make sure no column is entirely NA in count.type
    #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column
    #'   with all NA
    if (any(apply(data$count.type, 2, function(col) all(is.na(col))))) {
      errMsg <- "count.type contains a column with all NA."
      stop(errMsg)
    }
  }

  ## make sure family is either 'poisson', 'negbin', or 'gamma'
  #' @srrstats {G2.3,G2.3a,G2.3b} Permit only expected univariate
  #'   (case-insensitive) parameter values
  if(!c(tolower(family) %in% c('poisson','negbin','gamma'))){
    errMsg <- "Invalid family. Options include 'poisson', 'negbin', or 'gamma'."
    stop(errMsg)
  }

  ## the smallest count.type is 1
  if(q == TRUE && min(data$count.type,na.rm=TRUE) != 1){
    errMsg <- paste0("The first gear type should be referenced as 1 in ",
                     "count.type. Subsequent gear types should be ",
                     "referenced 2, 3, 4, etc.")
    stop(errMsg)
  }

  ## count are integers, if family is poisson or negbin
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #'   for distributional parameters (i.e., count data must be non-negative
  #'   integers if a poisson or negative binomial distribution is used),
  #'   implemented prior to analytic routines
  if(tolower(family) %in% c('poisson','negbin')){
    if(!all(data$count %% 1 %in% c(0,NA)) | any(data$count < 0,na.rm=TRUE)){
      errMsg <- paste0("All values in count should be non-negative integers. ",
                       "Use family = 'gamma' if count is continuous.")
      stop(errMsg)
    }
  }

  ## count.type are integers
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(q == TRUE && !all(data$count.type %% 1 %in% c(0,NA))){
    errMsg <- "All values in count.type should be integers."
    stop(errMsg)
  }

  ## phipriors is a vector of two numeric values
  #' @srrstats {G2.0,BS2.2,BS2.3,BS2.4,BS2.5} Checks of vector length and
  #'   appropriateness of distributional parameters (i.e., vector of length 2,
  #'   numeric values > 0), implemented prior to analytic routines
  if(family == 'negbin'){
    if(!is.numeric(phipriors) | length(phipriors) != 2 | any(phipriors <= 0)){
      errMsg <- paste0("phipriors should be a vector of two positive ",
                       "numeric values. ex. c(0.25,0.25)")
      stop(errMsg)
    }
  }

  ## check length and range of n.chain
  if(any(length(as.integer(n.chain)) > 1 | n.chain < 1)){
    errMsg <- "n.chain should be an integer > 0 and of length 1."
    stop(errMsg)
  }

  ## check length and range of n.iter
  if(any(length(as.integer(n.iter)) > 1 | n.iter < 1)){
    errMsg <- "n.iter should be an integer > 0 and of length 1."
    stop(errMsg)
  }

  ## check length and range of n.warmup
  if(any(length(as.integer(n.warmup)) > 1 | n.warmup < 1)){
    errMsg <- "n.warmup should be an integer > 0 and of length 1."
    stop(errMsg)
  }

  ## check length and range of thin
  if(any(length(as.integer(thin)) > 1 | thin < 1)){
    errMsg <- "thin should be an integer > 0 and of length 1."
    stop(errMsg)
  }

  ## check length and range of adapt_delta
  if(any(length(adapt_delta) > 1 | adapt_delta < 0 | adapt_delta > 1)){
    errMsg <- paste0("adapt_delta should be a numeric value > 0 and < 1 and ",
                     "of length 1.")
    stop(errMsg)
  }

  ## check length of seed
  if(!is.null(seed)){
    if(length(as.integer(seed)) > 1){
      errMsg <- "seed should be an integer of length 1."
      stop(errMsg)
    }
  }
}

# checks if initial values are provided
#' @noRd
initial_values_checks_trad <- function(initial_values,data,n.chain){

  ## length of initial values is equal to the number of chains
  if(length(initial_values) != n.chain){
    errMsg <- paste0("The length of the list of initial values should equal ",
                     "the number of chains (n.chain, default is 4).")
    stop(errMsg)
  }

  for(i in 1:n.chain){

    ## check mu input
    if('mu' %in% names(initial_values[[i]])){
      ## if mu is numeric
      if(any(!is.numeric(initial_values[[i]]$mu)) |
         any(initial_values[[i]]$mu < 0)){
        errMsg <- "Initial values for 'mu' should be numeric values > 0."
        stop(errMsg)
      }
      ## check mu length
      if(length(initial_values[[i]]$mu) != dim(data$count)[1]){
        errMsg <- paste0("The length of initial values for 'mu' should ",
                         "equal the number of sites.")
        stop(errMsg)
      }
    }

    ## check q input
    if('q' %in% names(initial_values[[i]])){
      ## if q is numeric
      if(any(!is.numeric(initial_values[[i]]$q)) |
         any(initial_values[[i]]$q < 0)){
        errMsg <- "Initial values for 'q' should be numeric."
        stop(errMsg)
      }
      ## check q length
      if(length(initial_values[[i]]$q) != (length(table(data$count.type))-1)){
        errMsg <- paste0("The length of initial values for 'q' should equal: ",
                         "# unique gear types - 1 (i.e., q for reference ",
                         "type = 1).")
        stop(errMsg)
      }
    }


  }
}



