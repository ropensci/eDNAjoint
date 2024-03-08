#' Specify and fit model using count data from traditional, non eDNA surveys
#'
#' This function implements a Bayesian model that estimates expected species catch rate using count data from traditional, non eDNA surveys. When multiple traditional gear types are used, an optional variation allows estimation of catchability coefficients, which scale the catchability of gear types relative to the expected catch rate of a reference gear type. Model is implemented using Bayesian inference using the `rstan` package, which uses Hamiltonian Monte Carlo to simulate the posterior distributions.
#'
#' @srrstats {G1.4} Roxygen function documentation begins here
#' @export
#' @srrstats {BS1.1,BS3.0} Descriptions of how to enter data, description of how NAs are handled
#' @param data A list containing data necessary for model fitting. Valid tags are `count` and `count.type`. `count` is a matrix or data frame of traditional survey count data, with first dimension equal to the number of sites (i) and second dimension equal to the maximum number of traditional survey replicates at a given site (j). `count.type` is an optional matrix or data frame of integers indicating gear type (k) used in corresponding count data, with first dimension equal to the number of sites (i) and second dimension equal to the maximum number of traditional survey replicates at a given site (j). Values should be integers beginning with 1 (referring to the first gear type) to n (last gear type). Empty cells should be NA and will be removed during processing. Sites, i, should be consistent in all count data.
#' @param family The distribution class used to model traditional survey count data. Options include poisson ('poisson'), negative binomial ('negbin'), and gamma ('gamma'). Default value is 'poisson'.
#' @param q A logical value indicating whether to estimate a catchability coefficient, q, for traditional survey gear types (TRUE) or to not estimate a catchability coefficient, q, for traditional survey gear types (FALSE). Default value is FALSE.
#' @srrstats {BS1.0,G2.1a,BS1.2} Description of hyperparameters and how to specify prior distributions, explicit documentation of vector input types
#' @param phipriors A numeric vector indicating gamma distribution hyperparameters (shape, rate) used as the prior distribution for phi, the overdispersion in the negative binomial distribution for traditional survey gear data. Used when family = 'negbin.' Default vector is c(0.25,0.25).
#' @param multicore A logical value indicating whether to parallelize chains with multiple cores. Default is TRUE.
#' @srrstats {BS1.3} Description of parameters used in the computational process begins here
#' @param n.chain Number of MCMC chains. Default value is 4.
#' @param n.iter.burn Number of warm-up MCMC iterations. Default value is 500.
#' @param n.iter.sample Number of sampling MCMC iterations. Default value is 2500.
#' @param thin A positive integer specifying the period for saving samples. Default value is 1.
#' @param adapt_delta Target average acceptance probability used in `rstan::sampling`. Default value is 0.9.
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @note  Before fitting the model, this function checks to ensure that the model specification is possible given the data files. These checks include:
#' \itemize{
#' \item  All tags in data are valid (i.e., include count and count.type).
#' \item  Number of sites in count and count type data are equal.
#' \item  All data are numeric (i.e., integer or NA).
#' \item  Empty data cells (NA) match in count and count.type.
#' @srrstats {G2.3,G2.3a,G2.3b} Permit only expected univariate (case-insensitive) parameter values
#' \item  family is 'poisson', 'negbin', or 'gamma'.
#' @srrstats {G2.0a} Explicit secondary documentation of any expectations on lengths of inputs
#' \item  phipriors (if used) is a vector of two numeric values.
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' @examples
#' \donttest{
#' # Load data
#' @srrstats {BS1.1} Descriptions of how to enter data
#' data(greencrabData)
#'
#' # Examine data in list
#' # This function uses only traditional survey count data and optionally the count type data
#' names(greencrabData)
#'
#' # Note that the surveyed sites (rows) should match in the data
#' dim(greencrabData$count)[1]
#' dim(greencrabData$count.type)[1]
#'
#' # Fit a model without estimating a catchability coefficient for traditional survey gear types.
#' # This model assumes all traditional survey methods have the same catchability.
#' # Count data is modeled using a poisson distribution.
#' fit.no.q = traditionalModel(data=greencrabData, family='poisson', q=FALSE)
#'
#'
#' # Fit a model estimating a catchability coefficient for traditional survey gear types.
#' # This model does not assume all traditional survey methods have the same catchability.
#' # Count data is modeled using a negative binomial distribution.
#' fit.q = traditionalModel(data=greencrabData, family='negbin', q=TRUE)
#' }
#'

traditionalModel <- function(data, family='poisson',
                             q=FALSE, phipriors=c(0.25,0.25), multicore=TRUE,
                             n.chain=4, n.iter.burn=500,
                             n.iter.sample=2500, thin=1,
                             adapt_delta=0.9) {

  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using this helper function
  traditionalModel_input_checks(data, family, q, phipriors)

  if (!requireNamespace("rstan", quietly = TRUE)){
    stop ("The 'rstan' package is not installed.", call. = FALSE)
  }

  # make character inputs case-insensitive
  #' @srrstats {G2.3b} Allow case-insensitive character parameter values
  family <- tolower(family)

  ###
  #convert data to long format
  '%>%' <- magrittr::`%>%`

  #convert count data to long format
  count_all <- as.data.frame(data$count) %>%
    dplyr::mutate(L=1:dim(data$count)[1]) %>%
    tidyr::pivot_longer(cols=!L,values_to='count') %>%
    tidyr::drop_na()

  #if q==TRUE, add count type data to count df
  if(q==TRUE){
    q_ref <- 1
    count.type_df <- as.data.frame(data$count.type) %>%
      dplyr::mutate(L=1:dim(data$count.type)[1]) %>%
      tidyr::pivot_longer(cols=!L,values_to='count.type') %>%
      tidyr::drop_na()
    count_all$count.type <- count.type_df$count.type

  #create vector of q coefficient names
  counttypes <- unique(count_all$count.type)
  names <- counttypes[!counttypes==q_ref]
  q_names <- paste0('q_',names)

  #add dummy variables for count type
  for(i in seq_along(q_names)){
    count_all[,q_names[i]] <- ifelse(count_all$count.type==names[i],1,0)
  }

  }

  # set up core number
  if(multicore == TRUE){
    cores <- parallel::detectCores()
  } else {
    cores <- 1
  }

  ##run model, catchability, pois/gamma
  if(q==TRUE&&family!='negbin'){
    model_index <- dplyr::case_when(family=='poisson'~ 1,
                                    family=='gamma' ~ 2)
    out <- rstan::sampling(c(stanmodels$traditional_catchability_pois,
                             stanmodels$traditional_catchability_gamma)[model_index][[1]],
                           data = list(
                             Nloc = length(unique(count_all$L)),
                             C = nrow(count_all),
                             R = count_all$L,
                             E = count_all$count,
                             nparams = length(q_names),
                             mat = as.matrix(count_all[,q_names]),
                             control = list(adapt_delta = adapt_delta)
                           ),
                           cores = cores,
                           chains = n.chain,
                           thin = thin,
                           warmup = n.iter.burn,
                           iter = n.iter.burn + n.iter.sample,
                           init = init_trad_catchability(n.chain, count_all, q_names)
    )
  } else if(q==TRUE&&family=='negbin'){
    ##run model, catchability, negbin
    out <- rstan::sampling(stanmodels$traditional_catchability_negbin,
                           data = list(
                             Nloc = length(unique(count_all$L)),
                             C = nrow(count_all),
                             R = count_all$L,
                             E = count_all$count,
                             nparams = length(q_names),
                             mat = as.matrix(count_all[,q_names]),
                             phipriors = phipriors,
                             control = list(adapt_delta = adapt_delta)
                           ),
                           cores = cores,
                           chains = n.chain,
                           thin = thin,
                           warmup = n.iter.burn,
                           iter = n.iter.burn + n.iter.sample,
                           init = init_trad_catchability(n.chain, count_all, q_names)
    )
  } else if(q==FALSE&&family!='negbin'){
    ##run model, no catchability, pois/gamma
    model_index <- dplyr::case_when(family=='poisson'~ 1,
                                    family=='gamma' ~ 2)
    out <- rstan::sampling(c(stanmodels$traditional_pois,
                             stanmodels$traditional_gamma)[model_index][[1]],
                           data = list(
                             Nloc = length(unique(count_all$L)),
                             C = nrow(count_all),
                             R = count_all$L,
                             E = count_all$count,
                             control = list(adapt_delta = adapt_delta,
                                            stepsize = 0.5)
                           ),
                           cores = cores,
                           chains = n.chain,
                           thin = thin,
                           warmup = n.iter.burn,
                           iter = n.iter.burn + n.iter.sample,
                           init = init_trad(n.chain, count_all)
    )
  } else if(q==FALSE&&family=='negbin'){
    ##run model, no catchability, negbin
    out <- rstan::sampling(stanmodels$traditional_negbin,
                           data = list(
                             Nloc = length(unique(count_all$L)),
                             C = nrow(count_all),
                             R = count_all$L,
                             E = count_all$count,
                             phipriors = phipriors,
                             control = list(adapt_delta = adapt_delta)
                           ),
                           cores = cores,
                           chains = n.chain,
                           thin = thin,
                           warmup = n.iter.burn,
                           iter = n.iter.burn + n.iter.sample,
                           init = init_trad(n.chain, count_all)
    )
  }

  # assert that the log likelihood is a double
  stopifnot(is.double(sum(colMeans(rstan::extract(out,par='log_lik')$log_lik))))

  return(out)
}

init_trad_catchability <- function(n.chain, count_all, q_names){
  #helper function
  #traditional model, catchability coefficient
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      mu = stats::runif(length(unique(count_all$L)), 0.01, 5),
      q = as.data.frame(stats::runif(length(q_names),0.01,1))
    )
  }
  return(A)
}

init_trad <- function(n.chain, count_all){
  #helper function
  #traditional model
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      mu = stats::runif(length(unique(count_all$L)), 0.01, 5)
    )
  }
  return(A)
}

# function for input checks
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique messages
traditionalModel_input_checks <- function(data, family, q, phipriors){

  ## make sure all data tags are valid -- if q == TRUE
  #' @srrstats {G2.13} Pre-processing routines to check for missing data
  if (q==TRUE && !all(c('count.type','count') %in% names(data))){
    errMsg = paste("Data should include 'count' and 'count.type'.")
    stop(errMsg)
  }

  ## make sure all data tags are valid -- if q == FALSE
  #' @srrstats {G2.13} Pre-processing routines to check for missing data
  if (q==FALSE && !all(c('count') %in% names(data))){
    errMsg = paste("Data should include 'count'.")
    stop(errMsg)
  }

  ## make sure count is not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for zero-length data
  if (dim(data$count)[1]==0) {
    errMsg = paste("count contains zero-length data.")
    stop(errMsg)
  }
  ## make sure no column is entirely NA in count
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with all NA
  if (any(apply(data$count, 2, function(col) all(is.na(col))))) {
    errMsg = paste("count contains a column with all NA.")
    stop(errMsg)
  }

  ## make sure dimensions of count and count.type are equal, if count.type is present
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data is dimensionally commensurate
  if (q==TRUE){
    if(dim(data$count)[1]!=dim(data$count.type)[1]|dim(data$count)[2]!=dim(data$count.type)[2]) {
      errMsg = paste("Dimensions of count and count.type do not match.")
      stop(errMsg)
    }
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  if(any(data$count==Inf) | any(data$count==-Inf)){
    errMsg = paste("count contains undefined values (i.e., Inf or -Inf)")
    stop(errMsg)
  }

  ## make sure all data is numeric -- if q == TRUE
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted for distributional parameters (i.e., count data must numeric), implemented prior to analytic routines
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of unsupported type
  if (q==TRUE) {
    if(is.numeric(data$count)==FALSE |
       is.numeric(data$count.type)==FALSE) {
      errMsg = paste("Data should be numeric.")
      stop(errMsg)
    }
  }

  ## make sure all data is numeric -- if q == FALSE
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted for distributional parameters (i.e., count data must positive and numeric), implemented prior to analytic routines
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of unsupported type
  if (q==FALSE) {
    if(is.numeric(data$count)==FALSE | any(data$count < 0)) {
      errMsg = paste("Data should be numeric.")
      stop(errMsg)
    }
  }

  if(q==TRUE){
    ## make sure locations of NAs in count data match locations of NAs in count.type data
    #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data is dimensionally commensurate
    if(any((which(is.na(data$count.type))==which(is.na(data$count)))==FALSE)){
      errMsg = paste("Empty data cells (NA) in count data should match empty data cells (NA) in count.type data.")
      stop(errMsg)
    }
    ## make sure count.type is not zero-length
    #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for zero-length data
    if (dim(data$count.type)[1]==0) {
      errMsg = paste("count.type contains zero-length data.")
      stop(errMsg)
    }
    ## make sure no column is entirely NA in count.type
    #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with all NA
    if (any(apply(data$count.type, 2, function(col) all(is.na(col))))) {
      errMsg = paste("count.type contains a column with all NA.")
      stop(errMsg)
    }
  }

  ## make sure family is either 'poisson', 'negbin', or 'gamma'
  #' @srrstats {G2.3,G2.3a,G2.3b} Permit only expected univariate (case-insensitive) parameter values
  if(!c(tolower(family) %in% c('poisson','negbin','gamma'))){
    errMsg = paste("Invalid family. Options include 'poisson', 'negbin', or 'gamma'.")
    stop(errMsg)
  }

  ## the smallest count.type is 1
  if(q==TRUE && min(data$count.type,na.rm=TRUE) != 1){
    errMsg = paste("The first gear type should be referenced as 1 in count.type. Subsequent gear types should be referenced 2, 3, 4, etc.")
    stop(errMsg)
  }

  ## count are integers, if family is poisson or negbin
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted for distributional parameters (i.e., count data must be non-negative integers if a poisson or negative binomial distribution is used), implemented prior to analytic routines
  if(!all(data$count %% 1 %in% c(0,NA)) && tolower(family) %in% c('poisson','negbin') | any(data$count < 0)){
    errMsg = paste("All values in count should be non-negative integers. Use family = 'gamma' if count is continuous.")
    stop(errMsg)
  }

  ## count.type are integers
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of unsupported type
  if(q==TRUE && !all(data$count.type %% 1 %in% c(0,NA))){
    errMsg = paste("All values in count.type should be integers.")
    stop(errMsg)
  }

  ## phipriors is a vector of two numeric values
  #' @srrstats {G2.0,BS2.2,BS2.3,BS2.4,BS2.5} Checks of vector length and appropriateness of distributional parameters (i.e., vector of length 2, numeric values > 0), implemented prior to analytic routines
  if(!is.numeric(phipriors) | length(phipriors)!=2 | any(phipriors<=0)){
    errMsg = paste("phipriors should be a vector of two positive numeric values. ex. c(0.25,0.25)")
    stop(errMsg)
  }
}
