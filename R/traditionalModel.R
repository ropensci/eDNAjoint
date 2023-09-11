#' Specify and fit model using count data from traditional (non-eDNA) surveys
#'
#' This function implements a Bayesian model that estimates expected species catch rate using count data from traditional (non-eDNA surveys). When multiple traditional gear types are used, an optional variation allows estimation of catchability coefficients, which scale the catchability of gear types relative to the expected catch rate of a reference gear type. Model is implemented using Bayesian inference using the `rstan` package, which uses Hamilton Monte Carlo to simulate the posterior distributions.
#'
#' @export
#' @param data A list containing data necessary for model fitting. Valid tags are `count` and `count.type`. `count` is a matrix or data frame of traditional survey count data, with first dimension equal to the number of sites (i) and second dimension equal to the maximum number of traditional survey replicates at a given site (j). `count.type` is an optional matrix or data frame of integers indicating gear type (k) used in corresponding count data, with first dimension equal to the number of sites (i) and second dimension equal to the maximum number of traditional survey replicates at a given site (j). Values should be integers beginning with 1 (referring to the first gear type) to n (last gear type). Empty cells should be NA. Sites, i, should be consistent in all count data.
#' @param family The distribution class used to model traditional survey count data. Options include poisson ('poisson') and negative binomial ('negbin'). Default value is 'poisson'.
#' @param q A logical value indicating whether to estimate a catchability coefficient, q, for traditional survey gear types (TRUE) or to not estimate a catchability coefficient, q, for traditional survey gear types (FALSE). Default value is FALSE.
#' @param n.chain Number of MCMC chains. Default value is 4.
#' @param n.iter.burn Number of warm-up MCMC iterations. Default value is 500.
#' @param n.iter.sample Number of sampling MCMC iterations. Default value is 2500.
#' @param adapt_delta Target average acceptance probability used in `rstan::sampling`. Default value is 0.9.
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @note  Before fitting the model, this function checks to ensure that the model specification is possible given the data files. These checks include:
#' \itemize{
#' \item  All tags in data are valid (i.e., include count and count.type).
#' \item  Number of sites in count and count type data are equal.
#' \item  All data are numeric (i.e., integer or NA).
#' \item  Empty data cells (NA) match in count and count.type.
#' \item  family is either 'poisson' or 'negbin'.
#' \item  q_ref is an integer >= 0.
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' @examples
#' \donttest{
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
                             q=FALSE,
                             n.chain=4, n.iter.burn=500,
                             n.iter.sample=2500,
                             adapt_delta=0.9) {

  ## #1. make sure all data tags are valid -- if q == TRUE
  if (q==TRUE && !all(c('count.type','count') %in% names(data))){
      errMsg = paste("Data should include 'count' and 'count.type'.")
      stop(errMsg)
  }

  ## #2. make sure all data tags are valid -- if q == FALSE
  if (q==TRUE && !all(c('count') %in% names(data))){
    errMsg = paste("Data should include 'count'.")
    stop(errMsg)
  }

  ## #3. make sure dimensions of count and count.type are equal, if count.type is present
  if (q==TRUE){
    if(dim(data$count)[1]!=dim(data$count.type)[1]|dim(data$count)[2]!=dim(data$count.type)[2]) {
      errMsg = paste("Dimensions of count and count.type do not match.")
      stop(errMsg)
    }
  }

  ## #4. make sure all data is numeric -- if q == TRUE
  if (q==TRUE) {
    if(is.numeric(data$count)==FALSE |
       is.numeric(data$count.type)==FALSE) {
      errMsg = paste("Data should be numeric (i.e. contain integers or NA).")
      stop(errMsg)
    }
  }

  ## #5. make sure all data is numeric -- if q == FALSE
  if (q==FALSE) {
    if(is.numeric(data$count)==FALSE ) {
      errMsg = paste("Data should be numeric (i.e. contain integers or NA).")
      stop(errMsg)
    }
  }

  ## #6. make sure locations of NAs in count data match locations of NAs in count.type data
  if(q==TRUE){
    if(any((which(is.na(data$count.type))==which(is.na(data$count)))==FALSE)){
      errMsg = paste("Empty data cells (NA) in count data should match empty data cells (NA) in count.type data.")
      stop(errMsg)
    }
  }

  ## #7. make sure family is either 'poisson' or 'negbin'
  if(!c(family %in% c('poisson','negbin'))){
    errMsg = paste("Invalid family. Options include 'poisson' and 'negbin'.")
    stop(errMsg)
  }

  ## #8. the smallest count.type is 1
  if(min(data$count.type,na.rm=TRUE) != 1){
    errMsg = paste("The first gear type should be referenced as 1 in count.type. Subsequent gear types should be referenced 2, 3, 4, etc.")
    stop(errMsg)
  }

  ## #9. count.type are integers
  if(!all(data$count.type %% 1 %in% c(0,NA))){
    errMsg = paste("All values in count.type should be integers.")
    stop(errMsg)
  }

  if (!requireNamespace("rstan", quietly = TRUE)){
    stop ("The 'rstan' package is not installed.", call. = FALSE)
  }

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

  ##run model, catchability, pois
  if(q==TRUE&&family=='poisson'){
    out <- rstan::sampling(stanmodels$traditional_catchability_pois,
                           data = list(
                             Nloc = length(unique(count_all$L)),
                             C = nrow(count_all),
                             R = count_all$L,
                             E = count_all$count,
                             nparams = length(q_names),
                             mat = as.matrix(count_all[,q_names]),
                             control = list(adapt_delta = adapt_delta)
                           ),
                           chains = n.chain,
                           thin = 1,
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
                              control = list(adapt_delta = adapt_delta)
                            ),
                            chains = n.chain,
                            thin = 1,
                            warmup = n.iter.burn,
                            iter = n.iter.burn + n.iter.sample,
                            init = init_trad_catchability(n.chain, count_all, q_names)
     )
   } else if(q==FALSE&&family=='poisson'){
     ##run model, no catchability, pois
     out <- rstan::sampling(stanmodels$traditional_pois,
                            data = list(
                              Nloc = length(unique(count_all$L)),
                              C = nrow(count_all),
                              R = count_all$L,
                              E = count_all$count,
                              control = list(adapt_delta = adapt_delta,
                                             stepsize = 0.5)
                            ),
                            chains = n.chain,
                            thin = 1,
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
                              control = list(adapt_delta = adapt_delta)
                            ),
                            chains = n.chain,
                            thin = 1,
                            warmup = n.iter.burn,
                            iter = n.iter.burn + n.iter.sample,
                            init = init_trad(n.chain, count_all)
     )
   }

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
