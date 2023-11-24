#' Specify and fit joint model using count data from traditional surveys and eDNA qPCR data
#'
#' This function implements a Bayesian model that integrates data from paired eDNA and traditional surveys. The model estimates parameters including the expected species catch rate and the probability of false positive eDNA detection. This function allows for optional model variations, like inclusion of site-level covariates that scale the sensitivity of eDNA sampling relative to traditional sampling, as well as estimation of catchability coefficients when multiple traditional gear types are used. Model is implemented using Bayesian inference using the `rstan` package, which uses Hamiltonian Monte Carlo to simulate the posterior distributions.
#'
#' @export
#' @param data A list containing data necessary for model fitting. Valid tags are `qPCR.N`, `qPCR.K`, `count`, `count.type`, and `site.cov`. `qPCR.N` and `qPCR.K` are matrices or data frames with first dimension equal to the number of sites (i) and second dimension equal to the maximum number of eDNA samples at a given site (m). `qPCR.N` contains the total number of qPCR replicates per site and eDNA sample, and `qPCR.K` contains the total number of positive qPCR detections per site and eDNA sample. `count` is a matrix or data frame of traditional survey count data, with first dimension equal to the number of sites (i) and second dimension equal to the maximum number of traditional survey replicates at a given site (j). `count.type` is an optional matrix or data frame of integers indicating gear type used in corresponding count data, with first dimension equal to the number of sites (i) and second dimension equal to the maximum number of traditional survey replicates at a given site. Values should be integers beginning with 1 (referring to the first gear type) to n (last gear type). `site.cov` is an optional matrix or data frame of site-level covariate data, with first dimension equal to the number of sites (i). `site.cov` should include informative column names. Empty cells should be NA. Sites, i, should be consistent in all qPCR, count, and site covariate data.
#' @param cov A character vector indicating the site-level covariates to include in the model. Default value is 'None'.
#' @param family The distribution class used to model traditional survey count data. Options include poisson ('poisson'), negative binomial ('negbin'), and gamma ('gamma'). Default value is 'poisson'.
#' @param p10priors A numeric vector indicating beta distribution parameters (alpha, beta) used as the prior distribution for the eDNA false positive probability (p10). Default vector is c(1,20).
#' @param q A logical value indicating whether to estimate a catchability coefficient, q, for traditional survey gear types (TRUE) or to not estimate a catchability coefficient, q, for traditional survey gear types (FALSE). Default value is FALSE.
#' @param phipriors A numeric vector indicating gamma distribution parameters (shape, rate) used as the prior distribution for phi, the overdispersion in the negative binomial distribution for traditional survey gear data. Used when family = 'negbin.' Default vector is c(0.25,0.25).
#' @param n.chain Number of MCMC chains. Default value is 4.
#' @param n.iter.burn Number of warm-up MCMC iterations. Default value is 500.
#' @param n.iter.sample Number of sampling MCMC iterations. Default value is 2500.
#' @param thin A positive integer specifying the period for saving samples. Default value is 1.
#' @param adapt_delta Target average acceptance probability used in `rstan::sampling`. Default value is 0.9.
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
#' @note  Before fitting the model, this function checks to ensure that the model specification is possible given the data files. These checks include:
#' \itemize{
#' \item  All tags in data are valid (i.e., include qPCR.N, qPCR.K, count, count.type, and site.cov).
#' \item  Dimensions of qPCR.N and qPCR.K are equal, and dimensions of count and count.type are equal (if count.type provided).
#' \item  Number of sites in qPCR and count data are equal.
#' \item  All data are numeric (i.e., integer or NA).
#' \item  Empty data cells (NA) match in qPCR.N and qPCR.K and in count and count.type.
#' \item  family is either 'poisson', 'negbin', or 'gamma'.
#' \item  p10priors is a vector of two numeric values.
#' \item  q_ref is an integer >= 0.
#' \item  site.cov has same number of rows as qPCR.N and count, if present
#' \item  site.cov is numeric, if present
#' \item  cov values match column names in site.cov, if present
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' @examples
#' \donttest{
#' # Ex. 1: Implementing the joint model with site-level covariates
#'
#' # Load data
#' data(gobyData)
#'
#' # Examine data in list
#' names(gobyData)
#'
#' # Note that the surveyed sites (rows) should match in all data
#' dim(gobyData$qPCR.N)[1]
#' dim(gobyData$count)[1]
#'
#' # Fit a model including 'Filter_time' and 'Salinity' site-level covariates
#' # These covariates will scale the sensitivity of eDNA sampling relative to traditional surveys
#' # Count data is modeled using a poisson distribution.
#' fit.cov = jointModel(data=gobyData, cov=c('Filter_time','Salinity'),
#'                      family='poisson', p10priors=c(1,20), q=FALSE)
#'
#'
#' # Ex. 2: Implementing the joint model with multiple traditional gear types
#'
#' # Load data
#' data(greencrabData)
#'
#' # Examine data in list
#' names(greencrabData)
#'
#' # Note that the surveyed sites (rows) should match in all data
#' dim(greencrabData$qPCR.N)[1]
#' dim(greencrabData$count)[1]
#'
#' # Fit a model estimating a catchability coefficient for traditional survey gear types.
#' # This model does not assume all traditional survey methods have the same catchability.
#' # Count data is modeled using a negative binomial distribution.
#' fit.q = jointModel(data=greencrabData, cov='None', family='negbin',
#'                    p10priors=c(1,20), q=TRUE)
#' }
#'

jointModel <- function(data, cov='None', family='poisson', p10priors=c(1,20), q=FALSE,
                       phipriors=c(0.25,0.25),
                       n.chain=4, n.iter.burn=500,
                       n.iter.sample=2500, thin=1, adapt_delta=0.9) {

  ####
  # input checks

  # all models
  all_checks(data,cov,family,p10priors)

  # model with catchability coefficients
  if (q==TRUE) {
    catchability_checks(data,cov)
  }

  # model without catchability coefficients
  if (q==FALSE) {
    no_catchability_checks(data,cov)
  }

  # model with covariates
  if (all(cov!='None')) {
    covariate_checks(data,cov)
  }

  if (!requireNamespace("rstan", quietly = TRUE)){
    stop ("The 'rstan' package is not installed.", call. = FALSE)
  }

  ###
  #convert data to long format

  '%>%' <- magrittr::`%>%`

  #convert qPCR data to long format
  qPCR_all <- as.data.frame(data$qPCR.N) %>%
    dplyr::mutate(L=1:dim(data$qPCR.N)[1]) %>%
    tidyr::pivot_longer(cols=!L,values_to='N') %>%
    tidyr::drop_na()
  qPCR.K_df <- as.data.frame(data$qPCR.K) %>%
    dplyr::mutate(L=1:dim(data$qPCR.K)[1]) %>%
    tidyr::pivot_longer(cols=!L,values_to='K') %>%
    tidyr::drop_na()
  qPCR_all$K <- qPCR.K_df$K

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

  #if present, prepare covariate data
  if(all(cov!='None')){
    site_mat <- as.data.frame(data$site.cov[,cov])
    site_mat <- cbind(as.data.frame(rep(1,length(site_mat[,1]))),site_mat)
    colnames(site_mat) <- c('int',cov)
  }

  #convert p10 prior
  #p10 prior: convert beta(1,20) to lognormal distribution
  #moment match from beta(alpha,beta) to normal(mu, sigma^2)
  alpha <- p10priors[1]
  beta <- p10priors[2]
  mu <- alpha/(alpha+beta)
  sigma_2 <- (alpha*beta)/((alpha+beta)^2*(alpha+beta+1))
  #convert normal(mu, sigma^2) to lognormal(mu, sigma^2)
  mu_ln <- log(mu^2/sqrt(mu^2+sigma_2))
  sigma_2_ln <- log(1+sigma_2/mu^2)
  sigma_ln <- sqrt(sigma_2_ln)

  # create data that will be present in all model variations
  data <- list(
    S = nrow(qPCR_all),
    L = qPCR_all$L,
    Nloc = length(unique(qPCR_all$L)),
    N = qPCR_all$N,
    K = qPCR_all$K,
    C = nrow(count_all),
    R = count_all$L,
    E = count_all$count,
    p10priors = c(mu_ln,sigma_ln),
    control = list(adapt_delta = adapt_delta)
    )


  ##run model, catchability, pois/gamma, no covariates
  if(q==TRUE&&family!='negbin'&&all(cov=='None')){
    model_index <- dplyr::case_when(family=='poisson'~ 1,
                                    family=='gamma' ~ 2)
    out <- rstan::sampling(c(stanmodels$joint_binary_catchability_pois,
                             stanmodels$joint_binary_catchability_gamma)[model_index][[1]],
                           data = rlist::list.append(
                             data,
                             nparams = length(q_names),
                             mat = as.matrix(count_all[,q_names])
                           ),
                           chains = n.chain,
                           thin = thin,
                           warmup = n.iter.burn,
                           iter = n.iter.burn + n.iter.sample,
                           init = init_joint_catchability(n.chain,count_all,q_names)
    )
  } else if(q==TRUE&&family=='negbin'&&all(cov=='None')){
    ##run model, catchability, negbin, no covariates
    out <- rstan::sampling(stanmodels$joint_binary_catchability_negbin,
                           data = rlist::list.append(
                             data,
                             phipriors = phipriors,
                             nparams = length(q_names),
                             mat = as.matrix(count_all[,q_names])
                           ),
                           chains = n.chain,
                           thin = thin,
                           warmup = n.iter.burn,
                           iter = n.iter.burn + n.iter.sample,
                           init = init_joint_catchability(n.chain,count_all,q_names)
    )
  } else if(q==FALSE&&family!='negbin'&&all(cov=='None')){
    ##run model, no catchability, pois/gamma, no covariates
      model_index <- dplyr::case_when(family=='poisson'~ 1,
                                      family=='gamma' ~ 2)
      out <- rstan::sampling(c(stanmodels$joint_binary_pois,
                               stanmodels$joint_binary_gamma)[model_index][[1]],
                           data = data,
                           chains = n.chain,
                           thin = thin,
                           warmup = n.iter.burn,
                           iter = n.iter.burn + n.iter.sample,
                           init = init_joint(n.chain,count_all)
    )
  } else if(q==FALSE&&family=='negbin'&&all(cov=='None')){
    ##run model, no catchability, negbin, no covariates
    out <- rstan::sampling(stanmodels$joint_binary_negbin,
                           data = rlist::list.append(
                             data,
                             phipriors = phipriors
                           ),
                           chains = n.chain,
                           thin = thin,
                           warmup = n.iter.burn,
                           iter = n.iter.burn + n.iter.sample,
                           init = init_joint(n.chain,count_all)
    )
  } else if(q==TRUE&&family=='negbin'&&all(cov!='None')){
    ##run model, catchability, negbin, covariates
    out <- rstan::sampling(stanmodels$joint_binary_cov_catchability_negbin,
                           data = rlist::list.append(
                             data,
                             phipriors = phipriors,
                             nparams = length(q_names),
                             mat = as.matrix(count_all[,q_names]),
                             nsitecov = length(cov)+1,
                             mat_site = as.matrix(site_mat),
                           ),
                           chains = n.chain,
                           thin = thin,
                           warmup = n.iter.burn,
                           iter = n.iter.burn + n.iter.sample,
                           init = init_joint_cov_catchability(n.chain,count_all,q_names,cov)
    )
  } else if(q==TRUE&&family!='negbin'&&all(cov!='None')){
    ##run model, catchability, pois/gamma, covariates
    model_index <- dplyr::case_when(family=='poisson'~ 1,
                                    family=='gamma' ~ 2)
    out <- rstan::sampling(c(stanmodels$joint_binary_cov_catchability_pois,
                             stanmodels$joint_binary_cov_catchability_gamma)[model_index][[1]],
                           data = rlist::list.append(
                             data,
                             nparams = length(q_names),
                             mat = as.matrix(count_all[,q_names]),
                             nsitecov = length(cov)+1,
                             mat_site = as.matrix(site_mat),
                           ),
                           chains = n.chain,
                           thin = thin,
                           warmup = n.iter.burn,
                           iter = n.iter.burn + n.iter.sample,
                           init = init_joint_cov_catchability(n.chain,count_all,q_names,cov)
    )
  } else if(q==FALSE&&family=='negbin'&&all(cov!='None')){
    ##run model, no catchability, negbin, covariates
    out <- rstan::sampling(stanmodels$joint_binary_cov_negbin,
                           data = rlist::list.append(
                             data,
                             phipriors = phipriors,
                             nsitecov = length(cov)+1,
                             mat_site = as.matrix(site_mat),
                           ),
                           chains = n.chain,
                           thin = thin,
                           warmup = n.iter.burn,
                           iter = n.iter.burn + n.iter.sample,
                           init = init_joint_cov(n.chain,count_all,cov)
    )
  } else if(q==FALSE&&family!='negbin'&&all(cov!='None')){
    ##run model, no catchability, pois/gamma, covariates
    model_index <- dplyr::case_when(family=='poisson'~ 1,
                                    family=='gamma' ~ 2)
    out <- rstan::sampling(c(stanmodels$joint_binary_cov_pois,
                             stanmodels$joint_binary_cov_gamma)[model_index][[1]],
                           data = rlist::list.append(
                             data,
                             nsitecov = length(cov)+1,
                             mat_site = as.matrix(site_mat),
                           ),
                           chains = n.chain,
                           thin = thin,
                           warmup = n.iter.burn,
                           iter = n.iter.burn + n.iter.sample,
                           init = init_joint_cov(n.chain,count_all,cov)
    )
  }

  return(out)
}

###########
#helper functions: initial values
###########

init_joint_cov <- function(n.chain,count_all,cov){
  #helper function
  #joint model, catchability coefficient, site covariates
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      mu = stats::runif(length(unique(count_all$L)), 0.01, 5),
      p10 = log(0.05),
      alpha = rep(0.1,length(cov)+1)
    )
  }
  return(A)
}

init_joint_cov_catchability <- function(n.chain,count_all,q_names,cov){
  #helper function
  #joint model, catchability coefficient, site covariates
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      mu = stats::runif(length(unique(count_all$L)), 0.01, 5),
      q = as.data.frame(stats::runif(length(q_names),0.01,1)),
      p10 = log(0.05),
      alpha = rep(0.1,length(cov)+1)
    )
  }
  return(A)
}

init_joint_catchability <- function(n.chain,count_all,q_names){
  #helper function
  #joint model, catchability coefficient, no site covariates
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      mu = stats::runif(length(unique(count_all$L)), 0.01, 5),
      q = as.data.frame(stats::runif(length(q_names),0.01,1)),
      p10 = log(0.05),
      beta = .5
    )
  }
  return(A)
}

init_joint <- function(n.chain,count_all){
  #helper function
  #joint model, no catchability coefficient, no site covariates
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      mu = stats::runif(length(unique(count_all$L)), 0.01, 5),
      p10 = log(0.05),
      beta = .5
    )
  }
  return(A)
}

################
#helper functions: input checks
################

# if q = TRUE
catchability_checks <- function(data,cov){
  ## All tags in data are valid (i.e., include qPCR.N, qPCR.K, count, count.type, and site.cov)
  #cov='None'
  if (all(cov=='None') && !all(c('qPCR.N', 'qPCR.K', 'count','count.type') %in% names(data))){
    errMsg = paste("Data should include 'qPCR.N', 'qPCR.K', 'count', and 'count.type'.")
    stop(errMsg)
  }
  #q=TRUE and cov!='None'
  if (all(cov!='None') && !all(c('qPCR.N', 'qPCR.K', 'count','count.type','site.cov') %in% names(data))){
    errMsg = paste("Data should include 'qPCR.N', 'qPCR.K', 'count', 'count.type', and 'site.cov'.")
    stop(errMsg)
  }
  ## make sure dimensions of count and count.type are equal, if count.type is present
  if(dim(data$count)[1]!=dim(data$count.type)[1]|dim(data$count)[2]!=dim(data$count.type)[2]) {
      errMsg = paste("Dimensions of count and count.type do not match.")
      stop(errMsg)
  }

  ## make sure all data is numeric -- if q == TRUE
  if(is.numeric(data$qPCR.K)==FALSE |
       is.numeric(data$qPCR.N)==FALSE |
       is.numeric(data$count)==FALSE |
       is.numeric(data$count.type)==FALSE) {
      errMsg = paste("Data should be numeric.")
      stop(errMsg)
  }
  ## make sure locations of NAs in count data match locations of NAs in count.type data
  if(sum(is.na(data$count.type))!=sum(is.na(data$count))){
      errMsg = paste("Empty data cells (NA) in count data should match empty data cells (NA) in count.type data.")
      stop(errMsg)
  }
  if(any((which(is.na(data$count))==which(is.na(data$count.type)))==FALSE)){
      errMsg = paste("Empty data cells (NA) in count data should match empty data cells (NA) in count.type data.")
      stop(errMsg)
  }
  ## the smallest count.type is 1
  if(min(data$count.type,na.rm=TRUE) != 1){
    errMsg = paste("The first gear type should be referenced as 1 in count.type. Subsequent gear types should be referenced 2, 3, 4, etc.")
    stop(errMsg)
  }

  ## count.type are integers
  if(!all(data$count.type %% 1 %in% c(0,NA))){
    errMsg = paste("All values in count.type should be integers.")
    stop(errMsg)
  }
}

no_catchability_checks <- function(data,cov){
  ## All tags in data are valid (i.e., include qPCR.N, qPCR.K, count, and site.cov)
  #cov='None'
  if (all(cov=='None') && !all(c('qPCR.N', 'qPCR.K', 'count') %in% names(data))){
    errMsg = paste("Data should include 'qPCR.N', 'qPCR.K', and 'count'.")
    stop(errMsg)
  }
  #cov!='None'
  if (all(cov!='None') && !all(c('qPCR.N', 'qPCR.K', 'count','site.cov') %in% names(data))){
    errMsg = paste("Data should include 'qPCR.N', 'qPCR.K', 'count', and 'site.cov'.")
    stop(errMsg)
  }

  ## make sure all data is numeric -- if q == FALSE
  if(is.numeric(data$qPCR.K)==FALSE |
       is.numeric(data$qPCR.N)==FALSE |
       is.numeric(data$count)==FALSE ) {
      errMsg = paste("Data should be numeric (i.e. contain integers or NA).")
      stop(errMsg)
  }
}

all_checks <- function(data,cov,family,p10priors){
  ## make sure dimensions of qPCR.N and qPCR.K are equal
  if (dim(data$qPCR.N)[1]!=dim(data$qPCR.K)[1]|dim(data$qPCR.N)[2]!=dim(data$qPCR.K)[2]) {
    errMsg = paste("Dimensions of qPCR.N and qPCR.K do not match.")
    stop(errMsg)
  }
  ## make sure number of rows in count = number of rows in qPCR.N and qPCR.K
  if (dim(data$qPCR.N)[1]!=dim(data$count)[1]) {
    errMsg = paste("Number of sites (rows) in qPCR data and traditional survey count data do not match.")
    stop(errMsg)
  }

  ## make sure locations of NAs in qPCR.N data match locations of NAs in qPCR.K data
  if(any((which(is.na(data$qPCR.N))==which(is.na(data$qPCR.K)))==FALSE)){
    errMsg = paste("Empty data cells (NA) in qPCR.N data should match empty data cells (NA) in qPCR.K data.")
    stop(errMsg)
  }

  ## make sure family is either 'poisson', 'negbin', or 'gamma'
  if(!c(family %in% c('poisson','negbin','gamma'))){
    errMsg = paste("Invalid family. Options include 'poisson', 'negbin', and 'gamma'.")
    stop(errMsg)
  }

  ## p10 priors is a vector of two integers
  if(!is.numeric(p10priors) | length(p10priors)!=2){
    errMsg = paste("p10priors should be a vector of two integers. ex. c(1,20)")
    stop(errMsg)
  }


  ## count are integers, if family is poisson or negbin
  if(!all(data$count %% 1 %in% c(0,NA)) && family %in% c('poisson','negbin')){
    errMsg = paste("All values in count should be integers. Use family = 'gamma' if count is continuous.")
    stop(errMsg)
  }

  ## qPCR.N are integers
  if(!all(data$qPCR.N %% 1 %in% c(0,NA))){
    errMsg = paste("All values in qPCR.N should be integers.")
    stop(errMsg)
  }

  ## qPCR.K are integers
  if(!all(data$qPCR.K %% 1 %in% c(0,NA))){
    errMsg = paste("All values in qPCR.K should be integers.")
    stop(errMsg)
  }

}

covariate_checks <- function(data,cov){
  ## site.cov is numeric, if present
  if(!is.numeric(data$site.cov)){
    errMsg = paste("site.cov should be numeric.")
    stop(errMsg)
  }

  ## cov values match column names in site.cov
  if(!all(cov %in% colnames(data$site.cov))){
    errMsg = paste("cov values should be listed in the column names of site.cov in the data.")
    stop(errMsg)
  }

  ## site.cov has same number of rows as qPCR.N and count, if present
  if(dim(data$qPCR.N)[1]!=dim(data$site.cov)[1]){
    errMsg = paste("The number of rows in site.cov matrix should match the number of rows in all other matrices.")
    stop(errMsg)
  }
}


