#' Specify and fit joint model using count data from traditional surveys and
#' eDNA qPCR data
#'
#' @srrstats {G1.0,G1.1} This software makes available a novel algorithm/model
#' that was previously published in the literature. The literature reference for
#' the joint model is provided here.
#' This function implements a Bayesian model that integrates data from paired
#' eDNA and traditional surveys, as described in Keller et al. (2022)
#' <https://doi.org/10.1002/eap.2561>. The model estimates parameters including
#' the expected species catch rate and the probability of false positive eDNA
#' detection. This function allows for optional model variations, like inclusion
#' of site-level covariates that scale the sensitivity of eDNA sampling relative
#' to traditional sampling, as well as estimation of catchability coefficients
#' when multiple traditional gear types are used. Model is implemented using
#' Bayesian inference using the `rstan` package, which uses Hamiltonian Monte
#' Carlo to simulate the posterior distributions.
#'
#' @srrstats {G1.4,G1.3} Roxygen function documentation begins here, with
#' definitions of statistical terminology and inputs
#' @export
#' @srrstats {BS1.1,BS3.0,G2.14} Descriptions of how to enter data, description
#' of how NAs are handled (which are informative and should be deliberately
#' included in input data, where necessary)
#' @param data A list containing data necessary for model fitting. Valid tags
#' are `qPCR.N`, `qPCR.K`, `count`, `count.type`, and `site.cov`. `qPCR.N` and
#' `qPCR.K` are matrices or data frames with first dimension equal to the number
#' of sites (i) and second dimension equal to the maximum number of eDNA samples
#' at a given site (m). `qPCR.N` contains the total number of qPCR replicates
#' per site and eDNA sample, and `qPCR.K` contains the total number of positive
#' qPCR detections per site and eDNA sample. `count` is a matrix or data frame
#' of traditional survey count data, with first dimension equal to the number of
#' sites (i) and second dimension equal to the maximum number of traditional
#' survey replicates at a given site (j). `count.type` is an optional matrix or
#' data frame of integers indicating gear type used in corresponding count data,
#' with first dimension equal to the number of sites (i) and second dimension
#' equal to the maximum number of traditional survey replicates at a given site.
#' Values should be integers beginning with 1 (referring to the first gear type)
#' to n (last gear type). `site.cov` is an optional matrix or data frame of
#' site-level covariate data, with first dimension equal to the number of sites
#' (i). `site.cov` should include informative column names. Empty cells should
#' be NA and will be removed during processing. Sites, i, should be consistent
#' in all qPCR, count, and site covariate data.
#' @srrstats {G2.1a} Here are explicit documentation of vector input types
#' @param cov A character vector indicating the site-level covariates to include
#' in the model. Default value is 'None'.
#' @param family The distribution class used to model traditional survey count
#' data. Options include poisson ('poisson'), negative binomial ('negbin'), and
#' gamma ('gamma'). Default value is 'poisson'.
#' @srrstats {BS1.0,G2.1a,BS1.2} Description of hyperparameters and how to
#' specify prior distributions, explicit documentation of vector input types
#' @param p10priors A numeric vector indicating beta distribution
#' hyperparameters (alpha, beta) used as the prior distribution for the eDNA
#' false positive probability (p10). Default vector is c(1,20).
#' @param q A logical value indicating whether to estimate a catchability
#' coefficient, q, for traditional survey gear types (TRUE) or to not estimate a
#' catchability coefficient, q, for traditional survey gear types (FALSE).
#' Default value is FALSE.
#' @srrstats {BS1.0,G2.1a,BS1.2} Description of hyperparameters and how to
#' specify prior distributions, explicit documentation of vector input types
#' @param phipriors A numeric vector indicating gamma distribution
#' hyperparameters (shape, rate) used as the prior distribution for phi, the
#' overdispersion in the negative binomial distribution for traditional survey
#' gear data. Used when family = 'negbin.' Default vector is c(0.25,0.25).
#' @param multicore A logical value indicating whether to parallelize chains
#' with multiple cores. Default is TRUE.
#' @srrstats {BS2.7,BS2.11} Option for user to provide initial values for each
#' chain
#' @param initial_values A list of lists of initial values to use in MCMC. The
#' length should equal the number of MCMC chains. Initial values can be provided
#'  for parameters: beta, p10 (log-scale), mu, q, alpha. If no initial values
#'  are provided, default random values are drawn.
#' @srrstats {BS1.3} Description of parameters used in the computational process
#' begins here
#' @param n.chain Number of MCMC chains. Default value is 4.
#' @param n.iter.burn Number of warm-up MCMC iterations. Default value is 500.
#' @param n.iter.sample Number of sampling MCMC iterations. Default value is
#' 2500.
#' @param thin A positive integer specifying the period for saving samples.
#' Default value is 1.
#' @param adapt_delta Target average acceptance probability used in
#' `rstan::sampling`. Default value is 0.9.
#' @srrstats {BS2.12} Parameter controlling the verbosity of output
#' @param verbose Logical value controlling the verbosity of output (i.e.,
#' warnings, messages, progress bar). Default is TRUE.
#' @srrstats {BS5.0} function returns initial values used in computation
#' @return A list of:
#' \itemize{
#' \item a model object of class `stanfit` returned by `rstan::sampling`
#' \item initial values used in MCMC
#' }
#'
#' @srrstats {G2.0a} Explicit secondary documentation of any expectations on
#' lengths of inputs
#' @note  Before fitting the model, this function checks to ensure that the
#' model specification is possible given the data files. These checks include:
#' \itemize{
#' \item  All tags in data are valid (i.e., include qPCR.N, qPCR.K, count,
#' count.type, and site.cov).
#' \item  Dimensions of qPCR.N and qPCR.K are equal, and dimensions of count and
#' count.type are equal (if count.type provided).
#' \item  Number of sites in qPCR and count data are equal.
#' \item  All data are numeric (i.e., integer or NA).
#' \item  Empty data cells (NA) match in qPCR.N and qPCR.K and in count and
#' count.type.
#' \item  family is either 'poisson', 'negbin', or 'gamma'.
#' \item  p10priors and phipriors (if used) is a vector of two numeric values.
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
#' # These covariates will scale the sensitivity of eDNA sampling relative to
#' # traditional surveys
#' # Count data is modeled using a poisson distribution.
#' fit.cov = jointModel(data=gobyData, cov=c('Filter_time','Salinity'),
#'                      family="poisson", p10priors=c(1,20), q=FALSE)
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
#' # Fit a model estimating a catchability coefficient for traditional survey
#' # gear types.
#' # This model does not assume all traditional survey methods have the same
#' # catchability.
#' # Count data is modeled using a negative binomial distribution.
#' fit.q = jointModel(data=greencrabData, cov="None", family="negbin",
#'                    p10priors=c(1,20), q=TRUE)
#' }
#'

jointModel <- function(data, cov='None', family='poisson', p10priors=c(1,20),
                       q=FALSE, phipriors=c(0.25,0.25), multicore=TRUE,
                       initial_values='None', n.chain=4, n.iter.burn=500,
                       n.iter.sample=2500, thin=1, adapt_delta=0.9,
                       verbose=TRUE) {

  ####
  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using the following
  #' helper functions

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

  # all models
  all_checks(data,cov,family,p10priors,phipriors)

  # initial value checks
  if(all(initial_values != 'None')){
    initial_values_checks(initial_values,data,cov,n.chain)
  }

  if (!requireNamespace("rstan", quietly = TRUE)){
    stop ("The 'rstan' package is not installed.", call. = FALSE)
  }

  # make character inputs case-insensitive
  #' @srrstats {G2.3b} Allow case-insensitive character parameter values
  family <- tolower(family)

  ###
  #convert data to long format

  '%>%' <- magrittr::`%>%`

  #convert qPCR data to long format
  #' @srrstats {G2.7} Use as.data.frame() to allow input list of any tabular
  #' form (i.e., matrix, etc.)
  qPCR_all <- as.data.frame(data$qPCR.N) %>%
    dplyr::mutate(L=1:dim(data$qPCR.N)[1]) %>%
    tidyr::pivot_longer(cols=!L,values_to='N') %>%
    #' @srrstats {G2.15} Software does not assume non-missingness and actually
    #' expects it if the number of observations across sites is unequal
    tidyr::drop_na()
  qPCR.K_df <- as.data.frame(data$qPCR.K) %>%
    dplyr::mutate(L=1:dim(data$qPCR.K)[1]) %>%
    tidyr::pivot_longer(cols=!L,values_to='K') %>%
    #' @srrstats {G2.15} Software does not assume non-missingness and actually
    #' expects it if the number of observations across sites is unequal
    tidyr::drop_na()
  qPCR_all$K <- qPCR.K_df$K

  #convert count data to long format
  #' @srrstats {G2.7} Use as.data.frame() to allow input list of any tabular
  #' form (i.e., matrix, etc.)
  count_all <- as.data.frame(data$count) %>%
    dplyr::mutate(L=1:dim(data$count)[1]) %>%
    tidyr::pivot_longer(cols=!L,values_to='count') %>%
    #' @srrstats {G2.15} Software does not assume non-missingness and actually
    #' expects it if the number of observations across sites is unequal
    tidyr::drop_na()

  #if q==TRUE, add count type data to count df
  if(q==TRUE){
    q_ref <- 1
    #' @srrstats {G2.7} Use as.data.frame() to allow input list of any tabular
    #' form (i.e., matrix, etc.)
    count.type_df <- as.data.frame(data$count.type) %>%
      dplyr::mutate(L=1:dim(data$count.type)[1]) %>%
      tidyr::pivot_longer(cols=!L,values_to='count.type') %>%
      tidyr::drop_na()
    count_all$count.type <- count.type_df$count.type

    #create vector of q coefficient names
    counttypes <- unique(count_all$count.type)
    names <- counttypes[!counttypes==q_ref]
    #' @srrstats {G2.4,G2.4c} Explicit conversion to character
    q_names <- as.character(paste0('q_',names))

    #add dummy variables for count type
    for(i in seq_along(q_names)){
      count_all[,q_names[i]] <- ifelse(count_all$count.type==names[i],1,0)
    }
  }

  #if present, prepare covariate data
  if(all(cov!='None')){
    #' @srrstats {G2.7,G2.10} Use as.data.frame() to allow input list of any
    #' tabular form (i.e., matrix, etc.) and converts before filtering columns
    #' based on input 'cov'
    site_mat <- as.data.frame(as.data.frame(data$site.cov)[,cov])
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

  # set up core number
  if(multicore == TRUE){
    cores <- parallel::detectCores()
  } else {
    cores <- 1
  }


  ##run model, catchability, pois/gamma, no covariates
  if(q==TRUE&&family!='negbin'&&all(cov=='None')){
    model_index <- dplyr::case_when(family=='poisson'~ 1,
                                    family=='gamma' ~ 2)
    inits <- init_joint_catchability(n.chain,count_all,q_names,initial_values)
    out <- rstan::sampling(
      c(stanmodels$joint_binary_catchability_pois,
        stanmodels$joint_binary_catchability_gamma)[model_index][[1]],
      data = rlist::list.append(
        data,
        nparams = length(q_names),
        mat = as.matrix(count_all[,q_names])
      ),
      cores = cores,
      #' @srrstats {G2.4,G2.4a} explicit conversion to
      #' integers for sampling arguments
      chains = as.integer(n.chain),
      thin = as.integer(thin),
      warmup = as.integer(n.iter.burn),
      iter = (
        as.integer(n.iter.burn) + as.integer(n.iter.sample)
      ),
      init = inits,
      refresh = ifelse(verbose==TRUE,500,0)
    )
  } else if(q==TRUE&&family=='negbin'&&all(cov=='None')){
    ##run model, catchability, negbin, no covariates
    inits <- init_joint_catchability(n.chain,count_all,q_names,initial_values)
    out <- rstan::sampling(stanmodels$joint_binary_catchability_negbin,
                           data = rlist::list.append(
                             data,
                             phipriors = phipriors,
                             nparams = length(q_names),
                             mat = as.matrix(count_all[,q_names])
                           ),
                           cores = cores,
                           #' @srrstats {G2.4,G2.4a} explicit conversion to
                           #' integers for sampling arguments
                           chains = as.integer(n.chain),
                           thin = as.integer(thin),
                           warmup = as.integer(n.iter.burn),
                           iter = (
                             as.integer(n.iter.burn) + as.integer(n.iter.sample)
                           ),
                           init = inits,
                           refresh = ifelse(verbose==TRUE,500,0)
    )
  } else if(q==FALSE&&family!='negbin'&&all(cov=='None')){
    ##run model, no catchability, pois/gamma, no covariates
    model_index <- dplyr::case_when(family=='poisson'~ 1,
                                    family=='gamma' ~ 2)
    inits <- init_joint(n.chain,count_all,initial_values)
    out <- rstan::sampling(c(stanmodels$joint_binary_pois,
                             stanmodels$joint_binary_gamma)[model_index][[1]],
                           data = data,
                           cores = cores,
                           #' @srrstats {G2.4,G2.4a} explicit conversion to
                           #' integers for sampling arguments
                           chains = as.integer(n.chain),
                           thin = as.integer(thin),
                           warmup = as.integer(n.iter.burn),
                           iter = (
                             as.integer(n.iter.burn) + as.integer(n.iter.sample)
                           ),
                           init = inits,
                           refresh = ifelse(verbose==TRUE,500,0)
    )
  } else if(q==FALSE&&family=='negbin'&&all(cov=='None')){
    ##run model, no catchability, negbin, no covariates
    inits <- init_joint(n.chain,count_all,initial_values)
    out <- rstan::sampling(stanmodels$joint_binary_negbin,
                           data = rlist::list.append(
                             data,
                             phipriors = phipriors
                           ),
                           cores = cores,
                           #' @srrstats {G2.4,G2.4a} explicit conversion to
                           #' integers for sampling arguments
                           chains = as.integer(n.chain),
                           thin = as.integer(thin),
                           warmup = as.integer(n.iter.burn),
                           iter = (
                             as.integer(n.iter.burn) + as.integer(n.iter.sample)
                           ),
                           init = inits,
                           refresh = ifelse(verbose==TRUE,500,0)
    )
  } else if(q==TRUE&&family=='negbin'&&all(cov!='None')){
    ##run model, catchability, negbin, covariates
    inits <- init_joint_cov_catchability(n.chain,count_all,
                                         q_names,cov,initial_values)
    out <- rstan::sampling(stanmodels$joint_binary_cov_catchability_negbin,
                           data = rlist::list.append(
                             data,
                             phipriors = phipriors,
                             nparams = length(q_names),
                             mat = as.matrix(count_all[,q_names]),
                             nsitecov = length(cov)+1,
                             mat_site = as.matrix(site_mat)
                           ),
                           cores = cores,
                           #' @srrstats {G2.4,G2.4a} explicit conversion to
                           #' integers for sampling arguments
                           chains = as.integer(n.chain),
                           thin = as.integer(thin),
                           warmup = as.integer(n.iter.burn),
                           iter = (
                             as.integer(n.iter.burn) + as.integer(n.iter.sample)
                           ),
                           init = inits,
                           refresh = ifelse(verbose==TRUE,500,0)
    )
  } else if(q==TRUE&&family!='negbin'&&all(cov!='None')){
    ##run model, catchability, pois/gamma, covariates
    model_index <- dplyr::case_when(family=='poisson'~ 1,
                                    family=='gamma' ~ 2)
    inits <- init_joint_cov_catchability(n.chain,count_all,
                                         q_names,cov,initial_values)
    out <- rstan::sampling(
      c(stanmodels$joint_binary_cov_catchability_pois,
        stanmodels$joint_binary_cov_catchability_gamma)[model_index][[1]],
      data = rlist::list.append(
        data,
        nparams = length(q_names),
        mat = as.matrix(count_all[,q_names]),
        nsitecov = length(cov)+1,
        mat_site = as.matrix(site_mat)
      ),
      cores = cores,
      #' @srrstats {G2.4,G2.4a} explicit conversion to
      #' integers for sampling arguments
      chains = as.integer(n.chain),
      thin = as.integer(thin),
      warmup = as.integer(n.iter.burn),
      iter = (
        as.integer(n.iter.burn) + as.integer(n.iter.sample)
      ),
      init = inits,
      refresh = ifelse(verbose==TRUE,500,0)
    )
  } else if(q==FALSE&&family=='negbin'&&all(cov!='None')){
    ##run model, no catchability, negbin, covariates
    inits <- init_joint_cov(n.chain,count_all,cov,initial_values)
    out <- rstan::sampling(stanmodels$joint_binary_cov_negbin,
                           data = rlist::list.append(
                             data,
                             phipriors = phipriors,
                             nsitecov = length(cov)+1,
                             mat_site = as.matrix(site_mat)
                           ),
                           cores = cores,
                           #' @srrstats {G2.4,G2.4a} explicit conversion to
                           #' integers for sampling arguments
                           chains = as.integer(n.chain),
                           thin = as.integer(thin),
                           warmup = as.integer(n.iter.burn),
                           iter = (
                             as.integer(n.iter.burn) + as.integer(n.iter.sample)
                           ),
                           init = inits,
                           refresh = ifelse(verbose==TRUE,500,0)
    )
  } else if(q==FALSE&&family!='negbin'&&all(cov!='None')){
    ##run model, no catchability, pois/gamma, covariates
    model_index <- dplyr::case_when(family=='poisson'~ 1,
                                    family=='gamma' ~ 2)
    inits <- init_joint_cov(n.chain,count_all,cov,initial_values)
    out <- rstan::sampling(c(stanmodels$joint_binary_cov_pois,
                             stanmodels$joint_binary_cov_gamma)[model_index][[1]],
                           data = rlist::list.append(
                             data,
                             nsitecov = length(cov)+1,
                             mat_site = as.matrix(site_mat)
                           ),
                           cores = cores,
                           #' @srrstats {G2.4,G2.4a} explicit conversion to
                           #' integers for sampling arguments
                           chains = as.integer(n.chain),
                           thin = as.integer(thin),
                           warmup = as.integer(n.iter.burn),
                           iter = (
                             as.integer(n.iter.burn) + as.integer(n.iter.sample)
                           ),
                           init = inits,
                           refresh = ifelse(verbose==TRUE,500,0)
    )
  }

  # assert that the log likelihood is a double
  #' @srrstats {G5.3} assert that model run worked and the log likelihood is
  #' valid (i.e., not NA)
  stopifnot(is.double(sum(colMeans(rstan::extract(out,par='log_lik')$log_lik))))

  # Create a list to store the results
  #' @srrstats {BS5.0} function returns initial values used in computation
  result_list <- list(model = out, inits = inits)

  #' @srrstats {BS5.5} the `model` return object is of class `stanfit`, which
  #' includes information about convergence
  return(result_list)
}

###########
#helper functions: initial values
###########
#' @srrstats {BS2.7,BS2.11} Option for user to provide initial values for each
#' chain

init_joint_cov <- function(n.chain,count_all,cov,initial_values){
  #helper function
  #joint model, catchability coefficient, site covariates
  A <- list()
  if(all(initial_values != 'None')){
    for(i in 1:n.chain){
      A[[i]] <- list(
        if('mu' %in% names(initial_values[[i]])){
          mu <- inits[[i]]$mu
        } else {
          mu <- stats::runif(length(unique(count_all$L)), 0.01, 5)
        },

        if('p10' %in% names(initial_values[[i]])){
          p10 <- inits[[i]]$p10
        } else {
          p10 <- stats::runif(1,log(0.0001),log(0.08))
        },

        if('alpha' %in% names(initial_values[[i]])){
          alpha <- inits[[i]]$alpha
        } else {
          alpha <- rep(0.1,length(cov)+1)
        }
      )
      names(A[[i]]) <- c('mu','p10','alpha')
    }
  } else {
    for(i in 1:n.chain){
      A[[i]] <- list(
        mu <- stats::runif(length(unique(count_all$L)), 0.01, 5),
        p10 <- stats::runif(1,log(0.0001),log(0.08)),
        alpha <- rep(0.1,length(cov)+1)
      )
    }
  }

  return(A)
}

init_joint_cov_catchability <- function(n.chain,count_all,q_names,cov,
                                        initial_values){
  #helper function
  #joint model, catchability coefficient, site covariates
  A <- list()
  if(all(initial_values != 'None')){
    for(i in 1:n.chain){
      A[[i]] <- list(
        if('mu' %in% names(initial_values[[i]])){
          mu <- inits[[i]]$mu
        } else {
          mu <- stats::runif(length(unique(count_all$L)), 0.01, 5)
        },

        if('q' %in% names(initial_values[[i]])){
          q <- as.data.frame(inits[[i]]$q)
        } else {
          q <- as.data.frame(stats::runif(length(q_names),0.01,1))
        },

        if('p10' %in% names(initial_values[[i]])){
          p10 <- inits[[i]]$p10
        } else {
          p10 <- stats::runif(1,log(0.0001),log(0.08))
        },

        if('alpha' %in% names(initial_values[[i]])){
          alpha <- inits[[i]]$alpha
        } else {
          alpha <- rep(0.1,length(cov)+1)
        }
      )
      names(A[[i]]) <- c('mu','q','p10','alpha')
    }
  } else {
    for(i in 1:n.chain){
      A[[i]] <- list(
        mu <- stats::runif(length(unique(count_all$L)), 0.01, 5),
        q <- as.data.frame(stats::runif(length(q_names),0.01,1)),
        p10 <- stats::runif(1,log(0.0001),log(0.08)),
        alpha <- rep(0.1,length(cov)+1)
      )
    }
  }

  return(A)
}

init_joint_catchability <- function(n.chain,count_all,q_names,initial_values){
  #helper function
  #joint model, catchability coefficient, no site covariates
  A <- list()
  if(all(initial_values != 'None')){
    for(i in 1:n.chain){
      A[[i]] <- list(
        if('mu' %in% names(initial_values[[i]])){
          mu <- inits[[i]]$mu
        } else {
          mu <- stats::runif(length(unique(count_all$L)), 0.01, 5)
        },

        if('q' %in% names(initial_values[[i]])){
          q <- as.data.frame(inits[[i]]$q)
        } else {
          q <- as.data.frame(stats::runif(length(q_names),0.01,1))
        },

        if('p10' %in% names(initial_values[[i]])){
          p10 <- inits[[i]]$p10
        } else {
          p10 <- stats::runif(1,log(0.0001),log(0.08))
        },

        if('beta' %in% names(initial_values[[i]])){
          beta <- inits[[i]]$beta
        } else {
          beta <- 0.5
        }
      )
      names(A[[i]]) <- c('mu','q','p10','beta')
    }
  } else {
    for(i in 1:n.chain){
      A[[i]] <- list(
        mu <- stats::runif(length(unique(count_all$L)), 0.01, 5),
        q <- as.data.frame(stats::runif(length(q_names),0.01,1)),
        p10 <- stats::runif(1,log(0.0001),log(0.08)),
        beta <- 0.5
      )
    }
  }

  return(A)
}

init_joint <- function(n.chain,count_all,initial_values){
  #helper function
  #joint model, no catchability coefficient, no site covariates
  A <- list()
  if(all(initial_values != 'None')){
    for(i in 1:n.chain){
      A[[i]] <- list(
        if('mu' %in% names(initial_values[[i]])){
          mu <- inits[[i]]$mu
        } else {
          mu <- stats::runif(length(unique(count_all$L)), 0.01, 5)
        },

        if('p10' %in% names(initial_values[[i]])){
          p10 <- inits[[i]]$p10
        } else {
          p10 <- stats::runif(1,log(0.0001),log(0.08))
        },

        if('beta' %in% names(initial_values[[i]])){
          beta <- inits[[i]]$beta
        } else {
          beta <- 0.5
        }
      )
      names(A[[i]]) <- c('mu','p10','beta')
    }
  } else {
    for(i in 1:n.chain){
      A[[i]] <- list(
        mu <- stats::runif(length(unique(count_all$L)), 0.01, 5),
        p10 <- stats::runif(1,log(0.0001),log(0.08)),
        beta <- 0.5
      )
    }
  }

  return(A)
}

################
#helper functions: input checks
################
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#'  messages

# input checks if catchabilty coefficients are used
catchability_checks <- function(data,cov){

  ## make sure count.type is not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for zero-length
  #' data
  if (dim(data$count.type)[1]==0) {
    errMsg <- "count.type contains zero-length data."
    stop(errMsg)
  }
  ## make sure no column is entirely NA in count.type
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with
  #' all NA
  if (any(apply(data$count.type, 2, function(col) all(is.na(col))))) {
    errMsg <- "count.type contains a column with all NA."
    stop(errMsg)
  }

  #' @srrstats {G2.13} Pre-processing routines to check for missing data
  ## All tags in data are valid (i.e., include qPCR.N, qPCR.K, count,
  ## count.type, and site.cov)
  #cov='None'
  if (all(cov=='None') && !all(c('qPCR.N', 'qPCR.K',
                                 'count','count.type') %in% names(data))){
    errMsg <- paste0("Data should include 'qPCR.N', 'qPCR.K', ",
                     "'count', and 'count.type'.")
    stop(errMsg)
  }
  #q=TRUE and cov!='None'
  if (all(cov!='None') && !all(c('qPCR.N', 'qPCR.K', 'count',
                                 'count.type','site.cov') %in% names(data))){
    errMsg <- paste0("Data should include 'qPCR.N', 'qPCR.K', ",
                     "'count', 'count.type', and 'site.cov'.")
    stop(errMsg)
  }

  ## make sure dimensions of count and count.type are equal, if
  ## count.type is present
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #' is dimensionally commensurate
  if(dim(data$count)[1]!=dim(data$count.type)[1]|
     dim(data$count)[2]!=dim(data$count.type)[2]) {
    errMsg <- "Dimensions of count and count.type do not match."
    stop(errMsg)
  }

  ## make sure all data is numeric -- if q == TRUE
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #' unsupported type
  if(is.numeric(data$qPCR.K)==FALSE |
     is.numeric(data$qPCR.N)==FALSE |
     is.numeric(data$count)==FALSE |
     is.numeric(data$count.type)==FALSE) {
    errMsg <- "Data should be numeric."
    stop(errMsg)
  }
  ## make sure locations of NAs in count data match locations of NAs in
  ## count.type data
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input
  #' data is dimensionally commensurate
  if(sum(is.na(data$count.type))!=sum(is.na(data$count))){
    errMsg <- paste0("Empty data cells (NA) in count data should match ",
                     "empty data cells (NA) in count.type data.")
    stop(errMsg)
  }
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input
  #' data is dimensionally commensurate
  if(any((which(is.na(data$count))==which(is.na(data$count.type)))==FALSE)){
    errMsg <- paste0("Empty data cells (NA) in count data should match ",
                     "empty data cells (NA) in count.type data.")
    stop(errMsg)
  }
  ## the smallest count.type is 1
  if(min(data$count.type,na.rm=TRUE) != 1){
    errMsg <- paste0("The first gear type should be referenced as 1 in ",
                     "count.type. Subsequent gear types should be referenced ",
                     "2, 3, 4, etc.")
    stop(errMsg)
  }

  ## count.type are integers
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #' unsupported type
  if(!all(data$count.type %% 1 %in% c(0,NA))){
    errMsg <- "All values in count.type should be integers."
    stop(errMsg)
  }
}

# input checks if no catchabilty coefficients are used
no_catchability_checks <- function(data,cov){

  #' @srrstats {G2.13} Pre-processing routines to check for missing data
  ## All tags in data are valid (i.e., include qPCR.N, qPCR.K, count,
  ## and site.cov)
  #cov='None'
  if (all(cov=='None') &&
      !all(c('qPCR.N', 'qPCR.K', 'count') %in% names(data))){
    errMsg <- "Data should include 'qPCR.N', 'qPCR.K', and 'count'."
    stop(errMsg)
  }
  #cov!='None'
  if (all(cov!='None') &&
      !all(c('qPCR.N', 'qPCR.K', 'count','site.cov') %in% names(data))){
    errMsg <- "Data should include 'qPCR.N', 'qPCR.K', 'count', and 'site.cov'."
    stop(errMsg)
  }

  ## make sure all data is numeric -- if q == FALSE
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'  unsupported type
  if(is.numeric(data$qPCR.K)==FALSE |
     is.numeric(data$qPCR.N)==FALSE |
     is.numeric(data$count)==FALSE ) {
    errMsg <- "Data should be numeric."
    stop(errMsg)
  }
}

# input checks for all variations
all_checks <- function(data,cov,family,p10priors,phipriors){


  ## make sure count, qPCR.N, and qPCR.K are not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for
  #' zero-length data
  if (dim(data$qPCR.N)[1]==0 | dim(data$qPCR.K)[1]==0 |
      dim(data$count)[1]==0) {
    errMsg <- "Input data contains zero-length data."
    stop(errMsg)
  }
  ## make sure no column is entirely NA in qPCR.N
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column
  #' with all NA
  if (any(apply(data$qPCR.N, 2, function(col) all(is.na(col))))) {
    errMsg <- "qPCR.N contains a column with all NA."
    stop(errMsg)
  }

  ## make sure no column is entirely NA in qPCR.K
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with
  #' all NA
  if (any(apply(data$qPCR.K, 2, function(col) all(is.na(col))))) {
    errMsg <- "qPCR.K contains a column with all NA."
    stop(errMsg)
  }

  ## make sure no column is entirely NA in count
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with
  #' all NA
  if (any(apply(data$count, 2, function(col) all(is.na(col))))) {
    errMsg <- "count contains a column with all NA."
    stop(errMsg)
  }

  ## make sure dimensions of qPCR.N and qPCR.K are equal
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #' is dimensionally commensurate
  if (dim(data$qPCR.N)[1]!=dim(data$qPCR.K)[1]|
      dim(data$qPCR.N)[2]!=dim(data$qPCR.K)[2]) {
    errMsg <- "Dimensions of qPCR.N and qPCR.K do not match."
    stop(errMsg)
  }
  ## make sure number of rows in count = number of rows in qPCR.N and qPCR.K
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #' is dimensionally commensurate
  if (dim(data$qPCR.N)[1]!=dim(data$count)[1]) {
    errMsg <- paste0("Number of sites (rows) in qPCR data and traditional ",
                     "survey count data do not match.")
    stop(errMsg)
  }

  ## make sure locations of NAs in qPCR.N data match locations of NAs in
  ## qPCR.K data
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #' is dimensionally commensurate
  if(any((which(is.na(data$qPCR.N))==which(is.na(data$qPCR.K)))==FALSE)){
    errMsg <- paste0("Empty data cells (NA) in qPCR.N data should match ",
                     "empty data cells (NA) in qPCR.K data.")
    stop(errMsg)
  }

  ## make sure family is either 'poisson', 'negbin', or 'gamma'
  #' @srrstats {G2.3,G2.3a,G2.3b} Permit only expected univariate
  #' (case-insensitive) parameter values
  if(!c(tolower(family) %in% c('poisson','negbin','gamma'))){
    errMsg <- paste0("Invalid family. Options include 'poisson', 'negbin', ",
                     "and 'gamma'.")
    stop(errMsg)
  }

  ## p10priors is a vector of two integers
  #' @srrstats {G2.0,BS2.2,BS2.3,BS2.4,BS2.5,BS2.6} Checks of vector length and
  #' appropriateness of distributional parameters (i.e., vector of length 2,
  #' numeric values > 0), implemented prior to analytic routines
  if(!is.numeric(p10priors) | length(p10priors)!=2 | any(p10priors<=0)){
    errMsg <- paste0("p10priors should be a vector of two positive numeric ",
                     "values. ex. c(1,20)")
    stop(errMsg)
  }

  ## phipriors is a vector of two numeric values
  #' @srrstats {G2.0,BS2.2,BS2.3,BS2.4,BS2.5,BS2.6} Checks of vector length
  #' and appropriateness of distributional parameters (i.e., vector of length
  #' 2, numeric values > 0), implemented prior to analytic routines
  if(!is.numeric(phipriors) | length(phipriors)!=2 | any(phipriors<=0)){
    errMsg <- paste0("phipriors should be a vector of two positive numeric ",
                     "values. ex. c(0.25,0.25)")
    stop(errMsg)
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  if(any(data$count==Inf,na.rm=TRUE) | any(data$count==-Inf,na.rm=TRUE)){
    errMsg <- "count contains undefined values (i.e., Inf or -Inf)"
    stop(errMsg)
  }

  ## count are integers, if family is poisson or negbin
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #' for distributional parameters (i.e., count data must be an integer if a
  #' poisson or negative binomial distribution is used), implemented prior to
  #' analytic routines
  if(tolower(family) %in% c('poisson','negbin')){
    if(!all(data$count %% 1 %in% c(0,NA)) | any(data$count < 0,na.rm=TRUE)){
      errMsg <- paste0("All values in count should be non-negative integers. ",
                       "Use family = 'gamma' if count is continuous.")
      stop(errMsg)
    }
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  if(any(data$qPCR.N==Inf,na.rm=TRUE) | any(data$qPCR.N==-Inf,na.rm=TRUE)){
    errMsg <- "qPCR.N contains undefined values (i.e., Inf or -Inf)"
    stop(errMsg)
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  if(any(data$qPCR.K==Inf,na.rm=TRUE) | any(data$qPCR.K==-Inf,na.rm=TRUE)){
    errMsg <- "qPCR.K contains undefined values (i.e., Inf or -Inf)"
    stop(errMsg)
  }

  ## qPCR.N are integers
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #' for distributional parameters (i.e., qPCR data are non-negative integers),
  #' implemented prior to analytic routines
  if(!all(data$qPCR.N %% 1 %in% c(0,NA)) | any(data$qPCR.N < 0,na.rm=TRUE)){
    errMsg <- "All values in qPCR.N should be non-negative integers."
    stop(errMsg)
  }

  ## qPCR.K are integers
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #' for distributional parameters (i.e., qPCR data are non-negative integers),
  #' implemented prior to analytic routines
  if(!all(data$qPCR.K %% 1 %in% c(0,NA)) | any(data$qPCR.K < 0,na.rm=TRUE)){
    errMsg <- "All values in qPCR.K should be non-negative integers."
    stop(errMsg)
  }



}

# input checks if site-level covariates are used
covariate_checks <- function(data,cov){

  ## make sure site.cov is not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for zero-length
  #' data
  if (dim(data$site.cov)[1]==0) {
    errMsg <- "site.cov contains zero-length data."
    stop(errMsg)
  }
  ## make sure no column is entirely NA in site.cov
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with
  #' all NA
  if (any(apply(data$site.cov, 2, function(col) all(is.na(col))))) {
    errMsg <- "site.cov contains a column with all NA."
    stop(errMsg)
  }

  ## site.cov is numeric, if present
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #' unsupported type
  if(!is.numeric(data$site.cov)){
    errMsg <- "site.cov should be numeric."
    stop(errMsg)
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  if(any(data$site.cov==Inf) | any(data$site.cov==-Inf)){
    errMsg <- "site.cov contains undefined values (i.e., Inf or -Inf)"
    stop(errMsg)
  }

  ## cov values match column names in site.cov
  if(!all(cov %in% colnames(data$site.cov))){
    errMsg <- paste0("cov values should be listed in the column names of ",
                     "site.cov in the data.")
    stop(errMsg)
  }

  ## site.cov has same number of rows as qPCR.N and count, if present
  #' @srrstats {BS2.1} Pre-processing routines to ensure all input data is
  #' dimensionally commensurate
  if(dim(data$qPCR.N)[1]!=dim(data$site.cov)[1]){
    errMsg <- paste0("The number of rows in site.cov matrix should match the ",
                     "number of rows in all other matrices.")
    stop(errMsg)
  }

  ## add warning if number of covariates is greater than the number of sites
  #' @srrstats {G5.8d} Pre-processing routines to check if data is outside
  #' scope of algorithm (i.e., # site-level covariates is greater than the
  #' number of sites)
  if(length(cov)>dim(data$site.cov)[2]){
    warnMsg = paste0("The number of site-level covariates exceeds the number ",
                     "of sites (i.e., n < p).")
    warning(warnMsg)
  }

  ## add warning if number of site-covariate data has perfect collinearity
  #' @srrstats {BS3.1} Pre-processing routines to check if site covariate
  #' data has perfect collinearity
  rank_mat <- qr(data$site.cov)$rank
  if(rank_mat < ncol(data$site.cov)){
    warnMsg = "Data in site.cov exhibits perfect collinearity."
    warning(warnMsg)
  }
}

# checks if initial values are provided
initial_values_checks <- function(initial_values,data,cov,n.chain){

  ## length of initial values is equal to the number of chains
  if(length(initial_values)!=n.chain){
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
      if(length(initial_values[[i]]$mu)!=dim(data$count)[1]){
        errMsg <- paste0("The length of initial values for 'mu' should equal ",
                         "the number of sites.")
        stop(errMsg)
      }
    }

    ## check p10 input
    if('p10' %in% names(initial_values[[i]])){
      ## if p10 is numeric
      if(!is.numeric(initial_values[[i]]$p10)){
        errMsg <- "Initial values for 'p10' should be numeric."
        stop(errMsg)
      }
      ## check p10 length
      if(length(initial_values[[i]]$p10)!=1){
        errMsg <- "The length of initial values for 'p10' should equal 1."
        stop(errMsg)
      }
    }

    ## check beta input
    if('beta' %in% names(initial_values[[i]])){
      ## if beta is numeric
      if(!is.numeric(initial_values[[i]]$beta)){
        errMsg <- "Initial values for 'beta' should be numeric."
        stop(errMsg)
      }
      ## check beta length
      if(length(initial_values[[i]]$beta)!=1){
        errMsg <- "The length of initial values for 'beta' should equal 1."
        stop(errMsg)
      }
    }

    ## check alpha input
    if('alpha' %in% names(initial_values[[i]])){
      ## if alpha is numeric
      if(any(!is.numeric(initial_values[[i]]$alpha)) |
         any(initial_values[[i]]$alpha < 0)){
        errMsg <- "Initial values for 'alpha' should be numeric."
        stop(errMsg)
      }
      ## check alpha length
      if(length(initial_values[[i]]$alpha)!=(length(cov)+1)){
        errMsg <- paste0("The length of initial values for 'alpha' should ",
                         "equal: # covariates + 1 (i.e., including intercept).")
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
      if(length(initial_values[[i]]$q)!=(length(table(data$count.type))-1)){
        errMsg <- paste0("The length of initial values for 'q' should equal:",
                         " # unique gear types - 1 (i.e., q for reference ",
                         "type = 1).")
        stop(errMsg)
      }
    }


  }
}


