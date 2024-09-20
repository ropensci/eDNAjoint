#' Specify and fit joint model using count data from traditional surveys and
#' eDNA qPCR data
#'
#' This function implements a Bayesian model that integrates data from paired
#' eDNA and traditional surveys, as described in Keller et al. (2022)
#' <https://doi.org/10.1002/eap.2561>. The model estimates parameters including
#' the expected species catch rate and the probability of false positive eDNA
#' detection. This function allows for optional model variations, like inclusion
#' of site-level covariates that scale the sensitivity of eDNA sampling relative
#' to traditional sampling, as well as estimation of catchability coefficients
#' when multiple traditional gear types are used. Model is implemented using
#' Bayesian inference using the `rstan` package, which uses Hamiltonian Monte
#' Carlo to simulate the posterior distributions. See more examples in the
#' \href{https://bookdown.org/abigailkeller/eDNAjoint_vignette/}{Package
#' Vignette}.
#'
#' @srrstats {G1.0,G1.1} This software makes available a novel algorithm/model
#'   that was previously published in the literature. The literature reference
#'   for the joint model is provided here.
#'
#' @srrstats {G1.4,G1.3} Roxygen function documentation begins here, with
#'   definitions of statistical terminology and inputs
#' @export
#' @srrstats {BS1.1,BS3.0,G2.14} Descriptions of how to enter data, description
#'   of how NAs are handled (which are informative and should be deliberately
#'   included in input data, where necessary)
#' @param data A list containing data necessary for model fitting. Valid tags
#'   are `qPCR.N`, `qPCR.K`, `count`, `count.type`, and `site.cov`. `qPCR.N` and
#'   `qPCR.K` are matrices or data frames with first dimension equal to the
#'   number of sites (i) and second dimension equal to the maximum number of
#'   eDNA samples at a given site (m). `qPCR.N` contains the total number of
#'   qPCR replicates per site and eDNA sample, and `qPCR.K` contains the total
#'   number of positive qPCR detections per site and eDNA sample. `count` is a
#'   matrix or data frame of traditional survey count data, with first
#'   dimension equal to the number of sites (i) and second dimension equal to
#'   the maximum number of traditional survey replicates at a given site (j).
#'   `count.type` is an optional matrix or data frame of integers indicating
#'   gear type used in corresponding count data, with first dimension equal to
#'   the number of sites (i) and second dimension equal to the maximum number
#'   of traditional survey replicates at a given site. Values should be
#'   integers beginning with 1 (referring to the first gear type) to n (last
#'   gear type). `site.cov` is an optional matrix or data frame of site-level
#'   covariate data, with first dimension equal to the number of sites (i).
#'   `site.cov` should include informative column names. Empty cells should
#'   be NA and will be removed during processing. Sites, i, should be consistent
#'   in all qPCR, count, and site covariate data.
#' @srrstats {G2.1a} Here are explicit documentation of vector input types
#' @param cov A character vector indicating the site-level covariates to include
#'   in the model. Default value is NULL.
#' @param family The distribution class used to model traditional survey count
#'   data. Options include poisson ('poisson'), negative binomial ('negbin'),
#'   and gamma ('gamma'). Default value is 'poisson'.
#' @srrstats {BS1.0,G2.1a,BS1.2} Description of hyperparameters and how to
#'   specify prior distributions, explicit documentation of vector input types
#' @param p10priors A numeric vector indicating beta distribution
#'   hyperparameters (alpha, beta) used as the prior distribution for the eDNA
#'   false positive probability (p10). Default vector is c(1,20).
#' @param q A logical value indicating whether to estimate a catchability
#'   coefficient, q, for traditional survey gear types (TRUE) or to not
#'   estimate a catchability coefficient, q, for traditional survey gear types
#'   (FALSE). Default value is FALSE.
#' @srrstats {BS1.0,G2.1a,BS1.2} Description of hyperparameters and how to
#'   specify prior distributions, explicit documentation of vector input types
#' @param phipriors A numeric vector indicating gamma distribution
#'   hyperparameters (shape, rate) used as the prior distribution for phi, the
#'   overdispersion in the negative binomial distribution for traditional survey
#'   gear data. Used when family = 'negbin.' If family = 'negbin', then
#'   default vector is c(0.25,0.25), otherwise, default is NULL.
#' @param multicore A logical value indicating whether to parallelize chains
#'   with multiple cores. Default is TRUE.
#' @srrstats {BS2.7,BS2.11} Option for user to provide initial values for each
#'   chain
#' @param initial_values A list of lists of initial values to use in MCMC. The
#'   length should equal the number of MCMC chains. Initial values can be
#'   provided for parameters: beta, p10 (log-scale), mu, q, alpha. If no
#'   initial values are provided, default random values are drawn.
#' @srrstats {BS1.3} Description of parameters used in the computational process
#'   begins here
#' @param n.chain Number of MCMC chains. Default value is 4.
#' @param n.iter.burn Number of warm-up MCMC iterations. Default value is 500.
#' @param n.iter.sample Number of sampling MCMC iterations. Default value is
#'   2500.
#' @param thin A positive integer specifying the period for saving samples.
#'   Default value is 1.
#' @param adapt_delta Numeric value between 0 and 1 indicating target average
#'   acceptance probability used in `rstan::sampling`. Default value is 0.9.
#' @srrstats {BS2.12} Parameter controlling the verbosity of output
#' @param verbose Logical value controlling the verbosity of output (i.e.,
#'   warnings, messages, progress bar). Default is TRUE.
#' @param seed A positive integer seed used for random number generation in
#'   MCMC. Default is NULL, which means the seed is generated from 1 to the
#'   maximum integer supported by R.
#' @srrstats {BS5.0} function returns initial values used in computation
#' @return A list of:
#' \itemize{
#' \item a model object of class `stanfit` returned by `rstan::sampling`
#' \item initial values used in MCMC
#' }
#'
#' @srrstats {G2.0a} Explicit secondary documentation of any expectations on
#'   lengths of inputs
#' @note  Before fitting the model, this function checks to ensure that the
#'   model specification is possible given the data files. These checks include:
#' \itemize{
#' \item  All tags in data are valid (i.e., include qPCR.N, qPCR.K, count,
#'   count.type, and site.cov).
#' \item  Dimensions of qPCR.N and qPCR.K are equal, and dimensions of count and
#'   count.type are equal (if count.type provided).
#' \item  Number of sites in qPCR and count data are equal.
#' \item  All data are numeric (i.e., integer or NA).
#' \item  Empty data cells (NA) match in qPCR.N and qPCR.K and in count and
#'   count.type.
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
#' # Ex. 1: Implementing the joint model
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
#' # Fit a basic model with paired eDNA and traditional survey data.
#' # Count data is modeled using a poisson distribution.
#' fit <- jointModel(data = gobyData, family = "poisson", p10priors = c(1,20),
#'                   multicore = FALSE)
#'
#' # Ex. 2: Implementing the joint model with site-level covariates
#'
#' # With the same data, fit a model including 'Filter_time' and 'Salinity'
#' # site-level covariates
#' # These covariates will scale the sensitivity of eDNA sampling relative to
#' # traditional surveys
#' # Count data is modeled using a poisson distribution.
#' fit.cov <- jointModel(data = gobyData, cov = c('Filter_time','Salinity'),
#'                       family = "poisson", p10priors = c(1,20),
#'                       multicore = FALSE)
#'
#'
#' # Ex. 3: Implementing the joint model with multiple traditional gear types
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
#' fit.q <- jointModel(data = greencrabData, cov = NULL, family = "negbin",
#'                     p10priors = c(1,20), q = TRUE, phipriors = c(0.25,0.25),
#'                     multicore = FALSE, initial_values = NULL,
#'                     n.chain = 4, n.iter.burn = 500,
#'                     n.iter.sample = 2500, thin = 1, adapt_delta = 0.9,
#'                     verbose = TRUE, seed = 123)
#' }
#'

jointModel <- function(data, cov = NULL, family = 'poisson',
                       p10priors = c(1,20),
                       q = FALSE, phipriors = NULL, multicore = TRUE,
                       initial_values = NULL, n.chain = 4, n.iter.burn = 500,
                       n.iter.sample = 2500, thin = 1, adapt_delta = 0.9,
                       verbose = TRUE, seed = NULL) {


  # make character inputs case-insensitive
  #' @srrstats {G2.3b} Allow case-insensitive character parameter values
  family <- tolower(family)


  # get phipriors
  if(family != 'negbin'){
    phipriors <- NULL
  } else if(family == 'negbin' && is.null(phipriors)){
    phipriors <- c(0.25,0.25)
  } else if(family == 'negbin' && !is.null(phipriors)){
    phipriors <- phipriors
  }

  ####
  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using the following
  #'   helper functions

  if (q == TRUE) {
    # model with catchability coefficients
    catchability_checks(data,cov)
  } else {
    # model without catchability coefficients
    no_catchability_checks(data,cov)
  }

  # model with covariates
  if (all(!is.null(cov))) {
    covariate_checks(data,cov)
  }

  # all models
  all_checks(data,cov,family,p10priors,phipriors,n.chain,n.iter.burn,
             n.iter.sample,thin,adapt_delta,seed)

  # initial value checks
  if(all(!is.null(initial_values))){
    initial_values_checks(initial_values,data,cov,n.chain)
  }

  if (!requireNamespace("rstan", quietly = TRUE)){
    stop ("The 'rstan' package is not installed.", call. = FALSE)
  }

  ###
  #convert data to long format

  '%>%' <- magrittr::`%>%`

  # convert qPCR data to long format
  #' @srrstats {G2.7} Use as.data.frame() to allow input list of any tabular
  #'   form (i.e., matrix, etc.)
  qPCR_all <- as.data.frame(data$qPCR.N) %>%
    dplyr::mutate(L = 1:dim(data$qPCR.N)[1]) %>%
    tidyr::pivot_longer(cols =! L,values_to = 'N') %>%
    #' @srrstats {G2.15} Software does not assume non-missingness and actually
    #'   expects it if the number of observations across sites is unequal
    tidyr::drop_na()
  qPCR.K_df <- as.data.frame(data$qPCR.K) %>%
    dplyr::mutate(L = 1:dim(data$qPCR.K)[1]) %>%
    tidyr::pivot_longer(cols =! L,values_to = 'K') %>%
    #' @srrstats {G2.15} Software does not assume non-missingness and actually
    #'   expects it if the number of observations across sites is unequal
    tidyr::drop_na()
  qPCR_all$K <- qPCR.K_df$K

  # convert count data to long format
  #' @srrstats {G2.7} Use as.data.frame() to allow input list of any tabular
  #'   form (i.e., matrix, etc.)
  count_all <- as.data.frame(data$count) %>%
    dplyr::mutate(L = 1:dim(data$count)[1]) %>%
    tidyr::pivot_longer(cols =! L,values_to = 'count') %>%
    #' @srrstats {G2.15} Software does not assume non-missingness and actually
    #' expects it if the number of observations across sites is unequal
    tidyr::drop_na()

  # subset count data to remove sites without traditional samples
  count_df <- as.data.frame(data$count)
  sub_count <- count_df[rowSums(is.na(count_df)) != ncol(count_df), ]

  # add site index to count data
  index_match <- as.data.frame(cbind(unique(count_all$L),
                                     1:dim(sub_count)[1]))
  colnames(index_match) <- c('L','R')
  count_all <- dplyr::left_join(count_all,index_match,by='L')

  # get site indices with paired and unpaired dna samples
  trad_ind <- index_match$L
  dna_ind <- unique(qPCR_all$L)[!unique(qPCR_all$L)%in% trad_ind]

  # subset qPCR data -- paired
  qPCR_all_trad <- qPCR_all[qPCR_all$L %in% trad_ind,]
  L_match_trad <- as.data.frame(cbind(unique(qPCR_all_trad$L),
                                      1:length(unique(qPCR_all_trad$L))))
  colnames(L_match_trad) <- c('L','L_unique')
  qPCR_all_trad <- dplyr::left_join(qPCR_all_trad,L_match_trad,by='L')

  # subset qPCR data -- unpaired
  if(length(dna_ind)>0){
    qPCR_all_dna <- qPCR_all[qPCR_all$L %in% dna_ind,]
    L_match_dna <- as.data.frame(cbind(unique(qPCR_all_dna$L),
                                       1:length(unique(qPCR_all_dna$L))))
    colnames(L_match_dna) <- c('L','L_unique')
    qPCR_all_dna <- dplyr::left_join(qPCR_all_dna,L_match_dna,by='L')
  } else {
    qPCR_all_dna <- as.data.frame(matrix(NA,nrow=0,ncol=4))
    colnames(qPCR_all_dna) <- c('L','N','K','L_unique')
    L_match_dna <- as.data.frame(matrix(NA,nrow=0,ncol=4))
  }

  # if q == TRUE, add count type data to count df
  if(isCatch_type(q)){
    q_ref <- 1
    #' @srrstats {G2.7} Use as.data.frame() to allow input list of any tabular
    #'   form (i.e., matrix, etc.)
    count.type_df <- as.data.frame(data$count.type) %>%
      dplyr::mutate(L = 1:dim(data$count.type)[1]) %>%
      tidyr::pivot_longer(cols =! L,values_to = 'count.type') %>%
      tidyr::drop_na()
    count_all$count.type <- count.type_df$count.type

    #create vector of q coefficient names
    counttypes <- unique(count_all$count.type)
    names <- counttypes[!counttypes == q_ref]
    #' @srrstats {G2.4,G2.4c} Explicit conversion to character
    q_names <- as.character(paste0('q_',names))

    #add dummy variables for count type
    for(i in seq_along(q_names)){
      count_all[,q_names[i]] <- ifelse(count_all$count.type == names[i],1,0)
    }
  }

  # if present, prepare covariate data
  if(all(!is.null(cov))){
    #' @srrstats {G2.7,G2.10} Use as.data.frame() to allow input list of any
    #'   tabular form (i.e., matrix, etc.) and converts before filtering columns
    #'   based on input 'cov'
    site_mat <- as.data.frame(as.data.frame(data$site.cov)[,cov])
    site_mat <- cbind(as.data.frame(rep(1,length(site_mat[,1]))),site_mat)
    colnames(site_mat) <- c('int',cov)
  }

  # convert p10 prior
  # p10 prior: convert beta(1,20) to lognormal distribution
  # moment match from beta(alpha,beta) to normal(mu, sigma^2)
  alpha <- p10priors[1]
  beta <- p10priors[2]
  mu <- alpha/(alpha+beta)
  sigma_2 <- (alpha*beta)/((alpha+beta)^2*(alpha+beta+1))
  # convert normal(mu, sigma^2) to lognormal(mu, sigma^2)
  mu_ln <- log(mu^2/sqrt(mu^2+sigma_2))
  sigma_2_ln <- log(1+sigma_2/mu^2)
  sigma_ln <- sqrt(sigma_2_ln)

  # create data that will be present in all model variations
  model_data <- list(
    S = nrow(qPCR_all_trad),
    S_dna = nrow(qPCR_all_dna),
    C = nrow(count_all),
    L = qPCR_all_trad$L_unique,
    L_dna = qPCR_all_dna$L_unique,
    R = count_all$R,
    Nloc_dna = length(unique(qPCR_all_dna$L)),
    Nloc_trad = length(unique(qPCR_all_trad$L)),
    trad_ind = trad_ind,
    dna_ind = as.array(dna_ind),
    E = count_all$count,
    N = qPCR_all_trad$N,
    K = qPCR_all_trad$K,
    N_dna = qPCR_all_dna$N,
    K_dna = qPCR_all_dna$K,
    p10priors = c(mu_ln,sigma_ln),
    control = list(adapt_delta = adapt_delta)
  )

  # append data if q == TRUE
  if(isCatch_type(q)){
    model_data <- rlist::list.append(
      model_data,
      nparams = length(q_names),
      mat = as.matrix(count_all[,q_names])
    )
  }
  # append data if family == negbin
  if(isNegbin_type(family)){
    model_data <- rlist::list.append(
      model_data,
      phipriors = phipriors
    )
  }
  # append data if cov != NULL
  if(isCov_type(cov)){
    model_data <- rlist::list.append(
      model_data,
      nsitecov = length(cov)+1,
      mat_site = as.matrix(site_mat)
    )
  }

  # set up core number
  cores <- ifelse(multicore == TRUE,
                  parallel::detectCores(),
                  1)

  # get seed
  SEED <- ifelse(!is.null(seed),
                 as.integer(seed),
                 sample.int(.Machine$integer.max, 1))

  # get stan model
  model_index <- get_stan_model(q, family, cov)

  # get initial values
  if(isCatch_type(q)){
    inits <- get_inits(n.chain,qPCR_all,initial_values,cov,L_match_trad,
                       L_match_dna,data,q_names)
  } else {
    inits <- get_inits(n.chain,qPCR_all,initial_values,cov,L_match_trad,
                       L_match_dna,data)
  }

  # run model
  out <- rstan::sampling(
    c(stanmodels$joint_binary_cov_catchability_pois,
      stanmodels$joint_binary_cov_catchability_negbin,
      stanmodels$joint_binary_cov_catchability_gamma,
      stanmodels$joint_binary_catchability_pois,
      stanmodels$joint_binary_catchability_negbin,
      stanmodels$joint_binary_catchability_gamma,
      stanmodels$joint_binary_cov_pois,
      stanmodels$joint_binary_cov_negbin,
      stanmodels$joint_binary_cov_gamma,
      stanmodels$joint_binary_pois,
      stanmodels$joint_binary_negbin,
      stanmodels$joint_binary_gamma)[model_index][[1]],
    data = model_data,
    cores = cores,
    seed = SEED,
    #' @srrstats {G2.4,G2.4a} explicit conversion to
    #'   integers for sampling arguments
    chains = as.integer(n.chain),
    thin = as.integer(thin),
    warmup = as.integer(n.iter.burn),
    iter = (
      as.integer(n.iter.burn) + as.integer(n.iter.sample)
    ),
    init = inits,
    refresh = ifelse(verbose == TRUE,500,0)
  )

  # assert that the log likelihood is a double
  #' @srrstats {G5.3} assert that model run worked and the log likelihood is
  #'   valid (i.e., not NA)
  stopifnot(is.double(sum(colMeans(rstan::extract(out,
                                                  par = 'log_lik')$log_lik))))

  # check for divergent transitions
  div_trans <- sum(lapply(rstan::get_sampler_params(out,
                                                    inc_warmup = FALSE),
                          div_check)[[1]])
  # print either troubleshooting or visualization tips
  if(div_trans>0){
    url <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette/',
                  'tips.html#troubleshooting-tips')
    message <- 'Refer to the eDNAjoint guide for troubleshooting tips: '
    cat(message, url, "\n")
  } else {
    url <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette/',
                  'tips.html#visualization-tips')
    message <- 'Refer to the eDNAjoint guide for visualization tips: '
    cat(message, url, "\n")
  }

  # add chain names to init list
  names(inits) <- paste0('chain',seq(1,n.chain,1))

  # Create a list to store the results
  #' @srrstats {BS5.0} function returns initial values used in computation
  result_list <- list(model = out, inits = inits)

  #' @srrstats {BS5.5} the `model` return object is of class `stanfit`, which
  #'   includes information about convergence
  return(result_list)
}


####################################
# helper functions: get model type #
####################################

isCatch_type <- function(q){
  out <- ifelse(q == TRUE,TRUE,FALSE)
  return(out)
}
isCov_type <- function(cov){
  out <- ifelse(all(!is.null(cov)),TRUE,FALSE)
  return(out)
}
isNegbin_type <- function(family){
  out <- ifelse(family == 'negbin',TRUE,FALSE)
  return(out)
}
get_family_index <- function(family){
  if(family == 'poisson'){
    index <- 1
  } else if(family == 'negbin'){
    index <- 2
  } else if(family == 'gamma'){
    index <- 3
  }
  return(index)
}

####################################
# helper functions: get stan model #
####################################

get_stan_model <- function(q, family, cov){

  index <- get_family_index(family)

  if(isCatch_type(q)){
    if(isCov_type(cov)){
      final_index <- index
    } else {
      final_index <- index + 3
    }
  } else {
    if(isCov_type(cov)){
      final_index <- index + 6
    } else {
      final_index <- index + 9
    }
  }
  return(final_index)
}


####################################
# helper functions: initial values #
####################################
#' @srrstats {BS2.7,BS2.11} Option for user to provide initial values for each
#'   chain


get_inits <- function(n.chain,qPCR_all,initial_values,cov,L_match_trad,
                      L_match_dna,data,q_names=NULL){
  if(!is.null(cov)){
    if(!is.null(q_names)){
      inits <- init_joint_cov_catchability(n.chain,qPCR_all,q_names,cov,
                                           initial_values,L_match_trad,
                                           L_match_dna,data)
    } else {
      inits <- init_joint_cov(n.chain,qPCR_all,cov,initial_values,
                              L_match_trad,L_match_dna,data)
    }
  } else {
    if(!is.null(q_names)){
      inits <- init_joint_catchability(n.chain,qPCR_all,q_names,
                                       initial_values,L_match_trad,
                                       L_match_dna,data)
    } else {
      inits <- init_joint(n.chain,qPCR_all,initial_values,L_match_trad,
                          L_match_dna,data)
    }
  }
  return(inits)
}

init_joint_cov <- function(n.chain,qPCR_all,cov,initial_values,
                           L_match_trad,L_match_dna,data){

  # get mu means
  mu_means_trad <- as.vector(na.omit(rowMeans(data$count,na.rm=TRUE)+0.01))
  mu_means_all <- rep(NA,dim(L_match_dna)[1]+dim(L_match_trad)[1])
  mu_means_all[L_match_trad$L] <- mu_means_trad
  if(dim(L_match_dna)[1]>0){
    mu_means_all[L_match_dna$L] <- rep(mean(mu_means_trad),dim(L_match_dna)[1])
  }

  # helper function
  # joint model, catchability coefficient, site covariates
  A <- list()
  if(all(!is.null(initial_values))){
    for(i in 1:n.chain){
      A[[i]] <- list(
        if('mu' %in% names(initial_values[[i]])){
          mu_trad <- initial_values[[i]]$mu[L_match_trad$L]
        } else {
          mu_trad <- mu_means_trad
        },

        if('mu' %in% names(initial_values[[i]])){
          mu <- initial_values[[i]]$mu
        } else {
          mu <- mu_means_all
        },

        if('p10' %in% names(initial_values[[i]])){
          log_p10 <- log(initial_values[[i]]$p10)
        } else {
          log_p10 <- stats::runif(1,log(0.0001),log(0.08))
        },

        if('alpha' %in% names(initial_values[[i]])){
          alpha <- initial_values[[i]]$alpha
        } else {
          alpha <- rep(0.1,length(cov)+1)
        }
      )
      names(A[[i]]) <- c('mu_trad','mu','log_p10','alpha')
    }
  } else {
    for(i in 1:n.chain){
      A[[i]] <- list(
        mu_trad <- mu_means_trad,
        mu <- mu_means_all,
        log_p10 <- stats::runif(1,log(0.0001),log(0.08)),
        alpha <- rep(0.1,length(cov)+1)
      )
      names(A[[i]]) <- c('mu_trad','mu','log_p10','alpha')
    }
  }

  return(A)
}

init_joint_cov_catchability <- function(n.chain,qPCR_all,q_names,cov,
                                        initial_values,L_match_trad,
                                        L_match_dna,data){

  # get mu means
  mu_means_trad <- as.vector(na.omit(rowMeans(data$count,na.rm=TRUE)+0.01))
  mu_means_all <- rep(NA,dim(L_match_dna)[1]+dim(L_match_trad)[1])
  mu_means_all[L_match_trad$L] <- mu_means_trad
  if(dim(L_match_dna)[1]>0){
    mu_means_all[L_match_dna$L] <- rep(mean(mu_means_trad),dim(L_match_dna)[1])
  }

  # helper function
  # joint model, catchability coefficient, site covariates
  A <- list()
  if(all(!is.null(initial_values))){
    for(i in 1:n.chain){
      A[[i]] <- list(
        if('mu' %in% names(initial_values[[i]])){
          mu_trad_1 <- initial_values[[i]]$mu[L_match_trad$L]
        } else {
          mu_trad_1 <- mu_means_trad
        },

        if('mu' %in% names(initial_values[[i]])){
          mu <- initial_values[[i]]$mu
        } else {
          mu <- mu_means_all
        },

        if('p10' %in% names(initial_values[[i]])){
          log_p10 <- log(initial_values[[i]]$p10)
        } else {
          log_p10 <- stats::runif(1,log(0.0001),log(0.08))
        },

        if('alpha' %in% names(initial_values[[i]])){
          alpha <- initial_values[[i]]$alpha
        } else {
          alpha <- rep(0.1,length(cov)+1)
        },

        if('q' %in% names(initial_values[[i]])){
          q_trans <- as.data.frame(initial_values[[i]]$q)
        } else {
          q_trans <- as.data.frame(stats::runif(length(q_names),0.01,1))
        }
      )
      names(A[[i]]) <- c('mu_trad_1','mu','log_p10','alpha','q_trans')

    }
  } else {
    for(i in 1:n.chain){
      A[[i]] <- list(
        mu_trad_1 <- mu_means_trad,
        mu <- mu_means_all,
        log_p10 <- stats::runif(1,log(0.0001),log(0.08)),
        alpha <- rep(0.1,length(cov)+1),
        q_trans <- as.data.frame(stats::runif(length(q_names),0.01,1))
      )
      names(A[[i]]) <- c('mu_trad_1','mu','log_p10','alpha','q_trans')
    }
  }

  return(A)
}

init_joint_catchability <- function(n.chain,qPCR_all,q_names,initial_values,
                                    L_match_trad,L_match_dna,data){

  # get mu means
  mu_means_trad <- as.vector(na.omit(rowMeans(data$count,na.rm=TRUE)+0.01))
  mu_means_all <- rep(NA,dim(L_match_dna)[1]+dim(L_match_trad)[1])
  mu_means_all[L_match_trad$L] <- mu_means_trad
  if(dim(L_match_dna)[1]>0){
    mu_means_all[L_match_dna$L] <- rep(mean(mu_means_trad),dim(L_match_dna)[1])
  }

  # helper function
  # joint model, catchability coefficient, no site covariates
  A <- list()
  if(all(!is.null(initial_values))){
    for(i in 1:n.chain){
      A[[i]] <- list(
        if('mu' %in% names(initial_values[[i]])){
          mu_trad_1 <- initial_values[[i]]$mu[L_match_trad$L]
        } else {
          mu_trad_1 <- mu_means_trad
        },

        if('mu' %in% names(initial_values[[i]])){
          mu <- initial_values[[i]]$mu
        } else {
          mu <- mu_means_all
        },

        if('p10' %in% names(initial_values[[i]])){
          log_p10 <- log(initial_values[[i]]$p10)
        } else {
          log_p10 <- stats::runif(1,log(0.0001),log(0.08))
        },

        if('beta' %in% names(initial_values[[i]])){
          beta <- initial_values[[i]]$beta
        } else {
          beta <- 0.5
        },

        if('q' %in% names(initial_values[[i]])){
          q_trans <- as.data.frame(initial_values[[i]]$q)
        } else {
          q_trans <- as.data.frame(stats::runif(length(q_names),0.01,1))
        }
      )
      names(A[[i]]) <- c('mu_trad_1','mu','log_p10','beta','q_trans')
    }
  } else {
    for(i in 1:n.chain){
      A[[i]] <- list(
        mu_trad_1 <- mu_means_trad,
        mu <- mu_means_all,
        log_p10 <- stats::runif(1,log(0.0001),log(0.08)),
        beta <- 0.5,
        q_trans <- as.data.frame(stats::runif(length(q_names),0.01,1))
      )
      names(A[[i]]) <- c('mu_trad_1','mu','log_p10','beta','q_trans')
    }
  }

  return(A)
}

init_joint <- function(n.chain,qPCR_all,initial_values,
                       L_match_trad,L_match_dna,data){

  # get mu means
  mu_means_trad <- as.vector(na.omit(rowMeans(data$count,na.rm=TRUE)+0.01))
  mu_means_all <- rep(NA,dim(L_match_dna)[1]+dim(L_match_trad)[1])
  mu_means_all[L_match_trad$L] <- mu_means_trad
  if(dim(L_match_dna)[1]>0){
    mu_means_all[L_match_dna$L] <- rep(mean(mu_means_trad),dim(L_match_dna)[1])
  }

  # helper function
  # joint model, no catchability coefficient, no site covariates
  A <- list()
  if(all(!is.null(initial_values))){
    for(i in 1:n.chain){
      A[[i]] <- list(
        if('mu' %in% names(initial_values[[i]])){
          mu_trad <- initial_values[[i]]$mu[L_match_trad$L]
        } else {
          mu_trad <- mu_means_trad
        },

        if('mu' %in% names(initial_values[[i]])){
          mu <- initial_values[[i]]$mu
        } else {
          mu <- mu_means_all
        },

        if('p10' %in% names(initial_values[[i]])){
          log_p10 <- log(initial_values[[i]]$p10)
        } else {
          log_p10 <- stats::runif(1,log(0.0001),log(0.08))
        },

        if('beta' %in% names(initial_values[[i]])){
          beta <- initial_values[[i]]$beta
        } else {
          beta <- 0.5
        }
      )
      names(A[[i]]) <- c('mu_trad','mu','log_p10','beta')
    }
  } else {
    for(i in 1:n.chain){
      A[[i]] <- list(
        mu_trad <- mu_means_trad,
        mu <- mu_means_all,
        log_p10 <- stats::runif(1,log(0.0001),log(0.08)),
        beta <- 0.5
      )
      names(A[[i]]) <- c('mu_trad','mu','log_p10','beta')
    }
  }

  return(A)
}

################
# input checks #
################
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#'   messages

# input checks if catchabilty coefficients are used
catchability_checks <- function(data,cov){

  #' @srrstats {G2.13} Pre-processing routines to check for missing data
  ## All tags in data are valid (i.e., include qPCR.N, qPCR.K, count,
  ## count.type, and site.cov)
  #cov='None'
  if (all(is.null(cov)) && !all(c('qPCR.N', 'qPCR.K',
                                  'count','count.type') %in% names(data))){
    errMsg1 <- paste0("Data should include 'qPCR.N', 'qPCR.K', ",
                      "'count', and 'count.type'.")
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase3.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
  #q=TRUE and cov != 'None'
  if (all(!is.null(cov)) && !all(c('qPCR.N', 'qPCR.K', 'count',
                                   'count.type','site.cov') %in% names(data))){
    errMsg1 <- paste0("Data should include 'qPCR.N', 'qPCR.K', ",
                      "'count', 'count.type', and 'site.cov'.")
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase3.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## make sure count.type is not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for zero-length
  #'   data
  if (dim(data$count.type)[1] == 0) {
    errMsg1 <- "count.type contains zero-length data."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase3.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
  ## make sure no column is entirely NA in count.type
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with
  #' all NA
  if (any(apply(data$count.type, 2, function(col) all(is.na(col))))) {
    errMsg1 <- "count.type contains a column with all NA."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase3.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## make sure dimensions of count and count.type are equal, if
  ## count.type is present
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #'   is dimensionally commensurate
  if(dim(data$count)[1] != dim(data$count.type)[1]|
     dim(data$count)[2] != dim(data$count.type)[2]) {
    errMsg1 <- "Dimensions of count and count.type do not match."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase3.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## make sure all data is numeric -- if q == TRUE
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(is.numeric(data$qPCR.K) == FALSE |
     is.numeric(data$qPCR.N) == FALSE |
     is.numeric(data$count) == FALSE |
     is.numeric(data$count.type) == FALSE) {
    errMsg1 <- "Data should be numeric."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase3.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
  ## make sure locations of NAs in count data match locations of NAs in
  ## count.type data
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input
  #'   data is dimensionally commensurate
  if(any((which(is.na(data$count)) == which(is.na(data$count.type))) == FALSE)){
    errMsg1 <- paste0("Empty data cells (NA) in count data should match ",
                      "empty data cells (NA) in count.type data.")
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase3.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
  ## the smallest count.type is 1
  if(min(data$count.type, na.rm =  TRUE) != 1){
    errMsg1 <- paste0("The first gear type should be referenced as 1 in ",
                      "count.type. Subsequent gear types should be referenced ",
                      "2, 3, 4, etc.")
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase3.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## count.type are integers
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(!all(data$count.type %% 1 %in% c(0,NA))){
    errMsg1 <- "All values in count.type should be integers."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase3.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
}

# input checks if no catchabilty coefficients are used
no_catchability_checks <- function(data,cov){

  #' @srrstats {G2.13} Pre-processing routines to check for missing data
  ## All tags in data are valid (i.e., include qPCR.N, qPCR.K, count,
  ## and site.cov)
  #cov='None'
  if (all(is.null(cov)) &&
      !all(c('qPCR.N', 'qPCR.K', 'count') %in% names(data))){
    errMsg1 <- "Data should include 'qPCR.N', 'qPCR.K', and 'count'."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
  #cov != 'None'
  if (all(!is.null(cov)) &&
      !all(c('qPCR.N', 'qPCR.K', 'count','site.cov') %in% names(data))){
    errMsg1 <- "Data should include 'qPCR.N', 'qPCR.K', 'count', and 'site.cov'."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase2.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## make sure all data is numeric -- if q == FALSE
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(is.numeric(data$qPCR.K) == FALSE |
     is.numeric(data$qPCR.N) == FALSE |
     is.numeric(data$count) == FALSE ) {
    errMsg1 <- "Data should be numeric."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
}

# input checks for all variations
all_checks <- function(data, cov, family, p10priors, phipriors, n.chain,
                       n.iter.burn, n.iter.sample, thin, adapt_delta, seed){


  ## make sure count, qPCR.N, and qPCR.K are not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for
  #'   zero-length data
  if (dim(data$qPCR.N)[1] == 0 | dim(data$qPCR.K)[1] == 0 |
      dim(data$count)[1] == 0) {
    errMsg1 <- "Input data contains zero-length data."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
  ## make sure no column is entirely NA in qPCR.N
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column
  #'   with all NA
  if (any(apply(data$qPCR.N, 2, function(col) all(is.na(col))))) {
    errMsg1 <- "qPCR.N contains a column with all NA."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## make sure no column is entirely NA in qPCR.K
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with
  #'   all NA
  if (any(apply(data$qPCR.K, 2, function(col) all(is.na(col))))) {
    errMsg1 <- "qPCR.K contains a column with all NA."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## make sure no column is entirely NA in count
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with
  #'   all NA
  if (any(apply(data$count, 2, function(col) all(is.na(col))))) {
    errMsg1 <- "count contains a column with all NA."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## make sure dimensions of qPCR.N and qPCR.K are equal
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #'   is dimensionally commensurate
  if (dim(data$qPCR.N)[1] != dim(data$qPCR.K)[1]|
      dim(data$qPCR.N)[2] != dim(data$qPCR.K)[2]) {
    errMsg1 <- "Dimensions of qPCR.N and qPCR.K do not match."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
  ## make sure number of rows in count = number of rows in qPCR.N and qPCR.K
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #'   is dimensionally commensurate
  if (dim(data$qPCR.N)[1] != dim(data$count)[1]) {
    errMsg1 <- paste0("Number of sites (rows) in qPCR data and traditional ",
                      "survey count data do not match.")
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## make sure locations of NAs in qPCR.N data match locations of NAs in
  ## qPCR.K data
  #' @srrstats {BS2.1,G2.13} Pre-processing routines to ensure all input data
  #'   is dimensionally commensurate
  if(any((which(is.na(data$qPCR.N)) == which(is.na(data$qPCR.K))) == FALSE)){
    errMsg1 <- paste0("Empty data cells (NA) in qPCR.N data should match ",
                      "empty data cells (NA) in qPCR.K data.")
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## make sure family is either 'poisson', 'negbin', or 'gamma'
  #' @srrstats {G2.3,G2.3a,G2.3b} Permit only expected univariate
  #'   (case-insensitive) parameter values
  if(!c(tolower(family) %in% c('poisson','negbin','gamma'))){
    errMsg <- paste0("Invalid family. Options include 'poisson', 'negbin', ",
                     "and 'gamma'.")
    stop(errMsg)
  }

  ## p10priors is a vector of two integers
  #' @srrstats {G2.0,BS2.2,BS2.3,BS2.4,BS2.5,BS2.6} Checks of vector length and
  #'   appropriateness of distributional parameters (i.e., vector of length 2,
  #'   numeric values > 0), implemented prior to analytic routines
  if(!is.numeric(p10priors) | length(p10priors) != 2 | any(p10priors<=0)){
    errMsg <- paste0("p10priors should be a vector of two positive numeric ",
                     "values. ex. c(1,20)")
    stop(errMsg)
  }

  ## phipriors is a vector of two numeric values
  #' @srrstats {G2.0,BS2.2,BS2.3,BS2.4,BS2.5,BS2.6} Checks of vector length
  #'   and appropriateness of distributional parameters (i.e., vector of length
  #'   2, numeric values > 0), implemented prior to analytic routines
  if(family == 'negbin'){
    if(!is.numeric(phipriors) | length(phipriors) != 2 | any(phipriors<=0)){
      errMsg <- paste0("phipriors should be a vector of two positive numeric ",
                       "values. ex. c(0.25,0.25)")
      stop(errMsg)
    }
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  if(any(data$count == Inf, na.rm =  TRUE) | any(data$count == -Inf,
                                                 na.rm =  TRUE)){
    errMsg1 <- "count contains undefined values (i.e., Inf or -Inf)"
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## count are integers, if family is poisson or negbin
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #'   for distributional parameters (i.e., count data must be an integer if a
  #'   poisson or negative binomial distribution is used), implemented prior to
  #'   analytic routines
  if(tolower(family) %in% c('poisson','negbin')){
    if(!all(data$count %% 1 %in% c(0,NA)) | any(data$count < 0, na.rm =  TRUE)){
      errMsg <- paste0("All values in count should be non-negative integers. ",
                       "Use family = 'gamma' if count is continuous.")
      stop(errMsg)
    }
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  if(any(data$qPCR.N == Inf, na.rm = TRUE) | any(data$qPCR.N == -Inf,
                                                 na.rm =  TRUE)){
    errMsg1 <- "qPCR.N contains undefined values (i.e., Inf or -Inf)"
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  if(any(data$qPCR.K == Inf, na.rm =  TRUE) | any(data$qPCR.K == -Inf,
                                                  na.rm =  TRUE)){
    errMsg1 <- "qPCR.K contains undefined values (i.e., Inf or -Inf)"
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## qPCR.N are integers
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #'   for distributional parameters (i.e., qPCR data are non-negative
  #'   integers), implemented prior to analytic routines
  if(!all(data$qPCR.N %% 1 %in% c(0,NA)) | any(data$qPCR.N < 0,
                                               na.rm =  TRUE)){
    errMsg1 <- "All values in qPCR.N should be non-negative integers."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## qPCR.K are integers
  #' @srrstats {BS2.5} Checks of appropriateness of numeric values submitted
  #'   for distributional parameters (i.e., qPCR data are non-negative
  #'   integers), implemented prior to analytic routines
  if(!all(data$qPCR.K %% 1 %in% c(0,NA)) | any(data$qPCR.K < 0,
                                               na.rm =  TRUE)){
    errMsg1 <- "All values in qPCR.K should be non-negative integers."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## check length and range of n.chain
  if(any(length(as.integer(n.chain)) > 1 | n.chain < 1)){
    errMsg <- "n.chain should be an integer > 0 and of length 1."
    stop(errMsg)
  }

  ## check length and range of n.iter.sample
  if(any(length(as.integer(n.iter.sample)) > 1 | n.iter.sample < 1)){
    errMsg <- "n.iter.sample should be an integer > 0 and of length 1."
    stop(errMsg)
  }

  ## check length and range of n.iter.burn
  if(any(length(as.integer(n.iter.burn)) > 1 | n.iter.burn < 1)){
    errMsg <- "n.iter.burn should be an integer > 0 and of length 1."
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

  ## check that N >= K
  if(any(data$qPCR.K > data$qPCR.N, na.rm =  TRUE)){
    errMsg1 <- paste0("N should be >= K in qPCR data. N is the number of qPCR ",
                      "replicates per sample, and K is the number of positive ",
                      "detections among replicates.")
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }




}

# input checks if site-level covariates are used
covariate_checks <- function(data,cov){

  ## make sure site.cov is not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for zero-length
  #'   data
  if (dim(data$site.cov)[1] == 0) {
    errMsg1 <- "site.cov contains zero-length data."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase2.html')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
  ## make sure no column is entirely NA in site.cov
  #' @srrstats {G5.8,G5.8c} Pre-processing routines to check for column with
  #'   all NA
  if (any(apply(data$site.cov, 2, function(col) all(is.na(col))))) {
    errMsg1 <- "site.cov contains a column with all NA."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase2.html')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## site.cov is numeric, if present
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(!is.numeric(data$site.cov)){
    errMsg1 <- "site.cov should be numeric."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase2.html')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  if(any(data$site.cov == Inf) | any(data$site.cov == -Inf)){
    errMsg1 <- "site.cov contains undefined values (i.e., Inf or -Inf)"
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase2.html')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## cov values match column names in site.cov
  if(!all(cov %in% colnames(data$site.cov))){
    errMsg1 <- paste0("cov values should be listed in the column names of ",
                      "site.cov in the data.")
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase2.html')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## site.cov has same number of rows as qPCR.N and count, if present
  #' @srrstats {BS2.1} Pre-processing routines to ensure all input data is
  #'   dimensionally commensurate
  if(dim(data$qPCR.N)[1] != dim(data$site.cov)[1]){
    errMsg1 <- paste0("The number of rows in site.cov matrix should match the ",
                      "number of rows in all other matrices.")
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_vignette',
                      '/usecase2.html')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## add warning if number of covariates is greater than the number of sites
  #' @srrstats {G5.8d} Pre-processing routines to check if data is outside
  #'   scope of algorithm (i.e., # site-level covariates is greater than the
  #'   number of sites)
  if(length(cov)>dim(data$site.cov)[2]){
    warnMsg <- paste0("The number of site-level covariates exceeds the number ",
                      "of sites (i.e., n < p).")
    warning(warnMsg)
  }

  ## add warning if number of site-covariate data has perfect collinearity
  #' @srrstats {BS3.1} Pre-processing routines to check if site covariate
  #'   data has perfect collinearity
  rank_mat <- qr(data$site.cov)$rank
  if(rank_mat < ncol(data$site.cov)){
    warnMsg <- "Data in site.cov exhibits perfect collinearity."
    warning(warnMsg)
  }
}

# checks if initial values are provided
initial_values_checks <- function(initial_values,data,cov,n.chain){

  ## length of initial values is equal to the number of chains
  if(length(initial_values) != n.chain){
    errMsg1 <- paste0("The length of the list of initial values should equal ",
                      "the number of chains (n.chain, default is 4).")
    errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                      'initial values: ')
    errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_',
                      'vignette/usecase1.html#initialvalues')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  for(i in 1:n.chain){

    ## check mu input
    if('mu' %in% names(initial_values[[i]])){
      ## if mu is numeric
      if(any(!is.numeric(initial_values[[i]]$mu)) |
         any(initial_values[[i]]$mu < 0)){
        errMsg1 <- "Initial values for 'mu' should be numeric values > 0."
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_',
                          'vignette/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
      ## check mu length
      if(length(initial_values[[i]]$mu) != dim(data$count)[1]){
        errMsg1 <- paste0("The length of initial values for 'mu' should equal ",
                          "the number of sites.")
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_',
                          'vignette/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
    }

    ## check p10 input
    if('p10' %in% names(initial_values[[i]])){
      ## if p10 is numeric
      if(!is.numeric(initial_values[[i]]$p10)){
        errMsg1 <- "Initial values for 'p10' should be numeric."
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_',
                          'vignette/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
      ## check p10 length
      if(length(initial_values[[i]]$p10) != 1){
        errMsg1 <- "The length of initial values for 'p10' should equal 1."
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_',
                          'vignette/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
    }

    ## check beta input
    if('beta' %in% names(initial_values[[i]])){
      ## if beta is numeric
      if(!is.numeric(initial_values[[i]]$beta)){
        errMsg1 <- "Initial values for 'beta' should be numeric."
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_',
                          'vignette/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
      ## check beta length
      if(length(initial_values[[i]]$beta) != 1){
        errMsg1 <- "The length of initial values for 'beta' should equal 1."
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_',
                          'vignette/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
    }

    ## check alpha input
    if('alpha' %in% names(initial_values[[i]])){
      ## if alpha is numeric
      if(any(!is.numeric(initial_values[[i]]$alpha))){
        errMsg1 <- "Initial values for 'alpha' should be numeric."
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_',
                          'vignette/usecase2.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
      ## check alpha length
      if(length(initial_values[[i]]$alpha) != (length(cov)+1)){
        errMsg1 <- paste0("The length of initial values for 'alpha' should ",
                          "equal: # covariates + 1 (i.e., including intercept).")
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_',
                          'vignette/usecase2.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
    }

    ## check q input
    if('q' %in% names(initial_values[[i]])){
      ## if q is numeric
      if(any(!is.numeric(initial_values[[i]]$q)) |
         any(initial_values[[i]]$q < 0)){
        errMsg1 <- "Initial values for 'q' should be numeric."
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_',
                          'vignette/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
      ## check q length
      if(length(initial_values[[i]]$q) != (length(table(data$count.type))-1)){
        errMsg1 <- paste0("The length of initial values for 'q' should equal:",
                          " # unique gear types - 1 (i.e., q for reference ",
                          "type = 1).")
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://bookdown.org/abigailkeller/eDNAjoint_',
                          'vignette/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
    }


  }
}



