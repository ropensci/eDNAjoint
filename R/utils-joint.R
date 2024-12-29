## helper functions for jointModel.R ##

####################################
# helper functions: get model type #
####################################

#' @noRd
isCatch_type <- function(q){
  out <- ifelse(q == TRUE,TRUE,FALSE)
  return(out)
}
#' @noRd
isCov_type <- function(cov){
  out <- ifelse(all(!is.null(cov)),TRUE,FALSE)
  return(out)
}
#' @noRd
isNegbin_type <- function(family){
  out <- ifelse(family == 'negbin',TRUE,FALSE)
  return(out)
}
#' @noRd
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

#' @noRd
get_stan_model <- function(family){

  index <- ifelse(family %in% c('poisson','negbin'),1,2)

  return(index)
}


####################################
# helper functions: initial values #
####################################
#' @srrstats {BS2.7,BS2.11} Option for user to provide initial values for each
#'   chain

#' @noRd
get_inits <- function(n.chain,qPCR_all,initial_values,cov,L_match_trad,
                      L_match_dna,data,q_names=NULL){
  if(!is.null(q_names)){
    inits <- init_joint_cov_catchability(n.chain,qPCR_all,q_names,cov,
                                         initial_values,L_match_trad,
                                         L_match_dna,data)
    } else {
      inits <- init_joint_cov(n.chain,qPCR_all,cov,initial_values,
                              L_match_trad,L_match_dna,data)
    }

  return(inits)
}

#' @noRd
init_joint_cov <- function(n.chain,qPCR_all,cov,initial_values,
                           L_match_trad,L_match_dna,data){

  # get mu means
  mu_means_trad <- as.vector(stats::na.omit(rowMeans(data$count,
                                                     na.rm=TRUE)+0.01))
  mu_means_all <- rep(NA,dim(L_match_dna)[1]+dim(L_match_trad)[1])
  mu_means_all[L_match_trad$L_ind] <- mu_means_trad
  if(dim(L_match_dna)[1]>0){
    mu_means_all[L_match_dna$L_ind] <- rep(mean(mu_means_trad),
                                           dim(L_match_dna)[1])
  }

  # helper function
  # joint model, catchability coefficient, site covariates
  A <- list()
  if(all(!is.null(initial_values))){
    for(i in 1:n.chain){
      A[[i]] <- list(
        if('mu' %in% names(initial_values[[i]])){
          mu_trad <- initial_values[[i]]$mu[L_match_trad$L_ind]
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
          log_p10 <- stats::runif(1,log(0.0001),log(0.01))
        },

        if('alpha' %in% names(initial_values[[i]])){
          alpha <- as.array(initial_values[[i]]$alpha)
        } else {
          alpha <- as.array(c(3.5,rep(0,length(cov))))
        },

        p_dna <- rep(0.4,dim(L_match_dna)[1]),
        p11_dna <- rep(0.4,dim(L_match_dna)[1]) - 0.01
      )
      names(A[[i]]) <- c('mu_trad','mu','log_p10','alpha','p_dna','p11_dna')
    }
  } else {
    for(i in 1:n.chain){
      A[[i]] <- list(
        mu_trad <- mu_means_trad,
        mu <- mu_means_all,
        log_p10 <- stats::runif(1,log(0.0001),log(0.01)),
        alpha <- as.array(c(3.5,rep(0,length(cov)))),
        p_dna <- rep(0.4,dim(L_match_dna)[1]),
        p11_dna <- rep(0.4,dim(L_match_dna)[1]) - 0.01
      )
      names(A[[i]]) <- c('mu_trad','mu','log_p10','alpha','p_dna','p11_dna')
    }
  }

  return(A)
}

#' @noRd
init_joint_cov_catchability <- function(n.chain,qPCR_all,q_names,cov,
                                        initial_values,L_match_trad,
                                        L_match_dna,data){

  # get mu means
  mu_means_trad <- as.vector(stats::na.omit(rowMeans(data$count,
                                                     na.rm=TRUE)+0.01))
  mu_means_all <- rep(NA,dim(L_match_dna)[1]+dim(L_match_trad)[1])
  mu_means_all[L_match_trad$L_ind] <- mu_means_trad
  if(dim(L_match_dna)[1]>0){
    mu_means_all[L_match_dna$L_ind] <- rep(mean(mu_means_trad),
                                           dim(L_match_dna)[1])
  }

  # helper function
  # joint model, catchability coefficient, site covariates
  A <- list()
  if(all(!is.null(initial_values))){
    for(i in 1:n.chain){
      A[[i]] <- list(
        if('mu' %in% names(initial_values[[i]])){
          mu_trad <- initial_values[[i]]$mu[L_match_trad$L_ind]
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
          log_p10 <- stats::runif(1,log(0.0001),log(0.01))
        },

        if('alpha' %in% names(initial_values[[i]])){
          alpha <- as.array(initial_values[[i]]$alpha)
        } else {
          alpha <- as.array(c(3.5,rep(0,length(cov))))
        },

        if('q' %in% names(initial_values[[i]])){
          q_trans <- as.data.frame(initial_values[[i]]$q)
        } else {
          q_trans <- as.data.frame(stats::runif(length(q_names),0.01,1))
        },
        p_dna <- rep(0.4,dim(L_match_dna)[1]),
        p11_dna <- rep(0.4,dim(L_match_dna)[1]) - 0.01
      )
      names(A[[i]]) <- c('mu_trad','mu','log_p10','alpha','q_trans',
                         'p_dna','p11_dna')

    }
  } else {
    for(i in 1:n.chain){
      A[[i]] <- list(
        mu_trad <- mu_means_trad,
        mu <- mu_means_all,
        log_p10 <- stats::runif(1,log(0.0001),log(0.01)),
        alpha <- as.array(c(3.5,rep(0,length(cov)))),
        q_trans <- as.data.frame(stats::runif(length(q_names),0.01,1)),
        p_dna <- rep(0.4,dim(L_match_dna)[1]),
        p11_dna <- rep(0.4,dim(L_match_dna)[1]) - 0.01
      )
      names(A[[i]]) <- c('mu_trad','mu','log_p10','alpha','q_trans',
                         'p_dna','p11_dna')
    }
  }

  return(A)
}



################
# input checks #
################
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#'   messages

#' @noRd
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
                      '/usecase3.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
}

#' @noRd
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
  #cov != 'None'
  if (all(!is.null(cov)) &&
      !all(c('qPCR.N', 'qPCR.K', 'count','site.cov') %in% names(data))){
    errMsg1 <- "Data should include 'qPCR.N', 'qPCR.K', 'count', and 'site.cov'."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }
}

#' @noRd
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
                      '/usecase1.html#prepare-the-data')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }




}

#' @noRd
# input checks if site-level covariates are used
covariate_checks <- function(data,cov){

  ## make sure site.cov is not zero-length
  #' @srrstats {G5.8,G5.8a} Pre-processing routines to check for zero-length
  #'   data
  if (dim(data$site.cov)[1] == 0) {
    errMsg1 <- "site.cov contains zero-length data."
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
                      '/usecase2.html')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## make sure no data are undefined
  #' @srrstats {G2.16} Pre-processing routines to check for undefined data
  if(any(data$site.cov == Inf) | any(data$site.cov == -Inf)){
    errMsg1 <- "site.cov contains undefined values (i.e., Inf or -Inf)"
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://ednajoint.netlify.app',
                      '/usecase2.html')
    errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
    stop(errMsg)
  }

  ## cov values match column names in site.cov
  if(!all(cov %in% colnames(data$site.cov))){
    errMsg1 <- paste0("cov values should be listed in the column names of ",
                      "site.cov in the data.")
    errMsg2 <- 'See the eDNAjoint guide for data formatting help: '
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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
    errMsg3 <- paste0('https://ednajoint.netlify.app',
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

#' @noRd
# checks if initial values are provided
initial_values_checks <- function(initial_values,data,cov,n.chain){

  ## length of initial values is equal to the number of chains
  if(length(initial_values) != n.chain){
    errMsg1 <- paste0("The length of the list of initial values should equal ",
                      "the number of chains (n.chain, default is 4).")
    errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                      'initial values: ')
    errMsg3 <- paste0('https://ednajoint.netlify.app',
                      '/usecase1.html#initialvalues')
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
        errMsg3 <- paste0('https://ednajoint.netlify.app',
                          '/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
      ## check mu length
      if(length(initial_values[[i]]$mu) != dim(data$count)[1]){
        errMsg1 <- paste0("The length of initial values for 'mu' should equal ",
                          "the number of sites.")
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://ednajoint.netlify.app',
                          '/usecase1.html#initialvalues')
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
        errMsg3 <- paste0('https://ednajoint.netlify.app',
                          '/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
      ## check p10 length
      if(length(initial_values[[i]]$p10) != 1){
        errMsg1 <- "The length of initial values for 'p10' should equal 1."
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://ednajoint.netlify.app',
                          '/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
    }

    ## check alpha input -- no covariates
    if('alpha' %in% names(initial_values[[i]]) && is.null(cov)){
      ## if alpha is numeric
      if(!is.numeric(initial_values[[i]]$alpha)){
        errMsg1 <- "Initial values for 'alpha' should be numeric."
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://ednajoint.netlify.app',
                          '/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
      ## check alpha length
      if(length(initial_values[[i]]$alpha) != 1){
        errMsg1 <- "The length of initial values for 'alpha' should equal 1."
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://ednajoint.netlify.app',
                          '/usecase1.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
    }

    ## check alpha input -- covariates
    if('alpha' %in% names(initial_values[[i]]) && !is.null(cov)){
      ## if alpha is numeric
      if(any(!is.numeric(initial_values[[i]]$alpha))){
        errMsg1 <- "Initial values for 'alpha' should be numeric."
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://ednajoint.netlify.app',
                          '/usecase2.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
      ## check alpha length
      if(length(initial_values[[i]]$alpha) != (length(cov)+1)){
        errMsg1 <- paste0("The length of initial values for 'alpha' should ",
                          "equal: # covariates + 1 (i.e., including intercept).")
        errMsg2 <- paste0('See the eDNAjoint guide for help formatting ',
                          'initial values: ')
        errMsg3 <- paste0('https://ednajoint.netlify.app',
                          '/usecase2.html#initialvalues')
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
        errMsg3 <- paste0('https://ednajoint.netlify.app',
                          '/usecase2.html#initialvalues')
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
        errMsg3 <- paste0('https://ednajoint.netlify.app',
                          '/usecase2.html#initialvalues')
        errMsg <- paste(errMsg1,errMsg2,errMsg3, sep = "\n")
        stop(errMsg)
      }
    }


  }
}



