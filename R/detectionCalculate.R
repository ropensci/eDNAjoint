#' Calculate the survey effort necessary to detect species presence, given the
#' species expected catch rate.
#'
#' This function calculates the number of survey effort units to necessary
#' detect species presence using median estimated parameter values from
#' jointModel(). Detecting species presence is defined
#' as producing at least one true positive eDNA detection or catching at least
#' one individual. See more examples in the
#' \href{https://bookdown.org/abigailkeller/eDNAjoint_vignette/}{Package
#' Vignette}.
#'
#' @srrstats {G1.4} Roxygen function documentation begins here
#' @export
#' @srrstats {G2.1a} Here are explicit documentation of vector input types
#' @param modelfit An object of class `stanfit`.
#' @param mu A numeric vector of species densities/capture rates. If multiple
#'   traditional gear types are represented in the model, mu is the catch rate
#'   of gear type 1.
#' @param cov.val A numeric vector indicating the values of site-level
#'   covariates to use for prediction. Default is NULL.
#' @param probability A numeric value indicating the probability of detecting
#'   presence. The default is 0.9.
#' @param qPCR.N An integer indicating the number of qPCR replicates per eDNA
#'   sample. The default is 3.
#' @return A summary table of survey efforts necessary to detect species
#'   presence, given mu, for each survey type.
#'
#' @srrstats {G2.0a,G2.2} Explicit secondary documentation of any expectations
#'   on lengths of inputs
#' @note  Before fitting the model, this function checks to ensure that the
#'   function is possible given the inputs. These checks include:
#' \itemize{
#' \item  Input model fit is an object of class 'stanfit'.
#' \item  Input mu is a numeric vector.
#' \item  Input probability is a univariate numeric value.
#' \item  If model fit contains alpha, cov.val must be provided.
#' \item  Input cov.val is numeric.
#' \item  Input cov.val is the same length as the number of estimated
#'   covariates.
#' \item  Input model fit has converged (i.e. no divergent transitions after
#'   warm-up).
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' @examples
#' \donttest{
#' # Ex. 1: Calculating necessary effort for detection with site-level
#' # covariates
#'
#' # Load data
#' data(gobyData)
#'
#' # Fit a model including 'Filter_time' and 'Salinity' site-level covariates
#' fit.cov <- jointModel(data = gobyData, cov = c('Filter_time','Salinity'),
#'                       family = "poisson", p10priors = c(1,20), q = FALSE,
#'                       multicore = FALSE)
#'
#' # Calculate at the mean covariate values
#' # (covariates are standardized, so mean = 0)
#' detectionCalculate(fit.cov$model, mu = seq(from = 0.1, to = 1, by = 0.1),
#'                    cov.val = c(0,0), qPCR.N = 3)
#'
#' # Calculate mu_critical at salinity 0.5 z-scores greater than the mean
#' detectionCalculate(fit.cov$model, mu = seq(from = 0.1, to = 1, by = 0.1),
#'                    cov.val = c(0,0.5), qPCR.N = 3)
#'
#' # Ex. 2: Calculating necessary effort for detection with multiple traditional
#' # gear types
#'
#' # Load data
#' data(greencrabData)
#'
#' # Fit a model with no site-level covariates
#' fit.q <- jointModel(data = greencrabData, cov = NULL, family = "negbin",
#'                     p10priors = c(1,20), q = TRUE, multicore = FALSE)
#'
#' # Calculate
#' detectionCalculate(fit.q$model, mu = seq(from = 0.1, to = 1, by = 0.1),
#'                    cov.val = NULL, qPCR.N = 3)
#'
#' # Change probability of detecting presence to 0.95
#' detectionCalculate(fit.q$model, mu = 0.1, cov.val = NULL,
#'                    probability = 0.95, qPCR.N = 3)
#' }
#'

detectionCalculate <- function(modelfit, mu, cov.val = NULL, probability = 0.9,
                               qPCR.N = 3){

  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using this helper
  #'   function
  detectionCalculate_input_checks(modelfit, mu, cov.val, probability, qPCR.N)

  if (!requireNamespace("rstan", quietly = TRUE)){
    stop ("The 'rstan' package is not installed.", call. = FALSE)
  }

  ## check to see if there are any divergent transitions
  #' @srrstats {BS4.5} Warning message if the input model fit has divergence
  #'   transitions
  if(sum(lapply(rstan::get_sampler_params(modelfit,inc_warmup = FALSE),
                div_check)[[1]]) > 0){

    sum <- sum(lapply(rstan::get_sampler_params(modelfit,inc_warmup = FALSE),
                      div_check)[[1]])

    warning <- paste0('Warning: There are ',sum,
                      ' divergent transitions in your model fit. ')
    print(warning)

  }

  # get n traditional samples
  if(isCatch(modelfit@model_pars)){
    ntrad_out <- get_ntrad_q(modelfit@model_pars, modelfit, mu, probability)
  } else {
    ntrad_out <- get_ntrad(modelfit@model_pars, modelfit, mu, probability)
  }

  # get n dna samples
  if(isJoint(modelfit@model_pars)){
    ndna_out <- get_ndna(modelfit@model_pars, modelfit, mu, qPCR.N,
                         probability, cov.val)

  }

  # combine into one df - joint model
  if(isJoint(modelfit@model_pars)){
    if(isCatch(modelfit@model_pars)){
      out <- cbind(mu,ntrad_out,ndna_out)
      #rename columns
      for(i in 1:modelfit@par_dims$q){
        trad_names <- paste('n_traditional_',i+1,sep = '')
      }
      colnames(out) <- c('mu','n_traditional_1',trad_names,'n_eDNA')
    } else {
      out <- cbind(mu,ntrad_out,ndna_out)
      colnames(out) <- c('mu','n_traditional','n_eDNA')
    }
  }

  # combine into one df - traditional model
  if(!isJoint(modelfit@model_pars)){
    if(isCatch(modelfit@model_pars)){
      out <- cbind(mu,ntrad_out)
      #rename columns
      for(i in 1:modelfit@par_dims$q){
        trad_names <- paste('n_traditional_',i+1,sep = '')
      }
      colnames(out) <- c('mu','n_traditional_1',trad_names)
    } else {
      out <- cbind(mu,ntrad_out)
      colnames(out) <- c('mu','n_traditional')
    }
  }

  return(out)
}

# functions to check model type
isJoint <- function(pars){
  out <- ifelse('p10' %in% pars,TRUE,FALSE)
  return(out)
}
isCatch <- function(pars){
  out <- ifelse('q' %in% pars,TRUE,FALSE)
  return(out)
}
isNegbin <- function(pars){
  out <- ifelse('phi' %in% pars,TRUE,FALSE)
  return(out)
}
isCov <- function(pars){
  out <- ifelse('alpha' %in% pars,TRUE,FALSE)
  return(out)
}

# function to get n traditional samples
get_ntrad <- function(pars, modelfit, mu, probability){

  # number of traditional survey effort
  ntrad_out <- vector(length = length(mu))

  for(i in seq_along(mu)){
    # number of traditional survey replicates
    ntrad <- seq(from = 0, to = 50000, by = 1)

    # P(X = 0 | mu) in one traditional survey trial
    if(isNegbin(pars)){
      phi <- stats::median(unlist(rstan::extract(modelfit, pars = 'phi')))
      pr <- stats::pnbinom(q = 0, mu = mu[i], size = phi)
    } else {
      pr <- stats::ppois(q = 0, lambda = mu[i])
    }

    # dnbinom: x = failed events; size = successful events
    # probability of catching >0 animals based on ntrad
    prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr)

    # find value
    value <- match(min(prob[prob >= probability]), prob)
    ntrad_out[i] <- ntrad[value]
  }
  return(ntrad_out)
}

# function to get n traditional samples -- catchability coefficient
get_ntrad_q <- function(pars, modelfit, mu, probability){

  # get catch coefficients
  # create empty q list
  q_list <- list()

  # fill in q list
  for(i in 1:modelfit@par_dims$q){
    name <- paste('q',i,sep = '')
    q_list[[name]] <- stats::median(rstan::extract(modelfit,
                                                   pars = 'q')$q[,i])
  }
  q_list <- c(1,unlist(q_list))

  # number of traditional survey effort
  ntrad_out <- matrix(NA,nrow = length(mu),ncol = length(q_list))

  for(j in seq_along(q_list)){
    for(i in seq_along(mu)){
      # number of traditional survey replicates
      ntrad <- seq(from = 0, to = 50000, by = 1)

      # P(X = 0 | mu) in one traditional survey trial
      if(isNegbin(pars)){
        phi <- stats::median(unlist(rstan::extract(modelfit, pars = 'phi')))
        pr <- stats::pnbinom(q = 0, mu = q_list[j]*mu[i], size = phi)
      } else {
        pr <- stats::ppois(q = 0, lambda = q_list[j]*mu[i])
      }

      # dnbinom: x = failed events; size = successful events
      # probability of catching >0 animals based on ntrad
      prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr)

      # find value
      value <- match(min(prob[prob >= probability]), prob)
      ntrad_out[i,j] <- ntrad[value]
    }
  }
  return(ntrad_out)
}

# function to get n eDNA samples
get_ndna <- function(pars, modelfit, mu, qPCR.N, probability, cov.val){

  # get beta
  if(isCov(pars)){
    alpha <- apply(rstan::extract(modelfit, pars = 'alpha')$alpha,2,'median')
    beta <- alpha %*% c(1,cov.val)
  } else {
    beta <- stats::median(unlist(rstan::extract(modelfit, pars = 'beta')))
  }

  # create output
  ndna_out <- vector(length = length(mu))

  for(i in seq_along(mu)){
    # number of eDNA samples
    ndna <- seq(from = 0, to = 50000, by = 1)

    p11 <- mu[i]/(mu[i]+exp(beta))
    # P(at least one detection|mu) out of total bottles/PCR replicates
    prob <- 1 - stats::dbinom(x = 0, size = ndna*qPCR.N, prob = p11)

    # find value
    value <- match(min(prob[prob >= probability]), prob)
    ndna_out[i] <- ndna[value]
  }

  return(ndna_out)
}

#function to check for divergent transitions
div_check <- function(x){
  divergent <- sum(x[,'divergent__'])
  return(divergent)
}

# function for input checks
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#'   messages
detectionCalculate_input_checks <- function(modelfit, mu, cov.val,
                                            probability, qPCR.N){
  ## #1. make sure model fit is of class stanfit
  #' @srrstats {G2.8} Makes sure input of sub-function is of class 'stanfit'
  #'   (i.e., output of jointModel())
  if(!is(modelfit,'stanfit')) {
    errMsg <- "modelfit must be of class 'stanfit'."
    stop(errMsg)
  }

  ## #2. make sure mu is a numeric vector
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(!any(is.numeric(mu)) | any(mu <= 0)) {
    errMsg <- "mu must be a numeric vector of positive values"
    stop(errMsg)
  }

  ## #3. make sure probability is a numeric value between 0 and 1
  #' @srrstats {G2.0,G2.2} Assertion on length of input, prohibit multivariate
  #'   input to parameters expected to be univariate
  if(!is.numeric(probability) | length(probability)>1 | any(probability < 0) |
     any(probability > 1)) {
    errMsg <- "probability must be a numeric value between 0 and 1"
    stop(errMsg)
  }

  ## #4. cov.val is numeric, if provided
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(all(!is.null(cov.val)) && !is.numeric(cov.val)) {
    errMsg <- "cov.val must be a numeric vector"
    stop(errMsg)
  }

  ## #5. Only include input cov.val if covariates are included in model
  if(all(!is.null(cov.val)) && !c('alpha') %in% modelfit@model_pars) {
    errMsg <- paste0("cov.val must be NULL if the model does not ",
                     "contain site-level covariates.")
    stop(errMsg)
  }

  ## #6. Input cov.val is the same length as the number of estimated covariates.
  #' @srrstats {G2.0} Assertion on length of input
  if(all(!is.null(cov.val)) &&
     length(cov.val) != (modelfit@par_dims$alpha-1)) {
    errMsg <- paste0("cov.val must be of the same length as the number of ",
                     "non-intercept site-level coefficients in the model.")
    stop(errMsg)
  }

  ## #7. If covariates are in model, cov.val must be provided
  if(all(c('alpha','p10') %in% modelfit@model_pars) && all(is.null(cov.val))) {
    errMsg <- paste0("cov.val must be provided if the model contains ",
                     "site-level covariates.")
    stop(errMsg)
  }

  ## #8. qPCR.N must be an integer
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(qPCR.N %% 1 != 0) {
    errMsg <- "qPCR.N should be an integer."
    stop(errMsg)
  }
}
