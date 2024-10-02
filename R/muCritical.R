#' Calculate mu_critical
#'
#' This function uses the full posterior distributions of parameters estimated
#' by `jointModel()` to calculate mu_critical, or the expected catch rate at
#' which the probabilities of a false positive eDNA detection and true positive
#' eDNA detection are equal. See more examples in the
#' \href{https://bookdown.org/abigailkeller/eDNAjoint_vignette/}{Package
#' Vignette}.
#'
#' @srrstats {G1.4} Roxygen function documentation begins here
#' @export
#' @srrstats {G2.1a} Here are explicit documentation of vector input types
#' @param modelfit An object of class `stanfit`
#' @param cov.val A numeric vector indicating the values of site-level
#'   covariates to use for prediction. Default is NULL.
#' @param ci Credible interval calculated using highest density interval (HDI).
#'   Default is 0.9 (i.e., 90% credible interval).
#' @return A list with median mu_critical and lower and upper bounds on the
#'   credible interval. If multiple gear types are used, a table of mu_critical
#'   and lower and upper credible interval bounds is returned with one column
#'   for each gear type.
#'
#' @srrstats {G2.0a,G2.2} Explicit secondary documentation of any expectations
#'   on lengths of inputs
#' @note  Before fitting the model, this function checks to ensure that the
#'   function is possible given the inputs. These checks include:
#' \itemize{
#' \item  Input model fit is an object of class 'stanfit'.
#' \item  Input credible interval is a univariate numeric value greater than 0
#'   and less than 1.
#' \item  Input model fit contains p10 parameter.
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
#' # Ex. 1: Calculating mu_critical with site-level covariates
#'
#' # Load data
#' data(gobyData)
#'
#' # Fit a model including 'Filter_time' and 'Salinity' site-level covariates
#' fit.cov <- jointModel(data = gobyData, cov = c('Filter_time','Salinity'),
#'                       family = "poisson", p10priors = c(1,20), q = FALSE,
#'                       multicore = FALSE)
#'
#' # Calculate mu_critical at the mean covariate values (covariates are
#' # standardized, so mean = 0)
#' muCritical(fit.cov$model, cov.val = c(0,0), ci = 0.9)
#'
#' # Calculate mu_critical at habitat size 0.5 z-scores greater than the mean
#' muCritical(fit.cov$model, cov.val = c(0,0.5), ci = 0.9)
#'
#' # Ex. 2: Calculating mu_critical with multiple traditional gear types
#'
#' # Load data
#' data(greencrabData)
#'
#' # Fit a model with no site-level covariates
#' fit.q <- jointModel(data = greencrabData, cov = NULL, family = "negbin",
#'                     p10priors = c(1,20), q = TRUE, multicore = FALSE)
#'
#' # Calculate mu_critical
#' muCritical(fit.q$model, cov.val = NULL, ci = 0.9)
#' }
#'

muCritical <- function(modelfit, cov.val = NULL, ci = 0.9) {

  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using this helper
  #'   function
  muCritical_input_checks(modelfit, cov.val, ci)

  if (!requireNamespace("rstan", quietly = TRUE)){
    stop ("The 'rstan' package is not installed.", call. = FALSE)
  }

  ## check to see if there are any divergent transitions
  #' @srrstats {BS4.5} Warning message if the input model fit has divergence
  #'   transitions
  if(sum(lapply(rstan::get_sampler_params(modelfit,
                                          inc_warmup = FALSE),
                div_check)[[1]]) > 0){

    sum <- sum(lapply(rstan::get_sampler_params(modelfit, inc_warmup = FALSE),
                      div_check)[[1]])

    warning_msg <- paste0('Warning: There are ', sum,
                          ' divergent transitions in your model fit. ')
    warning(warning_msg)

  }

  if(all(is.null(cov.val))){
    #extract posteriors for beta and p10 parameters
    ##beta
    posterior_beta <- unlist(rstan::extract(modelfit, pars = "beta"))
    #p10
    posterior_p10 <- unlist(rstan::extract(modelfit, pars = "p10"))
  } else {
    #extract posteriors for beta and p10 parameters
    ##alpha
    posterior_alpha <- rstan::extract(modelfit, pars = "alpha")$alpha
    ##beta
    posterior_beta <- posterior_alpha %*% c(1, cov.val)
    #p10
    posterior_p10 <- unlist(rstan::extract(modelfit, pars = "p10"))
  }

  if('q' %in% modelfit@model_pars){
    #extract q values
    posterior_q <- rstan::extract(modelfit, pars = "q")$q

    #create empty dataframe
    out <- as.data.frame(matrix(NA, nrow = 3, ncol = modelfit@par_dims$q+1))
    rownames(out) <- c('median','lower_ci','upper_ci')
    for(i in 1:modelfit@par_dims$q){
      gear_names <- paste('gear_', i+1, sep = '')
    }
    colnames(out) <- c('gear_1', gear_names)

    #calculate mu_critical -- gear type 1
    critical_mu_1 <- rep(NA, length(posterior_beta))
    for(i in seq_along(critical_mu_1)){
      critical_mu_1[i] <- (
        posterior_p10[i]*exp(posterior_beta[i])/(1-posterior_p10[i])
      )
    }
    out[,1] <- c(stats::median(critical_mu_1),
                 bayestestR::ci(critical_mu_1, method = 'HDI',
                                ci = ci)[2]$CI_low,
                 bayestestR::ci(critical_mu_1, method = 'HDI',
                                ci = ci)[3]$CI_high)

    #calculate mu_critical -- gear type 2+
    for(i in 1:modelfit@par_dims$q){
      out[,i+1] <- c(stats::median(critical_mu_1*posterior_q[,i]),
                     bayestestR::ci(critical_mu_1*posterior_q[,i],
                                    method = 'HDI', ci = ci)[2]$CI_low,
                     bayestestR::ci(critical_mu_1*posterior_q[,i],
                                    method = 'HDI', ci = ci)[3]$CI_high)
    }

  } else {
    #calculate mu_critical
    critical_mu <- rep(NA, length(posterior_beta))
    for(i in seq_along(critical_mu)){
      critical_mu[i] <- (
        posterior_p10[i]*exp(posterior_beta[i])/(1-posterior_p10[i])
      )
    }

    out <- list(median = stats::median(critical_mu),
                lower_ci = bayestestR::ci(critical_mu, method = 'HDI',
                                        ci = ci)[2],
                upper_ci = bayestestR::ci(critical_mu, method = 'HDI',
                                        ci = ci)[3])

  }

  return(out)

}

#function to check for divergent transitions
div_check <- function(x){
  divergent <- sum(x[,'divergent__'])
  return(divergent)
}

# function for input checks
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#'   messages
muCritical_input_checks <- function(modelfit, cov.val, ci){
  ## #1. make sure model fit is of class stanfit
  #' @srrstats {G2.8} Makes sure input of sub-function is of class 'stanfit'
  #'   (i.e., output of jointModel())
  if(!is(modelfit,'stanfit')) {
    errMsg <- "modelfit must be of class 'stanfit'."
    stop(errMsg)
  }

  ## #2. make sure ci is valid
  #' @srrstats {G2.0,G2.2} Assertion on length of input, prohibit multivariate
  #'   input to parameters expected to be univariate
  if(!is.numeric(ci) | ci <= 0 | ci >= 1 | length(ci) > 1) {
    errMsg <- "ci must be a numeric value >0 and <1."
    stop(errMsg)
  }

  ## #3. make sure model fit contains p10 parameter
  if(!("p10" %in% modelfit@model_pars)) {
    errMsg <- "modelfit must contain 'p10' parameter."
    stop(errMsg)
  }

  ## #4. if modelfit contains alpha, cov.val must be provided
  if('alpha' %in% modelfit@model_pars && all(is.null(cov.val))) {
    errMsg <- paste0("If modelfit contains site-level covariates, values ",
                     "must be provided for cov.val")
    stop(errMsg)
  }

  ## #5. cov.val is numeric, if provided
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(all(!is.null(cov.val)) && !is.numeric(cov.val)) {
    errMsg <- "cov.val must be a numeric vector"
    stop(errMsg)
  }

  ## #6. Only include input cov.val if covariates are included in model
  if(all(!is.null(cov.val)) && !c('alpha') %in% modelfit@model_pars) {
    errMsg <- paste0("cov.val must be NULL if the model does not contain ",
                     "site-level covariates.")
    stop(errMsg)
  }

  ## #7. Input cov.val is the same length as the number of estimated covariates.
  #' @srrstats {G2.0} Assertion on length of input
  if(all(!is.null(cov.val)) && length(cov.val) != (modelfit@par_dims$alpha-1)) {
    errMsg <- paste0("cov.val must be of the same length as the number of ",
                     "non-intercept site-level coefficients in the model.")
    stop(errMsg)
  }
}

