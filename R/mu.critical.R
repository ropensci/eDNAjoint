#' Calculate mu_critical
#'
#' This function uses the full posterior distributions of parameters estimated by `jointModel()` to calculate mu_critical, or the expected catch rate at which the probabilities of a false positive eDNA detection and true positive eDNA detection are equal.
#'
#' @export
#' @param modelfit An object of class `stanfit`
#' @param cov.val A numeric vector indicating the values of site-level covariates to use for prediction. Default is 'None'.
#' @param ci Credible interval calculated using highest density interval (HDI). Default is 0.9 (i.e., 90% credible interval).
#' @return A list with median mu_critical and lower and upper bounds on the credible interval.
#'
#' @note  Before fitting the model, this function checks to ensure that the function is possible given the inputs. These checks include:
#' \itemize{
#' \item  Input model fit is an object of class 'stanfit'.
#' \item  Input credible interval is a numeric value greater than 0 and less than 1.
#' \item  Input model fit contains p10 parameter.
#' \item  If model fit contains alpha, cov.val must be provided.
#' \item  Input cov.val is numeric.
#' \item  Input cov.val is the same length as the number of estimated covariates.
#' \item  Input model fit has converged (i.e. no divergent transitions after warm-up).
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
#' # Fit a model including 'Filter_time' and 'Hab_size' site-level covariates
#' fit.cov = jointModel(data=gobyData, cov=c('Filter_time','Hab_size'),
#'                      family='poisson', p10priors=c(1,20), q=FALSE)
#'
#' # Calculate mu_critical at the mean covariate values (covariates are standardized, so mean=0)
#' mu.critical(fit.cov, cov.val=c(0,0))
#'
#' # Calculate mu_critical at habitat size 0.5 z-scores greater than the mean
#' mu.critical(fit.cov, cov.val=c(0,0.5))
#'
#' # Ex. 2: Calculating mu_critical with multiple traditional gear types
#'
#' # Load data
#' data(greencrabData)
#'
#' # Fit a model with no site-level covariates
#' fit.q = jointModel(data=greencrabData, cov='None', family='negbin',
#'                    p10priors=c(1,20), q=TRUE, q_ref=1)
#'
#' # Calculate mu_critical
#' mu.critical(fit.q, cov.val='None')
#' }
#'

mu.critical <- function(modelfit, cov.val = 'None', ci = 0.9) {

  ## #1. make sure model fit is of class stanfit
  if(!is(modelfit,'stanfit')) {
      errMsg = paste("modelfit must be of class 'stanfit'.")
      stop(errMsg)
  }

  ## #2. make sure ci is valid
  if(!is.numeric(ci)|ci<=0|ci>=1) {
    errMsg = paste("ci must be a numeric value >0 and <1.")
    stop(errMsg)
  }

  ## #3. make sure model fit contains p10 parameter
  if(!("p10" %in% modelfit@model_pars)) {
    errMsg = paste("modelfit must be contain 'p10' parameter.")
    stop(errMsg)
  }

  ## #4. if modelfit contains alpha, cov.val must be provided
  if('alpha' %in% modelfit@model_pars && all(cov.val=='None')) {
    errMsg = paste("If modelfit contains site-level covariates, values must be provided for cov.val")
    stop(errMsg)
  }

  ## #5. cov.val is numeric, if provided
  if(all(cov.val != 'None') && !is.numeric(cov.val)) {
    errMsg = paste("cov.val must be a numeric vector")
    stop(errMsg)
  }

  ## #6. Only include input cov.val if covariates are included in model
  if(all(cov.val != 'None') && !c('alpha') %in% modelfit@model_pars) {
    errMsg = paste("cov.val must be 'None' if the model does not contain site-level covariates.")
    stop(errMsg)
  }

  ## #7. Input cov.val is the same length as the number of estimated covariates.
  if(all(cov.val != 'None') && length(cov.val)!=(modelfit@par_dims$alpha-1)) {
    errMsg = paste("cov.val must be of the same length as the number of non-intercept site-level coefficients in the model.")
    stop(errMsg)
  }

  if (!requireNamespace("rstan", quietly = TRUE)){
    stop ("The 'rstan' package is not installed.", call. = FALSE)
  }

  ## #8. check to see if there are any divergent transitions
  if(sum(rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[1]][,'divergent__'],
         rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[2]][,'divergent__'],
         rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[3]][,'divergent__'],
         rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[4]][,'divergent__']) > 0 ){

    sum <- sum(rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[1]][,'divergent__'],
               rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[2]][,'divergent__'],
               rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[3]][,'divergent__'],
               rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[4]][,'divergent__'])

    warning <- paste0('Warning: There are ',sum,' divergent transitions in your model fit. ')
    print(warning)

  }

  if(all(cov.val=='None')){
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
    posterior_beta <- posterior_alpha %*% c(1,cov.val)
    #p10
    posterior_p10 <- unlist(rstan::extract(modelfit, pars = "p10"))
  }

  #calculate mu_critical
  critical_mu <- rep(NA, length(posterior_beta))
  for(i in 1:length(critical_mu)){
    critical_mu[i] <- posterior_p10[i]*exp(posterior_beta[i])/(1-posterior_p10[i])
  }

  out <- list(median=stats::median(critical_mu),
              lower_ci=bayestestR::ci(critical_mu, method = 'HDI', ci = ci)[2],
              upper_ci=bayestestR::ci(critical_mu, method = 'HDI', ci = ci)[3])

  return(out)

}
