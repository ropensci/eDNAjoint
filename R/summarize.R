#' Summarize posterior distributions of model parameters.
#'
#' This function summarizes the posterior distributions of specified parameters from a model fit. Summary includes mean, sd, and specified quantiles, as well as effective sample size (n_eff) and Rhat for estimated parameters.
#'
#' @export
#' @param modelfit An object of class `stanfit`.
#' @param par A character vector of parameter names. The default is 'all'.
#' @param probs A numeric vector of quantiles of interest. The default is c(0.025,0.975).
#' @param digits An integer indicating the number of decimal places to round values in summary table. Default value is 3.
#' @return A summary table of parameter estimates.
#'
#' @note  Before fitting the model, this function checks to ensure that the function is possible given the inputs. These checks include:
#' \itemize{
#' \item  Input model fit is an object of class 'stanfit'.
#' \item  Input probs is a numeric vector.
#' \item  Input par is a character vector.
#' \item  Input par are present in fitted model.
#' \item  Input model fit has converged (i.e. no divergent transitions after warm-up).
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' @examples
#' \donttest{
#' data(greencrabData)
#'
#' # Fit a model
#' modelfit = jointModel(data=greencrabData, family='negbin', q=TRUE)
#'
#' # Create summary table of all parameters
#' summarize(modelfit)
#'
#' # Summarize just 'p10' parameter
#' summarize(modelfit, par = 'p10')
#' }
#'

summarize <- function(modelfit, par = 'all', probs = c(0.025,0.975), digits = 3) {

  ## #1. make sure model fit is of class stanfit
  if(!is(modelfit,'stanfit')) {
    errMsg = paste("modelfit must be of class 'stanfit'.")
    stop(errMsg)
  }

  ## #2. make sure probs is a numeric vector
  if(!is.numeric(probs)) {
    errMsg = paste("probs must be a numeric vector")
    stop(errMsg)
  }

  ## #3. make sure par is a character vector
  if(!is.character(par)) {
    errMsg = paste("par must be a character vector")
    stop(errMsg)
  }

  ## #4. make sure model fit contains all par input
  if(all(!(par %in% modelfit@model_pars)) && par != 'all') {
    errMsg = paste("modelfit must contain all selected parameters: ",par)
    stop(errMsg)
  }

  if (!requireNamespace("rstan", quietly = TRUE)){
    stop ("The 'rstan' package is not installed.", call. = FALSE)
  }

  ## #5. check to see if there are any divergent transitions
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

  if(all(c('p10','q','phi','beta') %in% modelfit@model_pars)==TRUE &&
     !c('alpha') %in% modelfit@model_pars &&
     all(par == 'all')){
    #joint, catchability, negbin
    out <- round(rstan::summary(modelfit,pars=c('p10','beta','q','phi','mu'),probs=probs,use_cache=FALSE)$summary,digits)

  } else if(all(c('p10','beta','q') %in% modelfit@model_pars)==TRUE &&
            all(!c('phi','alpha') %in% modelfit@model_pars) &&
            all(par == 'all')){
    #joint, catchability, pois
    out <- round(rstan::summary(modelfit,pars=c('p10','beta','q','mu'),probs=probs,use_cache=FALSE)$summary,digits)

  } else if(all(c('p10','phi','beta') %in% modelfit@model_pars)==TRUE &&
            all(!c('q','alpha') %in% modelfit@model_pars) &&
            all(par == 'all')){
    #joint, no catchability, negbin
    out <- round(rstan::summary(modelfit,pars=c('p10','beta','phi','mu'),probs=probs,use_cache=FALSE)$summary,digits)

  } else if(all(c('p10','beta') %in% modelfit@model_pars)==TRUE &&
            all(!c('q','phi','alpha') %in% modelfit@model_pars==TRUE) &&
            all(par == 'all')){
    #joint, no catchability, pois
    out <- round(rstan::summary(modelfit,pars=c('p10','beta','mu'),probs=probs,use_cache=FALSE)$summary,digits)

  } else if(all(c('p10','q','phi','alpha') %in% modelfit@model_pars)==TRUE &&
            all(par == 'all')){
    #joint, catchability, negbin, site covariates
    out <- round(rstan::summary(modelfit,pars=c('p10','alpha','beta','q','phi','mu'),probs=probs,use_cache=FALSE)$summary,digits)

  } else if(all(c('p10','alpha','q') %in% modelfit@model_pars)==TRUE &&
            !c('phi') %in% modelfit@model_pars &&
            all(par == 'all')){
    #joint, catchability, pois, site covariates
    out <- round(rstan::summary(modelfit,pars=c('p10','alpha','beta','q','mu'),probs=probs,use_cache=FALSE)$summary,digits)

  } else if(all(c('p10','phi','alpha') %in% modelfit@model_pars)==TRUE &&
            !c('q') %in% modelfit@model_pars &&
            all(par == 'all')){
    #joint, no catchability, negbin, site covariates
    out <- round(rstan::summary(modelfit,pars=c('p10','alpha','beta','phi','mu'),probs=probs,use_cache=FALSE)$summary,digits)

  } else if(all(c('p10','alpha') %in% modelfit@model_pars)==TRUE &&
            all(!c('q','phi') %in% modelfit@model_pars==TRUE) &&
            all(par == 'all')){
    #joint, no catchability, pois, site covariates
    out <- round(rstan::summary(modelfit,pars=c('p10','alpha','beta','mu'),probs=probs,use_cache=FALSE)$summary,digits)

  } else if(all(c('q','phi') %in% modelfit@model_pars)==TRUE &&
            all(!c('p10','beta') %in% modelfit@model_pars==TRUE) &&
            all(par == 'all')){
    #traditional, catchability, negbin
    out <- round(rstan::summary(modelfit,pars=c('q','phi','mu'),probs=probs,use_cache=FALSE)$summary,digits)

  } else if(all(c('q') %in% modelfit@model_pars)==TRUE &&
            all(!c('p10','beta','phi') %in% modelfit@model_pars==TRUE) &&
            all(par == 'all')){
    #traditional, catchability, pois
    out <- round(rstan::summary(modelfit,pars=c('q','mu'),probs=probs,use_cache=FALSE)$summary,digits)

  } else if(all(c('phi') %in% modelfit@model_pars)==TRUE &&
            all(!c('p10','beta','q') %in% modelfit@model_pars==TRUE) &&
            all(par == 'all')){
    #traditional, no catchability, negbin
    out <- round(rstan::summary(modelfit,pars=c('phi','mu'),probs=probs,use_cache=FALSE)$summary,digits)

  } else if(all(!c('p10','beta','q','phi') %in% modelfit@model_pars)==TRUE &&
            all(par == 'all')){
    #traditional, no catchability, pois
    out <- round(rstan::summary(modelfit,pars='mu',probs=probs,use_cache=FALSE)$summary,digits)

  } else {
    out <- round(rstan::summary(modelfit,pars=par,probs=probs,use_cache=FALSE)$summary,digits)
  }


  return(out)

}
