#' Summarize posterior distributions of model parameters.
#'
#' This function summarizes the posterior distributions of specified parameters
#' from a model fit. Summary includes mean, sd, and specified quantiles, as
#' well as effective sample size (n_eff) and Rhat for estimated parameters. See
#' more examples in the
#' \href{https://ednajoint.netlify.app}{Package
#' Vignette}.
#'
#' @srrstats {BS5.3,BS6.4} Function to summarize parameter samples and provide
#'   convergence statistics (using Stan functions), used to summarize parameter
#'   values relevant to the user
#'
#' @srrstats {G1.4} Roxygen function documentation begins here
#' @export
#' @srrstats {G2.1a} Here are explicit documentation of vector input types
#' @param modelfit An object of class `stanfit`.
#' @param par A character vector of parameter names. The default is 'all'.
#' @param probs A numeric vector of quantiles of interest. The default is
#'   c(0.025,0.975).
#' @param digits An integer indicating the number of decimal
#'   places to round values in summary table. Default value is 3.
#' @return A summary table of parameter estimates.
#'
#' @note  Before fitting the model, this function checks to ensure that the
#'   function is possible given the inputs. These checks include:
#' \itemize{
#' \item  Input model fit is an object of class 'stanfit'.
#' \item  Input probs is a numeric vector.
#' \item  Input par is a character vector.
#' \item  Input par are present in fitted model.
#' \item  Input model fit has converged (i.e. no divergent transitions after
#'   warm-up).
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' @examples
#' \donttest{
#' data(greencrabData)
#'
#' # Fit a model
#' modelfit <- jointModel(data = greencrabData, family = "negbin", q = TRUE,
#'                        multicore = FALSE)
#'
#' # Create summary table of all parameters
#' jointSummarize(modelfit$model)
#'
#' # Summarize just 'p10' parameter
#' jointSummarize(modelfit$model, par = "p10", probs = c(0.025, 0.975),
#'                digits = 3)
#' }
#'

jointSummarize <- function(modelfit, par = 'all', probs = c(0.025,0.975),
                           digits = 3) {

  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using this helper
  #'   function
  jointSummarize_input_checks(modelfit, par, probs)

  if (!requireNamespace("rstan", quietly = TRUE)){
    stop ("The 'rstan' package is not installed.", call. = FALSE)
  }

  ## check to see if there are any divergent transitions
  #' @srrstats {BS4.5} Warning message if the input model fit has
  #'   divergence transitions
  if(sum(lapply(rstan::get_sampler_params(modelfit,
                                          inc_warmup = FALSE),
                div_check)[[1]]) > 0){

    sum <- sum(lapply(rstan::get_sampler_params(modelfit,
                                                inc_warmup = FALSE),
                      div_check)[[1]])

    warning_msg <- paste0('Warning: There are ',sum,
                          ' divergent transitions in your model fit. ')
    warning(warning_msg)

  }

  # get summary
  if(all(par == 'all')){

    params <- get_all_params(modelfit@model_pars)
    #' @srrstats {G2.4,G2.4a} explicit conversion to integers for sampling
    #'   arguments
    out <- round(rstan::summary(modelfit, pars = params, probs = probs,
                                use_cache = FALSE)$summary, as.integer(digits))
  } else {
    #' @srrstats {G2.4,G2.4a} explicit conversion to integers for sampling
    #'   arguments
    out <- round(rstan::summary(modelfit, pars = par, probs = probs,
                                use_cache = FALSE)$summary, as.integer(digits))
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

# function to get vector of all param names
get_all_params <- function(pars){
  params <- c('mu')

  # catchability
  if(isCatch(pars)){
    params <- c(params,'q')
  }

  # joint
  if(isJoint(pars)){
    params <- c(params,'p10','beta','alpha')
  }

  # negbin
  if(isNegbin(pars)){
    params <- c(params,'phi')
  }

  return(params)
}

# function to check for divergent transitions
div_check <- function(x){
  divergent <- sum(x[,'divergent__'])
  return(divergent)
}

# function for input checks
#' @srrstats {G5.2a} Pre-processing routines to check inputs have
#'   unique messages
jointSummarize_input_checks <- function(modelfit, par, probs){
  ## #1. make sure model fit is of class stanfit
  #' @srrstats {G2.8} Makes sure input of sub-function is of class 'stanfit'
  #'   (i.e., output of jointModel())
  if(!is(modelfit,'stanfit')) {
    errMsg <- "modelfit must be of class 'stanfit'."
    stop(errMsg)
  }

  ## #2. make sure probs is a numeric vector
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(!is.numeric(probs)) {
    errMsg <- "probs must be a numeric vector."
    stop(errMsg)
  }

  ## #3. make sure all values of probs are between 0 and 1
  if(any(probs > 1 | probs < 0)) {
    errMsg <- "probs must be between 0 and 1."
    stop(errMsg)
  }

  ## #4. make sure par is a character vector
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(!is.character(par)) {
    errMsg <- "par must be a character vector."
    stop(errMsg)
  }

  ## #5. make sure model fit contains all par input
  if(all(!(par %in% modelfit@model_pars)) && par != 'all') {
    errMsg <- paste("modelfit must contain all selected parameters:",par)
    stop(errMsg)
  }
}
