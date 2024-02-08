#' Perform model selection using leave-one-out cross validation of model objects
#'
#' This function performs leave-one-out cross validation of a list of model fits using functions in the `loo` package. Compare models fit using `jointModel()` or models fits using `traditionalModel()`.
#'
#' @export
#' @param modelfits A list containing model fits of class `stanfit`.
#' @return A matrix of delta elpd (expected log pointwise predictive density) between model fits. Function is performed using the `loo` package.
#'
#' @note  Before model selection, this function makes the following check:
#' \itemize{
#' \item Input is a list of model fits of class 'stanfit'.
#' \item All models compared were fit wither either `jointModel()` or all with `traditionalModel().`
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' @examples
#' \donttest{
#' data(greencrabData)
#'
# Fit a model without estimating a catchability coefficient for traditional survey gear types.
#' # This model assumes all traditional survey methods have the same catchability.
#' # Count data is modeled using a poisson distribution.
#' fit.no.q = jointModel(data=greencrabData, family='poisson', p10priors=c(1,20), q=FALSE)
#'
#'
#' # Fit a model estimating a catchability coefficient for traditional survey gear types.
#' # This model does not assume all traditional survey methods have the same catchability.
#' # Gear type 1 is used as the reference gear type.
#' # Count data is modeled using a negative binomial distribution.
#' fit.q = jointModel(data=greencrabData, family='negbin', p10priors=c(1,20), q=TRUE)
#'
#' # Perform model selection
#' jointSelect(modelfits=list(fit.no.q, fit.q))
#' }
#'

jointSelect <- function(modelfits) {

  ## #1. make sure input is a list
  if (!is(modelfits,'list')){
      errMsg = paste("modelfits must be a list.")
      stop(errMsg)
    }

  ## #2. make sure all data objects are of class stanfit
  test <- function(modelfit){is(modelfit,'stanfit')}
  if (!all(lapply(modelfits,test)==TRUE)) {
      errMsg = paste("Model fits in modelfits input must be of class 'stanfit'.")
      stop(errMsg)
  }

  ## #3. make sure all models are of the same type
  type <- function(modelfit){'p10'%in%modelfit@model_pars}
  if (length(unique(lapply(modelfits,type)))!=1){
    errMsg = paste("All modelfits must be fit with either jointModel() or all with traditionalModel().")
    stop(errMsg)
  }

  #create empty loo list
  loo_list <- list()

  ##leave one out cross-validation
  for(i in 1:length(modelfits)){
    name <- paste('model',i,sep='')
    loo_list[[name]] <- loo::loo(modelfits[[i]])
  }

  #compare loo
  loo_out <- loo::loo_compare(loo_list)

  return(loo_out)
}

