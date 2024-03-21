#' Perform model selection using leave one out cross validation of model objects
#'
#' @srrstats {G1.0} The literature reference for leave one out cross validation
#' is provided here.
#' This function performs leave one out cross validation of a list of model
#' fits using functions in the `loo` package, as described in Vehtari, Gelman,
#' and Gabry (2017) <doi:10.1007/s11222-016-9696-4>. Compare models fit using
#' `jointModel()` or models fits using `traditionalModel()`.
#'
#' @srrstats {G1.4} Roxygen function documentation begins here
#' @export
#' @param modelfits A list containing model fits of class `stanfit`.
#' @return A matrix of delta elpd (expected log pointwise predictive density)
#' between model fits. Function is performed using the `loo` package.
#'
#' @note  Before model selection, this function makes the following check:
#' \itemize{
#' \item Input is a list of model fits of class 'stanfit'.
#' \item All models compared were fit wither either `jointModel()` or all with
#' `traditionalModel().`
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' @examples
#' \donttest{
#' data(greencrabData)
#'
#' # Fit a model without estimating a catchability coefficient for traditional
#' # survey gear types.
#' # This model assumes all traditional survey methods have the same
#' # catchability.
#' # Count data is modeled using a poisson distribution.
#' fit.no.q = jointModel(data=greencrabData, family='poisson',
#'                       p10priors=c(1,20), q=FALSE)
#'
#'
#' # Fit a model estimating a catchability coefficient for traditional
#' # survey gear types.
#' # This model does not assume all traditional survey methods have the
#' # same catchability.
#' # Gear type 1 is used as the reference gear type.
#' # Count data is modeled using a negative binomial distribution.
#' fit.q = jointModel(data=greencrabData, family='negbin',
#'                    p10priors=c(1,20), q=TRUE)
#'
#' # Perform model selection
#' jointSelect(modelfits=list(fit.no.q$model, fit.q$model))
#' }
#'

jointSelect <- function(modelfits) {

  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using this
  #' helper function
  jointSelect_input_checks(modelfits)


  #create empty loo list
  loo_list <- list()

  ##leave one out cross-validation
  for(i in seq_along(modelfits)){
    name <- paste('model',i,sep='')
    loo_list[[name]] <- loo::loo(modelfits[[i]])
  }

  #compare loo
  loo_out <- loo::loo_compare(loo_list)

  return(loo_out)
}

# function for input checks
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#' messages
jointSelect_input_checks <- function(modelfits){
  ## #1. make sure input is a list
  if (!is(modelfits,'list')){
    errMsg = "modelfits must be a list."
    stop(errMsg)
  }

  ## #2. make sure all data objects are of class stanfit
  #' @srrstats {G2.8} Makes sure input of sub-function is of class 'stanfit'
  #' (i.e., output of jointModel())
  test <- function(modelfit){is(modelfit,'stanfit')}
  if (!all(lapply(modelfits,test)==TRUE)) {
    errMsg <- "Model fits in modelfits input must be of class 'stanfit'."
    stop(errMsg)
  }

  ## #3. make sure all models are of the same type
  type <- function(modelfit){'p10'%in%modelfit@model_pars}
  if (length(unique(lapply(modelfits,type)))!=1){
    errMsg <- paste0("All modelfits must be fit with either jointModel() or ",
                     "all with traditionalModel().")
    stop(errMsg)
  }
}

