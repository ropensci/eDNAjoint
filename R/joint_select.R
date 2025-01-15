#' Perform model selection using leave one out cross validation of model objects
#'
#' This function performs leave one out cross validation of a list of model
#' fits using functions in the `loo` package, as described in Vehtari, Gelman,
#' and Gabry (2017) <doi:10.1007/s11222-016-9696-4>. Compare models fit using
#' `joint_model()` or models fits using `traditional_model()`. See more examples
#' in the \href{https://ednajoint.netlify.app}{Package
#' Vignette}.
#'
#' @srrstats {G1.0} The literature reference for leave one out cross validation
#'   is provided here.
#'
#' @srrstats {G1.4} Roxygen function documentation begins here
#' @export
#' @param model_fits A list containing model fits of class `stanfit`.
#' @return A matrix of delta elpd (expected log pointwise predictive density)
#'   between model fits. Function is performed using the `loo` package.
#'
#' @note  Before model selection, this function makes the following check:
#' \itemize{
#' \item Input is a list of model fits of class 'stanfit'.
#' \item All models compared were fit wither either `joint_model()` or all with
#'   `traditional_model().`
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' @examples
#' \donttest{
#' data(green_crab_data)
#'
#' # Fit a model without estimating a gear scaling coefficient for traditional
#' # survey gear types.
#' # This model assumes all traditional survey methods have the same
#' # catchability.
#' # Count data is modeled using a poisson distribution.
#' fit_no_q <- joint_model(data = green_crab_data, family = "poisson",
#'                         p10_priors = c(1,20), q = FALSE, multicore = FALSE)
#'
#'
#' # Fit a model estimating a gear scaling coefficient for traditional
#' # survey gear types.
#' # This model does not assume all traditional survey methods have the
#' # same catchability.
#' # Gear type 1 is used as the reference gear type.
#' # Count data is modeled using a negative binomial distribution.
#' fit_q <- joint_model(data = green_crab_data, family = "negbin",
#'                      p10_priors = c(1,20), q = TRUE, multicore = FALSE)
#'
#' # Perform model selection
#' joint_select(model_fits = list(fit_no_q$model, fit_q$model))
#' }
#'

joint_select <- function(model_fits) {

  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using this
  #'   helper function
  joint_select_input_checks(model_fits)


  #create empty loo list
  loo_list <- list()

  ##leave one out cross-validation
  for (i in seq_along(model_fits)) {
    name <- paste("model", i, sep = "")
    loo_list[[name]] <- loo::loo(model_fits[[i]])
  }

  #compare loo
  loo_out <- loo::loo_compare(loo_list)

  return(loo_out)
}

# function for input checks
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#'   messages
joint_select_input_checks <- function(model_fits) {
  ## #1. make sure input is a list
  if (!is(model_fits, "list")) {
    err_msg <- "model_fits must be a list."
    stop(err_msg)
  }

  ## #2. make sure all data objects are of class stanfit
  #' @srrstats {G2.8} Makes sure input of sub-function is of class 'stanfit'
  #'   (i.e., output of joint_model())
  test <- function(model_fit) {
    is(model_fit, "stanfit")
  }
  if (!all(lapply(model_fits, test) == TRUE)) {
    err_msg <- "Model fits in model_fits input must be of class 'stanfit'."
    stop(err_msg)
  }

  ## #3. make sure all models are of the same type
  type <- function(model_fit) {
    "p10" %in% model_fit@model_pars
  }
  if (length(unique(lapply(model_fits, type))) != 1) {
    err_msg <- paste0("All model_fits must be fit with either joint_model() ",
                      "or all with traditional_model().")
    stop(err_msg)
  }
}
