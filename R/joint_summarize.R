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
#' @param model_fit An object of class `stanfit`.
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
#' data(green_crab_data)
#'
#' # Fit a model
#' model_fit <- joint_model(data = green_crab_data, family = "negbin", q = TRUE,
#'                          multicore = FALSE)
#'
#' # Create summary table of all parameters
#' joint_summarize(model_fit$model)
#'
#' # Summarize just 'p10' parameter
#' joint_summarize(model_fit$model, par = "p10", probs = c(0.025, 0.975),
#'                 digits = 3)
#' }
#'

joint_summarize <- function(model_fit, par = "all", probs = c(0.025, 0.975),
                            digits = 3) {

  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using this helper
  #'   function
  joint_summarize_input_checks(model_fit, par, probs)

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("The 'rstan' package is not installed.", call. = FALSE)
  }

  ## check to see if there are any divergent transitions
  #' @srrstats {BS4.5} Warning message if the input model fit has
  #'   divergence transitions
  if (sum(lapply(rstan::get_sampler_params(model_fit,
                                           inc_warmup = FALSE),
                 div_check)[[1]]) > 0) {

    sum <- sum(lapply(rstan::get_sampler_params(model_fit,
                                                inc_warmup = FALSE),
                      div_check)[[1]])

    warning_msg <- paste0("Warning: There are ", sum,
                          " divergent transitions in your model fit. ")
    warning(warning_msg)

  }

  # get summary
  if (all(par == "all")) {

    params <- get_all_params(model_fit@model_pars, model_fit)
    #' @srrstats {G2.4,G2.4a} explicit conversion to integers for sampling
    #'   arguments
    out <- round(rstan::summary(model_fit, pars = params, probs = probs,
                                use_cache = FALSE)$summary, as.integer(digits))
  } else {
    #' @srrstats {G2.4,G2.4a} explicit conversion to integers for sampling
    #'   arguments
    out <- round(rstan::summary(model_fit, pars = par, probs = probs,
                                use_cache = FALSE)$summary, as.integer(digits))
  }

  # fix row name if phi present
  if ("phi[1]" %in% rownames(out)) {
    row_index <- which(rownames(out) == "phi[1]")
    rownames(out)[row_index] <- "phi"
  }

  # fix row names if q = FALSE
  if (model_fit@par_dims$q == 0) {
    names <- rownames(out)
    rownames(out) <- gsub("\\[([0-9]+),1\\]", "[\\1]", names)
  }

  return(out)

}

# functions to check model type
is_joint <- function(pars) {
  out <- ifelse("p10" %in% pars, TRUE, FALSE)
  return(out)
}
is_catch <- function(model_fit) {
  out <- ifelse(model_fit@par_dims$q > 0, TRUE, FALSE)
  return(out)
}
is_negbin <- function(model_fit) {
  if ("phi" %in% model_fit@model_pars) {
    out <- ifelse(model_fit@par_dims$phi == 1, TRUE, FALSE)
  } else {
    out <- FALSE
  }
  return(out)
}

# function to get vector of all param names
get_all_params <- function(pars, model_fit) {
  params <- c("mu")

  # catchability
  if (is_catch(model_fit)) {
    params <- c(params, "q")
  }

  # joint
  if (is_joint(pars)) {
    params <- c(params, "p10", "beta", "alpha")
  }

  # negbin
  if (is_negbin(model_fit)) {
    params <- c(params, "phi[1]")
  }

  return(params)
}

# function to check for divergent transitions
div_check <- function(x) {
  divergent <- sum(x[, "divergent__"])
  return(divergent)
}

# function for input checks
#' @srrstats {G5.2a} Pre-processing routines to check inputs have
#'   unique messages
joint_summarize_input_checks <- function(model_fit, par, probs) {
  ## #1. make sure model fit is of class stanfit
  #' @srrstats {G2.8} Makes sure input of sub-function is of class 'stanfit'
  #'   (i.e., output of joint_model())
  if (!is(model_fit, "stanfit")) {
    err_msg <- "model_fit must be of class 'stanfit'."
    stop(err_msg)
  }

  ## #2. make sure probs is a numeric vector
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if (!is.numeric(probs)) {
    err_msg <- "probs must be a numeric vector."
    stop(err_msg)
  }

  ## #3. make sure all values of probs are between 0 and 1
  if (any(probs > 1) || any(probs < 0)) {
    err_msg <- "probs must be between 0 and 1."
    stop(err_msg)
  }

  ## #4. make sure par is a character vector
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if (!is.character(par)) {
    err_msg <- "par must be a character vector."
    stop(err_msg)
  }

  ## #5. make sure model fit contains all par input
  if (any(!(par %in% model_fit@model_pars)) && all(par != "all")) {
    err_msg <- paste("model_fit must contain all selected parameters:", par)
    stop(err_msg)
  }
}
