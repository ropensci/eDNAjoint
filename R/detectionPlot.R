#' Plot the survey effort necessary to detect species presence, given the
#' species expected catch rate.
#'
#' This function plots the median number of survey effort units to necessary
#' detect species presence. Detecting species presence is defined as producing
#' at least one true positive eDNA detection or catching at least one
#' individual. See more examples in the
#' \href{https://bookdown.org/abigailkeller/eDNAjoint_vignette/}{Package
#' Vignette}.
#'
#' @srrstats {G1.4} Roxygen function documentation begins here
#' @export
#' @srrstats {G2.1a} Here are explicit documentation of vector input types
#' @param modelfit An object of class `stanfit`.
#' @param mu.min A value indicating the minimum expected species catch rate for
#'   plotting. If multiple traditional gear types are represented in the model,
#'   mu is the catch rate of gear type 1.
#' @param mu.max A value indicating the maximum expected species catch rate for
#'   plotting. If multiple traditional gear types are represented in the model,
#'   mu is the catch rate of gear type 1.
#' @param cov.val A numeric vector indicating the values of site-level
#'   covariates to use for prediction. Default is NULL.
#' @param qPCR.N An integer indicating the number of qPCR replicates per eDNA
#'   sample. The default is 3.
#' @param probability A numeric value indicating the probability of detecting
#'   presence. The default is 0.9.
#' @return A plot displaying survey efforts necessary to detect species
#'   presence, given mu, for each survey type.
#'
#' @srrstats {G2.0a,G2.2} Explicit secondary documentation of any expectations
#'   on lengths of inputs
#' @note  Before fitting the model, this function checks to ensure that the
#'   function is possible given the inputs. These checks include:
#' \itemize{
#' \item  Input model fit is an object of class 'stanfit'.
#' \item  Input mu.min is a numeric value greater than 0.
#' \item  Input mu.max is a numeric value.
#' \item  If model fit contains alpha, cov.val must be provided.
#' \item  Input cov.val is numeric.
#' \item  Input cov.val is the same length as the number of estimated
#'   covariates.
#' \item  Input probability is a univariate numeric value.
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
#' # Plot at the mean covariate values (covariates are standardized, so mean=0)
#' detectionPlot(fit.cov$model, mu.min = 0.1, mu.max = 1,
#'               cov.val = c(0,0), qPCR.N = 3)
#'
#' # Calculate mu_critical at salinity 0.5 z-scores greater than the mean
#' detectionPlot(fit.cov$model, mu.min = 0.1, mu.max = 1, cov.val = c(0,0.5),
#'               qPCR.N = 3)
#'
#' # Ex. 2: Calculating necessary effort for detection with multiple
#' # traditional gear types
#'
#' # Load data
#' data(greencrabData)
#'
#' # Fit a model with no site-level covariates
#' fit.q <- jointModel(data = greencrabData, cov = NULL, family = "negbin",
#'                     p10priors = c(1,20), q = TRUE,
#'                     multicore = FALSE)
#'
#' # Calculate
#' detectionPlot(fit.q$model, mu.min = 0.1, mu.max = 1,
#'               cov.val = NULL, qPCR.N = 3)
#'
#' # Change probability of detecting presence to 0.95
#' detectionPlot(fit.q$model, mu.min = 0.1, mu.max = 1, cov.val = NULL,
#'               probability = 0.95, qPCR.N = 3)
#' }
#'

detectionPlot <- function(modelfit, mu.min, mu.max, cov.val = NULL,
                          probability = 0.9, qPCR.N = 3){

  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using this helper
  #'   function
  detectionPlot_input_checks(modelfit, mu.min, mu.max, cov.val,
                             probability, qPCR.N)

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

  # create mu vector
  mu <- seq(from = mu.min, to = mu.max, by = (mu.max-mu.min)/1000)

  # get n samples df
  out <- detectionCalculate(modelfit, mu, cov.val, probability, qPCR.N)

  # convert to long df
  '%>%' <- magrittr::`%>%`
  out_long <- as.data.frame(out) %>%
    tidyr::pivot_longer(cols =! mu, names_to = 'survey_type')

  # get full names
  if(length(unique(out_long$survey_type))==1){
    names <- 'traditional'
  } else if(length(unique(out_long$survey_type))==2){
    names <- c('eDNA','traditional')
  } else{
    names <- 'eDNA'
    for(i in 1:(length(unique(out_long$survey_type))-1)){
      names <- c(names, paste0('traditional (gear type ',i,')'))
    }
  }

  plot <- ggplot2::ggplot(data = out_long)+
    ggplot2::geom_line(ggplot2::aes(x = mu, y = value, color = survey_type),
                       linewidth = 1)+
    ggplot2::labs(x = 'mu (expected catch rate)',y = '# survey units',
                  color = 'survey type')+
    ggplot2::scale_color_manual(values = scales::hue_pal()(length(names)),
                                labels = names)+
    ggplot2::theme_minimal()

  return(plot)
}

#function to check for divergent transitions
div_check <- function(x){
  divergent <- sum(x[,'divergent__'])
  return(divergent)
}

# function for input checks
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#'   messages
detectionPlot_input_checks <- function(modelfit, mu.min, mu.max, cov.val,
                                       probability, qPCR.N){
  ## #1. make sure model fit is of class stanfit
  #' @srrstats {G2.8} Makes sure input of sub-function is of class 'stanfit'
  #'   (i.e., output of jointModel())
  if(!is(modelfit,'stanfit')) {
    errMsg <- "modelfit must be of class 'stanfit'."
    stop(errMsg)
  }

  ## #2. make sure mu.min is a numeric value
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(!is.numeric(mu.min) | length(mu.min)>1 | mu.min <= 0) {
    errMsg <- "mu.min must be a numeric value greater than 0"
    stop(errMsg)
  }

  ## #3. make sure mu.max is a numeric value
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(!is.numeric(mu.max) | length(mu.max)>1 | mu.max <= mu.min) {
    errMsg <- "mu.max must be a numeric value greater than mu.min"
    stop(errMsg)
  }

  ## #4. make sure probability is a numeric value between 0 and 1
  #' @srrstats {G2.0} Assertion on length of input, prohibit multivariate input
  #'   to parameters expected to be univariate
  if(!is.numeric(probability) | length(probability)>1 | probability < 0 |
     probability > 1) {
    errMsg <- "probability must be a numeric value between 0 and 1"
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
  if(all(!is.null(cov.val)) && length(cov.val)!=(modelfit@par_dims$alpha-1)) {
    errMsg <- paste0("cov.val must be of the same length as the number of ",
                     "non-intercept site-level coefficients in the model.")
    stop(errMsg)
  }

  ## #8. If covariates are in model, cov.val must be provided
  if(all(c('alpha','p10') %in% modelfit@model_pars) && all(is.null(cov.val))) {
    errMsg <- paste0("cov.val must be provided if the model contains ",
                     "site-level covariates.")
    stop(errMsg)
  }

  ## #9. qPCR.N must be an integer
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  if(qPCR.N %% 1 != 0) {
    errMsg <- "qPCR.N should be an integer."
    stop(errMsg)
  }
}


