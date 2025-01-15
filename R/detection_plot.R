#' Plot the survey effort necessary to detect species presence, given the
#' species expected catch rate.
#'
#' This function plots the number of survey effort units to necessary detect
#' species presence, calculated using median estimated parameter values from
#' joint_model(). Detecting species presence is defined as producing
#' at least one true positive eDNA detection or catching at least one
#' individual. See more examples in the
#' \href{https://ednajoint.netlify.app/}{Package
#' Vignette}.
#'
#' @srrstats {G1.4} Roxygen function documentation begins here
#' @export
#' @srrstats {G2.1a} Here are explicit documentation of vector input types
#' @param model_fit An object of class `stanfit`.
#' @param mu_min A value indicating the minimum expected species catch rate for
#'   plotting. If multiple traditional gear types are represented in the model,
#'   mu is the catch rate of gear type 1.
#' @param mu_max A value indicating the maximum expected species catch rate for
#'   plotting. If multiple traditional gear types are represented in the model,
#'   mu is the catch rate of gear type 1.
#' @param cov_val A numeric vector indicating the values of site-level
#'   covariates to use for prediction. Default is NULL.
#' @param pcr_n An integer indicating the number of qPCR replicates per eDNA
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
#' \item  Input mu_min is a numeric value greater than 0.
#' \item  Input mu_max is a numeric value.
#' \item  If model fit contains alpha, cov_val must be provided.
#' \item  Input cov_val is numeric.
#' \item  Input cov_val is the same length as the number of estimated
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
#' data(goby_data)
#'
#' # Fit a model including 'Filter_time' and 'Salinity' site-level covariates
#' fit_cov <- joint_model(data = goby_data, cov = c('Filter_time','Salinity'),
#'                        family = "poisson", p10_priors = c(1,20), q = FALSE,
#'                        multicore = FALSE)
#'
#' # Plot at the mean covariate values (covariates are standardized, so mean=0)
#' detection_plot(fit_cov$model, mu_min = 0.1, mu_max = 1,
#'                cov_val = c(0,0), pcr_n = 3)
#'
#' # Calculate mu_critical at salinity 0.5 z-scores greater than the mean
#' detection_plot(fit_cov$model, mu_min = 0.1, mu_max = 1, cov_val = c(0,0.5),
#'                pcr_n = 3)
#'
#' # Ex. 2: Calculating necessary effort for detection with multiple
#' # traditional gear types
#'
#' # Load data
#' data(green_crab_data)
#'
#' # Fit a model with no site-level covariates
#' fit_q <- joint_model(data = green_crab_data, cov = NULL, family = "negbin",
#'                      p10_priors = c(1,20), q = TRUE,
#'                      multicore = FALSE)
#'
#' # Calculate
#' detection_plot(fit_q$model, mu_min = 0.1, mu_max = 1,
#'                cov_val = NULL, pcr_n = 3)
#'
#' # Change probability of detecting presence to 0.95
#' detection_plot(fit_q$model, mu_min = 0.1, mu_max = 1, cov_val = NULL,
#'                probability = 0.95, pcr_n = 3)
#' }
#'

detection_plot <- function(model_fit, mu_min, mu_max, cov_val = NULL,
                           probability = 0.9, pcr_n = 3) {

  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using this helper
  #'   function
  det_plot_input_checks_1(model_fit, mu_min, mu_max, cov_val,
                          probability, pcr_n)
  det_plot_input_checks_2(model_fit, mu_min, mu_max, cov_val,
                          probability, pcr_n)

  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("The 'rstan' package is not installed.", call. = FALSE)
  }

  ## check to see if there are any divergent transitions
  #' @srrstats {BS4.5} Warning message if the input model fit has divergence
  #'   transitions
  if (sum(lapply(rstan::get_sampler_params(model_fit, inc_warmup = FALSE),
                 div_check)[[1]]) > 0) {

    sum <- sum(lapply(rstan::get_sampler_params(model_fit, inc_warmup = FALSE),
                      div_check)[[1]])

    warning_msg <- paste0("Warning: There are ", sum,
                          " divergent transitions in your model fit. ")
    warning(warning_msg)

  }

  # create mu vector
  mu <- seq(from = mu_min, to = mu_max, by = (mu_max - mu_min) / 1000)

  # get n samples df
  out <- detection_calculate(model_fit, mu, cov_val, probability, pcr_n)

  # convert to long df
  out_long <- as.data.frame(out) |>
    tidyr::pivot_longer(cols = ! mu, names_to = "survey_type")

  # get full names
  if (length(unique(out_long$survey_type)) == 1) {
    names <- "traditional"
  } else if (length(unique(out_long$survey_type)) == 2) {
    names <- c("eDNA", "traditional")
  } else {
    names <- "eDNA"
    for (i in 1:(length(unique(out_long$survey_type)) - 1)) {
      names <- c(names, paste0("traditional (gear type ", i, ")"))
    }
  }

  plot <- ggplot2::ggplot(data = out_long) +
    ggplot2::geom_line(ggplot2::aes(x = mu, y = value, color = survey_type),
                       linewidth = 1) +
    ggplot2::labs(x = "mu (expected catch rate)", y = "# survey units",
                  color = "survey type") +
    ggplot2::scale_color_manual(values = scales::hue_pal()(length(names)),
                                labels = names) +
    ggplot2::theme_minimal()

  return(plot)
}

#function to check for divergent transitions
div_check <- function(x) {
  divergent <- sum(x[, "divergent__"])
  return(divergent)
}

# function for input checks
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#'   messages
#'
# Helper function to validate conditions
validate_condition <- function(condition, message) {
  if (!condition) stop(message)
}

det_plot_input_checks_1 <- function(model_fit, mu_min, mu_max, cov_val,
                                    probability, pcr_n) {
  ## #1. make sure model fit is of class stanfit
  #' @srrstats {G2.8} Makes sure input of sub-function is of class 'stanfit'
  #'   (i.e., output of joint_model())
  validate_condition(!is(model_fit, "stanfit"),
                     "model_fit must be of class 'stanfit'.")

  ## #2. make sure mu_min is a numeric value
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  validate_condition(!is.numeric(mu_min) || length(mu_min) > 1 || mu_min <= 0,
                     "mu_min must be a numeric value greater than 0")

  ## #3. make sure mu_max is a numeric value
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  validate_condition(!is.numeric(mu_max) || length(mu_max) > 1 ||
                       mu_max <= mu_min,
                     "mu_max must be a numeric value greater than mu_min")

  ## #4. make sure probability is a numeric value between 0 and 1
  #' @srrstats {G2.0} Assertion on length of input, prohibit multivariate input
  #'   to parameters expected to be univariate
  validate_condition(!is.numeric(probability) || length(probability) > 1 ||
                       probability < 0 || probability > 1,
                     "probability must be a numeric value between 0 and 1")

  ## #5. cov_val is numeric, if provided
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  validate_condition(all(!is.null(cov_val)) && !is.numeric(cov_val),
                     "cov_val must be a numeric vector")

}

det_plot_input_checks_2 <- function(model_fit, mu_min, mu_max, cov_val,
                                    probability, pcr_n) {

  ## #6. Only include input cov_val if covariates are included in model
  validate_condition(all(!is.null(cov_val)) &&
                       model_fit@par_dims$alpha == 1,
                     paste0("cov_val must be NULL if the model does not ",
                            "contain site-level covariates."))

  ## #7. Input cov_val is the same length as the number of estimated covariates.
  #' @srrstats {G2.0} Assertion on length of input
  validate_condition(all(!is.null(cov_val)) &&
                       length(cov_val) != (model_fit@par_dims$alpha - 1),
                     paste0("cov_val must be of the same length as the number ",
                            "of non-intercept site-level coefficients in the ",
                            "model."))

  ## #8. If covariates are in model, cov_val must be provided
  validate_condition(model_fit@par_dims$alpha > 1 &&
                       "p10" %in% model_fit@model_pars && all(is.null(cov_val)),
                     paste0("cov_val must be provided if the model contains ",
                            "site-level covariates."))

  ## #9. pcr_n must be an integer
  #' @srrstats {G5.8,G5.8b} Pre-processing routines to check for data of
  #'   unsupported type
  validate_condition(pcr_n %% 1 != 0, "pcr_n should be an integer.")
}
