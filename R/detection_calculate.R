#' Calculate the survey effort necessary to detect species presence, given the
#' species expected catch rate.
#'
#' This function calculates the number of survey effort units to necessary
#' detect species presence using median estimated parameter values from
#' joint_model(). Detecting species presence is defined
#' as producing at least one true positive eDNA detection or catching at least
#' one individual. See more examples in the
#' \href{https://ednajoint.netlify.app/}{Package
#' Vignette}.
#'
#' @srrstats {G1.4} Roxygen function documentation begins here
#' @export
#' @srrstats {G2.1a} Here are explicit documentation of vector input types
#' @param model_fit An object of class `stanfit`.
#' @param mu A numeric vector of species densities/capture rates. If multiple
#'   traditional gear types are represented in the model, mu is the catch rate
#'   of gear type 1.
#' @param cov_val A numeric vector indicating the values of site-level
#'   covariates to use for prediction. Default is NULL.
#' @param probability A numeric value indicating the probability of detecting
#'   presence. The default is 0.9.
#' @param pcr_n An integer indicating the number of qPCR replicates per eDNA
#'   sample. The default is 3.
#' @return A summary table of survey efforts necessary to detect species
#'   presence, given mu, for each survey type.
#'
#' @srrstats {G2.0a,G2.2} Explicit secondary documentation of any expectations
#'   on lengths of inputs
#' @note  Before fitting the model, this function checks to ensure that the
#'   function is possible given the inputs. These checks include:
#' \itemize{
#' \item  Input model fit is an object of class 'stanfit'.
#' \item  Input mu is a numeric vector.
#' \item  Input probability is a univariate numeric value.
#' \item  If model fit contains alpha, cov_val must be provided.
#' \item  Input cov_val is numeric.
#' \item  Input cov_val is the same length as the number of estimated
#'   covariates.
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
#' # Calculate at the mean covariate values
#' # (covariates are standardized, so mean = 0)
#' detection_calculate(fit_cov$model, mu = seq(from = 0.1, to = 1, by = 0.1),
#'                     cov_val = c(0,0), pcr_n = 3)
#'
#' # Calculate mu_critical at salinity 0.5 z-scores greater than the mean
#' detection_calculate(fit.cov$model, mu = seq(from = 0.1, to = 1, by = 0.1),
#'                     cov_val = c(0,0.5), pcr_n = 3)
#'
#' # Ex. 2: Calculating necessary effort for detection with multiple traditional
#' # gear types
#'
#' # Load data
#' data(green_crab_data)
#'
#' # Fit a model with no site-level covariates
#' fit_q <- joint_model(data = green_crab_data, cov = NULL, family = "negbin",
#'                      p10_priors = c(1,20), q = TRUE, multicore = FALSE)
#'
#' # Calculate
#' detection_calculate(fit_q$model, mu = seq(from = 0.1, to = 1, by = 0.1),
#'                     cov_val = NULL, pcr_n = 3)
#'
#' # Change probability of detecting presence to 0.95
#' detection_calculate(fit_q$model, mu = 0.1, cov_val = NULL,
#'                     probability = 0.95, pcr_n = 3)
#' }
#'

detection_calculate <- function(model_fit, mu, cov_val = NULL,
                                probability = 0.9, pcr_n = 3) {

  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using this helper
  #'   function
  det_calc_input_checks_1(model_fit, mu, cov_val, probability, pcr_n)
  det_calc_input_checks_2(model_fit, mu, cov_val, probability, pcr_n)

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

  # get n traditional samples
  if (is_catch(model_fit)) {
    ntrad_out <- get_ntrad_q(model_fit@model_pars, model_fit, mu, probability)
  } else {
    ntrad_out <- get_ntrad(model_fit@model_pars, model_fit, mu, probability)
  }

  # get n dna samples
  if (is_joint(model_fit@model_pars)) {
    ndna_out <- get_ndna(model_fit@model_pars, model_fit, mu, pcr_n,
                         probability, cov_val)

  }

  # combine into one df - joint model
  if (is_joint(model_fit@model_pars)) {
    if (is_catch(model_fit)) {
      out <- cbind(mu, ntrad_out, ndna_out)
      #rename columns
      trad_names <- c()
      for (i in 1:model_fit@par_dims$q) {
        trad_names <- c(trad_names, paste("n_traditional_", i + 1, sep = ""))
      }
      colnames(out) <- c("mu", "n_traditional_1", trad_names, "n_eDNA")
    } else {
      out <- cbind(mu, ntrad_out, ndna_out)
      colnames(out) <- c("mu", "n_traditional", "n_eDNA")
    }
  }

  # combine into one df - traditional model
  if (!is_joint(model_fit@model_pars)) {
    if (is_catch(model_fit)) {
      out <- cbind(mu, ntrad_out)
      #rename columns
      for (i in 1:model_fit@par_dims$q) {
        trad_names <- paste("n_traditional_", i + 1, sep = "")
      }
      colnames(out) <- c("mu", "n_traditional_1", trad_names)
    } else {
      out <- cbind(mu, ntrad_out)
      colnames(out) <- c("mu", "n_traditional")
    }
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
  out <- ifelse(model_fit@par_dims$phi == 1,
                TRUE, FALSE)
  return(out)
}

# function to get n traditional samples
get_ntrad <- function(pars, model_fit, mu, probability) {

  # number of traditional survey effort
  ntrad_out <- vector(length = length(mu))

  for (i in seq_along(mu)) {
    # number of traditional survey replicates
    ntrad <- seq(from = 0, to = 50000, by = 1)

    # P(X = 0 | mu) in one traditional survey trial
    if (is_negbin(model_fit)) {
      phi <- stats::median(unlist(rstan::extract(model_fit, pars = "phi[1]")))
      pr <- stats::pnbinom(q = 0, mu = mu[i], size = phi)
    } else {
      pr <- stats::ppois(q = 0, lambda = mu[i])
    }

    # dnbinom: x = failed events; size = successful events
    # probability of catching >0 animals based on ntrad
    prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr)

    # find value
    value <- match(min(prob[prob >= probability]), prob)
    ntrad_out[i] <- ntrad[value]
  }
  return(ntrad_out)
}

# function to get n traditional samples -- catchability coefficient
get_ntrad_q <- function(pars, model_fit, mu, probability) {

  # get catch coefficients
  # create empty q list
  q_list <- list()

  # fill in q list
  for (i in 1:model_fit@par_dims$q){
    name <- paste("q", i, sep = "")
    q_list[[name]] <- stats::median(rstan::extract(model_fit,
                                                   pars = "q")$q[, i])
  }
  q_list <- c(1, unlist(q_list))

  # number of traditional survey effort
  ntrad_out <- matrix(NA, nrow = length(mu), ncol = length(q_list))

  for (j in seq_along(q_list)) {
    for (i in seq_along(mu)) {
      # number of traditional survey replicates
      ntrad <- seq(from = 0, to = 50000, by = 1)

      # P(X = 0 | mu) in one traditional survey trial
      if (is_negbin(model_fit)) {
        phi <- stats::median(unlist(rstan::extract(model_fit, pars = "phi[1]")))
        pr <- stats::pnbinom(q = 0, mu = q_list[j] * mu[i], size = phi)
      } else {
        pr <- stats::ppois(q = 0, lambda = q_list[j] * mu[i])
      }

      # dnbinom: x = failed events; size = successful events
      # probability of catching >0 animals based on ntrad
      prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr)

      # find value
      value <- match(min(prob[prob >= probability]), prob)
      ntrad_out[i, j] <- ntrad[value]
    }
  }
  return(ntrad_out)
}

# function to get n eDNA samples
get_ndna <- function(pars, model_fit, mu, pcr_n, probability, cov_val) {

  # get beta
  alpha <- apply(rstan::extract(model_fit, pars = "alpha")$alpha, 2, "median")
  beta <- alpha %*% c(1, cov_val)

  # create output
  ndna_out <- vector(length = length(mu))

  for (i in seq_along(mu)) {
    # number of eDNA samples
    ndna <- seq(from = 0, to = 50000, by = 1)

    p11 <- mu[i] / (mu[i] + exp(beta))
    # P(at least one detection|mu) out of total bottles/PCR replicates
    prob <- 1 - stats::dbinom(x = 0, size = ndna * pcr_n, prob = p11)

    # find value
    value <- match(min(prob[prob >= probability]), prob)
    ndna_out[i] <- ndna[value]
  }

  return(ndna_out)
}

#function to check for divergent transitions
div_check <- function(x) {
  divergent <- sum(x[, "divergent__"])
  return(divergent)
}

# function for input checks
#' @srrstats {G5.2a} Pre-processing routines to check inputs have unique
#'   messages
# Helper function to validate conditions
validate_condition <- function(condition, message) {
  if (condition) stop(message)
}

det_calc_input_checks_1 <- function(model_fit, mu, cov_val,
                                    probability, pcr_n) {
  # 1. Validate model_fit is of class 'stanfit'
  validate_condition(!is(model_fit, "stanfit"),
                     "model_fit must be of class 'stanfit'.")

  # 2. Validate mu is a numeric vector of positive values
  validate_condition(!is.numeric(mu) || any(!mu > 0),
                     "mu must be a numeric vector of positive values.")

  # 3. Validate probability is a numeric value between 0 and 1
  validate_condition(!is.numeric(probability) || length(probability) != 1 ||
                       probability <= 0 || probability >= 1,
                     "probability must be a numeric value between 0 and 1.")
}

det_calc_input_checks_2 <- function(model_fit, mu, cov_val,
                                    probability, pcr_n) {

  # 4. Validate cov_val is numeric if provided
  if (!is.null(cov_val)) {
    validate_condition(!is.numeric(cov_val),
                       "cov_val must be a numeric vector.")

    # 5. Ensure cov_val matches model requirements
    if (model_fit@par_dims$alpha == 1) {
      stop(paste0("cov_val must be NULL if the model does not contain ",
                  "site-level covariates."))
    }
    validate_condition(length(cov_val) != (model_fit@par_dims$alpha - 1),
                       paste0("cov_val must be of the same length as the ",
                              "number of non-intercept site-level ",
                              "coefficients in the model."))
  }

  # 6. Validate cov_val is provided if covariates are in the model
  if (model_fit@par_dims$alpha > 1 && "p10" %in% model_fit@model_pars &&
        is.null(cov_val)) {
    stop(paste0("cov_val must be provided if the model contains site-level ",
                "covariates."))
  }

  # 7. Validate pcr_n is an integer
  validate_condition(pcr_n %% 1 != 0, "pcr_n should be an integer.")
}
