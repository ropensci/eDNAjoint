#' Specify and fit model using count data from traditional, non eDNA surveys
#'
#' This function implements a Bayesian model that estimates expected species
#' catch rate using count data from traditional, non eDNA surveys. When
#' multiple traditional gear types are used, an optional variation allows
#' estimation of gear scaling coefficients, which scale the catchability of
#' gear types relative to the expected catch rate of a reference gear type.
#' Model is implemented using Bayesian inference using the `rstan` package,
#' which uses Hamiltonian Monte Carlo to simulate the posterior distributions.
#' See more examples in the
#' \href{https://ednajoint.netlify.app}{Package
#' Vignette}.
#'
#' @srrstats {G1.4,G1.3} Roxygen function documentation begins here, with
#'   definitions of statistical terminology and inputs
#' @export
#' @srrstats {BS1.1,BS3.0,G2.14} Descriptions of how to enter data, description
#'   of how NAs are handled (which are informative and should be deliberately
#'   included in input data, where necessary)
#' @param data A list containing data necessary for model fitting. Valid tags
#'   are `count` and `count_type`. `count` is a matrix or data frame of
#'   traditional survey count data, with first dimension equal to the number of
#'   sites (i) and second dimension equal to the maximum number of traditional
#'   survey replicates at a given site (j). `count_type` is an optional matrix
#'   or data frame of integers indicating gear type (k) used in corresponding
#'   count data, with first dimension equal to the number of sites (i) and
#'   second dimension equal to the maximum number of traditional survey
#'   replicates at a given site (j). Values should be integers beginning with
#'   1 (referring to the first gear type) to n (last gear type). Empty cells
#'   should be NA and will be removed during processing. Sites, i, should be
#'   consistent in all count data.
#' @param family The distribution class used to model traditional survey count
#'   data. Options include poisson ('poisson'), negative binomial ('negbin'),
#'   and gamma ('gamma'). Default value is 'poisson'.
#' @param q A logical value indicating whether to estimate gear scaling
#'   coefficients, q, for traditional survey gear types (TRUE) or to not
#'   estimate gear scaling coefficients, q, for traditional survey gear types
#'   (FALSE). Default value is FALSE.
#' @srrstats {BS1.0,G2.1a,BS1.2} Description of hyperparameters and how to
#'   specify prior distributions, explicit documentation of vector input types
#' @param phi_priors A numeric vector indicating gamma distribution
#'   hyperparameters (shape, rate) used as the prior distribution for phi, the
#'   overdispersion in the negative binomial distribution for traditional survey
#'   gear data. Used when family = 'negbin.' If family = 'negbin', then
#'   default vector is c(0.25,0.25), otherwise, default is NULL.
#' @param multicore A logical value indicating whether to parallelize chains
#'   with multiple cores. Default is FALSE.
#' @srrstats {BS2.7,BS2.11} Option for user to provide initial values for each
#'   chain
#' @param initial_values A list of lists of initial values to use in MCMC. The
#'   length should equal the number of MCMC chains. Initial values can be
#'   provided for parameters: mu and q. If no initial values are provided,
#'   default random values are drawn.
#' @srrstats {BS1.3} Description of parameters used in the computational
#'   process begins here
#' @param n_chain Number of MCMC chains. Default value is 4.
#' @param n_warmup A positive integer specifying the number of warm-up MCMC
#'   iterations. Default value is 500.
#' @param n_iter A positive integer specifying the number of iterations for each
#'   chain (including warmup). Default value is 3000.
#' @param thin A positive integer specifying the period for saving samples.
#'   Default value is 1.
#' @param adapt_delta Numeric value between 0 and 1 indicating target average
#'   acceptance probability used in `rstan::sampling`. Default value is 0.9.
#' @srrstats {BS2.12} Parameter controlling the verbosity of output
#' @param verbose Logical value controlling the verbosity of output (i.e.,
#'   warnings, messages, progress bar). Default is TRUE.
#' @param seed A positive integer seed used for random number generation in
#'   MCMC. Default is NULL, which means the seed is generated from 1 to the
#'   maximum integer supported by R.
#' @return A list of:
#' \itemize{
#' \item a model object of class `stanfit` returned by `rstan::sampling`
#' \item initial values used in MCMC
#' }
#'
#' @srrstats {G2.0a} Explicit secondary documentation of any expectations on
#'   lengths of inputs
#' @note  Before fitting the model, this function checks to ensure that the
#'   model specification is possible given the data files. These checks include:
#' \itemize{
#' \item  All tags in data are valid (i.e., include count and count_type).
#' \item  Number of sites in count and count type data are equal.
#' \item  All data are numeric (i.e., integer or NA).
#' \item  Empty data cells (NA) match in count and count_type.
#' \item  family is 'poisson', 'negbin', or 'gamma'.
#' \item  phi_priors (if used) is a vector of two numeric values.
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' @examples
#' \donttest{
#' # Load data
#' data(green_crab_data)
#'
#' # Examine data in list
#' # This function uses only traditional survey count data and optionally
#' # the count type data
#' names(green_crab_data)
#'
#' # Note that the surveyed sites (rows) should match in the data
#' dim(green_crab_data$count)[1]
#' dim(green_crab_data$count_type)[1]
#'
#' # Fit a model without estimating a gear scaling coefficient for traditional
#' # survey gear types.
#' # This model assumes all traditional survey methods have the same
#' # catchability.
#' # Count data is modeled using a poisson distribution.
#' fit_no_q <- traditional_model(data = green_crab_data, family = "poisson",
#'                               q = FALSE, phi_priors = NULL,
#'                               multicore = FALSE, verbose = TRUE)
#'
#'
#' # Fit a model estimating a gear scaling coefficient for traditional survey
#' # gear types.
#' # This model does not assume all traditional survey methods have the same
#' # catchability.
#' # Count data is modeled using a negative binomial distribution.
#' fit_q <- traditional_model(data = green_crab_data, family = "negbin",
#'                            q = TRUE, phi_priors = c(0.25,0.25),
#'                            multicore = FALSE, initial_values = NULL,
#'                            n_chain = 4, n_warmup = 500, n_iter = 3000,
#'                            thin = 1, adapt_delta = 0.9, verbose = TRUE,
#'                            seed = 123)
#' }
#'

traditional_model <- function(data, family = "poisson",
                              q = FALSE, phi_priors = NULL,
                              multicore = FALSE, initial_values = NULL,
                              n_chain = 4, n_warmup = 500,
                              n_iter = 3000, thin = 1,
                              adapt_delta = 0.9, verbose = TRUE, seed = NULL) {


  # make character inputs case-insensitive
  #' @srrstats {G2.3b} Allow case-insensitive character parameter values
  family <- tolower(family)

  # get phi_priors
  if (family != "negbin") {
    phi_priors <- NULL
  } else if (family == "negbin" && is.null(phi_priors)) {
    phi_priors <- c(0.25, 0.25)
  } else if (family == "negbin" && !is.null(phi_priors)) {
    phi_priors <- phi_priors
  }

  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using this helper
  #'   function
  trad_model_input_checks_1(data, family, q, phi_priors, n_chain,
                            n_warmup, n_iter,
                            thin, adapt_delta, seed)
  trad_model_input_checks_2(data, family, q, phi_priors, n_chain,
                            n_warmup, n_iter,
                            thin, adapt_delta, seed)

  # initial value checks
  if (all(!is.null(initial_values))) {
    initial_values_checks_trad(initial_values, data, n_chain)
  }

  ###
  #convert count data to long format
  #' @srrstats {G2.7} Use as.data.frame() to allow input list of any tabular
  #'   form (i.e., matrix, etc.)
  count_all <- as.data.frame(data$count) |>
    dplyr::mutate(L_ind = seq_len(nrow(data$count))) |>
    tidyr::pivot_longer(cols = ! L_ind, values_to = "count") |>
    #' @srrstats {G2.15} Software does not assume non-missingness and actually
    #'   expects it if the number of observations across sites is unequal
    tidyr::drop_na()

  #if q == TRUE, add count type data to count df
  if (q == TRUE) {
    q_ref <- 1
    #' @srrstats {G2.7} Use as.data.frame() to allow input list of any tabular
    #'   form (i.e., matrix, etc.)
    count_type_df <- as.data.frame(data$count_type) |>
      dplyr::mutate(L_ind = seq_len(data$count_type)[1]) |>
      tidyr::pivot_longer(cols = ! L_ind, values_to = "count_type") |>
      #' @srrstats {G2.15} Software does not assume non-missingness and
      #'   actually expects it if the number of observations across sites is
      #'   unequal
      tidyr::drop_na()
    count_all$count_type <- count_type_df$count_type

    #create vector of q coefficient names
    counttypes <- unique(count_all$count_type)
    names <- counttypes[!counttypes == q_ref]
    #' @srrstats {G2.4,G2.4c} Explicit conversion to character
    q_names <- as.character(paste0("q_", names))

  }

  # set up core number
  cores <- ifelse(multicore == TRUE, parallel::detectCores(), 1)

  # create data that will be present in all model variations
  model_data <- list(
    Nloc = length(unique(count_all$L_ind)),
    n_C = nrow(count_all),
    R_ind = count_all$L_ind,
    n_E = count_all$count,
    control = list(adapt_delta = adapt_delta)
  )

  # append data based on catchability
  if (q == TRUE) {
    model_data <- rlist::list.append(
      model_data,
      nparams = length(q_names),
      mat = as.integer(count_all$count_type),
      ctch = 1
    )
  } else {
    model_data <- rlist::list.append(
      model_data,
      nparams = 0,
      mat = as.integer(rep(1, nrow(count_all))),
      ctch = 0
    )
  }

  # append data if family == negbin
  if (family == "negbin") {
    model_data <- rlist::list.append(
      model_data,
      phi_priors = phi_priors,
      negbin = 1
    )
  } else if (family == "poisson") {
    model_data <- rlist::list.append(
      model_data,
      phi_priors = c(1, 1),
      negbin = 0
    )
  } else if (family == "gamma") {
    model_data <- rlist::list.append(
      model_data,
      betapriors = c(0.25, 0.25),
      alphapriors = c(0.01, 0.01)
    )
  }

  # get stan model
  model_index <- get_stan_model(family)

  # get initial values
  if (q == TRUE) {
    inits <- init_trad_catchability(n_chain, count_all, q_names,
                                    initial_values)
  } else {
    inits <- init_trad(n_chain, count_all, initial_values)
  }

  # get seed
  seed_mod <- ifelse(!is.null(seed),
                     as.integer(seed),
                     sample.int(.Machine$integer.max, 1))

  # run model
  out <- rstan::sampling(
    c(stanmodels$traditional_count,
      stanmodels$traditional_continuous)[model_index][[1]],
    data = model_data,
    cores = cores,
    seed = seed_mod,
    #' @srrstats {G2.4,G2.4a} explicit conversion to
    #'   integers for sampling arguments
    chains = as.integer(n_chain),
    thin = as.integer(thin),
    warmup = as.integer(n_warmup),
    iter = as.integer(n_iter),
    init = inits,
    refresh = ifelse(verbose == TRUE, 500, 0)
  )



  # assert that the log likelihood is a double
  #' @srrstats {G5.3} assert that model run worked and the log likelihood is
  #'   valid (i.e., not NA)
  stopifnot(is.double(sum(colMeans(rstan::extract(out,
                                                  par = "log_lik")$log_lik))))

  # Create a list to store the results
  #' @srrstats {BS5.0} function returns initial values used in computation
  result_list <- list(model = out, inits = inits)

  #' @srrstats {BS5.5} the `model` return object is of class `stanfit`, which
  #'   includes information about convergence
  return(result_list)
}
