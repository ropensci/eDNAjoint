#' Specify and fit joint model using count data from traditional surveys and
#' eDNA PCR data
#'
#' This function implements a Bayesian model that integrates data from paired
#' eDNA and traditional surveys, as described in Keller et al. (2022)
#' <doi.org/10.1002/eap.2561>. The model estimates parameters including
#' the expected species catch rate and the probability of false positive eDNA
#' detection. This function allows for optional model variations, like inclusion
#' of site-level covariates that scale the sensitivity of eDNA sampling relative
#' to traditional sampling, as well as estimation of gear scaling coefficients
#' that scales the relative catchability of multiple traditional gear types.
#' Model is implemented using Bayesian inference using the `rstan` package,
#' which uses Hamiltonian Monte Carlo to simulate the posterior distributions.
#' See more examples in the
#' \href{https://ednajoint.netlify.app/}{Package
#' Vignette}.
#'
#' @srrstats {G1.0,G1.1} This software makes available a algorithm/model
#'   that was previously published in the literature. The literature reference
#'   for the joint model is provided here.
#'
#' @srrstats {G1.4,G1.3} Roxygen function documentation begins here, with
#'   definitions of statistical terminology and inputs
#' @export
#' @srrstats {BS1.1,BS3.0,G2.14} Descriptions of how to enter data, description
#'   of how NAs are handled (which are informative and should be deliberately
#'   included in input data, where necessary)
#' @param data A list containing data necessary for model fitting. Valid tags
#'   are `pcr_n`, `pcr_k`, `count`, `count_type`, and `site_cov`. `pcr_n` and
#'   `pcr_k` are matrices or data frames with first dimension equal to the
#'   number of sites (i) and second dimension equal to the maximum number of
#'   eDNA samples at a given site (m). `pcr_n` contains the total number of
#'   PCR replicates per site and eDNA sample, and `pcr_k` contains the total
#'   number of positive PCR detections per site and eDNA sample. `count` is a
#'   matrix or data frame of traditional survey count data, with first
#'   dimension equal to the number of sites (i) and second dimension equal to
#'   the maximum number of traditional survey replicates at a given site (j).
#'   `count_type` is an optional matrix or data frame of integers indicating
#'   gear type used in corresponding count data, with first dimension equal to
#'   the number of sites (i) and second dimension equal to the maximum number
#'   of traditional survey replicates at a given site. Values should be
#'   integers beginning with 1 (referring to the first gear type) to n (last
#'   gear type). `site_cov` is an optional matrix or data frame of site-level
#'   covariate data, with first dimension equal to the number of sites (i).
#'   `site_cov` should include informative column names. Empty cells should
#'   be NA and will be removed during processing. Sites, i, should be consistent
#'   in all PCR, count, and site covariate data.
#' @srrstats {G2.1a} Here are explicit documentation of vector input types
#' @param cov A character vector indicating the site-level covariates to include
#'   in the model. Default value is NULL.
#' @param family The distribution class used to model traditional survey count
#'   data. Options include poisson ('poisson'), negative binomial ('negbin'),
#'   and gamma ('gamma'). Default value is 'poisson'.
#' @srrstats {BS1.0,G2.1a,BS1.2} Description of hyperparameters and how to
#'   specify prior distributions, explicit documentation of vector input types
#' @param p10_priors A numeric vector indicating beta distribution
#'   hyperparameters (alpha, beta) used as the prior distribution for the eDNA
#'   false positive probability (p10). Default vector is c(1,20).
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
#'   provided for parameters: beta, p10 (log-scale), mu, q, alpha. If no
#'   initial values are provided, default random values are drawn.
#' @srrstats {BS1.3} Description of parameters used in the computational process
#'   begins here
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
#' @srrstats {BS5.0} function returns initial values used in computation
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
#' \item  All tags in data are valid (i.e., include pcr_n, pcr_k, count,
#'   count_type, and site_cov).
#' \item  Dimensions of pcr_n and pcr_k are equal, and dimensions of count and
#'   count_type are equal (if count_type provided).
#' \item  Number of sites in PCR and count data are equal.
#' \item  All data are numeric (i.e., integer or NA).
#' \item  Empty data cells (NA) match in pcr_n and pcr_k and in count and
#'   count_type.
#' \item  family is either 'poisson', 'negbin', or 'gamma'.
#' \item  p10_priors and phi_priors (if used) is a vector of two numeric values.
#' \item  site_cov has same number of rows as pcr_n and count, if present
#' \item  site_cov is numeric, if present
#' \item  cov values match column names in site_cov, if present
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' @examples
#' \donttest{
#' # Ex. 1: Implementing the joint model
#'
#' # Load data
#' data(goby_data)
#'
#' # Examine data in list
#' names(goby_data)
#'
#' # Note that the surveyed sites (rows) should match in all data
#' dim(goby_data$pcr_n)[1]
#' dim(goby_data$count)[1]
#'
#' # Fit a basic model with paired eDNA and traditional survey data.
#' # Count data is modeled using a poisson distribution.
#' fit <- joint_model(data = goby_data, family = "poisson",
#'                    p10_priors = c(1, 20),
#'                    multicore = FALSE)
#'
#' # Ex. 2: Implementing the joint model with site-level covariates
#'
#' # With the same data, fit a model including 'Filter_time' and 'Salinity'
#' # site-level covariates
#' # These covariates will scale the sensitivity of eDNA sampling relative to
#' # traditional surveys
#' # Count data is modeled using a poisson distribution.
#' fit_cov <- joint_model(data = goby_data, cov = c('Filter_time','Salinity'),
#'                        family = "poisson", p10_priors = c(1, 20),
#'                        multicore = FALSE)
#'
#'
#' # Ex. 3: Implementing the joint model with multiple traditional gear types
#'
#' # Load data
#' data(green_crab_data)
#'
#' # Examine data in list
#' names(green_crab_data)
#'
#' # Note that the surveyed sites (rows) should match in all data
#' dim(green_crab_data$pcr_n)[1]
#' dim(green_crab_data$count)[1]
#'
#' # Fit a model estimating a gear scaling coefficient for traditional survey
#' # gear types.
#' # This model does not assume all traditional survey methods have the same
#' # catchability.
#' # Count data is modeled using a negative binomial distribution.
#' fit_q <- joint_model(data = green_crab_data, cov = NULL, family = "negbin",
#'                      p10_priors = c(1, 20), q = TRUE,
#'                      phi_priors = c(0.25, 0.25),
#'                      multicore = FALSE, initial_values = NULL,
#'                      n_chain = 4, n_warmup = 500,
#'                      n_iter = 3000, thin = 1, adapt_delta = 0.9,
#'                      verbose = TRUE, seed = 123)
#' }
#'

joint_model <- function(data, cov = NULL, family = "poisson",
                        p10_priors = c(1, 20),
                        q = FALSE, phi_priors = NULL, multicore = FALSE,
                        initial_values = NULL, n_chain = 4, n_warmup = 500,
                        n_iter = 3000, thin = 1, adapt_delta = 0.9,
                        verbose = TRUE, seed = NULL) {


  # make character inputs case-insensitive
  #' @srrstats {G2.3b} Allow case-insensitive character parameter values
  family <- tolower(family)


  # get phi_priors
  if (family != "negbin") {
    phi_priors <- NULL
  } else if (family == "negbin") {
    if (is.null(phi_priors)) {
      phi_priors <- c(0.25, 0.25)
    }
  }

  ####
  # input checks
  #' @srrstats {G2.1} Types of inputs are checked/asserted using the following
  #'   helper functions

  if (q == TRUE) {
    # model with catchability coefficients
    catchability_checks_1(data, cov)
    catchability_checks_2(data, cov)
  } else {
    # model without catchability coefficients
    no_catchability_checks(data, cov)
  }

  # model with covariates
  if (all(!is.null(cov))) {
    covariate_checks(data, cov)
  }

  # all models
  all_checks_1(data, cov, family, p10_priors, phi_priors, n_chain, n_warmup,
               n_iter, thin, adapt_delta, seed)
  all_checks_2(data, cov, family, p10_priors, phi_priors, n_chain, n_warmup,
               n_iter, thin, adapt_delta, seed)

  # initial value checks
  initial_values_checks_1(initial_values, data, cov, n_chain)
  initial_values_checks_2(initial_values, data, cov, n_chain)

  ###
  #convert data to long format

  # convert PCR data to long format
  #' @srrstats {G2.7} Use as.data.frame() to allow input list of any tabular
  #'   form (i.e., matrix, etc.)
  pcr_all <- as.data.frame(data$pcr_n) |>
    dplyr::mutate(L_ind = seq_len(nrow(data$pcr_n))) |>
    tidyr::pivot_longer(cols = ! L_ind, values_to = "n_N") |>
    #' @srrstats {G2.15} Software does not assume non-missingness and actually
    #'   expects it if the number of observations across sites is unequal
    tidyr::drop_na()
  pcr_k_df <- as.data.frame(data$pcr_k) |>
    dplyr::mutate(L_ind = seq_len(nrow(data$pcr_k))) |>
    tidyr::pivot_longer(cols = ! L_ind, values_to = "n_K") |>
    #' @srrstats {G2.15} Software does not assume non-missingness and actually
    #'   expects it if the number of observations across sites is unequal
    tidyr::drop_na()
  pcr_all$n_K <- pcr_k_df$n_K

  # convert count data to long format
  #' @srrstats {G2.7} Use as.data.frame() to allow input list of any tabular
  #'   form (i.e., matrix, etc.)
  count_all <- as.data.frame(data$count) |>
    dplyr::mutate(L_ind = seq_len(nrow(data$count))) |>
    tidyr::pivot_longer(cols = ! L_ind, values_to = "count") |>
    #' @srrstats {G2.15} Software does not assume non-missingness and actually
    #' expects it if the number of observations across sites is unequal
    tidyr::drop_na()

  # subset count data to remove sites without traditional samples
  count_df <- as.data.frame(data$count)
  sub_count <- as.data.frame(
    count_df[rowSums(is.na(count_df)) != ncol(count_df), ]
  )

  # add site index to count data
  index_match <- as.data.frame(cbind(unique(count_all$L_ind),
                                     seq_len(nrow(sub_count))))
  colnames(index_match) <- c("L_ind", "R_ind")
  count_all <- dplyr::left_join(count_all, index_match, by = "L_ind")

  # get site indices with paired and unpaired dna samples
  trad_ind <- index_match$L_ind
  dna_ind <- unique(pcr_all$L_ind)[!unique(pcr_all$L_ind) %in% trad_ind]

  # subset PCR data -- paired
  pcr_all_trad <- pcr_all[pcr_all$L_ind %in% trad_ind, ]
  l_match_trad <- as.data.frame(cbind(unique(pcr_all_trad$L_ind),
                                      seq_along(unique(pcr_all_trad$L_ind))))
  colnames(l_match_trad) <- c("L_ind", "L_unique")
  pcr_all_trad <- dplyr::left_join(pcr_all_trad, l_match_trad, by = "L_ind")

  # subset PCR data -- unpaired
  if (length(dna_ind) > 0) {
    pcr_all_dna <- pcr_all[pcr_all$L_ind %in% dna_ind, ]
    l_match_dna <- as.data.frame(cbind(unique(pcr_all_dna$L_ind),
                                       seq_along(unique(pcr_all_dna$L_ind))))
    colnames(l_match_dna) <- c("L_ind", "L_unique")
    pcr_all_dna <- dplyr::left_join(pcr_all_dna, l_match_dna, by = "L_ind")
  } else {
    pcr_all_dna <- as.data.frame(matrix(NA, nrow = 0, ncol = 4))
    colnames(pcr_all_dna) <- c("L_ind", "n_N", "n_K", "L_unique")
    l_match_dna <- as.data.frame(matrix(NA, nrow = 0, ncol = 4))
  }

  # if present, prepare covariate data
  if (all(!is.null(cov))) {
    #' @srrstats {G2.7,G2.10} Use as.data.frame() to allow input list of any
    #'   tabular form (i.e., matrix, etc.) and converts before filtering columns
    #'   based on input 'cov'
    site_mat <- as.data.frame(as.data.frame(data$site_cov)[, cov])
    site_mat <- cbind(as.data.frame(rep(1, length(site_mat[, 1]))), site_mat)
    colnames(site_mat) <- c("int", cov)
  } else {
    # data structures when there are no covariates
    site_mat <- rep(1,
                    length(unique(pcr_all_trad$L_ind)) +
                      length(unique(pcr_all_dna$L_ind)))
  }

  # convert p10 prior
  # p10 prior: convert beta(1,20) to lognormal distribution
  # moment match from beta(alpha,beta) to normal(mu, sigma^2)
  alpha <- p10_priors[1]
  beta <- p10_priors[2]
  mu <- alpha / (alpha + beta)
  sigma_2 <- (alpha * beta) / ((alpha + beta) ^ 2 * (alpha + beta + 1))
  # convert normal(mu, sigma^2) to lognormal(mu, sigma^2)
  mu_ln <- log(mu ^ 2 / sqrt(mu ^ 2 + sigma_2))
  sigma_2_ln <- log(1 + sigma_2 / mu ^ 2)
  sigma_ln <- sqrt(sigma_2_ln)

  # create data that will be present in all model variations
  model_data <- list(
    n_S = nrow(pcr_all_trad),
    S_dna = nrow(pcr_all_dna),
    n_C = nrow(count_all),
    L_ind = pcr_all_trad$L_unique,
    L_dna = pcr_all_dna$L_unique,
    R_ind = count_all$R_ind,
    Nloc_dna = length(unique(pcr_all_dna$L_ind)),
    Nloc_trad = length(unique(pcr_all_trad$L_ind)),
    trad_ind = trad_ind,
    dna_ind = as.array(dna_ind),
    n_E = count_all$count,
    n_N = pcr_all_trad$n_N,
    n_K = pcr_all_trad$n_K,
    N_dna = pcr_all_dna$n_N,
    K_dna = pcr_all_dna$n_K,
    nsitecov = length(cov) + 1,
    mat_site = as.matrix(site_mat),
    p10_priors = c(mu_ln, sigma_ln),
    alphapriors = c(0, 10),
    control = list(adapt_delta = adapt_delta)
  )

  # append data based on catchability
  if (is_catch_type(q)) {
    # add count type data to count df
    q_ref <- 1
    #' @srrstats {G2.7} Use as.data.frame() to allow input list of any tabular
    #'   form (i.e., matrix, etc.)
    count_type_df <- as.data.frame(data$count_type) |>
      dplyr::mutate(L_ind = seq_len(nrow(data$count_type))) |>
      tidyr::pivot_longer(cols = ! L_ind, values_to = "count_type") |>
      tidyr::drop_na()
    count_all$count_type <- count_type_df$count_type

    #create vector of q coefficient names
    counttypes <- unique(count_all$count_type)
    names <- counttypes[!counttypes == q_ref]
    #' @srrstats {G2.4,G2.4c} Explicit conversion to character
    q_names <- as.character(paste0("q_", names))

    # append data
    model_data <- rlist::list.append(
      model_data,
      nparams = length(q_names),
      mat = as.integer(count_all$count_type),
      ctch = 1
    )
  } else {
    q_names <- NULL

    model_data <- rlist::list.append(
      model_data,
      nparams = 0,
      mat = as.integer(rep(1, nrow(count_all))),
      ctch = 0
    )
  }

  # append data if family == negbin
  if (get_family_index(family) == 2) {
    model_data <- rlist::list.append(
      model_data,
      phi_priors = phi_priors,
      negbin = 1
    )
  } else if (get_family_index(family) == 1) {
    model_data <- rlist::list.append(
      model_data,
      phi_priors = c(1, 1),
      negbin = 0
    )
  } else if (get_family_index(family) == 3) {
    model_data <- rlist::list.append(
      model_data,
      bgammapriors = c(0.25, 0.25),
      agammapriors = c(0.01, 0.01)
    )
  }

  # set up core number
  cores <- ifelse(multicore == TRUE,
                  parallel::detectCores(),
                  1)

  # get seed
  seed_mod <- ifelse(!is.null(seed),
                     as.integer(seed),
                     sample.int(.Machine$integer.max, 1))

  # get stan model
  model_index <- get_stan_model(family)

  # get initial values
  inits <- get_inits(n_chain, pcr_all, initial_values, cov, l_match_trad,
                     l_match_dna, data, q_names)

  # run model
  out <- rstan::sampling(
    c(stanmodels$joint_count,
      stanmodels$joint_continuous)[model_index][[1]],
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

  # check for divergent transitions
  div_trans <- sum(lapply(rstan::get_sampler_params(out,
                                                    inc_warmup = FALSE),
                          div_check)[[1]])
  # print either troubleshooting or visualization tips
  if (div_trans > 0 && verbose) {
    url <- "https://ednajoint.netlify.app/tips#troubleshooting-tips"
    message <- "Refer to the eDNAjoint guide for troubleshooting tips: "
  } else {
    url <- "https://ednajoint.netlify.app/tips#visualization-tips"
    message <- "Refer to the eDNAjoint guide for visualization tips: "
  }
  cat(message, url, "\n")

  # add chain names to init list
  names(inits) <- paste0("chain", seq(1, n_chain, 1))

  # Create a list to store the results
  #' @srrstats {BS5.0} function returns initial values used in computation
  result_list <- list(model = out, inits = inits)

  #' @srrstats {BS5.5} the `model` return object is of class `stanfit`, which
  #'   includes information about convergence
  return(result_list)
}
