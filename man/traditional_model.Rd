% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/traditional_model.R
\name{traditional_model}
\alias{traditional_model}
\title{Specify and fit model using count data from traditional, non eDNA surveys}
\usage{
traditional_model(
  data,
  family = "poisson",
  q = FALSE,
  phi_priors = NULL,
  multicore = FALSE,
  initial_values = NULL,
  n_chain = 4,
  n_warmup = 500,
  n_iter = 3000,
  thin = 1,
  adapt_delta = 0.9,
  verbose = TRUE,
  seed = NULL
)
}
\arguments{
\item{data}{A list containing data necessary for model fitting. Valid tags
are \code{count} and \code{count_type}. \code{count} is a matrix or data frame of
traditional survey count data, with first dimension equal to the number of
sites (i) and second dimension equal to the maximum number of traditional
survey replicates at a given site (j). \code{count_type} is an optional matrix
or data frame of integers indicating gear type (k) used in corresponding
count data, with first dimension equal to the number of sites (i) and
second dimension equal to the maximum number of traditional survey
replicates at a given site (j). Values should be integers beginning with
1 (referring to the first gear type) to n (last gear type). Empty cells
should be NA and will be removed during processing. Sites, i, should be
consistent in all count data.}

\item{family}{The distribution class used to model traditional survey count
data. Options include poisson ('poisson'), negative binomial ('negbin'),
and gamma ('gamma'). Default value is 'poisson'.}

\item{q}{A logical value indicating whether to estimate gear scaling
coefficients, q, for traditional survey gear types (TRUE) or to not
estimate gear scaling coefficients, q, for traditional survey gear types
(FALSE). Default value is FALSE.}

\item{phi_priors}{A numeric vector indicating gamma distribution
hyperparameters (shape, rate) used as the prior distribution for phi, the
overdispersion in the negative binomial distribution for traditional survey
gear data. Used when family = 'negbin.' If family = 'negbin', then
default vector is c(0.25,0.25), otherwise, default is NULL.}

\item{multicore}{A logical value indicating whether to parallelize chains
with multiple cores. Default is FALSE.}

\item{initial_values}{A list of lists of initial values to use in MCMC. The
length should equal the number of MCMC chains. Initial values can be
provided for parameters: mu and q. If no initial values are provided,
default random values are drawn.}

\item{n_chain}{Number of MCMC chains. Default value is 4.}

\item{n_warmup}{A positive integer specifying the number of warm-up MCMC
iterations. Default value is 500.}

\item{n_iter}{A positive integer specifying the number of iterations for each
chain (including warmup). Default value is 3000.}

\item{thin}{A positive integer specifying the period for saving samples.
Default value is 1.}

\item{adapt_delta}{Numeric value between 0 and 1 indicating target average
acceptance probability used in \code{rstan::sampling}. Default value is 0.9.}

\item{verbose}{Logical value controlling the verbosity of output (i.e.,
warnings, messages, progress bar). Default is TRUE.}

\item{seed}{A positive integer seed used for random number generation in
MCMC. Default is NULL, which means the seed is generated from 1 to the
maximum integer supported by R.}
}
\value{
A list of:
\itemize{
\item a model object of class \code{stanfit} returned by \code{rstan::sampling}
\item initial values used in MCMC
}
}
\description{
This function implements a Bayesian model that estimates expected species
catch rate using count data from traditional, non eDNA surveys. When
multiple traditional gear types are used, an optional variation allows
estimation of gear scaling coefficients, which scale the catchability of
gear types relative to the expected catch rate of a reference gear type.
Model is implemented using Bayesian inference using the \code{rstan} package,
which uses Hamiltonian Monte Carlo to simulate the posterior distributions.
See more examples in the
\href{https://ednajoint.netlify.app}{Package
Vignette}.
}
\note{
Before fitting the model, this function checks to ensure that the
model specification is possible given the data files. These checks include:
\itemize{
\item  All tags in data are valid (i.e., include count and count_type).
\item  Number of sites in count and count type data are equal.
\item  All data are numeric (i.e., integer or NA).
\item  Empty data cells (NA) match in count and count_type.
\item  family is 'poisson', 'negbin', or 'gamma'.
\item  phi_priors (if used) is a vector of two numeric values.
}

If any of these checks fail, the function returns an error message.
}
\examples{
\donttest{
# Load data
data(green_crab_data)

# Examine data in list
# This function uses only traditional survey count data and optionally
# the count type data
names(green_crab_data)

# Note that the surveyed sites (rows) should match in the data
dim(green_crab_data$count)[1]
dim(green_crab_data$count_type)[1]

# Fit a model without estimating a gear scaling coefficient for traditional
# survey gear types.
# This model assumes all traditional survey methods have the same
# catchability.
# Count data is modeled using a poisson distribution.
fit_no_q <- traditional_model(data = green_crab_data, family = "poisson",
                              q = FALSE, phi_priors = NULL,
                              multicore = FALSE, verbose = TRUE)


# Fit a model estimating a gear scaling coefficient for traditional survey
# gear types.
# This model does not assume all traditional survey methods have the same
# catchability.
# Count data is modeled using a negative binomial distribution.
fit_q <- traditional_model(data = green_crab_data, family = "negbin",
                           q = TRUE, phi_priors = c(0.25,0.25),
                           multicore = FALSE, initial_values = NULL,
                           n_chain = 4, n_warmup = 500, n_iter = 3000,
                           thin = 1, adapt_delta = 0.9, verbose = TRUE,
                           seed = 123)
}

}
