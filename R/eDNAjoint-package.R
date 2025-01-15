#' The 'eDNAjoint' package.
#'
#' @description Models integrate environmental DNA (eDNA) detection data and
#'   traditional survey data to jointly estimate species catch rate (see
#'   \href{https://bookdown.org/abigailkeller/eDNAjoint_vignette/}{Package
#'   Vignette}). Models can be used with count data via traditional survey
#'   methods (i.e., trapping, electrofishing, visual) and replicated eDNA
#'   detection/nondetection data via polymerase chain reaction (i.e., PCR or
#'   qPCR) from multiple survey locations. Estimated parameters include
#'   probability of a false positive eDNA detection, a site-level covariates
#'   that scale the sensitivity of eDNA surveys relative to traditional surveys,
#'   and catchability coefficients for traditional gear types. Models are
#'   implemented with a Bayesian framework (Markov chain Monte Carlo) using the
#'   'Stan' probabilistic programming language.
#'
#' @docType package
#' @name eDNAjoint-package
#' @aliases eDNAjoint
#' @useDynLib eDNAjoint, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan.
#'   https://mc-stan.org
#'
## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end
NULL
