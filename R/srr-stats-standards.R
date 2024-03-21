#' srr_stats
#'
#' All of the following standards initially have `@srrstatsTODO` tags.
#' These may be moved at any time to any other locations in your code.
#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,
#' or `@srrstatsNA`, ensuring that references to every one of the following
#' standards remain somewhere within your code.
#' (These comments may be deleted at any time.)
#'
#' @srrstatsVerbose TRUE
#' @noRd
NULL

#' NA_standards
#' @srrstatsNA {G1.2} I attempted to add lifecycle statements with the
#' lifecycle package but was getting ROxygen errors when doing
#' `usethis::use_lifecycle()`. It would be nice to get some guidance on the
#' best way to do this.
#' @srrstatsNA {G1.4a} There are no non-exported functions.
#' @srrstatsNA {G1.6} This software makes no performance claims.
#' @srrstatsNA {G2.4b} It didn't seems obvious where to put as.numeric, as many
#' of my function input checks tested if values that are expected to be numeric
#' are in fact numeric.
#' @srrstatsNA {G2.4d,G2.4e} I don't believe I have anything that should be a
#' factor in my functions.
#' @srrstatsNA {G2.5} There are no inputs of type 'factor' in this software.
#' @srrstatsNA {G2.6} I do not believe this applies but would be interested in
#' guidance as to where I can implement this.
#' @srrstatsNA {G2.9} I do not belive I am doing any type conversions that
#' could result in lost data.
#' @srrstatsNA {G2.11,G2.12} I do not believe this applies but would be
#' interested in guidance as to where I can implement this.
#' @srrstatsNA {G2.14a,G2.14b,G2.14c} NAs in data input are meaningful and
#' should not create an error, nor be imputed
#' @srrstatsNA {G3.1} No covariance methods are used in this software.
#' @srrstatsNA {G3.1a} No covariance methods are used in this software.
#' @srrstatsNA {G4.0} There are no functions that write to local files
#' @srrstatsNA {G5.0} This software is applicable to specific data (i.e., data
#' from environmental DNA and traditional surveys), so NIST datasets would not
#' be applicable
#' @srrstatsNA {G5.4a} This software is not proposing a new algorithm and
#' instead uses functions from the `rstan` package, so it is unclear to me if
#' and how comparing to a C++ would be useful.
#' @srrstatsNA {G5.4b} This is not a new implementation of an existing method
#' @srrstatsNA {G5.4c} This is not a new implementation of an existing method
#' @srrstatsNA {G5.9, G5.9a} It is unclear to me how to add trivial noise to
#' the data, since the expected data are all integers.
#' @srrstatsNA {G5.10,G5.11,G5.11a,G5.12} I have not implemented any extended
#' tests
#' @srrstatsNA {BS1.3b} This software uses Stan that implements only one type
#' of algorithm (HMC/NUTS)
#' @srrstatsNA {BS1.5} This software enables convergence checks through
#' interaction with other packages, but does not enable multiple types of
#' convergence checks itself necessary for comparison
#' @srrstatsNA {BS4.1} This software is not presenting a new sampler
#' @srrstatsNA {BS1.3a, BS2.8} Since all my sampling functions uses `sampling`
#' from `rstan`, it is unclear how the user could use results from previous
#' runs to be used as starting points for a new run, although there is a way
#' for the user to specify initial values (with an example in the vignette)
#' @srrstatsNA {BS2.9} This software uses `rstan`'s `sampling`, which ensures
#' different seeds are used (from rstan's sampling() documentation: "Even if
#' multiple chains are used, only one seed is needed, with other chains having
#' seeds derived from that of the first chain to avoid dependent samples")
#' @srrstatsNA {BS2.10} This software ensures that the same seed is not passed
#' to multiple computational chains
#' @srrstatsNA {BS2.13, BS2.14} This software uses `rstan`'s `sampling`, and
#' while it seems possible to suppress all messages/warnings/progress bars, it
#' is unclear how to suppress just progress bars or just warnings.
#' @srrstatsNA {BS3.2} I do not provide a distinct routine for processing
#' collinear data, although I give the user a warning if input data has
#' perfect collinearity
#' @srrstatsNA {BS4.4} This software uses `rstan`'s `sampling`, and it does
#' not appear to yet be a mechanism of stopping the chain upon convergence
#' @srrstatsNA {BS4.6,BS4.7,BS5.4} This software does not include any
#' convergence checkers, although I do provide examples of how to examine
#' convergence in the package vignette.
#' @srrstatsNA {BS5.1} If I'm interpreting this standard correctly, I think
#' the `rstan` `stanfit` object returns appropriate metadata? i.e., if you are
#' estimating parameters at 20 sites, the output mu vector would be of length 20
#' @srrstatsNA {BS5.2} Parameters for non-uniform priors are provided by the
#' user, so it is unclear to me if this would be relevant.
#' @srrstatsNA {BS6.0} The main function, `jointModel()`, returns an object of
#' class `stanfit`. I would prefer to not print this output, as this software
#' is focused on a user who may want the actual model output to be abstracted.
#' I have therefore created `jointSummarize()` so that the user can interpret
#' the output
#' in terms of the parameters described in the 'About the model' section of
#' the vignette. If a user prefers, they can inspect the original `stanfit`
#' object.
#' @srrstatsNA {BS6.1} It is unclear to me which information/return objects
#' should be plotted
#' @srrstatsNA {BS6.5} Providing joint plotting options seems redundant to
#' functions in packages like `bayesplot`, but I'm open to suggestions.
#' @srrstatsNA {BS7.0,BS7.1,BS7.2} This software uses `rstan`'s HMC/NUTS
#' algorithm, where I am assuming that the algorithm will recover expected
#' prior and posterior distributions. Although I am open to suggestions for
#' how to test this here.
#' @srrstatsNA {BS7.3} It is unclear to me how to evaluate algorithmic
#' efficiency, as the MCMC can't be stopped upon convergence with `rstan`'s
#' sampling, and using something like Sys.time() doesn't seem like it would
#' provide an accurate estimate of computation time, but I'm open to other
#' suggestions.
#'
#'
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#' @noRd
NULL
