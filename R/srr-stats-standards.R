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
#'
#' @srrstatsTODO {G1.1} *Statistical Software should document whether the algorithm(s) it implements are:* - *The first implementation of a novel algorithm*; or - *The first implementation within **R** of an algorithm which has previously been implemented in other languages or contexts*; or - *An improvement on other implementations of similar algorithms in **R***.
#' @srrstatsTODO {G1.2} *Statistical Software should include a* Life Cycle Statement *describing current and anticipated future states of development.*
#' @srrstatsTODO {G1.3} *All statistical terminology should be clarified and unambiguously defined.*
#' @srrstatsTODO {G1.6} *Software should include code necessary to compare performance claims with alternative implementations in other R packages.*
#' @srrstatsTODO {G2.4} *Provide appropriate mechanisms to convert between different data types, potentially including:*
#' @srrstatsTODO {G2.4a} *explicit conversion to `integer` via `as.integer()`*
#' @srrstatsTODO {G2.4b} *explicit conversion to continuous via `as.numeric()`*
#' @srrstatsTODO {G2.4c} *explicit conversion to character via `as.character()` (and not `paste` or `paste0`)*
#' @srrstatsTODO {G2.4d} *explicit conversion to factor via `as.factor()`*
#' @srrstatsTODO {G2.4e} *explicit conversion from factor via `as...()` functions*
#' @srrstatsTODO {G2.6} *Software which accepts one-dimensional input should ensure values are appropriately pre-processed regardless of class structures.*
#' @srrstatsTODO {G2.7} *Software should accept as input as many of the above standard tabular forms as possible, including extension to domain-specific forms.*
#' @srrstatsTODO {G2.9} *Software should issue diagnostic messages for type conversion in which information is lost (such as conversion of variables from factor to character; standardisation of variable names; or removal of meta-data such as those associated with [`sf`-format](https://r-spatial.github.io/sf/) data) or added (such as insertion of variable or column names where none were provided).*
#' @srrstatsTODO {G2.10} *Software should ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behaviour, and should ensure all column-extraction operations behave consistently regardless of the class of tabular data used as input.*
#' @srrstatsTODO {G2.11} *Software should ensure that `data.frame`-like tabular objects which have columns which do not themselves have standard class attributes (typically, `vector`) are appropriately processed, and do not error without reason. This behaviour should be tested. Again, columns created by the [`units` package](https://github.com/r-quantities/units/) provide a good test case.*
#' @srrstatsTODO {G2.12} *Software should ensure that `data.frame`-like tabular objects which have list columns should ensure that those columns are appropriately pre-processed either through being removed, converted to equivalent vector columns where appropriate, or some other appropriate treatment such as an informative error. This behaviour should be tested.*
#' @srrstatsTODO {G2.13} *Statistical Software should implement appropriate checks for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*
#' @srrstatsTODO {G2.14} *Where possible, all functions should provide options for users to specify how to handle missing (`NA`) data, with options minimally including:*
#' @srrstatsTODO {G2.14a} *error on missing data*
#' @srrstatsTODO {G2.14b} *ignore missing data with default warnings or messages issued*
#' @srrstatsTODO {G2.14c} *replace missing data with appropriately imputed values*
#' @srrstatsTODO {G2.15} *Functions should never assume non-missingness, and should never pass data with potential missing values to any base routines with default `na.rm = FALSE`-type parameters (such as [`mean()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/mean.html), [`sd()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/sd.html) or [`cor()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html)).*
#' @srrstatsTODO {G2.16} *All functions should also provide options to handle undefined values (e.g., `NaN`, `Inf` and `-Inf`), including potentially ignoring or removing such values.*
#' @srrstatsTODO {G3.0} *Statistical software should never compare floating point numbers for equality. All numeric equality comparisons should either ensure that they are made between integers, or use appropriate tolerances for approximate equality.*
#' @srrstatsTODO {G3.1} *Statistical software which relies on covariance calculations should enable users to choose between different algorithms for calculating covariances, and should not rely solely on covariances from the `stats::cov` function.*
#' @srrstatsTODO {G3.1a} *The ability to use arbitrarily specified covariance methods should be documented (typically in examples or vignettes).*
#' @srrstatsTODO {G5.0} *Where applicable or practicable, tests should use standard data sets with known properties (for example, the [NIST Standard Reference Datasets](https://www.itl.nist.gov/div898/strd/), or data sets provided by other widely-used R packages).*
#' @srrstatsTODO {G5.2} *Appropriate error and warning behaviour of all functions should be explicitly demonstrated through tests. In particular,*
#' @srrstatsTODO {G5.2a} *Every message produced within R code by `stop()`, `warning()`, `message()`, or equivalent should be unique*
#' @srrstatsTODO {G5.2b} *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
#' @srrstatsTODO {G5.3} *For functions which are expected to return objects containing no missing (`NA`) or undefined (`NaN`, `Inf`) values, the absence of any such values in return objects should be explicitly tested.*
#' @srrstatsTODO {G5.4a} *For new methods, it can be difficult to separate out correctness of the method from the correctness of the implementation, as there may not be reference for comparison. In this case, testing may be implemented against simple, trivial cases or against multiple implementations such as an initial R implementation compared with results from a C/C++ implementation.*
#' @srrstatsTODO {G5.4b} *For new implementations of existing methods, correctness tests should include tests against previous implementations. Such testing may explicitly call those implementations in testing, preferably from fixed-versions of other software, or use stored outputs from those where that is not possible.*
#' @srrstatsTODO {G5.4c} *Where applicable, stored values may be drawn from published paper outputs when applicable and where code from original implementations is not available*
#' @srrstatsTODO {G5.5} *Correctness tests should be run with a fixed random seed*
#' @srrstatsTODO {G5.6} **Parameter recovery tests** *to test that the implementation produce expected results given data with known properties. For instance, a linear regression algorithm should return expected coefficient values for a simulated data set generated from a linear model.*
#' @srrstatsTODO {G5.6a} *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
#' @srrstatsTODO {G5.6b} *Parameter recovery tests should be run with multiple random seeds when either data simulation or the algorithm contains a random component. (When long-running, such tests may be part of an extended, rather than regular, test suite; see G5.10-4.12, below).*
#' @srrstatsTODO {G5.7} **Algorithm performance tests** *to test that implementation performs as expected as properties of data change. For instance, a test may show that parameters approach correct estimates within tolerance as data size increases, or that convergence times decrease for higher convergence thresholds.*
#' @srrstatsTODO {G5.8} **Edge condition tests** *to test that these conditions produce expected behaviour such as clear warnings or errors when confronted with data with extreme properties including but not limited to:*
#' @srrstatsTODO {G5.8a} *Zero-length data*
#' @srrstatsTODO {G5.8b} *Data of unsupported types (e.g., character or complex numbers in for functions designed only for numeric data)*
#' @srrstatsTODO {G5.8c} *Data with all-`NA` fields or columns or all identical fields or columns*
#' @srrstatsTODO {G5.8d} *Data outside the scope of the algorithm (for example, data with more fields (columns) than observations (rows) for some regression algorithms)*
#' @srrstatsTODO {G5.9} **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
#' @srrstatsTODO {G5.9a} *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*
#' @srrstatsTODO {G5.9b} *Running under different random seeds or initial conditions does not meaningfully change results*
#' @srrstatsTODO {G5.10} *Extended tests should included and run under a common framework with other tests but be switched on by flags such as as a `<MYPKG>_EXTENDED_TESTS="true"` environment variable.* - The extended tests can be then run automatically by GitHub Actions for example by adding the following to the `env` section of the workflow:
#' @srrstatsTODO {G5.11} *Where extended tests require large data sets or other assets, these should be provided for downloading and fetched as part of the testing workflow.*
#' @srrstatsTODO {G5.11a} *When any downloads of additional data necessary for extended tests fail, the tests themselves should not fail, rather be skipped and implicitly succeed with an appropriate diagnostic message.*
#' @srrstatsTODO {G5.12} *Any conditions necessary to run extended tests such as platform requirements, memory, expected runtime, and artefacts produced that may need manual inspection, should be described in developer documentation such as a `CONTRIBUTING.md` or `tests/README.md` file.*
#' @srrstatsTODO {BS1.3a} *Bayesian Software should document, both in text and examples, how to use the output of previous simulations as starting points of subsequent simulations.*
#' @srrstatsTODO {BS1.4} *For Bayesian Software which implements or otherwise enables convergence checkers, documentation should explicitly describe and provide examples of use with and without convergence checkers.*
#' @srrstatsTODO {BS1.5} *For Bayesian Software which implements or otherwise enables multiple convergence checkers, differences between these should be explicitly tested.*
#' @srrstatsTODO {BS2.7} *Enable starting values to be explicitly controlled via one or more input parameters, including multiple values for software which implements or enables multiple computational "chains."*
#' @srrstatsTODO {BS2.8} *Enable results of previous runs to be used as starting points for subsequent runs.*
#' @srrstatsTODO {BS2.9} *Ensure each chain is started with a different seed by default.*
#' @srrstatsTODO {BS2.10} *Issue diagnostic messages when identical seeds are passed to distinct computational chains.*
#' @srrstatsTODO {BS2.11} *Software which accepts starting values as a vector should provide the parameter with a plural name: for example, "starting_values" and not "starting_value".*
#' @srrstatsTODO {BS2.12} *Bayesian Software should implement at least one parameter controlling the verbosity of output, defaulting to verbose output of all appropriate messages, warnings, errors, and progress indicators.*
#' @srrstatsTODO {BS2.13} *Bayesian Software should enable suppression of messages and progress indicators, while retaining verbosity of warnings and errors. This should be tested.*
#' @srrstatsTODO {BS2.14} *Bayesian Software should enable suppression of warnings where appropriate. This should be tested.*
#' @srrstatsTODO {BS2.15} *Bayesian Software should explicitly enable errors to be caught, and appropriately processed either through conversion to warnings, or otherwise captured in return values. This should be tested.*
#' @srrstatsTODO {BS3.1} *Implement pre-processing routines to diagnose perfect collinearity, and provide appropriate diagnostic messages or warnings*
#' @srrstatsTODO {BS3.2} *Provide distinct routines for processing perfectly collinear data, potentially bypassing sampling algorithms*
#' @srrstatsTODO {BS4.1} *Packages should provide explicit comparisons with external samplers which demonstrate intended advantage of implementation (generally via tests, vignettes, or both).*
#' @srrstatsTODO {BS4.4} *Enable computations to be stopped on convergence (although not necessarily by default).*
#' @srrstatsTODO {BS4.5} *Ensure that appropriate mechanisms are provided for models which do not converge.*
#' @srrstatsTODO {BS4.6} *Implement tests to confirm that results with convergence checker are statistically equivalent to results from equivalent fixed number of samples without convergence checking.*
#' @srrstatsTODO {BS4.7} *Where convergence checkers are themselves parametrised, the effects of such parameters should also be tested. For threshold parameters, for example, lower values should result in longer sequence lengths.*
#' @srrstatsTODO {BS5.0} *Return values should include starting value(s) or seed(s), including values for each sequence where multiple sequences are included*
#' @srrstatsTODO {BS5.1} *Return values should include appropriate metadata on types (or classes) and dimensions of input data*
#' @srrstatsTODO {BS5.2} *Bayesian Software should either return the input function or prior distributional specification in the return object; or enable direct access to such via additional functions which accept the return object as single argument.*
#' @srrstatsTODO {BS5.4} *Where multiple checkers are enabled, Bayesian Software should return details of convergence checker used*
#' @srrstatsTODO {BS5.5} *Appropriate diagnostic statistics to indicate absence of convergence should either be returned or immediately able to be accessed.*
#' @srrstatsTODO {BS6.0} *Software should implement a default `print` method for return objects*
#' @srrstatsTODO {BS6.1} *Software should implement a default `plot` method for return objects*
#' @srrstatsTODO {BS6.5} *Software may provide abilities to plot both sequences of posterior samples and distributional estimates together in single graphic*
#' @srrstatsTODO {BS7.0} *Software should demonstrate and confirm recovery of parametric estimates of a prior distribution*
#' @srrstatsTODO {BS7.1} *Software should demonstrate and confirm recovery of a prior distribution in the absence of any additional data or information*
#' @srrstatsTODO {BS7.2} *Software should demonstrate and confirm recovery of a expected posterior distribution given a specified prior and some input data*
#' @srrstatsTODO {BS7.3} *Bayesian software should include tests which demonstrate and confirm the scaling of algorithmic efficiency with sizes of input data.*
#' @srrstatsTODO {BS7.4} *Bayesian software should implement tests which confirm that predicted or fitted values are on (approximately) the same scale as input values.*
#' @srrstatsTODO {BS7.4a} *The implications of any assumptions on scales on input objects should be explicitly tested in this context; for example that the scales of inputs which do not have means of zero will not be able to be recovered.*
#' @noRd
NULL

#' NA_standards
#' @srrstatsNA {G1.4a} There are no non-exported functions.
#' @srrstatsNA {G2.5} There are no inputs of type 'factor' in this software.
#' @srrstatsNA {G4.0} There are no functions that write to local files
#' @srrstatsNA {BS1.3b} This software uses Stan that implements only one type of algorithm (HMC/NUTS)
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#' @noRd
NULL
