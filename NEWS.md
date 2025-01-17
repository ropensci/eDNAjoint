# eDNAjoint 0.3

### User-facing changes

-   Changed function names to snakemake (i.e., `jointModel()` to `joint_model()`)
-   Changed function arguments to snakemake (i.e., `n.chain` to `n_chain`)
-   Changed input data names to snakemake (i.e., `qPCR.N` to `pcr_n`)
-   Changed MCMC sampling arguments so the `n_iter` includes all MCMC iterations, including the number of warm-up iterations specified in `n_warmup`. Changed default values accordingly.

### Non-user facing changes
-   Collapsed Stan files with more conditional statements to accommodate more model variations within one file
-   Added Stan helper functions
-   Moved `joint_model` and `traditional_model` helper functions to new files
-   Linted package to match default `lintr` style
-   Updated compilation flags to allow installation from source on Windows RStudio

# eDNAjoint 0.2

-   Package source code now hosted at <https://github.com/ropensci/eDNAjoint>
-   Improved programming (i.e., reduced cyclomatic complexity, spacing, etc.)
-   Added functionality to accommodate semi-paired data (i.e., traditional samples not collected at all surveyed sites).
-   Improved integration of package print statements with user guide.
-   Updated default calculation of initial values for MCMC.

# eDNAjoint 0.1

-   Initial version of the package on GitHub.
