# Contributing to eDNAjoint

Thanks for considering contributing to eDNAjoint!

A few notes on contributing:

1. This package follows the [Tidyverse Style Guide](https://style.tidyverse.org/). Any contributed R code should pass `lintr::lint_package()`.
2. This package also contains [Stan code](https://mc-stan.org/) to specify the Bayesian models called by `joint_model()` and `traditional_model()`. 
These files are found in the [`inst/stan` folder](https://github.com/ropensci/eDNAjoint/tree/master/inst/stan) of this repo. 
Style of contributed Stan code can be checked using `lintr::lint(file.stan)`, which will follow a modification of the default lintr settings described [here](https://github.com/ropensci/eDNAjoint/tree/master/inst/stan/.lintr).

Feel free to email the author and maintainer (ednajoint@gmail.com) if you have any ideas or questions.
