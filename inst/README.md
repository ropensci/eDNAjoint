
The `/stan` folder in this folder contains Bayesian model specifications
written in the Stan probabalistic programming language. Each file
corresponds to a variation of a model (originally developed in Keller et
al., 2022) that uses environmental DNA (eDNA) data and “traditional”
survey data to jointly estimate parameters. These model variations are
accessed based on the type of input data and/or user-defined input
parameters, including distributional assumptions.

Probability distributions were chosen for the model specifications using
the model developed in Keller et al. 2022. These original models use:

1.  a binomial distribution to represent the probability of a qPCR
    detection (1) or non-detection (0) (Lahoz-Monfort et al., 2016)
2.  a poisson or negative binomial distribution to represent how
    traditional survey count data are generated from an expected catch
    rate (Lindén and Mäntyniemi, 2011).
3.  a normal distribution as the prior on the probability of a false
    positive eDNA detection. This prior is informative and specified by
    the user as an input argument when running `joint_model()`
    (Lahoz-Monfort et al., 2016).
4.  a gamma distribution as the prior for the overdispersion parameter
    for traditional survey count data, if a negative binomial
    distribution is used. This prior can be informative and can be
    specified by the user as an input argument when running
    `joint_model()`.

Other variations on this original model specification include:

1.  a gamma distribution to represent how continuous traditional survey
    data are generated from an expected catch rate.
2.  a gear scaling coefficient if multiple traditional survey gear types
    are used and have different catchabilities
3.  a regression to include site-level covariates that scale the
    sensitivity of eDNA sampling relative to traditional surveys. A
    “shrinkage” prior is used for these covariates as a form of Bayesian
    penalization (van Erp et al., 2019).

This folder also contains ‘traditional models’, which can be used to
model the traditional survey data in isolation. These models can be used
as a comparison with the joint model that adds eDNA survey data to
determine if and how the addition of eDNA data affects inference.

The four files in the `/stan` folder represent four model variations:

1.  `joint_continuous.stan`: joint model with continuous traditional
    survey data
2.  `joint_count.stan`: joint model with discrete count traditional
    survey data
3.  `traditional_continuous.stan`: traditional model with continuous
    survey data
4.  `traditional_count.stan`: traditional model with discrete count
    survey data

The `/stan/functions` folder contains helper function for the above
files.

## References

Keller, A.G., Grason, E.W., McDonald, P.S., Ramon-Laca, A., Kelly, R.P.
(2022). Tracking an invasion front with environmental DNA. *Ecological
Applications*. 32(4): e2561. <https://doi.org/10.1002/eap.2561>

Lahoz-Monfort, J., Guillera-Arroita, G., Tingley, R. (2016). Statistical
approaches to account for false-positive errors in environmental DNA
samples. *Molecular Ecology Resources*. 16(3): 673-685.
<https://doi.org/10.1111/1755-0998.12486>

Lindén, A., Mäntyniemi, S. (2011). Using the negative binomial
distribution to model overdispersion in ecological count data.
*Ecology*. 92(7): 1414-1421. <https://doi.org/10.1890/10-1831.1>

van Erp, S., Oberski, D.L., Mulder, J. (2019). Shrinkage priors for
Bayesian penalized regression. *Journal of Mathematical Psychology*. 89:
31-50. <https://doi.org/10.1016/j.jmp.2018.12.004>
