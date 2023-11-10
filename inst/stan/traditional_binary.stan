data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of samples
    array[C] int<lower=1> R;   // index of locations for traditional samples
    int<lower=1> Nloc;   // total number of locations
    array[C] int<lower=0> E;   // number of animals in sample C
    int<lower = 0, upper = 1> include_phi; // binary indicator of negbinomial

}

parameters{/////////////////////////////////////////////////////////////////////
    array[Nloc] real<lower=0> mu;  // expected density at each site
    vector<lower=0>[include_phi ? 1 : 0] phi;  // dispersion parameter, if include_phi = 1
    }

transformed parameters{/////////////////////////////////////////////////////////////////////


}

model{/////////////////////////////////////////////////////////////////////


    if (include_phi == 1)
       for (j in 1:C){
        E[j] ~ poisson(mu[R[j]]); // Eq. 1.1
       }
    else
       for (j in 1:C){
        E[j] ~ neg_binomial_2(mu[R[j]], phi); // Eq. 1.1
       }


}

generated quantities{
  vector[C] log_lik;

    if (include_phi == 1)
       for (j in 1:C){
        log_lik[j] = poisson_lpmf(E[j] | mu[R[j]]); //store log likelihood of traditional data given model
       }
    else
       for (j in 1:C){
        log_lik[j] = neg_binomial_2_lpmf(E[j] | mu[R[j]], phi); //store log likelihood of traditional data given model
       }

}

