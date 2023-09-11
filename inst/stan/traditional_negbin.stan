data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of trap samples
    array[C] int<lower=1> R;   // index of locations for traditional samples
    int<lower=1> Nloc;   // total number of locations
    array[C] int<lower=0> E;   // number of animals in sample C

}

parameters{/////////////////////////////////////////////////////////////////////
    array[Nloc] real<lower=0> mu;  // expected density at each site
    real<lower=0> phi;  // dispersion parameter
    }

transformed parameters{/////////////////////////////////////////////////////////////////////

}

model{/////////////////////////////////////////////////////////////////////


    for(j in 1:C){

      E[j] ~ neg_binomial_2(mu[R[j]],phi); // Eq. 1.1
    }


}

generated quantities{
  vector[C] log_lik;

    for(j in 1:C){
          log_lik[j] = neg_binomial_2_lpmf(E[j] | mu[R[j]], phi); //store log likelihood of traditional data given model
      }

}

