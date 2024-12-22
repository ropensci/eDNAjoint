data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of samples
    array[C] int<lower=1> R;   // index of locations for traditional samples
    int<lower=1> Nloc;   // total number of locations
    array[C] int<lower=0> E;   // number of animals in sample C

}

parameters{/////////////////////////////////////////////////////////////////////
    array[Nloc] real<lower=0> mu;  // expected density at each site
    }

transformed parameters{/////////////////////////////////////////////////////////////////////


}

model{/////////////////////////////////////////////////////////////////////


    for(j in 1:C){

      E[j] ~ poisson(mu[R[j]]); // Eq. 1.1
    }


}

generated quantities{
  vector[C] log_lik;

  ////////////////////////////////
  // get point-wise log likelihood
  for(j in 1:C){
      log_lik[j] = poisson_lpmf(E[j] | mu[R[j]]); //store log likelihood of traditional data given model
  }

}

