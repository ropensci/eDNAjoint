data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of trap samples
    int<lower=1> R[C];   // index of locations for trap samples
    int<lower=1> Nloc;   // total number of locations 
    int<lower=0> E[C];   // number of crabs trapped in trap sample C
    
}

parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0> mu[Nloc];  // expected density at each site 
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

