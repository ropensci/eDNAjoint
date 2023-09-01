data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of samples
    int<lower=1> R[C];   // index of locations for samples
    int<lower=1> Nloc;   // total number of locations 
    int<lower=0> E[C];   // number of animals trapped in sample C
    
}

parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0> mu[Nloc];  // expected density at each site
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
  
    for(j in 1:C){
          log_lik[j] = poisson_lpmf(E[j] | mu[R[j]]); //store log likelihood of traditional data given model
      }
  
}

