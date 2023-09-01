data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of samples
    int<lower=1> R[C];   // index of locations for samples
    int<lower=1> Nloc;   // total number of locations 
    int<lower=0> E[C];   // number of animals trapped in sample C
    int<lower=0> nparams;  // number of gear types
    matrix[C,nparams] mat;  // matrix of gear type integers
    
}

parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0> mu[Nloc];  // expected density at each site
    vector<lower=-0.99999>[nparams] q_trans; // catchability coefficients
    }

transformed parameters{/////////////////////////////////////////////////////////////////////
    vector<lower=0>[C] coef;
    
    for(i in 1:C){
      coef[i] = 1 + dot_product(mat[i],q_trans);
    }
    
  
}

model{/////////////////////////////////////////////////////////////////////

  
    for(j in 1:C){
      
      E[j] ~ poisson(coef[j]*mu[R[j]]); // Eq. 1.1
    }
    
  
}

generated quantities{
  vector[nparams] q;
  vector[C] log_lik;
  
  q = q_trans + 1;
  
    for(j in 1:C){
          log_lik[j] = poisson_lpmf(E[j] | coef[j]*mu[R[j]]); //store log likelihood of traditional data given model
      }
  
}

