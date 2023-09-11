data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of samples
    array[C] int<lower=1> R;   // index of locations for traditional samples
    int<lower=1> Nloc;   // total number of locations
    array[C] int<lower=0> E;   // number of animals in sample C
    int<lower=0> nparams;  // number of gear types
    matrix[C,nparams] mat;  // matrix of gear type integers

}

parameters{/////////////////////////////////////////////////////////////////////
    array[Nloc] real<lower=0> mu_1;  // expected density at each site
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

      E[j] ~ poisson(coef[j]*mu_1[R[j]]); // Eq. 1.1
    }


}

generated quantities{
  vector[nparams] q;
  vector[C] log_lik;
  matrix[Nloc,nparams+1] mu;  // matrix of catch rates

  q = q_trans + 1;

  mu[,1] = mu_1

  for(i in 1:length(nparams)){
    mu[,i+1] = mu_1*q[i]
  }

    for(j in 1:C){
          log_lik[j] = poisson_lpmf(E[j] | coef[j]*mu_1[R[j]]); //store log likelihood of traditional data given model
      }

}

