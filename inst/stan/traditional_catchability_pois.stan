data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of samples
    array[C] int<lower=1> R;   // index of locations for traditional samples
    int<lower=1> Nloc;   // total number of locations
    array[C] int<lower=0> E;   // number of animals in sample C
    int<lower=0> nparams;  // number of gear types
    array[C] int<lower=1> mat;  // vector of gear type integers

}

parameters{/////////////////////////////////////////////////////////////////////
    array[Nloc] real<lower=0> mu_1;  // expected density at each site
    vector<lower=-0.99999>[nparams] q_trans; // catchability coefficients
    }

transformed parameters{/////////////////////////////////////////////////////////////////////
    vector<lower=0>[nparams+1] coef;

    coef = append_row(1, 1+q_trans);


}

model{/////////////////////////////////////////////////////////////////////


    for(j in 1:C){

      E[j] ~ poisson(coef[mat[j]]*mu_1[R[j]]); // Eq. 1.1
    }


}

generated quantities{
  vector[nparams] q;
  vector[C] log_lik;
  matrix[Nloc,nparams+1] mu;  // matrix of catch rates

  ////////////////////////////////////
  // transform to interpretable params
  q = q_trans + 1;

  mu[, 1] = to_vector(mu_1);

  mu[, 2:(nparams + 1)] = to_vector(mu_1) * q';

  ////////////////////////////////
  // get point-wise log likelihood

  for(j in 1:C){
      log_lik[j] = poisson_lpmf(E[j] | coef[mat[j]]*mu_1[R[j]]); //store log likelihood of traditional data given model
  }

}

