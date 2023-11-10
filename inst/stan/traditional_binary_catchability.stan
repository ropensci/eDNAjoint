data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of samples
    array[C] int<lower=1> R;   // index of locations for traditional samples
    int<lower=1> Nloc;   // total number of locations
    array[C] int<lower=0> E;   // number of animals in sample C
    int<lower=0> nparams;  // number of gear types
    matrix[C,nparams] mat;  // matrix of gear type integers
    int<lower = 0, upper = 1> include_phi; // binary indicator of negbinomial

}

parameters{/////////////////////////////////////////////////////////////////////
    array[Nloc] real<lower=0> mu_1;  // expected density at each site
    vector<lower=-0.99999>[nparams] q_trans; // catchability coefficients
    vector<lower=0>[include_phi ? 1 : 0] phi;  // dispersion parameter, if include_phi = 1
    }

transformed parameters{/////////////////////////////////////////////////////////////////////
    vector<lower=0>[C] coef;

    for(i in 1:C){
      coef[i] = 1 + dot_product(mat[i],q_trans);
    }


}

model{/////////////////////////////////////////////////////////////////////


    if (include_phi == 1)
       for (j in 1:C){
        E[j] ~ poisson(coef[j]*mu_1[R[j]]); // Eq. 1.1
       }
    else
       for (j in 1:C){
        E[j] ~ neg_binomial_2(coef[j]*mu_1[R[j]], phi); // Eq. 1.1
       }


}

generated quantities{
  vector[nparams] q;
  vector[C] log_lik;
  matrix[Nloc,nparams+1] mu;  // matrix of catch rates

  q = q_trans + 1;

  mu[,1] = to_vector(mu_1);

  for(i in 1:nparams){
    mu[,i+1] = to_vector(mu_1)*q[i];
  }

    if (include_phi == 1)
       for (j in 1:C){
        log_lik[j] = poisson_lpmf(E[j] | coef[j]*mu_1[R[j]]); //store log likelihood of traditional data given model
       }
    else
       for (j in 1:C){
        log_lik[j] = neg_binomial_2_lpmf(E[j] | coef[j]*mu_1[R[j]], phi); //store log likelihood of traditional data given model
       }

}

