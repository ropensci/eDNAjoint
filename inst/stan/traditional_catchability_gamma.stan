data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of trap samples
    array[C] int<lower=1> R;   // index of locations for traditional samples
    int<lower=1> Nloc;   // total number of locations
    array[C] real<lower=0> E;   // number of animals in sample C
    int<lower=0> nparams;  // number of gear types
    matrix[C,nparams] mat;  // matrix of gear type integers

}

parameters{/////////////////////////////////////////////////////////////////////
    array[Nloc] real<lower=0> alpha;  // alpha param for gamma distribution
    array[Nloc] real <lower=0.01> beta;  // beta param for gamma distribution
    vector<lower=-0.99999>[nparams] q_trans; // catchability coefficients
    }

transformed parameters{/////////////////////////////////////////////////////////////////////
    vector<lower=0>[C] coef;
    array[C] real<lower=0> E_trans;

    for(i in 1:C){
      coef[i] = 1 + dot_product(mat[i],q_trans);
    }

    for(j in 1:C){
      E_trans[j] = E[j] + 0.0000000000001;
    }
}

model{/////////////////////////////////////////////////////////////////////


    for(j in 1:C){

      E_trans[j] ~ gamma(coef[j]*alpha[R[j]],beta[R[j]]); // Eq. 1.1
    }

    beta ~ gamma(0.25,0.25);
    alpha ~ gamma(0.01,0.01);

}

generated quantities{
  vector[nparams] q;
  vector[C] log_lik;
  matrix[Nloc,nparams+1] mu;  // matrix of catch rates

  q = q_trans + 1;

  for(j in 1:Nloc){
    mu[j,1] = alpha[j]/beta[j];
  }

  for(i in 1:nparams){
    mu[,i+1] = to_vector(mu[,1])*q[i];
  }

    for(j in 1:C){
          log_lik[j] = gamma_lpdf(E_trans[j] | coef[j]*alpha[R[j]], beta[R[j]]); //store log likelihood of traditional data given model
      }

}

