data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of trap samples
    array[C] int<lower=1> R;   // index of locations for traditional samples
    int<lower=1> Nloc;   // total number of locations
    array[C] real<lower=0> E;   // number of animals in sample C

}

parameters{/////////////////////////////////////////////////////////////////////
    array[Nloc] real<lower=0> alpha;  // alpha param for gamma distribution
    array[Nloc] real <lower=0> beta;  // beta param for gamma distribution
    }

transformed parameters{/////////////////////////////////////////////////////////////////////
    array[C] real<lower=0> E_trans;  //

    for(j in 1:C){
      E_trans[j] = E[j] + 0.0000000000001;
    }
}

model{/////////////////////////////////////////////////////////////////////



    for(j in 1:C){

      E_trans[j] ~ gamma(alpha[R[j]],beta[R[j]]); // Eq. 1.1
    }

    beta ~ gamma(0.25,0.25);
    alpha ~ gamma(0.25,0.25);

}

generated quantities{
  vector[C] log_lik;
  vector[Nloc] mu;


  for(j in 1:Nloc){
    mu[j] = alpha[j]/beta[j];
  }

    for(j in 1:C){
          log_lik[j] = gamma_lpdf(E_trans[j] | alpha[R[j]], beta[R[j]]); //store log likelihood of traditional data given model
      }

}

