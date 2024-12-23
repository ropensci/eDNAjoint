data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of trap samples
    array[C] int<lower=1> R;   // index of locations for traditional samples
    int<lower=1> Nloc;   // total number of locations
    array[C] real<lower=0> E;   // number of animals in sample C
    int<lower=0> nparams;  // number of gear types
    array[C] int<lower=1> mat;  // vector of gear type integers
    int<lower=0,upper=1> ctch; // binary indicator of presence of catchability coefficient

}

parameters{/////////////////////////////////////////////////////////////////////
    array[Nloc] real<lower=0> alpha;  // alpha param for gamma distribution
    array[Nloc] real <lower=0.01> beta;  // beta param for gamma distribution
    vector<lower=-0.99999>[nparams] q_trans; // catchability coefficients
    }

transformed parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0> coef[(ctch == 1) ? nparams+1 :  0];
    array[C] real<lower=0> E_trans;

    if(ctch == 1)
      coef = to_array_1d(append_row(1, 1+q_trans));

    for(j in 1:C){
      E_trans[j] = E[j] + 0.0000000000001;
    }
}

model{/////////////////////////////////////////////////////////////////////


    for (j in 1:C) {
      real lambda = (ctch == 1) ? coef[mat[j]]*alpha[R[j]] : alpha[R[j]];
      E_trans[j] ~ gamma(lambda, beta[R[j]]);  // Eq. 1.1
    }

    beta ~ gamma(0.25,0.25);
    alpha ~ gamma(0.01,0.01);

}

generated quantities{
  vector[nparams] q;
  vector[C] log_lik;
  matrix[Nloc,nparams+1] mu;  // matrix of catch rates

  ////////////////////////////////////
  // transform to interpretable params
  for(j in 1:Nloc){
    mu[j,1] = alpha[j]/beta[j];
  }

  if(ctch == 1)
    q = q_trans + 1;
    for(i in 1:nparams){
      mu[,i+1] = to_vector(mu[,1])*q[i];
    }

  ////////////////////////////////
  // get point-wise log likelihood

  //store log likelihood of traditional data given model
  for (j in 1:C) {
    real lambda = (ctch == 1) ? coef[mat[j]]*alpha[R[j]] : alpha[R[j]];
    log_lik[j] = gamma_lpdf(E_trans[j] | lambda, beta[R[j]]);
  }

}

