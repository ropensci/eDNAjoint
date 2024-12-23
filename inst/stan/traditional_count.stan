data{/////////////////////////////////////////////////////////////////////
    int<lower=1> C;    // number of trap samples
    array[C] int<lower=1> R;   // index of locations for traditional samples
    int<lower=1> Nloc;   // total number of locations
    array[C] int<lower=0> E;   // number of animals in sample C
    int<lower=0> nparams;  // number of gear types
    array[C] int<lower=1> mat;  // vector of gear type integers
    array[2] real phipriors; // priors for gamma distrib on phi
    int<lower=0,upper=1> negbin; // binary indicator of negative binomial
    int<lower=0,upper=1> ctch; // binary indicator of presence of catchability coefficient

}

parameters{/////////////////////////////////////////////////////////////////////
    array[Nloc] real<lower=0> mu_1;  // expected density at each site
    real<lower=0> phi[(negbin == 1) ? 1 :  0]; // dispersion parameter
    vector<lower=-0.99999>[nparams] q_trans; // catchability coefficients
    }

transformed parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0> coef[(ctch == 1) ? nparams+1 :  0];

    if(ctch == 1)
      coef = to_array_1d(append_row(1, 1+q_trans));

}

model{/////////////////////////////////////////////////////////////////////


    if (negbin == 1) {
      for (j in 1:C) {
        real lambda = (ctch == 1) ? coef[mat[j]] * mu_1[R[j]] : mu_1[R[j]];
        E[j] ~ neg_binomial_2(lambda, phi);  // Eq. 1.1
      }
    } else {
      for (j in 1:C) {
        real lambda = (ctch == 1) ? coef[mat[j]] * mu_1[R[j]] : mu_1[R[j]];
        E[j] ~ poisson(lambda);  // Eq. 1.1
      }
    }

    if(negbin == 1)
      phi ~ gamma(phipriors[1], phipriors[2]); // phi prior

}

generated quantities{
  vector[nparams] q;
  vector[C] log_lik;
  matrix[Nloc,nparams+1] mu;  // matrix of catch rates

  ////////////////////////////////////
  // transform to interpretable params
  mu[, 1] = to_vector(mu_1);

  if(ctch == 1)
    q = q_trans + 1;
    mu[, 2:(nparams + 1)] = to_vector(mu_1) * q';

  ////////////////////////////////
  // get point-wise log likelihood

  //store log likelihood of traditional data given model
  if (negbin == 1) {
    for (j in 1:C) {
      real lambda = (ctch == 1) ? coef[mat[j]] * mu_1[R[j]] : mu_1[R[j]];
      log_lik[j] = neg_binomial_2_lpmf(E[j] | lambda, phi);
    }
  } else {
    for (j in 1:C) {
      real lambda = (ctch == 1) ? coef[mat[j]] * mu_1[R[j]] : mu_1[R[j]];
      log_lik[j] = poisson_lpmf(E[j] | lambda);
    }
  }

}

