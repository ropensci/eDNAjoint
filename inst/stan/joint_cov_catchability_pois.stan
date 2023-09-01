data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of qPCR samples
    int<lower=1> C;    // number of traditional samples
    int<lower=1> L[S];   // index of locations for qPCR samples
    int<lower=1> R[C];   // index of locations for traditional samples
    int<lower=1> Nloc;   // number of locations
    int<lower=0> E[C];   // number of animals in sample C
    int<lower=1> N[S];   // number of qPCR replicates per site
    int<lower=0> K[S]; // number of qPCR detections among these replicates
    int<lower=0> nsitecov;  // number of site-level covariates
    real p10priors[2]; // priors for normal distrib on p10
    matrix[Nloc,nsitecov] mat_site;  // matrix of site-level covariates
    int<lower=0> nparams;  // number of gear types
    matrix[C,nparams] mat;  // matrix of gear type integers

}

parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0> mu[Nloc];  // expected catch at each site
    real<lower=0,upper=1> p10;  // p10, false-positive rate.
    vector[nsitecov] alpha; // site-level beta covariates
    vector<lower=-0.99999>[nparams] q_trans; // catchability coefficients
}

transformed parameters{/////////////////////////////////////////////////////////////////////
  real<lower=0, upper = 1> p11[Nloc]; // true-positive detection probability
  real<lower=0, upper = 1> p[Nloc];   // total detection probability
  vector<lower=0>[C] coef;

  for (i in 1:Nloc){
    p11[i] = mu[i] / (mu[i] + exp(dot_product(mat_site[i],alpha))); // Eq. 1.2
    p[i] = p11[i] + p10; // Eq. 1.3
  }

  for(k in 1:C){
      coef[k] = 1 + dot_product(mat[k],q_trans);
    }
}

model{/////////////////////////////////////////////////////////////////////


    for (j in 1:C){
        E[j] ~ poisson(coef[j]*mu[R[j]]); // Eq. 1.1
    }

    for (i in 1:S){
        K[i] ~ binomial(N[i], p[L[i]]); // Eq. 1.4
    }


  //priors
  p10 ~ beta(p10priors[1], p10priors[2]); // p10 prior
  alpha ~ normal(0,10); // sitecov shrinkage priors

}

generated quantities{
  vector[nparams] q;
  vector[C] log_lik;
  vector[Nloc] beta;

  for (i in 1:Nloc){
    beta[i] = dot_product(mat_site[i],alpha);
  }

  q = q_trans + 1;

    for(j in 1:C){
          log_lik[j] = poisson_lpmf(E[j] | coef[j]*mu[R[j]]); //store log likelihood of traditional data given model
      }

}

