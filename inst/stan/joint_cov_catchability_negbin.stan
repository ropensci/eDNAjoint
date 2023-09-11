data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of qPCR samples
    int<lower=1> C;    // number of traditional samples
    array[S] int<lower=1> L;   // index of locations for qPCR samples
    array[C] int<lower=1> R;   // index of locations for traditional samples
    int<lower=1> Nloc;   // number of locations
    array[C] int<lower=0> E;   // number of animals in sample C
    array[S] int<lower=1> N;   // number of qPCR replicates per site
    array[S] int<lower=0> K; // number of qPCR detections among these replicates
    array[2] real p10priors; // priors for normal distrib on p10
    int<lower=0> nsitecov;  // number of site-level covariates
    matrix[Nloc,nsitecov] mat_site;  // matrix of site-level covariates
    int<lower=0> nparams;  // number of gear types
    matrix[C,nparams] mat;  // matrix of gear type integers

}

parameters{/////////////////////////////////////////////////////////////////////
    array[Nloc] real<lower=0> mu_1;  // expected catch at each site of gear type 1
    real<lower=0> phi;  // dispersion parameter
    real<lower=0,upper=1> p10;  // p10, false-positive rate.
    vector[nsitecov] alpha; // site-level beta covariates
    vector<lower=-0.99999>[nparams] q_trans; // catchability coefficients
}

transformed parameters{/////////////////////////////////////////////////////////////////////
  array[Nloc] real<lower=0, upper = 1> p11; // true-positive detection probability
  array[Nloc] real<lower=0, upper = 1> p;   // total detection probability
  vector<lower=0>[C] coef;

  for (i in 1:Nloc){
    p11[i] = mu_1[i] / (mu_1[i] + exp(dot_product(mat_site[i],alpha))); // Eq. 1.2
    p[i] = p11[i] + p10; // Eq. 1.3
  }

  for(k in 1:C){
      coef[k] = 1 + dot_product(mat[k],q_trans);
    }
}

model{/////////////////////////////////////////////////////////////////////


    for (j in 1:C){
        E[j] ~ neg_binomial_2(coef[j]*mu_1[R[j]], phi); // Eq. 1.1
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
  matrix[Nloc,nparams+1] mu;  // matrix of catch rates

  for (i in 1:Nloc){
    beta[i] = dot_product(mat_site[i],alpha);
  }

  q = q_trans + 1;

  mu[,1] = to_vector(mu_1);

  for(i in 1:nparams){
    mu[,i+1] = to_vector(mu_1)*q[i];
  }

    for(j in 1:C){
          log_lik[j] = neg_binomial_2_lpmf(E[j] | coef[j]*mu_1[R[j]], phi); //store log likelihood of traditional data given model
      }

}

