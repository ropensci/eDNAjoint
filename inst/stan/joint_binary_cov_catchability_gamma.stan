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
    array[Nloc] real<lower=0> alpha_gamma;  // alpha param for gamma distribution
    array[Nloc] real <lower=0.01> beta_gamma;  // beta param for gamma distribution
    real<upper=0> log_p10;  // p10, false-positive rate.
    vector[nsitecov] alpha; // site-level beta covariates
    vector<lower=-0.99999>[nparams] q_trans; // catchability coefficients
}

transformed parameters{/////////////////////////////////////////////////////////////////////
  array[Nloc] real<lower=0, upper = 1> p11; // true-positive detection probability
  array[Nloc] real<lower=0, upper = 1> p;   // total detection probability
  vector<lower=0>[C] coef;
  array[Nloc] real<lower=0> mu_1;  // expected catch at each site
  array[C] real<lower=0> E_trans;

  for (i in 1:Nloc){
    mu_1[i] = alpha_gamma[i]/beta_gamma[i];
    p11[i] = mu_1[i] / (mu_1[i] + exp(dot_product(mat_site[i],alpha))); // Eq. 1.2
    p[i] = p11[i] + exp(log_p10); // Eq. 1.3
  }

  for(k in 1:C){
      coef[k] = 1 + dot_product(mat[k],q_trans);
    }

    for(j in 1:C){
      E_trans[j] = E[j] + 0.0000000000001;
    }
}

model{/////////////////////////////////////////////////////////////////////


    for (j in 1:C){
      E_trans[j] ~ gamma(coef[j]*alpha_gamma[R[j]],beta_gamma[R[j]]); // Eq. 1.1
    }

    for (i in 1:S){
        K[i] ~ binomial(N[i], p[L[i]]); // Eq. 1.4
    }


  //priors
  log_p10 ~ normal(p10priors[1], p10priors[2]); // p10 prior
  alpha ~ normal(0,10); // sitecov shrinkage priors
  beta_gamma ~ gamma(0.25,0.25);
  alpha_gamma ~ gamma(0.01,0.01);

}

generated quantities{
  vector[nparams] q;
  vector[C+S] log_lik;
  vector[Nloc] beta;
  matrix[Nloc,nparams+1] mu;  // matrix of catch rates
  real p10;

  p10 = exp(log_p10);

  for (i in 1:Nloc){
    beta[i] = dot_product(mat_site[i],alpha);
  }

  q = q_trans + 1;

  mu[,1] = to_vector(mu_1);

  for(i in 1:nparams){
    mu[,i+1] = to_vector(mu_1)*q[i];
  }

    for(j in 1:C){
          log_lik[j] = gamma_lpdf(E_trans[j] | coef[j]*alpha_gamma[R[j]], beta_gamma[R[j]]); //store log likelihood of traditional data given model
      }

      for(i in 1:S){
          log_lik[C+i] = binomial_lpmf(K[i] | N[i], p[L[i]]); //store log likelihood of eDNA data given model
      }

}

