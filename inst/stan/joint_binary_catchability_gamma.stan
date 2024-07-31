data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of paired qPCR samples
    int<lower=0> S_dna;    // number of unpaired qPCR samples
    int<lower=1> C;    // number of traditional samples
    array[S] int<lower=1> L;   // index of locations for qPCR samples w/traditional pair
    array[S_dna] int<lower=1> L_dna;   // index of locations for unpaired qPCR samples
    array[C] int<lower=1> R;   // index of locations for traditional samples
    int<lower=0> Nloc_dna;   // number of locations with unpaired dna samples
    int<lower=1> Nloc_trad;   // number of locations with paired samples
    array[Nloc_trad] int<lower=0> trad_ind;   // index of traditional samples
    array[Nloc_dna] int<lower=0> dna_ind;   // index of unpaired dna samples
    array[C] real<lower=0> E;   // number of animals in sample C
    array[S] int<lower=1> N;   // number of qPCR replicates per site
    array[S] int<lower=0> K; // number of qPCR detections among these replicates
    array[S_dna] int<lower=1> N_dna;   // number of qPCR replicates per site of unpaired samples
    array[S_dna] int<lower=0> K_dna; // number of qPCR detections among these replicates of unpaired samples
    array[2] real p10priors; // priors for normal distrib on p10
    int<lower=0> nparams;  // number of gear types
    matrix[C,nparams] mat;  // matrix of gear type integers

}

parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0> beta; // scaling coefficient in saturation function
    real<upper=0> log_p10;  // p10, false-positive rate.
    array[Nloc_dna] real<lower=0, upper = 1> p_dna;   // total detection probability
    vector<lower=-0.99999>[nparams] q_trans; // catchability coefficients
    array[Nloc_trad] real<lower=0> alpha_gamma;  // alpha param for gamma distribution
    array[Nloc_trad] real <lower=0.01> beta_gamma;  // beta param for gamma distribution
}

transformed parameters{/////////////////////////////////////////////////////////////////////
  array[Nloc_trad] real<lower=0> mu_trad_1;  // expected catch at each site for sites with traditional samples
  array[Nloc_trad] real<lower=0, upper = 1> p11_trad; // true-positive detection probability
  array[Nloc_trad] real<lower=0, upper = 1> p_trad;   // total detection probability
  array[C] real<lower=0> E_trans;
  vector<lower=0>[C] coef;

  for (i in 1:Nloc_trad){
    mu_trad_1[i] = alpha_gamma[i]/beta_gamma[i];
    p11_trad[i] = mu_trad_1[i] / (mu_trad_1[i] + exp(beta)); // Eq. 1.2
    p_trad[i] = p11_trad[i] + exp(log_p10); // Eq. 1.3
  }

  for(k in 1:C){
      coef[k] = 1 + dot_product(mat[k],q_trans);
    }

  for(j in 1:C){
      E_trans[j] = E[j] + 0.0000000000001;
    }
}

model{/////////////////////////////////////////////////////////////////////


    for(j in 1:C){
      E_trans[j] ~ gamma(coef[j]*alpha_gamma[R[j]],beta_gamma[R[j]]); // Eq. 1.1
    }

    for (i in 1:S){
        K[i] ~ binomial(N[i], p_trad[L[i]]); // Eq. 1.4
    }

    if(Nloc_dna > 0)
       for (i in 1:S_dna){
           K_dna[i] ~ binomial(N_dna[i], p_dna[L_dna[i]]); // Eq. 1.4
       }


  //priors
  log_p10 ~ normal(p10priors[1], p10priors[2]); // p10 prior
  beta ~ normal(0,10);
  beta_gamma ~ gamma(0.25,0.25);
  alpha_gamma ~ gamma(0.01,0.01);

}

generated quantities{
  vector[nparams] q;
  vector[C+S+S_dna] log_lik;
  real p10;
  matrix[Nloc_dna+Nloc_trad,nparams+1] mu;  // matrix of catch rates
  array[Nloc_dna] real<lower=0, upper = 1> p11_dna; // true-positive detection probability

  p10 = exp(log_p10);

  q = q_trans + 1;

  for(i in 1:Nloc_trad){
    mu[trad_ind[i],1] = mu_trad_1[i];
    for(j in 1:nparams){
       mu[trad_ind[i],j+1] = mu_trad_1[i]*q[j];
    }
  }


  if(Nloc_dna > 0)
     for (i in 1:Nloc_dna){
       p11_dna[i] = p_dna[i] - p10;
       mu[dna_ind[i],1] = p11_dna[i]*exp(beta)/(1-p11_dna[i]);
       for(j in 1:nparams){
           mu[dna_ind[i],j+1] = mu[dna_ind[i],1]*q[j];
        }
     }

  for(j in 1:C){
    log_lik[j] = gamma_lpdf(E_trans[j] | coef[j]*alpha_gamma[R[j]], beta_gamma[R[j]]);
  }

  for(i in 1:S){
    log_lik[C+i] = binomial_lpmf(K[i] | N[i], p_trad[L[i]]); //store log likelihood of eDNA data given model
  }

  if(Nloc_dna > 0)
     for(i in 1:S_dna){
       log_lik[C+S+i] = binomial_lpmf(K_dna[i] | N_dna[i], p_dna[L_dna[i]]); //store log likelihood of eDNA data given model
     }

}

