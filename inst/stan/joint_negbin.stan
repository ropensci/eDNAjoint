data{/////////////////////////////////////////////////////////////////////
    int<lower=1> S;    // number of qPCR samples
    int<lower=1> C;    // number of traditional samples
    int<lower=1> L[S];   // index of locations for qPCR samples
    int<lower=1> R[C];   // index of locations for traditional samples
    int<lower=1> Nloc;   // number of locations
    int<lower=0> E[C];   // number of animals in sample C
    int<lower=1> N[S];   // number of qPCR replicates per site
    int<lower=0> K[S]; // number of qPCR detections among these replicates
    real p10priors[2]; // priors for normal distrib on p10

}

parameters{/////////////////////////////////////////////////////////////////////
    real<lower=0> mu[Nloc];  // expected catch at each site
    real<lower=0> phi;  // dispersion parameter
    real<lower=0> beta; // scaling coefficient in saturation function
    real<lower=0,upper=1> p10;  // p10, false-positive rate.
}

transformed parameters{/////////////////////////////////////////////////////////////////////
  real<lower=0, upper = 1> p11[Nloc]; // true-positive detection probability
  real<lower=0, upper = 1> p[Nloc];   // total detection probability
  
  for (i in 1:Nloc){
    p11[i] = mu[i] / (mu[i] + exp(beta)); // Eq. 1.2
    p[i] = p11[i] + p10; // Eq. 1.3
  }
}

model{/////////////////////////////////////////////////////////////////////


    for (j in 1:C){
        E[j] ~ neg_binomial_2(mu[R[j]], phi); // Eq. 1.1
    }
    
    for (i in 1:S){
        K[i] ~ binomial(N[i], p[L[i]]); // Eq. 1.4
    }
    
    
  //priors
  p10 ~ beta(p10priors[1], p10priors[2]); // p10 prior
  beta ~ normal(0,10); // beta shrinkage priors
  
}

generated quantities{
  vector[C] log_lik;
  
    for(j in 1:C){
          log_lik[j] = neg_binomial_2_lpmf(E[j] | mu[R[j]], phi); //store log likelihood of traditional data given model
      }
  
}

