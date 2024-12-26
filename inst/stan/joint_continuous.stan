functions {
  #include /functions/calc_loglik.stan
}


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
    array[C] int<lower=1> mat;  // vector of gear type integers
    int<lower=0> nsitecov;  // number of site-level covariates
    matrix[Nloc_trad+Nloc_dna,nsitecov] mat_site;  // matrix of site-level covariates
    int<lower=0,upper=1> ctch; // binary indicator of presence of catchability coefficient

}

parameters{/////////////////////////////////////////////////////////////////////
    real<upper=0> log_p10;  // p10, false-positive rate.
    array[Nloc_dna] real<lower=0, upper = 1> p_dna;   // total detection probability
    vector<lower=-0.99999>[nparams] q_trans; // catchability coefficients
    vector[nsitecov] alpha; // site-level beta covariates
    vector<lower=0>[Nloc_trad] alpha_gamma;  // alpha param for gamma distribution
    vector<lower=0.01>[Nloc_trad] beta_gamma;  // beta param for gamma distribution
}

transformed parameters{/////////////////////////////////////////////////////////////////////
  vector<lower=0, upper = 1>[Nloc_trad] p11_trad; // true-positive detection probability
  vector<lower=0, upper = 1>[Nloc_trad] p_trad;   // total detection probability
  real<lower=0> coef[(ctch == 1) ? nparams+1 :  0];
  vector<lower=0>[Nloc_trad] mu_trad;  // expected catch at each site for sites with traditional samples
  array[C] real<lower=0> E_trans;

  mu_trad = alpha_gamma ./ beta_gamma;
  p11_trad = mu_trad ./ (mu_trad + exp(mat_site[trad_ind, ] * alpha)); // Eq. 1.2
  p_trad = p11_trad + exp(log_p10); // Eq. 1.3

  if(ctch == 1)
    coef = to_array_1d(append_row(1, 1+q_trans));

  for(j in 1:C){
      E_trans[j] = E[j] + 0.0000000000001;
    }
}

model{/////////////////////////////////////////////////////////////////////


    // get lambda
    real lambda[C];
    lambda = get_lambda_continuous(ctch, coef, mat, alpha_gamma, R, C);

    for (j in 1:C) {
      E_trans[j] ~ gamma(lambda[j], beta_gamma[R[j]]);  // Eq. 1.1
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
  alpha ~ normal(0,10); // sitecov shrinkage priors
  beta_gamma ~ gamma(0.25,0.25);
  alpha_gamma ~ gamma(0.01,0.01);

}

generated quantities{
  vector[nparams] q;
  vector[C+S+S_dna] log_lik;
  real p10;
  matrix[Nloc_dna+Nloc_trad,nparams+1] mu;  // matrix of catch rates
  array[Nloc_dna] real<lower=0, upper = 1> p11_dna; // true-positive detection probability
  vector[Nloc_trad] beta;

  ////////////////////////////////////
  // transform to interpretable params
  p10 = exp(log_p10);

  mu[trad_ind, 1] = mu_trad;
  if(ctch == 1)
    q = q_trans + 1;
    mu[trad_ind, 2:(nparams + 1)] = mu_trad * q';

  beta = mat_site[trad_ind] * alpha;


  if(Nloc_dna > 0)
     for (i in 1:Nloc_dna){
       p11_dna[i] = p_dna[i] - p10;
       mu[dna_ind[i],1] = p11_dna[i]*exp(dot_product(mat_site[dna_ind[i]],alpha))/(1-p11_dna[i]);
     }
     if(ctch == 1)
       for (i in 1:Nloc_dna){
         mu[dna_ind[i], 2:(nparams + 1)] = mu[dna_ind[i], 1] * q';
       }

  ////////////////////////////////
  // get point-wise log likelihood

  log_lik = calc_loglik_continuous(ctch, coef, mat, alpha_gamma, beta_gamma, R,
                                   E_trans, K, N, p_trad, L, C, S, S_dna,
                                   Nloc_dna, K_dna, N_dna, p_dna, L_dna);


}

