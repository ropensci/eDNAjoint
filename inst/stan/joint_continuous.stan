functions {
  #include /functions/calc_loglik.stan
  #include /functions/calc_mu.stan
  #include /functions/calc_p11.stan
}


data{/////////////////////////////////////////////////////////////////////
    int<lower=1> n_S;    // number of paired qPCR samples
    int<lower=0> S_dna;    // number of unpaired qPCR samples
    int<lower=1> n_C;    // number of traditional samples
    array[n_S] int<lower=1> L_ind;   // index of locations for qPCR samples w/traditional pair
    array[S_dna] int<lower=1> L_dna;   // index of locations for unpaired qPCR samples
    array[n_C] int<lower=1> R_ind;   // index of locations for traditional samples
    int<lower=0> Nloc_dna;   // number of locations with unpaired dna samples
    int<lower=1> Nloc_trad;   // number of locations with paired samples
    array[Nloc_trad] int<lower=0> trad_ind;   // index of traditional samples
    array[Nloc_dna] int<lower=0> dna_ind;   // index of unpaired dna samples
    array[n_C] real<lower=0> n_E;   // number of animals in sample C
    array[n_S] int<lower=1> n_N;   // number of qPCR replicates per site
    array[n_S] int<lower=0> n_K; // number of qPCR detections among these replicates
    array[S_dna] int<lower=1> N_dna;   // number of qPCR replicates per site of unpaired samples
    array[S_dna] int<lower=0> K_dna; // number of qPCR detections among these replicates of unpaired samples
    array[2] real p10priors; // priors for normal distrib on p10
    array[2] real alphapriors; // priors for site covariate shrinkage params
    array[2] real bgammapriors; // priors for beta param in expected catch rate
    array[2] real agammapriors; // priors for alpha param in expected catch rate
    int<lower=0> nparams;  // number of gear types
    array[n_C] int<lower=1> mat;  // vector of gear type integers
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
  array[n_C] real<lower=0> E_trans;

  mu_trad = alpha_gamma ./ beta_gamma;
  p11_trad = calc_p11(Nloc_trad, mu_trad, mat_site, trad_ind, alpha); // Eq. 1.2
  p_trad = p11_trad + exp(log_p10); // Eq. 1.3

  if(ctch == 1)
    coef = to_array_1d(append_row(1, 1+q_trans));

  for(j in 1:n_C){
      E_trans[j] = n_E[j] + 0.0000000000001;
    }
}

model{/////////////////////////////////////////////////////////////////////


    // get lambda
    real lambda[n_C];
    lambda = get_lambda_continuous(ctch, coef, mat, alpha_gamma, R_ind, n_C);

    for (j in 1:n_C) {
      E_trans[j] ~ gamma(lambda[j], beta_gamma[R_ind[j]]);  // Eq. 1.1
    }

    for (i in 1:n_S){
      n_K[i] ~ binomial(n_N[i], p_trad[L_ind[i]]); // Eq. 1.4
    }

    if(Nloc_dna > 0)
       for (i in 1:S_dna){
         K_dna[i] ~ binomial(N_dna[i], p_dna[L_dna[i]]); // Eq. 1.4
       }


  //priors
  log_p10 ~ normal(p10priors[1], p10priors[2]); // p10 prior
  alpha ~ normal(alphapriors[1], alphapriors[2]); // sitecov shrinkage priors
  beta_gamma ~ gamma(bgammapriors[1], bgammapriors[2]);
  alpha_gamma ~ gamma(agammapriors[1], agammapriors[2]);

}

generated quantities{
  vector[nparams] q;
  vector[n_C+n_S+S_dna] log_lik;
  real p10;
  matrix[Nloc_dna+Nloc_trad,nparams+1] mu;  // matrix of catch rates
  vector[Nloc_trad] beta;

  ////////////////////////////////////
  // transform to interpretable params
  p10 = exp(log_p10);

  if(ctch == 1)
    q = q_trans + 1;

  beta = mat_site[trad_ind] * alpha;

  mu = calc_mu(trad_ind, dna_ind, mu_trad, ctch, nparams, q,
               Nloc_dna, Nloc_trad, p_dna, p10,
               mat_site, alpha);

  ////////////////////////////////
  // get point-wise log likelihood

  log_lik = calc_loglik_continuous(ctch, coef, mat, alpha_gamma, beta_gamma,
                                   R_ind, E_trans, n_K, n_N, p_trad, L_ind, n_C,
                                   n_S, S_dna, Nloc_dna, K_dna, N_dna, p_dna,
                                   L_dna);


}

