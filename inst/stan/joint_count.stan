functions {
  #include /functions/calc_loglik.stan
  #include /functions/calc_mu.stan
  #include /functions/calc_p11.stan
}

data {
  // number of paired qPCR samples
  int<lower = 1> n_S;
  // number of unpaired qPCR samples
  int<lower = 0> S_dna;
  // number of traditional samples
  int<lower = 1> n_C;
  // index of locations for qPCR samples w/traditional pair
  array[n_S] int<lower = 1> L_ind;
  // index of locations for unpaired qPCR samples
  array[S_dna] int<lower = 1> L_dna;
  // index of locations for traditional samples
  array[n_C] int<lower = 1> R_ind;
  // number of locations with unpaired dna samples
  int<lower = 0> Nloc_dna;
  // number of locations with paired samples
  int<lower = 1> Nloc_trad;
  // index of traditional samples
  array[Nloc_trad] int<lower = 0> trad_ind;
  // index of unpaired dna samples
  array[Nloc_dna] int<lower = 0> dna_ind;
  // number of animals in sample C
  array[n_C] int<lower = 0> n_E;
  // number of qPCR replicates per site
  array[n_S] int<lower = 1> n_N;
  // number of qPCR detections among these replicates
  array[n_S] int<lower = 0> n_K;
  // number of qPCR replicates per site of unpaired samples
  array[S_dna] int<lower = 1> N_dna;
  // number of qPCR detections among these replicates of unpaired samples
  array[S_dna] int<lower = 0> K_dna;
  // priors for normal distrib on p10
  array[2] real p10_priors;
  // priors for site covariate shrinkage params
  array[2] real alphapriors;
  // number of gear types
  int<lower = 0> nparams;
  // vector of gear type integers
  array[n_C] int<lower = 1> mat;
  // number of site-level covariates
  int<lower = 0> nsitecov;
  // matrix of site-level covariates
  matrix[Nloc_trad + Nloc_dna, nsitecov] mat_site;
  // priors for gamma distrib on phi
  array[2] real phi_priors;
  // binary indicator of negative binomial
  int<lower = 0, upper = 1> negbin;
  // binary indicator of presence of catchability coefficient
  int<lower = 0, upper = 1> ctch;
}

parameters {
  // expected catch at each site for sites with traditional samples
  vector<lower = 0>[Nloc_trad] mu_trad;
  // p10, false-positive rate
  real<upper = 0> log_p10;
  // total detection probability
  array[Nloc_dna] real<lower = 0, upper = 1> p_dna;
  // catchability coefficients
  vector[nparams] q_log;
  // site-level beta covariates
  vector[nsitecov] alpha;
  // dispersion parameter
  array[negbin] real<lower = 0> phi;
}

transformed parameters {
  // true-positive detection probability
  vector<lower = 0, upper = 1>[Nloc_trad] p11_trad;
  // total detection probability
  vector<lower = 0, upper = 1>[Nloc_trad] p_trad;
  // traditional sample-specific catchability coefficient
  array[(ctch == 1) ? nparams + 1 : 0] real coef;

  p11_trad = calc_p11(Nloc_trad, mu_trad, mat_site, trad_ind, alpha); // Eq. 1.2
  p_trad = p11_trad + exp(log_p10); // Eq. 1.3

  if (ctch) {
    coef = to_array_1d(append_row(0, q_log));
  }
}

model {
  // get lambda
  array[n_C] real lambda;
  lambda = get_lambda_count(ctch, coef, mat, mu_trad, R_ind, n_C);

  if (negbin) {
    for (j in 1:n_C) {
      n_E[j] ~ neg_binomial_2(lambda[j], phi);  // Eq. 1.1
    }
  } else {
    n_E ~ poisson(lambda);  // Eq. 1.1
  }

  n_K ~ binomial(n_N, p_trad[L_ind]); // Eq. 1.4

  if (Nloc_dna > 0) {
    K_dna ~ binomial(N_dna, p_dna[L_dna]); // Eq. 1.4
  }

  //priors
  log_p10 ~ normal(p10_priors[1], p10_priors[2]); // p10 prior
  alpha ~ normal(alphapriors[1], alphapriors[2]); // sitecov shrinkage priors
  if (negbin) {
    phi ~ gamma(phi_priors[1], phi_priors[2]); // phi prior
  }
}

generated quantities {
  vector[nparams] q;
  vector[n_C + n_S + S_dna] log_lik;
  real p10;
  matrix[Nloc_dna + Nloc_trad, nparams + 1] mu;
  vector[Nloc_trad] beta;

  ////////////////////////////////////
  // transform to interpretable params

  p10 = exp(log_p10);

  if (ctch) {
    q = exp(q_log);
  }

  beta = mat_site[trad_ind] * alpha;

  mu = calc_mu(
    trad_ind, dna_ind, mu_trad, ctch, nparams, q,
    Nloc_dna, Nloc_trad, p_dna, p10, mat_site, alpha
  );

  ////////////////////////////////
  // get point-wise log likelihood

  log_lik = calc_loglik_count(
    ctch, coef, mat, mu_trad, R_ind, negbin, phi,
    n_E, n_K, n_N, p_trad, L_ind, n_C, n_S, S_dna,
    Nloc_dna, K_dna, N_dna, p_dna, L_dna
  );
}
