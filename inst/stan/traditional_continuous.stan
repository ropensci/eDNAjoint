functions {
  #include /functions/calc_loglik.stan
  #include /functions/calc_mu.stan
}

data {
  // number of trap samples
  int<lower = 1> n_C;
  // index of locations for traditional samples
  array[n_C] int<lower = 1> R_ind;
  // total number of locations
  int<lower = 1> Nloc;
  // number of animals in sample C
  array[n_C] real<lower = 0> n_E;
  // number of gear types
  int<lower = 0> nparams;
  // vector of gear type integers
  array[n_C] int<lower = 1> mat;
  // binary indicator of presence of catchability coefficient
  int<lower = 0, upper = 1> ctch;
  // priors for beta param in expected catch rate
  array[2] real betapriors;
  // priors for alpha param in expected catch rate
  array[2] real alphapriors;
}

parameters {
  // alpha param for gamma distribution
  vector<lower = 0>[Nloc] alpha;
  // beta param for gamma distribution
  vector<lower = 0.01>[Nloc] beta;
  // catchability coefficients
  vector<lower = -0.99999>[nparams] q_trans;
}

transformed parameters {
  // traditional sample-specific catchability coefficient
  array[(ctch == 1) ? nparams + 1 :  0] real<lower = 0> coef;
  // transformed traditional data so that E > 0
  array[n_C] real<lower = 0> E_trans;

  if (ctch == 1) {
    coef = to_array_1d(append_row(1, 1 + q_trans));
  }

  for (j in 1:n_C) {
    E_trans[j] = n_E[j] + 0.0000000000001;
  }
}

model {
  // get lambda
  array[n_C] real lambda;
  lambda = get_lambda_continuous(ctch, coef, mat, alpha, R_ind, n_C);

  for (j in 1:n_C) {
    E_trans[j] ~ gamma(lambda[j], beta[R_ind[j]]);  // Eq. 1.1
  }

  alpha ~ gamma(alphapriors[1], alphapriors[2]);
  beta ~ gamma(betapriors[1], betapriors[2]);
}

generated quantities{
  vector[nparams] q;
  vector[n_C] log_lik;
  matrix[Nloc, nparams + 1] mu;  // matrix of catch rates

  ////////////////////////////////////
  // transform to interpretable params
  if (ctch == 1) {
    q = q_trans + 1;
  }

  mu = calc_mu_trad_continuous(Nloc, nparams, alpha, beta, q, ctch);

  ////////////////////////////////
  // get point-wise log likelihood

  // get lambda
  log_lik = calc_loglik_tradmod_continuous(
    beta, E_trans, R_ind, n_C, ctch, coef, mat, alpha
  );
}
