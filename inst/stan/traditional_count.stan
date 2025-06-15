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
  array[n_C] int<lower = 0> n_E;
  // number of gear types
  int<lower = 0> nparams;
  // vector of gear type integers
  array[n_C] int<lower = 1> mat;
  // priors for gamma distrib on phi
  array[2] real phi_priors;
  // binary indicator of negative binomial
  int<lower = 0, upper=1> negbin;
  // binary indicator of presence of catchability coefficient
  int<lower = 0, upper=1> ctch;
}

parameters {
  // expected density at each site
  vector<lower = 0>[Nloc] mu_1;
  // dispersion parameter
  array[negbin] real<lower = 0> phi;
  // catchability coefficients
  vector[nparams] q_log;
}

transformed parameters {
  array[(ctch == 1) ? nparams + 1 : 0] real<lower = 0> coef;

  if (ctch) {
    coef = to_array_1d(append_row(0, q_log));
  }

}

model {
  // get lambda
  array[n_C] real lambda;
  lambda = get_lambda_count(ctch, coef, mat, mu_1, R_ind, n_C);

  if (negbin) {
    for (j in 1:n_C) {
      n_E[j] ~ neg_binomial_2(lambda[j], phi);  // Eq. 1.1
    }
  } else {
    n_E ~ poisson(lambda);  // Eq. 1.1
  }

  if (negbin) {
    phi ~ gamma(phi_priors[1], phi_priors[2]); // phi prior
  }
}

generated quantities{
  vector[nparams] q;
  vector[n_C] log_lik;
  matrix[Nloc, nparams + 1] mu;  // matrix of catch rates

  ////////////////////////////////////
  // transform to interpretable params

  if (ctch) {
    q = exp(q_log);
  }

  mu = calc_mu_trad_count(Nloc, nparams, mu_1, q, ctch);

  ////////////////////////////////
  // get point-wise log likelihood
  log_lik = calc_loglik_tradmod_count(
    negbin, phi, n_E, n_C, ctch, coef, mat, mu_1, R_ind
  );
}
