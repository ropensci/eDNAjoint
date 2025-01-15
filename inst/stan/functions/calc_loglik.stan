/* functions for calculating log likelihood */

// calculate log likelihood of eDNA data
vector calc_loglik_dna(
  int n_S,
  int S_dna,
  int Nloc_dna,
  int[] n_K,
  int[] n_N,
  vector p_trad,
  int[] L_ind,
  int[] K_dna,
  int[] N_dna,
  real[] p_dna,
  int[] L_dna){

    vector[n_S + S_dna] log_lik;

    for (i in 1:n_S) {
      log_lik[i] = binomial_lpmf(n_K[i] | n_N[i], p_trad[L_ind[i]]);
    }

    if (Nloc_dna > 0) {
      //store log likelihood of eDNA data given model
      for (i in 1:S_dna) {
        log_lik[n_S+i] = binomial_lpmf(K_dna[i] | N_dna[i], p_dna[L_dna[i]]);
      }
    }

    return log_lik;
  }

// calculate log likelihood of traditional count data
vector calc_loglik_trad_count(
  real[] lambda,
  int negbin,
  real[] phi,
  int[] n_E,
  int n_C){

    vector[n_C] log_lik;

    //store log likelihood of traditional data given model
    if (negbin == 1) {
      for (j in 1:n_C) {
        log_lik[j] = neg_binomial_2_lpmf(n_E[j] | lambda[j], phi);
      }
    } else {
      for (j in 1:n_C) {
        log_lik[j] = poisson_lpmf(n_E[j] | lambda[j]);
      }
    }
    return log_lik;
  }

// calculate log likelihood of traditional continuous data
vector calc_loglik_trad_continuous(
  real[] lambda,
  vector beta_gamma,
  real[] E_trans,
  int[] R_ind,
  int n_C){

    vector[n_C] log_lik;

    //store log likelihood of traditional data given model
    for (j in 1:n_C) {
      log_lik[j] = gamma_lpdf(E_trans[j] | lambda, beta_gamma[R_ind[j]]);
    }

    return log_lik;
  }

// calculate log likelihood of data in joint count model
vector calc_loglik_count(
  int ctch,
  real[] coef,
  int[] mat,
  vector mu_trad,
  int[] R_ind,
  int negbin,
  real[] phi,
  int[] n_E,
  int[] n_K,
  int[] n_N,
  vector p_trad,
  int[] L_ind,
  int n_C,
  int n_S,
  int S_dna,
  int Nloc_dna,
  int[] K_dna,
  int[] N_dna,
  real[] p_dna,
  int[] L_dna){

    vector[n_C + n_S + S_dna] log_lik;

    // get lambda
    real lambda[n_C];
    lambda = get_lambda_count(ctch, coef, mat, mu_trad, R_ind, n_C);

    // traditional data
    log_lik[1:n_C] = calc_loglik_trad_count(lambda, negbin, phi, n_E, n_C);

    // eDNA data
    int end;
    end = n_C + n_S + S_dna;
    log_lik[n_C + 1:end] = calc_loglik_dna(n_S, S_dna, Nloc_dna, n_K,
                                           n_N, p_trad, L_ind, K_dna,
                                           N_dna, p_dna, L_dna);

    return log_lik;

  }

// calculate log likelihood of data in joint continuous model
vector calc_loglik_continuous(
  int ctch,
  real[] coef,
  int[] mat,
  vector alpha_gamma,
  vector beta_gamma,
  int[] R_ind,
  real[] E_trans,
  int[] n_K,
  int[] n_N,
  vector p_trad,
  int[] L_ind,
  int n_C,
  int n_S,
  int S_dna,
  int Nloc_dna,
  int[] K_dna,
  int[] N_dna,
  real[] p_dna,
  int[] L_dna){

    vector[n_C+n_S+S_dna] log_lik;

    // get lambda
    real lambda[n_C];
    lambda = get_lambda_continuous(ctch, coef, mat, alpha_gamma, R_ind, n_C);

    // traditional data
    log_lik[1:n_C] = calc_loglik_trad_continuous(lambda, beta_gamma, E_trans,
                                                 R_ind, n_C);

    // eDNA data
    int end;
    end = n_C + n_S + S_dna;
    log_lik[n_C + 1:end] = calc_loglik_dna(n_S, S_dna, Nloc_dna, n_K,
                                           n_N, p_trad, L_ind, K_dna,
                                           N_dna, p_dna, L_dna);

    return log_lik;

  }

// calculate lambda for count data
real[] get_lambda_count(
  int ctch,
  real[] coef,
  int[] mat,
  vector mu_trad,
  int[] R_ind,
  int n_C){

    real lambda[n_C];

    for (j in 1:n_C) {
      lambda[j] = (
        (ctch == 1) ? coef[mat[j]] * mu_trad[R_ind[j]] : mu_trad[R_ind[j]]
      );
    }

    return lambda;
}

// calculate lambda for continuous data
real[] get_lambda_continuous(
  int ctch,
  real[] coef,
  int[] mat,
  vector alpha_gamma,
  int[] R_ind,
  int n_C){

    real lambda[n_C];

    for (j in 1:n_C) {
      lambda[j] = (
        (ctch == 1) ? coef[mat[j]] *
        alpha_gamma[R_ind[j]] : alpha_gamma[R_ind[j]]
      );
    }

    return lambda;
}

// calculate log likelihood of data in traditional count model
vector calc_loglik_tradmod_count(
  int negbin,
  real[] phi,
  int[] n_E,
  int n_C,
  int ctch,
  real[] coef,
  int[] mat,
  vector mu_1,
  int[] R_ind){

    vector[n_C] log_lik;

    //get lambda
    real lambda[n_C];
    lambda = get_lambda_count(ctch, coef, mat, mu_1, R_ind, n_C);

    //store log likelihood of traditional data given model
    log_lik = calc_loglik_trad_count(lambda, negbin, phi, n_E, n_C);

    return log_lik;
  }

// calculate log likelihood of data in traditional continuous model
vector calc_loglik_tradmod_continuous(
  vector beta,
  real[] E_trans,
  int[] R_ind,
  int n_C,
  int ctch,
  real[] coef,
  int[] mat,
  vector alpha){

    vector[n_C] log_lik;

    //get lambda
    real lambda[n_C];
    lambda = get_lambda_continuous(ctch, coef, mat, alpha, R_ind, n_C);

    //store log likelihood of traditional data given model
    log_lik = calc_loglik_trad_continuous(lambda, beta, E_trans, R_ind, n_C);

    return log_lik;
  }
