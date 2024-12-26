/* functions for calculating log likelihood */

vector calc_loglik_dna(
  int S,
  int S_dna,
  int Nloc_dna,
  int[] K,
  int[] N,
  vector p_trad,
  int[] L,
  int[] K_dna,
  int[] N_dna,
  real[] p_dna,
  int[] L_dna){

    vector[S+S_dna] log_lik;

    for(i in 1:S){
      log_lik[i] = binomial_lpmf(K[i] | N[i], p_trad[L[i]]);
      }

    if(Nloc_dna > 0)
       //store log likelihood of eDNA data given model
       for(i in 1:S_dna){
         log_lik[S+i] = binomial_lpmf(K_dna[i] | N_dna[i], p_dna[L_dna[i]]);
         }

    return log_lik;
  }

vector calc_loglik_trad_count(
  real[] lambda,
  int negbin,
  real[] phi,
  int[] E,
  int C){

    vector[C] log_lik;

    //store log likelihood of traditional data given model
    if (negbin == 1) {
      for (j in 1:C) {
        log_lik[j] = neg_binomial_2_lpmf(E[j] | lambda[j], phi);
      }
    } else {
      for (j in 1:C) {
        log_lik[j] = poisson_lpmf(E[j] | lambda[j]);
      }
    }
    return log_lik;
  }

vector calc_loglik_trad_continuous(
  real[] lambda,
  vector beta_gamma,
  real[] E_trans,
  int[] R,
  int C){

    vector[C] log_lik;

    //store log likelihood of traditional data given model
    for (j in 1:C) {
      log_lik[j] = gamma_lpdf(E_trans[j] | lambda, beta_gamma[R[j]]);
    }

    return log_lik;
  }


vector calc_loglik_count(
  int ctch,
  real[] coef,
  int[] mat,
  vector mu_trad,
  int[] R,
  int negbin,
  real[] phi,
  int[] E,
  int[] K,
  int[] N,
  vector p_trad,
  int[] L,
  int C,
  int S,
  int S_dna,
  int Nloc_dna,
  int[] K_dna,
  int[] N_dna,
  real[] p_dna,
  int[] L_dna){

    vector[C+S+S_dna] log_lik;

    // get lambda
    real lambda[C];
    lambda = get_lambda_count(ctch, coef, mat, mu_trad, R, C);

    // traditional data
    log_lik[1:C] = calc_loglik_trad_count(lambda, negbin, phi, E, C);

    // eDNA data
    log_lik[C+1:C+S+S_dna] = calc_loglik_dna(S, S_dna, Nloc_dna, K, N,
                                             p_trad, L, K_dna, N_dna,
                                             p_dna, L_dna);

    return log_lik;

  }

vector calc_loglik_continuous(
  int ctch,
  real[] coef,
  int[] mat,
  vector alpha_gamma,
  vector beta_gamma,
  int[] R,
  real[] E_trans,
  int[] K,
  int[] N,
  vector p_trad,
  int[] L,
  int C,
  int S,
  int S_dna,
  int Nloc_dna,
  int[] K_dna,
  int[] N_dna,
  real[] p_dna,
  int[] L_dna){

    vector[C+S+S_dna] log_lik;

    // get lambda
    real lambda[C];
    lambda = get_lambda_continuous(ctch, coef, mat, alpha_gamma, R, C);

    // traditional data
    log_lik[1:C] = calc_loglik_trad_continuous(lambda, beta_gamma, E_trans, R, C);

    // eDNA data
    log_lik[C+1:C+S+S_dna] = calc_loglik_dna(S, S_dna, Nloc_dna, K, N,
                                             p_trad, L, K_dna, N_dna,
                                             p_dna, L_dna);

    return log_lik;

  }

real[] get_lambda_count(
  int ctch,
  real[] coef,
  int[] mat,
  vector mu_trad,
  int[] R,
  int C){

    real lambda[C];

    for (j in 1:C) {
      lambda[j] = (ctch == 1) ? coef[mat[j]] * mu_trad[R[j]] : mu_trad[R[j]];
    }

    return lambda;
}

real[] get_lambda_continuous(
  int ctch,
  real[] coef,
  int[] mat,
  vector alpha_gamma,
  int[] R,
  int C){

    real lambda[C];

    for (j in 1:C) {
      lambda[j] = (ctch == 1) ? coef[mat[j]]*alpha_gamma[R[j]] : alpha_gamma[R[j]];
    }

    return lambda;
}

