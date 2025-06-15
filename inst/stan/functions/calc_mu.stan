/* functions for creating mu as a generated quantity */

// calculate mu for joint models
matrix calc_mu(
  array[] int trad_ind,
  array[] int dna_ind,
  vector mu_trad,
  int ctch,
  int nparams,
  vector q,
  int Nloc_dna,
  int Nloc_trad,
  array[] real p_dna,
  real p10,
  matrix mat_site,
  vector alpha){

    matrix[Nloc_dna + Nloc_trad, nparams + 1] mu;

    mu[trad_ind, 1] = mu_trad;
    if (ctch == 1) {
      mu[trad_ind, 2:(nparams + 1)] = mu_trad * q';
    }

    if (Nloc_dna > 0) {
      array[Nloc_dna] real p11_dna;
      for (i in 1:Nloc_dna) {
        p11_dna[i] = p_dna[i] - p10;
        mu[dna_ind[i], 1] = (
          p11_dna[i]*exp(dot_product(to_vector(mat_site[dna_ind[i]]), alpha)) /
          (1 - p11_dna[i])
          );
       }
      if (ctch == 1) {
        for (i in 1:Nloc_dna) {
          mu[dna_ind[i], 2:(nparams + 1)] = mu[dna_ind[i], 1] * q';
        }
      }
    }

    return mu;
  }


// calculate mu for traditional count model
matrix calc_mu_trad_count(
  int Nloc,
  int nparams,
  vector mu_1,
  vector q,
  int ctch){

    matrix[Nloc, nparams + 1] mu;

    mu[, 1] = mu_1;

    if (ctch == 1) {
      mu[, 2:(nparams + 1)] = mu_1 * q';
    }

    return mu;
  }

// calculate mu for traditional continuous model
matrix calc_mu_trad_continuous(
  int Nloc,
  int nparams,
  vector alpha,
  vector beta,
  vector q,
  int ctch){

    matrix[Nloc, nparams + 1] mu;

    for (j in 1:Nloc) {
      mu[j, 1] = alpha[j] / beta[j];
    }

    if (ctch == 1) {
      for (i in 1:nparams) {
        mu[, i + 1] = to_vector(mu[, 1]) * q[i];
      }
    }

    return mu;
  }
