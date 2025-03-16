/* function for calculating p11 */

// calculate p11 for joint models
vector calc_p11(
  int Nloc_trad,
  vector mu_trad,
  matrix mat_site,
  array[] int trad_ind,
  vector alpha){

    vector[Nloc_trad] p11_trad;

    // Eq. 1.2
    p11_trad = mu_trad ./ (mu_trad + exp(mat_site[trad_ind, ] * alpha));

    return p11_trad;
  }
