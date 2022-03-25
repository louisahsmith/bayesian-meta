functions {
  // count_elem() and which_elem() adapted slightly from https://discourse.mc-stan.org/t/dealing-with-data-subsetting-in-stan/9842/7
  // count number times elem appears in test set
  int count_elem(real[] test, real elem) {
    int count;
    count = 0;
    for(i in 1:num_elements(test))
    if(test[i] == elem)
    count = count + 1;
    return(count);
  }
  
  // find elements in test which are equal to elem
  int[] which_elem(real[] test, real elem) {
    int res[count_elem(test, elem)];
    int ci;
    ci = 1;
    for(i in 1:num_elements(test))
    if(test[i] == elem) {
      res[ci] = i;
      ci = ci + 1;
    }
    return(res);
  }
  
  // interpolate from a sequence of values from log-likelihood function
  // doesn't deal with values outside the grid (proposal will be rejected with warning)
  // assumes they are ordered, with no duplicates
  real grid_lp(vector theta, vector beta, vector[] ll, real[] vals_evaled, real[] ll_indices, int M, int n_vals, real min_val, real max_val) {
    real x; // theta + beta proposal
    real floor_x; // index of closest value below x
    int floor_ind[1]; // index as integer
    int ceil_ind; // index of value above x
    real val_lo; // actual values evaluated at those indices
    real val_up;
    real ll_lo; // log-likelihoods at those indices
    real ll_up;
    real lprob = 0; // collect the sum of the log-likelihoods

    for (i in 1:M) {
      x = theta[i] + beta[i];
      if (x < min_val || x > max_val)
          reject("proposal not in range");
      // find the index of the value right below x
      floor_x = floor((x-min_val)*n_vals/(max_val-min_val)); 
      // this needs to be converted to an integer
      floor_ind = which_elem(ll_indices, floor_x);
      // and then the index above x
      ceil_ind = floor_ind[1] + 1;
      // LL at loside grid value
      ll_lo = ll[floor_ind[1]][i];
      // LL at upside grid value
      ll_up = ll[ceil_ind][i];
      val_lo = vals_evaled[floor_ind[1]];
      val_up = vals_evaled[ceil_ind];
      // estimate LL at x based on this linear slope within interval around x
      lprob = ll_lo + (x - val_lo)*(ll_up - ll_lo)/(val_up - val_lo) + lprob;
    }
    return lprob;
  }
}

data {
  int<lower=1> M; // number of sites
  int<lower=1> n_vals; // number of likelihood values per site (grid size)
  vector[M] ll[n_vals]; // likelihood values (array of size M containing vectors with G elements)
  real vals_evaled[n_vals]; // values at which likelihood was evaluated
  real<lower=1> ll_indices[n_vals]; // 1:G; just to avoid creating it
  real min_val; // min & max values at which likelihood evaluated
  real max_val;
  
  int<lower=1> N; // number of negative controls (overall)
  int<lower=1,upper=N> site[N]; // which site a negative control belongs to
  real x[N]; // log RRs
  real<lower=0> s_j[N]; // standard errors
}

parameters {
  // parameters for the overall effect distribution
  real mu;
  real<lower=0> tau;
  
  //data source-specific parameters for the bias distribution
  real delta[M];
  real<lower=0> gamma[M];
  
  // true data source-specific effect
  vector[M] theta_0;
  
  // true biases for the effects of interest and for negative controls
  vector[M] beta_0;
  vector[N] beta_j;
}

model {
  // priors for overall effect
  mu ~ normal(0, 10);
  tau ~ normal(0, 5);
  
  // data-source specific bias priors
  delta ~ normal(0, 10);
  gamma ~ normal(0, 5);
  
  // negative controls
  for (j in 1:N){
    beta_j[j] ~ normal(delta[site[j]], gamma[site[j]]);
  }
  x ~ normal(beta_j, s_j);
  
  // effect of interest
  theta_0 ~ normal(mu, tau);
  beta_0 ~ normal(delta, gamma);
  target += grid_lp(theta_0, beta_0, ll, vals_evaled, ll_indices, M, n_vals, min_val, max_val);
}
