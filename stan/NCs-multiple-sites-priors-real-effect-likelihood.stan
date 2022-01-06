functions {

   // count_elem() and which_elem() from https://discourse.mc-stan.org/t/dealing-with-data-subsetting-in-stan/9842/7
   // count number times elem appears in test set
  int count_elem(vector test, real elem) {
    int count;
    count = 0;
    for(i in 1:num_elements(test))
    if(test[i] == elem)
    count = count + 1;
    return(count);
  }
  
  // find elements in test which are equal to elem
  int[] which_elem(vector test, real elem) {
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
  // doesn't deal with values outside the grid
  // assumes they are ordered, with no duplicates
  real grid_lpdf(real x, vector ll, vector vals_evaled, vector ll_indices, int n_val, real min_val, real max_val) {
    real lprob;
    real ll_lo;
    real ll_up;
    real floor_x;
    int floor_ind[1];
    int ceil_ind;
    real slope;
    
    // need to add these min and max vals to the data, as well as extend beyond just in case
    // find the index of the value right below x
    floor_x = floor((x-min_val)*n_val/(max_val-min_val)); 
    // this needs to be converted to an integer
    floor_ind = which_elem(ll_indices, floor_x);
    // and then the index above x
    ceil_ind = floor_ind[1] + 1;

    // LL at loside grid value
    ll_lo = ll[floor_ind[1]];
    // LL at upside grid value
    ll_up = ll[ceil_ind];
    
    // slope of LL function (assumed linear) in interval around x
    slope = (ll_up - ll_lo)/(vals_evaled[ceil_ind] - vals_evaled[floor_ind[1]]);
    // estimate LL at x based on this slope
    lprob = ll_lo + (x - vals_evaled[floor_ind[1]])*slope;
    
    return lprob;
  }
  
}
data {
  int<lower=0> J;          // number of negative controls over all sites
  real b[J];               // estimated effect on negative control
  real<lower=0> se_b[J];   // std err of negative control estimate
  int<lower=0> M;           // number of sites
  int<lower=0,upper=J> site[J];    // indicator of which site a given negative control belongs to
  
  int<lower=0> G; // grid size x M if all have the same grid size
  vector[G] ll;
  vector[G] vals_evaled; 
  vector[G] ll_indices;
  int<lower=0> grid_id[M,2]; // start and end indices of the sites along the likelihood grid
  int n_vals[M]; 
  vector[M] min_vals; 
  vector[M] max_vals;

  // parameters for priors
  // will be half-normal priors for the constrained positive params
  real mean_prior_THETA;
  real<lower=0> sd_prior_THETA;
  real<lower=0> mean_prior_GAMMA;
  real<lower=0> sd_prior_GAMMA;
  real mean_prior_BETA ;
  real<lower=0> sd_prior_BETA;
  real<lower=0> mean_prior_TAU ;
  real<lower=0> sd_prior_TAU;
  real mean_prior_beta[M];
  real<lower=0> sd_prior_beta[M];
  real<lower=0> mean_prior_tau[M];
  real<lower=0> sd_prior_tau[M];
}
parameters {
  real THETA;             // overall population true effect of interest
  real<lower=0> GAMMA;    // standard deviation in overall effects

  real theta_tilde_i[M];  // site-specific true biased effect
  real b_i[M];            // site specific true bias on true effect (so that true effect is theta_tilde_i - b_i)
  
  real true_b[J];         // true bias
  real beta[M];           // mean bias at site M
  real<lower=0> tau[M];   // standard deviation of bias distribution at site M
  real BETA;              // mean overall bias
  real<lower=0> TAU;      // standard deviation of mean bias across all sites
}
transformed parameters {
  real theta_i[M]; // true effect of interest
  real theta_b[M]; // center of the distribution from which theta_i is drawn
  for (i in 1:M) {
    theta_b[i] = THETA + b_i[i];
    theta_i[i] = theta_tilde_i[i] - b_i[i];
  }
}
model {
  // priors
  BETA ~ normal(mean_prior_BETA, sd_prior_BETA);
  TAU ~ normal(mean_prior_TAU, sd_prior_TAU);
  beta ~ normal(BETA, TAU); // distribution of mean bias across sites

  THETA ~ normal(mean_prior_THETA, sd_prior_THETA);
  GAMMA ~ normal(mean_prior_GAMMA, sd_prior_GAMMA);

  for (m in 1:M){
    beta[m] ~ normal(mean_prior_beta[m], sd_prior_beta[m]);
    tau[m] ~ normal(mean_prior_tau[m], sd_prior_tau[m]);
  }

  // negative control distribution
  for (j in 1:J){
    true_b[j] ~ normal(beta[site[j]], tau[site[j]]); // site-specific bias distribution
    b[j] ~ normal(true_b[j], se_b[j]); // sampling distribution
  }
  
  //effect of interest
  for (m in 1:M){
    theta_tilde_i[m] ~ normal(theta_b[m], GAMMA); // site-specific true effect
    theta_tilde_i[m] ~ grid(ll[grid_id[m,1]:grid_id[m,2]], // likelihood "function"
                            vals_evaled[grid_id[m,1]:grid_id[m,2]], 
                            ll_indices[grid_id[m,1]:grid_id[m,2]], 
                            n_vals[m], 
                            min_vals[m], 
                            max_vals[m]); 
    b_i[m] ~ normal(beta[m], tau[m]); // bias draw
  }
}
