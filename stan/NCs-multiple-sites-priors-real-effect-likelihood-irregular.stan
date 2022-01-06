functions {

  int bisect_left(real x, vector vals, int hi_ind) {
    int mid;
    int lo = 1;
    int hi = hi_ind;
    
    while(lo < hi) {
      mid = (lo + hi) / 2;
      if (vals[mid] < x) lo = mid + 1;
      else hi = mid;
    }
    return lo;
  }
  
  // interpolate from a sequence of values from log-likelihood function
  // doesn't deal with values outside the grid
  // assumes they are ordered, with no duplicates
  real grid_lpdf(real x, vector ll, vector vals_evaled, int hi_ind) {
    real lprob;
    real ll_lo;
    real ll_up;
    int floor_ind;
    int ceil_ind;
    real slope;
    
    // index below
    floor_ind = bisect_left(x, vals_evaled, hi_ind);
    // and then the index above x
    ceil_ind = floor_ind + 1;

    // LL at loside grid value
    ll_lo = ll[floor_ind];
    // LL at upside grid value
    ll_up = ll[ceil_ind];
    
    // slope of LL function (assumed linear) in interval around x
    slope = (ll_up - ll_lo)/(vals_evaled[ceil_ind] - vals_evaled[floor_ind]);
    // estimate LL at x based on this slope
    lprob = ll_lo + (x - vals_evaled[floor_ind])*slope;
    
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
  // vector[G] ll_indices; 
  int<lower=0> grid_id[M,2]; // start and end indices of the sites along the likelihood grid
  int<lower=0> n_vals[M];

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
                            vals_evaled[grid_id[m,1]:grid_id[m,2]], n_vals[m]); 
    b_i[m] ~ normal(beta[m], tau[m]); // bias draw
  }
}
