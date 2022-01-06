data {
  int<lower=0> J;          // number of negative controls over all sites
  real b[J];               // estimated effect on negative control
  real<lower=0> se_b[J];   // std err of effect
  int<lower=0> M;           // number of sites
  int<lower=0,upper=J> site[J];    // indicator of which site a given negative control belongs to
  
  // parameters for priors
  // will be half-normal priors for the constrained positive params
  real mean_prior_BETA ;
  real<lower=0> sd_prior_BETA;
  real<lower=0> mean_prior_TAU ;
  real<lower=0> sd_prior_TAU;
  real mean_prior_beta[M] ;
  real<lower=0> sd_prior_beta[M];
  real<lower=0> mean_prior_tau[M] ;
  real<lower=0> sd_prior_tau[M];
}
parameters {
  real true_b[J];           // true bias
  real beta[M];                // mean bias at site M
  real<lower=0> tau[M];        // standard deviation of bias distribution at site M
  real BETA;               // mean overall bias
  real<lower=0> TAU;      // standard deviation of mean bias across all sites
} 
model {
  BETA ~ normal(mean_prior_BETA, sd_prior_BETA);          // priors
  TAU ~ normal(mean_prior_TAU, sd_prior_TAU);        
  for (m in 1:M){
        beta[m] ~ normal(mean_prior_beta[m], sd_prior_beta[m]); 
        tau[m] ~ normal(mean_prior_tau[m], sd_prior_tau[m]);
  }

  beta ~ normal(BETA, TAU); // distribution of bias across sites
  for (j in 1:J){
      true_b[j] ~ normal(beta[site[j]], tau[site[j]]); // site-specific bias distribution
      b[j] ~ normal(true_b[j], se_b[j]); // sampling distribution
  }
}
