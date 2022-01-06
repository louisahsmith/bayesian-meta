data {
  int<lower=0> J;          // number of negative controls over all sites
  real b[J];               // estimated effect on negative control
  real<lower=0> se_b[J];   // std err of negative control estimate
  int<lower=0> M;           // number of sites
  int<lower=0,upper=J> site[J];    // indicator of which site a given negative control belongs to
  
  real theta_hat[M]; //biased estimate of effect of interst
  real se_theta[M]; // std err of effect of interest
  
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
  real THETA;                 // overall population true effect of interest
  real<lower=0> GAMMA;        // standard deviation in overall effects

  real theta_tilde_i[M];      // site-specific true biased effect
  real b_i[M];                // site specific true bias on true effect (so that true effect is theta_tilde_i - b_i)
  
  real true_b[J];           // true bias
  real beta[M];             // mean bias at site M
  real<lower=0> tau[M];    // standard deviation of bias distribution at site M
  real BETA;               // mean overall bias
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
    theta_hat[m] ~ normal(theta_tilde_i[m], se_theta[m]); // sampling distribution
    b_i[m] ~ normal(beta[m], tau[m]); // bias draw
  }
}
