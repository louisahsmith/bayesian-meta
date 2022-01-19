functions {

  real log_RR_lpdf(real x, real counterfactualOutcomes, real counterfactualDays, real exposureOutcomes, real exposureDays) {
    real mu0;
    real lprob;
    
    mu0 = log(counterfactualOutcomes/counterfactualDays);
    
    lprob = exposureOutcomes*(mu0 + x + log(exposureDays)) - 
      exp(mu0 + x + log(exposureDays)) +
      counterfactualOutcomes*(mu0 + log(counterfactualDays)) - 
      exp(mu0 + log(counterfactualDays));
    
    return lprob;
  }
  
}
data {
  int<lower=0> J;          // number of negative controls over all sites
  int<lower=0,upper=J> site[J]; // indicator of which site a given negative control belongs to

  vector<lower=0>[J] counterfactualOutcomesNCs;
  vector<lower=0>[J] counterfactualDaysNCs;
  vector<lower=0>[J] exposureOutcomesNCs;
  vector<lower=0>[J] exposureDaysNCs;
  
  int<lower=0> M;           // number of sites
  vector<lower=0>[M] counterfactualOutcomesInterest;
  vector<lower=0>[M] counterfactualDaysInterest;
  vector<lower=0>[M] exposureOutcomesInterest;
  vector<lower=0>[M] exposureDaysInterest;

  // parameters for priors
  // will be half-normal priors for the constrained positive params
  real mean_prior_THETA;
  real<lower=0> sd_prior_THETA;
  real<lower=0> mean_prior_GAMMA;
  real<lower=0> sd_prior_GAMMA;
  vector[M] mean_prior_beta;
  vector<lower=0>[M] sd_prior_beta;
  vector<lower=0>[M] mean_prior_tau;
  vector<lower=0>[M] sd_prior_tau;
}
parameters {
  real THETA;             // overall population true effect of interest
  real<lower=0> GAMMA;    // standard deviation in overall effects

  real theta_tilde_i[M];  // site-specific true biased effect
  real b_i[M];            // site specific true bias on true effect (so that true effect is theta_tilde_i - b_i)
  
  real beta[M];           // mean bias at site M
  real<lower=0> tau[M];   // standard deviation of bias distribution at site M
  real b_j[J]; // bias of specific negative control 

}
transformed parameters {
  real theta_i[M]; // true effect of interest
  for (i in 1:M) {
    theta_i[i] = theta_tilde_i[i] - b_i[i];
  }
}
model {
  // priors for overall effect of interest
  THETA ~ normal(mean_prior_THETA, sd_prior_THETA);
  GAMMA ~ normal(mean_prior_GAMMA, sd_prior_GAMMA);

  // priors for the negative control distribution
  for (m in 1:M){
    beta[m] ~ normal(mean_prior_beta[m], sd_prior_beta[m]);
    tau[m] ~ normal(mean_prior_tau[m], sd_prior_tau[m]);
  }

  // negative control distribution
  for (j in 1:J){
    b_j[j] ~ normal(beta[site[j]], tau[site[j]]); // site-specific bias distribution
    b_j[j] ~ log_RR(counterfactualOutcomesNCs[j], counterfactualDaysNCs[j], exposureOutcomesNCs[j], exposureDaysNCs[j]); // likelihood
  }
  
  // effect of interest
  for (m in 1:M){
    theta_i[m] ~ normal(THETA, GAMMA); // site-specific true effect
    theta_tilde_i[m] ~ log_RR(counterfactualOutcomesInterest[m], counterfactualDaysInterest[m], //likelihood
                              exposureOutcomesInterest[m], exposureDaysInterest[m]);
    b_i[m] ~ normal(beta[m], tau[m]); // bias draw
  }
}
