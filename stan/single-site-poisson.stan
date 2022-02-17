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
  real<lower=0> counterfactualOutcomesInterest;
  real<lower=0> counterfactualDaysInterest;
  real<lower=0> exposureOutcomesInterest;
  real<lower=0> exposureDaysInterest;

  // parameters for priors
  // will be half-normal priors for the constrained positive params
  real mean_prior_THETA;
  real<lower=0> sd_prior_THETA;
}
parameters {
  real THETA;             // overall population true effect of interest
}

model {
  // priors for overall effect of interest
  THETA ~ normal(mean_prior_THETA, sd_prior_THETA);
  THETA ~ log_RR(counterfactualOutcomesInterest, counterfactualDaysInterest, //likelihood
                exposureOutcomesInterest, exposureDaysInterest);
}
