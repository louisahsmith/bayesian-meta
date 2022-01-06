data {
  int<lower=0> J;          // # negative controls
  real b[J];               // estimated effect
  real<lower=0> se_b[J];   // std err of effect
}
parameters {
  real true_b[J];           // true bias
  real beta;                // mean bias
  real<lower=0> tau;        // standard deviation of bias distribution
} 
model {
  true_b ~ normal(beta, tau);
  b ~ normal(true_b, se_b); // sampling distribution
}
