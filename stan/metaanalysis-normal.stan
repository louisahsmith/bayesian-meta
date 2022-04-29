functions {
  real biased_normal_lpdf(real[] y, real[] theta, real[] beta, real[] s, int M) {
    real theta_tilde[M];
    for (i in 1:M){
      theta_tilde[i] = theta[i] + beta[i];
    }
    return normal_lpdf(y | theta_tilde, s);
  }
}

data {
  int<lower=1> M; // number of sites
  real y[M]; // log RRs
  real<lower=0> s_0[M]; // standard errors
  
  int<lower=1> N; // number of negative controls (overall)
  int<lower=1,upper=N> site[N]; // which site a negative control belongs to
  real x[N]; // log RRs
  real<lower=0> s_j[N]; // standard errors
  
  // prior distributions
  real mu_mean;
  real<lower=0> mu_sd;
  real tau_mean;
  real<lower=0> tau_sd;
  real lambda_mean;
  real<lower=0> lambda_sd;
  real eta_mean;
  real<lower=0> eta_sd;
  // can have separate priors for the sd
  // of the site-specific bias if you want
  real gamma_mean[M];
  real<lower=0> gamma_sd[M];
}

parameters {
  // parameters for the overall effect distribution
  real mu;
  real<lower=0> tau;
  
  // parameters for the overall bias distribution
  real lambda;
  real<lower=0> eta;
  
  //data source-specific parameters for the bias distribution
  real delta[M];
  real<lower=0> gamma[M];
  
  // true data source-specific effect
  real theta_0[M];
  
  // true biases for the effects of interest and for negative controls
  real beta_0[M];
  real beta_j[N];
}

model {
  // priors for overall effect
  mu ~ normal(mu_mean, mu_sd);
  tau ~ normal(tau_mean, tau_sd);
  
  // priors for bias distribution
  lambda ~ normal(lambda_mean, lambda_sd);
  eta ~ normal(eta_mean, eta_sd);
  
  // data-source specific bias distribution
  delta ~ normal(lambda, eta);
  gamma ~ normal(gamma_mean, gamma_sd);
  
  // negative controls
  for (j in 1:N){
    beta_j[j] ~ normal(delta[site[j]], gamma[site[j]]);
  }
  x ~ normal(beta_j, s_j);
  
  // effect of interest
  theta_0 ~ normal(mu, tau);
  beta_0 ~ normal(delta, gamma);
  y ~ biased_normal(theta_0, beta_0, s_0, M);
}
