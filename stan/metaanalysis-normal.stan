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
}

parameters {
  // parameters for the overall effect distribution
  real mu;
  real<lower=0> tau;
  
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
  y ~ biased_normal(theta_0, beta_0, s_0, M);
}
