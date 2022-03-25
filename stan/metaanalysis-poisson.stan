functions {
  // Poisson likelihood for the number of outcomes in terms of the log rate ratio theta
  // and the bias beta given a fixed expected number of outcomes in the comparison group
  real poisson_RR_lpmf(int[] y, real[] theta, real[] beta, real[] y_star) {
    int N = num_elements(theta);
    real theta_tilde[N];
    real lambda[N];
    for (i in 1:N) {
      theta_tilde[i] = theta[i] + beta[i];
      lambda[i] = exp(theta_tilde[i]) * y_star[i];
    }
    return poisson_lpmf(y | lambda);
  }
}

data {
  int<lower=1> M; // number of sites
  int y[M]; // y_i
  real<lower=0> y_star[M]; // y*

  int<lower=1> N; // number of negative controls (overall)
  int<lower=1,upper=N> site[N]; // which site a negative control belongs to
  int x[N]; // x_ij
  real<lower=0> x_star[N]; // x*_ij
  real zeros[N]; //literally just 0s
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
  x ~ poisson_RR(zeros, beta_j, x_star);

  // effect of interest
  theta_0 ~ normal(mu, tau);
  beta_0 ~ normal(delta, gamma);
  y ~ poisson_RR(theta_0, beta_0, y_star);
}
