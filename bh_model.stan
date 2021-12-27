//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data{

  int<lower = 0> n; // number of observations

  vector[n] ssb; // vector of observed ssb

  vector[n] r; // vector of recruits

  real max_r; // max observed recruitment

}

transformed data{

  vector[n] log_r; // log recruitment

  log_r = log(r);

}

parameters {

  real<lower = 0.2, upper = 1> h; // steepness
  
  real<lower = 0> alpha; // max recruitment
  
  real<lower = 0> sigma; // recruitment standard deviation

}

transformed parameters{

  vector[n] rhat;
  
  vector[n] log_rhat;
  
  rhat = (0.8 * alpha * h * ssb) ./ (0.2 * alpha * (1 - h) +(h - 0.2) * ssb); // beverton holt model
  
  log_rhat = log(rhat);
  
}

model{

  log_r ~ normal(log_rhat - 0.5 * sigma^2, sigma); // account for retransformation bias
  
  sigma ~ cauchy(0, 2.5);
  
  alpha ~ normal(2 * max_r, 0.2 * max_r);

}

generated quantities{

  vector[n] pp_rhat;
  
  for (i in 1:n) {
  
    pp_rhat[i] = exp(normal_rng(log_rhat[i] - 0.5 * sigma^2, sigma)); // generate posterior predictives
    
  }
  
}

