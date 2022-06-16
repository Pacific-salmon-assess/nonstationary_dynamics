data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R; //recruits
  vector[N] S; //spawners in time T
}
parameters {
  real a;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this
  real<lower=0> c; // depensation parameter

 //variance components  
  real<lower = 0> sigma_e;
}

transformed parameters{
  vector[N] mu;
  real b;
  
  b=exp(log_b);
  mu = a*S ./(1 + b*pow(S,c));
  
}
model{
  //priors
  a ~ normal(0,2.5); //initial productivity - wide prior
  log_b ~ normal(-12,3); //initial productivity - wide prior
  c ~ uniform(0,5);
  
  //variance terms
  sigma_e ~ inv_gamma(1, 1);

  R ~ lognormal(mu, sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = lognormal_lpdf(R|mu[t], sigma_e);
  } 