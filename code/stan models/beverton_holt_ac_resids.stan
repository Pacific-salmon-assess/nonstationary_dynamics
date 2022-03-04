data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R; //recruits
  vector[N] S; //spawners in time T
}
 transformed data{
 vector[N] log_R; //recruits in the log scale
 log_R = log(R);
 }
parameters {
  real log_a;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = -1, upper = 1> phi;
}

transformed parameters{
  real a = exp(log_a);
  real b = exp(log_b);
  vector[N] mu;
  vector[N] log_mu;
  vector[N] epsilon; //residuals

  mu = a*S ./(1 + b*S);
  log_mu = log(mu);
  
  epsilon[1] = log_R[1] - log_mu[1];
  for(t in 2:N){
    epsilon[t] =(log_R[t] - log_mu[t]);
    log_mu[t] = log_mu[t] + (phi*epsilon[t-1]);
  }

}
model{
  //priors
  log_a ~ normal(0,10); //initial productivity - wide prior
  log_b ~ normal(0,10); //initial productivity - wide prior
  phi ~ uniform(-1,1);
  
  //variance terms
  sigma_e ~ cauchy(0,5);

  log_R ~ normal(log_mu, sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(log_R|log_mu[t], sigma_e);
  } 