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
  real<lower = 0> a;// initial productivity (on log scale)
  real<lower = 0> b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
}
transformed parameters{
  vector[N] log_mu;

  
  for(n in 1:N){
    log_mu[n] = log(S)[n]+a-b*S[n];
  }
}
model{
  //priors
  a ~ cauchy(0, 2); //initial productivity - wide prior
  b ~ cauchy(0, 10); //initial productivity - wide prior
  
  //variance terms
  sigma_e ~ inv_gamma(1, 1);

  log_R ~ normal(log_mu, sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(log_R|log_mu[t], sigma_e);
}