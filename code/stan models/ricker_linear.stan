data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real log_a;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;

}
transformed parameters{
	real b;
	
	b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
}
model{
  //priors
  log_a ~ normal(0,5); //initial productivity - wide prior
  log_b ~ normal(0,5); //initial productivity - wide prior
  
  //variance terms
  sigma_e ~ inv_gamma(2, 1);

  R_S ~ normal(log_a - b*S, sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(R_S[t]|log_a - S[t]*b, sigma_e);
  }
   