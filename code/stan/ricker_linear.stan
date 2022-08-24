data{
  int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;

}
transformed parameters{
	real b;
	
	b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
}
model{
  //priors
  log_a ~ normal(0,2.5); //intrinsic productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  //variance terms
  sigma_e ~ gamma(2,3);
  
  R_S ~ normal(log_a - S*b, sigma_e);
}

   