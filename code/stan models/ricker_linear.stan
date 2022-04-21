data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
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
  log_a ~ normal(0,5); //intrinsic productivity - wide prior
  log_b ~ student_t(5,-9,1); //per capita capacity parameter - wide prior
  
  //variance terms
  sigma_e ~ student_t(7,0,1);

  R_S ~ normal(log_a - b*S, sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  vector[N] y_rep;
  for (t in 1:N){ log_lik[t] = normal_lpdf(R_S[t]|log_a - S[t]*b, sigma_e);
				  y_rep[t] = normal_rng(log_a - S[t]*b, sigma_e);
  }
  }
   