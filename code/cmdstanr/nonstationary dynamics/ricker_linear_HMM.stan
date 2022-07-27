data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int N_reg; //number of regimes
 }
parameters {
  real<lower = 0> log_a[N_reg];// initial productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  
  simplex[N_reg] theta[N]; //probability estimates for each regime by year
	
 
}
transformed parameters{
	real b;
	vector[N] reg[N]; //regimes by year
	
	for (t in 1:N){
	reg[t] ~ categorical_rng(theta[t]);
	}
	
	b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
}
model{
  //priors
  log_a ~ normal(0,2.5); //intrinsic productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  //variance terms
  sigma_e ~ gamma(2,3);
  
  theta[1] ~ dirichlet(rep_vector(1, N_reg));
	for(t in 2:N){
		theta[t] ~ dirichlet(theta[t-1]);
	}
  
 for(t in 1:N){
  R_S ~ normal(log_a[reg[t]] - S[t]*b, sigma_e);
  }
}
generated quantities{
  vector[N] log_lik;
  vector[N] y_rep;
  for (t in 1:N){ log_lik[t] = normal_lpdf(R_S[t]|log_a[reg[t]] - S[t]*b, sigma_e);
				  y_rep[t] = normal_rng(log_a[reg[t]] - S[t]*b, sigma_e);			
  }
 }
   
   