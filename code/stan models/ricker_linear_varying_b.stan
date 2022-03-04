data{
 int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real log_a;// initial productivity (on log scale) - fixed in this
  real b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[N-1] b_dev; //year-to-year deviations in a

}

transformed parameters{
  vector[N] log_b; //b in each year (log scale)
  vector[N] b; //b in each year
  
  log_b[1] = b0;
  for(t in 2:N){
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a ~ normal(0,10); //initial productivity - wide prior
  b0 ~ normal(0,2); //covariates - reef
  
  //variance terms
  sigma_e ~ inv_gamma(2, 1);
  sigma_b ~ inv_gamma(2, 1);
   
  b_dev ~ std_normal();
  for(t in 1:N)   R_S[t] ~ normal(log_a-b[t]*S[t], sigma_e);
}
  generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(R_S[t]|log_a-b[t]*S[t], sigma_e); 
  }