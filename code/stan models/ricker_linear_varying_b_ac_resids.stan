data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real a;// initial productivity (on log scale) - fixed in this
  real b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[N-1] b_dev; //year-to-year deviations in a
  vector[N] b; //b in each year 
}

transformed parameters{
  vector[N] mu;
  vector[N] epsilon;
  
  mu[1] = a + b0*S;
  b[1] = b0;
  epsilon[1] = R_S[1] - mu[1]; 
  for(t in 2:N){
    b[t] = b[t-1] + b_dev[t]*sigma_b;
    mu[t] = a + b[t]*S + phi*epsilon[t-1];
  }
}  

model{
  //priors
  a ~ normal(0,10); //initial productivity - wide prior
  beta0 ~ normal(0,2); //covariates - reef
  
  //variance terms
  sigma_e ~ inv_gamma(2, 1);
  sigma_a ~ inv_gamma(2, 1);
   
  b_dev ~ std_normal();
  
  R_S ~ normal(mu, sigma_e);
}
  generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpfg(R_S[t]|mu[t], sigma_e) 
  } 