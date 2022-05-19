data{
 int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale) - fixed in this
  real<lower = 0> b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[N-1] b_dev; //year-to-year deviations in a

}

transformed parameters{
  vector<lower = 0>[N] b; //b in each year
  
  b[1] = b0;
  for(t in 2:N){
    b[t] = b[t-1] + b_dev[t-1]*sigma_b;
  } 
}  

model{
  //priors
  log_a ~ normal(0,2.5); //initial productivity - wide prior
  b0 ~ normal(0,0.1); //
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_b ~ gamma(2,3);
   
  b_dev ~ std_normal();
  
  for(t in 1:N)   R_S[t] ~ normal(log_a-b[t]*S[t], sigma_e);
}
  generated quantities{
 vector[N] log_lik;
  vector[N] y_rep;
  for (t in 1:N){
  log_lik[t] = normal_lpdf(R_S[t]|log_a-b[t]*S[t], sigma_e);
   y_rep[t] = normal_rng(log_a-b[t]*S[t], sigma_e);
 }
 }
 
 
 