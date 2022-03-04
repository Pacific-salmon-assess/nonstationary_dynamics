data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters{
  real log_a0;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;

  //time-varying parameters
  vector[N] a_dev; //year-to-year deviations in a
  
}
transformed parameters{
  real b;
  vector[N] log_a; //a in each year (on log scale)
  
  b=exp(log_b);
  
  log_a[1] = log_a0; //initial value
  for(t in 2:N){
    log_a[t] = log_a[t-1] + a_dev[t]*sigma_a; //random walk of log_a
  }
}  
model{
  //priors
  log_a0 ~ normal(0,10); //initial productivity - wide prior
  log_b ~ normal(0, 2); //capacity rate
  a_dev ~ std_normal(); //standardized (z-scales) deviances
  
  //variance terms
  sigma_e ~ inv_gamma(2, 1);
  sigma_a ~ inv_gamma(2, 1);
   
  to_vector(a_dev) ~ std_normal();
  for(t in 1:N){
    R_S[t] ~ normal(log_a[t] - S[t]*b, sigma_e);
  }
  
}
  generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(R_S[t]|log_a[t] - S[t]*b, sigma_e);
 }
 
 