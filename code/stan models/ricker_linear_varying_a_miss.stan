data{
  int<lower=1> N;//number years in the data series(time-series length)
  int<lower=1> Nobs;//number of annual samples (N-missisng)
  int ii[Nobs];//index of years with data
  int TT[N];//index of years
  vector[Nobs] R_S; //log(recruits per spawner)
  vector[Nobs] S; //spawners in time T
 }
parameters{
  real<lower = 0> log_a0;// initial productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;

  //time-varying parameters
  vector[N-1] a_dev; //year-to-year deviations in a
  
}
transformed parameters{
  real b;
  vector[N] log_a; //a in each year (on log scale)
  
  b=exp(log_b);
  
  log_a[1] = log_a0; //initial value
  for(t in 2:N){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a; //random walk of log_a
  }
}  
model{
  //priors
  log_a0 ~ normal(0,2.5); //initial productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  a_dev ~ std_normal(); //standardized (z-scales) deviances
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_a ~ gamma(2,3);
   
  for(t in 1:Nobs){
    R_S[t] ~ normal(log_a[t] - S[t]*b, sigma_e);
  }
  
}
  generated quantities{
  vector[N] log_lik;
   vector[N] y_rep;
  for (t in 1:N){
  log_lik[t] = normal_lpdf(R_S[t]|log_a[ii[t]] - S[t]*b, sigma_e);
   y_rep[t] = normal_rng(log_a[ii[t]] - S[t]*b, sigma_e);
 }
 }
 
 
 