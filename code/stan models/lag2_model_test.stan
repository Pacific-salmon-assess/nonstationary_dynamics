data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters{
  real<upper = 0> log_b; // rate capacity - fixed in this
  
  vector<lower=0>[N] log_a; //a in each year (on log scale)
  vector[N-1] v; //lagged deviations in a
  
 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_v;
}
transformed parameters{
  real b;
 
  b=exp(log_b);
  
}  
model{
  //priors
   //initial productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  v[1] ~ normal(0, sigma_v);
   for(t in 2:N-1)
    v[t] ~ normal(v[t-1], sigma_v);
	
  log_a[1] ~ normal(0,2.5);
  for(t in 2:N)
    log_a[t] ~ normal(log_a[t-1]+ v[t-1], sigma_a); //random walk of log_a
 
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_v ~ gamma(2,3);
  sigma_a ~ gamma(2,3);
  
  for(t in 1:N){
    R_S[t] ~ normal(log_a[t] - S[t]*b, sigma_e);
  }
  
}
  generated quantities{
  vector[N] log_lik;
  vector[N] y_rep;
  for (t in 1:N){
  log_lik[t] = normal_lpdf(R_S[t]|log_a[t] - S[t]*b, sigma_e);
   y_rep[t] = normal_rng(log_a[t] - S[t]*b, sigma_e);
 }
 }
 
 
 