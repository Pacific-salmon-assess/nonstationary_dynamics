data{
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
  real<lower=0> log_a0;// initial productivity (on log scale) - fixed in this
  real<upper=0> log_b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[N-1] a_dev; //year-to-year deviations in a
  vector[N-1] b_dev; //year-to-year deviations in a
}

transformed parameters{
  vector[N] log_a; //a in each year (log scale)
  vector[N] log_b; //b in each year (log scale)
  vector[N] b; //b in each year
  
  log_a[1] = log_a0;
  log_b[1] = log_b0;
  for(t in 2:N){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a;
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a0 ~ normal(0,2.5); //initial productivity - wide prior
  log_b0 ~ normal(-12,3); //covariates - reef
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_a ~ gamma(2,3);
  sigma_b ~ gamma(2,3);
  
  a_dev ~ std_normal();
  b_dev ~ std_normal();
  
  for(t in 1:N)   R_S[t] ~ normal(log_a[t]-b[t]*S[t], sigma_e);
}
generated quantities{
  vector[N] log_lik;
  real b_3b;
  real b_5b;
  real log_a_3b;
  real log_a_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
  
  for (t in 1:N) log_lik[t] = normal_lpdf(R_S[t]|log_a[t]-b[t]*S[t], sigma_e);
  
  log_a_3b = (log_a[N]+log_a[N-1]+log_a[N-2])/3;
  log_a_5b = (log_a[N]+log_a[N-1]+log_a[N-2]+log_a[N-3]+log_a[N-4])/5;
  
  b_3b = (b[N]+b[N-1]+b[N-2])/3;
  b_5b = (b[N]+b[N-1]+b[N-2]+b[N-3]+b[N-4])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a[N] - x_oos*b[N], sigma_e);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b_3b, sigma_e);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b_5b, sigma_e);
 }
 
 
 
 