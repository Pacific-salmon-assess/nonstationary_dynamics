data{
  int<lower=1> N;//number of annual samples (time-series length)
  int L; //total years covered by time-series
  int ii[N];//index of years with data
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
  vector[L-1] a_dev; //year-to-year deviations in a
  vector[L-1] b_dev; //year-to-year deviations in a
}

transformed parameters{
  vector[L] log_a; //a in each year (log scale)
  vector[L] log_b; //b in each year (log scale)
  vector[L] b; //b in each year
  
  log_a[1] = log_a0;
  log_b[1] = log_b0;
  for(t in 2:L){
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
  
  for(n in 1:N) R_S ~ normal(log_a[ii[n]]-b[ii[n]]*S[n], sigma_e);
}
generated quantities{
  real b_3b;
  real b_5b;
  real log_a_3b;
  real log_a_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
   
  log_a_3b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]])/3;
  log_a_5b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]]+log_a[ii[N-3]]+log_a[ii[N-4]])/5;
  
  b_3b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]])/3;
  b_5b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]]+b[ii[N-3]]+b[ii[N-4]])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a[ii[N]] - x_oos*b[ii[N]], sigma_e);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b_3b, sigma_e);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b_5b, sigma_e);
 }
 
 
 
 