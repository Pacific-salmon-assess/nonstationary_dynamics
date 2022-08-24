data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale) - fixed in this
  real<upper = 0> b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] b_dev; //year-to-year deviations in a

}

transformed parameters{
  vector[L] log_b; //b in each year
  vector[L] b; //b in each year
  
  log_b[1] = b0;
  for(t in 2:L){
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a ~ normal(0,2.5); //initial productivity - wide prior
  b0 ~ normal(-12,3); //covariates - reef
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_b ~ gamma(2,3);
   
  b_dev ~ std_normal();
  
  for(n in 1:N) R_S[n] ~ normal(log_a-S[n]*b[ii[n]], sigma_e);
}
generated quantities{
  real b_3b;
  real b_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
  
  b_3b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]])/3;
  b_5b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]]+b[ii[N-3]]+b[ii[N-4]])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a - x_oos*b[ii[N]], sigma_e);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a - x_oos*b_3b, sigma_e);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a - x_oos*b_5b, sigma_e);
 }
 
 
 