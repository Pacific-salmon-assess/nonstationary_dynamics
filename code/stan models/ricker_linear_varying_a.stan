data{
  int<lower=1> N;//number of annual samples
  int<lower=1> TT;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters{
  real<lower = 0> log_a0;// initial productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;

  //time-varying parameters
  vector[TT-1] a_dev; //year-to-year deviations in a
  
}
transformed parameters{
  real b;
  vector[TT] log_a; //a in each year (on log scale)
  
  b=exp(log_b);
  
  log_a[1] = log_a0; //initial value
  for(t in 2:TT){
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
   
  for(t in 1:N){
    R_S[t] ~ normal(log_a[ii[t]] - S[t]*b, sigma_e);
  }
  
}