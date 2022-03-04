data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real a0;// initial productivity (on log scale)
  real b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;
  
  //time-varying parameters
  vector[N-1] a_dev; //year-to-year deviations in a
  
  //correlation structures
  real<lower = -1, upper = 1> phi;
}
transformed parameters{
  vector[N] a; //a in each year (on log scale)
  vector[N] mu; //expected mean
  vector[N] epsilon; //residuals
   
   a[1] = a0;
   mu = a0-b*S;
   epsilon[1] = R_S[1] - mu[1];
  
  for(t in 2:N){
    a[t] = a[t-1] + a_dev[t-1]*sigma_a;
	mu[t] = a[t]-b*S[t] + (phi*epsilon[t-1]);
  }
}  

model{
  //priors
  a0 ~ normal(0,10); //initial productivity - wide prior
  b ~ normal(0, 2); //capacity rate
  a_dev ~ std_normal();
  phi ~ uniform(-1,1);
  
  //variance terms
  sigma_e ~ inv_gamma(2, 1);
  sigma_a ~ inv_gamma(2, 1);
 
  for(t in 1:N){
    R_S[t] ~ normal(mu[t], sigma_e);
  }
}

generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(R_S[t]|mu[t], sigma_e);
 }
 