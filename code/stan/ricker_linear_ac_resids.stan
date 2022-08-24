data{
  int<lower=1> N;//number of annual samples
  int<lower=1> TT;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters{
  real log_a;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = -1, upper = 1> phi;

}
transformed parameters{
real b;
vector[N] mu;
vector[N] epsilon; //residuals

b = exp(log_b);
mu = log_a-b*S;

epsilon[1] = R_S[1] - mu[1];
  for(t in 2:N){
    epsilon[t] =(R_S[t] - mu[t]);
    mu[t] = mu[t] + (phi^(ii[t]-ii[t-1])*epsilon[t-1]); //phi raised the power of the number of time-steps between successive productivity estimates
  }
}
model{
  //priors
  log_a ~ normal(0,2.5); //initial productivity - wide prior
  log_b ~ normal(-12,3); //initial productivity - wide prior
  phi ~ uniform(-1,1);
  
  //variance terms
  sigma_e ~ gamma(2,3);

  R_S ~ normal(mu, sigma_e);
  
}

