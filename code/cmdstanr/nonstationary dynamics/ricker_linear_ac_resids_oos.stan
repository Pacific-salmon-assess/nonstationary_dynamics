data{
  int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T

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
    mu[t] = mu[t] + (phi*epsilon[t-1]);
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
generated quantities{
  vector[N] log_lik;
  real ep_3b;
  real ep_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
  
  for (t in 1:N) log_lik[t] = normal_lpdf(R_S[t]|mu[t], sigma_e);
  
  ep_3b = (epsilon[N]+epsilon[N-1]+epsilon[N-2])/3;
  ep_5b = (epsilon[N]+epsilon[N-1]+epsilon[N-2]+epsilon[N-3]+epsilon[N-4])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a - x_oos*b+phi*epsilon[N], sigma_e);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a - x_oos*b+phi*ep_3b, sigma_e);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a - x_oos*b+phi*ep_5b, sigma_e);
 } 


