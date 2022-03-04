data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R; //recruits
  vector[N] S; //spawners in time T
  real max_R; //maximum observed recruitment
}
 transformed data{
 vector[N] log_R; //recruits in the log scale
 log_R = log(R);
 }
parameters {
  real log_a;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this
  real<lower = 0.2, upper = 1> h; //steepness


 //variance components  
  real<lower = 0> sigma_e;
}
transformed parameters{
  real a = exp(log_a)
  real b = exp(log_b)
  vector[N] R_est;
  vector[N] log_R_est;
  
  for(n in 1:N){
    R_est[n] = (0.8*a*h*S[n])/(0.2*beta*(1-h) + (h-0.2)*S[n]);
  }
  
  log_R_est = log(R_est)
  
}
model{
  //priors
  a ~ normal(0,100); //initial productivity - wide prior
  b ~ normal(0,100); //initial productivity - wide prior
  
  //variance terms
  sigma_e ~ inv_gamma(2, 1);

  log_R ~ normal(log_R_est, sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(log_R|R_est, sigma_e);
  }