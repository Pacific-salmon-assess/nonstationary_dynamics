data{
 int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real a0;// initial productivity (on log scale) - fixed in this
  real b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  vector<lower = 0>[2] sigma_a_b;
  
  //time-varying parameters
  cholesky_factor_corr[2] Lcorr;
  matrix[N-1,2] z_ab; //matrix of correlated a/b 
}

transformed parameters{
  matrix[N-1,2] z_mat; //matrix of correlated a/b 
  vector[N] log_a; //a in each year (log scale)
  vector[N] log_b; //b in each year (log scale)
  vector[N] b; //b in each year
  
  log_a[1] = a0;
  log_b[1] = b0;
	
	for(i in 1:N-1) 
    
  for(t in 2:N){
   z_mat[t-1,] = (diag_pre_multiply(sigma_a_b, Lcorr) * to_vector(z_ab[t-1,]))';
	log_a[t] = log_a[t-1] + z_mat[t-1,1];
    log_b[t] = log_b[t-1] + z_mat[t-1,2];
  } 
  b=exp(log_b);
}  

model{
  //priors
  a0 ~ normal(0,10); //initial productivity - wide prior
  b0 ~ normal(0,2); //covariates - reef
  Lcorr ~ lkj_corr_cholesky(1.0); //prior for cholesky factor
  
  //variance terms
  sigma_e ~ inv_gamma(2, 1);
  sigma_a_b ~ inv_gamma(2, 1);
  
  to_vector(z_ab) ~ std_normal();
  
  for(t in 1:N)   R_S[t] ~ normal(log_a[t]-b[t]*S[t], sigma_e);
}
  generated quantities{
   corr_matrix[2] Cor_1 = multiply_lower_tri_self_transpose(Lcorr);
   vector[N] log_lik;
  
  for (t in 1:N) log_lik[t] = normal_lpdf(R_S[t]|log_a[t]-b[t]*S[t], sigma_e); 
  }
 