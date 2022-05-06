data{
  int<lower=1> N;//number of annual samples (time-series length)
  real TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real<lower = 0> log_a0;// average productivity (on log scale)
  real<upper = 0> log_b0; // average rate capacity (on log scale)
  vector[N] z_yr_a; //annual deviations in slope
  vector[N] z_yr_b; //annual deviations in slope
  
 //variance components  
  real<lower=0> sigma_e;
  real<lower=0> gp_tau_a;
  real<lower=0> gp_rho_a;
   real<lower=0> gp_tau_b;
  real<lower=0> gp_rho_b;

}
transformed parameters{
	vector[N] log_a;
	vector[N] log_b;
	vector[N] b;
	
	matrix[N,N] L_Tmat1;
	matrix[N,N] Tmat1;
	Tmat1 = cov_exp_quad(TT, gp_tau_a, gp_rho_b);
	for(n in 1:N) Tmat1[n, n] += 1e-9;
	L_Tmat1 = cholesky_decompose(Tmat1);
	
	matrix[N,N] L_Tmat2;
	matrix[N,N] Tmat2;
	Tmat2 = cov_exp_quad(TT, gp_tau_b, gp_rho_b);
	for(n in 1:N) Tmat2[n, n] += 1e-9;
	L_Tmat2 = cholesky_decompose(Tmat2);

	log_a =log_a0+L_Tmat1*z_yr_a;
	log_b =log_b0+L_Tmat2*z_yr_b;
	
	b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
}
model{
  //priors
  log_a0 ~ normal(0,5); //intrinsic productivity - wide prior
  log_b0 ~ student_t(5,-9,1); //average per capita capacity parameter
  
  //variance terms
  sigma_e ~ gamma(2,3);
  gp_tau_a ~ gamma(2,3);
  gp_rho_a ~ inv_gamma(5, 5);
  gp_tau_b ~ gamma(2,3);
  gp_rho_b ~ inv_gamma(5, 5);
  
  //scaled effects of time
  z_yr_a ~ std_normal();
  z_yr_b ~ std_normal();

  for(i in 1:N) R_S[i] ~ normal(log_a[i] - b[i]*S[i], sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  vector[N] y_rep;
  for (t in 1:N){ log_lik[t] = normal_lpdf(R_S[t]|log_a[t] - b[t]*S[t], sigma_e);
				  y_rep[t] = normal_rng(log_a[t] - b[t]*S[t], sigma_e);
  }
  }
   

