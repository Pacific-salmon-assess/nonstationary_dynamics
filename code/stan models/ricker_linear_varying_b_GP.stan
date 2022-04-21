data{
  int<lower=1> N;//number of annual samples (time-series length)
  real TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale)
  real<upper = 0> b0; // average rate capacity
  vector[N] z_yr; //annual deviations in slope
  
 //variance components  
  real<lower = 0> sigma_e;
  real<lower=0> gp_tau;
  real<lower=0> gp_rho;

}
transformed parameters{
	vector[N] log_b;
	vector[N] b;
	matrix[N,N] L_Tmat;
	matrix[N,N] Tmat;
	Tmat = cov_exp_quad(TT, gp_tau, gp_rho);
	for(n in 1:N) Tmat[n, n] += 1e-9;
	L_Tmat = cholesky_decompose(Tmat);

	log_b =b0+L_Tmat*z_yr;
	b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
}
model{
  //priors
  log_a ~ normal(0,5); //intrinsic productivity - wide prior
  b0 ~ student_t(5,-9,1); //average per capita capacity parameter
  
  //variance terms
  sigma_e ~ gamma(2,5);
  gp_tau ~ gamma(2,5);
  gp_rho ~ inv_gamma(5, 5);
  
  //scaled effects of time
  z_yr ~ std_normal();

  for(i in 1:N) R_S[i] ~ normal(log_a - b[i]*S[i], sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  vector[N] y_rep;
  for (t in 1:N){ log_lik[t] = normal_lpdf(R_S[t]|log_a - S[t]*b[t], sigma_e);
				  y_rep[t] = normal_rng(log_a - S[t]*b[t], sigma_e);
  }
  }
   

