data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int<lower=1> R;//number of ocean regions
  int<lower=1> SR; //mumber of subregions
  array[N] int J_i;//index of stocks
  array[SR] int sr_or;//index of ocean-region by sub-region
  array[J] int j_or;//index of ocean region by stock  
  array[J] int j_sr;//index of sub-region by stock  
  array[N] real R_S; //log(recruits per spawner)
  array[N] real S; //design matrix of spawners in time T
  int<lower=1> Nb; //number of B-splines
  matrix[N,Nb] B_year; //B-spline matrix
  matrix[L,Nb] B_pred_year; //B-spline prediction matrix to get the mean function
  array[J] real pSmax_mean; //priors on smax - based on observed spawner abundance
  array[J] real pSmax_sig;
}
parameters{
//long-term alpha
 real<lower=0> log_a0;// global productivity mean
  vector[R] z_a_or; //ocean region-level (WC/GOA/BS) productivity deviation
  vector[SR] z_a_sr;//sub-region productivity deviation
  vector[J] z_a_j;//stock-specific productivity deviation
  real<lower=0> sd_a_or; //variance in ocean-region productivity
  real<lower=0> sd_a_sr; //variance in sub-region productivity by ocean region
  real<lower=0> sd_a; //variance in stock productivity by sub-region

//density dependence 
vector<lower = 0>[J] Smax; // spawners where recruitment is maximized
//variance components  
 vector<lower = 0>[J] sigma; //total variance
//Hierarchical splines
vector[Nb] z_w;
real<lower=0> tau_w;
matrix[J,Nb] z_w_j; //stock-specific z-scores for each basis function
matrix[R,Nb] z_w_r; //ocean region-specific z-scores for each basis function
vector<lower=0>[J] tau_wj; //sigma for stock-specific taus
vector<lower=0>[R] tau_wr; //sigma for stock-specific taus - degree to which 
}
transformed parameters{
vector[J] b=1.0./Smax; //capacity rate
vector[Nb] w; //b-spline global weight
matrix[R,Nb] w_r; //stock-specific b-spline weights productivity
matrix[J,Nb] w_j; //stock-specific b-spline weights productivity
vector[N] mu; //expectation for each observation
 
 //productivity estimates
  vector[R] log_a_or; //ocean region specific productivity
  vector[SR] log_a_sr; //sub-region specific productivity
  vector[J] log_a_j; //stock-specific productivity
   
   //productivity terms
   log_a_or = log_a0 + z_a_or*sd_a_or; //ocean-region-specific productivity
   log_a_sr = log_a_or[sr_or] + z_a_sr*sd_a_sr; //subregion-specific productivity
   log_a_j = log_a_sr[j_sr] + z_a_j*sd_a; //stock-specific productivity

 //global GAM weights as a function of z-scores x shared variance tau_w
 w = z_w*tau_w; 
 
//stock-specific weight deviates
for(k in 1:Nb){
  w_r[,k]=w[k] + z_w_r[,k].*tau_wr;
  w_j[,k]= w_r[j_or,k] + z_w_j[,k].*tau_wj;
}

for(n in 1:N) mu[n]=log_a_j[J_i[n]]+B_year[n,]*to_vector(w_j[J_i[n],])-b[J_i[n]]*S[n];
}
model{
  
  //long-term prod
  log_a0 ~ normal(1.5,2.5); //initial productivity - wide prior
  z_a_or ~ std_normal(); //N(0,1)
  z_a_sr ~ std_normal(); //N(0,1)
  z_a_j ~ std_normal(); //N(0,1)
  sd_a_or ~ normal(0,1);
  sd_a_sr ~ normal(0,1);
  sd_a ~ normal(0,1);

   //capacity for each stock
  for(j in 1:J) Smax[j] ~ normal(pSmax_mean[j],pSmax_sig[j]);
  
  //GAM terms
  z_w ~ normal(0,1);
  to_vector(z_w_r) ~ normal(0,1);
  to_vector(z_w_j) ~ normal(0,1);
  tau_w ~ normal(0,3);
  tau_wr ~ normal(0,3);
  tau_wj ~ normal(0,3);
 
  //variance terms
  sigma ~ gamma(2,1);
 
  R_S ~ normal(mu,sigma[J_i]); 
  
}
generated quantities{
vector[L] log_a_t; //predicted global productivity per year
matrix[L,R] log_a_rt; //non linear prediction
matrix[L,J] log_a_jt; //non linear prediction
matrix[L,J] Smsy_jt;
matrix[L,J] Umsy_jt;
matrix[L,R] Umsy_rt;
vector[L] Umsy_t;
vector[N] log_lik;

for(l in 1:L){
    log_a_t[l]=log_a0+B_pred_year[l,]*to_vector(w);
    Umsy_t[l] = 1-lambert_w0(exp(1-log_a_t[l]));
	
	for(r in 1:R){
	log_a_rt[l,r]=log_a_or[r]+B_pred_year[l,]*to_vector(w_r[r,]);
    Umsy_rt[l,r]= 1-lambert_w0(exp(1-log_a_rt[l,r]));
	}
	for(j in 1:J){
 	log_a_jt[l,j]=log_a_j[j]+B_pred_year[l,]*to_vector(w_j[j,]);
    	Umsy_jt[l,j] = 1-lambert_w0(exp(1-log_a_jt[l,j]));
    	Smsy_jt[l,j] = (1-lambert_w0(exp(1-log_a_jt[l,j])))/b[j];
	}
}

for(n in 1:N)log_lik[n]=normal_lpdf(R_S[n]|mu[n],sigma[J_i[n]]);
}