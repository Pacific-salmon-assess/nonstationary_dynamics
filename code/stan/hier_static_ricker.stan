data{
  int<lower=1> N;//number of annual samples 
  int<lower=1> J;// number of time series
  int<lower=1> R;//number of ocean regions
  int<lower=1> SR; //mumber of subregions
  int L; //years covered by time-series
  array[N] int ii; //index of brood years
  array[N] int ij;//index of stocks
  array[SR] int sr_or;//index of ocean-region by sub-region
  array[J] int j_sr;//index of sub-region by stock  
  array[J] int  start_y; // ragged start point for observations (N)
  array[J] int end_y;   // ragged end points for observations (N)
  array[N] real R_S; //log(recruits per spawner)
  array[N] real S; //vector of spawners
  array[J] real pSmax_mean; //priors on smax - based on observed spawner abundance
  array[J] real pSmax_sig;
}
parameters{
 //productivity
  real<lower=0> log_a0;// global productivity mean
  vector[R] z_a_or; //ocean region-level (WC/GOA/BS) productivity deviation
  vector[SR] z_a_sr;//sub-region productivity deviation
  vector[J] z_a_j;//stock-specific productivity deviation
  real<lower=0> sd_a_or; //variance in ocean-region productivity
  real<lower=0> sd_a_sr; //variance in sub-region productivity by ocean region
  real<lower=0> sd_a; //variance in stock productivity by sub-region
 
//capacity
 vector<lower = 0>[J] Smax; //spawners at max. recruitment - independent for each stock
 
 //variance components
 vector<lower = 0>[J] sigma; //mean variance among all stocks
 
//autocorrelation
 vector<lower=-1,upper=1>[J] rho; //non-hierarchical - independently estimated for each stock
 
 }
transformed parameters{
  //expectation
  vector[N] mu;
  vector[N] epsilon;
  vector[N] sigma_t;
  
  //per capita density dependence
  vector[J] b;
  
  //productivity estimates
  vector[R] log_a_or; //ocean region specific productivity
  vector[SR] log_a_sr; //sub-region specific productivity
  vector[J] log_a; //stock-specific productivity
  
  //sigma estimates
  vector<lower=0>[J] sigmaAR; //sigma - adjusted for autocorrelation
  
   //productivity terms
   log_a_or = log_a0 + z_a_or*sd_a_or; //ocean-region-specific productivity
   log_a_sr = log_a_or[sr_or] + z_a_sr*sd_a_sr; //subregion-specific productivity
   log_a = log_a_sr[j_sr] + z_a_j*sd_a; //stock-specific productivity

 //variance terms
   
 
    //expectation, autocorrelation, residuals
   for(j in 1:J){ //for every stock
        b[j]=1/Smax[j];
		mu[start_y[j]]=log_a[j]- S[start_y[j]]*b[j]; //initial expectation
		epsilon[start_y[j]]=R_S[start_y[j]]-mu[start_y[j]];//initial residual
		sigma_t[start_y[j]]=sigma[j];
		
		for(t in (start_y[j]+1):(end_y[j])){
		mu[t]=log_a[j] - S[t]*b[j] + (rho[j]^(ii[t]-ii[t-1]))*epsilon[t-1]; //autocorrelated expectation
		epsilon[t]=R_S[t]-(mu[t]-(rho[j]^(ii[t]-ii[t-1]))*epsilon[t-1]);
		//residual
		sigma_t[t] = sigma[j]*sqrt(1-rho[j]^(2*(ii[t]-ii[t-1])));
		//sigma - adjusted for autocorrelation and gaps in time-series (hence why observation-specific sigma);
		}
	sigmaAR[j] = mean(sigma_t[start_y[j]:end_y[j]]); //average AR(1) adjusted sigma for stock j
   }		
}  
model{
  //priors
  
  //prod
  log_a0 ~ normal(1.5,2.5); //initial productivity - wide prior
  z_a_or ~ std_normal(); //N(0,1)
  z_a_sr ~ std_normal(); //N(0,1)
  z_a_j ~ std_normal(); //N(0,1)
  sd_a_or ~ normal(0,1);
  sd_a_sr ~ normal(0,1);
  sd_a ~ normal(0,1);
  
  //capacity for each stock
  for(j in 1:J) Smax[j] ~ normal(pSmax_mean[j],pSmax_sig[j]);
  
  //variance terms
  sigma ~ gamma(2,1); //gamma prior on variance (lower limit of zero)
  rho~ uniform(-1,1);
  
  //likelihood statement
  R_S ~ normal(mu, sigma_t);   
}
 generated quantities{
     vector[N] log_lik;
     vector[J] Umsy;
     vector[J] Smsy;
     
    for(n in 1:N) log_lik[n] = normal_lpdf(R_S[n]|mu[n], sigmaAR[ij[n]]); 
   
	for(j in 1:J){
    Umsy[j] = 1-lambert_w0(exp(1-(log_a[j])));
    Smsy[j] = (1-lambert_w0(exp(1-(log_a[j]))))/b[j];
    }
}
