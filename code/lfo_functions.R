source(here('code','dlm-wrapper.R'))

##Functions####
##Functions for LFO-CV
log_sum_exp <- function(x) {
  max_x <- max(x)  
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

#
sr_mod<- function(type=c('static','tv','regime'),ac=FALSE,par=c('a','b','both'),loglik=FALSE){
  if(type=='static'&ac==F){
    if(loglik==FALSE){
    m="data{
      int<lower=1> N;//number of annual samples (time-series length)
      vector[N] R_S; //log(recruits per spawner)
      vector[N] S; //spawners in time T
     }
    parameters {
      real<lower = 0> log_a;// initial productivity (on log scale)
      real<upper = 0> log_b; // rate capacity - fixed in this
    
     //variance components  
      real<lower = 0> sigma_e;
    
    }
    transformed parameters{
    	real b;
    	
    	b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
    }
    model{
      //priors
      log_a ~ normal(0,2.5); //intrinsic productivity - wide prior
      log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
      
      //variance terms
      sigma_e ~ gamma(2,3);
      
      R_S ~ normal(log_a - S*b, sigma_e);
    }
    "
    }
if(loglik==TRUE){
    m ="data{
      int<lower=1> N;//number of annual samples (time-series length)
      vector[N] R_S; //log(recruits per spawner)
      vector[N] S; //spawners in time T
      real y_oos; //log(recruits per spawner)
      real x_oos; //spawners in time T
     }
    parameters {
      real<lower = 0> log_a;// initial productivity (on log scale)
      real<upper = 0> log_b; // rate capacity - fixed in this
    
     //variance components  
      real<lower = 0> sigma_e;
    
    }
    transformed parameters{
    	real b;
    	
    	b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
    }
    model{
      //priors
      log_a ~ normal(0,2.5); //intrinsic productivity - wide prior
      log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
      
      //variance terms
      sigma_e ~ gamma(2,3);
      
      R_S ~ normal(log_a - S*b, sigma_e);
    }
    generated quantities{
     vector[N] log_lik;
     real log_lik_oos;
    	for(n in 1:N)log_lik[n] = normal_lpdf(R_S[n]|log_a - S[n]*b, sigma_e);
    	log_lik_oos = normal_lpdf(y_oos|log_a - x_oos*b, sigma_e);
    }
    "}
  }
if(type=='tv'&par=='a'){
    if(loglik==FALSE){
      m="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
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
  vector[L-1] a_dev; //year-to-year deviations in a
  
}
transformed parameters{
  real b;
  vector[L] log_a; //a in each year (on log scale)
  
  b=exp(log_b);
  
  log_a[1] = log_a0; //initial value
  for(t in 2:L){
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
   
 
  for(n in 1:N) R_S[n] ~ normal(log_a[ii[n]] - S[n]*b, sigma_e); 
  
} "
    }
    if(loglik==TRUE){
      m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters{
  real<lower = 0> log_a0;// initial productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;

  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  
}
transformed parameters{
  real b;
  vector[L] log_a; //a in each year (on log scale)
  
  b=exp(log_b);
  
  log_a[1] = log_a0; //initial value
  for(t in 2:L){
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
   
  for(n in 1:N) R_S[n] ~ normal(log_a[ii[n]] - S[n]*b, sigma_e);
}
  generated quantities{
  real log_a_3b;
  real log_a_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
  
  log_a_3b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]])/3;
  log_a_5b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]]+log_a[ii[N-3]]+log_a[ii[N-4]])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a[ii[N]] - x_oos*b, sigma_e);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b, sigma_e);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b, sigma_e);
 }
    "}
  }  
if(type=='tv'&par=='b'){
  if(loglik==FALSE){
    m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale) - fixed in this
  real<upper = 0> b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] b_dev; //year-to-year deviations in a

}

transformed parameters{
  vector[L] log_b; //b in each year
  vector[L] b; //b in each year
  
  log_b[1] = b0;
  for(t in 2:L){
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a ~ normal(0,2.5); //initial productivity - wide prior
  b0 ~ normal(-12,3); //covariates - reef
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_b ~ gamma(2,3);
   
  b_dev ~ std_normal();
 for(n in 1:N) R_S[n] ~ normal(log_a-b[ii[n]]*S[n], sigma_e);
}
 "
  }
if(loglik==TRUE){
  m="data{
  int<lower=1> N;//number of annual samples
  int<lower=1> L;//number years in the data series(time-series length)
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale) - fixed in this
  real<upper = 0> b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] b_dev; //year-to-year deviations in a

}

transformed parameters{
  vector[L] log_b; //b in each year
  vector[L] b; //b in each year
  
  log_b[1] = b0;
  for(t in 2:L){
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a ~ normal(0,2.5); //initial productivity - wide prior
  b0 ~ normal(-12,3); //covariates - reef
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_b ~ gamma(2,3);
   
  b_dev ~ std_normal();
  
  for(n in 1:N) R_S[n] ~ normal(log_a-S[n]*b[ii[n]], sigma_e);
}
generated quantities{
  real b_3b;
  real b_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
  
  b_3b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]])/3;
  b_5b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]]+b[ii[N-3]]+b[ii[N-4]])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a - x_oos*b[ii[N]], sigma_e);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a - x_oos*b_3b, sigma_e);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a - x_oos*b_5b, sigma_e);
 }
 "}
}
if(type=='tv'&par=='both'){
  if(loglik==FALSE){
    m="data{
  int<lower=1> N;//number of annual samples (time-series length)
  int L; //total years covered by time-series
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real<lower=0> log_a0;// initial productivity (on log scale) - fixed in this
  real<upper=0> log_b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  vector[L-1] b_dev; //year-to-year deviations in a
}

transformed parameters{
  vector[L] log_a; //a in each year (log scale)
  vector[L] log_b; //b in each year (log scale)
  vector[L] b; //b in each year
  
  log_a[1] = log_a0;
  log_b[1] = log_b0;
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a;
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a0 ~ normal(0,2.5); //initial productivity - wide prior
  log_b0 ~ normal(-12,3); //covariates - reef
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_a ~ gamma(2,3);
  sigma_b ~ gamma(2,3);
  
  a_dev ~ std_normal();
  b_dev ~ std_normal();
  
  for(n in 1:N) R_S[n] ~ normal(log_a[ii[n]]-b[ii[n]]*S[n], sigma_e);
}"
  }
if(loglik==TRUE){
  m="data{
  int<lower=1> N;//number of annual samples (time-series length)
  int L; //total years covered by time-series
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
  real<lower=0> log_a0;// initial productivity (on log scale) - fixed in this
  real<upper=0> log_b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  vector[L-1] b_dev; //year-to-year deviations in a
}

transformed parameters{
  vector[L] log_a; //a in each year (log scale)
  vector[L] log_b; //b in each year (log scale)
  vector[L] b; //b in each year
  
  log_a[1] = log_a0;
  log_b[1] = log_b0;
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a;
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a0 ~ normal(0,2.5); //initial productivity - wide prior
  log_b0 ~ normal(-12,3); //covariates - reef
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_a ~ gamma(2,3);
  sigma_b ~ gamma(2,3);
  
  a_dev ~ std_normal();
  b_dev ~ std_normal();
  
  for(n in 1:N) R_S ~ normal(log_a[ii[n]]-b[ii[n]]*S[n], sigma_e);
}
generated quantities{
  real b_3b;
  real b_5b;
  real log_a_3b;
  real log_a_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
   
  log_a_3b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]])/3;
  log_a_5b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]]+log_a[ii[N-3]]+log_a[ii[N-4]])/5;
  
  b_3b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]])/3;
  b_5b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]]+b[ii[N-3]]+b[ii[N-4]])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a[ii[N]] - x_oos*b[ii[N]], sigma_e);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b_3b, sigma_e);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b_5b, sigma_e);
 }
 "}
}
if(type=='regime'&par=='a'){
  if(loglik==FALSE){
    m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
ordered[K] log_a; // max. productivity
real<upper = 0> log_b; // rate capacity - fixed in this
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
real b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator1[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator1[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator1);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities {
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];

{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward

{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator2[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator2[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] | log_a[i] - b*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator2);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 
}
"
  }
if(loglik==TRUE){
  m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
  
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
ordered[K] log_a; // max. productivity
real<upper = 0> log_b; // rate capacity - fixed in this
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
real b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator1[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator1[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator1);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities{
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];

//out of sample log-likelihoods
real log_lik_oos_1b; //OOS log likelihood - non weighted
real log_lik_oos_1bw;//OOS log likelihood - weighted
real log_lik_oos_3b; //OOS log likelihood - non weighted
real log_lik_oos_3bw;//OOS log likelihood - weighted
real log_lik_oos_5b; //OOS log likelihood - non weighted
real log_lik_oos_5bw;//OOS log likelihood - weighted

//slope based on regime in year N
real log_a_1b;
real sigma_1b;
real log_a_3b;
real sigma_3b;
real log_a_5b;
real sigma_5b;

//slope weighted by probability of each regime 
vector[K] log_a_1bw_k;
vector[K] sigma_1bw_k;
vector[K] log_a_3bw_k;
vector[K] sigma_3bw_k;
vector[K] log_a_5bw_k;
vector[K] sigma_5bw_k;

real log_a_1bw;
real sigma_1bw;
real log_a_3bw;
real sigma_3bw;
real log_a_5bw;
real sigma_5bw;

{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward

{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator2[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator2[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] | log_a[i] - b*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator2);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 

log_a_1b = log_a[zstar[N]]; //intercept
sigma_1b = sigma[zstar[N]]; //sigma error
log_a_3b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]])/3; //intercept
sigma_3b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]])/3; //intercept
log_a_5b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]]+log_a[zstar[N-3]]+log_a[zstar[N-4]])/5; //intercept
sigma_5b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]]+sigma[zstar[N-3]]+sigma[zstar[N-4]])/5; //intercept

//slope weighted by probability of each regime 
for(k in 1:K) log_a_1bw_k[k]=gamma[N,k]*log_a[k]; //prob of each regime x productivity for each regime
for(k in 1:K) sigma_1bw_k[k] = gamma[N,k]*sigma[k]; //prob of each regime x residual error for each regime
for(k in 1:K) log_a_3bw_k[k]=(gamma[N,k]*log_a[k]+gamma[N-1,k]*log_a[k]+gamma[N-2,k]*log_a[k])/3; //prob of each regime x productivity for each regime
for(k in 1:K) sigma_3bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k])/3; //prob of each regime x residual error for each regime
for(k in 1:K) log_a_5bw_k[k]=(gamma[N,k]*log_a[k]+gamma[N-1,k]*log_a[k]+gamma[N-2,k]*log_a[k]+gamma[N-3,k]*log_a[k]+gamma[N-4,k]*log_a[k])/5; //prob of each regime x productivity for each regime
for(k in 1:K) sigma_5bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k]+gamma[N-3,k]*sigma[k]+gamma[N-4,k]*sigma[k])/5; //prob of each regime x residual error for each regime

log_a_1bw=sum(log_a_1bw_k); //weighted productivity
sigma_1bw=sum(sigma_1bw_k); //weighted sigma
log_a_3bw=sum(log_a_3bw_k); //weighted productivity
sigma_3bw=sum(sigma_3bw_k); //weighted sigma
log_a_5bw=sum(log_a_5bw_k); //weighted productivity
sigma_5bw=sum(sigma_5bw_k); //weighted sigma

log_lik_oos_1b = normal_lpdf(y_oos|log_a_1b - x_oos*b, sigma_1b);
log_lik_oos_1bw = normal_lpdf(y_oos|log_a_1bw - x_oos*b, sigma_1bw);
log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b, sigma_3b);
log_lik_oos_3bw = normal_lpdf(y_oos|log_a_3bw - x_oos*b, sigma_3bw);
log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b, sigma_5b);
log_lik_oos_5bw = normal_lpdf(y_oos|log_a_5bw - x_oos*b, sigma_5bw);
}
 "}
}
if(type=='regime'&par=='b'){
  if(loglik==FALSE){
    m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
real<lower = 0>  log_a; // max. productivity
ordered[K]<upper = 0> log_b; // rate capacity - fixed in this
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
ordered[K] b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a - b[j]*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,2.5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities{
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];
{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward
{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a - b[i]*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a - b[j]*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a - b[j]*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 
}
"
  }
if(loglik==TRUE){
  m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
  real y_oos; //out of sample (1-year ahead) log(R/S)
  real x_oos; //spawners 1-year ahead
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
real<lower = 0>  log_a; // max. productivity
ordered[K] log_b; // rate capacity - fixed in this
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
ordered[K] b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a - b[j]*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,2.5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities {
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];

//out of sample log-likelihoods
real log_lik_oos_1b; //OOS log likelihood - non weighted
real log_lik_oos_1bw;//OOS log likelihood - weighted
real log_lik_oos_3b; //OOS log likelihood - non weighted
real log_lik_oos_3bw;//OOS log likelihood - weighted
real log_lik_oos_5b; //OOS log likelihood - non weighted
real log_lik_oos_5bw;//OOS log likelihood - weighted

//slope based on regime in year N
real b_1b;
real sigma_1b;
real b_3b;
real sigma_3b;
real b_5b;
real sigma_5b;

//slope weighted by probability of each regime 
vector[K] b_1bw_k;
vector[K] sigma_1bw_k;
vector[K] b_3bw_k;
vector[K] sigma_3bw_k;
vector[K] b_5bw_k;
vector[K] sigma_5bw_k;

real b_1bw;
real sigma_1bw;
real b_3bw;
real sigma_3bw;
real b_5bw;
real sigma_5bw;


{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward
{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a - b[i]*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a - b[j]*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a - b[j]*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 

b_1b = b[zstar[N]]; //slope based on most probable state in sample N
sigma_1b = sigma[zstar[N]]; //sigma error of most probable state in sample N
b_3b = (b[zstar[N]]+b[zstar[N-1]]+b[zstar[N-2]])/3; //intercept
sigma_3b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]])/3; //intercept
b_5b = (b[zstar[N]]+b[zstar[N-1]]+b[zstar[N-2]]+b[zstar[N-3]]+b[zstar[N-4]])/5; //intercept
sigma_5b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]]+sigma[zstar[N-3]]+sigma[zstar[N-4]])/5; //intercept

for(k in 1:K) b_1bw_k[k]=gamma[N,k]*b[k]; //prob of each regime x productivity for each regime
for(k in 1:K) sigma_1bw_k[k] = gamma[N,k]*sigma[k]; //prob of each regime x residual error for each regime
b_3bw_k[k]=(gamma[N,k]*b[k]+gamma[N-1,k]*b[k]+gamma[N-2,k]*b[k])/3; //prob of each regime x productivity for each regime
for(k in 1:K) sigma_3bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k])/3; //prob of each regime x residual error for each regime
for(k in 1:K) b_5bw_k[k]=(gamma[N,k]*b[k]+gamma[N-1,k]*b[k]+gamma[N-2,k]*b[k]+gamma[N-3,k]*b[k]+gamma[N-4,k]*b[k])/5; //prob of each regime x productivity for each regime
for(k in 1:K) sigma_5bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k]+gamma[N-3,k]*sigma[k]+gamma[N-4,k]*sigma[k])/5; //prob of each regime x residual error for each regime

b_1bw=sum(b_1bw_k); //weighted capacity
sigma_1bw=sum(sigma_1bw_k); //weighted sigma
b_3bw=sum(b_3bw_k); //weighted productivity
sigma_3bw=sum(sigma_3bw_k); //weighted sigma
b_5bw=sum(b_5bw_k); //weighted productivity
sigma_5bw=sum(sigma_5bw_k); //weighted sigma

//LL for each prediction
log_lik_oos_1b = normal_lpdf(y_oos|log_a - x_oos*b_1b, sigma_1b);
log_lik_oos_1bw = normal_lpdf(y_oos|log_a - x_oos*b_1bw, sigma_1bw);
log_lik_oos_3b = normal_lpdf(y_oos|log_a - x_oos*b_3b, sigma_3b);
log_lik_oos_3bw = normal_lpdf(y_oos|log_a - x_oos*b_3bw, sigma_3bw);
log_lik_oos_5b = normal_lpdf(y_oos|log_a - x_oos*b_5b, sigma_5b);
log_lik_oos_5bw = normal_lpdf(y_oos|log_a - x_oos*b_5bw, sigma_5bw);
}"}
}
if(type=='regime'&par=='both'){
  if(loglik==FALSE){
    m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
ordered[K]<lower = 0>  log_a; // regime max. productivity
ordered[K]<upper = 0> log_b; // regime rate capacity 
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
ordered[K] b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b[j]*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities {
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];
{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward
{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a[i] - b[i]*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b[j]*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b[j]*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 
}
"
  }
if(loglik==TRUE){
  m="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
  real y_oos; //out of sample (1-year ahead) log(R/S)
  real x_oos; //spawners 1-year ahead
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
ordered[K]  log_a; // regime max. productivity
ordered[K] log_b; // regime rate capacity 
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
ordered[K] b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b[j]*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities {
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];

//out of sample log-likelihoods
real log_lik_oos_1b; //OOS log likelihood - non weighted
real log_lik_oos_1bw;//OOS log likelihood - weighted
real log_lik_oos_3b; //OOS log likelihood - non weighted
real log_lik_oos_3bw;//OOS log likelihood - weighted
real log_lik_oos_5b; //OOS log likelihood - non weighted
real log_lik_oos_5bw;//OOS log likelihood - weighted

//slope based on regime in year N
real log_a_1b;
real log_a_3b;
real log_a_5b;
real b_1b;
real sigma_1b;
real b_3b;
real sigma_3b;
real b_5b;
real sigma_5b;

//slope weighted by probability of each regime 
vector[K] log_a_1bw_k;
vector[K] sigma_1bw_k;
vector[K] log_a_3bw_k;
vector[K] sigma_3bw_k;
vector[K] log_a_5bw_k;
vector[K] sigma_5bw_k;
vector[K] b_1bw_k;
vector[K] b_3bw_k;
vector[K] b_5bw_k;

real b_1bw;
real b_3bw;
real b_5bw;
real log_a_1bw;
real sigma_1bw;
real log_a_3bw;
real sigma_3bw;
real log_a_5bw;
real sigma_5bw;


{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward
{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] |log_a[i] - b[i]*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b[j]*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b[j]*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 

log_a_1b = log_a[zstar[N]]; //intercept
b_1b = b[zstar[N]]; //slope based on most probable state in sample N
sigma_1b = sigma[zstar[N]]; //sigma error
log_a_3b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]])/3; //intercept 3-y back
b_3b = (b[zstar[N]]+b[zstar[N-1]]+b[zstar[N-2]])/3; //intercept
sigma_3b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]])/3; //intercept
log_a_5b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]]+log_a[zstar[N-3]]+log_a[zstar[N-4]])/5; //intercept
b_5b = (b[zstar[N]]+b[zstar[N-1]]+b[zstar[N-2]]+b[zstar[N-3]]+b[zstar[N-4]])/5; 
sigma_5b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]]+sigma[zstar[N-3]]+sigma[zstar[N-4]])/5; //intercept

//slope weighted by probability of each regime
//slope weighted by probability of each regime 
for(k in 1:K){
log_a_1bw_k[k]=gamma[N,k]*log_a[k]; //prob of each regime x productivity for each regime
b_1bw_k[k]=gamma[N,k]*b[k]; //prob of each regime x productivity for each regime
sigma_1bw_k[k] = gamma[N,k]*sigma[k]; //prob of each regime x residual error for each regime

log_a_3bw_k[k]=(gamma[N,k]*log_a[k]+gamma[N-1,k]*log_a[k]+gamma[N-2,k]*log_a[k])/3; //prob of each regime x productivity for each regime
b_3bw_k[k]=(gamma[N,k]*b[k]+gamma[N-1,k]*b[k]+gamma[N-2,k]*b[k])/3; //prob of each regime x productivity for each regime
sigma_3bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k])/3; //prob of each regime x residual error for each regime

log_a_5bw_k[k]=(gamma[N,k]*log_a[k]+gamma[N-1,k]*log_a[k]+gamma[N-2,k]*log_a[k]+gamma[N-3,k]*log_a[k]+gamma[N-4,k]*log_a[k])/5; //prob of each regime x productivity for each regime
b_5bw_k[k]=(gamma[N,k]*b[k]+gamma[N-1,k]*b[k]+gamma[N-2,k]*b[k]+gamma[N-3,k]*b[k]+gamma[N-4,k]*b[k])/5; //prob of each regime x productivity for each regime
sigma_5bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k]+gamma[N-3,k]*sigma[k]+gamma[N-4,k]*sigma[k])/5; //prob of each regime x residual error for each regime
}
log_a_1bw=sum(log_a_1bw_k); //weighted productivity
b_1bw=sum(b_1bw_k); //weighted capacity - 1 year previous
sigma_1bw=sum(sigma_1bw_k); //weighted sigma
log_a_3bw=sum(log_a_3bw_k); //weighted productivity
b_3bw=sum(b_3bw_k); //weighted capacity - 3 year previous average
sigma_3bw=sum(sigma_3bw_k); //weighted sigma
log_a_5bw=sum(log_a_5bw_k); //weighted productivity
b_5bw=sum(b_5bw_k); //weighted capacity - 5 year previous average
sigma_5bw=sum(sigma_5bw_k); //weighted sigma

//LL for each prediction
log_lik_oos_1b = normal_lpdf(y_oos|log_a_1b - x_oos*b_1b, sigma_1b);
log_lik_oos_1bw = normal_lpdf(y_oos|log_a_1bw - x_oos*b_1bw, sigma_1bw);
log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b_3b, sigma_3b);
log_lik_oos_3bw = normal_lpdf(y_oos|log_a_3bw - x_oos*b_1bw, sigma_3bw);
log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b_5b, sigma_5b);
log_lik_oos_5bw = normal_lpdf(y_oos|log_a_5bw - x_oos*b_5bw, sigma_5bw);

}

"}
}
return(m)  
}

#Refit rstan
stan_refit<- function(sm,newdata,oos,regime=FALSE,K=NULL){
  #mod = model file name - eg. 'ricker_linear_oos.stan'
  #newdata = data to train model
  #oosdata = data to predict onto
  #regime = TRUE or FALSE for regime shift models (have different data inputs)
  #K = number of potential regimes (2 or 3)
  
    oosdata=newdata[oos,]
    newdata=newdata[-oos,]
    if(regime==FALSE){
      r = rstan::sampling(sm, 
                      data = list(N=nrow(newdata),
                                  L=max(newdata$broodyear)-min(newdata$broodyear)+1,
                                  ii=newdata$broodyear-min(newdata$broodyear)+1,
                                  R_S =newdata$logR_S,
                                  S=newdata$spawners,
                                  y_oos=oosdata$logR_S,
                                  x_oos=oosdata$spawners),
                      control = list(adapt_delta = 0.99), warmup = 200, chains = 6, iter = 700)
    }
    if(regime==TRUE){
      r = rstan::sampling(sm, 
                      data = list(N=nrow(newdata),
                                  R_S=newdata$logR_S,
                                  S=newdata$spawners,
                                  K=K,
                                  alpha_dirichlet=rep(1,K),
                                  y_oos=oosdata$logR_S,
                                  x_oos=oosdata$spawners), #prior for state transition probabilities (this makes them equal)
                      control = list(adapt_delta = 0.99), warmup = 200, chains = 6, iter = 700)
    }
  
  return(r)
}

#Refit cmdstanr model
stan_mod_refit<- function(mod,newdata,oos){
  if(is.null(oos)==T){
    r=mod$sample(
      data = list(R_S = newdata$y,
                  N=nrow(newdata),
                  TT=as.numeric(factor(newdata$year)),
                  S=c((newdata$spawners))),
      seed = 123, 
      chains = 6, 
      parallel_chains = 6,
      iter_warmup = 200,
      iter_sampling = 500,
      refresh = 200,
      adapt_delta = 0.99,
      max_treedepth = 20 # print update every 500 iters
    )
  }
  
  if(is.null(oos)==F){
    oosdata=newdata[oos,]
    newdata=newdata[-oos,]
    r=mod$sample(
      data = list(R_S = newdata$y,
                  N=nrow(newdata),
                  TT=as.numeric(factor(newdata$year)),
                  S=c((newdata$spawners)),
                  y_oos=oosdata$y,
                  x_oos=oosdata$spawners),
      seed = 123, 
      chains = 6, 
      parallel_chains = 6,
      iter_warmup = 200,
      iter_sampling = 500,
      refresh = 200,
      adapt_delta = 0.99,
      max_treedepth = 20 # print update every 500 iters
    )
  } 
  
  return(r)
}

#Leave-future-out cross-validation function
stan_mod_lfo_cv=function(mod,type=c('static','tv','regime'),df,L){
  #mod = model to fit (model name for cmdstanr)
  #tv = 0 for static model; 1 for time-varying (for calculating elpds)
  #df = full data frame
  #L = starting point for LFO-CV (min. 10)
  loglik_exact <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for static model
  loglik_exact_1b <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for 1-year back estimates of productivity/capacity
  loglik_exact_3b <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for average of last 3-years of productivity/capacity
  loglik_exact_5b <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for average of last 5-years of productivity/capacity
  
  for (i in L:(nrow(df) - 1)) {
    past <- 1:i
    oos <- i + 1
    df_past <- df[past, , drop = FALSE]
    df_oos <- df[c(past, oos), , drop = FALSE]
    fit_past<- stan_mod_refit(mod=mod,newdata=df_oos,oos=i+1)
    if(tv==0){
      loglik_exact[,i+1]<- fit_past$draws(variables=c('log_lik_oos'),format='draws_matrix')
    }else{
      loglik_exact_1b[, i + 1] <- fit_past$draws(variables=c('log_lik_oos_1b'),format='draws_matrix')
      loglik_exact_3b[, i + 1] <- fit_past$draws(variables=c('log_lik_oos_3b'),format='draws_matrix')
      loglik_exact_5b[, i + 1] <- fit_past$draws(variables=c('log_lik_oos_5b'),format='draws_matrix')
    }
  }
  if(tv==0){
    exact_elpds<- apply(loglik_exact, 2, log_mean_exp); exact_elpds=exact_elpds[-(1:L)]
    return(exact_elpds)
  }else{
    exact_elpds_1b <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b=exact_elpds_5b[-(1:L)]
    
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }
}


stan_lfo_cv=function(mod,type=c('static','tv','regime'),df,L=10,K=NULL){
  #mod = model to fit (model name for cmdstanr)
  #tv = 0 for static model; 1 for time-varying (for calculating elpds)
  #df = full data frame
  #L = starting point for LFO-CV (default 10)
  sm <- stan_model(model_code = mod)
  
  loglik_exact <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for static model
  loglik_exact_1b <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for 1-year back estimates of productivity/capacity
  loglik_exact_3b <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for average of last 3-years of productivity/capacity
  loglik_exact_5b <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for average of last 5-years of productivity/capacity
  if(type=='regime'){
    loglik_exact_1bw <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    loglik_exact_3bw <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for average of last 3-years of productivity/capacity
    loglik_exact_5bw <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for average of last 5-years of productivity/capacity
  }
  for (i in L:(nrow(df) - 1)){
    past <- 1:i
    oos <- i + 1
    df_past <- df[past, , drop = FALSE]
    df_oos <- df[c(past, oos), , drop = FALSE]
    if(type=='static'){
      fit_past<- stan_refit(sm=sm,newdata=df_oos,oos=i+1)
      ll=extract(fit_past,pars=c('log_lik_oos'))
      loglik_exact[,i+1]<- ll$log_lik_oos
      
    }
    if(type=='tv'){
      fit_past<- stan_refit(sm=sm,newdata=df_oos,oos=i+1)
      ll=extract(fit_past,pars=c('log_lik_oos_1b','log_lik_oos_3b','log_lik_oos_5b'))
      loglik_exact_1b[, i + 1] <-ll$log_lik_oos_1b
      loglik_exact_3b[, i + 1] <-ll$log_lik_oos_3b
      loglik_exact_5b[, i + 1] <-ll$log_lik_oos_5b
    }
    if(type=='regime'){
      fit_past<- stan_refit(sm=sm,newdata=df_oos,oos=i+1,regime=TRUE,K=K)
      ll=extract(fit_past,pars=c('log_lik_oos_1b','log_lik_oos_3b','log_lik_oos_5b','log_lik_oos_1bw','log_lik_oos_3bw','log_lik_oos_5bw'))
      loglik_exact_1b[, i + 1] <- ll$log_lik_oos_1b
      loglik_exact_3b[, i + 1] <- ll$log_lik_oos_3b
      loglik_exact_5b[, i + 1] <- ll$log_lik_oos_5b
      loglik_exact_1bw[, i + 1] <- ll$log_lik_oos_1bw
      loglik_exact_3bw[, i + 1] <- ll$log_lik_oos_3bw
      loglik_exact_5bw[, i + 1] <- ll$log_lik_oos_5bw
    }
  }
  
  if(type=='static'){
    exact_elpds<- apply(loglik_exact, 2, log_mean_exp); exact_elpds=exact_elpds[-(1:L)]
    return(exact_elpds)
  }
  if(type=='tv'){
    exact_elpds_1b <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b=exact_elpds_5b[-(1:L)]
    
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }
  if(type=='regime'){
    exact_elpds_1b <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b=exact_elpds_5b[-(1:L)]
    exact_elpds_1bw <- apply(loglik_exact_1bw, 2, log_mean_exp); exact_elpds_1bw=exact_elpds_1b[-(1:L)]
    exact_elpds_3bw <- apply(loglik_exact_3bw, 2, log_mean_exp); exact_elpds_3bw=exact_elpds_3b[-(1:L)]
    exact_elpds_5bw <- apply(loglik_exact_5bw, 2, log_mean_exp); exact_elpds_5bw=exact_elpds_5b[-(1:L)]
    
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b,exact_elpds_1bw,exact_elpds_3bw,exact_elpds_5bw))
  }
}

#Returns the pointwise log-likelihoods for L+1:N
dlm_mod_refit<- function(mod,newdat){
  require(dlm)
  source(here('code','dlm-wrapper.R'))
  if(mod==2){
    avary<-fitDLM(newdat, alpha_vary = TRUE, beta_vary = FALSE)
    return(avary)
  }
  if(mod==3){
    bvary<-fitDLM(newdat, alpha_vary = FALSE, beta_vary = TRUE)
    return(bvary)
  }
  if(mod==4){
    abvary<-fitDLM(newdat, alpha_vary = TRUE, beta_vary = TRUE)
    return(abvary)
  }
}

dlm_mod_lfo_cv=function(mod,df,L){
  #mod = model to fit (model name for cmdstanr)
  #tv = 0 for static model; 1 for time-varying (for calculating elpds)
  #df = full data frame
  #L = starting point for LFO-CV (min. 10)
  exact_elpds_1b<- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
  exact_elpds_3b<- numeric(nrow(df)) #loglik for average of last 3-years of productivity/capacity
  exact_elpds_5b<- numeric(nrow(df))#loglik for average of last 5-years of productivity/capacity
  for (i in L:(nrow(df) - 1)) {
    past <- 1:i
    oos <- i + 1
    df_past <- df[past, , drop = FALSE]
    df_oos <- df[c(past, oos), , drop = FALSE]
    fit_past<- dlm_mod_refit(mod=mod,newdat=df_past)
    rs_pred_1b=fit_past$results$alpha[i]+fit_past$results$beta[i]*df_oos$spwn[i + 1]
    rs_pred_3b=mean(fit_past$results$alpha[(i-2):i])+mean(fit_past$results$beta[(i-2):i])*df_oos$spwn[i + 1]
    rs_pred_5b=mean(fit_past$results$alpha[(i-4):i])+mean(fit_past$results$beta[(i-4):i])*df_oos$spwn[i + 1]
    
    exact_elpds_1b[i+1] <- log(dnorm(log(df_oos$rec[i + 1]/df_oos$spwn[i + 1]),mean=rs_pred_1b,sd=fit_past$sd.est[1]))
    exact_elpds_3b[i+1] <- log(dnorm(log(df_oos$rec[i + 1]/df_oos$spwn[i + 1]),mean=rs_pred_3b,sd=fit_past$sd.est[1]))
    exact_elpds_5b[i+1] <- log(dnorm(log(df_oos$rec[i + 1]/df_oos$spwn[i + 1]),mean=rs_pred_5b,sd=fit_past$sd.est[1]))
  }
  exact_elpds_1b=exact_elpds_1b[-(1:L)]
  exact_elpds_3b=exact_elpds_3b[-(1:L)]
  exact_elpds_5b=exact_elpds_5b[-(1:L)]
  
  return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
}

#
lm_mod_refit<- function(y,x,ac=FALSE){
  require(nlme)
  if(ac==FALSE){
    lm_rf<- lm(y~x)
    return(lm_rf)
  }
  if(ac==TRUE){
    lm_ac_rf<- nlme::gls(y~x,correlation=corAR1(form = ~ 1))
    return(lm_ac_rf)
  }
}

lm_mod_lfo_cv=function(df,ac=F,L){
  #df = full data frame
  #ac = autocorrelation, if ac=T then implement AR-1
  #L = starting point for LFO-CV (min. 10)
  if(ac==F){
    exact_elpds_1b <- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    for (i in L:(nrow(df) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      fit_past_lm<- lm_mod_refit(y=df_past$y,x=df_past$spawners,ac=FALSE)
      rs_pred_1b=fit_past_lm$coefficients[1]+fit_past_lm$coefficients[2]*df_oos$spawners[i + 1]
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_1b,sd=sigma(fit_past_lm)))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    return(exact_elpds_1b)
  }
  if(ac==T){
    exact_elpds_ac_1b <- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_ac_3b <- numeric(nrow(df)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_ac_5b <- numeric(nrow(df))#loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(df) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      
      fit_past_ac_lm<- lm_mod_refit(y=df_past$y,x=df_past$spawners,ac=TRUE)
      #AR-1 model - extract phi and residuals 1-year back, 3-years back, and 5-years back
      phi=intervals(fit_past_ac_lm)$corStruct[2]
      ac_resids_1b=residuals(fit_past_ac_lm)[i]
      ac_resids_3b=mean(residuals(fit_past_ac_lm)[(i-2):i])
      ac_resids_5b=mean(residuals(fit_past_ac_lm)[(i-4):i])
      
      rs_pred_ac_1b=fit_past_ac_lm$coefficients[1]+fit_past_ac_lm$coefficients[2]*df_oos$spawners[i + 1]+phi*ac_resids_1b
      rs_pred_ac_3b=fit_past_ac_lm$coefficients[1]+fit_past_ac_lm$coefficients[2]*df_oos$spawners[i + 1]+phi*ac_resids_3b
      rs_pred_ac_5b=fit_past_ac_lm$coefficients[1]+fit_past_ac_lm$coefficients[2]*df_oos$spawners[i + 1]+phi*ac_resids_5b
      
      exact_elpds_ac_1b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_ac_1b,sd=fit_past_ac_lm$sigma))
      exact_elpds_ac_3b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_ac_3b,sd=fit_past_ac_lm$sigma))
      exact_elpds_ac_5b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_ac_5b,sd=fit_past_ac_lm$sigma))
    }
    exact_elpds_ac_1b=exact_elpds_ac_1b[-(1:L)]
    exact_elpds_ac_3b=exact_elpds_ac_3b[-(1:L)]
    exact_elpds_ac_5b=exact_elpds_ac_5b[-(1:L)]
    return(list(exact_elpds_ac_1b,exact_elpds_ac_3b,exact_elpds_ac_5b))
  }
}

tmb_ricker_refit<- function(df){
  srm <- lm(df$y~ df$spawners)
  SRdata<-list(obs_logRS=df$y,obs_S=df$spawners)
  
  parameters_simple<- list(
    alpha=srm$coefficients[1],
    logbeta = log(ifelse(-srm$coefficients[2]<0,1e-08,-srm$coefficients[2])),
    logsigobs=log(.4)
  )
  obj<- TMB::MakeADFun(SRdata,parameters_simple,DLL="Ricker_simple")
  opt_simple <- nlminb(obj$par,obj$fn,obj$gr)
  if(opt_simple$convergence==0){
    return(list(obj$report(),opt_simple))  
  }
}

tmb_mod_refit<- function(df,tv.par=c('alpha','beta','both')){
  srm <- lm(df$y~df$spawners)
  SRdata<-list(obs_logRS=df$y,obs_S=df$spawners)
  if(tv.par=='alpha'){
    parameters<- list(
      alphao=srm$coefficients[1],
      logbeta = log(ifelse(-srm$coefficients[2]<0,1e-08,-srm$coefficients[2])),
      #rho=.5,
      #logvarphi=0,
      logsigobs=log(.4),
      logsiga=log(.4),
      alpha=rep(srm$coefficients[1],length(df$y))
    )
    
    obja<- MakeADFun(SRdata,parameters,DLL="Ricker_tva_Smax",random="alpha")#,lower = -Inf, upper = Inf)
    obja$fn()
    obja$gr()
    opta=nlminb(obja$par,obja$fn,obja$gr)
    if(opta$convergence==0){
      return(list(obja$report(),opta))  
    }
  }
  
  if(tv.par=='beta'){
    parameters<- list(
      logbetao = log(ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[[2]])),
      alpha=srm$coefficients[[1]],
      logsigobs=log(.4),
      logsigb=log(.4),
      #rho=.4,
      #logvarphi= 0.5,
      logbeta=log(rep(-srm$coefficients[2],length(df$spawners)))
    )
    objb <- MakeADFun(SRdata,parametersb,DLL="Ricker_tvb",random="logbeta",lower=c(-10,-20,-6,-6),
                      upper=c(10,0,2,2))
    objb$fn()
    objb$gr()
    optb=nlminb(objb$par,objb$fn,objb$gr)
    if(optb$convergence==0){
      return(list(objb$report(),optb))  
    }
  }
  
  if(tv.par=='both'){
    parametersab<- list(
      logbetao = log(ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[[2]])),
      alphao=srm$coefficients[[1]],
      logsigobs=log(.5),
      logsiga=log(.1),
      logsigb=log(.1),
      logbeta=rep(log(ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[[2]])),length(df$spawners)),
      alpha=rep(srm$coefficients[[1]],length(df$spawners))
    )
    
    objab <- MakeADFun(SRdata,parametersab,DLL="Ricker_tva_tvb",random=c("logbeta","alpha"))
    objab$fn()
    objab$gr()
    optab=nlminb(objab$par,objab$fn,objab$gr)
    if(optab$convergence==0){
      return(list(objab$report(),optab))  
    }
  }
}

tmb_mod_lfo_cv=function(df,tv.par=c('static','alpha','beta','both'),L){
  #df = full data frame
  #ac = autocorrelation, if ac=T then implement AR-1
  #L = starting point for LFO-CV (min. 10)
  if(tv.par=='static'){
    exact_elpds_1b <- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    for (i in L:(nrow(df) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      
      fit_past_tmb<- tmb_ricker_refit(df=df_past)
      rs_pred_1b=fit_past_tmb[[1]]$alpha-fit_past_tmb[[1]]$beta*df_oos$spawners[i + 1]
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_1b,sd=exp(fit_past_tmb[[2]]$par[3])))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    return(exact_elpds_1b)
  }
  if(tv.par=='alpha'){
    exact_elpds_1b <- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(df)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(df))#loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(df) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      
      fit_past_tv_a_tmb<- tmb_mod_refit(df=df_past,tv.par='alpha')
      
      rs_pred_1b=fit_past_tv_a_tmb[[1]]$alpha[i]-fit_past_tv_a_tmb[[1]]$beta*df_oos$spawners[i + 1]
      rs_pred_3b=mean(fit_past_tv_a_tmb[[1]]$alpha[(i-2):i])-fit_past_tv_a_tmb[[1]]$beta*df_oos$spawners[i + 1]
      rs_pred_5b=mean(fit_past_tv_a_tmb[[1]]$alpha[(i-4):i])-fit_past_tv_a_tmb[[1]]$beta*df_oos$spawners[i + 1]
      
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_1b,sd=exp(fit_past_tv_a_tmb[[2]]$par[3])))
      exact_elpds_3b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_3b,sd=exp(fit_past_tv_a_tmb[[2]]$par[3])))
      exact_elpds_5b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_5b,sd=exp(fit_past_tv_a_tmb[[2]]$par[3])))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }
  if(tv.par=='beta'){
    exact_elpds_1b <- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(df)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(df))#loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(df) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      
      fit_past_tv_b_tmb<- tmb_mod_refit(df=df_past,tv.par='beta')
      
      rs_pred_1b=fit_past_tv_b_tmb[[1]]$alpha-fit_past_tv_b_tmb[[1]]$beta[i]*df_oos$spawners[i + 1]
      rs_pred_3b=fit_past_tv_b_tmb[[1]]$alpha-mean(fit_past_tv_b_tmb[[1]]$beta[(i-2):i])*df_oos$spawners[i + 1]
      rs_pred_5b=fit_past_tv_b_tmb[[1]]$alpha-mean(fit_past_tv_b_tmb[[1]]$beta[(i-4):i])*df_oos$spawners[i + 1]
      
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_1b,sd=exp(fit_past_tv_b_tmb[[2]]$par[3])))
      exact_elpds_3b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_3b,sd=exp(fit_past_tv_b_tmb[[2]]$par[3])))
      exact_elpds_5b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_5b,sd=exp(fit_past_tv_b_tmb[[2]]$par[3])))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }
  if(tv.par=='both'){
    exact_elpds_1b <- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(df)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(df))#loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(df) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      
      fit_past_tv_ab_tmb<- tmb_mod_refit(df=df_past,tv.par='both')
      
      rs_pred_1b=fit_past_tv_ab_tmb[[1]]$alpha[i]-fit_past_tv_ab_tmb[[1]]$beta[i]*df_oos$spawners[i + 1]
      rs_pred_3b=mean(fit_past_tv_ab_tmb[[1]]$alpha[(i-2):i])-mean(fit_past_tv_ab_tmb[[1]]$beta[(i-2):i])*df_oos$spawners[i + 1]
      rs_pred_5b=mean(fit_past_tv_ab_tmb[[1]]$alpha[(i-4):i])-mean(fit_past_tv_ab_tmb[[1]]$beta[(i-4):i])*df_oos$spawners[i + 1]
      
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_1b,sd=exp(fit_past_tv_ab_tmb[[2]]$par[3])))
      exact_elpds_3b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_3b,sd=exp(fit_past_tv_ab_tmb[[2]]$par[3])))
      exact_elpds_5b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_5b,sd=exp(fit_past_tv_ab_tmb[[2]]$par[3])))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }
}



##Stan Models####
#static linearized Ricker stock-recruitment
#m1: Static####
m1="data{
  int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;

}
transformed parameters{
	real b;
	
	b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
}
model{
  //priors
  log_a ~ normal(0,2.5); //intrinsic productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  //variance terms
  sigma_e ~ gamma(2,3);
  
  R_S ~ normal(log_a - S*b, sigma_e);
}
"
#m1oos: Static (LFO)####
m1oos="data{
  int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
  real<lower = 0> log_a;// initial productivity (on log scale)
  real<upper = 0> log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;

}
transformed parameters{
	real b;
	
	b = exp(log_b); //prevents b (density dependence) from being negative (ie. positive)
}
model{
  //priors
  log_a ~ normal(0,2.5); //intrinsic productivity - wide prior
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  //variance terms
  sigma_e ~ gamma(2,3);
  
  R_S ~ normal(log_a - S*b, sigma_e);
}
generated quantities{
 vector[N] log_lik;
 real log_lik_oos;
	for(n in 1:N)log_lik[n] = normal_lpdf(R_S[n]|log_a - S[n]*b, sigma_e);
	log_lik_oos = normal_lpdf(y_oos|log_a - x_oos*b, sigma_e);
}
"
#m5oos: TV a&b (LFO)####
m5oos="data{
  int<lower=1> N;//number of annual samples (time-series length)
  int L; //total years covered by time-series
  int ii[N];//index of years with data
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
  real<lower=0> log_a0;// initial productivity (on log scale) - fixed in this
  real<upper=0> log_b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[L-1] a_dev; //year-to-year deviations in a
  vector[L-1] b_dev; //year-to-year deviations in a
}

transformed parameters{
  vector[L] log_a; //a in each year (log scale)
  vector[L] log_b; //b in each year (log scale)
  vector[L] b; //b in each year
  
  log_a[1] = log_a0;
  log_b[1] = log_b0;
  for(t in 2:L){
    log_a[t] = log_a[t-1] + a_dev[t-1]*sigma_a;
    log_b[t] = log_b[t-1] + b_dev[t-1]*sigma_b;
  } 
  b=exp(log_b);
}  

model{
  //priors
  log_a0 ~ normal(0,2.5); //initial productivity - wide prior
  log_b0 ~ normal(-12,3); //covariates - reef
  
  //variance terms
  sigma_e ~ gamma(2,3);
  sigma_a ~ gamma(2,3);
  sigma_b ~ gamma(2,3);
  
  a_dev ~ std_normal();
  b_dev ~ std_normal();
  
  for(n in 1:N) R_S ~ normal(log_a[ii[n]]-b[ii[n]]*S[n], sigma_e);
}
generated quantities{
  real b_3b;
  real b_5b;
  real log_a_3b;
  real log_a_5b;
  real log_lik_oos_1b;
  real log_lik_oos_3b;
  real log_lik_oos_5b;
   
  log_a_3b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]])/3;
  log_a_5b = (log_a[ii[N]]+log_a[ii[N-1]]+log_a[ii[N-2]]+log_a[ii[N-3]]+log_a[ii[N-4]])/5;
  
  b_3b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]])/3;
  b_5b = (b[ii[N]]+b[ii[N-1]]+b[ii[N-2]]+b[ii[N-3]]+b[ii[N-4]])/5;
  
  log_lik_oos_1b = normal_lpdf(y_oos|log_a[ii[N]] - x_oos*b[ii[N]], sigma_e);
  log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b_3b, sigma_e);
  log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b_5b, sigma_e);
 }
"

#mod 6oos: Regime shift a####
m6oos="functions {
vector normalize(vector x) {
return x / sum(x);
}
}
data {
 int<lower=1> N;//number of annual samples (time-series length)
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
  int<lower=1> K; //number of hidden regime states
  vector[K] alpha_dirichlet; //prior inputs for dirichlet 
  
  real y_oos; //log(recruits per spawner)
  real x_oos; //spawners in time T
 }
parameters {
// Discrete state model
simplex[K] pi1; // initial state probabilities
simplex[K] A[K]; // transition probabilities

// A[i][j] = p(z_t = j | z_{t-1} = i)
// Continuous observation model
ordered[K] log_a; // max. productivity
real<upper = 0> log_b; // rate capacity - fixed in this
real<lower=0> sigma[K]; // observation standard deviations
}

transformed parameters {
vector[K] logalpha[N];
real b; //

b=exp(log_b);
 
{ // Forward algorithm log p(z_t = j | y_{1:t})
real accumulator1[K];

logalpha[1] = log(pi1) + normal_lpdf(R_S[1] |log_a - b*S[1], sigma);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
for (i in 1:K) { // i = previous (t-1)
// Murphy (2012) p. 609 eq. 17.48
// belief state + transition prob + local evidence at t
accumulator1[i] = logalpha[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] |log_a[j] - b*S[t], sigma[j]);
}
logalpha[t, j] = log_sum_exp(accumulator1);
}
}
} // Forward
}
model{
sigma ~ gamma(2,3);
log_a ~ normal(0,5);
log_b ~ normal(-12,3);

pi1 ~ dirichlet(rep_vector(1, K));

for(k in 1:K){
A[k,] ~ dirichlet(alpha_dirichlet);
}

target += log_sum_exp(logalpha[N]);
}
generated quantities{
int<lower=1, upper=K> zstar[N];
real logp_zstar;
vector[K] alpha[N];
vector[K] logbeta[N];
vector[K] loggamma[N];
vector[K] beta[N];
vector[K] gamma[N];

//out of sample log-likelihoods
real log_lik_oos_1b; //OOS log likelihood - non weighted
real log_lik_oos_1bw;//OOS log likelihood - weighted
real log_lik_oos_3b; //OOS log likelihood - non weighted
real log_lik_oos_3bw;//OOS log likelihood - weighted
real log_lik_oos_5b; //OOS log likelihood - non weighted
real log_lik_oos_5bw;//OOS log likelihood - weighted

//slope based on regime in year N
real log_a_1b;
real sigma_1b;
real log_a_3b;
real sigma_3b;
real log_a_5b;
real sigma_5b;

//slope weighted by probability of each regime 
vector[K] log_a_1bw_k;
vector[K] sigma_1bw_k;
vector[K] log_a_3bw_k;
vector[K] sigma_3bw_k;
vector[K] log_a_5bw_k;
vector[K] sigma_5bw_k;

real log_a_1bw;
real sigma_1bw;
real log_a_3bw;
real sigma_3bw;
real log_a_5bw;
real sigma_5bw;

{ // Forward algortihm
for (t in 1:N)
alpha[t] = softmax(logalpha[t]);
} // Forward

{ // Backward algorithm log p(y_{t+1:T} | z_t = j)
real accumulator2[K];
for (j in 1:K)
logbeta[N, j] = 1;
for (tforward in 0:(N-2)) {
int t;
t = N - tforward;
for (j in 1:K) { // j = previous (t-1)
for (i in 1:K) { // i = next (t)
// Murphy (2012) Eq. 17.58
// backwards t + transition prob + local evidence at t
accumulator2[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(R_S[t] | log_a[i] - b*S[t], sigma[i]);
}
logbeta[t-1, j] = log_sum_exp(accumulator2);
}
}
for (t in 1:N)
beta[t] = softmax(logbeta[t]);
} // Backward
{ // forward-backward algorithm log p(z_t = j | y_{1:N})
for(t in 1:N) {
loggamma[t] = alpha[t] .* beta[t];
}
for(t in 1:N)
gamma[t] = normalize(loggamma[t]);
} // forward-backward

{ // Viterbi algorithm
int bpointer[N, K]; // backpointer to the most likely previous state on the most probable path
real delta[N, K]; // max prob for the sequence up to t
// that ends with an emission from state k
for (j in 1:K)
delta[1, K] = normal_lpdf(R_S[1] | log_a[j] - b*S[1], sigma[j]);
for (t in 2:N) {
for (j in 1:K) { // j = current (t)
delta[t, j] = negative_infinity();
for (i in 1:K) { // i = previous (t-1)
real logp;
logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(R_S[t] | log_a[j] - b*S[t], sigma[j]);
if (logp > delta[t, j]) {
bpointer[t, j] = i;
delta[t, j] = logp;
}
}
}
}
logp_zstar = max(delta[N]);
for (j in 1:K)
if (delta[N, j] == logp_zstar)
zstar[N] = j;
for (t in 1:(N - 1)) {
zstar[N - t] = bpointer[N - t + 1, zstar[N - t + 1]];
}
} 

log_a_1b = log_a[zstar[N]]; //intercept
sigma_1b = sigma[zstar[N]]; //sigma error
log_a_3b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]])/3; //intercept
sigma_3b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]])/3; //intercept
log_a_5b = (log_a[zstar[N]]+log_a[zstar[N-1]]+log_a[zstar[N-2]]+log_a[zstar[N-3]]+log_a[zstar[N-4]])/5; //intercept
sigma_5b = (sigma[zstar[N]]+sigma[zstar[N-1]]+sigma[zstar[N-2]]+sigma[zstar[N-3]]+sigma[zstar[N-4]])/5; //intercept

//slope weighted by probability of each regime 
for(k in 1:K) log_a_1bw_k[k]=gamma[N,k]*log_a[k]; //prob of each regime x productivity for each regime
for(k in 1:K) sigma_1bw_k[k] = gamma[N,k]*sigma[k]; //prob of each regime x residual error for each regime
for(k in 1:K) log_a_3bw_k[k]=(gamma[N,k]*log_a[k]+gamma[N-1,k]*log_a[k]+gamma[N-2,k]*log_a[k])/3; //prob of each regime x productivity for each regime
for(k in 1:K) sigma_3bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k])/3; //prob of each regime x residual error for each regime
for(k in 1:K) log_a_5bw_k[k]=(gamma[N,k]*log_a[k]+gamma[N-1,k]*log_a[k]+gamma[N-2,k]*log_a[k]+gamma[N-3,k]*log_a[k]+gamma[N-4,k]*log_a[k])/5; //prob of each regime x productivity for each regime
for(k in 1:K) sigma_5bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k]+gamma[N-3,k]*sigma[k]+gamma[N-4,k]*sigma[k])/5; //prob of each regime x residual error for each regime

log_a_1bw=sum(log_a_1bw_k); //weighted productivity
sigma_1bw=sum(sigma_1bw_k); //weighted sigma
log_a_3bw=sum(log_a_3bw_k); //weighted productivity
sigma_3bw=sum(sigma_3bw_k); //weighted sigma
log_a_5bw=sum(log_a_5bw_k); //weighted productivity
sigma_5bw=sum(sigma_5bw_k); //weighted sigma

log_lik_oos_1b = normal_lpdf(y_oos|log_a_1b - x_oos*b, sigma_1b);
log_lik_oos_1bw = normal_lpdf(y_oos|log_a_1bw - x_oos*b, sigma_1bw);
log_lik_oos_3b = normal_lpdf(y_oos|log_a_3b - x_oos*b, sigma_3b);
log_lik_oos_3bw = normal_lpdf(y_oos|log_a_3bw - x_oos*b, sigma_3bw);
log_lik_oos_5b = normal_lpdf(y_oos|log_a_5b - x_oos*b, sigma_5b);
log_lik_oos_5bw = normal_lpdf(y_oos|log_a_5bw - x_oos*b, sigma_5bw);
}
"
