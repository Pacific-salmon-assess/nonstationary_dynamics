functions {
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
for(k in 1:K){
log_a_1bw_k[k]=gamma[N,k]*log_a[k]; //prob of each regime x productivity for each regime

sigma_1bw_k[k] = gamma[N,k]*sigma[k]; //prob of each regime x residual error for each regime

log_a_3bw_k[k]=(gamma[N,k]*log_a[k]+gamma[N-1,k]*log_a[k]+gamma[N-2,k]*log_a[k])/3; //prob of each regime x productivity for each regime

sigma_3bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k])/3; //prob of each regime x residual error for each regime

log_a_5bw_k[k]=(gamma[N,k]*log_a[k]+gamma[N-1,k]*log_a[k]+gamma[N-2,k]*log_a[k]+gamma[N-3,k]*log_a[k]+gamma[N-4,k]*log_a[k])/5; //prob of each regime x productivity for each regime

sigma_5bw_k[k] = (gamma[N,k]*sigma[k]+gamma[N-1,k]*sigma[k]+gamma[N-2,k]*sigma[k]+gamma[N-3,k]*sigma[k]+gamma[N-4,k]*sigma[k])/5; //prob of each regime x residual error for each regime
}
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
