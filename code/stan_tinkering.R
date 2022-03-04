#Stan models
Ricker_linear<- "data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real a;// initial productivity (on log scale)
  real b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;

}
model{
  //priors
  a ~ normal(0,5); //initial productivity - wide prior
  b ~ normal(0,5); //initial productivity - wide prior
  
  //variance terms
  sigma_e ~ inv_gamma(2, 1);

  R_S ~ normal(a - b*S, sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(R_S|a - S[t]*b, sigma_e);
  } 
"

Ricker_linear_ac<- "data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters{
  real a;// initial productivity (on log scale)
  real b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = -1, upper = 1> phi;

}
transformed parameters{
vector[N] mu;
vector[N] epsilon; //residuals
mu = a-b*S;

epsilon[1] = R_S[1] - mu[1];
  for(t in 2:N){
    epsilon[t] =(R_S[t] - mu[t]);
    mu[t] = mu[t] + (phi*epsilon[t-1]);
  }

}
model{
  //priors
  a ~ normal(0,5); //initial productivity - wide prior
  b ~ normal(0,5); //initial productivity - wide prior
  phi ~ uniform(-1,1);
  
  //variance terms
  sigma_e ~ cauchy(0,5);

  R_S ~ normal(mu, sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(R_S|mu[t], sigma_e);
  } 
"

Ricker_nl<- "data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R; //recruits
  vector[N] S; //spawners in time T
}
 transformed data{
 vector[N] log_R; //recruits in the log scale
 log_R = log(R);
 }
parameters {
  real<lower = 0> a;// initial productivity (on log scale)
  real<lower = 0> b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
}
transformed parameters{
  vector[N] log_mu;

  
  for(n in 1:N){
    log_mu[n] = log(S)[n]+a-b*S[n];
  }
}
model{
  //priors
  a ~ cauchy(0, 2); //initial productivity - wide prior
  b ~ cauchy(0, 10); //initial productivity - wide prior
  
  //variance terms
  sigma_e ~ inv_gamma(1, 1);

  log_R ~ normal(log_mu, sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(log_R|log_mu[t], sigma_e);
  } 
"

Ricker_nl2<- "data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R; //recruits
  vector[N] S; //spawners in time T
}
 transformed data{
 vector[N] log_R; //recruits in the log scale
 log_R = log(R);
 }
parameters {
  real<lower = 0> a;// initial productivity (on log scale)
  real<lower = 0> b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
}
transformed parameters{
  vector[N] mu;
  vector[N] log_mu;
  
  for(n in 1:N){
    mu[n] = a*S[n]*exp(b*S[n]);
  }
  log_mu = log(mu);
}
model{
  //priors
  a ~ cauchy(0, 10); //initial productivity - wide prior
  b ~ cauchy(0, 5); //initial productivity - wide prior
  
  //variance terms
  sigma_e ~ inv_gamma(1, 1);

  log_R ~ normal(log_mu, sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(log_R|log_mu[t], sigma_e);
  } 
"

BevHolt<- "data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R; //recruits
  vector[N] S; //spawners in time T
}
 transformed data{
 vector[N] log_R; //recruits in the log scale
 log_R = log(R);
 }
parameters {
  real log_a;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
}

transformed parameters{
  real a = exp(log_a);
  real b = exp(log_b);
  vector[N] mu;
  vector[N] log_mu;
  
  mu = a*S ./(1 + b*S);
  
  log_mu = log(mu);
}
model{
  //priors
  log_a ~ normal(0,10); //initial productivity - wide prior
  log_b ~ normal(0,10); //initial productivity - wide prior
  
  //variance terms
  sigma_e ~ inv_gamma(2, 1);

  log_R ~ normal(log_mu, sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(log_R|log_mu[t], sigma_e);
  } 
"

BevHolt_ac<- "data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R; //recruits
  vector[N] S; //spawners in time T
}
 transformed data{
 vector[N] log_R; //recruits in the log scale
 log_R = log(R);
 }
parameters {
  real log_a;// initial productivity (on log scale)
  real log_b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = -1, upper = 1> phi;
}

transformed parameters{
  real a = exp(log_a);
  real b = exp(log_b);
  vector[N] mu;
  vector[N] log_mu;
  vector[N] epsilon; //residuals

  mu = a*S ./(1 + b*S);
  log_mu = log(mu);
  
  epsilon[1] = log_R[1] - log_mu[1];
  for(t in 2:N){
    epsilon[t] =(log_R[t] - log_mu[t]);
    log_mu[t] = log_mu[t] + (phi*epsilon[t-1]);
  }

}
model{
  //priors
  log_a ~ normal(0,10); //initial productivity - wide prior
  log_b ~ normal(0,10); //initial productivity - wide prior
  phi ~ uniform(-1,1);
  
  //variance terms
  sigma_e ~ cauchy(0,5);

  log_R ~ normal(log_mu, sigma_e);
  
}
generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(log_R|log_mu[t], sigma_e);
  } 
"

BH_2<- "data{
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
"

SR_fit_var_a<- "data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  vector[N] R_S; //log(recruits per spawner)
  vector[N] S; //spawners in time T
 }
parameters {
  real a0;// initial productivity (on log scale)
  real b; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_a;

  //time-varying parameters
  vector[N-1] a_dev; //year-to-year deviations in a
  
}
transformed parameters{
  vector[N] a; //a in each year (on log scale)
   a[1] = a0;
   
  for(t in 2:N){
    a[t] = a[t-1] + a_dev[t-1]*sigma_a;
  }

}  

model{
  //priors
  a0 ~ normal(0,5); //initial productivity - wide prior
  b ~ normal(0, 2); //capacity rate
  a_dev ~ std_normal();
  
  //variance terms
  sigma_e ~ inv_gamma(2, 1);
  sigma_a ~ inv_gamma(2, 1);
   
  to_vector(a_dev) ~ std_normal();
  for(t in 1:N){
    R_S[t] ~ normal(a[t] - S[t]*b, sigma_e);
  }
  
}
  generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal_lpdf(R_S|a[t] - S[t]*b, sigma_e);
  } 
"

SR_fit_var_b<- "data{
  int<lower=1> N;//number of annual samples (time-series length)
  int TT[N];//index of years
  int R_S[TT]; //log(recruits per spawner)
  int S[TT]; //spawners in time T
 }
parameters {
  real a;// initial productivity (on log scale) - fixed in this
  real b0; // rate capacity - fixed in this

 //variance components  
  real<lower = 0> sigma_e;
  real<lower = 0> sigma_b;
  
  //time-varying parameters
  vector[N-1] b_dev; //year-to-year deviations in a
  vector[N] b; //b in each year 
}

transformed parameters{
  vector[N] mu;
  vector[N] epsilon;
  mu[1] = a + b0*S;
  b[1] = b0;
  epsilon[1] = R_S[1] - mu[1];
   
  for(t in 2:N){
    b[t] = b[t-1] + b_dev[t]*sigma_b;
    mu[t] = a + b[t]*S + phi*epsilon[t-1];
    epsilon[t] = (R_S[t] - mu[t]);
  }
}  

model{
  //priors
  a ~ normal(0,5); //initial productivity - wide prior
  beta2 ~ normal(0,2); //covariates - reef
  
  //variance terms
  sigma_e ~ inv_gamma(2, 1);
  sigma_a ~ inv_gamma(2, 1);
   
  b_dev ~ std_normal();
  
  R_S ~ normal(mu, sigma_e);
}
  generated quantities{
  vector[N] log_lik;
  for (t in 1:N) log_lik[t] = normal(mu[t], sigma_e) 
  } 
"

#Load data from above
chinook_dat<- read.csv(here('data','filtered datasets','chinook_final.csv'))

library(rstan);library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

test1<- subset(chinook_dat,stock.id==unique(chinook_dat$stock.id)[1])
test2<- subset(chinook_dat,stock.id==unique(chinook_dat$stock.id)[2])
test3<- subset(chinook_dat,stock.id==unique(chinook_dat$stock.id)[3])

sr_static1<- rstan::stan(model_code = Ricker_linear, data = list(R_S = test1$logR_S,
                                                         N=nrow(test1),
                                                         TT=as.numeric(factor(test1$broodyear)),
                                                         S=c((test1$spawners)/1e3)),
                        pars = c('a','b','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_var_a1<- rstan::stan(model_code = SR_fit_var_a, data = list(R_S = test1$logR_S,
                                                                 N=nrow(test1),
                                                                 TT=as.numeric(factor(test1$broodyear)),
                                                                 S=c((test1$spawners)/1e3)),
                         pars = c('a','b','sigma_a','sigma_e','log_lik'),
                         control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

psis_loo_static1<- loo(extract_log_lik(sr_static1))
psis_loo_var_a1<- loo(extract_log_lik(sr_var_a1))
loo::loo_compare(psis_loo_static1,psis_loo_var_a1)

sr_static2<- rstan::stan(model_code = Ricker_linear, data = list(R_S = test2$logR_S,
                                                                 N=nrow(test2),
                                                                 TT=as.numeric(factor(test2$broodyear)),
                                                                 S=c((test2$spawners)/1e3)),
                         pars = c('a','b','sigma_e','log_lik'),
                         control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_var_a2<- rstan::stan(model_code = SR_fit_var_a, data = list(R_S = test2$logR_S,
                                                               N=nrow(test2),
                                                               TT=as.numeric(factor(test2$broodyear)),
                                                               S=c((test2$spawners)/1e3)),
                        pars = c('a','b','sigma_a','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

psis_loo_static2<- loo(extract_log_lik(sr_static2))
psis_loo_var_a2<- loo(extract_log_lik(sr_var_a2))
loo::loo_compare(psis_loo_static2,psis_loo_var_a2)

sr_static3<- rstan::stan(model_code = Ricker_linear, data = list(R_S = test3$logR_S,
                                                                 N=nrow(test3),
                                                                 TT=as.numeric(factor(test3$broodyear)),
                                                                 S=c((test3$spawners)/1e3)),
                         pars = c('a','b','sigma_e','log_lik'),
                         control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_var_a3<- rstan::stan(model_code = SR_fit_var_a, data = list(R_S = test3$logR_S,
                                                               N=nrow(test3),
                                                               TT=as.numeric(factor(test3$broodyear)),
                                                               S=c((test3$spawners)/1e3)),
                        pars = c('a','b','sigma_a','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

psis_loo_static3<- loo(extract_log_lik(sr_static3))
psis_loo_var_a3<- loo(extract_log_lik(sr_var_a3))
loo::loo_compare(psis_loo_static3,psis_loo_var_a3)


sr_static_ac<- rstan::stan(model_code = Ricker_linear_ac, data = list(R_S = test1$logR_S,
                                                         N=nrow(test1),
                                                         TT=as.numeric(factor(test1$broodyear)),
                                                         S=c((test1$spawners)/1e3)),
                        pars = c('a','b','phi','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_static_nl<- rstan::stan(model_code = Ricker_nl, data = list(R = test1$recruits,
                                                                N=nrow(test1),
                                                                TT=as.numeric(factor(test1$broodyear)),
                                                                S=c((test1$spawners)/1e3)),
                        pars = c('a','log_a','b','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_static_nl2<- rstan::stan(model_code = Ricker_nl2, data = list(R = test1$recruits,
                                                               N=nrow(test1),
                                                               TT=as.numeric(factor(test1$broodyear)),
                                                               S=c((test1$spawners)/1e3)),
                           pars = c('a','b','sigma_e','log_lik'),
                           control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)


bh_static<- rstan::stan(model_code = BevHolt, data = list(R = test1$recruits,
                                                               N=nrow(test1),
                                                               TT=as.numeric(factor(test1$broodyear)),
                                                               S=c((test1$spawners))),
                           pars = c('a','b','sigma_e','log_lik'),
                           control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

bh_static_ac<- rstan::stan(model_code = BevHolt_ac, data = list(R = test1$recruits,
                                                          N=nrow(test1),
                                                          TT=as.numeric(factor(test1$broodyear)),
                                                          S=c((test1$spawners))),
                        pars = c('a','b','phi','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

#Testing out logR_S vs. R nl parameterizations
plot(recruits~spawners,data=test_dat,bty='l',cex=1.5,pch=21)
params_sr_s<- rstan::extract(sr_static)
params_sr_s_nl<- rstan::extract(sr_static_nl)

r_s_pred<- exp(median(params_sr_s$a)-median(params_sr_s$b)*c(seq(0,25000)/1e3))
r_pred_nl<- exp(median(params_sr_s$a))*c(seq(0,25000))*exp(-median(params_sr_s_nl$b)/1e3*c(seq(0,25000)))
r_pred<- r_s_pred*seq(0,25000)
lines(r_pred~seq(0,25000))

r_pred_nl<- median(params_sr_s_nl$a)*c(seq(0,25000))*exp(-median(params_sr_s_nl$b)/1e3*c(seq(0,25000)))
lines(r_pred_nl~seq(0,25000),col='darkred')

lines(median(params_sr_s$a)*c(seq(0,25000)/1e3)*exp(median(-params_sr_s$b)*c(seq(0,25000)/1e3)))

bh_static_ac<- rstan::stan(model_code = BevHolt_ac, data = list(R = test_dat$recruits,
                                                          N=nrow(test_dat),
                                                          TT=as.numeric(factor(test_dat$broodyear)),
                                                          S=c((test_dat$spawners))),
                        pars = c('a','b','phi','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

mod_test2_nl<- rstan::stan(model_code = BH, data = list(R = test_dat$recruits,
                                                        N=nrow(test_dat),
                                                        TT=as.numeric(factor(test_dat$broodyear)),
                                                        S=c((test_dat$spawners))),
                           pars = c('a','b','sigma_e','log_lik'),
                           control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

shinystan::launch_shinystan(mod_test1_nl)


lm_m_test<- lm(logR_S~c(scale(spawners)),data=test_dat)
vcov(lm_m_test)
summary(lm_m_test)

plot(logR_S~spawners,data=test_dat)
lines(lm_m_test$coefficients[1]+lm_m_test$coefficients[2]*c(seq(0,25000)/1e3))
lines()


mod_test2<- rstan::stan(model_code = SR_fit_var_a, data = list(R_S = test_dat$logR_S,
                                                               N=nrow(test_dat),
                                                               TT=as.numeric(factor(test_dat$broodyear)),
                                                               S=c(scale(test_dat$spawners))),
                        pars = c('a','b','sigma_a','sigma_e','sigma_a','a_dev','log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1, init=0)

params_2<- extract(mod_test2)

elpd1= loo::loo(extract_log_lik(mod_test1))
elpd2= loo::loo(extract_log_lik(mod_test2))
elpd_comp<- loo::loo_compare(elpd1,elpd2)


shinystan::launch_shinystan(mod_test)


###BASIC SHIT I DONT GET LETS DO IT###
##Ricker functions
a=4
b=1e-3
S=seq(0,5000)
R=a*S*exp(-b*S)
plot(R~S,type='l')
A=log(a)

R2=S*exp(A-b*S)
lines(R2~S,lty=5)

plot(R~R2)

R_S = exp(A - b*S)
R3 = R_S*S

log_R<- log(S)+A-b*S
exp(log_R)
plot(exp(log_R)~R2)
log_R2<- log(a)+log(S)-b*S

#all of these line up - phew

#with real data
library(nlme)
R=test1$recruits
S=test1$spawners
r1<- nls(R~S*exp(A-b*S),start=list(A=4,b=1e-4))
r2<- nls(R~a*S*exp(-b*S),start=list(a=1,b=1e-4))
r1$m$getPars()[1]
log(r2$m$getPars()[1]) #check

log_r1<- nls(log(R)~log(S/1e3)+A-b*S/1e3,start=list(A=1,b=1e-4))
log_r2<- nls(log(R)~log(a)+log(S)-b*S,start=list(a=1,b=1e-4))
log_r1$m$getPars()[1]
log(log_r2$m$getPars()[1]) #check
