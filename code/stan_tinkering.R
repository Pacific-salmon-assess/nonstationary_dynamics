library(here)

#Load data from above
sock_dat<- read.csv(here('data','filtered datasets','sockeye_final.csv'))

library(rstan);library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

test1<- subset(sock_dat,stock.id==unique(sock_dat$stock.id)[1])
test2<- subset(sock_dat,stock.id==unique(sock_dat$stock.id)[2])
test3<- subset(sock_dat,stock.id==unique(sock_dat$stock.id)[3])

sr_static1<- rstan::stan(file = here('code','stan models','ricker_linear.stan'), data = list(R_S = test1$logR_S,
                                                                                                       N=nrow(test1),
                                                                                                       TT=as.numeric(factor(test1$broodyear)),
                                                                                                       S=c((test1$spawners/1e5))),
                         pars = c('log_a','log_b','b','sigma_e','log_lik'),
                         control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_static1a<- rstan::stan(file = here('code','stan models','ricker_linear_ac_resids.stan'), data = list(R_S = test1$logR_S,
                                                         N=nrow(test1),
                                                         TT=as.numeric(factor(test1$broodyear)),
                                                         S=c((test1$spawners/1e5))),
                        pars = c('log_a','log_b','b','phi','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_var_a1<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a.stan'), 
                                                    data = list(R_S = test1$logR_S,
                                                                 N=nrow(test1),
                                                                 TT=as.numeric(factor(test1$broodyear)),
                                                                 S=c((test1$spawners/1e5))),
                         pars = c('log_a','b','sigma_a','sigma_e','log_lik'),
                         control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)

sr_var_b1<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_b.stan'), 
                        data = list(R_S = test1$logR_S,
                                    N=nrow(test1),
                                    TT=as.numeric(factor(test1$broodyear)),
                                    S=c((test1$spawners/1e5))),
                        pars = c('log_a','log_b','b','sigma_b','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)

sr_var_a_b1<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a_and_b.stan'), 
                        data = list(R_S = test1$logR_S,
                                    N=nrow(test1),
                                    TT=as.numeric(factor(test1$broodyear)),
                                    S=c((test1$spawners/1e5))),
                        pars = c('log_a','log_b','b','sigma_a','sigma_b','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)

sr_var_a_b2<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a_and_b2.stan'), 
                          data = list(R_S = test1$logR_S,
                                      N=nrow(test1),
                                      TT=as.numeric(factor(test1$broodyear)),
                                      S=c((test1$spawners/1e5))),
                          pars = c('sigma_a_b','sigma_e','log_a','log_b','Cor_1','b','z_mat','log_lik'),
                          control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)

psis_loo_static1<- loo(extract_log_lik(sr_static1))
psis_loo_static1a<- loo(extract_log_lik(sr_static1a))
psis_loo_var_a1<- loo(extract_log_lik(sr_var_a1))
psis_loo_var_b1<- loo(extract_log_lik(sr_var_b1))
psis_loo_var_a_b1<- loo(extract_log_lik(sr_var_a_b1))
psis_loo_var_a_b2<- loo(extract_log_lik(sr_var_a_b2))
loo::loo_compare(psis_loo_static1,psis_loo_static1a,psis_loo_var_a1,psis_loo_var_b1,psis_loo_var_a_b1,psis_loo_var_a_b2)


sr_static2<- rstan::stan(file = here('code','stan models','ricker_linear.stan'), data = list(R_S = test2$logR_S,
                                                                                             N=nrow(test2),
                                                                                             TT=as.numeric(factor(test2$broodyear)),
                                                                                             S=c((test2$spawners/1e5))),
                         pars = c('log_a','log_b','b','sigma_e','log_lik'),
                         control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_static2a<- rstan::stan(file = here('code','stan models','ricker_linear_ac_resids.stan'), data = list(R_S = test2$logR_S,
                                                                                                        N=nrow(test2),
                                                                                                        TT=as.numeric(factor(test2$broodyear)),
                                                                                                        S=c((test2$spawners/1e5))),
                          pars = c('log_a','log_b','b','phi','sigma_e','log_lik'),
                          control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_var_a2<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a.stan'), 
                        data = list(R_S = test2$logR_S,
                                    N=nrow(test2),
                                    TT=as.numeric(factor(test2$broodyear)),
                                    S=c((test2$spawners/1e5))),
                        pars = c('log_a','b','sigma_a','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)

sr_var_b2<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_b.stan'), 
                        data = list(R_S = test2$logR_S,
                                    N=nrow(test2),
                                    TT=as.numeric(factor(test2$broodyear)),
                                    S=c((test2$spawners/1e5))),
                        pars = c('log_a','log_b','b','sigma_b','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)

sr_var_a_b2_1<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a_and_b.stan'), 
                          data = list(R_S = test2$logR_S,
                                      N=nrow(test2),
                                      TT=as.numeric(factor(test2$broodyear)),
                                      S=c((test2$spawners/1e5))),
                          pars = c('log_a','log_b','b','sigma_a','sigma_b','sigma_e','log_lik'),
                          control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)

sr_var_a_b2_2<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a_and_b2.stan'), 
                          data = list(R_S = test2$logR_S,
                                      N=nrow(test2),
                                      TT=as.numeric(factor(test2$broodyear)),
                                      S=c((test2$spawners/1e5))),
                          pars = c('sigma_a_b','sigma_e','log_a','log_b','Cor_1','b','z_mat','log_lik'),
                          control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)

loo1<- loo(extract_log_lik(sr_static2))
psis_loo_static2a<- loo(extract_log_lik(sr_static2a))
loo2<- loo(extract_log_lik(sr_var_a2))
loo3<- loo(extract_log_lik(sr_var_b2))
loo4<- loo(extract_log_lik(sr_var_a_b2_1))
psis_loo_var_a_b2_2<- loo(extract_log_lik(sr_var_a_b2_2))
loo::loo_compare(psis_loo_static2,psis_loo_static2a,psis_loo_var_a2,psis_loo_var_b2,psis_loo_var_a_b2_1,psis_loo_var_a_b2_2)


lpd_point <- cbind(
  loo1$pointwise[,"elpd_loo"],
  loo2$pointwise[,"elpd_loo"],
  loo3$pointwise[,"elpd_loo"],
  loo4$pointwise[,"elpd_loo"]
)



sr_static3<- rstan::stan(file = here('code','stan models','ricker_linear.stan'), data = list(R_S = test3$logR_S,
                                                                                             N=nrow(test3),
                                                                                             TT=as.numeric(factor(test3$broodyear)),
                                                                                             S=c((test3$spawners/1e5))),
                         pars = c('log_a','log_b','b','sigma_e','log_lik'),
                         control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_static3a<- rstan::stan(file = here('code','stan models','ricker_linear_ac_resids.stan'), data = list(R_S = test3$logR_S,
                                                                                                        N=nrow(test3),
                                                                                                        TT=as.numeric(factor(test3$broodyear)),
                                                                                                        S=c((test3$spawners/1e5))),
                          pars = c('log_a','log_b','b','phi','sigma_e','log_lik'),
                          control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_var_a3<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a.stan'), 
                        data = list(R_S = test3$logR_S,
                                    N=nrow(test3),
                                    TT=as.numeric(factor(test3$broodyear)),
                                    S=c((test3$spawners/1e5))),
                        pars = c('log_a','b','sigma_a','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)

sr_var_b3<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_b.stan'), 
                        data = list(R_S = test3$logR_S,
                                    N=nrow(test3),
                                    TT=as.numeric(factor(test3$broodyear)),
                                    S=c((test3$spawners/1e5))),
                        pars = c('log_a','log_b','b','sigma_b','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)

sr_var_a_b3_1<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a_and_b.stan'), 
                            data = list(R_S = test3$logR_S,
                                        N=nrow(test3),
                                        TT=as.numeric(factor(test3$broodyear)),
                                        S=c((test3$spawners/1e5))),
                            pars = c('log_a','log_b','b','sigma_a','sigma_b','sigma_e','log_lik'),
                            control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)

sr_var_a_b3_2<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a_and_b2.stan'), 
                            data = list(R_S = test3$logR_S,
                                        N=nrow(test3),
                                        TT=as.numeric(factor(test3$broodyear)),
                                        S=c((test3$spawners/1e5))),
                            pars = c('sigma_a_b','sigma_e','log_a','log_b','Cor_1','b','z_mat','log_lik'),
                            control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)

psis_loo_static3<- loo(extract_log_lik(sr_static2))
psis_loo_static3a<- loo(extract_log_lik(sr_static2a))
psis_loo_var_a3<- loo(extract_log_lik(sr_var_a2))
psis_loo_var_b3<- loo(extract_log_lik(sr_var_b2))
psis_loo_var_a_b3_1<- loo(extract_log_lik(sr_var_a_b3_1))
psis_loo_var_a_b3_2<- loo(extract_log_lik(sr_var_a_b3_2))
loo::loo_compare(psis_loo_static3,psis_loo_static3a,psis_loo_var_a3,psis_loo_var_b3,psis_loo_var_a_b3_1,psis_loo_var_a_b3_2)


sr_var_a1_ac<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a_ac_resids.stan'), 
                        data = list(R_S = test1$logR_S,
                                    N=nrow(test1),
                                    TT=as.numeric(factor(test1$broodyear)),
                                    S=c((test1$spawners))/1e3),
                        pars = c('a','b','phi','sigma_a','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 1, iter = 800, thin = 1)

sr_var_a1<- rstan::stan(model_code = SR_fit_var_a, data = list(R_S = test1$logR_S,
                                                               N=nrow(test1),
                                                               TT=as.numeric(factor(test1$broodyear)),
                                                               S=c((test1$spawners)/1e3)),
                        pars = c('a','b','sigma_a','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)


psis_loo_static1<- loo(extract_log_lik(sr_static1))
psis_loo_var_a1<- loo(extract_log_lik(sr_var_a1))
loo::loo_compare(psis_loo_static1,psis_loo_var_a1)

sr_static2<- rstan::stan(file = here('code','stan models','ricker_linear.stan'), data = list(R_S = test2$logR_S,
                                                                 N=nrow(test2),
                                                                 TT=as.numeric(factor(test2$broodyear)),
                                                                 S=c((test2$spawners)/1e5)),
                         pars = c('a','b','sigma_e','log_lik'),
                         control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_var_a2<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a.stan'), data = list(R_S = test2$logR_S,
                                                               N=nrow(test2),
                                                               TT=as.numeric(factor(test2$broodyear)),
                                                               S=c((test2$spawners)/1e5)),
                        pars = c('a','b','sigma_a','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

sr_var_b2<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_b.stan'), 
                        data = list(R_S = test2$logR_S,
                                    N=nrow(test2),
                                    TT=as.numeric(factor(test2$broodyear)),
                                    S=c((test2$spawners/1e5))),
                        pars = c('a','b','sigma_b','sigma_e','log_lik'),
                        control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1,init=0)


psis_loo_static2<- loo(extract_log_lik(sr_static2))
psis_loo_var_a2<- loo(extract_log_lik(sr_var_a2))
psis_loo_var_b2<- loo(extract_log_lik(sr_var_b2))
loo::loo_compare(psis_loo_static2,psis_loo_var_a2,psis_loo_var_b2)

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


log_bh1<- nls(R~(a*S)/(1+b*S),start=list(a=0,b=0))
log_r2<- nls(log(R)~log(a)+log(S)-b*S,start=list(a=1,b=1e-4))

