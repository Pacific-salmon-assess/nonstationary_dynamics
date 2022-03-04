#Preliminary batched model support for varying dynamics in sockeye
library(here)
sock_dat<- read.csv(here('data','filtered datasets','sockeye_final.csv'))
sock_info<- read.csv(here('data','sockeye','sockeye_info.csv'))
library(rstan);library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

sock_info<- subset(sock_info, Stock.ID %in% sock_dat$stock.id)

for(i in 1:nrow(sock_info)){
  s<- subset(sock_dat,stock.id==sock_info$Stock.ID[i])
  
  #Model 1: static a & b
  mod1<- rstan::stan(file = here('code','stan models','ricker_linear.stan'), data = list(R_S = s$logR_S,
                                                                                               N=nrow(s),
                                                                                               TT=as.numeric(factor(s$broodyear)),
                                                                                               S=c((s$spawners/1e5))),
                           pars = c('log_a','b','log_b','sigma_e','log_lik'),
                           control = list(adapt_delta = 0.999,max_treedepth = 15), warmup = 1000, chains = 10, iter = 2000, thin = 1)
  
  mod2<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a.stan'), 
                          data = list(R_S = s$logR_S,
                                      N=nrow(s),
                                      TT=as.numeric(factor(s$broodyear)),
                                      S=c((s$spawners/1e5))),
                          pars = c('log_a','b','log_b','sigma_a','sigma_e','log_lik'),
                          control = list(adapt_delta = 0.999,max_treedepth = 15), warmup = 1000, chains = 10, iter = 2000, thin = 1,init=0)
  
  mod3<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_b.stan'), 
                     data = list(R_S = s$logR_S,
                                 N=nrow(s),
                                 TT=as.numeric(factor(s$broodyear)),
                                 S=c((s$spawners/1e5))),
                     pars = c('log_a','b','log_b','sigma_b','sigma_e','log_lik'),
                     control = list(adapt_delta = 0.999,max_treedepth = 15), warmup = 1000, chains = 10, iter = 2000, thin = 1,init=0)
 
   mod4<- rstan::stan(file =  here('code','stan models','ricker_linear_varying_a_and_b.stan'), 
                     data = list(R_S = s$logR_S,
                                 N=nrow(s),
                                 TT=as.numeric(factor(s$broodyear)),
                                 S=c((s$spawners/1e5))),
                     pars = c('log_a','b','log_b','sigma_a','sigma_b','sigma_e','log_lik'),
                     control = list(adapt_delta = 0.999,max_treedepth = 15), warmup = 1000, chains = 10, iter = 2000, thin = 1,init=0)
   
   loo1<- loo(extract_log_lik(mod1))
   loo2<- loo(extract_log_lik(mod2))
   loo3<- loo(extract_log_lik(mod3))
   loo4<- loo(extract_log_lik(mod4))
   
   lpd_point <- cbind(
     loo1$pointwise[,"elpd_loo"],
     loo2$pointwise[,"elpd_loo"],
     loo3$pointwise[,"elpd_loo"],
     loo4$pointwise[,"elpd_loo"]
   )
   
   weights=stacking_weights(lpd_point)
   sock_info$w1[i]=weights[1]
   sock_info$w2[i]=weights[2]
   sock_info$w3[i]=weights[3]
   sock_info$w4[i]=weights[4]
   
   params1<- rstan::extract(mod1)
   params2<- rstan::extract(mod2)
   params3<- rstan::extract(mod3)
   params4<- rstan::extract(mod4)
   
   pars_mod1=c('log_a','b','log_b','sigma_e')
   pars_mod2=c('log_a','b','log_b','sigma_e','sigma_a')
   pars_mod3=c('log_a','b','log_b','sigma_e','sigma_b')
   pars_mod4=c('log_a','b','log_b','sigma_e','sigma_a','sigma_b')
   
   mod_par_path_1<- here('outputs','initial stan runs','sockeye','1 - static')
   mod_par_path_2<- here('outputs','initial stan runs','sockeye','2 - a')
   mod_par_path_3<- here('outputs','initial stan runs','sockeye','3 - b')
   mod_par_path_4<- here('outputs','initial stan runs','sockeye','4 - a and b')
   
   write.csv(as.data.frame(params1),file.path(mod_par_path_1,paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model1','.csv',sep='')))
   write.csv(as.data.frame(params2),file.path(mod_par_path_2,paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model2','.csv',sep='')))
   write.csv(as.data.frame(params3),file.path(mod_par_path_3,paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model3','.csv',sep='')))
   write.csv(as.data.frame(params4),file.path(mod_par_path_4,paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model4','.csv',sep='')))
   
   mod_par_path_sum<- here('outputs','initial stan runs','sockeye','model summaries')
   write.csv(as.data.frame(summary(mod1,pars = pars_mod1)$summary),file.path(mod_par_path_sum,paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model1_summary','.csv',sep='')))
   write.csv(as.data.frame(summary(mod2,pars = pars_mod2)$summary),file.path(mod_par_path_sum,paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model2_summary','.csv',sep='')))
   write.csv(as.data.frame(summary(mod3,pars = pars_mod3)$summary),file.path(mod_par_path_sum,paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model3_summary','.csv',sep='')))
   write.csv(as.data.frame(summary(mod4,pars = pars_mod4)$summary),file.path(mod_par_path_sum,paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model4_summary','.csv',sep='')))
   
}
