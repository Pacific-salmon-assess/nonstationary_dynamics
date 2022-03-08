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

sock_info$mod_sec=apply(sock_info[,14:17],1,which.max)
summary(factor(sock_info$mod_sec))

for(i in 1:nrow(sock_info)){
  if(sock_info$mod_sec[i]==1){
    params=read.csv(here('outputs','initial stan runs','sockeye','1 - static',paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model1','.csv',sep='')))
  
    log_a= params[,(gsub('\\..*','',colnames(params))=='log_a')]
    b= params[,(gsub('\\..*','',colnames(params))=='b')]
    a.med=median(log_a)
    a.l80=quantile(log_a,0.1)
    a.u80=quantile(log_a,0.9)
    
    b.med=median(b)
    b.l80=quantile(b,0.1)
    b.u80=quantile(b,0.9)
    }
  if(sock_info$mod_sec[i]==2){
    params=read.csv(here('outputs','initial stan runs','sockeye','2 - a',paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model2','.csv',sep='')))
  
    log_a= params[,(gsub('\\..*','',colnames(params))=='log_a')]
    b= params[,(gsub('\\..*','',colnames(params))=='b')]
    a.med=apply(log_a,2,median)
    a.l80=apply(log_a,2, quantile, probs=c(0.1))
    a.u80=apply(log_a,2, quantile, probs=c(0.9))
    
    b.med=median(b)
    b.l80=quantile(b,0.1)
    b.u80=quantile(b,0.9)
    
    }
  if(sock_info$mod_sec[i]==3){
    params=read.csv(here('outputs','initial stan runs','sockeye','3 - b',paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model3','.csv',sep='')))
  
    log_a= params[,(gsub('\\..*','',colnames(params))=='log_a')]
    b= params[,(gsub('\\..*','',colnames(params))=='b')]
    a.med=median(log_a)
    a.l80=quantile(log_a,0.1)
    a.u80=quantile(log_a,0.9)
    
    b.med=apply(b,2,median)
    b.l80=apply(b,2, quantile, probs=c(0.1))
    b.u80=apply(b,2, quantile, probs=c(0.9))
    
    
    }
  if(sock_info$mod_sec[i]==4){
    params=read.csv(here('outputs','initial stan runs','sockeye','4 - a and b',paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model4','.csv',sep='')))
  
    log_a= params[,(gsub('\\..*','',colnames(params))=='log_a')]
    b= params[,(gsub('\\..*','',colnames(params))=='b')]
    a.med=apply(log_a,2,median)
    a.l80=apply(log_a,2, quantile, probs=c(0.1))
    a.u80=apply(log_a,2, quantile, probs=c(0.9))
    
    b.med=apply(b,2,median)
    b.l80=apply(b,2, quantile, probs=c(0.1))
    b.u80=apply(b,2, quantile, probs=c(0.9))
    }
 
  
  pdf(file.path(here('outputs','initial stan runs','sockeye','plots'),paste('Best fit',sock_info$Stock[i],sock_info$Species[i],sep='_','.pdf')),width=14,height=8.5)
  par(mfrow=c(1,2))
  plot(a.med,type='l',ylim=c(min(a.l80),max(a.u80)),ylab='',main=paste('Productivity -',sock_info$Stock[i],sock_info$Species[i],sep=' '))
  lines(a.l80,lty=5);lines(a.u80,lty=5)
  plot(b.med,type='l',ylim=c(min(b.l80),max(b.u80)),ylab='',main=paste('Capacity -',sock_info$Stock[i],sock_info$Species[i],sep=' '))
  lines(b.l80,lty=5);lines(b.u80,lty=5)
  dev.off()

}
