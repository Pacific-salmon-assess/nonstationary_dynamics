#Preliminary batched model support for varying dynamics in sockeye
library(here);library(dplyr)
sock_dat<- read.csv(here('data','filtered datasets','sockeye_final.csv'))
sock_info<- read.csv(here('data','filtered datasets','filtered_sockeye_info.csv'))
library(cmdstanr);library(loo);library(dlm)

source(here('code','plot_functions.R'))
source(here('code','dlm-wrapper.R'))

sock_info<- subset(sock_info, stock.id %in% sock_dat$stock.id)
weights=data.frame(w1=NA,w2=NA,w2gp=NA,w3=NA,w3gp=NA,w4=NA,w4gp=NA)
sock_info<- cbind(sock_info,weights)
####Stan models####
set_cmdstan_path()

file1 <- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear.stan")
#file1_nl <- file.path(cmdstan_path(), "ricker_linear_nolik.stan") #for testing priors no longer needed
file2 <- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_varying_a.stan") #random walk log_a model
file2gp <- file.path(cmdstan_path(), 'nonstationary dynamics',"ricker_linear_varying_a_GP.stan")
file2ou <- file.path(cmdstan_path(), 'nonstationary dynamics',"ricker_linear_varying_a_OU.stan") #Ornstein-Uhlbeck process
file3 <- file.path(cmdstan_path(),'nonstationary dynamics', "ricker_linear_varying_b.stan")
file3t<- file.path(cmdstan_path(), 'nonstationary dynamics',"ricker_linear_varying_b2.stan")
file3gp <- file.path(cmdstan_path(),'nonstationary dynamics', "ricker_linear_varying_b_GP.stan")
file4 <- file.path(cmdstan_path(),'nonstationary dynamics', "ricker_linear_varying_a_and_b.stan")
file4.2 <- file.path(cmdstan_path(),'nonstationary dynamics', "ricker_linear_varying_a_and_b2.stan")
file4gp <- file.path(cmdstan_path(),'nonstationary dynamics', "ricker_linear_varying_a_and_b_GP.stan")
#test <- file.path(cmdstan_path(), "nonstationary dynamics","test_mod2.stan") #2-lag model, doesn't work well at all

mod1 <- cmdstan_model(file1)
mod2 <- cmdstan_model(file2)
mod2gp <- cmdstan_model(file2gp)
mod2ou <- cmdstan_model(file2ou)
mod3 <- cmdstan_model(file3)
mod3_t<- cmdstan_model(file3t)
mod3gp <- cmdstan_model(file3gp)
mod4 <- cmdstan_model(file4)
mod4.2 <- cmdstan_model(file4.2)
mod4gp <- cmdstan_model(file4gp)
#mtest<- cmdstan_model(test)

#testing out some alternative model fits
for(i in 1:nrow(sock_info)){
  s<- subset(sock_dat,stock.id==sock_info$stock.id[i])
  
  data=list(R_S = s$logR_S,
            N=nrow(s),
            TT=as.numeric(factor(s$broodyear)),
            S=c((s$spawners)))
  SRdata <- data.frame(byr=s$broodyear,
                       spwn=s$spawners,
                       rec=s$recruits)
  
  
  #Static SR model
  fit1<- mod1$sample(
    data = data,
    seed = 123, 
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 1000,
    iter_sampling = 2000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
  
  #TV productivity
  fit2<- mod2$sample(
    data = data,
    seed = 123, 
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
  #TV productivity - GP
  fit2gp<- mod2gp$sample(
    data = data,
    init=0,
    seed = 123, 
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
  #TV productivity - GP
  fit2ou<- mod2ou$sample(
    data = data,
    seed = 123, 
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
  avary<-fitDLM(SRdata, alpha_vary = TRUE, beta_vary = FALSE)
  
  
 # fit2test<- mtest$sample(
 #   data = data,
#    init=0,
#    seed = 123, 
 #   chains = 6, 
 #   parallel_chains = 6,
 #   iter_warmup = 500,
 #   iter_sampling = 1000,
#    refresh = 500,
 #   adapt_delta = 0.99,
 #   max_treedepth = 20 # print update every 500 iters
 # )
  
 # params2test<- fit2test$draws(format='df',variables=c('log_a','b','log_b','sigma_a','sigma_e','sigma_v'))
  #fit2test$summary(variables=c('sigma_a','sigma_e','sigma_v'))
  
  #TV capacity
  fit3<- mod3$sample(
    data = data,
    seed = 123, 
    init=0,
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
  bvary<-fitDLM(SRdata, alpha_vary = FALSE, beta_vary = TRUE)
  
 # fit3_test<- mod3_t$sample(
  #  data = data,
  #  seed = 123, 
   # init=0,
  #chains = 6, 
   # parallel_chains = 6,
    #iter_warmup = 500,
    #iter_sampling = 1000,
    #refresh = 500,
    #adapt_delta = 0.99,
    #max_treedepth = 20 # print update every 500 iters
  #)
  
  #TV capacity - GP
  fit3gp<- mod3gp$sample(
    data = data,
    seed = 123,
    init=0,
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
  #TV productivity & capacity
  fit4<- mod4$sample(
    data = data,
    init=0,
    seed = 123, 
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
  #correlated walks with a & b
#  fit4.2<- mod4.2$sample(
#    data = data,
#    init=0,
 #   seed = 123, 
#    chains = 6, 
#    parallel_chains = 6,
 #   iter_warmup = 500,
 #   iter_sampling = 1000,
 #   refresh = 500,
  #  adapt_delta = 0.99,
  #  max_treedepth = 20 # print update every 500 iters
 # )
  #TV productivity & capacity - GP
  fit4gp<- mod4gp$sample(
    data = data,
    init=0,
    seed = 123, 
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
  abvary<-fitDLM(SRdata, alpha_vary = TRUE, beta_vary = TRUE)
  
 # elpd1= fit1$loo(cores=2)
#  elpd2= fit2$loo(cores=2)
 # elpd2gp= fit2gp$loo(cores=2)
#  elpd3= fit3$loo(cores=2)
 # elpd3gp= fit3gp$loo(cores=2)
#  elpd4= fit4$loo(cores=2)
 # elpd4gp= fit4gp$loo(cores=2)
  #elpd_comp<- loo::loo_compare(elpd1,elpd2,elpd2gp,elpd3,elpd3gp,elpd4,elpd4gp)
  
 # lpd_point <- cbind(
#    elpd1$pointwise[,"elpd_loo"],
#    elpd2$pointwise[,"elpd_loo"],
 #   elpd2gp$pointwise[,"elpd_loo"],
#    elpd3$pointwise[,"elpd_loo"],
 #   elpd3gp$pointwise[,"elpd_loo"],
#    elpd4$pointwise[,"elpd_loo"],
 #   elpd4gp$pointwise[,"elpd_loo"]
#  )
#  mod_weights<-stacking_weights(lpd_point)
  
 #  weights=stacking_weights(lpd_point)
#   sock_info$w1[i]=weights[1]
#   sock_info$w2[i]=weights[2]
#   sock_info$w2gp[i]=weights[3]
 #  sock_info$w3[i]=weights[4]
#   sock_info$w3gp[i]=weights[5]
#  sock_info$w4[i]=weights[6]
 #  sock_info$w4gp[i]=weights[7]
   
   params1<- fit1$draws(format='df',variables=c('log_a','b','log_b'))
   params2<- fit2$draws(format='df',variables=c('log_a','b','log_b'))
   params2gp<- fit2gp$draws(format='df',variables=c('log_a','b','log_b'))
   params2ou<- fit2ou$draws(format='df',variables=c('log_a'))
   params3<- fit3$draws(format='df',variables=c('log_a','b','log_b'))
   params3gp<- fit3gp$draws(format='df',variables=c('log_a','b','log_b'))
   params4<- fit4$draws(format='df',variables=c('log_a','b','log_b'))
   params4.2<- fit4.2$draws(format='df',variables=c('log_a','b','log_b'))
   params4gp<- fit4gp$draws(format='df',variables=c('log_a','b','log_b'))
   
   pars_mod1=c('log_a','b','log_b','sigma_e')
   pars_mod2=c('log_a','b','log_b','sigma_e','sigma_a')
   pars_mod2gp=c('log_a','b','log_b','sigma_e','gp_tau','gp_rho')
   pars_mod3=c('log_a','b','log_b','sigma_e','sigma_b')
   pars_mod4=c('log_a','b','log_b','sigma_e','sigma_a','sigma_b')
   pars_mod4gp=c('log_a','b','log_b','sigma_e','gp_tau_a','gp_rho_a','gp_tau_b','gp_rho_b')
   
 #  mod_par_path_1<- here('outputs','initial stan runs','sockeye','1 - static')
#   mod_par_path_2<- here('outputs','initial stan runs','sockeye','2 - a')
 #  mod_par_path_3<- here('outputs','initial stan runs','sockeye','3 - b')
  # mod_par_path_4<- here('outputs','initial stan runs','sockeye','4 - a and b')
   
  # write.csv(fit1$draws(format='df',variables=pars_mod1),file.path(mod_par_path_1,paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model1','.csv',sep='')))
  # write.csv(fit2$draws(format='df',variables=pars_mod1),file.path(mod_par_path_2,paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model2','.csv',sep='')))
  # write.csv(fit3$draws(format='df',variables=pars_mod1),file.path(mod_par_path_3,paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model3','.csv',sep='')))
   #write.csv(fit4$draws(format='df',variables=pars_mod1),file.path(mod_par_path_4,paste(sprintf("%02d",i),'_',sock_info$Stock[i],sock_info$Species[i],'_model4','.csv',sep='')))
   
  # mod_par_path_sum<- here('outputs','initial stan runs','sockeye','model summaries')
  # write.csv(as.data.frame(fit1$summary(variables=pars_mod1)),file.path(mod_par_path_sum,paste(sprintf("%02d",i),'_',s$stock[i],s$species[i],'_model1_summary','.csv',sep='')))
  # write.csv(as.data.frame(fit2$summary(variables=pars_mod2)),file.path(mod_par_path_sum,paste(sprintf("%02d",i),'_',s$stock[i],s$species[i],'_model2_summary','.csv',sep='')))
  # write.csv(as.data.frame(fit2gp$summary(variables=pars_mod2gp)),file.path(mod_par_path_sum,paste(sprintf("%02d",i),'_',s$stock[i],s$species[i],'_model2gp_summary','.csv',sep='')))
  # write.csv(as.data.frame(fit3$summary(variables=pars_mod3)),file.path(mod_par_path_sum,paste(sprintf("%02d",i),'_',s$stock[i],s$species[i],'_model3_summary','.csv',sep='')))
  # write.csv(as.data.frame(fit3gp$summary(variables=pars_mod2gp)),file.path(mod_par_path_sum,paste(sprintf("%02d",i),'_',s$stock[i],s$species[i],'_model3gp_summary','.csv',sep='')))
  # write.csv(as.data.frame(fit4$summary(variables=pars_mod4)),file.path(mod_par_path_sum,paste(sprintf("%02d",i),'_',s$stock[i],s$species[i],'_model4_summary','.csv',sep='')))
  # write.csv(as.data.frame(fit4gp$summary(variables=pars_mod4gp)),file.path(mod_par_path_sum,paste(sprintf("%02d",i),'_',s$stock[i],s$species[i],'_model4gp_summary','.csv',sep='')))
   
   rw_gp_dlm_ou_plot_comp(params=params2,params_gp = params2gp,params_dlm=avary$results,type=1,x=s,pdf=1)
   rw_gp_dlm_plot_comp(params=params2,params_gp = params2gp,params_dlm=avary$results,type=1,x=s,pdf=1)
   rw_gp_dlm_plot_comp(params=params3,params_gp = params3gp,params_dlm=bvary$results,type=2,x=s,pdf=1)
   rw_gp_dlm_plot_comp(params=params4,params_gp = params4gp,params_dlm=abvary$results,type=3,x=s,pdf=1)
}


for(i in 1:nrow(sock_info)){
  s<- subset(sock_dat,stock.id==sock_info$stock.id[i])
  
  data=list(R_S = s$logR_S,
            N=nrow(s),
            TT=as.numeric(factor(s$broodyear)),
            S=c((s$spawners)))
  SRdata <- data.frame(byr=s$broodyear,
                       spwn=s$spawners,
                       rec=s$recruits)
  
  #TV productivity
  fit2<- mod2$sample(
    data = data,
    seed = 123, 
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
  #TV productivity - GP
  fit2gp<- mod2gp$sample(
    data = data,
    init=0,
    seed = 123, 
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
  #TV productivity - GP
  fit2ou<- mod2ou$sample(
    data = data,
    seed = 123, 
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
  avary<-fitDLM(SRdata, alpha_vary = TRUE, beta_vary = FALSE)
  
  params2<- fit2$draws(format='df',variables=c('log_a','b','log_b'))
  params2gp<- fit2gp$draws(format='df',variables=c('log_a','b','log_b'))
  params2ou<- fit2ou$draws(format='df',variables=c('log_a'))
  
  rw_gp_dlm_ou_plot_comp(params=params2,params_gp = params2gp,params_dlm=avary$results,params_ou=params2ou,x=s,pdf=0)

}

#correlated vs uncorrelated random walks
summ_corr=data.frame(stock=sock_info$stock,elpd1=NA,elpd2=NA,elpd_diff=NA,elpd_diff_se=NA,w1=NA,w2=NA,cor.median=NA,cor.l90=NA,cor.u90=NA)
for(i in 2:nrow(sock_info)){
  s<- subset(sock_dat,stock.id==sock_info$stock.id[i])
  
  data=list(R_S = s$logR_S,
            N=nrow(s),
            TT=as.numeric(factor(s$broodyear)),
            S=c((s$spawners)))
  SRdata <- data.frame(byr=s$broodyear,
                       spwn=s$spawners,
                       rec=s$recruits)
  
  #TV productivity & capacity
  fit4<- mod4$sample(
    data = data,
    init=0,
    seed = 123, 
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
#  correlated walks with a & b
    fit4.2<- mod4.2$sample(
      data = data,
      init=0,
     seed = 123, 
      chains = 6, 
      parallel_chains = 6,
     iter_warmup = 500,
     iter_sampling = 1000,
     refresh = 500,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
   )
  
  elpd4= fit4$loo(cores=2)
  elpd4.2= fit4.2$loo(cores=2)
  elpd_comp<- loo::loo_compare(elpd4,elpd4.2)
  
  lpd_point <- cbind(
      elpd4$pointwise[,"elpd_loo"],
     elpd4.2$pointwise[,"elpd_loo"]
    )
  weights=stacking_weights(lpd_point)
  
  summ_corr[i,2]=elpd4$estimates[1,1]
  summ_corr[i,3]=elpd4.2$estimates[1,1]
  summ_corr[i,4]=elpd_comp[2,1]
  summ_corr[i,5]=elpd_comp[2,2]
  summ_corr[i,6]=weights[1]
  summ_corr[i,7]=weights[2]
  
 
  params4<- fit4$draws(format='df',variables=c('log_a','b','log_b'))
  params4.2<- fit4.2$draws(format='df',variables=c('log_a','b','log_b','Cor_1'))
  
  summ_corr[i,8]=median(params4.2$`Cor_1[1,2]`)
  summ_corr[i,9]=quantile(params4.2$`Cor_1[1,2]`,0.05)
  summ_corr[i,10]=quantile(params4.2$`Cor_1[1,2]`,0.95)
  
  prod_cap_corr_comp(params1=params4,params2 = params4.2,x=s,pdf=1)
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
  plot(a.med~seq(sock_info$ts.start[i],sock_info$ts.end[i]),type='l',ylim=c(min(a.l80),max(a.u80)),ylab='',main=paste('Productivity -',sock_info$Stock[i],sock_info$Species[i],sep=' '))
  lines(a.l80,lty=5);lines(a.u80,lty=5)
  plot(b.med~seq(sock_info$ts.start[i],sock_info$ts.end[i]),type='l',ylim=c(min(b.l80),max(b.u80)),ylab='',main=paste('Capacity -',sock_info$Stock[i],sock_info$Species[i],sep=' '))
  lines(b.l80,lty=5);lines(b.u80,lty=5)
  dev.off()

}








log_a = rnorm(1000,0,2.5)
b = rnorm(1000,-12,3)
sigma=rgamma(100,2,3)


fit1_nl<- mod1_nl$sample(
  data = data,
  seed = 123, 
  chains = 6, 
  parallel_chains = 6,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 500,
  adapt_delta = 0.99,
  max_treedepth = 20 # print update every 500 iters
)
bayesplot::ppc_dens_overlay(data$R_S, as.matrix(fit1_nl$draws(format='df',variables='y_rep')[,1:44]))

nl_dat=fit1_nl$draws(format='df',variables='y_rep')[,1:44]
plot(data$R_S,type='n',bty='l',ylim=c(-10,10),ylab='log(R/S)')
nl_dat_sub=nl_dat[sample(nrow(nl_dat),1000),]
for(i in 1:1000){
  lines(as.numeric(nl_dat_sub[i,])~seq(1:ncol(nl_dat_sub)),lwd=0.5,col='lightblue')
}
lines(data$R_S,lwd=2,col='navy')

par(mfrow=c(2,1))
bayesplot::ppc_dens_overlay(data$R_S, as.matrix(fit1$draws(format='df',variables='y_rep')[,1:44]))
bayesplot::mcmc_hist(fit1$draws("y_rep"), binwidth = 0.025)

