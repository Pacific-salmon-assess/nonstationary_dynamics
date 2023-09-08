library(here);library(dplyr);library(ggplot2)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_may2023.csv'))
stock_info<- read.csv(here('data','filtered datasets','all_stocks_info_may2023.csv'))

#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)

library(samEst)
options(mc.cores = parallel::detectCores())
#1


###Load in data####
#Remove stocks with less than 15 years of recruitment data
stock_info_filtered=subset(stock_info,n.years>=16) #252 stocks
stock_info_filtered$stock.name=gsub('/','_',stock_info_filtered$stock.name)
stock_info_filtered$stock.name=gsub('&','and',stock_info_filtered$stock.name)

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock.id)) #267
stock_info_filtered$stock.id2=seq(1:nrow(stock_info_filtered))

if(any(stock_dat2$spawners==0)){stock_dat2$spawners=stock_dat2$spawners+1;stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)}
if(any(stock_dat2$recruits==0)){stock_dat2$recruits=stock_dat2$recruits+1;stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)}
stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)
stock_dat2=stock_dat2[complete.cases(stock_dat2$logR_S),]

#Model fits and figures
for(i in 1:nrow(stock_info_filtered)){
  
  
}



###LFO-CV####
#Define models (helps prevent crashing)
m1=samEst::sr_mod(type='static',ac = FALSE,par='n',lfo=TRUE)
m2=samEst::sr_mod(type='static',ac = TRUE,par='n',lfo=TRUE)
m3=samEst::sr_mod(type='rw',par='a',lfo=TRUE)
m4=samEst::sr_mod(type='rw',par='b',lfo=TRUE)
m5=samEst::sr_mod(type='rw',par='both',lfo=TRUE)
m6=samEst::sr_mod(type='hmm',par='a',lfo=TRUE)
m7=samEst::sr_mod(type='hmm',par='b',lfo=TRUE)
m8=samEst::sr_mod(type='hmm',par='both',lfo=TRUE)


#Summary data frame and store for pointwise LLs
lfo_summary=data.frame(stock=stock_info_filtered$stock.name,LL_m1=NA,LL_m2=NA,LL_m3.1=NA,LL_m3.3=NA,
                          LL_m3.5=NA,LL_m4.1=NA,LL_m4.3=NA,LL_m4.5=NA,LL_m5.1=NA,LL_m5.3=NA,LL_m5.5=NA,LL_m6.1=NA,LL_m6.3=NA,LL_m6.5=NA,
                         LL_m7.1=NA,LL_m7.3=NA,LL_m7.5=NA,
                          LL_m8.1=NA,LL_m8.3=NA,LL_m8.5=NA)
full_lfo=list() #full pointwise log likelihoods for each year by model type
mod_lfo=list() #top variant for each model type
for(i in 2:10){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  #Assess model fits for each model type
  #model 1 - static Ricker
  ll1<- samEst::stan_lfo_cv(mod=m1,type='static',df=df,L=10)
  #model 2 - static autocorrelated Ricker
  ll2<- samEst::stan_lfo_cv(mod=m2,type='static',df=df,L=10)
  #model 3 - dynamic productivity Ricker
  ll3<- samEst::stan_lfo_cv(mod=m3,type='tv',df=df,L=10)
  
  #model 4 - dynamic capacity Ricker
  ll4<- samEst::stan_lfo_cv(mod=m4,type='tv',df=df,L=10)
  
  #model 5 - dynamic productivity & capacity Ricker
  ll5<- samEst::stan_lfo_cv(mod=m5,type='tv',df=df,L=10)
  
  #model 6 - productivity regime shift - 2 regimes
  ll6<- samEst::stan_lfo_cv(mod=m6,type='regime',df=df,L=10,K=2)
  
  #model 7 - capacity regime shift
  ll7<- samEst::stan_lfo_cv(mod=m7,type='regime',df=df,L=10,K=2)
  
  #model 8 - productivity and capacity regime shift
  ll8<- samEst::stan_lfo_cv(mod=m8,type='regime',df=df,L=10,K=2)
 #Take the sum of the estimated pointwise likelihood estimates (=elpd_loo)
  lfo_summary[i,2]=sum(ll1)
  lfo_summary[i,3]=sum(ll2)
  lfo_summary[i,4:6]=apply(ll3,1,sum)
  lfo_summary[i,7:9]=apply(ll4,1,sum)
  lfo_summary[i,10:12]=apply(ll5,1,sum)
  lfo_summary[i,13:15]=apply(ll6,1,sum)
  lfo_summary[i,16:18]=apply(ll7,1,sum)
  lfo_summary[i,19:21]=apply(ll8,1,sum)

  full_lfo[[i]]=rbind(ll1,ll2,ll3,ll4,ll5,ll6,ll7,ll8)
  rownames(full_lfo[[i]])=c('m1','m2',paste('m3',c(1,3,5),sep="_"),paste('m4',c(1,3,5),sep="_"),paste('m5',c(1,3,5),sep="_"),paste('m6',c(1,3,5),sep="_"),paste('m7',c(1,3,5),sep="_"),paste('m8',c(1,3,5),sep="_"))
  wm3=which.max(apply(ll3,1,sum))#best fit for model 3
  wm4=which.max(apply(ll4,1,sum))#best fit for model 4
  wm5=which.max(apply(ll5,1,sum)) #best fit for model 5
  wm6=which.max(apply(ll6,1,sum)) #best fit for model 6
  wm7=which.max(apply(ll7,1,sum)) #best fit for model 7
  wm8=which.max(apply(ll8,1,sum)) #best fit for model 8 (two combinations of higher alpha, lower cap; higher alpha, higher cap)

  mod_lfo[[i]]=rbind(ll1,ll2,ll3[wm3,],ll4[wm4,],ll5[wm5,],ll6[wm6,],ll7[wm7,],ll8[wm8,])
  rownames(mod_lfo[[i]])[1:8]=paste('m',seq(1:8),sep='');rownames(mod_lfo[[i]])[8]='m8'
  
  write.csv(full_lfo[[i]],here('outputs','full_LFO_loglik','l10',paste(sprintf("%02d",i),'_',gsub(' ','_',stock_info_filtered$stock.name[i]),'full_lfo.csv',sep='')))
  write.csv(mod_lfo[[i]],here('outputs','model_LFO_loglik','l10',paste(sprintf("%02d",i),'_',stock_info_filtered$stock.name[i],'mod_lfo.csv',sep='')))
  
}

#Summarize results
f=list.files(here('outputs','model_LFO_lfo'))
s_i=as.numeric(gsub("[^0-9]", "", substr(f, 1, 3)))
m_ll=list()
full_mw=matrix(ncol=9,nrow=length(s_i))
full_mw[,1]=s_i
d90_mw=matrix(ncol=9,nrow=length(s_i))
d90_mw[,1]=s_i
d80_mw=matrix(ncol=9,nrow=length(s_i))
d80_mw[,1]=s_i
for(i in 1:length(f)){
  m_ll[[i]]<- read.csv(here('outputs','model_LFO_lfo',f[i]))
  full_mw[i,2:9]=model_weights(m_ll[[i]][,2:ncol(m_ll[[i]])],type='full')
  d90_mw[i,2:9]=model_weights(m_ll[[i]][,2:ncol(m_ll[[i]])],type='d90')
  d80_mw[i,2:9]=model_weights(m_ll[[i]][,2:ncol(m_ll[[i]])],type='d80')
  }

full_mw=full_mw[match(stock_info_filtered$stock.id2,full_mw[,1]),]
best_mod=apply(full_mw[,2:9],1,which.max) #this isn't working for some reason :{
sum_mod=summary(factor(best_mod))

library(formattable)
lfo_ms_df=data.frame(stock_info_filtered[,3:8],round(full_mw[,2:9],2),best_mod)
lfo_ms_df=lfo_ms_df[with(lfo_ms_df,order(species,lat)),]
write.csv(lfo_ms_df,here('outputs','ms_rmd','lfo_table.csv'))

formattable(lfo_ms_df,
            list(area(col = 7:14) ~ color_tile("transparent", "red")))

d90_mw=d90_mw[order(s_i),]
d80_mw=d80_mw[order(s_i),]

model_weights(m_ll[[1]])

#AIC/BIC/LOO stan####
#Define models (helps prevent crashing)
m1f=samEst::sr_mod(type='static',ac = FALSE,par='n',lfo =F)
m2f=samEst::sr_mod(type='static',ac = TRUE,par='n',lfo=F)
m3f=samEst::sr_mod(type='rw',par='a',lfo=F)
m4f=samEst::sr_mod(type='rw',par='b',lfo=F)
m5f=samEst::sr_mod(type='rw',par='both',lfo=F)
m6f=samEst::sr_mod(type='hmm',par='a',lfo=F)
m7f=samEst::sr_mod(type='hmm',par='b',lfo=F)
m8f=samEst::sr_mod(type='hmm',par='both',lfo=F)

#plots
for(i in 1:nrow(stock_info_filtered)){
  dir.create(here('outputs','figures','sr plots',paste(i,stock_info_filtered$stock.name[i],sep='-')))
}

aic_weights=matrix(ncol=8,nrow=nrow(stock_info_filtered))
bic_weights=matrix(ncol=8,nrow=nrow(stock_info_filtered))
loo_pbma_weights=matrix(ncol=8,nrow=nrow(stock_info_filtered))
loo_stack_weights=matrix(ncol=8,nrow=nrow(stock_info_filtered))
aic_weights_d90=matrix(ncol=8,nrow=nrow(stock_info_filtered))
bic_weights_d90=matrix(ncol=8,nrow=nrow(stock_info_filtered))
aic_weights_d80=matrix(ncol=8,nrow=nrow(stock_info_filtered))
bic_weights_d80=matrix(ncol=8,nrow=nrow(stock_info_filtered))
loo_pbma_weights_d90=matrix(ncol=8,nrow=nrow(stock_info_filtered))
loo_stack_weights_d90=matrix(ncol=8,nrow=nrow(stock_info_filtered))
loo_pbma_weights_d80=matrix(ncol=8,nrow=nrow(stock_info_filtered))
loo_stack_weights_d80=matrix(ncol=8,nrow=nrow(stock_info_filtered))
for(i in 2:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  s<- s[complete.cases(s$spawners),]
  
  df=data.frame(S=s$spawners,R=s$recruits,by=s$broodyear)
  #Fit each model
  #model 1 - static Ricker
  f1 = rstan::sampling(m1f, 
                       data = list(N=nrow(s),
                                   L=max(s$broodyear)-min(s$broodyear)+1,
                                   ii=s$broodyear-min(s$broodyear)+1,
                                   R_S=s$logR_S,
                                   S=s$spawners),
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 4, iter = 700)
  
  samEst::sr_plot(df=df,mod=f1,title=stock_info_filtered$stock.name[i],make.pdf=T,path=here('outputs','figures','sr plots',paste(i,stock_info_filtered$stock.name[i],sep='-')),type='static',form='stan',par='n')
dev.off()
  
#model 2 - static autocorrelated Ricker
  f2 = rstan::sampling(m2f, 
                       data = list(N=nrow(s),
                                   L=max(s$broodyear)-min(s$broodyear)+1,
                                   ii=s$broodyear-min(s$broodyear)+1,
                                   R_S=s$logR_S,
                                   S=s$spawners),
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 4, iter = 700)
  
  sr_plot(df=df,mod=f2,title=stock_info_filtered$stock.name[i],make.pdf=T,path=here('outputs','figures','sr plots',paste(i,stock_info_filtered$stock.name[i],sep='-')),type='static',form='stan',par='n',ac=TRUE)
  dev.off()
  
  #model 3 - dynamic productivity Ricker
  f3 = rstan::sampling(m3f, 
                       data = list(N=nrow(s),
                                   L=max(s$broodyear)-min(s$broodyear)+1,
                                   ii=s$broodyear-min(s$broodyear)+1,
                                   R_S=s$logR_S,
                                   S=s$spawners),
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 4, iter = 700)
  
 
  
  sr_plot(df=df,mod=f3,title=stock_info_filtered$stock.name[i],make.pdf=T,path=here('outputs','figures','sr plots',paste(i,stock_info_filtered$stock.name[i],sep='-')),type='rw',form='stan',par='a')
  dev.off()
  #model 4 - dynamic capacity Ricker
  f4 = rstan::sampling(m4f, 
                       data = list(N=nrow(s),
                                   L=max(s$broodyear)-min(s$broodyear)+1,
                                   ii=s$broodyear-min(s$broodyear)+1,
                                   R_S=s$logR_S,
                                   S=s$spawners),
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 4, iter = 700)
  
  sr_plot(df=df,mod=f4,title=stock_info_filtered$stock.name[i],make.pdf=T,path=here('outputs','figures','sr plots',paste(i,stock_info_filtered$stock.name[i],sep='-')),type='rw',form='stan',par='b')
  dev.off()
  
  #model 5 - dynamic productivity & capacity Ricker
  f5 = rstan::sampling(m5f, 
                       data = list(N=nrow(s),
                                   L=max(s$broodyear)-min(s$broodyear)+1,
                                   ii=s$broodyear-min(s$broodyear)+1,
                                   R_S=s$logR_S,
                                   S=s$spawners),
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 4, iter = 700)
  
  sr_plot(df=df,mod=f5,title=stock_info_filtered$stock.name[i],make.pdf=T,path=here('outputs','figures','sr plots',paste(i,stock_info_filtered$stock.name[i],sep='-')),type='rw',form='stan',par='both')
  dev.off()
  
  #model 6 - productivity regime shift - 2 regimes
  f6 = rstan::sampling(m6f, 
                       data = list(N=nrow(s),
                                   R_S=s$logR_S,
                                   S=s$spawners,
                                   K=2,
                                   alpha_dirichlet=rep(1,2)), #prior for state transition probabilities (this makes them equal)
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 4, iter = 700)
 
  sr_plot(df=df,mod=f6,title=stock_info_filtered$stock.name[i],make.pdf=T,path=here('outputs','figures','sr plots',paste(i,stock_info_filtered$stock.name[i],sep='-')),type='hmm',form='stan',par='a')
  dev.off()
  
  #model 7 - capacity regime shift
  f7 = rstan::sampling(m7f, 
                       data = list(N=nrow(s),
                                   R_S=s$logR_S,
                                   S=s$spawners,
                                   K=2,
                                   alpha_dirichlet=rep(1,2)), #prior for state transition probabilities (this makes them equal)
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 4, iter = 700)
  
 sr_plot(df=df,mod=f7,title=stock_info_filtered$stock.name[i],make.pdf=T,path=here('outputs','figures','sr plots',paste(i,stock_info_filtered$stock.name[i],sep='-')),type='hmm',form='stan',par='b')
 dev.off()
  
  #model 8 - productivity and capacity regime shift
  f8 = rstan::sampling(m8f, 
                       data = list(N=nrow(s),
                                   R_S=s$logR_S,
                                   S=s$spawners,
                                   K=2,
                                   alpha_dirichlet=rep(1,2)), #prior for state transition probabilities (this makes them equal)
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 4, iter = 700)
  
 sr_plot(df=df,mod=f8,title=stock_info_filtered$stock.name[i],make.pdf=T,path=here('outputs','figures','sr plots',paste(i,stock_info_filtered$stock.name[i],sep='-')),type='hmm',form='stan',par='both')
 dev.off()
   
   d1=rstan::extract(f1)
   d2=rstan::extract(f2)
   d3=rstan::extract(f3)
   d4=rstan::extract(f4)
   d5=rstan::extract(f5)
   d6=rstan::extract(f6)
   d7=rstan::extract(f7)
   d8=rstan::extract(f8)
   
   loo1=loo::loo(f1)
   loo2=loo::loo(f2)
   loo3=loo::loo(f3)
   loo4=loo::loo(f4)
   loo5=loo::loo(f5)
   loo6=loo::loo(f6)
   loo7=loo::loo(f7)
   loo8=loo::loo(f8)
   
   
   lpd_point <- cbind(
     loo1$pointwise[,"elpd_loo"], 
     loo2$pointwise[,"elpd_loo"],
     loo3$pointwise[,"elpd_loo"], 
     loo4$pointwise[,"elpd_loo"],
     loo5$pointwise[,"elpd_loo"],
     loo6$pointwise[,"elpd_loo"],
     loo7$pointwise[,"elpd_loo"],
     loo8$pointwise[,"elpd_loo"]
   )
   
   dl=list(d1$log_lik,d2$log_lik,d3$log_lik,d4$log_lik,d5$log_lik,d6$log_lik,d7$log_lik,d8$log_lik)
   
   
   aic_weights[i,]=stan_aic(dl,form='aic',type='full',k=c(3,4,4,4,5,5,5,6))
   bic_weights[i,]=stan_aic(dl,form='bic',type='full',k=c(3,4,4,4,5,5,5,6))
   aic_weights_d90[i,]=stan_aic(dl,form='aic',type='d90',k=c(3,4,5,5,6,5,5,6))
   bic_weights_d90[i,]=stan_aic(dl,form='bic',type='d90',k=c(3,4,5,5,6,5,5,6))
   aic_weights_d80[i,]=stan_aic(dl,form='aic',type='d80',k=c(3,4,5,5,6,5,5,6))
   bic_weights_d80[i,]=stan_aic(dl,form='bic',type='d80',k=c(3,4,5,5,6,5,5,6))
   loo_pbma_weights[i,]=loo::pseudobma_weights(lpd_point)
   loo_stack_weights[i,]=loo::stacking_weights(lpd_point)
   loo_pbma_weights_d90[i,]=loo::pseudobma_weights(lpd_point[apply(lpd_point,1,mean)>=quantile(apply(lpd_point,1,mean),0.1),])
   loo_stack_weights_d90[i,]=loo::stacking_weights(lpd_point[apply(lpd_point,1,mean)>=quantile(apply(lpd_point,1,mean),0.1),])
   loo_pbma_weights_d90[i,]=loo::pseudobma_weights(lpd_point[apply(lpd_point,1,mean)>=quantile(apply(lpd_point,1,mean),0.2),])
   loo_stack_weights_d90[i,]=loo::stacking_weights(lpd_point[apply(lpd_point,1,mean)>=quantile(apply(lpd_point,1,mean),0.2),])
   
   print(i)
}


best_mod=apply(aic_weights,1,which.max) #this isn't working for some reason :{
stan_aic_df=data.frame(stock_info_filtered[,3:8],round(aic_weights,2),best_mod)
write.csv(stan_aic_df,here('outputs','ms_rmd','stan_aic_table_weights.csv'))

best_mod=apply(bic_weights,1,which.max) #this isn't working for some reason :{
stan_bic_df=data.frame(stock_info_filtered[,3:8],round(bic_weights,2),best_mod)
write.csv(stan_bic_df,here('outputs','ms_rmd','stan_bic_table_weights.csv'))

best_mod=apply(loo_stack_weights,1,which.max) #this isn't working for some reason :{
stan_looic_df=data.frame(stock_info_filtered[,3:8],round(loo_stack_weights,2),best_mod)
write.csv(stan_looic_df,here('outputs','ms_rmd','stan_looic_table_stackweight.csv'))


best_mod=apply(aic_weights_d80,1,which.max) #this isn't working for some reason :{
sum_mod=summary(factor(best_mod))
sum_mod

stan_looic_df_bmaw=data.frame(stock_info_filtered[,3:8],round(bma_weights,2),best_mod)
write.csv(stan_looic_df_bmaw,here('outputs','ms_rmd','stan_looic_table_bmaweight.csv'))

best_mod=apply(bma_weightsl30,1,which.max) #this isn't working for some reason :{
sum_mod=summary(factor(best_mod))
sum_mod

stan_looic_df_stackl30=data.frame(stock_info_filtered[,3:8],round(stack_weightsl30,2),best_mod)
write.csv(stan_looic_df_stackl30,here('outputs','ms_rmd','stan_looic_table_stackweightl30.csv'))


best_mod=apply(stack_weightsl30,1,which.max) #this isn't working for some reason :{
sum_mod=summary(factor(best_mod))
sum_mod

stan_looic_df_stackl30=data.frame(stock_info_filtered[,3:8],round(stack_weightsl30,2),best_mod)
write.csv(stan_looic_df_stackl30,here('outputs','ms_rmd','stan_looic_table_bmaweightl30.csv'))





#TMB runs - LFO/AIC/BIC####
for(u in 1:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  #lfo comparison
  lfostatic<-samEst::tmb_mod_lfo_cv(data=df,model='static', L=10)
  lfoac <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='staticAC', L=10),error = function(e) {lfoac=list(lastparam=rep(-999,length(lfoac$lastparam)))})
  lfoalpha <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='rw_a', siglfo="obs", L=10),error = function(e) {lfoalpha=list(lastparam=rep(-999,length(lfoac$lastparam)), 
                                                                                                                   last3param=rep(-999,length(lfoac$lastparam)), 
                                                                                                                   last5param=rep(-999,length(lfoac$lastparam)))})
  lfobeta <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='rw_b', siglfo="obs", L=10),error = function(e) {lfobeta=list(lastparam=rep(-999,length(lfoac$lastparam)), 
                                                                                                                 last3param=rep(-999,length(lfoac$lastparam)), 
                                                                                                                 last5param=rep(-999,length(lfoac$lastparam)))})
  lfoalphabeta <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='rw_both', siglfo="obs", L=10),error = function(e) {lfoalphabeta=list(lastparam=rep(-999,length(lfoac$lastparam)), 
                                                                                                                              last3param=rep(-999,length(lfoac$lastparam)), 
                                                                                                                              last5param=rep(-999,length(lfoac$lastparam)))})
  lfohmma <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='HMM_a', L=10),error = function(e) {lfohmma=list(lastregime_pick=rep(-999,length(lfoac$lastparam)), 
                                                                                                    last3regime_pick=rep(-999,length(lfoac$lastparam)), 
                                                                                                    last5regime_pick=rep(-999,length(lfoac$lastparam)))})
  lfohmmb <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='HMM_b', L=10),error = function(e) {lfohmmb=list(lastregime_pick=rep(-999,length(lfoac$lastparam)), 
                                                                                                                   last3regime_pick=rep(-999,length(lfoac$lastparam)), 
                                                                                                                   last5regime_pick=rep(-999,length(lfoac$lastparam)))})
  lfohmm <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='HMM', L=10),error = function(e) {lfohmm=list(lastregime_pick=rep(-999,length(lfoac$lastparam)), 
                                                                                                               last3regime_pick=rep(-999,length(lfoac$lastparam)), 
                                                                                                               last5regime_pick=rep(-999,length(lfoac$lastparam)))})
  TMBstatic <- ricker_TMB(data=df)
  TMBac <- ricker_TMB(data=df, AC=TRUE)
  TMBtva <- tryCatch(samEst::ricker_rw_TMB(data=df,tv.par='a'),error = function(e) {TMBtva=list(conv_problem=TRUE)})
  TMBtvb <- tryCatch(samEst::ricker_rw_TMB(data=df, tv.par='b'),error = function(e){TMBtvb=list(conv_problem=TRUE)})
  TMBtvab <- tryCatch(samEst::ricker_rw_TMB(data=df, tv.par='both'),error = function(e){TMBtvab=list(conv_problem=TRUE)})
  TMBhmma <- tryCatch(samEst::ricker_hmm_TMB(data=df, tv.par='a'),error = function(e){TMBhmma=list(conv_problem=TRUE)})
  TMBhmmb <- tryCatch(samEst::ricker_hmm_TMB(data=df, tv.par='b'),error = function(e){TMBhmmb=list(conv_problem=TRUE)})
  TMBhmm  <- tryCatch(samEst::ricker_hmm_TMB(data=df, tv.par='both'),error = function(e) {TMBhmm=list(conv_problem=TRUE)})
  
  sr_plot(df=df,mod=TMBstatic,title=stock_info_filtered$stock.name[u],make.pdf=F,path=here('outputs','figures','sr plots',stock_info_filtered$stock.name[u]),type='static',form='tmb',par='n')
  sr_plot(df=df,mod=TMBac,title=stock_info_filtered$stock.name[i],make.pdf=F,path=here('outputs','figures','sr plots',stock_info_filtered$stock.name[u]),type='static',form='tmb',par='n')
  sr_plot(df=df,mod=TMBtva,title=stock_info_filtered$stock.name[i],make.pdf=F,path=here('outputs','figures','sr plots',stock_info_filtered$stock.name[u]),type='rw',form='tmb',par='a')
  sr_plot(df=df,mod=TMBtvb,title=stock_info_filtered$stock.name[i],make.pdf=F,path=here('outputs','figures','sr plots',stock_info_filtered$stock.name[u]),type='rw',form='tmb',par='b')
  sr_plot(df=df,mod=TMBtvab,title=stock_info_filtered$stock.name[i],make.pdf=F,path=here('outputs','figures','sr plots',stock_info_filtered$stock.name[u]),type='rw',form='tmb',par='both')
  sr_plot(df=df,mod=TMBhmma,title=stock_info_filtered$stock.name[i],make.pdf=F,path=here('outputs','figures','sr plots',stock_info_filtered$stock.name[u]),type='hmm',form='tmb',par='a')
  sr_plot(df=df,mod=TMBhmmb,title=stock_info_filtered$stock.name[i],make.pdf=F,path=here('outputs','figures','sr plots',stock_info_filtered$stock.name[u]),type='hmm',form='tmb',par='b')
  sr_plot(df=df,mod=TMBhmm,title=stock_info_filtered$stock.name[i],make.pdf=F,path=here('outputs','figures','sr plots',stock_info_filtered$stock.name[u]),type='hmm',form='tmb',par='both')
  
  LLdf<-rbind(lfostatic$lastparam,lfoac$lastparam,
              lfoalpha$lastparam,lfoalpha$last3param,lfoalpha$last5param,
              lfobeta$lastparam,lfobeta$last3param,lfobeta$last5param,
              lfoalphabeta$lastparam,lfoalphabeta$last3param,lfoalphabeta$last5param,
              lfohmma$lastregime_pick,lfohmma$last3regime_pick,lfohmma$last5regime_pick,
              lfohmmb$lastregime_pick,lfohmmb$last3regime_pick,lfohmmb$last5regime_pick,
              lfohmm$lastregime_pick,lfohmm$last3regime_pick,lfohmm$last5regime_pick
  )
  
  
  llfdf_b=rbind(
  LLdf[1:2,],
  LLdf[3:5,][which.max(apply(LLdf[3:5,],1,sum)),],#best fit for model 3
  LLdf[6:8,][which.max(apply(LLdf[6:8,],1,sum)),],#best fit for model 4
  LLdf[9:11,][which.max(apply(LLdf[9:11,],1,sum)),], #best fit for model 5
  LLdf[12:14,][which.max(apply(LLdf[12:14,],1,sum)),], #best fit for model 6
  LLdf[15:17,][which.max(apply(LLdf[15:17,],1,sum)),], #best fit for model 7
  LLdf[18:20,][which.max(apply(LLdf[18:20,],1,sum)),] #best fit for model 8 (two combinations of higher alpha, lower cap; higher alpha, higher cap)
  )
  
  
  LLdf_1b<-rbind(lfostatic$lastparam,lfoac$lastparam,
              lfoalpha$lastparam,
              lfobeta$lastparam,
              lfoalphabeta$lastparam,
              lfohmma$lastregime_pick,
              lfohmmb$lastregime_pick,
              lfohmm$lastregime_pick)
              
  
  LLdf_3b<-rbind(lfostatic$lastparam,lfoac$lastparam,
                 lfoalpha$last3param,
                 lfobeta$last3param,
                 lfoalphabeta$last3param,
                 lfohmma$last3regime_pick,
                 lfohmmb$last3regime_pick,
                 lfohmm$last3regime_pick
  )
  
  LLdf_5b<-rbind(lfostatic$lastparam,lfoac$lastparam,
                 lfoalpha$lastparam,lfoalpha$last3param,lfoalpha$last5param,
                 lfobeta$lastparam,lfobeta$last3param,lfobeta$last5param,
                 lfoalphabeta$lastparam,lfoalphabeta$last3param,lfoalphabeta$last5param,
                 lfohmma$lastregime_pick,lfohmma$last3regime_pick,lfohmma$last5regime_pick,
                 lfohmma$lastregime_average,lfohmma$last3regime_average,lfohmma$last5regime_average,
                 lfohmmb$lastregime_pick,lfohmmb$last3regime_pick,lfohmmb$last5regime_pick,
                 lfohmmb$lastregime_average,lfohmmb$last3regime_average,lfohmmb$last5regime_average,
                 lfohmm$lastregime_pick,lfohmm$last3regime_pick,lfohmm$last5regime_pick,
                 lfohmm$lastregime_average,lfohmm$last3regime_average,lfohmm$last5regime_average
  )
  
  AICdf<-c(ifelse(TMBstatic$conv_problem,999,TMBstatic$AICc),
               ifelse(TMBac$conv_problem,999,TMBac$AICc),
               ifelse(TMBtva$conv_problem,999,TMBtva$AICc),
               ifelse(TMBtvb$conv_problem,999,TMBtvb$AICc),
               ifelse(TMBtvab$conv_problem,999,TMBtvab$AICc),
               ifelse(TMBhmma$conv_problem,999,TMBhmma$AICc),
               ifelse(TMBhmmb$conv_problem,999,TMBhmmb$AICc),
               ifelse(TMBhmm$conv_problem,999,TMBhmm$AICc))
  
  BICdf<-c(ifelse(TMBstatic$conv_problem,999,TMBstatic$BIC),
               ifelse(TMBac$conv_problem,999,TMBac$BIC),
               ifelse(TMBtva$conv_problem,999,TMBtva$BIC),
               ifelse(TMBtvb$conv_problem,999,TMBtvb$BIC),
               ifelse(TMBtvab$conv_problem,999,TMBtvab$BIC),
               ifelse(TMBhmma$conv_problem,999,TMBhmma$BIC),
               ifelse(TMBhmmb$conv_problem,999,TMBhmmb$BIC),
               ifelse(TMBhmm$conv_problem,999,TMBhmm$BIC))
  
  write.csv(llfdf_b,here('outputs','TMB LFO','top model set',paste(u,'..',gsub(' ','_',stock_info_filtered$stock.name[u]),'_top_lfo.csv',sep='')))
  write.csv(LLdf,here('outputs','TMB LFO','full model set',paste(u,'..',gsub(' ','_',stock_info_filtered$stock.name[u]),'_top_lfo.csv',sep='')))
  write.csv(LLdf_1b,here('outputs','TMB LFO','1year back set',paste(u,'..',gsub(' ','_',stock_info_filtered$stock.name[u]),'_top_lfo.csv',sep='')))
  write.csv(LLdf_3b,here('outputs','TMB LFO','3year back set',paste(u,'..',gsub(' ','_',stock_info_filtered$stock.name[u]),'_top_lfo.csv',sep='')))
  write.csv(LLdf_5b,here('outputs','TMB LFO','5year back set',paste(u,'..',gsub(' ','_',stock_info_filtered$stock.name[u]),'_top_lfo.csv',sep='')))
  
  write.csv(AICdf,here('outputs','TMB AIC',paste(u,'..',gsub(' ','_',stock_info_filtered$stock.name[u]),'_top_lfo.csv',sep='')))
  write.csv(BICdf,here('outputs','TMB BIC',paste(u,'..',gsub(' ','_',stock_info_filtered$stock.name[u]),'_top_lfo.csv',sep='')))
}

#sig_a/sig_obs
sig_a=NA
sig=NA
for(u in 1:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  TMBac <- samEst::ricker_TMB(data=df, AC=TRUE)
  TMBtva <- tryCatch(samEst::ricker_rw_TMB(data=df,tv.par='a'),error = function(e) {TMBtva=list(conv_problem=TRUE)})
  
  sig[u]=TMBtva$sig
  sig_a[u]=TMBac$sigar
  
}

hist(sig,breaks=30,col='darkgray',border='white',xlim=c(0,2),main='',xlab=expression(sigma),cex.lab=1.2)
abline(v=quantile(sig,0.5),lty=5)
abline(v=quantile(sig,0.25),lty=5)
abline(v=quantile(sig,0.75),lty=5)

hist(sig_a,breaks=30,col='darkgray',border='white',xlim=c(0,2),main='',xlab=expression(sigma),cex.lab=1.2)
abline(v=quantile(sig_a,0.5),lty=5)
abline(v=quantile(sig_a,0.25),lty=5)
abline(v=quantile(sig_a,0.75),lty=5)


#Summarize results - LFO####
f=list.files(here('outputs','TMB LFO','top model set'))
s_i=as.numeric(gsub("[^0-9]", "", substr(f, 1, 3)))
m_ll=list()
full_mw=matrix(ncol=9,nrow=length(s_i))
full_mw[,1]=s_i
d90_mw=matrix(ncol=9,nrow=length(s_i))
d90_mw[,1]=s_i
d80_mw=matrix(ncol=9,nrow=length(s_i))
d80_mw[,1]=s_i

na_f= function(x){
  ifelse(is.na(x)==TRUE,-999,x)
}
inf_f= function(x){
  ifelse(is.finite(x)==FALSE,-999,x)
}

for(i in 1:length(f)){
  m_ll[[i]]<- read.csv(here('outputs','TMB LFO','full model set',f[i]))
  m_ll[[i]]<- apply(m_ll[[i]],2,na_f)
  m_ll[[i]]<- apply(m_ll[[i]],2,inf_f)
  llfdf_b=rbind(
    m_ll[[i]][1:2,],
    m_ll[[i]][3:5,][which.max(apply(m_ll[[i]][3:5,],1,sum)),],#best fit for model 3
    m_ll[[i]][6:8,][which.max(apply(m_ll[[i]][6:8,],1,sum)),],#best fit for model 4
    m_ll[[i]][9:11,][which.max(apply(m_ll[[i]][9:11,],1,sum)),], #best fit for model 5
    m_ll[[i]][12:14,][which.max(apply(m_ll[[i]][12:14,],1,sum)),], #best fit for model 6
    m_ll[[i]][15:17,][which.max(apply(m_ll[[i]][15:17,],1,sum)),], #best fit for model 7
    m_ll[[i]][18:20,][which.max(apply(m_ll[[i]][18:20,],1,sum)),] #best fit for model 8 (two combinations of higher alpha, lower cap; higher alpha, higher cap)
  )
  
  full_mw[i,2:9]=model_weights(llfdf_b[,2:ncol(llfdf_b)],form='PBMA',type='full')
  d90_mw[i,2:9]=model_weights(llfdf_b[,2:ncol(llfdf_b)],form='PBMA',type='d90')
  d80_mw[i,2:9]=model_weights(llfdf_b[,2:ncol(llfdf_b)],form='PBMA',type='d80')
}

full_mw=full_mw[match(stock_info_filtered$stock.id2,full_mw[,1]),]
best_mod=apply(full_mw[,2:9],1,which.max) #this isn't working for some reason :{
sum_mod=summary(factor(best_mod))
sum_mod

tmb_lfo_df=data.frame(stock_info_filtered[,3:8],round(full_mw[,2:9],2),best_mod)
tmb_lfo_df=tmb_lfo_df[with(tmb_lfo_df,order(species,lat)),]
write.csv(tmb_lfo_df,here('outputs','ms_rmd','tmb_lfo_table.csv'))

d90_mw=d90_mw[match(stock_info_filtered$stock.id2,d90_mw[,1]),]
best_mod=apply(d90_mw[,2:9],1,which.max) #this isn't working for some reason :{
sum_mod=summary(factor(best_mod))
sum_mod

tmb_lfo_df_d90=data.frame(stock_info_filtered[,3:8],round(d90_mw[,2:9],2),best_mod)
tmb_lfo_df_d90=tmb_lfo_df_d90[with(tmb_lfo_df_d90,order(species,lat)),]
write.csv(tmb_lfo_df_d90,here('outputs','ms_rmd','tmb_lfo_table_d90.csv'))

d80_mw=d80_mw[match(stock_info_filtered$stock.id2,d80_mw[,1]),]
best_mod=apply(d80_mw[,2:9],1,which.max) #this isn't working for some reason :{
sum_mod=summary(factor(best_mod))
sum_mod

tmb_lfo_df_d80=data.frame(stock_info_filtered[,3:8],round(d80_mw[,2:9],2),best_mod)
tmb_lfo_df_d80=tmb_lfo_df_d80[with(tmb_lfo_df_d80,order(species,lat)),]
write.csv(tmb_lfo_df_d80,here('outputs','ms_rmd','tmb_lfo_table_d80.csv'))

#Summarize results - AIC####
f_aic=list.files(here('outputs','TMB AIC'))
s_i_aic=as.numeric(gsub("[^0-9]", "", substr(f, 1, 3)))
m_ll=list()
full_mw=matrix(ncol=9,nrow=length(s_i))
full_mw[,1]=s_i


for(i in 1:length(f)){
  m_ll[[i]]<- read.csv(here('outputs','TMB AIC',f[i]))
 
  full_mw[i,2:9]=model_weights(m_ll[[i]][,2],form='AIC')
}

full_mw=full_mw[match(stock_info_filtered$stock.id2,full_mw[,1]),]
best_mod=apply(full_mw[,2:9],1,which.max) #this isn't working for some reason :{
sum_mod=summary(factor(best_mod))
sum_mod

aic_df=data.frame(stock_info_filtered[,3:8],round(full_mw[,2:9],2),best_mod)
aic_df=aic_df[with(aic_df,order(species,lat)),]
write.csv(aic_df,here('outputs','ms_rmd','aic_table.csv'))


#Summarize results - BIC####
f_bic=list.files(here('outputs','TMB BIC'))
s_i_bic=as.numeric(gsub("[^0-9]", "", substr(f, 1, 3)))
m_ll=list()
full_mw=matrix(ncol=9,nrow=length(s_i_bic))
full_mw[,1]=s_i_bic

for(i in 1:length(f)){
  m_ll[[i]]<- read.csv(here('outputs','TMB BIC',f[i]))
  
  full_mw[i,2:9]=samEst::model_weights(m_ll[[i]][,2],form='AIC')
}
stock_info_filtered$stock.id2=seq(1:nrow(stock_info_filtered))

full_mw=full_mw[match(stock_info_filtered$stock.id2,full_mw[,1]),]
best_mod=apply(full_mw[,2:9],1,which.max) #this isn't working for some reason :{
sum_mod=summary(factor(best_mod))
sum_mod

bic_df=data.frame(stock_info_filtered[,3:8],round(full_mw[,2:9],2),best_mod)
bic_df=bic_df[with(bic_df,order(species,lat)),]
write.csv(bic_df,here('outputs','ms_rmd','bic_table.csv'))


#Pseudo-BMA+
#only evaluate likelihood for one version of model 8:
if(wm8<=6){ll8=ll8_1}else{ll8=ll8_2;wm8=wm8-6}

##Explorations....
#Revised model weights
model_weights(mod_lfo[[1]])


#Mean (log lik vs. sum)
m_ll1=apply(pw_lfo[[1]],1,mean)
pw_ll=pw_lfo[[5]]
hist(model_weights(pw_lfo[[1]],type='full'))
hist(model_weights(pw_ll,type='d90'))
hist(model_weights(pw_ll,type='d80'))


mod_weights=do.call(rbind.data.frame,lapply(pw_lfo,model_weights))
names(mod_weights)=paste('mod',seq(1:8),sep='_')
mod_weights$top_model=apply(mod_weights,1,which.max)
summary(factor(mod_weights$top_model))
summary(mod_weights)

mod_weights_d90=do.call(rbind.data.frame,lapply(pw_lfo,model_weights,type='d90'))
names(mod_weights_d90)=paste('mod',seq(1:8),sep='_')
mod_weights_d90$top_model=apply(mod_weights_d90,1,which.max)
summary(factor(mod_weights_d90$top_model))

mod_weights_d80=do.call(rbind.data.frame,lapply(pw_lfo,model_weights,type='d80'))
names(mod_weights_d80)=paste('mod',seq(1:8),sep='_')
mod_weights_d80$top_model=apply(mod_weights_d80,1,which.max)
summary(factor(mod_weights_d80$top_model))

top_model=apply(modelweight_summary,1,which.max)
top.model3=do.call(rbind.data.frame, top_model)
top.model3[90:nrow(modelweight_summary),1]=NA
modelweight_summary$top.model=top.model3[,1]

alpha_vary_w=mod_weights[,3]+mod_weights[,6]
summary(alpha_vary_w)
beta_vary_w=mod_weights[,4]+mod_weights[,7]
summary(beta_vary_w)
static_w=mod_weights[,1]+mod_weights[,2]
summary(static_w)
both_w=mod_weights[,5]+mod_weights[,8]
summary(both_w)



###Sensitivity analysis & minimum data retained ###
#harvest history
#Fraser sockeye - really high ER, to really low ER when alpha changed
#Exploitation rates - plot them with ER/U_msy (time-invariant U_msy)
#Use the simulation to get stock statistics that may be an indicator for genuine change
#Summary data frame and store for pointwise LLs
stock_info_filtered_25=subset(stock_info,n.years>=25) #242 stocks

stock_dat3=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat3$stock.id)) #242

stock_dat3$logR_S=log(stock_dat3$recruits/stock_dat3$spawners)


lfo_summary_L10=data.frame(stock=stock_info_filtered$stock.name,LL_m1=NA,LL_m2.1=NA,LL_m2.3=NA,LL_m2.5=NA,LL_m3.1=NA,LL_m3.3=NA,
                          LL_m3.5=NA,LL_m4.1=NA,LL_m4.3=NA,LL_m4.5=NA,LL_m5.1=NA,LL_m5.3=NA,LL_m5.5=NA,LL_m6.1=NA,LL_m6.3=NA,LL_m6.5=NA,
                          LL_m6.1w=NA,LL_m6.3w=NA,LL_m6.5w=NA,LL_m7.1=NA,LL_m7.3=NA,LL_m7.5=NA,LL_m7.1w=NA,LL_m7.3w=NA,LL_m7.5w=NA,
                          LL_m8.1.1=NA,LL_m8.1.3=NA,LL_m8.1.5=NA,LL_m8.1.1w=NA,LL_m8.1.3w=NA,LL_m8.1.5w=NA,LL_m8.2.1=NA,LL_m8.2.3=NA,LL_m8.2.5=NA,LL_m8.2.1w=NA,LL_m8.2.3w=NA,LL_m8.2.5w=NA)
lfo_summary_L15=data.frame(stock=stock_info_filtered$stock.name,LL_m1=NA,LL_m2.1=NA,LL_m2.3=NA,LL_m2.5=NA,LL_m3.1=NA,LL_m3.3=NA,
                              LL_m3.5=NA,LL_m4.1=NA,LL_m4.3=NA,LL_m4.5=NA,LL_m5.1=NA,LL_m5.3=NA,LL_m5.5=NA,LL_m6.1=NA,LL_m6.3=NA,LL_m6.5=NA,
                              LL_m6.1w=NA,LL_m6.3w=NA,LL_m6.5w=NA,LL_m7.1=NA,LL_m7.3=NA,LL_m7.5=NA,LL_m7.1w=NA,LL_m7.3w=NA,LL_m7.5w=NA,
                              LL_m8.1.1=NA,LL_m8.1.3=NA,LL_m8.1.5=NA,LL_m8.1.1w=NA,LL_m8.1.3w=NA,LL_m8.1.5w=NA,LL_m8.2.1=NA,LL_m8.2.3=NA,LL_m8.2.5=NA,LL_m8.2.1w=NA,LL_m8.2.3w=NA,LL_m8.2.5w=NA)
lfo_summary_L20=data.frame(stock=stock_info_filtered$stock.name,LL_m1=NA,LL_m2.1=NA,LL_m2.3=NA,LL_m2.5=NA,LL_m3.1=NA,LL_m3.3=NA,
                              LL_m3.5=NA,LL_m4.1=NA,LL_m4.3=NA,LL_m4.5=NA,LL_m5.1=NA,LL_m5.3=NA,LL_m5.5=NA,LL_m6.1=NA,LL_m6.3=NA,LL_m6.5=NA,
                              LL_m6.1w=NA,LL_m6.3w=NA,LL_m6.5w=NA,LL_m7.1=NA,LL_m7.3=NA,LL_m7.5=NA,LL_m7.1w=NA,LL_m7.3w=NA,LL_m7.5w=NA,
                              LL_m8.1.1=NA,LL_m8.1.3=NA,LL_m8.1.5=NA,LL_m8.1.1w=NA,LL_m8.1.3w=NA,LL_m8.1.5w=NA,LL_m8.2.1=NA,LL_m8.2.3=NA,LL_m8.2.5=NA,LL_m8.2.1w=NA,LL_m8.2.3w=NA,LL_m8.2.5w=NA)

L10_lfo=list() #pointwise lfoelihood
L15_lfo=list() #pointwise lfoelihood
L20_lfo=list() #pointwise lfoelihood
for(i in 1:5){
  s<- subset(stock_dat3,stock.id==stock_info_filtered_25$stock.id[i])
  s<- s[complete.cases(s$spawners),]
  
  #Assess model fits for each model type
  #model 1 - static Ricker
  ll1_10<- stan_lfo_cv(mod=m1,type='static',df=s,L=10)
  ll1_15<- stan_lfo_cv(mod=m1,type='static',df=s,L=15)
  ll1_20<- stan_lfo_cv(mod=m1,type='static',df=s,L=20)
  
  #model 2 - static autocorrelated Ricker
  ll2_10<- stan_lfo_cv(mod=m2,type='tv',df=s,L=10)
  ll2_15<- stan_lfo_cv(mod=m2,type='tv',df=s,L=15)
  ll2_20<- stan_lfo_cv(mod=m2,type='tv',df=s,L=20)
  
  #model 3 - dynamic productivity Ricker
  ll3_10<- stan_lfo_cv(mod=m3,type='tv',df=s,L=10)
  ll3_15<- stan_lfo_cv(mod=m3,type='tv',df=s,L=15)
  ll3_20<- stan_lfo_cv(mod=m3,type='tv',df=s,L=20)
  
  #model 4 - dynamic capacity Ricker
  ll4_10<- stan_lfo_cv(mod=m4,type='tv',df=s,L=10)
  ll4_15<- stan_lfo_cv(mod=m4,type='tv',df=s,L=15)
  ll4_20<- stan_lfo_cv(mod=m4,type='tv',df=s,L=20)
  
  #model 5 - dynamic productivity & capacity Ricker
  ll5_10<- stan_lfo_cv(mod=m5,type='tv',df=s,L=10)
  ll5_15<- stan_lfo_cv(mod=m5,type='tv',df=s,L=15)
  ll5_20<- stan_lfo_cv(mod=m5,type='tv',df=s,L=20)
  
  #model 6 - productivity regime shift - 2 regimes
  ll6_10<- stan_lfo_cv(mod=m6,type='regime',df=s,L=10,K=2)
  ll6_15<- stan_lfo_cv(mod=m6,type='regime',df=s,L=15,K=2)
  ll6_20<- stan_lfo_cv(mod=m6,type='regime',df=s,L=20,K=2)
  
  #model 7 - capacity regime shift
  ll7_10<- stan_lfo_cv(mod=m7,type='regime',df=s,L=10,K=2)
  ll7_15<- stan_lfo_cv(mod=m7,type='regime',df=s,L=15,K=2)
  ll7_20<- stan_lfo_cv(mod=m7,type='regime',df=s,L=20,K=2)
  
  #model 8 - productivity and capacity regime shift
  ll8_1_10<- stan_lfo_cv(mod=m8_1,type='regime',df=s,L=10,K=2)
  ll8_1_15<- stan_lfo_cv(mod=m8_1,type='regime',df=s,L=15,K=2)
  ll8_1_20<- stan_lfo_cv(mod=m8_1,type='regime',df=s,L=20,K=2)
  
  #model 8 - productivity and capacity regime shift
  ll8_2_10<- stan_lfo_cv(mod=m8_2,type='regime',df=s,L=10,K=2)
  ll8_2_15<- stan_lfo_cv(mod=m8_2,type='regime',df=s,L=15,K=2)
  ll8_2_20<- stan_lfo_cv(mod=m8_2,type='regime',df=s,L=20,K=2)
  
  #Take the sum of the estimated pointwise likelihood estimates (=elpd_loo)
  L10_lfos[[i]]=rbind(ll1_10,ll2_10,ll3_10,ll4_10,ll5_10,ll6_10,ll7_10,ll8_1_10,ll8_2_10) #pointwise lfoelihood
  L15_lfo[[i]]=rbind(ll1_15,ll2_15,ll3_15,ll4_15,ll5_15,ll6_15,ll7_15,ll8_1_15,ll8_2_15) #pointwise lfoelihood
  L20_lfo[[i]]=rbind(ll1_20,ll2_20,ll3_20,ll4_20,ll5_20,ll6_20,ll7_20,ll8_1_20,ll8_2_20) #pointwise lfoelihood
}

wm_L10=list()
which.max(apply(L10_lfo[[1]][2:4,],1,sum))
model_weights(L10_lfo[[1]][2:4,])
which.max(apply(L10_lfo[[1]][5:7,],1,sum))
model_weights(L10_lfo[[1]][5:7,])
which.max(apply(L10_lfo[[1]][8:10,],1,sum))

which.max(apply(L10_lfo[[1]][11:13,],1,sum))
which.max(apply(L10_lfo[[1]][14:19,],1,sum))
which.max(apply(L10_lfo[[1]][20:25,],1,sum))
which.max(apply(L10_lfo[[1]][26:31,],1,sum))
which.max(apply(L10_lfo[[1]][32:37,],1,sum))

##comp of model weights
if(sum(L10_lfo[[1]][26,])>sum(L10_lfo[[1]][32,])){
  L10_comp1=rbind(L10_lfo[[1]][1,],L10_lfo[[1]][2,],L10_lfo[[1]][5,],L10_lfo[[1]][8,],L10_lfo[[1]][11,],L10_lfo[[1]][14,],L10_lfo[[1]][20,],L10_lfo[[1]][26,])
}else{
  L10_comp1=rbind(L10_lfo[[1]][1,],L10_lfo[[1]][2,],L10_lfo[[1]][5,],L10_lfo[[1]][8,],L10_lfo[[1]][11,],L10_lfo[[1]][14,],L10_lfo[[1]][20,],L10_lfo[[1]][32,])    
  }
model_weights(L10_comp1)

if(sum(L15_lfo[[1]][26,])>sum(L15_lfo[[1]][32,])){
  L15_comp1=rbind(L15_lfo[[1]][1,],L15_lfo[[1]][2,],L15_lfo[[1]][5,],L15_lfo[[1]][8,],L15_lfo[[1]][11,],L15_lfo[[1]][14,],L15_lfo[[1]][20,],L15_lfo[[1]][26,])
}else{
  L15_comp1=rbind(L15_lfo[[1]][1,],L15_lfo[[1]][2,],L15_lfo[[1]][5,],L15_lfo[[1]][8,],L15_lfo[[1]][11,],L15_lfo[[1]][14,],L15_lfo[[1]][20,],L15_lfo[[1]][32,])    
}
model_weights(L15_comp1)

if(sum(L20_lfo[[2]][26,])>sum(L20_lfo[[2]][32,])){
  L20_comp1=rbind(L20_lfo[[2]][1,],L20_lfo[[2]][2,],L20_lfo[[2]][5,],L20_lfo[[2]][8,],L20_lfo[[2]][11,],L20_lfo[[2]][14,],L20_lfo[[2]][20,],L20_lfo[[2]][26,])
}else{
  L20_comp1=rbind(L20_lfo[[2]][1,],L20_lfo[[2]][2,],L20_lfo[[2]][5,],L20_lfo[[2]][8,],L20_lfo[[2]][11,],L20_lfo[[2]][14,],L20_lfo[[2]][20,],L20_lfo[[2]][32,])    
}
model_weights(L20_comp1)

if(sum(L10_lfo[[2]][26,])>sum(L10_lfo[[2]][32,])){
  L10_comp1=rbind(L10_lfo[[2]][1,],L10_lfo[[2]][2,],L10_lfo[[2]][5,],L10_lfo[[2]][8,],L10_lfo[[2]][11,],L10_lfo[[2]][14,],L10_lfo[[2]][20,],L10_lfo[[2]][26,])
}else{
  L10_comp1=rbind(L10_lfo[[2]][1,],L10_lfo[[2]][2,],L10_lfo[[2]][5,],L10_lfo[[2]][8,],L10_lfo[[2]][11,],L10_lfo[[2]][14,],L10_lfo[[2]][20,],L10_lfo[[2]][32,])    
}
model_weights(L10_comp1)

if(sum(L15_lfo[[2]][26,])>sum(L15_lfo[[2]][32,])){
  L15_comp1=rbind(L15_lfo[[2]][1,],L15_lfo[[2]][2,],L15_lfo[[2]][5,],L15_lfo[[2]][8,],L15_lfo[[2]][11,],L15_lfo[[2]][14,],L15_lfo[[2]][20,],L15_lfo[[2]][26,])
}else{
  L15_comp1=rbind(L15_lfo[[2]][1,],L15_lfo[[2]][2,],L15_lfo[[2]][5,],L15_lfo[[2]][8,],L15_lfo[[2]][11,],L15_lfo[[2]][14,],L15_lfo[[2]][20,],L15_lfo[[2]][32,])    
}
model_weights(L15_comp1)

if(sum(L20_lfo[[2]][26,])>sum(L20_lfo[[2]][32,])){
  L20_comp1=rbind(L20_lfo[[2]][1,],L20_lfo[[2]][2,],L20_lfo[[2]][5,],L20_lfo[[2]][8,],L20_lfo[[2]][11,],L20_lfo[[2]][14,],L20_lfo[[2]][20,],L20_lfo[[2]][26,])
}else{
  L20_comp1=rbind(L20_lfo[[2]][1,],L20_lfo[[2]][2,],L20_lfo[[2]][5,],L20_lfo[[2]][8,],L20_lfo[[2]][11,],L20_lfo[[2]][14,],L20_lfo[[2]][20,],L20_lfo[[2]][32,])    
}
model_weights(L20_comp1)

if(sum(L10_lfo[[3]][26,])>sum(L10_lfo[[3]][32,])){
  L10_comp1=rbind(L10_lfo[[3]][1,],L10_lfo[[3]][2,],L10_lfo[[3]][5,],L10_lfo[[3]][8,],L10_lfo[[3]][11,],L10_lfo[[3]][14,],L10_lfo[[3]][20,],L10_lfo[[3]][26,])
}else{
  L10_comp1=rbind(L10_lfo[[3]][1,],L10_lfo[[3]][2,],L10_lfo[[3]][5,],L10_lfo[[3]][8,],L10_lfo[[3]][11,],L10_lfo[[3]][14,],L10_lfo[[3]][20,],L10_lfo[[3]][32,])    
}
model_weights(L10_comp1)

if(sum(L15_lfo[[3]][26,])>sum(L15_lfo[[3]][32,])){
  L15_comp1=rbind(L15_lfo[[3]][1,],L15_lfo[[3]][2,],L15_lfo[[3]][5,],L15_lfo[[3]][8,],L15_lfo[[3]][11,],L15_lfo[[3]][14,],L15_lfo[[3]][20,],L15_lfo[[3]][26,])
}else{
  L15_comp1=rbind(L15_lfo[[3]][1,],L15_lfo[[3]][2,],L15_lfo[[3]][5,],L15_lfo[[3]][8,],L15_lfo[[3]][11,],L15_lfo[[3]][14,],L15_lfo[[3]][20,],L15_lfo[[3]][32,])    
}
model_weights(L15_comp1)

if(sum(L20_lfo[[3]][26,])>sum(L20_lfo[[3]][32,])){
  L20_comp1=rbind(L20_lfo[[3]][1,],L20_lfo[[3]][2,],L20_lfo[[3]][5,],L20_lfo[[3]][8,],L20_lfo[[3]][11,],L20_lfo[[3]][14,],L20_lfo[[3]][20,],L20_lfo[[3]][26,])
}else{
  L20_comp1=rbind(L20_lfo[[3]][1,],L20_lfo[[3]][2,],L20_lfo[[3]][5,],L20_lfo[[3]][8,],L20_lfo[[3]][11,],L20_lfo[[3]][14,],L20_lfo[[3]][20,],L20_lfo[[3]][32,])    
}
model_weights(L20_comp1)

if(sum(L10_lfo[[4]][26,])>sum(L10_lfo[[4]][42,])){
  L10_comp1=rbind(L10_lfo[[4]][1,],L10_lfo[[4]][2,],L10_lfo[[4]][5,],L10_lfo[[4]][8,],L10_lfo[[4]][11,],L10_lfo[[4]][14,],L10_lfo[[4]][20,],L10_lfo[[4]][26,])
}else{
  L10_comp1=rbind(L10_lfo[[4]][1,],L10_lfo[[4]][2,],L10_lfo[[4]][5,],L10_lfo[[4]][8,],L10_lfo[[4]][11,],L10_lfo[[4]][14,],L10_lfo[[4]][20,],L10_lfo[[4]][42,])    
}
model_weights(L10_comp1)

if(sum(L15_lfo[[4]][26,])>sum(L15_lfo[[4]][42,])){
  L15_comp1=rbind(L15_lfo[[4]][1,],L15_lfo[[4]][2,],L15_lfo[[4]][5,],L15_lfo[[4]][8,],L15_lfo[[4]][11,],L15_lfo[[4]][14,],L15_lfo[[4]][20,],L15_lfo[[4]][26,])
}else{
  L15_comp1=rbind(L15_lfo[[4]][1,],L15_lfo[[4]][2,],L15_lfo[[4]][5,],L15_lfo[[4]][8,],L15_lfo[[4]][11,],L15_lfo[[4]][14,],L15_lfo[[4]][20,],L15_lfo[[4]][42,])    
}
model_weights(L15_comp1)

if(sum(L20_lfo[[4]][26,])>sum(L20_lfo[[4]][42,])){
  L20_comp1=rbind(L20_lfo[[4]][1,],L20_lfo[[4]][2,],L20_lfo[[4]][5,],L20_lfo[[4]][8,],L20_lfo[[4]][11,],L20_lfo[[4]][14,],L20_lfo[[4]][20,],L20_lfo[[4]][26,])
}else{
  L20_comp1=rbind(L20_lfo[[4]][1,],L20_lfo[[4]][2,],L20_lfo[[4]][5,],L20_lfo[[4]][8,],L20_lfo[[4]][11,],L20_lfo[[4]][14,],L20_lfo[[4]][20,],L20_lfo[[4]][42,])    
}
model_weights(L20_comp1)

which.max(apply(L15_lfo[[1]][2:4,],1,sum))
which.max(apply(L15_lfo[[1]][5:7,],1,sum))
which.max(apply(L15_lfo[[1]][8:10,],1,sum))
model_weights(L15_lfo[[1]][8:10,])
which.max(apply(L15_lfo[[1]][11:13,],1,sum))
which.max(apply(L15_lfo[[1]][14:19,],1,sum))
which.max(apply(L15_lfo[[1]][20:25,],1,sum))
which.max(apply(L15_lfo[[1]][26:31,],1,sum))
which.max(apply(L15_lfo[[1]][32:37,],1,sum))

which.max(apply(L20_lfo[[1]][2:4,],1,sum))
which.max(apply(L20_lfo[[1]][5:7,],1,sum))
which.max(apply(L20_lfo[[1]][8:10,],1,sum))
which.max(apply(L20_lfo[[1]][11:13,],1,sum))
which.max(apply(L20_lfo[[1]][14:19,],1,sum))
which.max(apply(L20_lfo[[1]][20:25,],1,sum))
which.max(apply(L20_lfo[[1]][26:31,],1,sum))
which.max(apply(L20_lfo[[1]][32:37,],1,sum))


which.max(apply(L10_lfo[[2]][2:4,],1,sum))
which.max(apply(L10_lfo[[2]][5:7,],1,sum))
which.max(apply(L10_lfo[[2]][8:10,],1,sum))
which.max(apply(L10_lfo[[2]][11:13,],1,sum))
which.max(apply(L10_lfo[[2]][14:19,],1,sum))
which.max(apply(L10_lfo[[2]][20:25,],1,sum))
which.max(apply(L10_lfo[[2]][26:31,],1,sum))
which.max(apply(L10_lfo[[2]][32:37,],1,sum))

which.max(apply(L15_lfo[[2]][2:4,],1,sum))
which.max(apply(L15_lfo[[2]][5:7,],1,sum))
which.max(apply(L15_lfo[[2]][8:10,],1,sum))
which.max(apply(L15_lfo[[2]][11:13,],1,sum))
which.max(apply(L15_lfo[[2]][14:19,],1,sum))
which.max(apply(L15_lfo[[2]][20:25,],1,sum))
which.max(apply(L15_lfo[[2]][26:31,],1,sum))
which.max(apply(L15_lfo[[2]][32:37,],1,sum))

which.max(apply(L20_lfo[[2]][2:4,],1,sum))
which.max(apply(L20_lfo[[2]][5:7,],1,sum))
which.max(apply(L20_lfo[[2]][8:10,],1,sum))
which.max(apply(L20_lfo[[2]][11:13,],1,sum))
which.max(apply(L20_lfo[[2]][14:19,],1,sum))
which.max(apply(L20_lfo[[2]][20:25,],1,sum))
which.max(apply(L20_lfo[[2]][26:31,],1,sum))
which.max(apply(L20_lfo[[2]][32:37,],1,sum))

which.max(apply(L10_lfo[[3]][2:4,],1,sum))
which.max(apply(L10_lfo[[3]][5:7,],1,sum))
which.max(apply(L10_lfo[[3]][8:10,],1,sum))
which.max(apply(L10_lfo[[3]][11:13,],1,sum))
which.max(apply(L10_lfo[[3]][14:19,],1,sum))
which.max(apply(L10_lfo[[3]][20:25,],1,sum))
which.max(apply(L10_lfo[[3]][26:31,],1,sum))
which.max(apply(L10_lfo[[3]][32:37,],1,sum))

which.max(apply(L10_lfo[[4]][2:4,],1,sum))
which.max(apply(L10_lfo[[4]][5:7,],1,sum))
which.max(apply(L10_lfo[[4]][8:10,],1,sum))
which.max(apply(L10_lfo[[4]][11:13,],1,sum))
which.max(apply(L10_lfo[[4]][14:19,],1,sum))
which.max(apply(L10_lfo[[4]][20:25,],1,sum))
which.max(apply(L10_lfo[[4]][26:31,],1,sum))
which.max(apply(L10_lfo[[4]][32:37,],1,sum))


ggplot(sp_sum, aes(fill=species, y=n, x=region))+
  scale_fill_manual(values=c("darkgray", "darkgreen", "darkblue","darksalmon","darkred"))+ 
  geom_bar(position="stack", stat="identity") + theme_minimal() +
  xlab('')+ylab('No. time-series')

ggplot(stock_info_filtered, aes(x=n.years, color=species)) +
  geom_histogram()



##Cowichan chinook plots####
sr_dat=subset(stock_dat2,stock.id==49)
sr_dat$logR

plot(recruits~spawners,data=sr_dat,ylab='Recruits',xlab='Spawners',xlim=c(0,max(sr_dat$spawners)))


#TMB runs - LFO/AIC/BIC####
AICdf=matrix(nrow=nrow(stock_info_filtered),ncol=3)
BICdf=matrix(nrow=nrow(stock_info_filtered),ncol=3)
for(u in 1:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  TMBstatic <- samEst::ricker_TMB(data=df)
  TMBac <- samEst::ricker_TMB(data=df, AC=TRUE)
  TMBtva <- tryCatch(samEst::ricker_rw_TMB(data=df,tv.par='a'),error = function(e) {TMBtva=list(conv_problem=TRUE)})
  
  AICdf[u,]<-c(TMBstatic$AICc,
               TMBac$AICc,
               TMBtva$AICc)
  BICdf[u,]<-c(TMBstatic$BIC,
               TMBac$BIC,
               TMBtva$BIC)
}
AICdf=AICdf[complete.cases(AICdf),]
write.csv(AICdf,'AICdf.csv')
BICdf=BICdf[complete.cases(BICdf),]
write.csv(BICdf,'BICdf.csv')


aic=read.csv(here('AICdf.csv'))
bic=read.csv(here('BICdf.csv'))

stock_inf2=stock_info_filtered[-134,]
stock_inf2$state2=stock_inf2$state
stock_inf2=stock_inf2 %>% mutate(state2=ifelse(state=='OR','OR-WA',state2),
                                 state2=ifelse(state=='WA','OR-WA',state2))

w1=matrix(nrow=nrow(stock_inf2),ncol=3)
w2=matrix(nrow=nrow(stock_inf2),ncol=3)
for(i in 1:nrow(AICdf)){
  w1[i,]=samEst::model_weights(AICdf[i,],form='AIC')
  w2[i,]=samEst::model_weights(BICdf[i,],form='AIC')
}
apply(w1,1,which.max)
apply(w2,1,which.max)
stock_inf2$bm.aic=apply(w1,1,which.max)
stock_inf2$bm.aic=ifelse(stock_inf2$bm.aic==1,'stationary',stock_inf2$bm.aic)
stock_inf2$bm.aic=ifelse(stock_inf2$bm.aic==2,'autocorr',stock_inf2$bm.aic)
stock_inf2$bm.aic=ifelse(stock_inf2$bm.aic==3,'dynamic',stock_inf2$bm.aic)
stock_inf2$bm.aic=factor(stock_inf2$bm.aic,levels=c('stationary','autocorr','dynamic'))
stock_inf2$bm.bic=apply(w2,1,which.max)
stock_inf2$bm.bic=ifelse(stock_inf2$bm.bic==1,'stationary',stock_inf2$bm.bic)
stock_inf2$bm.bic=ifelse(stock_inf2$bm.bic==2,'autocorr',stock_inf2$bm.bic)
stock_inf2$bm.bic=ifelse(stock_inf2$bm.bic==3,'dynamic',stock_inf2$bm.bic)
stock_inf2$bm.bic=factor(stock_inf2$bm.bic,levels=c('stationary','autocorr','dynamic'))

sp = stock_inf2 %>% group_by(species,bm.aic) %>% summarize(n=n(),.groups='keep')
sp1 = stock_inf2 %>% group_by(species) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$species,sp1$species)]

mod_col=cbind(c("azure4", "cyan4","deepskyblue4"),mod=c('stationary','autocorr','dynamic'))

p=ggplot(sp, aes(fill=factor(bm.aic), y=prop, x=species))+
  scale_fill_manual(values=mod_col[match(levels(factor(sp$bm.aic)),mod_col[,2])],name='Model')+ 
  ggtitle("Top-ranked model (AIC)")+
  geom_bar(position="stack", stat="identity",color='white') +
  geom_text(label=sp$n, position = position_stack(vjust = 0.5),colour='white',size=2.2)+
  ylim(0,1)+
  theme_minimal() +
  xlab('')+
  ylab('Prop. populations')+ 
  theme(text=element_text(size=14),axis.text=element_text(size=14),axis.line = element_line(colour = "black", 
                                                                                            size = 1, linetype = "solid"),plot.title = element_text(size=14))+
  annotation_custom(grid::textGrob(sp1$n[1]), xmin = 1, xmax = 1, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[2]), xmin = 2, xmax = 2, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[3]), xmin = 3, xmax = 3, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[4]), xmin = 4, xmax = 4, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)
pdf(file='aic.pdf',height=6,width=8)
p
dev.off()


sp = stock_inf2 %>% group_by(species,bm.bic) %>% summarize(n=n(),.groups='keep')
sp1 = stock_inf2 %>% group_by(species) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$species,sp1$species)]

b=ggplot(sp, aes(fill=factor(bm.bic), y=prop, x=species))+
  scale_fill_manual(values=mod_col[match(levels(factor(sp$bm.bic)),mod_col[,2])],name='Model')+ 
  ggtitle("Top-ranked model (BIC)")+
  geom_bar(position="stack", stat="identity",color='white') +
  geom_text(label=sp$n, position = position_stack(vjust = 0.5),colour='white',size=2.2)+
  ylim(0,1)+
  theme_minimal() +
  xlab('')+
  ylab('Prop. populations')+ 
  theme(text=element_text(size=14),axis.text=element_text(size=14),axis.line = element_line(colour = "black", 
                                                                                            size = 1, linetype = "solid"),plot.title = element_text(size=14))+
  annotation_custom(grid::textGrob(sp1$n[1]), xmin = 1, xmax = 1, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[2]), xmin = 2, xmax = 2, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[3]), xmin = 3, xmax = 3, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[4]), xmin = 4, xmax = 4, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

pdf(file='bic.pdf',height=6,width=8)
b
dev.off()




p=ggplot(sp, aes(fill=factor(bm.aic), y=prop, x=species))+
  scale_fill_manual(values=mod_col[match(levels(factor(sp$bm.aic)),mod_col[,2])],name='Model')+ 
  ggtitle("Top-ranked model (AIC)")+
  geom_bar(position="stack", stat="identity",color='white') +
  geom_text(label=sp$n, position = position_stack(vjust = 0.5),colour='white',size=2.2)+
  ylim(0,1)+
  theme_minimal() +
  xlab('')+
  ylab('Prop. populations')+ 
  theme(text=element_text(size=14),axis.text=element_text(size=14),axis.line = element_line(colour = "black", 
                                                                                            size = 1, linetype = "solid"),plot.title = element_text(size=14))+
  annotation_custom(grid::textGrob(sp1$n[1]), xmin = 1, xmax = 1, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[2]), xmin = 2, xmax = 2, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[3]), xmin = 3, xmax = 3, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[4]), xmin = 4, xmax = 4, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)
pdf(file='test.pdf',height=6,width=8)
p
dev.off()
