rm(list=ls())
library(here);library(dplyr);library(rstan)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_aug2022.csv'))
stock_info<- read.csv(here('data','filtered datasets','all_stocks_info_aug2022.csv'))

source(here('code','samEst code','stan_functions.R'))
source(here('code','samEst code','lfo_stan_functions.R'))

#Remove stocks with less than 15 years of recruitment data
stock_info_filtered=subset(stock_info,n.years>=18) #242 stocks

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock.id)) #242

stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)

#Define models (helps prevent crashing)
m1=sr_mod(type='static',ac = FALSE,par='n',loglik=T)
m2=sr_mod(type='static',ac = TRUE,par='n',loglik=T)
m3=sr_mod(type='rw',par='a',loglik=T)
m4=sr_mod(type='rw',par='b',loglik=T)
m5=sr_mod(type='rw',par='both',loglik=T)
m5f=sr_mod(type='rw',par='both',loglik=F)
m6=sr_mod(type='hmm',par='a',loglik=T)
m6f=sr_mod(type='hmm',par='a',loglik=F)
m7=sr_mod(type='hmm',par='b',loglik=T)
m7f=sr_mod(type='hmm',par='b',loglik=F)
m8_1=sr_mod(type='hmm',par='both',loglik=T,caphigh=F)
m8_2=sr_mod(type='hmm',par='both',loglik=T,caphigh=T)
m8_1f=sr_mod(type='hmm',par='both',loglik=F,caphigh=F)
m8_2f=sr_mod(type='hmm',par='both',loglik=F,caphigh=T)

#Summary data frame and store for pointwise LLs
loglik_summary=data.frame(stock=stock_info_filtered$stock.name,LL_m1=NA,LL_m2.1=NA,LL_m2.3=NA,LL_m2.5=NA,LL_m3.1=NA,LL_m3.3=NA,
                          LL_m3.5=NA,LL_m4.1=NA,LL_m4.3=NA,LL_m4.5=NA,LL_m5.1=NA,LL_m5.3=NA,LL_m5.5=NA,LL_m6.1=NA,LL_m6.3=NA,LL_m6.5=NA,
                          LL_m6.1w=NA,LL_m6.3w=NA,LL_m6.5w=NA,LL_m7.1=NA,LL_m7.3=NA,LL_m7.5=NA,LL_m7.1w=NA,LL_m7.3w=NA,LL_m7.5w=NA,
                          LL_m8.1=NA,LL_m8.3=NA,LL_m8.5=NA,LL_m8.1w=NA,LL_m8.3w=NA,LL_m8.5w=NA)
pw_loglik=list() #pointwise loglikelihood
se_elpd_loo=list()
modelweight_summary=data.frame(stock=stock_info_filtered$stock.name,w_m1=NA,w_m2=NA,w_m3=NA,w_m4=NA,w_m5=NA,w_m6=NA,
                               w_m7=NA,w_m8=NA)
for(i in 90:100){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  s<- s[complete.cases(s$spawners),]
  
  #Assess model fits for each model type
  #model 1 - static Ricker
  ll1<- stan_lfo_cv(mod=m1,type='static',df=s,L=10)
  #model 2 - static autocorrelated Ricker
  ll2<- stan_lfo_cv(mod=m2,type='tv',df=s,L=10)
  #model 3 - dynamic productivity Ricker
  ll3<- stan_lfo_cv(mod=m3,type='tv',df=s,L=10)
  
  #model 4 - dynamic capacity Ricker
  ll4<- stan_lfo_cv(mod=m4,type='tv',df=s,L=10)
  
  #model 5 - dynamic productivity & capacity Ricker
  ll5<- stan_lfo_cv(mod=m5,type='tv',df=s,L=10)
  
  #model 6 - productivity regime shift - 2 regimes
  ll6<- stan_lfo_cv(mod=m6,type='regime',df=s,L=10,K=2)
  
  #model 7 - capacity regime shift
  ll7<- stan_lfo_cv(mod=m7,type='regime',df=s,L=10,K=2)
  
  #model 8 - productivity and capacity regime shift
  ll8_1<- stan_lfo_cv(mod=m8_1,type='regime',df=s,L=10,K=2)
  #model 8 - productivity and capacity regime shift
  ll8_2<- stan_lfo_cv(mod=m8_2,type='regime',df=s,L=10,K=2)
  
  #Take the sum of the estimated pointwise likelihood estimates (=elpd_loo)
  loglik_summary[i,2]=sum(ll1)
  loglik_summary[i,3]=sum(ll2[[1]])
  loglik_summary[i,4]=sum(ll2[[2]])
  loglik_summary[i,5]=sum(ll2[[3]])
  loglik_summary[i,6]=sum(ll3[[1]])
  loglik_summary[i,7]=sum(ll3[[2]])
  loglik_summary[i,8]=sum(ll3[[3]])
  loglik_summary[i,9]=sum(ll4[[1]])
  loglik_summary[i,10]=sum(ll4[[2]])
  loglik_summary[i,11]=sum(ll4[[3]])
  loglik_summary[i,12]=sum(ll5[[1]])
  loglik_summary[i,13]=sum(ll5[[2]])
  loglik_summary[i,14]=sum(ll5[[3]])
  loglik_summary[i,15]=sum(ll6[[1]])
  loglik_summary[i,16]=sum(ll6[[2]])
  loglik_summary[i,17]=sum(ll6[[3]])
  loglik_summary[i,18]=sum(ll6[[4]])
  loglik_summary[i,19]=sum(ll6[[5]])
  loglik_summary[i,20]=sum(ll6[[6]])
  loglik_summary[i,21]=sum(ll7[[1]])
  loglik_summary[i,22]=sum(ll7[[2]])
  loglik_summary[i,23]=sum(ll7[[3]])
  loglik_summary[i,24]=sum(ll7[[4]])
  loglik_summary[i,25]=sum(ll7[[5]])
  loglik_summary[i,26]=sum(ll7[[6]])
  loglik_summary[i,27]=sum(ll8[[1]])
  loglik_summary[i,28]=sum(ll8[[2]])
  loglik_summary[i,29]=sum(ll8[[3]])
  loglik_summary[i,30]=sum(ll8[[4]])
  loglik_summary[i,31]=sum(ll8[[5]])
  loglik_summary[i,32]=sum(ll8[[6]])
  
  #Pseudo-BMA+
  wm2=which.max(loglik_summary[i,3:5]) #select best likelihood from different timeframes (1-y back, 3-y back, 5-y back)
  wm3=which.max(loglik_summary[i,6:8]) #best fit for model 2
  wm4=which.max(loglik_summary[i,9:11]) #best fit for model 2
  wm5=which.max(loglik_summary[i,12:14]) #best fit for model 2
  wm6=which.max(loglik_summary[i,15:20]) #best fit for model 2
  wm7=which.max(loglik_summary[i,21:26]) #best fit for model 2
  wm8=which.max(loglik_summary[i,27:31]) #best fit for model 2

  pw_loglik[[i]]=rbind(ll1,ll2[[wm2]],ll3[[wm3]],ll4[[wm4]],ll5[[wm5]],ll6[[wm6]],ll7[[wm7]],ll8[[wm8]])
  
  modelweight_summary[i,2:ncol(modelweight_summary)]= pseudobma_weights(pw_loglik[[i]])
}

##Explorations....

#Mean (log lik vs. sum)
m_ll1=apply(pw_loglik[[1]],1,mean)
pw_ll=pw_loglik[[5]]
hist(model_weights(pw_loglik[[1]],type='full'))
hist(model_weights(pw_ll,type='d90'))
hist(model_weights(pw_ll,type='d80'))

mod_weights=do.call(rbind.data.frame,lapply(pw_loglik,model_weights))
names(mod_weights)=paste('mod',seq(1:8),sep='_')
mod_weights$top_model=apply(mod_weights,1,which.max)
summary(factor(mod_weights$top_model))

mod_weights_d90=do.call(rbind.data.frame,lapply(pw_loglik,model_weights,type='d90'))
names(mod_weights_d90)=paste('mod',seq(1:8),sep='_')
mod_weights_d90$top_model=apply(mod_weights_d90,1,which.max)
summary(factor(mod_weights_d90$top_model))

mod_weights_d80=do.call(rbind.data.frame,lapply(pw_loglik,model_weights,type='d80'))
names(mod_weights_d80)=paste('mod',seq(1:8),sep='_')
mod_weights_d80$top_model=apply(mod_weights_d80,1,which.max)
summary(factor(mod_weights_d80$top_model))

top_model=apply(modelweight_summary,1,which.max)
top.model3=do.call(rbind.data.frame, top_model)
top.model3[90:nrow(modelweight_summary),1]=NA
modelweight_summary$top.model=top.model3[,1]

summary(factor(modelweight_summary$top.model))
summary(modelweight_summary[,2:9])
alpha_vary_w=modelweight_summary[,4]+modelweight_summary[,7]
summary(alpha_vary_w)
beta_vary_w=modelweight_summary[,5]+modelweight_summary[,8]
summary(beta_vary_w)
static_w=modelweight_summary[,2]+modelweight_summary[,3]
summary(static_w)
both_w=modelweight_summary[,6]+modelweight_summary[,9]
summary(both_w)
