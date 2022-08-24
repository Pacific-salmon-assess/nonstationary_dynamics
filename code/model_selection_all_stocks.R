rm(list=ls())
library(here);library(dplyr);library(rstan)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_aug2022.csv'))
stock_info<- read.csv(here('data','filtered datasets','all_stocks_info_aug2022.csv'))

source(here('code','functions.R'))
source(here('code','lfo_functions.R'))


#Remove stocks with less than 15 years of recruitment data
stock_info_filtered=subset(stock_info,n.years>=18) #242 stocks

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock.id)) #242

stock_dat2$logR_S=stock_dat2$recruits/stock_dat2$spawners

#Goal - fit each model and assess model fit using LFO-CV likelihood
s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])

#test...
stan_lfo_cv(mod=sr_mod(type='regime',par='both',loglik=T),type='regime',df=s,K=3)


#Define models (helps prevent crashing)
m1=sr_mod(type='static',ac = FALSE,loglik=T)
m2=sr_mod(type='static',ac = TRUE,loglik=T)
m3=sr_mod(type='tv',par='a',loglik=T)
m4=sr_mod(type='tv',par='b',loglik=T)
m5=sr_mod(type='tv',par='both',loglik=T)
m6=sr_mod(type='regime',par='a',loglik=T)
m7=sr_mod(type='regime',par='b',loglik=T)
m8=sr_mod(type='regime',par='both',loglik=T)

#Summary data frame and store for pointwise LLs
loglik_summary=data.frame(stock=stock_info_filtered$stock.name,LL_m1=NA,LL_m2.1=NA,LL_m2.3=NA,LL_m2.5=NA,LL_m3.1=NA,LL_m3.3=NA,
                          LL_m3.5=NA,LL_m4.1=NA,LL_m4.3=NA,LL_m4.5=NA,LL_m5.1=NA,LL_m5.3=NA,LL_m5.5=NA,LL_m6.1=NA,LL_m6.3=NA,LL_m6.5=NA,
                          LL_m6.1w=NA,LL_m6.3w=NA,LL_m6.5w=NA,LL_m7.1=NA,LL_m7.3=NA,LL_m7.5=NA,LL_m7.1w=NA,LL_m7.3w=NA,LL_m7.5w=NA,
                          LL_m8.1=NA,LL_m8.3=NA,LL_m8.5=NA,LL_m8.1w=NA,LL_m8.3w=NA,LL_m8.5w=NA)
pw_loglik=list()
for(i in 1:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  
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
  ll8<- stan_lfo_cv(mod=m8,type='regime',df=s,L=10,K=2)
  
  
  
}

plot_resid_t(resid_trends,m.col=3, l95.col = 4,u95.col=5,sp='Sockeye')
