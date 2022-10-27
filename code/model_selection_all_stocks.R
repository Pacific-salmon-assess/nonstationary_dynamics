rm(list=ls())
library(here);library(dplyr);library(rstan)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_aug2022.csv'))
stock_info<- read.csv(here('data','filtered datasets','all_stocks_info_aug2022.csv'))

source(here('code','samEst code','stan_functions.R'))
source(here('code','samEst code','lfo_stan_functions.R'))
source(here('code','samEst code','util_functions.R'))

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
m6=sr_mod(type='hmm',par='a',loglik=T)
m7=sr_mod(type='hmm',par='b',loglik=T)
m8_1=sr_mod(type='hmm',par='both',loglik=T,caphigh=F)
m8_2=sr_mod(type='hmm',par='both',loglik=T,caphigh=T)


#Summary data frame and store for pointwise LLs
loglik_summary=data.frame(stock=stock_info_filtered$stock.name,LL_m1=NA,LL_m2.1=NA,LL_m2.3=NA,LL_m2.5=NA,LL_m3.1=NA,LL_m3.3=NA,
                          LL_m3.5=NA,LL_m4.1=NA,LL_m4.3=NA,LL_m4.5=NA,LL_m5.1=NA,LL_m5.3=NA,LL_m5.5=NA,LL_m6.1=NA,LL_m6.3=NA,LL_m6.5=NA,
                          LL_m6.1w=NA,LL_m6.3w=NA,LL_m6.5w=NA,LL_m7.1=NA,LL_m7.3=NA,LL_m7.5=NA,LL_m7.1w=NA,LL_m7.3w=NA,LL_m7.5w=NA,
                          LL_m8.1.1=NA,LL_m8.1.3=NA,LL_m8.1.5=NA,LL_m8.1.1w=NA,LL_m8.1.3w=NA,LL_m8.1.5w=NA,LL_m8.2.1=NA,LL_m8.2.3=NA,LL_m8.2.5=NA,LL_m8.2.1w=NA,LL_m8.2.3w=NA,LL_m8.2.5w=NA)
full_loglik=list() #full pointwise log likelihoods for each year by model type
mod_loglik=list() #top variant for each model type
for(i in 30:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  s<- s[complete.cases(s$spawners),]
  
  #Assess model fits for each model type
  #model 1 - static Ricker
  ll1<- stan_lfo_cv(mod=m1,type='static',df=s,L=round((2/3)*stock_info_filtered$n.years[i]))
  #model 2 - static autocorrelated Ricker
  ll2<- stan_lfo_cv(mod=m2,type='tv',df=s,L=round((2/3)*stock_info_filtered$n.years[i]))
  #model 3 - dynamic productivity Ricker
  ll3<- stan_lfo_cv(mod=m3,type='tv',df=s,L=round((2/3)*stock_info_filtered$n.years[i]))
  
  #model 4 - dynamic capacity Ricker
  ll4<- stan_lfo_cv(mod=m4,type='tv',df=s,L=round((2/3)*stock_info_filtered$n.years[i]))
  
  #model 5 - dynamic productivity & capacity Ricker
  ll5<- stan_lfo_cv(mod=m5,type='tv',df=s,L=round((2/3)*stock_info_filtered$n.years[i]))
  
  #model 6 - productivity regime shift - 2 regimes
  ll6<- stan_lfo_cv(mod=m6,type='regime',df=s,L=round((2/3)*stock_info_filtered$n.years[i]),K=2)
  
  #model 7 - capacity regime shift
  ll7<- stan_lfo_cv(mod=m7,type='regime',df=s,L=round((2/3)*stock_info_filtered$n.years[i]),K=2)
  
  #model 8 - productivity and capacity regime shift
  ll8_1<- stan_lfo_cv(mod=m8_1,type='regime',df=s,L=round((2/3)*stock_info_filtered$n.years[i]),K=2)
  #model 8 - productivity and capacity regime shift
  ll8_2<- stan_lfo_cv(mod=m8_2,type='regime',df=s,L=round((2/3)*stock_info_filtered$n.years[i]),K=2)
  
  #Take the sum of the estimated pointwise likelihood estimates (=elpd_loo)
  loglik_summary[i,2]=sum(ll1)
  loglik_summary[i,3:5]=apply(ll2,1,sum)
  loglik_summary[i,6:8]=apply(ll3,1,sum)
  loglik_summary[i,9:11]=apply(ll4,1,sum)
  loglik_summary[i,12:14]=apply(ll5,1,sum)
  loglik_summary[i,15:20]=apply(ll6,1,sum)
  loglik_summary[i,21:26]=apply(ll7,1,sum)
  loglik_summary[i,27:32]=apply(ll8_1,1,sum)
  loglik_summary[i,33:38]=apply(ll8_2,1,sum)

  
  full_loglik[[i]]=rbind(ll1,ll2,ll3,ll4,ll5,ll6,ll7,ll8_1,ll8_2)
  rownames(full_loglik[[i]])=c('m1',paste('m2',c(1,3,5),sep="_"),paste('m3',c(1,3,5),sep="_"),paste('m4',c(1,3,5),sep="_"),paste('m5',c(1,3,5),sep="_"),paste('m6',c(1,3,5,'1w','3w','5w'),sep="_"),paste('m7',c(1,3,5,'1w','3w','5w'),sep="_"),paste('m8_1',c(1,3,5,'1w','3w','5w'),sep="_"),paste('m8_2',c(1,3,5,'1w','3w','5w'),sep="_"))
  wm2=which.max(apply(ll2,1,sum)) #select best likelihood from different timeframes (1-y back, 3-y back, 5-y back)
  wm3=which.max(apply(ll3,1,sum))#best fit for model 3
  wm4=which.max(apply(ll4,1,sum))#best fit for model 4
  wm5=which.max(apply(ll5,1,sum)) #best fit for model 5
  wm6=which.max(apply(ll6,1,sum)) #best fit for model 6
  wm7=which.max(apply(ll7,1,sum)) #best fit for model 7
  wm8=which.max(apply(rbind(ll8_1,ll8_2),1,sum)) #best fit for model 8 (two combinations of higher alpha, lower cap; higher alpha, higher cap)
  if(wm8<4){
    mod_loglik[[i]]=rbind(ll1,ll2[wm2,],ll3[wm3,],ll4[wm4,],ll5[wm5,],ll6[wm6,],ll7[wm7,],ll8_1[wm8,])
    rownames(mod_loglik[[i]])[1:8]=paste('m',seq(1:8),sep='');rownames(mod_loglik[[i]])[8]='m8_1'
  }
  if(wm8>4){
    mod_loglik[[i]]=rbind(ll1,ll2[wm2,],ll3[wm3,],ll4[wm4,],ll5[wm5,],ll6[wm6,],ll7[wm7,],ll8_2[wm8-6,])
    rownames(mod_loglik[[i]])[1:8]=paste('m',seq(1:8),sep='');rownames(mod_loglik[[i]])[8]='m8_2'
  }
  write.csv(full_loglik[[i]],here('outputs','full LFO loglik',paste(i,stock_info_filtered$stock.name[i],'full_loglik.csv',sep='')))
  write.csv(mod_loglik[[i]],here('outputs','model LFO loglik',paste(i,stock_info_filtered$stock.name[i],'mod_loglik.csv',sep='')))
  
}

#Pseudo-BMA+
#only evaluate likelihood for one version of model 8:
if(wm8<=6){ll8=ll8_1}else{ll8=ll8_2;wm8=wm8-6}

##Explorations....
#Revised model weights
model_weights(mod_loglik[[1]])


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
summary(mod_weights)

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


loglik_summary_L10=data.frame(stock=stock_info_filtered$stock.name,LL_m1=NA,LL_m2.1=NA,LL_m2.3=NA,LL_m2.5=NA,LL_m3.1=NA,LL_m3.3=NA,
                          LL_m3.5=NA,LL_m4.1=NA,LL_m4.3=NA,LL_m4.5=NA,LL_m5.1=NA,LL_m5.3=NA,LL_m5.5=NA,LL_m6.1=NA,LL_m6.3=NA,LL_m6.5=NA,
                          LL_m6.1w=NA,LL_m6.3w=NA,LL_m6.5w=NA,LL_m7.1=NA,LL_m7.3=NA,LL_m7.5=NA,LL_m7.1w=NA,LL_m7.3w=NA,LL_m7.5w=NA,
                          LL_m8.1.1=NA,LL_m8.1.3=NA,LL_m8.1.5=NA,LL_m8.1.1w=NA,LL_m8.1.3w=NA,LL_m8.1.5w=NA,LL_m8.2.1=NA,LL_m8.2.3=NA,LL_m8.2.5=NA,LL_m8.2.1w=NA,LL_m8.2.3w=NA,LL_m8.2.5w=NA)
loglik_summary_L15=data.frame(stock=stock_info_filtered$stock.name,LL_m1=NA,LL_m2.1=NA,LL_m2.3=NA,LL_m2.5=NA,LL_m3.1=NA,LL_m3.3=NA,
                              LL_m3.5=NA,LL_m4.1=NA,LL_m4.3=NA,LL_m4.5=NA,LL_m5.1=NA,LL_m5.3=NA,LL_m5.5=NA,LL_m6.1=NA,LL_m6.3=NA,LL_m6.5=NA,
                              LL_m6.1w=NA,LL_m6.3w=NA,LL_m6.5w=NA,LL_m7.1=NA,LL_m7.3=NA,LL_m7.5=NA,LL_m7.1w=NA,LL_m7.3w=NA,LL_m7.5w=NA,
                              LL_m8.1.1=NA,LL_m8.1.3=NA,LL_m8.1.5=NA,LL_m8.1.1w=NA,LL_m8.1.3w=NA,LL_m8.1.5w=NA,LL_m8.2.1=NA,LL_m8.2.3=NA,LL_m8.2.5=NA,LL_m8.2.1w=NA,LL_m8.2.3w=NA,LL_m8.2.5w=NA)
loglik_summary_L20=data.frame(stock=stock_info_filtered$stock.name,LL_m1=NA,LL_m2.1=NA,LL_m2.3=NA,LL_m2.5=NA,LL_m3.1=NA,LL_m3.3=NA,
                              LL_m3.5=NA,LL_m4.1=NA,LL_m4.3=NA,LL_m4.5=NA,LL_m5.1=NA,LL_m5.3=NA,LL_m5.5=NA,LL_m6.1=NA,LL_m6.3=NA,LL_m6.5=NA,
                              LL_m6.1w=NA,LL_m6.3w=NA,LL_m6.5w=NA,LL_m7.1=NA,LL_m7.3=NA,LL_m7.5=NA,LL_m7.1w=NA,LL_m7.3w=NA,LL_m7.5w=NA,
                              LL_m8.1.1=NA,LL_m8.1.3=NA,LL_m8.1.5=NA,LL_m8.1.1w=NA,LL_m8.1.3w=NA,LL_m8.1.5w=NA,LL_m8.2.1=NA,LL_m8.2.3=NA,LL_m8.2.5=NA,LL_m8.2.1w=NA,LL_m8.2.3w=NA,LL_m8.2.5w=NA)

L10_loglik=list() #pointwise loglikelihood
L15_loglik=list() #pointwise loglikelihood
L20_loglik=list() #pointwise loglikelihood
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
  L10_logliks[[i]]=rbind(ll1_10,ll2_10,ll3_10,ll4_10,ll5_10,ll6_10,ll7_10,ll8_1_10,ll8_2_10) #pointwise loglikelihood
  L15_loglik[[i]]=rbind(ll1_15,ll2_15,ll3_15,ll4_15,ll5_15,ll6_15,ll7_15,ll8_1_15,ll8_2_15) #pointwise loglikelihood
  L20_loglik[[i]]=rbind(ll1_20,ll2_20,ll3_20,ll4_20,ll5_20,ll6_20,ll7_20,ll8_1_20,ll8_2_20) #pointwise loglikelihood
}

wm_L10=list()
which.max(apply(L10_loglik[[1]][2:4,],1,sum))
model_weights(L10_loglik[[1]][2:4,])
which.max(apply(L10_loglik[[1]][5:7,],1,sum))
model_weights(L10_loglik[[1]][5:7,])
which.max(apply(L10_loglik[[1]][8:10,],1,sum))

which.max(apply(L10_loglik[[1]][11:13,],1,sum))
which.max(apply(L10_loglik[[1]][14:19,],1,sum))
which.max(apply(L10_loglik[[1]][20:25,],1,sum))
which.max(apply(L10_loglik[[1]][26:31,],1,sum))
which.max(apply(L10_loglik[[1]][32:37,],1,sum))

##comp of model weights
if(sum(L10_loglik[[1]][26,])>sum(L10_loglik[[1]][32,])){
  L10_comp1=rbind(L10_loglik[[1]][1,],L10_loglik[[1]][2,],L10_loglik[[1]][5,],L10_loglik[[1]][8,],L10_loglik[[1]][11,],L10_loglik[[1]][14,],L10_loglik[[1]][20,],L10_loglik[[1]][26,])
}else{
  L10_comp1=rbind(L10_loglik[[1]][1,],L10_loglik[[1]][2,],L10_loglik[[1]][5,],L10_loglik[[1]][8,],L10_loglik[[1]][11,],L10_loglik[[1]][14,],L10_loglik[[1]][20,],L10_loglik[[1]][32,])    
  }
model_weights(L10_comp1)

if(sum(L15_loglik[[1]][26,])>sum(L15_loglik[[1]][32,])){
  L15_comp1=rbind(L15_loglik[[1]][1,],L15_loglik[[1]][2,],L15_loglik[[1]][5,],L15_loglik[[1]][8,],L15_loglik[[1]][11,],L15_loglik[[1]][14,],L15_loglik[[1]][20,],L15_loglik[[1]][26,])
}else{
  L15_comp1=rbind(L15_loglik[[1]][1,],L15_loglik[[1]][2,],L15_loglik[[1]][5,],L15_loglik[[1]][8,],L15_loglik[[1]][11,],L15_loglik[[1]][14,],L15_loglik[[1]][20,],L15_loglik[[1]][32,])    
}
model_weights(L15_comp1)

if(sum(L20_loglik[[2]][26,])>sum(L20_loglik[[2]][32,])){
  L20_comp1=rbind(L20_loglik[[2]][1,],L20_loglik[[2]][2,],L20_loglik[[2]][5,],L20_loglik[[2]][8,],L20_loglik[[2]][11,],L20_loglik[[2]][14,],L20_loglik[[2]][20,],L20_loglik[[2]][26,])
}else{
  L20_comp1=rbind(L20_loglik[[2]][1,],L20_loglik[[2]][2,],L20_loglik[[2]][5,],L20_loglik[[2]][8,],L20_loglik[[2]][11,],L20_loglik[[2]][14,],L20_loglik[[2]][20,],L20_loglik[[2]][32,])    
}
model_weights(L20_comp1)

if(sum(L10_loglik[[2]][26,])>sum(L10_loglik[[2]][32,])){
  L10_comp1=rbind(L10_loglik[[2]][1,],L10_loglik[[2]][2,],L10_loglik[[2]][5,],L10_loglik[[2]][8,],L10_loglik[[2]][11,],L10_loglik[[2]][14,],L10_loglik[[2]][20,],L10_loglik[[2]][26,])
}else{
  L10_comp1=rbind(L10_loglik[[2]][1,],L10_loglik[[2]][2,],L10_loglik[[2]][5,],L10_loglik[[2]][8,],L10_loglik[[2]][11,],L10_loglik[[2]][14,],L10_loglik[[2]][20,],L10_loglik[[2]][32,])    
}
model_weights(L10_comp1)

if(sum(L15_loglik[[2]][26,])>sum(L15_loglik[[2]][32,])){
  L15_comp1=rbind(L15_loglik[[2]][1,],L15_loglik[[2]][2,],L15_loglik[[2]][5,],L15_loglik[[2]][8,],L15_loglik[[2]][11,],L15_loglik[[2]][14,],L15_loglik[[2]][20,],L15_loglik[[2]][26,])
}else{
  L15_comp1=rbind(L15_loglik[[2]][1,],L15_loglik[[2]][2,],L15_loglik[[2]][5,],L15_loglik[[2]][8,],L15_loglik[[2]][11,],L15_loglik[[2]][14,],L15_loglik[[2]][20,],L15_loglik[[2]][32,])    
}
model_weights(L15_comp1)

if(sum(L20_loglik[[2]][26,])>sum(L20_loglik[[2]][32,])){
  L20_comp1=rbind(L20_loglik[[2]][1,],L20_loglik[[2]][2,],L20_loglik[[2]][5,],L20_loglik[[2]][8,],L20_loglik[[2]][11,],L20_loglik[[2]][14,],L20_loglik[[2]][20,],L20_loglik[[2]][26,])
}else{
  L20_comp1=rbind(L20_loglik[[2]][1,],L20_loglik[[2]][2,],L20_loglik[[2]][5,],L20_loglik[[2]][8,],L20_loglik[[2]][11,],L20_loglik[[2]][14,],L20_loglik[[2]][20,],L20_loglik[[2]][32,])    
}
model_weights(L20_comp1)

if(sum(L10_loglik[[3]][26,])>sum(L10_loglik[[3]][32,])){
  L10_comp1=rbind(L10_loglik[[3]][1,],L10_loglik[[3]][2,],L10_loglik[[3]][5,],L10_loglik[[3]][8,],L10_loglik[[3]][11,],L10_loglik[[3]][14,],L10_loglik[[3]][20,],L10_loglik[[3]][26,])
}else{
  L10_comp1=rbind(L10_loglik[[3]][1,],L10_loglik[[3]][2,],L10_loglik[[3]][5,],L10_loglik[[3]][8,],L10_loglik[[3]][11,],L10_loglik[[3]][14,],L10_loglik[[3]][20,],L10_loglik[[3]][32,])    
}
model_weights(L10_comp1)

if(sum(L15_loglik[[3]][26,])>sum(L15_loglik[[3]][32,])){
  L15_comp1=rbind(L15_loglik[[3]][1,],L15_loglik[[3]][2,],L15_loglik[[3]][5,],L15_loglik[[3]][8,],L15_loglik[[3]][11,],L15_loglik[[3]][14,],L15_loglik[[3]][20,],L15_loglik[[3]][26,])
}else{
  L15_comp1=rbind(L15_loglik[[3]][1,],L15_loglik[[3]][2,],L15_loglik[[3]][5,],L15_loglik[[3]][8,],L15_loglik[[3]][11,],L15_loglik[[3]][14,],L15_loglik[[3]][20,],L15_loglik[[3]][32,])    
}
model_weights(L15_comp1)

if(sum(L20_loglik[[3]][26,])>sum(L20_loglik[[3]][32,])){
  L20_comp1=rbind(L20_loglik[[3]][1,],L20_loglik[[3]][2,],L20_loglik[[3]][5,],L20_loglik[[3]][8,],L20_loglik[[3]][11,],L20_loglik[[3]][14,],L20_loglik[[3]][20,],L20_loglik[[3]][26,])
}else{
  L20_comp1=rbind(L20_loglik[[3]][1,],L20_loglik[[3]][2,],L20_loglik[[3]][5,],L20_loglik[[3]][8,],L20_loglik[[3]][11,],L20_loglik[[3]][14,],L20_loglik[[3]][20,],L20_loglik[[3]][32,])    
}
model_weights(L20_comp1)

if(sum(L10_loglik[[4]][26,])>sum(L10_loglik[[4]][42,])){
  L10_comp1=rbind(L10_loglik[[4]][1,],L10_loglik[[4]][2,],L10_loglik[[4]][5,],L10_loglik[[4]][8,],L10_loglik[[4]][11,],L10_loglik[[4]][14,],L10_loglik[[4]][20,],L10_loglik[[4]][26,])
}else{
  L10_comp1=rbind(L10_loglik[[4]][1,],L10_loglik[[4]][2,],L10_loglik[[4]][5,],L10_loglik[[4]][8,],L10_loglik[[4]][11,],L10_loglik[[4]][14,],L10_loglik[[4]][20,],L10_loglik[[4]][42,])    
}
model_weights(L10_comp1)

if(sum(L15_loglik[[4]][26,])>sum(L15_loglik[[4]][42,])){
  L15_comp1=rbind(L15_loglik[[4]][1,],L15_loglik[[4]][2,],L15_loglik[[4]][5,],L15_loglik[[4]][8,],L15_loglik[[4]][11,],L15_loglik[[4]][14,],L15_loglik[[4]][20,],L15_loglik[[4]][26,])
}else{
  L15_comp1=rbind(L15_loglik[[4]][1,],L15_loglik[[4]][2,],L15_loglik[[4]][5,],L15_loglik[[4]][8,],L15_loglik[[4]][11,],L15_loglik[[4]][14,],L15_loglik[[4]][20,],L15_loglik[[4]][42,])    
}
model_weights(L15_comp1)

if(sum(L20_loglik[[4]][26,])>sum(L20_loglik[[4]][42,])){
  L20_comp1=rbind(L20_loglik[[4]][1,],L20_loglik[[4]][2,],L20_loglik[[4]][5,],L20_loglik[[4]][8,],L20_loglik[[4]][11,],L20_loglik[[4]][14,],L20_loglik[[4]][20,],L20_loglik[[4]][26,])
}else{
  L20_comp1=rbind(L20_loglik[[4]][1,],L20_loglik[[4]][2,],L20_loglik[[4]][5,],L20_loglik[[4]][8,],L20_loglik[[4]][11,],L20_loglik[[4]][14,],L20_loglik[[4]][20,],L20_loglik[[4]][42,])    
}
model_weights(L20_comp1)

which.max(apply(L15_loglik[[1]][2:4,],1,sum))
which.max(apply(L15_loglik[[1]][5:7,],1,sum))
which.max(apply(L15_loglik[[1]][8:10,],1,sum))
model_weights(L15_loglik[[1]][8:10,])
which.max(apply(L15_loglik[[1]][11:13,],1,sum))
which.max(apply(L15_loglik[[1]][14:19,],1,sum))
which.max(apply(L15_loglik[[1]][20:25,],1,sum))
which.max(apply(L15_loglik[[1]][26:31,],1,sum))
which.max(apply(L15_loglik[[1]][32:37,],1,sum))

which.max(apply(L20_loglik[[1]][2:4,],1,sum))
which.max(apply(L20_loglik[[1]][5:7,],1,sum))
which.max(apply(L20_loglik[[1]][8:10,],1,sum))
which.max(apply(L20_loglik[[1]][11:13,],1,sum))
which.max(apply(L20_loglik[[1]][14:19,],1,sum))
which.max(apply(L20_loglik[[1]][20:25,],1,sum))
which.max(apply(L20_loglik[[1]][26:31,],1,sum))
which.max(apply(L20_loglik[[1]][32:37,],1,sum))


which.max(apply(L10_loglik[[2]][2:4,],1,sum))
which.max(apply(L10_loglik[[2]][5:7,],1,sum))
which.max(apply(L10_loglik[[2]][8:10,],1,sum))
which.max(apply(L10_loglik[[2]][11:13,],1,sum))
which.max(apply(L10_loglik[[2]][14:19,],1,sum))
which.max(apply(L10_loglik[[2]][20:25,],1,sum))
which.max(apply(L10_loglik[[2]][26:31,],1,sum))
which.max(apply(L10_loglik[[2]][32:37,],1,sum))

which.max(apply(L15_loglik[[2]][2:4,],1,sum))
which.max(apply(L15_loglik[[2]][5:7,],1,sum))
which.max(apply(L15_loglik[[2]][8:10,],1,sum))
which.max(apply(L15_loglik[[2]][11:13,],1,sum))
which.max(apply(L15_loglik[[2]][14:19,],1,sum))
which.max(apply(L15_loglik[[2]][20:25,],1,sum))
which.max(apply(L15_loglik[[2]][26:31,],1,sum))
which.max(apply(L15_loglik[[2]][32:37,],1,sum))

which.max(apply(L20_loglik[[2]][2:4,],1,sum))
which.max(apply(L20_loglik[[2]][5:7,],1,sum))
which.max(apply(L20_loglik[[2]][8:10,],1,sum))
which.max(apply(L20_loglik[[2]][11:13,],1,sum))
which.max(apply(L20_loglik[[2]][14:19,],1,sum))
which.max(apply(L20_loglik[[2]][20:25,],1,sum))
which.max(apply(L20_loglik[[2]][26:31,],1,sum))
which.max(apply(L20_loglik[[2]][32:37,],1,sum))

which.max(apply(L10_loglik[[3]][2:4,],1,sum))
which.max(apply(L10_loglik[[3]][5:7,],1,sum))
which.max(apply(L10_loglik[[3]][8:10,],1,sum))
which.max(apply(L10_loglik[[3]][11:13,],1,sum))
which.max(apply(L10_loglik[[3]][14:19,],1,sum))
which.max(apply(L10_loglik[[3]][20:25,],1,sum))
which.max(apply(L10_loglik[[3]][26:31,],1,sum))
which.max(apply(L10_loglik[[3]][32:37,],1,sum))

which.max(apply(L10_loglik[[4]][2:4,],1,sum))
which.max(apply(L10_loglik[[4]][5:7,],1,sum))
which.max(apply(L10_loglik[[4]][8:10,],1,sum))
which.max(apply(L10_loglik[[4]][11:13,],1,sum))
which.max(apply(L10_loglik[[4]][14:19,],1,sum))
which.max(apply(L10_loglik[[4]][20:25,],1,sum))
which.max(apply(L10_loglik[[4]][26:31,],1,sum))
which.max(apply(L10_loglik[[4]][32:37,],1,sum))


ggplot(sp_sum, aes(fill=species, y=n, x=region))+
  scale_fill_manual(values=c("darkgray", "darkgreen", "darkblue","darksalmon","darkred"))+ 
  geom_bar(position="stack", stat="identity") + theme_minimal() +
  xlab('')+ylab('No. time-series')

ggplot(stock_info_filtered, aes(x=n.years, color=species)) +
  geom_histogram()

