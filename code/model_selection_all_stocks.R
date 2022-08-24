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
L=10
i=L

s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])

#test...
stan_lfo_cv(mod=sr_mod(type='regime',par='both',loglik=T),type='regime',df=s,K=3)


for(i in 1:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  

  #model 1 - static Ricker
  m1<- rstan::stan(file = here('code','stan','ricker_linear_oos.stan'), 
                   data = list(N=nrow(s),
                               R_S =s$logR_S,
                               S=s$spawners,
                               y),
                   pars = c('c','sd_site','sd_dv','sd_dmy','sd_r','sd_q','x','a_yr','beta'),
                   control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 150, chains = 4, iter = 450, thin = 1)
  
  
  #model 2 - static autocorrelated Ricker
  
  #model 3 - dynamic productivity Ricker
  
  #model 4 - dynamic capacity Ricker
  
  #model 5 - dynamic productivity & capacity Ricker
  
  #model 6 - productivity regime shift
  
  #model 7 - capacity regime shift
  
  #model 8 - productivity and capacity regime shift
  
}

plot_resid_t(resid_trends,m.col=3, l95.col = 4,u95.col=5,sp='Sockeye')
