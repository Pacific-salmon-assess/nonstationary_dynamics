rm(list=ls())
library(here);library(dplyr);library(cmdstanr)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_aug2022.csv'))
stock_info<- read.csv(here('data','filtered datasets','all_stocks_info_aug2022.csv'))

source(here('code','functions.R'))
source(here('code','lfo_functions.R'))


#Remove stocks with less than 15 years of recruitment data
stock_info_filtered=subset(stock_info,n.years>=18) #242 stocks

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock.id)) #242

stock_dat2$logR_S=stock_dat2$recruits/stock_dat2$spawners

#Load cmdstanr models
file1<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_oos.stan")
file2<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_ac_resids_oos.stan")
file3<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_varying_a_oos.stan")
file4<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_varying_b_oos.stan")
file5<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_varying_a_and_b_oos.stan")
file6<- file.path(cmdstan_path(),'nonstationary dynamics',"regime_shift_mod_alpha_oos.stan")
file7<- file.path(cmdstan_path(),'nonstationary dynamics',"regime_shift_mod_beta_oos.stan")
file8<- file.path(cmdstan_path(),'nonstationary dynamics',"regime_shift_mod_alphabeta_oos.stan")

m1oos<- cmdstan_model(file1) #model 1 - static Ricker
m2oos<- cmdstan_model(file2) #model 2 - static autocorrelated Ricker
m3oos<- cmdstan_model(file3) #model 3 - dynamic productivity Ricker
m4oos<- cmdstan_model(file4) #model 4 - dynamic capacity Ricker
m5oos<- cmdstan_model(file5) #model 5 - dynamic productivity & capacity Ricker
m6oos<- cmdstan_model(file6) #model 6 - productivity regime shift
m7oos<- cmdstan_model(file7) #model 7 - capacity regime shift
m8oos<- cmdstan_model(file8) #model 8 - productivity and capacity regime shift

#Goal - fit each model and assess model fit using LFO-CV likelihood
for(i in 1:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  
  #model 1 - static Ricker
  m1<- mod1$sample(
    data = list(R_S = df$y,
                N=nrow(df),
                S=c((df$spawners))),
    seed = 123, 
    chains = 6, 
    parallel_chains = 6,
    iter_warmup = 200,
    iter_sampling = 500,
    refresh = 200,
    adapt_delta = 0.99,
    max_treedepth = 20 # print update every 500 iters
  )
  
  
  #model 2 - static autocorrelated Ricker
  
  #model 3 - dynamic productivity Ricker
  
  #model 4 - dynamic capacity Ricker
  
  #model 5 - dynamic productivity & capacity Ricker
  
  #model 6 - productivity regime shift
  
  #model 7 - capacity regime shift
  
  #model 8 - productivity and capacity regime shift
  
}

plot_resid_t(resid_trends,m.col=3, l95.col = 4,u95.col=5,sp='Sockeye')
