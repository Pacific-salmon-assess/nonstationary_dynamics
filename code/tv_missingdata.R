#Test stan scripts for missing data


#log(b) stan
#log(b) TMB
#log(b) TMB centralized
#(b) dlm

#Check estimates and convergence in real data 


source(here('code','dlm-wrapper.R'))
library(here)
#library(TMB)
library(cmdstanr);
library(loo);
#library(dlm)
library(ggplot2)
library(gridExtra)

tvamiss <- file.path(cmdstan_path(),'timevarmodels', "ricker_linear_varying_a_miss.stan")
mod_tvamiss <- cmdstan_model(tvamiss)


sal_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_jun2022.csv'))

sal_dat[sal_dat$spawners==0,]

s<-sal_dat[sal_dat$stock.id==197,]
s$spawners[s$spawners==0] <- NA
s$logR_S <- log(s$recruits/s$spawners)
s$logR_S[is.infinite(s$logR_S)]<-NA



data=list(R_S = s$logR_S[!is.na(s$logR_S)],
            N = nrow(s),
            Nobs = sum(!is.na(s$logR_S)),
            ii= which(!is.na(s$logR_S)),
            TT = as.numeric(factor(s$broodyear)),
            S = c(s$spawnersavg))

  
    fit3 <- mod_tvamiss$sample(
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
  
    params<- fit3$draws(format='df',variables=c('log_a','b','log_b','sigma_b','sigma_e'))
  

    parssummary<-summary(params)

