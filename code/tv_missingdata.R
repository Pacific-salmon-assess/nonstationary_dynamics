#Test stan scripts for missing data


#log(b) stan
#log(b) TMB
#log(b) TMB centralized
#(b) dlm

#Check estimates and convergence in real data 
library(here)
#library(TMB)
library(cmdstanr);
library(loo);
#library(dlm)
library(ggplot2)
library(gridExtra)

tvamiss <- file.path(cmdstan_path(),'timevarmodels', "ricker_linear_varying_a_miss.stan")
mod_tvamiss <- cmdstan_model(tvamiss)

tva <- file.path(cmdstan_path(),'timevarmodels', "ricker_linear_varying_a.stan")
mod_tva <- cmdstan_model(tva)


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
            S = c(s$spawners[!is.na(s$logR_S)]))

  
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
  
params<- fit3$draws(format='df',variables=c('log_a','b','log_b','sigma_a','sigma_e'))
parssummary<-summary(params)
apply(params[,grep("log_a",names(params))]
head(parssummary)


#----------------------------------------------------------------

s<-sal_dat[sal_dat$stock.id==197,]
s$spawners[s$spawners==0] <- NA
s$logR_S <- log(s$recruits/s$spawners)
s$logR_S[is.infinite(s$logR_S)]<-NA

sc<-s[!is.na(s$logR_S),]


datasc=list(N = nrow(sc),
          TT = as.numeric(factor(sc$broodyear)),
          R_S = sc$logR_S,
          S = c(sc$spawners))

  
fit <- mod_tva$sample(
      data = datasc,
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
  
paramsc<- fit$draws(format='df',variables=c('log_a','b','log_b','sigma_a','sigma_e'))
parssummarysc<-summary(paramsc)



dt<- data.frame( value= parssummary[grep("log_a",parssummary$variable),]$median,
  low=parssummary[grep("log_a",parssummary$variable),]$q5,
  high=parssummary[grep("log_a",parssummary$variable),]$q95,
    yr=s$broodyear, 
    model="with NAs")


dtm<- data.frame( value=parssummarysc[grep("log_a",parssummarysc$variable),]$median,
  low=parssummarysc[grep("log_a",parssummarysc$variable),]$q5,
  high=parssummarysc[grep("log_a",parssummarysc$variable),]$q95,
  yr=sc$broodyear, 
  model="no NAs")

nadtm<-data.frame(value=c(NA,NA),
  low=c(NA,NA),
  high=c(NA,NA),
  yr=s$broodyear[is.na(s$logR_S)],
  model="no NAs")

dd<-rbind(dt,dtm,nadtm)

p<-ggplot(dd)+
 geom_line(aes(x=yr,y=value,color=model))+
 geom_ribbon(aes(x=yr,ymin=low,ymax=high,fill=model),alpha=.4)+
 theme_bw(16)+
 scale_color_viridis_d() +scale_fill_viridis_d()
 p

