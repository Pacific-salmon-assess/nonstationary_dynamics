rm(list=ls())
library(here);library(cmdstanr)
here()

##Functions####
##Functions for LFO-CV
log_sum_exp <- function(x) {
  max_x <- max(x)  
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

#Refit cmdstanr model
mod_refit_cmd<- function(mod,newdata,oos=NULL){
  if(is.null(oos)==T){
    r=mod$sample(
      data = list(R_S = newdata$y,
                  N=nrow(newdata),
                  TT=as.numeric(factor(newdata$year)),
                  S=c((newdata$spawners))),
      seed = 123, 
      chains = 6, 
      parallel_chains = 6,
      iter_warmup = 200,
      iter_sampling = 500,
      refresh = 200,
      adapt_delta = 0.99,
      max_treedepth = 20 # print update every 500 iters
    )
  }
  
  if(is.null(oos)==F){
    oosdata=newdata[oos,]
    newdata=newdata[-oos,]
    r=mod$sample(
      data = list(R_S = newdata$y,
                  N=nrow(newdata),
                  TT=as.numeric(factor(newdata$year)),
                  S=c((newdata$spawners)),
                  y_oos=oosdata$y,
                  x_oos=oosdata$spawners),
      seed = 123, 
      chains = 6, 
      parallel_chains = 6,
      iter_warmup = 200,
      iter_sampling = 500,
      refresh = 200,
      adapt_delta = 0.99,
      max_treedepth = 20 # print update every 500 iters
    )
  } 
  
  return(r)
}

#Leave-future-out cross-validation function
mod_lfo_cv=function(mod,tv=1,df,L){
  #mod = model to fit (model name for cmdstanr)
  #tv = 0 for static model; 1 for time-varying (for calculating elpds)
  #df = full data frame
  #L = starting point for LFO-CV (min. 10)
  loglik_exact <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for static model
  loglik_exact_1b <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for 1-year back estimates of productivity/capacity
  loglik_exact_3b <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for average of last 3-years of productivity/capacity
  loglik_exact_5b <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for average of last 5-years of productivity/capacity
  
  for (i in L:(nrow(df) - 1)) {
    past <- 1:i
    oos <- i + 1
    df_past <- df[past, , drop = FALSE]
    df_oos <- df[c(past, oos), , drop = FALSE]
    fit_past<- mod_refit_cmd(mod=mod,newdata=df_oos,oos=i+1)
    if(tv==0){
      loglik_exact[,i+1]<- fit_past$draws(variables=c('log_lik_oos'),format='draws_matrix')
    }else{
      loglik_exact_1b[, i + 1] <- fit_past$draws(variables=c('log_lik_oos_1b'),format='draws_matrix')
      loglik_exact_3b[, i + 1] <- fit_past$draws(variables=c('log_lik_oos_3b'),format='draws_matrix')
      loglik_exact_5b[, i + 1] <- fit_past$draws(variables=c('log_lik_oos_5b'),format='draws_matrix')
    }
  }
  if(tv==0){
    exact_elpds<- apply(loglik_exact, 2, log_mean_exp); exact_elpds=exact_elpds[-(1:L)]
    return(exact_elpds)
  }else{
    exact_elpds_1b_rs <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b_rs=exact_elpds_1b_rs[-(1:L)]
    exact_elpds_3b_rs <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b_rs=exact_elpds_3b_rs[-(1:L)]
    exact_elpds_5b_rs <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b_rs=exact_elpds_5b_rs[-(1:L)]
    
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }
}
#Returns the pointwise log-likelihoods for L+1:N

##Data####
fulldata<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_jun2022.csv'))

#example stock
stock1<- subset(fulldata,stock.id==unique(stock.id)[12])

stock1$logRS=log(stock1$recruits/stock1$spawners)
plot(recruits~spawners,data=stock1)
text(stock1$spawners,stock1$recruits,stock1$broodyear)

N <- length(stock1$logRS)
df <- data.frame(
  y = as.numeric(stock1$logRS),
  spawners=as.numeric(stock1$spawners),
  year = as.numeric(stock1$broodyear),
  time = 1:N
)


#CmdstanR models####
set_cmdstan_path() #copy these models to local cmdstan folder under 'nonstationary dynamics'

#identify files for each model
file1<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear.stan")
file1oos<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_oos.stan")
#static S-R model
file1.1<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_ac_resids.stan")
file1.1oos<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_ac_resids_oos.stan")
#1.1 = static S-R with autocorrelated residuals
file2<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_varying_a.stan")
file2oos<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_varying_a_oos.stan")
#2 = time-varying productivity
file3<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_varying_b.stan")
file3oos <- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_varying_b_oos.stan")
#3 = time-varying capacity
file4<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_varying_a_and_b.stan")
file4oos<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_varying_a_and_b_oos.stan")
#4 = time-varying productivity and capacity

#models - fit to full data:
mod1<- cmdstan_model(file1)
mod1.1<- cmdstan_model(file1.1)
mod2<- cmdstan_model(file2)
mod3<- cmdstan_model(file3)
mod4<- cmdstan_model(file4)

#models - with out of sample prediction:
mod1oos<- cmdstan_model(file1oos)
mod1.1oos<- cmdstan_model(file1.1oos)
mod2oos<- cmdstan_model(file2oos)
mod3oos<- cmdstan_model(file3oos)
mod4oos<- cmdstan_model(file4oos)

#fit models
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

m1.1<- mod1.1$sample(
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

m2<- mod2$sample(
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

m3<- mod3$sample(
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

m4<- mod4$sample(
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

#LFO-CV of models
###Model assessment setup
ll1=mod_lfo_cv(mod=mod1oos,tv=0,df=df,L=10)
ll1.1=mod_lfo_cv(mod=mod1.1oos,tv=1,df=df,L=10)
ll2=mod_lfo_cv(mod=mod2oos,df=df,L=10)
ll3=mod_lfo_cv(mod=mod3oos,df=df,L=10)
ll4=mod_lfo_cv(mod=mod4oos,df=df,L=10)

ELPDS=c(mod1=sum(ll1),mod1.1=sum(ll1.1),mod2=sum(ll2),mod3=sum(ll3),mod4=sum(ll4))

#Model summary list
mod_list<- list()

mod_list[[1]]=data.frame(stock=NA,species=NA,alpha=NA,alpha.l95=NA,alpha.u95=NA,Smax=NA,Smax.l95=NA,Smax,u95=NA,sigma=NA,sigma.l95=NA,sigma.u95=NA)
mod_list[[2]]=data.frame(stock=NA,species=NA,alpha=NA,alpha.l95=NA,alpha.u95=NA,Smax=NA,Smax.l95=NA,Smax,u95=NA,phi=NA,phi.l95=NA,phi.u95=NA,sigma=NA,sigma.l95=NA,sigma.u95=NA)
mod_list[[3]]=data.frame(stock=NA,species=NA,alpha=NA,alpha.l95=NA,alpha.u95=NA,Smax=NA,Smax.l95=NA,Smax,u95=NA,sigma=NA,sigma.l95=NA,sigma.u95=NA)
mod_list[[4]]=data.frame(stock=NA,species=NA,alpha=NA,alpha.l95=NA,alpha.u95=NA,Smax=NA,Smax.l95=NA,Smax,u95=NA,sigma=NA,sigma.l95=NA,sigma.u95=NA)
mod_list[[5]]=data.frame(stock=NA,species=NA,alpha=NA,alpha.l95=NA,alpha.u95=NA,Smax=NA,Smax.l95=NA,Smax,u95=NA,sigma=NA,sigma.l95=NA,sigma.u95=NA)
