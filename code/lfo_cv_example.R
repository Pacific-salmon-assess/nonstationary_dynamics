rm(list=ls())
library(here);library(cmdstanr)
here()

source(here('code','dlm-wrapper.R'))

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
stan_mod_refit<- function(mod,newdata,oos=NULL){
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
stan_mod_lfo_cv=function(mod,tv=1,df,L){
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
    fit_past<- stan_mod_refit(mod=mod,newdata=df_oos,oos=i+1)
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
    exact_elpds_1b <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b=exact_elpds_5b[-(1:L)]
    
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }
}
#Returns the pointwise log-likelihoods for L+1:N
dlm_mod_refit<- function(mod,newdat){
  require(dlm)
  source(here('code','dlm-wrapper.R'))
  if(mod==2){
    avary<-fitDLM(newdat, alpha_vary = TRUE, beta_vary = FALSE)
    return(avary)
  }
  if(mod==3){
    bvary<-fitDLM(newdat, alpha_vary = FALSE, beta_vary = TRUE)
    return(bvary)
  }
  if(mod==4){
    abvary<-fitDLM(newdat, alpha_vary = TRUE, beta_vary = TRUE)
    return(abvary)
  }
}

dlm_mod_lfo_cv=function(mod,df,L){
  #mod = model to fit (model name for cmdstanr)
  #tv = 0 for static model; 1 for time-varying (for calculating elpds)
  #df = full data frame
  #L = starting point for LFO-CV (min. 10)
  exact_elpds_1b<- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
  exact_elpds_3b<- numeric(nrow(df)) #loglik for average of last 3-years of productivity/capacity
  exact_elpds_5b<- numeric(nrow(df))#loglik for average of last 5-years of productivity/capacity
  for (i in L:(nrow(df) - 1)) {
    past <- 1:i
    oos <- i + 1
    df_past <- df[past, , drop = FALSE]
    df_oos <- df[c(past, oos), , drop = FALSE]
    fit_past<- dlm_mod_refit(mod=mod,newdat=df_past)
    rs_pred_1b=fit_past$results$alpha[i]+fit_past$results$beta[i]*df_oos$spwn[i + 1]
    rs_pred_3b=mean(fit_past$results$alpha[(i-2):i])+mean(fit_past$results$beta[(i-2):i])*df_oos$spwn[i + 1]
    rs_pred_5b=mean(fit_past$results$alpha[(i-4):i])+mean(fit_past$results$beta[(i-4):i])*df_oos$spwn[i + 1]
    
    exact_elpds_1b[i+1] <- log(dnorm(log(df_oos$rec[i + 1]/df_oos$spwn[i + 1]),mean=rs_pred_1b,sd=fit_past$sd.est[1]))
    exact_elpds_3b[i+1] <- log(dnorm(log(df_oos$rec[i + 1]/df_oos$spwn[i + 1]),mean=rs_pred_3b,sd=fit_past$sd.est[1]))
    exact_elpds_5b[i+1] <- log(dnorm(log(df_oos$rec[i + 1]/df_oos$spwn[i + 1]),mean=rs_pred_5b,sd=fit_past$sd.est[1]))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
}

#
lm_mod_refit<- function(y,x,ac=FALSE){
  require(nlme)
  if(ac==FALSE){
    lm_rf<- lm(y~x)
    return(lm_rf)
  }
  if(ac==TRUE){
    lm_ac_rf<- nlme::gls(y~x,correlation=corAR1(form = ~ 1))
    return(lm_ac_rf)
  }
}

lm_mod_lfo_cv=function(df,ac=F,L){
  #df = full data frame
  #ac = autocorrelation, if ac=T then implement AR-1
  #L = starting point for LFO-CV (min. 10)
  if(ac==F){
    exact_elpds_1b <- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    for (i in L:(nrow(df) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      fit_past_lm<- lm_mod_refit(y=df_past$y,x=df_past$spawners,ac=FALSE)
      rs_pred_1b=fit_past_lm$coefficients[1]+fit_past_lm$coefficients[2]*df_oos$spawners[i + 1]
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_1b,sd=sigma(fit_past_lm)))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    return(exact_elpds_1b)
    }
  if(ac==T){
    exact_elpds_ac_1b <- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_ac_3b <- numeric(nrow(df)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_ac_5b <- numeric(nrow(df))#loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(df) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      
      fit_past_ac_lm<- lm_mod_refit(y=df_past$y,x=df_past$spawners,ac=TRUE)
      #AR-1 model - extract phi and residuals 1-year back, 3-years back, and 5-years back
      phi=intervals(fit_past_ac_lm)$corStruct[2]
      ac_resids_1b=residuals(fit_past_ac_lm)[i]
      ac_resids_3b=mean(residuals(fit_past_ac_lm)[(i-2):i])
      ac_resids_5b=mean(residuals(fit_past_ac_lm)[(i-4):i])
      
      rs_pred_ac_1b=fit_past_ac_lm$coefficients[1]+fit_past_ac_lm$coefficients[2]*df_oos$spawners[i + 1]+phi*ac_resids_1b
      rs_pred_ac_3b=fit_past_ac_lm$coefficients[1]+fit_past_ac_lm$coefficients[2]*df_oos$spawners[i + 1]+phi*ac_resids_3b
      rs_pred_ac_5b=fit_past_ac_lm$coefficients[1]+fit_past_ac_lm$coefficients[2]*df_oos$spawners[i + 1]+phi*ac_resids_5b
      
      exact_elpds_ac_1b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_ac_1b,sd=fit_past_ac_lm$sigma))
      exact_elpds_ac_3b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_ac_3b,sd=fit_past_ac_lm$sigma))
      exact_elpds_ac_5b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_ac_5b,sd=fit_past_ac_lm$sigma))
    }
    exact_elpds_ac_1b=exact_elpds_ac_1b[-(1:L)]
    exact_elpds_ac_3b=exact_elpds_ac_3b[-(1:L)]
    exact_elpds_ac_5b=exact_elpds_ac_5b[-(1:L)]
    return(list(exact_elpds_ac_1b,exact_elpds_ac_3b,exact_elpds_ac_5b))
  }
}

tmb_ricker_refit<- function(df){
  srm <- lm(df$y~ df$spawners)
  SRdata<-list(obs_logRS=df$y,obs_S=df$spawners)
  
  parameters_simple<- list(
    alpha=srm$coefficients[1],
    logbeta = log(ifelse(-srm$coefficients[2]<0,1e-08,-srm$coefficients[2])),
    logsigobs=log(.4)
  )
  obj<- TMB::MakeADFun(SRdata,parameters_simple,DLL="Ricker_simple")
  opt_simple <- nlminb(obj$par,obj$fn,obj$gr)
  if(opt_simple$convergence==0){
    return(list(obj$report(),opt_simple))  
  }
}

tmb_mod_refit<- function(df,tv.par=c('alpha','beta','both')){
  srm <- lm(df$y~df$spawners)
  SRdata<-list(obs_logRS=df$y,obs_S=df$spawners)
  if(tv.par=='alpha'){
    parameters<- list(
      alphao=srm$coefficients[1],
      logbeta = log(ifelse(-srm$coefficients[2]<0,1e-08,-srm$coefficients[2])),
      #rho=.5,
      #logvarphi=0,
      logsigobs=log(.4),
      logsiga=log(.4),
      alpha=rep(srm$coefficients[1],length(df$y))
    )
    
    obja<- MakeADFun(SRdata,parameters,DLL="Ricker_tva_Smax",random="alpha")#,lower = -Inf, upper = Inf)
    obja$fn()
    obja$gr()
    opta=nlminb(obja$par,obja$fn,obja$gr)
    if(opta$convergence==0){
      return(list(obja$report(),opta))  
    }
  }
  
  if(tv.par=='beta'){
    parameters<- list(
      logbetao = log(ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[[2]])),
      alpha=srm$coefficients[[1]],
      logsigobs=log(.4),
      logsigb=log(.4),
      #rho=.4,
      #logvarphi= 0.5,
      logbeta=log(rep(-srm$coefficients[2],length(df$spawners)))
    )
    objb <- MakeADFun(SRdata,parametersb,DLL="Ricker_tvb",random="logbeta",lower=c(-10,-20,-6,-6),
                      upper=c(10,0,2,2))
    objb$fn()
    objb$gr()
    optb=nlminb(objb$par,objb$fn,objb$gr)
    if(optb$convergence==0){
      return(list(objb$report(),optb))  
    }
  }
  
  if(tv.par=='both'){
    parametersab<- list(
      logbetao = log(ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[[2]])),
      alphao=srm$coefficients[[1]],
      logsigobs=log(.5),
      logsiga=log(.1),
      logsigb=log(.1),
      logbeta=rep(log(ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[[2]])),length(df$spawners)),
      alpha=rep(srm$coefficients[[1]],length(df$spawners))
    )
    
    objab <- MakeADFun(SRdata,parametersab,DLL="Ricker_tva_tvb",random=c("logbeta","alpha"))
    objab$fn()
    objab$gr()
    optab=nlminb(objab$par,objab$fn,objab$gr)
    if(optab$convergence==0){
      return(list(objab$report(),optab))  
    }
  }
}

tmb_mod_lfo_cv=function(df,tv.par=c('static','alpha','beta','both'),L){
  #df = full data frame
  #ac = autocorrelation, if ac=T then implement AR-1
  #L = starting point for LFO-CV (min. 10)
  if(tv.par=='static'){
    exact_elpds_1b <- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    for (i in L:(nrow(df) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      
      fit_past_tmb<- tmb_ricker_refit(df=df_past)
      rs_pred_1b=fit_past_tmb[[1]]$alpha-fit_past_tmb[[1]]$beta*df_oos$spawners[i + 1]
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_1b,sd=exp(fit_past_tmb[[2]]$par[3])))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    return(exact_elpds_1b)
  }
  if(tv.par=='alpha'){
    exact_elpds_1b <- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(df)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(df))#loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(df) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      
      fit_past_tv_a_tmb<- tmb_mod_refit(df=df_past,tv.par='alpha')
    
      rs_pred_1b=fit_past_tv_a_tmb[[1]]$alpha[i]-fit_past_tv_a_tmb[[1]]$beta*df_oos$spawners[i + 1]
      rs_pred_3b=mean(fit_past_tv_a_tmb[[1]]$alpha[(i-2):i])-fit_past_tv_a_tmb[[1]]$beta*df_oos$spawners[i + 1]
      rs_pred_5b=mean(fit_past_tv_a_tmb[[1]]$alpha[(i-4):i])-fit_past_tv_a_tmb[[1]]$beta*df_oos$spawners[i + 1]
      
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_1b,sd=exp(fit_past_tv_a_tmb[[2]]$par[3])))
      exact_elpds_3b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_3b,sd=exp(fit_past_tv_a_tmb[[2]]$par[3])))
      exact_elpds_5b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_5b,sd=exp(fit_past_tv_a_tmb[[2]]$par[3])))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }
  if(tv.par=='beta'){
    exact_elpds_1b <- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(df)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(df))#loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(df) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      
      fit_past_tv_b_tmb<- tmb_mod_refit(df=df_past,tv.par='beta')
      
      rs_pred_1b=fit_past_tv_b_tmb[[1]]$alpha-fit_past_tv_b_tmb[[1]]$beta[i]*df_oos$spawners[i + 1]
      rs_pred_3b=fit_past_tv_b_tmb[[1]]$alpha-mean(fit_past_tv_b_tmb[[1]]$beta[(i-2):i])*df_oos$spawners[i + 1]
      rs_pred_5b=fit_past_tv_b_tmb[[1]]$alpha-mean(fit_past_tv_b_tmb[[1]]$beta[(i-4):i])*df_oos$spawners[i + 1]
      
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_1b,sd=exp(fit_past_tv_b_tmb[[2]]$par[3])))
      exact_elpds_3b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_3b,sd=exp(fit_past_tv_b_tmb[[2]]$par[3])))
      exact_elpds_5b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_5b,sd=exp(fit_past_tv_b_tmb[[2]]$par[3])))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }
  if(tv.par=='both'){
    exact_elpds_1b <- numeric(nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(df)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(df))#loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(df) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- df[past, , drop = FALSE]
      df_oos <- df[c(past, oos), , drop = FALSE]
      
      fit_past_tv_ab_tmb<- tmb_mod_refit(df=df_past,tv.par='both')
      
      rs_pred_1b=fit_past_tv_ab_tmb[[1]]$alpha[i]-fit_past_tv_ab_tmb[[1]]$beta[i]*df_oos$spawners[i + 1]
      rs_pred_3b=mean(fit_past_tv_ab_tmb[[1]]$alpha[(i-2):i])-mean(fit_past_tv_ab_tmb[[1]]$beta[(i-2):i])*df_oos$spawners[i + 1]
      rs_pred_5b=mean(fit_past_tv_ab_tmb[[1]]$alpha[(i-4):i])-mean(fit_past_tv_ab_tmb[[1]]$beta[(i-4):i])*df_oos$spawners[i + 1]
      
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_1b,sd=exp(fit_past_tv_ab_tmb[[2]]$par[3])))
      exact_elpds_3b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_3b,sd=exp(fit_past_tv_ab_tmb[[2]]$par[3])))
      exact_elpds_5b[i+1] <- log(dnorm(df_oos$y[i],mean=rs_pred_5b,sd=exp(fit_past_tv_ab_tmb[[2]]$par[3])))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b))
  }
}

##Data####
fulldata<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_jun2022.csv'))

#example stock
stock1<- subset(fulldata,stock.id==unique(stock.id)[6])

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
file2miss<- file.path(cmdstan_path(),'nonstationary dynamics',"ricker_linear_varying_a_miss.stan")

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
mod2m<- cmdstan_model(file2miss)

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
              TT=max(df$year)-min(df$year)+1,
              ii=df$year-min(df$year)+1,
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

m2m<- mod2m$sample(
  data = list(R_S = s$logrps,
              TT=max(s$broodyear)-min(s$broodyear)+1,
              ii=s$broodyear-min(s$broodyear)+1,
              N=nrow(s),
              S=c((s$spawners))),
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
ll1.1=stanmod_lfo_cv(mod=mod1.1oos,tv=1,df=df,L=10)
ll2=mod_lfo_cv(mod=mod2oos,df=df,L=10)
ll3=mod_lfo_cv(mod=mod3oos,df=df,L=10)
ll4=mod_lfo_cv(mod=mod4oos,df=df,L=10)

ELPDS=c(mod1=sum(ll1),mod1.1=sum(ll1.1),mod2=sum(ll2),mod3=sum(ll3),mod4=sum(ll4))

#Model summary list - for running through all stocks
mod_list<- list()

mod_list[[1]]=data.frame(stock=NA,species=NA,alpha=NA,alpha.l95=NA,alpha.u95=NA,Smax=NA,Smax.l95=NA,Smax,u95=NA,sigma=NA,sigma.l95=NA,sigma.u95=NA)
mod_list[[2]]=data.frame(stock=NA,species=NA,alpha=NA,alpha.l95=NA,alpha.u95=NA,Smax=NA,Smax.l95=NA,Smax,u95=NA,phi=NA,phi.l95=NA,phi.u95=NA,sigma=NA,sigma.l95=NA,sigma.u95=NA)
mod_list[[3]]=data.frame(stock=NA,species=NA,alpha=NA,alpha.l95=NA,alpha.u95=NA,Smax=NA,Smax.l95=NA,Smax,u95=NA,sigma=NA,sigma.l95=NA,sigma.u95=NA)
mod_list[[4]]=data.frame(stock=NA,species=NA,alpha=NA,alpha.l95=NA,alpha.u95=NA,Smax=NA,Smax.l95=NA,Smax,u95=NA,sigma=NA,sigma.l95=NA,sigma.u95=NA)
mod_list[[5]]=data.frame(stock=NA,species=NA,alpha=NA,alpha.l95=NA,alpha.u95=NA,Smax=NA,Smax.l95=NA,Smax,u95=NA,sigma=NA,sigma.l95=NA,sigma.u95=NA)


#DLM assessment####
df <- data.frame(byr=df$byr,
                     spwn=df$spwn,
                     rec=df$rec)

dlm_avary=dlm_mod_lfo_cv(mod=2,df=df,L=10)


#LM assessment####
lm_mod_lfo_cv(df=df,ac=F,L=10)
lm_mod_lfo_cv(df=df,ac=T,L=10)


#TMB assessment####
library(TMB)

compile(here("code","TMBmodels","Ricker_simple.cpp"))
dyn.load(dynlib(here("code","TMBmodels","Ricker_simple")))

compile(here("code","TMBmodels","Ricker_tva_Smax.cpp"))
dyn.load(dynlib(here("code","TMBmodels","Ricker_tva_Smax")))

#Model 4 tv  b
compile(here("code","TMBmodels","Ricker_tvb.cpp"))
dyn.load(dynlib(here("code","TMBmodels","Ricker_tvb")))

compile(here("code","TMBmodels","Ricker_tvlogb.cpp"))
dyn.load(dynlib(here("code","TMBmodels","Ricker_tvlogb")))

compile(here("code","TMBmodels","Ricker_tva_tvb.cpp"))
dyn.load(dynlib(here("code","TMBmodels","Ricker_tva_tvb")))


tmb_mod_lfo_cv(df=df,tv.par='static',L=10)
