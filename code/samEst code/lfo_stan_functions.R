#' stan_refit function
#'
#' This function refits a stan model (via rstan) to iteratively estimate the out-of-sample likelihood
#' @param sm Compiled stan model, which can be generated with sr_mod function.
#' @param newdata new data to feed to refit models.
#' @param oos What row of the new data should be treated as the out-of-sample test
#' @param regime TRUE or FALSE statement - is this a regime shift model or not
#' @param K Number of potential regime shifts - default NULL
#' @return returns the model fit
#' @export
#' @examples
#' r=stan_refit(sm=mod3,newdata=df,oos=12)
stan_refit<- function(sm,newdata,oos,regime=FALSE,K=NULL){
  #mod = model file name - eg. 'ricker_linear_oos.stan'
  #newdata = data to train model
  #oosdata = data to predict onto
  #regime = TRUE or FALSE for regime shift models (have different data inputs)
  #K = number of potential regimes (2 or 3)
  
  oosdata=newdata[oos,]
  newdata=newdata[-oos,]
  if(regime==FALSE){
    r = rstan::sampling(sm, 
                        data = list(N=nrow(newdata),
                                    L=max(newdata$by)-min(newdata$by)+1,
                                    ii=newdata$by-min(newdata$by)+1,
                                    R_S =newdata$logRS,
                                    S=newdata$S,
                                    y_oos=oosdata$logRS,
                                    x_oos=oosdata$S),
                        control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
  }
  if(regime==TRUE){
    r = rstan::sampling(sm, 
                        data = list(N=nrow(newdata),
                                    R_S=newdata$logRS,
                                    S=newdata$S,
                                    K=K,
                                    alpha_dirichlet=rep(1,K),
                                    y_oos=oosdata$logRS,
                                    x_oos=oosdata$S), #prior for state transition probabilities (this makes them equal)
                        control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
  }
  
  return(r)
}


#' stan_lfo_cv function
#'
#' This function implements exact Leave-Future-Out Cross-Validation (LFO-CV) (see https://arxiv.org/abs/1902.06281), 
#' which will iteratively estimate model likelihood in an out-of-sample context using past data as a training set.
#' @param mod Compiled stan model, which can be generated with sr_mod function.
#' @param type The type of stock-recruitment model that is being fit. Either 'static' for the time-invariant model (model 1: static Ricker), 
#' 'tv' for time-varying models (including autocorrelated residual model), or 'regime' for regime shift models
#' @param df Dataframe to use
#' @param L Minimum dataset to retain, default is to 10 years of data
#' @param K Number of potential regime shifts - default NULL
#' @return returns the pointwise out-of-sample log likelihoods. For 'tv' or 'regime' models returns a vector that uses either 1. the last year's time-varying parameters, 
#' 2. an average of the last 3 years of the time-varying parameters, or 5. an average of the last 5 years of the time-varying parameters. For 'regime' models also includes
#' these estimates that are also probability-weighted (ie. prob regime x regime parameter) in addition to the parameters for the most likely state 1, 3 or 5 years back.
#' @export
#' @examples
#' r=stan_refit(sm=mod3,newdata=df,oos=12)
stan_lfo_cv=function(mod,type=c('static','tv','regime'),df,L=10,K=NULL){
  #mod = model to fit (model name for cmdstanr)
  #tv = 0 for static model; 1 for time-varying (for calculating elpds)
  #df = full data frame
  #L = starting point for LFO-CV (default 10)
  # K = number of regimes
  
  loglik_exact <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for static model
  loglik_exact_1b <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for 1-year back estimates of productivity/capacity
  loglik_exact_3b <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for average of last 3-years of productivity/capacity
  loglik_exact_5b <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for average of last 5-years of productivity/capacity
  if(type=='regime'){
    loglik_exact_1bw <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for 1-year back estimates of productivity/capacity
    loglik_exact_3bw <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for average of last 3-years of productivity/capacity
    loglik_exact_5bw <- matrix(nrow = 3000, ncol = nrow(df)) #loglik for average of last 5-years of productivity/capacity
  }
  for (i in L:(nrow(df) - 1)){
    past <- 1:i
    oos <- i + 1
    df_past <- df[past, , drop = FALSE]
    df_oos <- df[c(past, oos), , drop = FALSE]
    if(type=='static'){
      fit_past<- stan_refit(sm=mod,newdata=df_oos,oos=i+1)
      ll=rstan::extract(fit_past,pars=c('log_lik_oos'))
      loglik_exact[,i+1]<- ll$log_lik_oos
      
    }
    if(type=='tv'){
      fit_past<- stan_refit(sm=mod,newdata=df_oos,oos=i+1)
      ll=rstan::extract(fit_past,pars=c('log_lik_oos_1b','log_lik_oos_3b','log_lik_oos_5b'))
      loglik_exact_1b[, i + 1] <-ll$log_lik_oos_1b
      loglik_exact_3b[, i + 1] <-ll$log_lik_oos_3b
      loglik_exact_5b[, i + 1] <-ll$log_lik_oos_5b
    }
    if(type=='regime'){
      fit_past<- stan_refit(sm=mod,newdata=df_oos,oos=i+1,regime=TRUE,K=K)
      ll=rstan::extract(fit_past,pars=c('log_lik_oos_1b','log_lik_oos_3b','log_lik_oos_5b','log_lik_oos_1bw','log_lik_oos_3bw','log_lik_oos_5bw'))
      loglik_exact_1b[, i + 1] <- ll$log_lik_oos_1b
      loglik_exact_3b[, i + 1] <- ll$log_lik_oos_3b
      loglik_exact_5b[, i + 1] <- ll$log_lik_oos_5b
      loglik_exact_1bw[, i + 1] <- ll$log_lik_oos_1bw
      loglik_exact_3bw[, i + 1] <- ll$log_lik_oos_3bw
      loglik_exact_5bw[, i + 1] <- ll$log_lik_oos_5bw
    }
  }
  
  if(type=='static'){
    exact_elpds<- apply(loglik_exact, 2, log_mean_exp); exact_elpds=exact_elpds[-(1:L)]
    r=exact_elpds
  }
  if(type=='tv'){
    exact_elpds_1b <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b=exact_elpds_5b[-(1:L)]
    
    r=rbind(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b)
  }
  if(type=='regime'){
    exact_elpds_1b <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b=exact_elpds_5b[-(1:L)]
    exact_elpds_1bw <- apply(loglik_exact_1bw, 2, log_mean_exp); exact_elpds_1bw=exact_elpds_1bw[-(1:L)]
    exact_elpds_3bw <- apply(loglik_exact_3bw, 2, log_mean_exp); exact_elpds_3bw=exact_elpds_3bw[-(1:L)]
    exact_elpds_5bw <- apply(loglik_exact_5bw, 2, log_mean_exp); exact_elpds_5bw=exact_elpds_5bw[-(1:L)]
    
    r=rbind(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b,exact_elpds_1bw,exact_elpds_3bw,exact_elpds_5bw)
  }
  return(r)
}

#' model_weights function
#'
#' This function implements estimates Pseudo Bayesian Model Averaging (Pseudo-BMA) weights based on Yao et al. 2018  (eq. 8, see http://www.stat.columbia.edu/~gelman/research/published/stacking_paper_discussion_rejoinder.pdf)
#' @param x a dataframe of pointwise out-of-sample loglikelihoods, where each row represents model likelihood predictions - these can be estimated with the stan_lfo_cv function
#' @return returns the relative weights for each of the included models
#' @export
#' @examples
#' model_weights(rbind(ll1,ll2))
model_weights<- function(x){
  #x = dataframe of pointwise log likelihoods
  elpd_1=apply(x,1,sum) #
  elpd_2=NA
  w=NA
  se_elpd=NA
  for(i in 1:nrow(x)){
    se_elpd[i]=sqrt(sum((x[i,]-(sum(x[i,])/ncol(x)))^2))
    elpd_2[i]=exp(elpd_1[i]-0.5*se_elpd[i])
  }
  for(i in 1:nrow(x)){
    w[i]=elpd_2[i]/sum(elpd_2) 
  }
  return(w)
}

