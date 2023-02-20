#TMB lfo function
# Dan greenberg and modified by Catarina for the package




#' Simple Ricker model estimated with TMB
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param model string. Which parameters are time-varying? Options are c('static','rw_a','rw_b',
#' 'rw_both', 'HMM', 'HMM_a', 'HMM_b')
#' @param L starting point for LFO-CV (minimum value is 10)
#' @param siglfo string. Incating whether full variance should be used for lfo of models with random walks in parameters
#' "obs" incates that only observation variance is considered for lfo calculations, "total" indicates that sum of 
#' process and observation variances are used. Option valud only for 'alpha' tv par. 
#' 
#' 
#' @returns vector of lfo by year
#' 
#' @importFrom stats nlminb 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' tmb_mod_lfo_cv(data=harck, model=c('static'))
#' 
tmb_mod_lfo_cv=function(data, model=c('static','staticAC','rw_a','rw_b','rw_both', 'HMM', 'HMM_a','HMM_b'), 
  L=10, siglfo=c("obs","total")){
  #df = full data frame
  #ac = autocorrelation, if ac=T then implement AR-1
  #L = starting point for LFO-CV (min. 10)
  

  if(L<10){
  	warning("L values < 10 are not recommended, results may be unstable")
  }
  if(sum(c("S","logRS") %in% names(data))!=2){
  	stop("data must contain 'S' and 'logRS', please revise names(data)") 
  }

  conv_problem<-numeric(length(L:(nrow(data) - 1)))

  if(model=='static'){
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_tmb<- ricker_TMB(data=df_past)
      conv_problem[i-(L-1)] <- fit_past_tmb$conv_problem
      rs_pred_1b=fit_past_tmb$alpha-fit_past_tmb$beta*df_oos$S[i + 1]
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=fit_past_tmb$sig))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    return(list(lastparam=exact_elpds_1b,
       conv_problem=conv_problem)
      )
  }else if(model=='staticAC'){
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_tmb<- ricker_TMB(data=df_past,AC=TRUE)
      conv_problem[i-(L-1)] <- fit_past_tmb$conv_problem
      rs_pred_1b<-fit_past_tmb$alpha-fit_past_tmb$beta*df_oos$S[i + 1] + fit_past_tmb$residuals[i] * fit_past_tmb$rho
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=fit_past_tmb$sigar))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    return(list(lastparam=exact_elpds_1b,
      conv_problem=conv_problem))
  }else if(model=='rw_a'){
  
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(data)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(data)) #loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_tv_a_tmb<- ricker_rw_TMB(data=df_past,tv.par='a')
      conv_problem[i-(L-1)] <- fit_past_tv_a_tmb$conv_problem
    
      rs_pred_1b=fit_past_tv_a_tmb$alpha[i]-fit_past_tv_a_tmb$beta*df_oos$S[i + 1]
      rs_pred_3b=mean(fit_past_tv_a_tmb$alpha[(i-2):i])-fit_past_tv_a_tmb$beta*df_oos$S[i + 1]
      rs_pred_5b=mean(fit_past_tv_a_tmb$alpha[(i-4):i])-fit_past_tv_a_tmb$beta*df_oos$S[i + 1]
      
      if(siglfo=="obs"){
        exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=fit_past_tv_a_tmb$sig))
        exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3b,sd=fit_past_tv_a_tmb$sig))
        exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5b,sd=fit_past_tv_a_tmb$sig))
      }else if(siglfo=="total"){

        sigtot <-sqrt(fit_past_tv_a_tmb$sig^2+fit_past_tv_a_tmb$siga^2)
        exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=sigtot))
        exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3b,sd=sigtot))
        exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5b,sd=sigtot))
      }else{
        stop("siglfo incorrectly defined options are `total` or `obs`")
      }
      
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(lastparam=exact_elpds_1b,
    	last3paramavg=exact_elpds_3b,
    	last5paramavg=exact_elpds_5b,
      conv_problem=conv_problem))
  }else if(model=='rw_b'){
  	#stop("not defined")
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(data)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(data))#loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_tv_b_tmb<- ricker_rw_TMB(data=df_past,tv.par='b')
      conv_problem[i-(L-1)] <- fit_past_tv_b_tmb$conv_problem

      rs_pred_1b=fit_past_tv_b_tmb$alpha-fit_past_tv_b_tmb$beta[i]*df_oos$S[i + 1]
      rs_pred_3b=fit_past_tv_b_tmb$alpha-mean(fit_past_tv_b_tmb$beta[(i-2):i])*df_oos$S[i + 1]
      rs_pred_5b=fit_past_tv_b_tmb$alpha-mean(fit_past_tv_b_tmb$beta[(i-4):i])*df_oos$S[i + 1]

      if(siglfo=="obs"){
        exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=fit_past_tv_b_tmb$sig))
        exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3b,sd=fit_past_tv_b_tmb$sig))
        exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5b,sd=fit_past_tv_b_tmb$sig))
      }else if(siglfo=="total"){

        sigtot <-sqrt(fit_past_tv_b_tmb$sig^2+fit_past_tv_b_tmb$sigb^2)
        exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=sigtot))
        exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3b,sd=sigtot))
        exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5b,sd=sigtot))
      }else{
        stop("siglfo incorrectly defined options are `total` or `obs`")
      }
      

    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(lastparam=exact_elpds_1b,
    	last3paramavg=exact_elpds_3b,
    	last5paramavg=exact_elpds_5b,
      conv_problem=conv_problem))
  }else if(model=='rw_both'){
  	
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(data)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(data)) #loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_tv_ab_tmb<- ricker_rw_TMB(data=df_past,tv.par='both')
      conv_problem[i-(L-1)] <- fit_past_tv_ab_tmb$conv_problem
      
      rs_pred_1b=fit_past_tv_ab_tmb$alpha[i]-fit_past_tv_ab_tmb$beta[i]*df_oos$S[i + 1]
      rs_pred_3b=mean(fit_past_tv_ab_tmb$alpha[(i-2):i])-mean(fit_past_tv_ab_tmb$beta[(i-2):i])*df_oos$S[i + 1]
      rs_pred_5b=mean(fit_past_tv_ab_tmb$alpha[(i-4):i])-mean(fit_past_tv_ab_tmb$beta[(i-4):i])*df_oos$S[i + 1]
      
      exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i],mean=rs_pred_1b,sd=exp(fit_past_tv_ab_tmb$sig)))
      exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i],mean=rs_pred_3b,sd=exp(fit_past_tv_ab_tmb$sig)))
      exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i],mean=rs_pred_5b,sd=exp(fit_past_tv_ab_tmb$sig)))
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(lastparam=exact_elpds_1b,
    	last3paramavg=exact_elpds_3b,
    	last5paramavg=exact_elpds_5b,
      conv_problem=conv_problem))

  }else if(model=='HMM'){
    #stop("not defined")
    exact_elpds_1k <- numeric(nrow(data)) #loglik choosing a specific regime in a given year
    exact_elpds_3k  <- numeric(nrow(data))
    exact_elpds_5k  <- numeric(nrow(data))
   
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_hmm_tmb<- ricker_hmm_TMB(data=df_past,tv.par='both')
      conv_problem[i-(L-1)] <- fit_past_hmm_tmb$conv_problem
     
      alpha <- fit_past_hmm_tmb$alpha[fit_past_hmm_tmb$regime]
      beta <- fit_past_hmm_tmb$beta[fit_past_hmm_tmb$regime]
      sigma <- fit_past_hmm_tmb$sigma
      #sigmaw <- sqrt((fit_past_hmm_tmb$sigma^2)%*%fit_past_hmm_tmb$probregime)

      rs_pred_1k=alpha[i]-beta[i]*df_oos$S[i + 1]
      rs_pred_3k=mean(alpha[(i-2):i])-mean(beta[(i-2):i])*df_oos$S[i + 1]
      rs_pred_5k=mean(alpha[(i-4):i])-mean(beta[(i-4):i])*df_oos$S[i + 1]
      
      
      exact_elpds_1k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1k,sd=sigma))
      exact_elpds_3k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3k,sd=sigma))
      exact_elpds_5k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5k,sd=sigma))
    }
    exact_elpds_1k=exact_elpds_1k[-(1:L)]
    exact_elpds_3k=exact_elpds_3k[-(1:L)]
    exact_elpds_5k=exact_elpds_5k[-(1:L)]
    
    return(list(lastregime_pick=exact_elpds_1k, 
      last3regime_pick=exact_elpds_3k, 
      last5regime_pick=exact_elpds_5k,
      conv_problem=conv_problem))

    }else if(model=='HMM_a'){
    #stop("not defined")
    exact_elpds_1k <- numeric(nrow(data)) #loglik choosing a specific regime in a given year
    exact_elpds_3k  <- numeric(nrow(data))
    exact_elpds_5k  <- numeric(nrow(data))
  
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_hmm_tmb<- ricker_hmm_TMB(data=df_past,tv.par='a')
      conv_problem[i-(L-1)] <- fit_past_hmm_tmb$conv_problem
     
      alpha <- fit_past_hmm_tmb$alpha[fit_past_hmm_tmb$regime]
      beta <- fit_past_hmm_tmb$beta
      sigma <- fit_past_hmm_tmb$sigma
      #sigmaw <- sqrt((fit_past_hmm_tmb$sigma^2)%*%fit_past_hmm_tmb$probregime)

      rs_pred_1k=alpha[i]-beta*df_oos$S[i + 1]
      rs_pred_3k=mean(alpha[(i-2):i])-mean(beta)*df_oos$S[i + 1]
      rs_pred_5k=mean(alpha[(i-4):i])-mean(beta)*df_oos$S[i + 1]
      
      exact_elpds_1k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1k,sd=sigma))
      exact_elpds_3k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3k,sd=sigma))
      exact_elpds_5k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5k,sd=sigma))
        }
    exact_elpds_1k=exact_elpds_1k[-(1:L)]
    exact_elpds_3k=exact_elpds_3k[-(1:L)]
    exact_elpds_5k=exact_elpds_5k[-(1:L)]
    
    return(list(lastregime_pick=exact_elpds_1k, 
      last3regime_pick=exact_elpds_3k, 
      last5regime_pick=exact_elpds_5k,
      conv_problem=conv_problem))

  }else if(model=='HMM_b'){
    #stop("not defined")
    exact_elpds_1k <- numeric(nrow(data)) #loglik choosing a specific regime in a given year
    exact_elpds_3k  <- numeric(nrow(data))
    exact_elpds_5k  <- numeric(nrow(data))

    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_hmm_tmb<- ricker_hmm_TMB(data=df_past,tv.par='b')
      fit_past_hmm_tmb$conv_problem
     
      alpha <- fit_past_hmm_tmb$alpha
      beta <- fit_past_hmm_tmb$beta[fit_past_hmm_tmb$regime]
      sigma <- fit_past_hmm_tmb$sigma
      #sigmaw <- sqrt((fit_past_hmm_tmb$sigma^2)%*%fit_past_hmm_tmb$probregime)

      rs_pred_1k<-alpha-beta[i]*df_oos$S[i + 1]
      rs_pred_3k=mean(alpha)-mean(beta[(i-2):i])*df_oos$S[i + 1]
      rs_pred_5k=mean(alpha)-mean(beta[(i-4):i])*df_oos$S[i + 1]
    
      exact_elpds_1k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1k,sd=sigma))
      exact_elpds_3k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3k,sd=sigma))
      exact_elpds_5k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5k,sd=sigma))
  }
    exact_elpds_1k=exact_elpds_1k[-(1:L)]
    exact_elpds_3k=exact_elpds_3k[-(1:L)]
    exact_elpds_5k=exact_elpds_5k[-(1:L)]
    
    return(list(lastregime_pick=exact_elpds_1k, 
      last3regime_pick=exact_elpds_3k, 
      last5regime_pick=exact_elpds_5k,
      conv_problem=conv_problem))

  }else{
  	stop(paste("model", model,"not defined, valid options are c('static','rw_a','rw_b','rw_both', 'HMM', 'HMM_a','HMM_b')"))
  }
}

#' Generate AIC estimate from different quantiles of data for a set of TMB models
#'
#' @param data A data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param fit A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param model string. Which parameters are time-varying? Options are c('static','rw_a','rw_b',
#' 'rw_both', 'HMM', 'HMM_a', 'HMM_b')
#' @returns 
#' 
#' @importFrom stats nlminb 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' tmb_mod_lfo_cv(data=harck, model=c('static'))
#' 
tmb_aic=function(data,fit,type=c('full','d90','d80'),form=c('aic','bic'),k){

  L_mat=matrix(data=NA,nrow=length(fit),ncol=nrow(data))
  preds=matrix(data=NA,nrow=length(fit),ncol=nrow(data))
  preds[1,]=fit[[1]]$alpha-fit[[1]]$beta*data$S
  L_mat[1,]=log(dnorm(data$logRS,mean=preds[1,],sd=fit[[1]]$sig))
  preds[2,]=fit[[2]]$alpha-fit[[2]]$beta*data$S+fit[[2]]$residuals*fit[[2]]$rho
  L_mat[2,]=log(dnorm(data$logRS,mean=preds[2,],sd=fit[[2]]$sigar))
  preds[3,]=fit[[3]]$alpha-fit[[3]]$beta*data$S
  L_mat[3,]=log(dnorm(data$logRS,mean=preds[3,],sd=fit[[3]]$sig))
  preds[4,]=fit[[4]]$alpha-fit[[4]]$beta*data$S
  L_mat[4,]=log(dnorm(data$logRS,mean=preds[4,],sd=fit[[4]]$sig))
  preds[5,]=fit[[5]]$alpha-fit[[5]]$beta*data$S
  L_mat[5,]=log(dnorm(data$logRS,mean=preds[5,],sd=fit[[5]]$sig))
  preds[6,]=fit[[6]]$alpha[fit[[6]]$regime]-fit[[6]]$beta*data$S
  L_mat[6,]=log(dnorm(data$logRS,mean=preds[6,],sd=fit[[6]]$sigma))
  preds[7,]=fit[[7]]$alpha-fit[[7]]$beta[fit[[7]]$regime]*data$S
  L_mat[7,]=log(dnorm(data$logRS,mean=preds[7,],sd=fit[[7]]$sigma))
  preds[8,]=fit[[8]]$alpha[fit[[8]]$regime]-fit[[8]]$beta[fit[[8]]$regime]*data$S
  L_mat[8,]=log(dnorm(data$logRS,mean=preds[8,],sd=fit[[8]]$sigma))
 
if(type=='full'){
  LL=apply(L_mat,1,sum)
  AIC=-2*LL+2*k+(2*k^2+2*k)/(nrow(data)-k-1)
  BIC=-2*LL+k*log(nrow(data))
  
  dAIC=AIC-min(AIC)
  dBIC=BIC-min(BIC)
  w_aic=NA
  w_bic=NA
  for(i in 1:length(fit)){w_aic[i]=exp(-0.5*dAIC[i])/sum(exp(-0.5*dBIC))}
  for(i in 1:length(fit)){w_bic[i]=exp(-0.5*dBIC[i])/sum(exp(-0.5*dBIC))}
}
if(type=='d90'){
  L_mat=L_mat[,apply(L_mat,2,log_mean_exp)>=quantile(apply(L_mat,2,log_mean_exp),0.1)]
  
  LL=apply(L_mat,1,sum)
  AIC=-2*LL+2*k+(2*k^2+2*k)/(nrow(data)-k-1)
  BIC=-2*LL+k*log(nrow(data))
  
  dAIC=AIC-min(AIC)
  dBIC=BIC-min(BIC)
  w_aic=NA
  w_bic=NA
  for(i in 1:length(fit)){w_aic[i]=exp(-0.5*dAIC[i])/sum(exp(-0.5*dBIC))}
  for(i in 1:length(fit)){w_bic[i]=exp(-0.5*dBIC[i])/sum(exp(-0.5*dBIC))}
}
if(type=='d80'){
  L_mat=L_mat[,apply(L_mat,2,log_mean_exp)>=quantile(apply(L_mat,2,log_mean_exp),0.2)]
  
  LL=apply(L_mat,1,sum)
  AIC=-2*LL+2*k+(2*k^2+2*k)/(nrow(data)-k-1)
  BIC=-2*LL+k*log(nrow(data))
  
  dAIC=AIC-min(AIC)
  dBIC=BIC-min(BIC)
  w_aic=NA
  w_bic=NA
  for(i in 1:length(fit)){w_aic[i]=exp(-0.5*dAIC[i])/sum(exp(-0.5*dBIC))}
  for(i in 1:length(fit)){w_bic[i]=exp(-0.5*dBIC[i])/sum(exp(-0.5*dBIC))}
}
if(form=='aic'){
  return(w_aic)
}
if(form=='bic'){
  return(w_bic)
}
  
  
  }


