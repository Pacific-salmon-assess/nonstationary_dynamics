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
    exact_elpds_wk <- numeric(nrow(data)) #loglik using weighted average of regimes
    exact_elpds_3k  <- numeric(nrow(data))
    exact_elpds_5k  <- numeric(nrow(data))
    exact_elpds_3wk <- numeric(nrow(data))
    exact_elpds_5wk <- numeric(nrow(data))

    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_hmm_tmb<- ricker_hmm_TMB(data=df_past,tv.par='both')
      conv_problem[i-(L-1)] <- fit_past_hmm_tmb$conv_problem
     
      alpha <- fit_past_hmm_tmb$alpha[fit_past_hmm_tmb$regime]
      alphaw <- fit_past_hmm_tmb$alpha%*%fit_past_hmm_tmb$probregime
      beta <- fit_past_hmm_tmb$beta[fit_past_hmm_tmb$regime]
      betaw <- fit_past_hmm_tmb$beta%*%fit_past_hmm_tmb$probregime
      sigma <- fit_past_hmm_tmb$sigma
      #sigmaw <- sqrt((fit_past_hmm_tmb$sigma^2)%*%fit_past_hmm_tmb$probregime)

      rs_pred_1k=alpha[i]-beta[i]*df_oos$S[i + 1]
      rs_pred_wk=alphaw[i]-betaw[i]*df_oos$S[i + 1]
      
      rs_pred_3k=mean(alpha[(i-2):i])-mean(beta[(i-2):i])*df_oos$S[i + 1]
      rs_pred_5k=mean(alpha[(i-4):i])-mean(beta[(i-4):i])*df_oos$S[i + 1]
      
      rs_pred_3wk=mean(alphaw[(i-2):i])-mean(betaw[(i-2):i])*df_oos$S[i + 1]
      rs_pred_5wk=mean(alphaw[(i-4):i])-mean(betaw[(i-4):i])*df_oos$S[i + 1]
      
      exact_elpds_1k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1k,sd=sigma))
      exact_elpds_wk[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_wk,sd=sigma))
      exact_elpds_3k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3k,sd=sigma))
      exact_elpds_5k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5k,sd=sigma))
      exact_elpds_3wk[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3wk,sd=sigma))
      exact_elpds_5wk[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5wk,sd=sigma))
    }
    exact_elpds_1k=exact_elpds_1k[-(1:L)]
    exact_elpds_wk=exact_elpds_wk[-(1:L)]
    exact_elpds_3k=exact_elpds_3k[-(1:L)]
    exact_elpds_5k=exact_elpds_5k[-(1:L)]
    exact_elpds_3wk=exact_elpds_3wk[-(1:L)]
    exact_elpds_5wk=exact_elpds_5wk[-(1:L)]
    
    return(list(lastregime_pick=exact_elpds_1k, 
      last3regime_pick=exact_elpds_3k, 
      last5regime_pick=exact_elpds_5k,
      lastregime_average=exact_elpds_wk,
      last3regime_average=exact_elpds_3wk,
      last5regime_average=exact_elpds_5wk,
      conv_problem=conv_problem))

    }else if(model=='HMM_a'){
    #stop("not defined")
    exact_elpds_1k <- numeric(nrow(data)) #loglik choosing a specific regime in a given year
    exact_elpds_wk <- numeric(nrow(data)) #loglik using weighted average of regimes
    exact_elpds_3k  <- numeric(nrow(data))
    exact_elpds_5k  <- numeric(nrow(data))
    exact_elpds_3wk <- numeric(nrow(data))
    exact_elpds_5wk <- numeric(nrow(data))

    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_hmm_tmb<- ricker_hmm_TMB(data=df_past,tv.par='a')
      conv_problem[i-(L-1)] <- fit_past_hmm_tmb$conv_problem
     
      alpha <- fit_past_hmm_tmb$alpha[fit_past_hmm_tmb$regime]
      alphaw <- fit_past_hmm_tmb$alpha%*%fit_past_hmm_tmb$probregime
      beta <- fit_past_hmm_tmb$beta
      betaw <- fit_past_hmm_tmb$beta
      sigma <- fit_past_hmm_tmb$sigma
      #sigmaw <- sqrt((fit_past_hmm_tmb$sigma^2)%*%fit_past_hmm_tmb$probregime)

      rs_pred_1k=alpha[i]-beta*df_oos$S[i + 1]
      rs_pred_wk=alphaw[i]-betaw*df_oos$S[i + 1]
      
      rs_pred_3k=mean(alpha[(i-2):i])-mean(beta)*df_oos$S[i + 1]
      rs_pred_5k=mean(alpha[(i-4):i])-mean(beta)*df_oos$S[i + 1]
      
      rs_pred_3wk=mean(alphaw[(i-2):i])-mean(betaw)*df_oos$S[i + 1]
      rs_pred_5wk=mean(alphaw[(i-4):i])-mean(betaw)*df_oos$S[i + 1]

      
      exact_elpds_1k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1k,sd=sigma))
      exact_elpds_wk[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_wk,sd=sigma))
      exact_elpds_3k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3k,sd=sigma))
      exact_elpds_5k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5k,sd=sigma))
      exact_elpds_3wk[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3wk,sd=sigma))
      exact_elpds_5wk[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5wk,sd=sigma))
    }
    exact_elpds_1k=exact_elpds_1k[-(1:L)]
    exact_elpds_wk=exact_elpds_wk[-(1:L)]
    exact_elpds_3k=exact_elpds_3k[-(1:L)]
    exact_elpds_5k=exact_elpds_5k[-(1:L)]
    exact_elpds_3wk=exact_elpds_3wk[-(1:L)]
    exact_elpds_5wk=exact_elpds_5wk[-(1:L)]
    
    return(list(lastregime_pick=exact_elpds_1k, 
      last3regime_pick=exact_elpds_3k, 
      last5regime_pick=exact_elpds_5k,
      lastregime_average=exact_elpds_wk,
      last3regime_average=exact_elpds_3wk,
      last5regime_average=exact_elpds_5wk,
      conv_problem=conv_problem))

  }else if(model=='HMM_b'){
    #stop("not defined")
    exact_elpds_1k <- numeric(nrow(data)) #loglik choosing a specific regime in a given year
    exact_elpds_wk <- numeric(nrow(data)) #loglik using weighted average of regimes
    exact_elpds_3k  <- numeric(nrow(data))
    exact_elpds_5k  <- numeric(nrow(data))
    exact_elpds_3wk <- numeric(nrow(data))
    exact_elpds_5wk <- numeric(nrow(data))


    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      fit_past_hmm_tmb<- ricker_hmm_TMB(data=df_past,tv.par='b')
      fit_past_hmm_tmb$conv_problem
     
      alpha <- fit_past_hmm_tmb$alpha
      alphaw <- fit_past_hmm_tmb$alpha
      beta <- fit_past_hmm_tmb$beta[fit_past_hmm_tmb$regime]
      betaw <- fit_past_hmm_tmb$beta%*%fit_past_hmm_tmb$probregime
      sigma <- fit_past_hmm_tmb$sigma
      #sigmaw <- sqrt((fit_past_hmm_tmb$sigma^2)%*%fit_past_hmm_tmb$probregime)

      rs_pred_1k<-alpha-beta[i]*df_oos$S[i + 1]
      rs_pred_wk<-alphaw-betaw[i]*df_oos$S[i + 1]

      rs_pred_3k=mean(alpha)-mean(beta[(i-2):i])*df_oos$S[i + 1]
      rs_pred_5k=mean(alpha)-mean(beta[(i-4):i])*df_oos$S[i + 1]
      
      rs_pred_3wk=mean(alphaw)-mean(betaw[(i-2):i])*df_oos$S[i + 1]
      rs_pred_5wk=mean(alphaw)-mean(betaw[(i-4):i])*df_oos$S[i + 1]
      
      
      exact_elpds_1k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1k,sd=sigma))
      exact_elpds_wk[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_wk,sd=sigma))
      exact_elpds_3k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3k,sd=sigma))
      exact_elpds_5k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5k,sd=sigma))
      exact_elpds_3wk[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3wk,sd=sigma))
      exact_elpds_5wk[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5wk,sd=sigma))
    }
    exact_elpds_1k=exact_elpds_1k[-(1:L)]
    exact_elpds_wk=exact_elpds_wk[-(1:L)]
    exact_elpds_3k=exact_elpds_3k[-(1:L)]
    exact_elpds_5k=exact_elpds_5k[-(1:L)]
    exact_elpds_3wk=exact_elpds_3wk[-(1:L)]
    exact_elpds_5wk=exact_elpds_5wk[-(1:L)]
    
    return(list(lastregime_pick=exact_elpds_1k, 
      last3regime_pick=exact_elpds_3k, 
      last5regime_pick=exact_elpds_5k,
      lastregime_average=exact_elpds_wk,
      last3regime_average=exact_elpds_3wk,
      last5regime_average=exact_elpds_5wk,
      conv_problem=conv_problem))

  }else{
  	stop(paste("model", model,"not defined, valid options are c('static','rw_a','rw_b','rw_both', 'HMM', 'HMM_a','HMM_b')"))
  }
}

