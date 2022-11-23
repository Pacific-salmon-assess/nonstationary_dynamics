#----------------------------------------------------
#Catarina put this functions together to be used in simulation
# 
#----------------------------------------------------







#' compile stan models
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param AC Logical. Are residuals autocorrelated? Default is FALSE
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
#' @param ... Anything else that would be passed to rstan::sampling
#' 
#' 
#' @returns a list containing the 
#' * alpha - median estimates for the alpha parameter vector
#' * beta - median estimates for the beta parameter 
#' * sigobs - median estimates for the observation error sigma         
#' * stanfit - a stanfit model object
#' * mcmcsummary - summary over kept samples
#' * c_mcmcsummary - chain specific summary 
#' * list of samples
#' 
#' 
#' @importFrom rstan stan extract summary 
#' 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' rickerstan(data=harck)
#' 
compile_code<-function(type=c('static','rw','hmm'), ac=FALSE, par=c('n','a','b','both'), caphigh=FALSE) {
 
   sm <- sr_mod(type=type, ac=ac, par=par, caphigh=caphigh, loglik=FALSE, modelcode=TRUE)
  
   mod <- stan_model(model_name="stanmod",model_code=sm)

   return(mod)
}

#' Simple Ricker model estimated with stan
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param AC Logical. Are residuals autocorrelated? Default is FALSE
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
#' @param ... Anything else that would be passed to rstan::sampling
#' 
#' 
#' @returns a list containing the 
#' * alpha - median estimates for the alpha parameter vector
#' * beta - median estimates for the beta parameter 
#' * sigobs - median estimates for the observation error sigma         
#' * stanfit - a stanfit model object
#' * mcmcsummary - summary over kept samples
#' * c_mcmcsummary - chain specific summary 
#' * list of samples
#' 
#' 
#' @importFrom rstan stan extract summary 
#' 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' rickerstan(data=harck)
#' 
ricker_stan <- function(data,  AC=FALSE, control = stancontrol(), mod=NULL, warmup=300,  chains = 6, iter = 1000,...) {
 

 if(is.null(mod)){
   sm <- compile_code(type='static',ac=AC,par='n',caphigh=FALSE)
  }else{
   sm <-mod
  }
  

  if(AC){
    datm = list(N=nrow(data),
                L=max(data$by)-min(data$by)+1,
                ii=seq_len(nrow(data)),
                R_S =data$logRS,
                S=data$S)
  }else{
    datm = list(N=nrow(data),
                R_S =data$logRS,
                S=data$S)
  

  }
  
  fit<-rstan::sampling(sm, data=datm,
                      control = control, warmup = warmup, 
                      chains = chains, iter = iter )
  #fit <- rstan::stan(model_code = sm, 
  #                      data = datm,
  #                      control = control, warmup = warmup, chains = chains, iter = iter)
  #
  
  mc <- rstan::extract(fit, 
                inc_warmup=FALSE, permuted=FALSE)
    
    
  aa <- rstan::summary(fit)
  


  return(list(alpha=aa$summary["log_a","50%"],
   beta=aa$summary["b","50%"],
   Smax=aa$summary["S_max","50%"],
   sigobs=aa$summary["sigma_e","50%"], 
   Smsy=aa$summary["S_msy","50%"],
   umsy=aa$summary["U_msy","50%"],
   stanfit=fit, 
   mcmcsummary=aa$summary,
   c_mcmcsummary=aa$c_summary, 
   samples=mc ) )

}








#' random walks Ricker model estimated with stan
#'
#' @param data A data frame containing Spawners (S) and log(Recruits/Spawners) (R_S) time series. Use sr_format for the correct column names. 
#' @param par Which parameter should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
#' @param ... Anything else that would be passed to rstan::sampling
#' @param sm_ext Default is null, external stan rw model. Implemented to speed up simulation evaluation, use with caution.
#' 
#' @returns a list containing the 
#' * alpha - median estimates for the alpha parameter vector
#' * beta - median estimates for the beta parameter 
#' * sigobs - median estimates for the observation error sigma         
#' * stanfit - a stanfit model object
#' * mcmcsummary - summary over kept samples
#' * c_mcmcsummary - chain specific summary 
#' * list of samples
#' 
#' 
#' @importFrom rstan stan extract summary 
#' 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' ricker_rw_stan(data=harck)
#' 
ricker_rw_stan <- function(data, par=c('a','b','both'),  control = stancontrol(), mod=NULL,
  warmup=300,  chains = 6, iter = 1000,...) {
  #par='b'

  if(is.null(mod)){
   sm <- compile_code(type='rw',ac=FALSE,par=par,caphigh=FALSE)
  }else{
   sm <- mod
  }
  #if(is.null(sm_ext)){
  #  sm <- sr_mod(type='rw',ac=FALSE,par=par,loglik=FALSE, modelcode=TRUE)
  #}else{
  #  sm <-sm_ext
  #}

  
  fit <- rstan::sampling(sm, data = list(N=nrow(data),
                                    L=max(data$by)-min(data$by)+1,
                                    ii=seq_along(data$logRS),
                                    R_S =data$logRS,
                                    S=data$S),
                        control = control, warmup = warmup, chains = chains, iter = iter)
  
  mc <- rstan::extract(fit, 
                inc_warmup=FALSE, permuted=FALSE)
    
    
  aa <- rstan::summary(fit)
  


  return(list(alpha=aa$summary[grep("log_a",row.names(aa$summary)),"50%"],
   beta=aa$summary[grep("^b",row.names(aa$summary)),"50%"],
   Smax=aa$summary[grep("S_max",row.names(aa$summary)),"50%"],
   sigobs=aa$summary["sigma_e","50%"],
   siga=ifelse(par=="a"|par=="both",aa$summary["sigma_a","50%"],NA),
   sigb=ifelse(par=="b"|par=="both",aa$summary["sigma_b","50%"],NA),
   Smsy=aa$summary[grep("S_msy",row.names(aa$summary)),"50%"],
   umsy=aa$summary[grep("U_msy",row.names(aa$summary)),"50%"],
   stanfit=fit, 
   mcmcsummary=aa$summary,
   c_mcmcsummary=aa$c_summary, 
   samples=mc ) )

}



#' random walks Ricker model estimated with stan
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param par Which parameter should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
#' @param ... Anything else that would be passed to rstan::sampling
#' @param sm_ext Default is null, external stan hmm model. Implemented to speed up simulation evaluation, use with caution.
#' 
#' @returns a list containing the 
#' * alpha - median estimates for the alpha parameter vector
#' * beta - median estimates for the beta parameter 
#' * sigobs - median estimates for the observation error sigma         
#' * stanfit - a stanfit model object
#' * mcmcsummary - summary over kept samples
#' * c_mcmcsummary - chain specific summary 
#' * list of samples 
#' 
#' 
#' @importFrom rstan stan extract summary 
#'
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' ricker_hmm_stan(data=harck)
#' 
ricker_hmm_stan <- function(data, par=c('a','b','both'), k_regime=2, 
  control = stancontrol(), warmup=300,  chains = 6, iter = 1000,..., mod=NULL) {
  #par='both'
  
  if(is.null(mod)){
   sm <- compile_code(type='hmm',ac=FALSE,par=par,caphigh=FALSE)
  }else{
   sm <-mod
  }
  #if(is.null(sm_ext)){
  #  sm <- sr_mod(type='hmm',ac=FALSE,par=par,loglik=FALSE, modelcode=TRUE)
  #}else{
  #  sm <- sm_ext
  #}

  
  fit <- rstan::sampling(sm, 
                        data = list(N=nrow(data),
                                    R_S =data$logRS,
                                    S=data$S,
                                    K=k_regime,
                                    alpha_dirichlet=rep(1,k_regime)
                                    ),
                        control = control, warmup = warmup, chains = chains, iter = iter)
  

  mc <- rstan::extract(fit, 
          inc_warmup=FALSE, permuted=FALSE)
    

  aa <- rstan::summary(fit)

  #extract time-series of parameters
  parts <-stan_regime_rps(m=fit,par=par)

#row.names(aa$summary)
#[grep("^b\\[",row.names(aa$summary))]

#aa$summary[grep("zstar",row.names(aa$summary)),"50%"]
# these are meaningless need to calc for each mcmc draw check with Dan
#need regime and weighted alpha, beta, sigma  
   return(list(

   alpha_regime=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$log_a_t,NA),
   alpha_wgt=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$log_a_wt,NA),
   Smax_regime=ifelse(rep(par=="b"|par=="both",nrow(data)),parts$S_max_t,NA),
   Smax_wgt=ifelse(rep(par=="b"|par=="both",nrow(data)),parts$S_max_wt,NA),
   Smsy_regime=parts$S_msy_t,
   Smsy_wgt=parts$S_msy_wt,
   umsy_regime=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$U_msy_t,NA),
   umsy_wgt=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$U_msy_wt,NA),
   alpha=ifelse(par=="a"|par=="both",list(aa$summary[grep("log_a\\[",row.names(aa$summary)),"50%"]),aa$summary["log_a","50%"])[[1]],
   beta=ifelse(par=="b"|par=="both",list(aa$summary[grep("b\\[",row.names(aa$summary)),"50%"]),aa$summary["b","50%"])[[1]],
   Smax=ifelse(par=="b"|par=="both",list(aa$summary[grep("S_max\\[",row.names(aa$summary)),"50%"]),aa$summary["S_max","50%"])[[1]],
   sigobs=aa$summary["sigma","50%"],
   pi=aa$summary[grep("pi1",row.names(aa$summary)),"50%"],
   A=aa$summary[grep("A",row.names(aa$summary)),"50%"],
   probregime =matrix(aa$summary[grep("gamma\\[",row.names(aa$summary)),"50%"],ncol=k_regime, byrow=T),
   regime = aa$summary[grep("^zstar",row.names(aa$summary)),"50%"],
   Smsy=aa$summary[grep("S_msy\\[",row.names(aa$summary)),"50%"],
   umsy=ifelse(par=="a"|par=="both",list(aa$summary[grep("U_msy\\[",row.names(aa$summary)),"50%"]),aa$summary["U_msy","50%"])[[1]],
   stanfit=fit, 
   mcmcsummary=aa$summary,
   c_mcmcsummary=aa$c_summary, 
   samples=mc ) )

}







#' Sampling control options. 
#'
#' Any arguments to pass to [rstan::sampling].
#'
#' @param eval.max Maximum number of evaluations of the objective function
#'   allowed.
#' @param iter.max Maximum number of iterations allowed.
#' @param ... Anything else. See the 'Control parameters' section of
#'   [rstan::sampling].
#'
#' @export
stancontrol <- function(adapt_delta = 0.99,  ...) {
  list(adapt_delta = 0.99, ...)
}



