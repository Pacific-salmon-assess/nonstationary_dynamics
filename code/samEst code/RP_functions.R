#Function to calculate reference points once SR parameters are computed. 






#sgen functions from samSim



#' compute Sgen likelihood function
#'
#' @param S Spawner numbers set to the interval between 0 and Smsy in sGenSolver.
#' @param a alpha parameter in Ricker function: R=S*exp(a-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(a-b*S)
#' @param Smsy estimate of Smsy based on the alpha and beta parameters above.
#' 
#' 
#' 
#' @returns Sgen likelihood 
#' 
#' 
#' 
Sgencompute <- function(S, a,b, Smsy ) {
	#modified from samsim sGenOptimSmsyum
  
  prt <- S * exp(a - b * S)
  epsilon <- log(Smsy) - log(prt)
  nLogLike <- sum(dnorm(epsilon, 0, 0.01, log = T))
  return(list(nSS = sum(nLogLike)))
  
}





#' compute Sgen estimate
#'
#' @param a alpha parameter in Ricker function: R=S*exp(a-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(a-b*S)
#' @param Smsy estimate of Smsy based on the alpha and beta parameters above.
#' 
#' 
#' 
#' @returns Sgen estimate 
#' 
#' 
#' 
sGenSolver <- function(a,b, Smsy) {
  #gives the min Ricker log-likelihood
  if(a>0){
    fnSGen <- function(S, a, b, Smsy) -1.0 * Sgencompute(S, a, b, Smsy)$nSS
    fit <- optimize(f = fnSGen, interval = c(0, ((a / b) * (0.5 - 0.07 * a))),
                 a=a,b=b, Smsy = Smsy)
  }else{
    fit <-list(minimum=NA)
  }

  return(list(fit = fit$minimum))
}






#' compute Smsy estimate based on Scheuerell 2016. An explicit solution for 
#' calculating optimum spawning stock size from Ricker’s stock recruitment model
#'
#' @param a alpha parameter in Ricker function: R=S*exp(a-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(a-b*S)
#' @param Smsy estimate of Smsy based on the alpha and beta parameters above.
#' 
#' @importFrom gsl lambert_W0
#' 
#' @returns Smsy estimate 
#' 
#' 
#' 
smsySolver <- function(a,b) {
  #gives the min Ricker log-likelihood
  Smsy <- (1 - gsl::lambert_W0(exp(1 - a))) /b

  return(Smsy)
}



#' compute Umsy estimate based on Scheuerell 2016. An explicit solution for 
#' calculating optimum spawning stock size from Ricker’s stock recruitment model
#'
#' @param a alpha parameter in Ricker function: R=S*exp(a-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(a-b*S)
#' @param Smsy estimate of Smsy based on the alpha and beta parameters above.
#' 
#' @importFrom gsl lambert_W0
#' 
#' @returns Smsy estimate 
#' 
#' 
#' 
umsySolver <- function(a) {
  #gives the min Ricker log-likelihood
  umsy <- (1 - gsl::lambert_W0(exp(1 - a)))

  return(umsy)
}

#' Compute reference points from parameters of Hidden Markov Model (regime shift) Stan model
#'
#' @param m estimated fit for rstan model
#' @returns Estimates for Smax, Smsy, and Umsy over time, either unweighted (ie. returns parameters for the most likely regime sequence) or weighted (ie. probability of each regime x regime parameters) 
#' 
#' 
#' 
stan_regime_rps<- function(m,par=c('a','b','both')){
  d=extract(m)
  if(par=='a'){
    log_a=apply(d$log_a,2,median)
    beta=median(d$b)
    S_max=median(d$S_max)
    U_msy=apply(d$U_msy,2,median)
    S_msy=apply(d$S_msy,2,median)
    gamma=cbind(apply(d$gamma[,,1],2,median),apply(d$gamma[,,2],2,median))
    zstar=apply(d$zstar,2,median)
    
    #Most likely sequence of parameters:
    log_a_t=log_a[zstar]
    U_msy_t=U_msy[zstar]
    S_msy_t=S_msy[zstar]
    #Weighted parameters (probability of regime state X regime parameter)
    log_a_wt=gamma%*%log_a
    U_msy_wt = gamma%*%U_msy
    S_msy_wt = gamma%*%S_msy
    return(data.frame(log_a_t,log_a_wt,U_msy_t,U_msy_wt,S_msy_t,S_msy_wt))
  }
  if(par=='b'){
    log_a=median(d$log_a)
    beta=apply(d$b,2,median)
    S_max=apply(d$S_max,2,median)
    U_msy=median(d$U_msy)
    S_msy=apply(d$S_msy,2,median)
    gamma=cbind(apply(d$gamma[,,1],2,median),apply(d$gamma[,,2],2,median))
    zstar=apply(d$zstar,2,median)
    
    #Most likely sequence of parameters:
    beta_t=beta[zstar]
    S_max_t=S_max[zstar]
    S_msy_t=S_msy[zstar]
    #Weighted parameters (probability of regime state X regime parameter)
    beta_wt=gamma%*%beta
    S_max_wt = gamma%*%S_max
    S_msy_wt = gamma%*%S_msy
    return(data.frame(beta_t,beta_wt,S_max_t,S_max_wt,S_msy_t,S_msy_wt))
  }
  if(par=='both'){
    log_a=apply(d$log_a,2,median)
    beta=apply(d$b,2,median)
    S_max=apply(d$S_max,2,median)
    U_msy=apply(d$U_msy,2,median)
    S_msy=apply(d$S_msy,2,median)
    gamma=cbind(apply(d$gamma[,,1],2,median),apply(d$gamma[,,2],2,median))
    zstar=apply(d$zstar,2,median)
    
    #Most likely sequence of parameters:
    log_a_t=log_a[zstar]
    beta_t=beta[zstar]
    S_max_t=S_max[zstar]
    U_msy_t=U_msy[zstar]
    S_msy_t=S_msy[zstar]
    #Weighted parameters (probability of regime state X regime parameter)
    log_a_wt=gamma%*%log_a
    beta_wt=gamma%*%beta
    S_max_wt = gamma%*%S_max
    U_msy_wt = gamma%*%U_msy
    S_msy_wt = gamma%*%S_msy
    return(data.frame(log_a_t,log_a_wt,beta_t,beta_wt,S_max_t,S_max_wt,U_msy_t,U_msy_wt,S_msy_t,S_msy_wt))
  }
    
  
}


