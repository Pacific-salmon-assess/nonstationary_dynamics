# util functions
#===============================================
#Structure of the code copied from the sdmTMB package :
#https://github.com/pbs-assess/sdmTMB


#' Find order of magnitude
#'
#' @param x a number. 
#' @export
#' 
#' @returns order of magnitude
#' 
#' 
log10_ceiling <- function(x) {
    10^(ceiling(log10(x)))
}



#' Find lalpha intial value for logistic transformed TMB HMM model parameters
#'
#' @param U a number. upper limit for parameter 
#' @param L a number. lower limit for parameter 
#' @param initpar Inital value of parameter to be transformed 
#' 
#' @export
#' 
#' @returns logistic transformed parameter guess
#' 
#' 
find_linit <- function(U, L, initpar) {
    linitpar <- -log((U-L)/(initpar-L)-1)

    return( linitpar)
}

#'log sum of exponentials
#'
#' @param x a vector 
#' @export
#' 
#' @returns log sum of exponentials
#' 
#' 
log_sum_exp <- function(x) {
  max_x <- max(x)  
  max_x + log(sum(exp(x - max_x)))
}

#' mean of log-sum-exp
#'
#' @param x a vector 
#' @export
#' 
#' @returns mean of exponentials
#' 
#' 
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

#' Function to format the names for a stock-recruitment dataset to avoid compatibility problems with other samEst functions
#'
#' @param S a vector of spawners
#' @param R a vector of recruits
#' @param R a vector of years for each cohort
#' @export
#' 
#' @returns mean of exponentials
#' 
#' 
sr_format<- function(S,R,by){
  x=data.frame(S=S,R=R,by=by,logRS=log(R/S))
  return(x)
}




# END
#***********************************************************************************
