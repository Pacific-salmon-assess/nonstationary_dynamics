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



log_sum_exp <- function(x) {
  max_x <- max(x)  
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}






# END
#***********************************************************************************
