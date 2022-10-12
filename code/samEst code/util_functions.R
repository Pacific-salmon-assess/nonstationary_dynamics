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










# END
#***********************************************************************************
