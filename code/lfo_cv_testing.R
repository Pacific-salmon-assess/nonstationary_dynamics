rm(list=ls())
library(here);library(brms);library(rstan);library(loo)
here()


##Functions for LFO-CV

log_sum_exp <- function(x) {
  max_x <- max(x)  
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

sum_log_ratios <- function(loglik, ids = NULL) {
  if (!is.null(ids)) loglik <- loglik[, ids, drop = FALSE]
  rowSums(loglik)
}

rbind_print <- function(...) {
  round(rbind(...), digits = 2)
}

mod_refit<- function(mod,newdata){
  rstan::stan(file =  here('code','stan models',mod), 
              data = list(R_S = newdata$y,
                          N=nrow(newdata),
                          TT=as.numeric(factor(newdata$year)),
                          S=c((newdata$spawners))),
              pars = c(),
              control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 300, chains = 6, iter = 800, thin = 1,init=0)
  
}

mod_refit2<- function(mod,newdata,oosdata,oos){
  rstan::stan(file =  here('code','stan models',mod), 
              data = list(R_S = newdata$y,
                          N=nrow(newdata),
                          TT=as.numeric(factor(newdata$year)),
                          S=c((newdata$spawners)),
                          N_oos=1,
                          y_oos=oosdata$y[oos],
                          x_oos=oosdata$spawners[oos]),
              pars = c(),
              control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 300, chains = 6, iter = 800, thin = 1,init=0)
  
}

##Function to get out-of-sample log likelihoods from models
log_lik_oos<- function(mod=1,params,newdata,oos){ 
  #mod = 1; static stock-recruitment model
  #mod = 2; time-varying productivity model
  #mod = 3; time-varying capacity model
  #mod = 4; time-varying productivity and capacity model
  
  #each model has its own pointwise likelihood generating model (eg. loglik_static, etc.)
  if(mod==1){ #static ricker model
    ll_mod<- rstan::stan(file =  here('code','stan models','loglik_static.stan'), 
                data = list(y_test = newdata$y,
                            x_test = newdata$spawners,
                            N=nrow(newdata),
                            N_samples=length(params$log_a),
                            alpha=params$log_a,
                            beta=params$b,
                            sigma=params$sigma_e),
                pars=c('log_lik'),
                chains = 1, iter = 1,algorithm="Fixed_param")
    #fits a stan model using the posterior estimates to generate likelihoods for all the data, including out of sample
    #does not sample parameters = 'Fixed_param') and only needs to be run once
    ll_param<- extract(ll_mod) #extract parameters
    ll<- ll_param$log_lik[,,] #out of sample log likelihood
    return(ll)
  }
  if(mod==2){ #TV alpha model
    ll_mod<- rstan::stan(file =  here('code','stan models',mod), 
                         data = list(y_test = newdata$y,
                                     x_test = newdata$spawners,
                                     N=nrow(newdata),
                                     N_samples=length(params$log_a),
                                     alpha=params$log_a,
                                     beta=params$b,
                                     sigma=params$sigma_e,
                                     S_oos=newdata$spawners[oos]),
                         pars = c(),
                         control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 300, chains = 6, iter = 800, thin = 1,init=0)
    ll<- loo::extract_log_lik(ll_mod)
    return(ll)
  }
  
}

#Data####
fulldata<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_jun2022.csv'))

stock1<- subset(fulldata,stock.id==unique(stock.id)[18])

stock1$logRS=log(stock1$recruits/stock1$spawners)


N <- length(stock1$logRS)
df <- data.frame(
  y = as.numeric(stock1$logRS),
  spawners=as.numeric(stock1$spawners),
  year = as.numeric(stock1$broodyear),
  time = 1:N
)



#Test with BRMS####

CHAINS <- 4
SEED <- 5838296
set.seed(SEED)

fit <- brm(
  y ~ spawners, 
  data = df,
  control = list(adapt_delta = 0.99), 
  seed = SEED, 
  chains = CHAINS
)

L=10
#
loglik_exact <- matrix(nrow = ndraws(fit), ncol = N)
for (i in L:(N - 1)) {
  past <- 1:i
  oos <- i + 1
  df_past <- df[past, , drop = FALSE]
  df_oos <- df[c(past, oos), , drop = FALSE]
  fit_i <- update(fit, newdata = df_past, recompile = FALSE)
  loglik_exact[, i + 1] <- log_lik(fit_i, newdata = df_oos, oos = oos)[, oos]
}

exact_elpds_1sap <- apply(loglik_exact, 2, log_mean_exp)
exact_elpd_1sap <- c(ELPD = sum(exact_elpds_1sap[-(1:L)]))


k_thres <- 0.7

approx_elpds_1sap <- rep(NA, N)

# initialize the process for i = L
past <- 1:L
oos <- L + 1
df_past <- df[past, , drop = FALSE]
df_oos <- df[c(past, oos), , drop = FALSE]
fit_past <- update(fit, newdata = df_past, recompile = FALSE)
loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
approx_elpds_1sap[L + 1] <- log_mean_exp(loglik[, oos])

# iterate over i > L
i_refit <- L
refits <- L
ks <- NULL
for (i in (L + 1):(N - 1)) {
  past <- 1:i
  oos <- i + 1
  df_past <- df[past, , drop = FALSE]
  df_oos <- df[c(past, oos), , drop = FALSE]
  loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
  
  logratio <- sum_log_ratios(loglik, (i_refit + 1):i)
  psis_obj <- suppressWarnings(psis(logratio))
  k <- pareto_k_values(psis_obj)
  ks <- c(ks, k)
  if (k > k_thres) {
    # refit the model based on the first i observations
    i_refit <- i
    refits <- c(refits, i)
    fit_past <- update(fit_past, newdata = df_past, recompile = FALSE)
    loglik <- log_lik(fit_past, newdata = df_oos, oos = oos)
    approx_elpds_1sap[i + 1] <- log_mean_exp(loglik[, oos])
  } else {
    lw <- weights(psis_obj, normalize = TRUE)[, 1]
    approx_elpds_1sap[i + 1] <- log_sum_exp(lw + loglik[, oos])
  }
} 

exact_elpds_1sap
approx_elpds_1sap
sum(na.omit(approx_elpds_1sap))
sum(na.omit(exact_elpds_1sap))

#Recreate with rstan####
L=10
past <- 1:L
oos <- L + 1
df_past <- df[past, , drop = FALSE]
df_oos <- df[c(past, oos), , drop = FALSE]
#full model
sr_static1<- rstan::stan(file = here('code','stan models','ricker_linear.stan'), data = list(R_S = df_past$y,
                                                                                             N=nrow(df_past),
                                                                                             TT=as.numeric(factor(df_past$year)),
                                                                                             S=c((df_past$spawners))),
                         pars = c('log_a','log_b','b','sigma_e'),
                         control = list(adapt_delta = 0.99,max_treedepth = 15), warmup = 200, chains = 4, iter = 800, thin = 1)

k_thres <- 0.7

approx_elpds_1sap_rstan<- rep(NA, N)

# initialize the process for i = L
past <- 1:L
oos <- L + 1
df_past <- df[past, , drop = FALSE]
df_oos <- df[c(past, oos), , drop = FALSE]

fit_past<- mod_refit(mod='ricker_linear.stan',newdata=df_past)

ll_oos<-log_lik_oos(mod=1,params=extract(fit_past),newdata = df_oos)

approx_elpds_1sap_rstan[L + 1] <- log_mean_exp(ll_oos[,oos])

i_refit <- L
refits <- L
ks <- NULL
for (i in (L+1):(N - 1)) {
  past <- 1:i
  oos <- i + 1
  df_past <- df[past, , drop = FALSE]
  df_oos <- df[c(past, oos), , drop = FALSE]
  loglik <- log_lik_oos(mod=1,params=extract(fit_past), newdata = df_oos)
  
  logratio <- sum_log_ratios(loglik, (i_refit + 1):i)
  psis_obj <- suppressWarnings(psis(logratio))
  k <- pareto_k_values(psis_obj)
  ks <- c(ks, k)
  if (k > k_thres) {
    # refit the model based on the first i observations
    i_refit <- i
    refits <- c(refits, i)
    fit_past<- mod_refit(mod='ricker_linear.stan',newdata=df_past)
    loglik <- log_lik_oos(mod=1,params=extract(fit_past), newdata = df_oos)
    approx_elpds_1sap_rstan[i + 1] <- log_mean_exp(loglik[, oos])
  } else {
    lw <- weights(psis_obj, normalize = TRUE)[, 1]
    approx_elpds_1sap_rstan[i + 1] <- log_sum_exp(lw + loglik[, oos])
  }
} 


loglik_exact_rs <- matrix(nrow = 3000, ncol = N)
for (i in L:(N - 1)) {
  past <- 1:i
  oos <- i + 1
  df_past <- df[past, , drop = FALSE]
  df_oos <- df[c(past, oos), , drop = FALSE]
  fit_i <- mod_refit2(mod='ricker_linear_oos.stan',newdata=df_past,oosdata=df_oos,oos=i+1)
  loglik_exact_rs[, i + 1] <- extract(fit_i)$log_lik_oos
}

exact_elpds_1sap_rs <- apply(loglik_exact_rs, 2, log_mean_exp)
exact_elpd_1sap_rs <- c(ELPD = sum(exact_elpds_1sap_rs[-(1:L)]))
exact_elpd_1sap_rs




