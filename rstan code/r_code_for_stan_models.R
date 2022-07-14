#Not sure what iterations/chains you were using for TMB (if MCMC?) - should probably make it consistent.

stan_rb=rstan::stan(file=here('stan code','ricker_linear_varying_a.stan'), data=list(R_S = s$logR_S,
                                                 N=nrow(s),
                                                 TT=as.numeric(factor(s$broodyear)),
                                                 S=c((s$spawners))),
                                                                                                      
            pars = c('log_a','b','sigma_e','sigma_a'),
            control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 500, chains = 4, iter = 2000, thin = 1)

params_stan_rb=rstan::extract(stan_rb)


stan_GP=rstan::stan(file=here('stan code','ricker_linear_varying_a_GP.stan'),data=list(R_S = s$logR_S,
                                                                             N=nrow(s),
                                                                             TT=as.numeric(factor(s$broodyear)),
                                                                             S=c((s$spawners))),
            pars = c('log_a','b','sigma_e','sigma_a'),
            control = list(adapt_delta = 0.95,max_treedepth = 15), warmup = 500, chains = 4, iter = 2000, thin = 1)

params_stan_rb=rstan::extract(stan_GP)