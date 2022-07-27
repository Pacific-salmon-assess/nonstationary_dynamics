

library(here)
library(TMB)
library(cmdstanr)

here()


#TMB model from Tang et al 2021
compile("TMBmodels/SR_HMM.cpp")
dyn.load(dynlib("TMBmodels/SR_HMM"))

#HMM testing

fulldata<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_jun2022.csv'))

stock1<- subset(fulldata,stock.id==unique(stock.id)[12])
stock1$spawnersavg<-stock1$spawners/mean(stock1$spawners,na.rm=T)
stock1$logRS=log(stock1$recruits/stock1$spawnersavg)

srm<-lm(stock1$logRS~stock1$spawnersavg)
plot(stock1$spawnersavg,stock1$logRS)


k_regime <- 2

tmb.data<-list(yt=stock1$logRS,
  st=stock1$spawnersavg, 
  alpha_u= 20,
  alpha_l=-0.5,
  beta_u=1,
  sigma_u = 2
  )

parameters <- list( 
         lalpha = -log(k_regime+1-1:k_regime),
         lbeta = rep(1/max(s$spawnersavg),k_regime),
         lsigma = rep(log(1),k_regime),
         pi1_tran = rep(0,k_regime-1),
         qij_tran = matrix(0,nrow=k_regime,ncol=k_regime-1)          
       ) 


obj<- MakeADFun(tmb.data,parameters,DLL="SR_HMM")
  

newtonOption(obj, smartsearch=FALSE)
opt <- nlminb(obj$par,obj$fn,obj$gr)
obj$gr(opt$par)
rep <- obj$report() 


HMM2_rt<-apply(rep$r_pred, 2,which.max) # give regime for each year


p <- tmb_obj$par
r <- tmb_obj$report(p)
####Stan models####
set_cmdstan_path()
#system(paste("cp", here("stanmodels","hmm_test.stan"), paste0(cmdstan_path(),"/timevarmodels")))

file_hmm<- file.path(cmdstan_path(),'timevarmodels',"hmm_test.stan")

mod_hmm<- cmdstan_model(file_hmm)

data=list(R_S = stock1$logRS,
          T=nrow(stock1),
          TT=as.numeric(factor(stock1$broodyear)),
          S=c((stock1$spawnersavg)),
          K=2,
          alpha_dirichlet=c(1,1))

test_hmm<- mod_hmm$sample(
  data = data,
  seed = 123, 
  chains = 6, 
  parallel_chains = 6,
  iter_warmup = 1000,
  iter_sampling = 2000,
  refresh = 500,
  adapt_delta = 0.99,
  max_treedepth = 20 # print update every 500 iters
)
test_hmm$print(max_rows=50)
as.data.frame(test_hmm$summary(variables=c('zstar')))$median
HMM2_rt

stansummary<-test_hmm$summary()

rep$alpha
stansummary[grep("a_max",stansummary$variable ),]$median