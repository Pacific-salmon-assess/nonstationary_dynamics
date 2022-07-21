rm(list=ls())
library(here);
library(cmdstanr)
here()
#HMM testing
fulldata<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_jun2022.csv'))

stock1<- subset(fulldata,stock.id==unique(stock.id)[18])

stock1$logRS=log(stock1$recruits/stock1$spawners)



####Stan models####
set_cmdstan_path()

file_hmm<- file.path(cmdstan_path(),'nonstationary dynamics',"hmm_test.stan")

mod_hmm<- cmdstan_model(file_hmm)

data=list(R_S = stock1$logRS,
          T=nrow(stock1),
          TT=as.numeric(factor(stock1$broodyear)),
          S=c((stock1$spawners)),
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
as.data.frame(test_hmm$summary(variables=c('zstar')))

