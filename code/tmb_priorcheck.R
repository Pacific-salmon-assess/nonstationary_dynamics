# Pre data prior check of TMB models
# Author: Catarina Wor
# Date April 2022
#======================================================

#Preliminary batched model support for varying dynamics in sockeye
library(here)
library(TMB)
library(tmbstan)
library("bayesplot")


#compile and load TMB models
#Model 1 static model
#not implemented yet
compile("TMBmodels/Ricker_simple_predata.cpp")
dyn.load(dynlib("TMBmodels/Ricker_simple_predata"))

#Model 2 tv a and static b (vary Srep)
compile("TMBmodels/Ricker_tva_Smax_predata.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_Smax_predata"))

#Model 3 tv a and static Srep (vary b)
compile("TMBmodels/Ricker_tva_Srep.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_Srep"))

#Model 4 tv  b
compile("TMBmodels/Ricker_tvb.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tvb"))


#Model 5 tv a and b
compile("TMBmodels/Ricker_tva_tvb.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_tvb"))


#compile("TMBmodels/Ricker_tvbdois.cpp")
#dyn.load(dynlib("TMBmodels/Ricker_tvbdois"))


#==================================
#load in data as placeholder
sock_dat<- read.csv(here('data','filtered datasets','sockeye_final.csv'))
sock_info<- read.csv(here('data','sockeye','sockeye_info.csv'))

sock_info<- subset(sock_info, Stock.ID %in% sock_dat$stock.id)
s <- subset(sock_dat,stock.id==sock_info$Stock.ID[2])



#==================================
#pre data checks


#simple
SRdata<-list(obs_logRS=s$logR_S,obs_S=s$spawners)
  
#Model 1 - static a & b
#to be implemented
parameters_simple<- list(
    alpha=1,
    logbeta = log(1e-08),
    logsigobs=log(.4)
    )

objsimple <- MakeADFun(SRdata,parameters_simple,DLL="Ricker_simple_predata")
fitmcmcsimple <- tmbstan(objsimple, chains=2,
              iter=10000, init="random",
              lower=c(-15,-50,-6),
               upper=c(15,0,6))

df <- reshape::melt(as.array(fitmcmcsimple))


pm <- ggplot(df) +
   geom_density(aes(x=value, color=chains)) +
   facet_wrap(~parameters, scales="free") +
   theme_bw(14)
pm

    
# tva

parameters<- list(
    alphao=2.7,
    logbeta = log(1e-08),
    logsigobs=log(.4),
    logsiga=log(.4),
    alpha=rep(2.0,length(s$recruits))
    )

obj <- MakeADFun(SRdata,parameters,DLL="Ricker_tva_Smax_predata",random="alpha")#,lower = -Inf, upper = Inf)
   newtonOption(obj, smartsearch=FALSE)
  #obj$fn()
  #obj$gr()

fitmcmctva <- tmbstan(obj, chains=2,
              iter=10000, init="random",
              lower=c(-15,-20,-6,-6),
               upper=c(15,0,6,6))

dftva <- reshape::melt(as.array(fitmcmctva))


pm <- ggplot(dftva) +
   geom_density(aes(x=value, color=chains)) +
   facet_wrap(~parameters, scales="free") +
   theme_bw(14)
pm

#these look wonky
mcmc_dens(fitmcmctva, pars = c("logsigobs", "logsiga"), transformations =list("logsigobs" = "exp", "logsiga" = "exp"))
mcmc_dens(fitmcmctva, pars = c("logbeta")) #+coord_cartesian(xlim=c(0.00000000001,0.01))







  
  