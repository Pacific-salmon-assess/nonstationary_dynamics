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
s <- sock_dat #subset(sock_dat,stock.id==sock_info$Stock.ID[2])



#==================================
#pre data checks


#simple
compile("TMBmodels/Ricker_simple_predata.cpp")
dyn.load(dynlib("TMBmodels/Ricker_simple_predata"))

SRdata<-list(obs_logRS=s$logR_S,obs_S=s$spawners)
  
#Model 1 - static a & b
#to be implemented
parameters_simple<- list(
    alpha=1,
    logbeta = log(1e-08),
    logsigobs=log(.4)
    )

objsimple <- MakeADFun(SRdata,parameters_simple,DLL="Ricker_simple_predata",lower=c(-10,-20,-6),
               upper=c(10,0,6))
fitmcmcsimple <- tmbstan(objsimple, chains=2,
              iter=10000, init="random",
              lower=c(-10,-20,-6),
               upper=c(10,0,6))

df <- reshape::melt(as.array(fitmcmcsimple))


pm <- ggplot(df) +
   geom_density(aes(x=value, color=chains)) +
   facet_wrap(~parameters, scales="free") +
   theme_bw(14)
pm



R_Projs <- matrix(nrow = 100, ncol = nrow(s))
for(i in 1:100){
  R_Projs[i, ] <- objsimple$simulate()$R_Proj
}

simR<-reshape::melt(R_Projs,varnames = c("ind","sim"))

recpj <- data.frame(recruits = c(simR$value,s$recruits),
  stock.id=s$stock.id,
  sim=c(rep("sim",nrow(simR)),rep("obs",nrow(s))),
  spawners=s$spawners)


pr <- ggplot(recpj) +
   geom_point(aes(x=spawners, y=recruits,color=sim),alpha=.6)+
   theme_bw(14)
pr


    
# tva

compile("TMBmodels/Ricker_tva_Smax_predata.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_Smax_predata"))


parameters<- list(
    alphao=2.7,
    logbeta = log(1e-08),
    logsigobs=log(.4),
    logsiga=log(.4),
    alpha=rep(2.0,length(s$recruits))
    )

obj <- MakeADFun(SRdata,parameters,DLL="Ricker_tva_Smax_predata",random="alpha", lower=c(-10,-20,-6,-6),
               upper=c(10,0,6,6))#,lower = -Inf, upper = Inf)
   newtonOption(obj, smartsearch=FALSE)
  #obj$fn()
  #obj$gr()

fitmcmctva <- tmbstan(obj, chains=2,
              iter=1000, init="random",
              lower=c(-10,-20,-6,-6),
               upper=c(10,0,6,6))

dftva <- reshape::melt(as.array(fitmcmctva))


pm <- ggplot(dftva) +
   geom_density(aes(x=value, color=chains)) +
   facet_wrap(~parameters, scales="free") +
   theme_bw(14)
pm




#these look wonky
mcmc_dens(fitmcmctva, pars = c("logsigobs", "logsiga"), transformations =list("logsigobs" = "exp", "logsiga" = "exp"))
mcmc_dens(fitmcmctva, pars = c("logbeta")) #+coord_cartesian(xlim=c(0.00000000001,0.01))

recpj_all<- list()

for(i in seq_len(nrow(sock_info))){

s <- subset(sock_dat,stock.id==sock_info$Stock.ID[i])

SRdata<-list(obs_logRS=s$logR_S,obs_S=s$spawners)
  

obj <- MakeADFun(SRdata,parameters,DLL="Ricker_tva_Smax_predata",random="alpha", lower=c(-10,-20,-6,-6),
               upper=c(10,0,2,2))#,lower = -Inf, upper = Inf)
   newtonOption(obj, smartsearch=FALSE)


R_Projs <- matrix(nrow = 100, ncol = nrow(s))
for(y in 1:100){
  R_Projs[y, ] <- obj$simulate()$R_Proj
}

simR<-reshape::melt(R_Projs,varnames = c("ind","sim"))

recpj_all[[i]] <- data.frame(recruits = c(simR$value,s$recruits),
  stock.id=s$stock.id,
  sim=c(rep("sim",nrow(simR)),rep("obs",nrow(s))),
  spawners=s$spawners)


}

recpj <- do.call(rbind.data.frame, recpj_all)


pr <- ggplot(recpj) +
   geom_point(aes(x=spawners, y=recruits,color=sim),alpha=.6)+
   theme_bw(14)
pr







  
  