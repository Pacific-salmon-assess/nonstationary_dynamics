#Preliminary batched model support for varying dynamics in sockeye

library(here)
library(TMB)

sock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_jun2022.csv'))
sock_info<- read.csv(here('data','filtered datasets','all_stocks_info_jun2022.csv'))

head(sock_dat)
head(sock_info)


sock_info<- subset(sock_info, stock.id %in% sock_dat$stock.id)

#compile and load TMB models
#Model 1 static model
#not implemented yet
compile("TMBmodels/Ricker_simple.cpp")
dyn.load(dynlib("TMBmodels/Ricker_simple"))

#Model 2 tv a and static b (vary Srep)
compile("TMBmodels/Ricker_tva_Smax.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_Smax"))

#Model 3 tv a and static Srep (vary b)
#compile("TMBmodels/Ricker_tva_Srep.cpp")
#dyn.load(dynlib("TMBmodels/Ricker_tva_Srep"))

#Model 4 tv  b
#compile("TMBmodels/Ricker_tvb.cpp")
#dyn.load(dynlib("TMBmodels/Ricker_tvb"))

#Model 4.2 tv  b
compile("TMBmodels/Ricker_tvlogb.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tvlogb"))

#Model 5 tv a and b
compile("TMBmodels/Ricker_tva_tvb.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_tvb"))


#compile("TMBmodels/Ricker_tvbdois.cpp")
#dyn.load(dynlib("TMBmodels/Ricker_tvbdois"))

SmaxAIC<-NULL
#SrepAIC<-NULL
simpleAIC<-NULL
#bvaryAIC<-NULL
logbvaryAIC<-NULL
abvaryAIC<-NULL


for(i in seq_len(nrow(sock_info))){
  #i<-2
  s <- subset(sock_dat,stock.id==sock_info$stock.id[i])
  s$spawners[s$spawners==0]<-NA

  s$spawnersavg <- s$spawners#/1000# mean(s$spawners,na.rm=T)
  s$logR_S <- log(s$recruits/s$spawnersavg)
  #SRdata<-list(obs_logRS=s$logR_S,obs_S=s$spawnersavg)
  

  #s$spawnersavg <- s$spawners/mean(s$spawners,na.rm=T)
  #s$logR_S <- log(s$recruits/s$spawnersavg)

  srm <- lm(s$logR_S~ s$spawnersavg, na.action=na.omit)
  
  SRdata<-list(obs_logRS=s$logR_S,obs_S=s$spawnersavg)
  
  #Model 1 - static a & b
  parameters_simple<- list(
    alpha=srm$coefficients[[1]],
    logbeta = log(ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[[2]])),
    logsigobs=log(.4)
    )

  obj_simple <- MakeADFun(SRdata,parameters_simple,DLL="Ricker_simple")
  
  newtonOption(obj_simple, smartsearch=FALSE)

  opt_simple <- nlminb(obj_simple$par,obj_simple$fn,obj_simple$gr)
  #rep_simple <- obj_simple$report()

  simpleAIC[i] <- 2*3-2*-opt_simple$objective


  
  #Model 2 - tv a and static b

  parameters<- list(
    alphao=srm$coefficients[1],
    logbeta = log(ifelse(-srm$coefficients[2]<0,1e-08,-srm$coefficients[2])),
    #rho=.5,
    #logvarphi=0,
    logsigobs=log(.4),
    logsiga=log(.4),
    alpha=rep(srm$coefficients[1],length(s$recruits))
    )

  obj <- MakeADFun(SRdata,parameters,DLL="Ricker_tva_Smax",random="alpha")#,lower = -Inf, upper = Inf)
   newtonOption(obj, smartsearch=FALSE)
  #obj$fn()
  #obj$gr()

  opt <- nlminb(obj$par,obj$fn,obj$gr)
  #plot <- obj$report()

  SmaxAIC[i]<-2*4-2*-opt$objective

  #Model 3 tv a and static Srep (vary b)
  #parametersSrep<- list(
  #  alphao=srm$coefficients[1],
  #  logSrep = log(ifelse(srm$coefficients[1]/-srm$coefficients[2]<0,10000,srm$coefficients[1]/-srm$coefficients[2])),
  #  #rho=.2,
  #  #logvarphi= 0.1,
  #  logsigobs=log(.4),
  #  logsiga=log(.4),
  #  alpha=rep(srm$coefficients[1],length(s$recruits))
  #  )
  #
  #
  #objSrep <- MakeADFun(SRdata,parametersSrep,DLL="Ricker_tva_Srep",random="alpha")
  #newtonOption(objSrep, smartsearch=FALSE)
  ##objSrep$fn()
  ##objSrep$gr()
  #optSrep <- nlminb(objSrep$par,objSrep$fn,objSrep$gr)
  ##repSrep <- objSrep$report()

  #SrepAIC[i]<-2*4-2*-optSrep$objective

  #Model 4.2 tv logb
 
  parametersb<- list(
    logbetao = log(ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[[2]])),
    alpha=srm$coefficients[[1]],
    logsigobs=log(.4),
    logsigb=log(.4),
    #rho=.4,
    #logvarphi= 0.5,
    logbeta=log(rep(-srm$coefficients[2],length(s$spawners)))
    )

  objlogb <- MakeADFun(SRdata,parametersb,DLL="Ricker_tvlogb",random="logbeta")
  newtonOption( objlogb, smartsearch=FALSE)
  objlogb$fn()
  objlogb$gr()
  skip_to_next<-FALSE
  tryCatch(
    {optlogb <- nlminb( objlogb$par, objlogb$fn, objlogb$gr)},
    error =function(e){ skip_to_next <<- TRUE}) 
  
 
  if(skip_to_next) { 
    logbvaryAIC[i]<-NA
    next 
  }else{
   logbvaryAIC[i]<-2*4-2*-optlogb$objective
  }   


#Model 4 tv b 
# compile("TMBmodels/Ricker_tvb.cpp")
#dyn.load(dynlib("TMBmodels/Ricker_tvb"))


  #parametersb<- list(
  #  logbetao = log(ifelse(-srm$coefficients[2]<0,1e-08,-srm$coefficients[2])),
  #  alpha=srm$coefficients[1],
  #  logsigobs=log(.4),
  #  logsigb=log(.4),
  #  #rho=.4,
  #  #logvarphi= 0.5,
  #  logbeta=(rep(-srm$coefficients[2],length(s$spawners)))
  #  )
  #
  #objb <- MakeADFun(SRdata,parametersb,DLL="Ricker_tvb",random="logbeta",lower=c(-10,-20,-6,-6),
  #             upper=c(10,0,2,2))
  #newtonOption(objb, smartsearch=FALSE)
  ##  objb$env$inner.control$tol10 <- 0
  #objb$fn()
  #objb$gr()
  #optb <- nlminb(objb$par,objb$fn,objb$gr)
  #
  #skip_to_next<-FALSE
  #tryCatch(
  #  {optb <- nlminb(objb$par,objb$fn,objb$gr)},
  #  error =function(e){ skip_to_next <<- TRUE}) 
  #
  #if(skip_to_next) { 
  #  bvaryAIC[i]<-NA
  #  next 
  #}else{
  #  bvaryAIC[i]<-2*4-2*-optb$objective
  #}   


  #Model 5 tv a and b 
  #compile("TMBmodels/Ricker_tva_tvb.cpp")
  #dyn.load(dynlib("TMBmodels/Ricker_tva_tvb"))

  
  parametersab<- list(
    logbetao = log(ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[[2]])),
    alphao=srm$coefficients[[1]],
    logsigobs=log(.5),
    logsiga=log(.1),
    logsigb=log(.1),
    logbeta=rep(log(ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[[2]])),length(s$spawners)),
    alpha=rep(srm$coefficients[[1]],length(s$spawners))
    )

  objab <- MakeADFun(SRdata,parametersab,DLL="Ricker_tva_tvb",random=c("logbeta","alpha"))
  newtonOption(objab, smartsearch=FALSE)
  objab$fn()
  objab$gr()
  skip_to_next<-FALSE
  tryCatch(
    {optab <- nlminb(objab$par,objab$fn,objab$gr)},
    error =function(e){ skip_to_next <<- TRUE}) 
  
 
  if(skip_to_next) {
    abvaryAIC[i]<-NA 
    next 
  }else{
    abvaryAIC[i]<-2*5-2*-optab$objective
  }   

}


deltaSmaxAIC<-SmaxAIC-pmin(SmaxAIC,simpleAIC, logbvaryAIC,abvaryAIC, na.rm=T)
deltasimpleAIC<-simpleAIC-pmin(SmaxAIC,simpleAIC, logbvaryAIC,abvaryAIC, na.rm=T)
deltalogbAIC<-logbvaryAIC-pmin(SmaxAIC,simpleAIC, logbvaryAIC,abvaryAIC, na.rm=T)
deltaabAIC<-abvaryAIC-pmin(SmaxAIC,simpleAIC, logbvaryAIC,abvaryAIC, na.rm=T)

sum(!is.na(deltasimpleAIC)&deltasimpleAIC==0)/length(deltasimpleAIC)
sum(!is.na(deltaSmaxAIC)&deltaSmaxAIC==0)/length(deltaSmaxAIC)
sum(!is.na(deltalogbAIC)&deltabAIC==0)/length(deltabAIC)
sum(!is.na(deltaabAIC)&deltaabAIC==0)/length(deltaabAIC)



#deltaSmaxAIC<-SmaxAIC-pmin(SmaxAIC,SrepAIC,simpleAIC,bvaryAIC, logbvaryAIC,abvaryAIC, na.rm=T)
#deltaSrepAIC<-SrepAIC-pmin(SmaxAIC,SrepAIC,simpleAIC,bvaryAIC,logbvaryAIC,abvaryAIC, na.rm=T)
#deltasimpleAIC<-simpleAIC-pmin(SmaxAIC,SrepAIC,simpleAIC,bvaryAIC,logbvaryAIC,abvaryAIC, na.rm=T)
#deltabAIC<-bvaryAIC-pmin(SmaxAIC,SrepAIC,simpleAIC,bvaryAIC,logbvaryAIC,abvaryAIC, na.rm=T)
#deltalogbAIC<-logbvaryAIC-pmin(SmaxAIC,SrepAIC,simpleAIC,bvaryAIC,logbvaryAIC,abvaryAIC, na.rm=T)
#deltaabAIC<-abvaryAIC-pmin(SmaxAIC,SrepAIC,simpleAIC,bvaryAIC,logbvaryAIC,abvaryAIC, na.rm=T)

#sum(!is.na(deltasimpleAIC)&deltasimpleAIC==0)/length(deltasimpleAIC)
#sum(!is.na(deltaSmaxAIC)&deltaSmaxAIC==0)/length(deltaSmaxAIC)
#sum(!is.na(deltaSrepAIC)&deltaSrepAIC==0)/length(deltaSrepAIC)
#sum(!is.na(deltabAIC)&deltabAIC==0)/length(deltabAIC)
#sum(!is.na(deltaabAIC)&deltaabAIC==0)/length(deltaabAIC)










  
  