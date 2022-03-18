#Preliminary batched model support for varying dynamics in sockeye
library(here)
library(TMB)

sock_dat<- read.csv(here('data','filtered datasets','sockeye_final.csv'))
sock_info<- read.csv(here('data','sockeye','sockeye_info.csv'))

sock_info<- subset(sock_info, Stock.ID %in% sock_dat$stock.id)

#compile and load TMB models
#Model 1 static model
#not implemented yet
compile("TMBmodels/Ricker_simple.cpp")
dyn.load(dynlib("TMBmodels/Ricker_simple"))

#Model 2 tv a and static b (vary Srep)
compile("TMBmodels/Ricker_tva_Smax.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_Smax"))

#Model 3 tv a and static Srep (vary b)

compile("TMBmodels/Ricker_tva_Srep.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tva_Srep"))

#Model 4 tv  b
compile("TMBmodels/Ricker_tvb.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tvb"))


compile("TMBmodels/Ricker_tvbdois.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tvbdois"))

SmaxAIC<-NULL
SrepAIC<-NULL
simpleAIC<-NULL
bvaryAIC<-NULL



for(i in seq_len(nrow(sock_info))){
  s <- subset(sock_dat,stock.id==sock_info$Stock.ID[i])

  srm <- lm(s$logR_S~ s$spawners)
  
  SRdatas<-list(obs_logR=log(s$recruits),obs_S=s$spawners)
  
  #Model 1 - static a & b
  #to be implemented
  parameters_simple<- list(
    alpha=srm$coefficients[1],
    logbeta = log(ifelse(-srm$coefficients[2]<0,1e-08,-srm$coefficients[2])),
    logSigObs=log(.4)
    )

  obj_simple <- MakeADFun(SRdatas,parameters_simple,DLL="Ricker_simple")
  newtonOption(obj_simple, smartsearch=FALSE)

  opt_simple <- nlminb(obj_simple$par,obj_simple$fn,obj_simple$gr)
  #rep_simple <- obj_simple$report()

  simpleAIC[i]<-2*3-2*-opt_simple$objective


  SRdata<-list(obs_logR=log(s$recruits),obs_S=s$spawners, prbeta1=3,prbeta2=3)
  #Model 2 - tv a and static b
  parameters<- list(
    alphao=srm$coefficients[1],
    logSmax = log(ifelse(1/-srm$coefficients[2]<0,1/1e-08,1/-srm$coefficients[2])),
    rho=.2,
    logvarphi= 0.1,
    alpha=rep(0.9,length(s$recruits))
    )

  obj <- MakeADFun(SRdata,parameters,DLL="Ricker_tva_Smax",random="alpha")
   newtonOption(obj, smartsearch=FALSE)

  opt <- nlminb(obj$par,obj$fn,obj$gr)
  #rep <- obj$report()

  SmaxAIC[i]<-2*4-2*-opt$objective

  #Model 3 tv a and static Srep (vary b)

  parametersSrep<- list(
    alphao=srm$coefficients[1],
    logSrep = log(ifelse(srm$coefficients[1]/-srm$coefficients[2]<0,10000,srm$coefficients[1]/-srm$coefficients[2])),
    rho=.2,
    logvarphi= 0.1,
    alpha=rep(srm$coefficients[1],length(s$recruits))
    )


  objSrep <- MakeADFun(SRdata,parametersSrep,DLL="Ricker_tva_Srep",random="alpha")
  newtonOption(objSrep, smartsearch=FALSE)

  optSrep <- nlminb(objSrep$par,objSrep$fn,objSrep$gr)
  #repSrep <- objSrep$report()

  SrepAIC[i]<-2*4-2*-optSrep$objective

  #Model 4 tv b -- not working need to recode.

  SRdata<-list(obs_logR=log(s$recruits),obs_S=s$spawners, prbeta1=3,prbeta2=3)
  parametersb<- list(
    logbetao=-11,
    alpha=srm$coefficients[1],
    rho=.4,
    logvarphi= 0.5,
    logbeta=rep(-11,length(s$recruits))
    )

  objb <- MakeADFun(SRdata,parametersb,DLL="Ricker_tvb",random="logbeta")
  newtonOption(objb, smartsearch=FALSE)
  #objb$fn()
  #objb$gr()
  skip_to_next<-FALSE
  tryCatch(
    {optb <- nlminb(objb$par,objb$fn,objb$gr)},
    error =function(e){ skip_to_next <<- TRUE}) 
  
 
  if(skip_to_next) { 
    next 
  }else{
   bvaryAIC[i]<-2*4-2*-optb$objective
  }   

  #if(i>1 ){
  #  if(!is.na(bvaryAIC[i-1])){
  #    if(bvaryAIC[i-1]==2*4-2*-optb$objective){
  #      bvaryAIC[i]<-NA
  #    }else{
  #      bvaryAIC[i]<-2*4-2*-optb$objective
  #    }}
  #}else{
  #  bvaryAIC[i]<-2*4-2*-optb$objective
  #}
  
  
  
  

  #Model 4 tv b -- not working need to recode.
  #SRdatab<-list(obs_logRS=s$logR_S,obs_S=s$spawners)


  #parametersbvary<- list(
  #  logbetao=log(-srm$coefficients[[2]]),#log(ifelse(-srm$coefficients[2]<0, 1e-08,-srm$coefficients[2])),
  #  alpha=srm$coefficients[1],
  #  rho=.2,
  #  logvarphi= 0.1,
  #  logbeta_dev=rep(0,length(s$recruits)) #log(rep(ifelse(-srm$coefficients[2]<0, 1e-08,-srm$coefficients[2]),length(s$recruits)))
  #  )

  #objbvary <- MakeADFun(SRdatab,parametersbvary,DLL="Ricker_tvbdois",random="logbeta_dev")
  #newtonOption(objbvary, smartsearch=FALSE)
  #objbvary$fn()
  #objbvary$gr()

  #optbvary <- nlminb(objbvary$par,objbvary$fn,objbvary$gr)
  #repbvary <- objbvary$report()

  #bvaryAIC[i]<-2*4-2*-optbvary$objective


}



SmaxAIC-SrepAIC

deltaSmaxAIC<-SmaxAIC-pmin(SmaxAIC,SrepAIC,simpleAIC,bvaryAIC, na.rm=T)
deltaSrepAIC<-SrepAIC-pmin(SmaxAIC,SrepAIC,simpleAIC,bvaryAIC, na.rm=T)
deltasimpleAIC<-simpleAIC-pmin(SmaxAIC,SrepAIC,simpleAIC,bvaryAIC, na.rm=T)
deltabAIC<-bvaryAIC-pmin(SmaxAIC,SrepAIC,simpleAIC,bvaryAIC, na.rm=T)

sum(!is.na(deltaSmaxAIC)&deltaSmaxAIC==0)/length(deltaSmaxAIC)
sum(!is.na(deltaSrepAIC)&deltaSrepAIC==0)/length(deltaSrepAIC)
sum(!is.na(deltasimpleAIC)&deltasimpleAIC==0)/length(deltasimpleAIC)
sum(!is.na(deltabAIC)&deltabAIC==0)/length(deltabAIC)









  
  