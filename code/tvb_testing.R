#Test and compare various Scripts for time varying b


#log(b) stan
#log(b) TMB
#log(b) TMB centralized
#(b) dlm

#Check estimates and convergence in real data 


source(here('code','dlm-wrapper.R'))
library(here)
library(TMB)
library(cmdstanr);
library(loo);
library(dlm)
library(ggplot2)
library(gridExtra)

sal_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_jun2022.csv'))
sal_info<- read.csv(here('data','filtered datasets','all_stocks_info_jun2022.csv'))

head(sal_dat)
sal_info[sal_info$stock.id==158,]
sal_dat[sal_dat$stock.id==158,]

sal_info <- subset(sal_info, stock.id %in% sal_dat$stock.id)

set_cmdstan_path()
tvbstan <- file.path(cmdstan_path(),'timevarmodels', "ricker_linear_varying_b.stan")
mod_tvbstan <- cmdstan_model(tvbstan)

#Model 4.2 tv  b
compile("TMBmodels/Ricker_tvlogb.cpp")
dyn.load(dynlib("TMBmodels/Ricker_tvlogb"))

#Model
#compile("TMBmodels/Ricker_tvlogbdev.cpp")
#dyn.load(dynlib("TMBmodels/Ricker_tvlogbdev"))

tmb_beta <- list()
stan_beta <- list()
dlm_beta <- list()

tmb_Smax <- list()
stan_Smax <- list()
dlm_Smax <- list()

tmb_conv <- list()
stan_conv <- list()
dlm_conv <- list()

spw_scaler<-rep(NA,length(unique(sal_dat$stock.id)))


for(i in seq_along(unique(sal_dat$stock.id))){
  print(i)
  #i<-2
  s <- subset(sal_dat, stock.id==unique(sal_dat$stock.id)[i])
  s$spawners[s$spawners==0] <- NA
  
  s$spawnersavg <- s$spawners#/10^(round(log10(min(s$spawners))-1))
  s$recruitsavg <- s$recruits#/10^(round(log10(min(s$spawners))-1))

  s$logR_S <- log(s$recruitsavg/s$spawnersavg)
  s$logR_S[is.infinite(s$logR_S)]<-NA
  SRdata<-list(obs_logRS=s$logR_S,obs_S=s$spawnersavg)
  SRdatadlm <- data.frame(byr=s$broodyear,
                       spwn=s$spawnersavg,
                       rec=s$recruits)

  srm <- lm(s$logR_S~ s$spawnersavg, na.action=na.omit)

  #Model tv logb TMB

  parametersb<- list(
    logbetao = log(ifelse(-srm$coefficients[[2]]<0,1e-08,-srm$coefficients[[2]])),
    alpha=3,
    logsigobs=log(.1),
    logsigb=log(.1),
    logbeta=log(rep(-srm$coefficients[2],length(s$spawners)))
    )
  
  objlogb <- MakeADFun(SRdata,parametersb,DLL="Ricker_tvlogb",random="logbeta")
  
  newtonOption( objlogb, smartsearch=FALSE)
  
  tryCatch({

    optlogb <- nlminb( objlogb$par, objlogb$fn, objlogb$gr)
    summary(sdreport(objlogb))
    p <- objlogb$par
    tmb_beta[[i]] <- objlogb$report()$beta
    tmb_Smax[[i]] <- objlogb$report()$Smax
    tmb_conv[[i]] <-  rep(optlogb$convergence,nrow(s))
  },  error = function(e) {
        print(paste("logb failed in iteration",i))
        tmb_conv[[i]] <<-  rep(1,nrow(s))
        tmb_beta[[i]] <<- rep(NA,nrow(s))
        tmb_Smax[[i]] <<- rep(NA,nrow(s))
  })

  
  #stan

  data=list(R_S = s$logR_S,
            N = nrow(s),
            TT = as.numeric(factor(s$broodyear)),
            S = c(s$spawnersavg))

  tryCatch({
    fit3 <- mod_tvbstan$sample(
      data = data,
      seed = 123, 
      init=0,
      chains = 6, 
      parallel_chains = 6,
      iter_warmup = 500,
      iter_sampling = 1000,
      refresh = 500,
      adapt_delta = 0.99,
      max_treedepth = 20 # print update every 500 iters
    )
  
    params<- fit3$draws(format='df',variables=c('log_a','b','log_b','sigma_b','sigma_e'))
  

    parssummary<-summary(params)
    stan_beta[[i]] <- parssummary[grep("^b",parssummary$variable),]$median
    stan_Smax[[i]] <- apply(1/params[,grep("^b",names(params))],2,median)
    stan_conv[[i]] <-  as.numeric(abs(1-parssummary[grep("^b",parssummary$variable),]$rhat)>.1)  
  },  error = function(e) {
        print(paste("logb failed in iteration",i))
        stan_conv[[i]] <<-  rep(1,nrow(s))
        stan_beta[[i]] <<- rep(NA,nrow(s))
        stan_Smax[[i]] <<- rep(NA,nrow(s))
  })

  #dlm

  tryCatch({
    bvary<-fitDLM(SRdatadlm, alpha_vary = FALSE, beta_vary = TRUE)
    dlm_beta[[i]] <- -bvary$results$beta
    dlm_Smax[[i]] <- 1/-bvary$results$beta
    dlm_conv[[i]] <-  rep(bvary$convergence,nrow(s))
  },  error = function(e) {
        print(paste("logb failed in iteration",i))
        dlm_conv[[i]] <<-  rep(1,nrow(s))
        dlm_beta[[i]] <<- rep(NA,nrow(s))
        dlm_Smax[[i]] <<- rep(NA,nrow(s))
  })



}


#plot estimates of b
# % of convergence for each model

save(dlm_Smax,dlm_beta,dlm_conv,
  stan_Smax,stan_beta,stan_conv,
  tmb_Smax,tmb_beta,tmb_conv,
 file = "C:/Users/worc/Documents/timevarproject/Rdta/tvbcomp.RData")




sum(unlist(lapply(tmb_conv,sum))==0)/length(tmb_conv) # 0.5121951


sum(unlist(lapply(stan_conv,sum))==0)/length(tmb_conv) #0.8878049

sum(unlist(lapply(dlm_conv,sum))==0)/length(tmb_conv) #0.9853659



length(unlist(dlm_Smax))
length(unlist(stan_Smax))
length(unlist(tmb_Smax))

length(unlist(dlm_conv))
nrow(sal_dat)

unlist(sapply(sapply(dlm_conv,length),seq_len))

df<- data.frame(smax=c(unlist(dlm_Smax),unlist(stan_Smax),unlist(tmb_Smax)),
  conv=c(unlist(dlm_conv),unlist(stan_conv),unlist(tmb_conv)),
  beta=c(unlist(dlm_beta),unlist(stan_beta),unlist(tmb_beta)),
  model=rep(c("dlm","stan","tmb"),each=length(unlist(dlm_Smax))),
  stock=unlist(mapply(rep,unique(sal_dat$stock.id),sapply(dlm_conv,length), USE.NAMES = F)),
  yr=unlist(sapply(sapply(dlm_conv,length),seq_len)))
  

df2<-df[df$conv==0,]

unique(df2$stock)

pps <- list()

for(n in seq_along(unique(df2$stock))){
  dff <- df2[df2$stock==unique(df2$stock)[n],]
  pp <- ggplot(dff) +
  geom_line(aes(x=yr, y=smax, color=model), size=2)+
  ggtitle(unique(dff$stock))+
  theme_bw(16)+
  scale_color_viridis_d()

  pps[[n]]<-pp
}


ggsave(
   filename = "../outputs/figures/bestcomp.pdf", 
   plot = marrangeGrob(pps, nrow=1, ncol=1),  
   width = 15, height = 9
)


unlist(sapply(dlm_conv,seq_along))

nrow(sal_dat)

sapply(dlm_conv,length)

(sal_info$stock.name, each=sapply(dlm_conv,length))


length(unlist(mapply(rep,sal_info$stock.name,sapply(dlm_conv,length))))

?rep

#try same thing without standardizing # spawners





