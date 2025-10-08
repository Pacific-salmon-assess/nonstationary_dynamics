library(here);library(dplyr)
source(here('code','util_func.R'))


#models####
cmdstanr::set_cmdstan_path(path='C:/Users/greenbergda/Documents/.cmdstan/cmdstan-2.29.2')
library(cmdstanr)
options(mc.cores=6)

file=file.path(cmdstanr::cmdstan_path(),'sr models', "hier_static_ricker.stan")
m.static=cmdstanr::cmdstan_model(file)

file.tv=file.path(cmdstanr::cmdstan_path(),'sr models', "ind_tvalpha_ricker.stan")
m.tv2=cmdstanr::cmdstan_model(file.tv)

file.gam=file.path(cmdstanr::cmdstan_path(),'sr models', "gam_tvalpha_ricker.stan")
m.tv3=cmdstanr::cmdstan_model(file.gam)

file.gam2=file.path(cmdstanr::cmdstan_path(),'sr models', "hier_gam_tvalpha_ricker.stan")
m.tv32=cmdstanr::cmdstan_model(file.gam2)

file.gam3=file.path(cmdstanr::cmdstan_path(),'sr models', "gam_tvalpha_hier.stan")
m.tv33=cmdstanr::cmdstan_model(file.gam3)

#stock data ####
stk.dat<- read.csv(here('data','salmon_productivity_compilation2025-09-25.csv'))
stk.info=read.csv(here('data','stock_info2025-09-25.csv'))

#Static model set####
##chinook####
ck.info=subset(stk.info,species=='Chinook')
ck.dat=stk.dat[stk.dat$stock.id%in%ck.info$stock.id,]

N.s.ck=rag_n(ck.dat$stock.id)

smax.ck=ck.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

sr.ck=dplyr::distinct(ck.info,region,.keep_all=T)
sr.ck=sr.ck[order(sr.ck$region),] #dataframe of ordered subregions


dl1=list(N=nrow(ck.dat),
        J=length(unique(ck.dat$stock.id)),
        SR=length(unique(ck.info$region)),
        R=length(unique(ck.info$ocean.basin)),
        L=max(ck.dat$broodyear)-min(ck.dat$broodyear)+1,
        ii=ck.dat$broodyear-min(ck.dat$broodyear)+1,
        ij=as.numeric(factor(ck.dat$stock.id)),
        sr_or=as.numeric(factor(sr.ck$ocean.basin)),
        j_sr=as.numeric(factor(ck.info$region)),
        start_y=N.s.ck[,1],
        end_y=N.s.ck[,2],
        R_S=ck.dat$logRS,
        S=ck.dat$spawners,
        pSmax_mean=0.5*smax.ck$max.S,
        pSmax_sig=2*smax.ck$max.S)

stc.ck=     m.static$sample(data=dl1,
                  chains = 6, 
                  iter_warmup = 500,
                  iter_sampling =500,
                  refresh = 100,
                  adapt_delta = 0.99,
                  max_treedepth = 20)

stc.ck$save_object('./outputs/model fits/static_fit_chinook.RDS')

##chum####
cm.info=subset(stk.info,species=='Chum')
cm.dat=stk.dat[stk.dat$stock.id%in%cm.info$stock.id,]

N.s.cm=rag_n(cm.dat$stock.id)

smax.cm=cm.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

sr.cm=dplyr::distinct(cm.info,region,.keep_all=T)
sr.cm=sr.cm[order(sr.cm$region),] #dataframe of ordered subregions


dl1=list(N=nrow(cm.dat),
         J=length(unique(cm.dat$stock.id)),
         SR=length(unique(cm.info$region)),
         R=length(unique(cm.info$ocean.basin)),
         L=max(cm.dat$broodyear)-min(cm.dat$broodyear)+1,
         ii=cm.dat$broodyear-min(cm.dat$broodyear)+1,
         ij=as.numeric(factor(cm.dat$stock.id)),
         sr_or=as.numeric(factor(sr.cm$ocean.basin)),
         j_sr=as.numeric(factor(cm.info$region)),
         start_y=N.s.cm[,1],
         end_y=N.s.cm[,2],
         R_S=cm.dat$logRS,
         S=cm.dat$spawners,
         pSmax_mean=0.5*smax.cm$max.S,
         pSmax_sig=2*smax.cm$max.S)

stc.cm=     m.static$sample(data=dl1,
                            chains = 6, 
                            iter_warmup = 500,
                            iter_sampling =500,
                            refresh = 100,
                            adapt_delta = 0.995,
                            max_treedepth = 20)

stc.cm$save_object('./outputs/model fits/static_fit_chum.RDS')

##coho####
cho.info=subset(stk.info,species=='Coho')
cho.dat=stk.dat[stk.dat$stock.id%in%cho.info$stock.id,]

N.s.cho=rag_n(cho.dat$stock.id)

smax.cho=cho.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

sr.cho=dplyr::distinct(cho.info,region,.keep_all=T)
sr.cho=sr.cho[order(sr.cho$region),] #dataframe of ordered subregions


dl1=list(N=nrow(cho.dat),
         J=length(unique(cho.dat$stock.id)),
         SR=length(unique(cho.info$region)),
         R=length(unique(cho.info$ocean.basin)),
         L=max(cho.dat$broodyear)-min(cho.dat$broodyear)+1,
         ii=cho.dat$broodyear-min(cho.dat$broodyear)+1,
         ij=as.numeric(factor(cho.dat$stock.id)),
         sr_or=as.numeric(factor(sr.cho$ocean.basin)),
         j_sr=as.numeric(factor(cho.info$region)),
         start_y=N.s.cho[,1],
         end_y=N.s.cho[,2],
         R_S=cho.dat$logRS,
         S=cho.dat$spawners,
         pSmax_mean=0.5*smax.cho$max.S,
         pSmax_sig=2*smax.cho$max.S)

stc.cho=     m.static$sample(data=dl1,
                             chains = 6, 
                             iter_warmup = 500,
                             iter_sampling =500,
                             refresh = 100,
                             adapt_delta = 0.99,
                             max_treedepth = 20)

stc.cho$save_object('./outputs/model fits/static_fit_coho.RDS')

##Pinks####
###Even####
pke.info=subset(stk.info,species=='Pink-Even')
pke.dat=stk.dat[stk.dat$stock.id%in%pke.info$stock.id,]

N.s.pke=rag_n(pke.dat$stock.id)

smax.pke=pke.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

sr.pke=dplyr::distinct(pke.info,region,.keep_all=T)
sr.pke=sr.pke[order(sr.pke$region),] #dataframe of ordered subregions


dl1=list(N=nrow(pke.dat),
         J=length(unique(pke.dat$stock.id)),
         SR=length(unique(pke.info$region)),
         R=length(unique(pke.info$ocean.basin)),
         L=max(pke.dat$broodyear)-min(pke.dat$broodyear)+1,
         ii=pke.dat$broodyear-min(pke.dat$broodyear)+1,
         ij=as.numeric(factor(pke.dat$stock.id)),
         sr_or=as.numeric(factor(sr.pke$ocean.basin)),
         j_sr=as.numeric(factor(pke.info$region)),
         start_y=N.s.pke[,1],
         end_y=N.s.pke[,2],
         R_S=pke.dat$logRS,
         S=pke.dat$spawners,
         pSmax_mean=0.5*smax.pke$max.S,
         pSmax_sig=2*smax.pke$max.S)

stc.pke=     m.static$sample(data=dl1,
                             chains = 6, 
                             iter_warmup = 500,
                             iter_sampling =500,
                             refresh = 100,
                             adapt_delta = 0.999,
                             max_treedepth = 20)

stc.pke$save_object('./outputs/model fits/static_fit_pinke.RDS')

###Odd####
pko.info=subset(stk.info,species=='Pink-Odd')
pko.dat=stk.dat[stk.dat$stock.id%in%pko.info$stock.id,]

N.s.pko=rag_n(pko.dat$stock.id)

smax.pko=pko.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

sr.pko=dplyr::distinct(pko.info,region,.keep_all=T)
sr.pko=sr.pko[order(sr.pko$region),] #dataframe of ordered subregions


dl1=list(N=nrow(pko.dat),
         J=length(unique(pko.dat$stock.id)),
         SR=length(unique(pko.info$region)),
         R=length(unique(pko.info$ocean.basin)),
         L=max(pko.dat$broodyear)-min(pko.dat$broodyear)+1,
         ii=pko.dat$broodyear-min(pko.dat$broodyear)+1,
         ij=as.numeric(factor(pko.dat$stock.id)),
         sr_or=as.numeric(factor(sr.pko$ocean.basin)),
         j_sr=as.numeric(factor(pko.info$region)),
         start_y=N.s.pko[,1],
         end_y=N.s.pko[,2],
         R_S=pko.dat$logRS,
         S=pko.dat$spawners,
         pSmax_mean=0.5*smax.pko$max.S,
         pSmax_sig=2*smax.pko$max.S)

stc.pko=     m.static$sample(data=dl1,
                             chains = 6, 
                             iter_warmup = 500,
                             iter_sampling =500,
                             refresh = 100,
                             adapt_delta = 0.999,
                             max_treedepth = 20)

stc.pko$save_object('./outputs/model fits/static_fit_pinko.RDS')

##Sockeye####
scy.info=subset(stk.info,species=='Sockeye')
scy.dat=stk.dat[stk.dat$stock.id%in%scy.info$stock.id,]

N.s.scy=rag_n(scy.dat$stock.id)

smax.scy=scy.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

sr.scy=dplyr::distinct(scy.info,region,.keep_all=T)
sr.scy=sr.scy[order(sr.scy$region),] #dataframe of ordered subregions


dl1=list(N=nrow(scy.dat),
         J=length(unique(scy.dat$stock.id)),
         SR=length(unique(scy.info$region)),
         R=length(unique(scy.info$ocean.basin)),
         L=max(scy.dat$broodyear)-min(scy.dat$broodyear)+1,
         ii=scy.dat$broodyear-min(scy.dat$broodyear)+1,
         ij=as.numeric(factor(scy.dat$stock.id)),
         sr_or=as.numeric(factor(sr.scy$ocean.basin)),
         j_sr=as.numeric(factor(scy.info$region)),
         start_y=N.s.scy[,1],
         end_y=N.s.scy[,2],
         R_S=scy.dat$logRS,
         S=scy.dat$spawners,
         pSmax_mean=0.5*smax.scy$max.S,
         pSmax_sig=2*smax.scy$max.S)

stc.scy=     m.static$sample(data=dl1,
                             chains = 6, 
                             iter_warmup = 500,
                             iter_sampling =500,
                             refresh = 100,
                             adapt_delta = 0.995,
                             max_treedepth = 20)



stc.scy$save_object('./outputs/model fits/static_fit_sockeye.RDS')
#Time-varying model set####

#Chinook####
chk.info=subset(stk.info,species=='Chinook')
chk.dat=stk.dat[stk.dat$stock.id%in%chk.info$stock.id,]

smax.chk=chk.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

sr.chk=dplyr::distinct(chk.info,region,.keep_all=T)
sr.chk=sr.chk[order(sr.chk$region),] #dataframe of ordered subregions

bspline_year=splines::bs(chk.dat$broodyear,knots=seq(min(chk.dat$broodyear),max(chk.dat$broodyear),by=2),degree=3,intercept=F)
predspline_year=splines::bs(seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),knots=seq(min(chk.dat$broodyear),max(chk.dat$broodyear),by=2),degree=3,intercept=F)

dl2=list(N=nrow(chk.dat),
         J=length(unique(chk.dat$stock.id)),
         L=max(chk.dat$broodyear)-min(chk.dat$broodyear)+1,
         J_i=as.numeric(factor(chk.dat$stock.id)),
         sr_or=as.numeric(factor(sr.chk$ocean.basin)),
         j_or=as.numeric(factor(chk.info$ocean.basin)),
         j_sr=as.numeric(factor(chk.info$region)),
         R=length(unique(chk.info$ocean.basin)),
         SR=length(unique(chk.info$region)),
         R_S=chk.dat$logRS,
         S=chk.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.chk$max.S,
         pSmax_sig=2*smax.chk$max.S)

tv.chk_k2=     m.tv33$sample(data=dl2,
                            chains = 6, 
                            iter_warmup = 500,
                            iter_sampling =500,
                            refresh = 250,
                            adapt_delta = 0.999,
                            max_treedepth = 20)

tv.chk_k2$save_object('./outputs/model fits/ns_fit_chinook_k2.RDS')

log.a.t.chk=tv.chk_k2$draws(variable='log_a_t',format='draws_matrix')
m.log.a.t.chk=apply(log.a.t.chk,2,median)
l90.log.a.t.chk=apply(log.a.t.chk,2,quantile,0.025)
u90.log.a.t.chk=apply(log.a.t.chk,2,quantile,0.975)
log.a.jt.chk=tv.chk_k2$draws(variable='log_a_jt',format='draws_matrix')

plot(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.chk[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chk))],2,median)
  lines(logaj~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=chk.info$begin[j]-min(chk.dat$broodyear)+1
  end=chk.info$end[j]-min(chk.dat$broodyear)+1
  lines(logaj[st:end]~seq(chk.info$begin[j],chk.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,col='darkred')


log.a.r.chk=tv.chk_k2$draws(variable='log_a_rt',format='draws_matrix')
log.a.r.bs=log.a.r.chk[,grepl(',1]',colnames(log.a.r.chk))]
log.a.r.goa=log.a.r.chk[,grepl(',2]',colnames(log.a.r.chk))]
log.a.r.nc=log.a.r.chk[,grepl(',3]',colnames(log.a.r.chk))]
log.a.r.sc=log.a.r.chk[,grepl(',4]',colnames(log.a.r.chk))]

plot(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity log(alpha)',xlab='Brood Cohort Year',type='n',bty='l')
lines(apply(log.a.r.bs,2,median)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),col='goldenrod',lwd=2)
lines(apply(log.a.r.goa,2,median)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),col='darkgreen',lwd=2)
lines(apply(log.a.r.nc,2,median)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),col='navy',lwd=2)
lines(apply(log.a.r.sc,2,median)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),col='orange',lwd=2)


smax.chk=chk.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat.chk=make_design_matrix(chk.dat$spawners,grp=chk.dat$stock.id)

stk.year=expand.grid(levels(factor(chk.dat$stock.id)),seq(min(chk.dat$broodyear),max(chk.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
chk.dat$stk.year=match(paste(chk.dat$stock.id,chk.dat$broodyear,sep='.'),stk.year[,3])

dl.rw=list(N=nrow(chk.dat),
         J=length(unique(chk.dat$stock.id)),
         L=max(chk.dat$broodyear)-min(chk.dat$broodyear)+1,
         R=length(unique(chk.info$ocean.basin)),
         J_i=as.numeric(factor(chk.dat$stock.id)),
         J_ii=chk.dat$stk.year,
         J_or=as.numeric(factor(chk.info$ocean.basin)),
         R_S=chk.dat$logRS,
         S=S.mat.chk,
         pSmax_mean=0.5*smax.chk$max.S,
         pSmax_sig=2*smax.chk$max.S)

tv.chk_rw=     m.tv2$sample(data=dl.rw,
                             chains = 6, 
                             iter_warmup = 500,
                             iter_sampling =500,
                             refresh = 250,
                             adapt_delta = 0.999,
                             max_treedepth = 20)

tv.chk_rw$save_object('./outputs/model fits/ns_fit_chinook.RDS')

log.a.t.chk=tv.chk_rw$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.chk=apply(exp(log.a.t.chk),2,median)
l90.log.a.t.chk=apply(exp(log.a.t.chk),2,quantile,0.025)
u90.log.a.t.chk=apply(exp(log.a.t.chk),2,quantile,0.975)
log.a.jt.chk=tv.chk_rw$draws(variable='log_a_t',format='draws_matrix')

plot(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),ylim=c(0,12),lwd=3,col='navy',ylab='Productivity - max. R/S',xlab='Brood Cohort Year',type='n',bty='l',main='Chinook')
for(j in 1:dl.rw$J){
  logaj=apply(log.a.jt.chk[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chk))],2,median)
  lines(exp(logaj)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=3,col=adjustcolor('black',alpha.f=0.4))
  st=chk.info$begin[j]-min(chk.dat$broodyear)+1
  end=chk.info$end[j]-min(chk.dat$broodyear)+1
  lines(exp(logaj[st:end])~seq(chk.info$begin[j],chk.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lwd=4,col='navy')
lines(l90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,lwd=2,col='navy')
lines(u90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,lwd=2,col='navy')


log.a.r.chk=tv.chk_rw$draws(variable='mu_log_a_r',format='draws_matrix')
log.a.r.bs=log.a.r.chk[,grepl(',1]',colnames(log.a.r.chk))]
log.a.r.goa=log.a.r.chk[,grepl(',2]',colnames(log.a.r.chk))]
log.a.r.nc=log.a.r.chk[,grepl(',3]',colnames(log.a.r.chk))]
log.a.r.sc=log.a.r.chk[,grepl(',4]',colnames(log.a.r.chk))]

plot(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity log(alpha)',xlab='Brood Cohort Year',type='n',bty='l')
lines(apply(log.a.r.bs,2,median)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),col='goldenrod',lwd=2)
lines(apply(log.a.r.goa,2,median)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),col='darkgreen',lwd=2)
lines(apply(log.a.r.nc,2,median)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),col='navy',lwd=2)
lines(apply(log.a.r.sc,2,median)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),col='orange',lwd=2)

bspline_year=splines::bs(chk.dat$broodyear,knots=seq(min(chk.dat$broodyear),max(chk.dat$broodyear),by=3),degree=3,intercept=F)
predspline_year=splines::bs(seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),knots=seq(min(chk.dat$broodyear),max(chk.dat$broodyear),by=3),degree=3,intercept=F)

dl3=list(N=nrow(chk.dat),
         J=length(unique(chk.dat$stock.id)),
         L=max(chk.dat$broodyear)-min(chk.dat$broodyear)+1,
         J_i=as.numeric(factor(chk.dat$stock.id)),
         sr_or=as.numeric(factor(sr.chk$ocean.basin)),
         j_or=as.numeric(factor(chk.info$ocean.basin)),
         j_sr=as.numeric(factor(chk.info$region)),
         R=length(unique(chk.info$ocean.basin)),
         SR=length(unique(chk.info$region)),
         R_S=chk.dat$logRS,
         S=chk.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.chk$max.S,
         pSmax_sig=2*smax.chk$max.S)

tv.chk_k3=     m.tv33$sample(data=dl3,
                             chains = 6, 
                             iter_warmup = 500,
                             iter_sampling =500,
                             refresh = 250,
                             adapt_delta = 0.999,
                             max_treedepth = 20)


plot(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
lines(apply(log.a.r.bs,2,median)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),col='goldenrod',lwd=2)
lines(apply(log.a.r.goa,2,median)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),col='darkgreen',lwd=2)
lines(apply(log.a.r.nc,2,median)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),col='navy',lwd=2)
lines(apply(log.a.r.sc,2,median)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),col='orange',lwd=2)

bspline_year=splines::bs(chk.dat$broodyear,knots=seq(min(chk.dat$broodyear),max(chk.dat$broodyear),by=3),degree=5,intercept=F)
predspline_year=splines::bs(seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),knots=seq(min(chk.dat$broodyear),max(chk.dat$broodyear),by=3),degree=5,intercept=F)

sr.chk=dplyr::distinct(chk.info,region,.keep_all=T)
sr.chk=sr.chk[order(sr.chk$region),] #dataframe of ordered subregions

dl2=list(N=nrow(chk.dat),
         J=length(unique(chk.dat$stock.id)),
         L=max(chk.dat$broodyear)-min(chk.dat$broodyear)+1,
         J_i=as.numeric(factor(chk.dat$stock.id)),
         sr_or=as.numeric(factor(sr.chk$ocean.basin)),
         j_sr=as.numeric(factor(chk.info$region)),
         R=length(unique(chk.info$ocean.basin)),
         SR=length(unique(chk.info$region)),
         R_S=chk.dat$logRS,
         S=chk.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.chk$max.S,
         pSmax_sig=2*smax.chk$max.S)

tv.chk_k3=     m.tv3$sample(data=dl3,
                            chains = 5, 
                            init=0,
                            iter_warmup = 200,
                            iter_sampling =600,
                            refresh = 250,
                            adapt_delta = 0.99,
                            max_treedepth = 30)

log.a.t.chk=tv.chk_k3$draws(variable='log_a_t',format='draws_matrix')
m.log.a.t.chk=apply(exp(log.a.t.chk),2,median)
l90.log.a.t.chk=apply(exp(log.a.t.chk),2,quantile,0.025)
u90.log.a.t.chk=apply(exp(log.a.t.chk),2,quantile,0.975)
log.a.jt.chk=tv.chk_k3$draws(variable='log_a_jt',format='draws_matrix')

plot(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),ylim=c(0,20),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.chk[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chk))],2,median)
  lines(exp(logaj)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=chk.info$begin[j]-min(chk.dat$broodyear)+1
  end=chk.info$end[j]-min(chk.dat$broodyear)+1
  lines(exp(logaj[st:end])~seq(chk.info$begin[j],chk.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,col='darkred')

bspline_year=splines::bs(chk.dat$broodyear,knots=seq(min(chk.dat$broodyear),max(chk.dat$broodyear),by=4),degree=3,intercept=F)
predspline_year=splines::bs(seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),knots=seq(min(chk.dat$broodyear),max(chk.dat$broodyear),by=4),degree=3,intercept=F)

dl4=list(N=nrow(chk.dat),
         J=length(unique(chk.dat$stock.id)),
         L=max(chk.dat$broodyear)-min(chk.dat$broodyear)+1,
         J_i=as.numeric(factor(chk.dat$stock.id)),
         R_S=chk.dat$logRS,
         S=chk.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.chk$max.S,
         pSmax_sig=2*smax.chk$max.S)

tv.chk_k4=     m.tv3$sample(data=dl4,
                            chains = 5, 
                            init=0,
                            iter_warmup = 200,
                            iter_sampling =600,
                            refresh = 250,
                            adapt_delta = 0.99,
                            max_treedepth = 30)

log.a.t.chk=tv.chk_k4$draws(variable='log_a_t',format='draws_matrix')
m.log.a.t.chk=apply(exp(log.a.t.chk),2,median)
l90.log.a.t.chk=apply(exp(log.a.t.chk),2,quantile,0.025)
u90.log.a.t.chk=apply(exp(log.a.t.chk),2,quantile,0.975)
log.a.jt.chk=tv.chk_k4$draws(variable='log_a_jt',format='draws_matrix')

plot(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),ylim=c(0,20),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.chk[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chk))],2,median)
  lines(exp(logaj)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=chk.info$begin[j]-min(chk.dat$broodyear)+1
  end=chk.info$end[j]-min(chk.dat$broodyear)+1
  lines(exp(logaj[st:end])~seq(chk.info$begin[j],chk.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.6))
}
lines(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,lwd=2,col='darkred')
lines(u90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,lwd=2,col='darkred')


llk2=apply(tv.chk_k2$draws(variables='log_lik',format='draws_matrix'),1,sum)
llk3=apply(tv.chk_k3$draws(variables='log_lik',format='draws_matrix'),1,sum)
llk4=apply(tv.chk_k4$draws(variables='log_lik',format='draws_matrix'),1,sum)

#tv.chk_rw=readRDS('./outputs/model fits/ns_fit_chinook.RDS')

#Chum####
chm.info=subset(stk.info,species=='Chum')
chm.dat=stk.dat[stk.dat$stock.id%in%chm.info$stock.id,]

smax.chm=chm.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

bspline_year=splines::bs(chm.dat$broodyear,knots=seq(min(chm.dat$broodyear),max(chm.dat$broodyear),by=2),degree=3,intercept=F)
predspline_year=splines::bs(seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),knots=seq(min(chm.dat$broodyear),max(chm.dat$broodyear),by=2),degree=3,intercept=F)

sr.chm=dplyr::distinct(chm.info,region,.keep_all=T)
sr.chm=sr.chm[order(sr.chm$region),] #dataframe of ordered subregions

dl2=list(N=nrow(chm.dat),
         J=length(unique(chm.dat$stock.id)),
         L=max(chm.dat$broodyear)-min(chm.dat$broodyear)+1,
         J_i=as.numeric(factor(chm.dat$stock.id)),
         sr_or=as.numeric(factor(sr.chm$ocean.basin)),
         j_or=as.numeric(factor(chm.info$ocean.basin)),
         j_sr=as.numeric(factor(chm.info$region)),
         R=length(unique(chm.info$ocean.basin)),
         SR=length(unique(chm.info$region)),
         R_S=chm.dat$logRS,
         S=chm.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.chm$max.S,
         pSmax_sig=2*smax.chm$max.S)

tv.chm_k2=     m.tv33$sample(data=dl2,
                             chains = 5, 
                             iter_warmup = 500,
                             iter_sampling =500,
                             refresh = 250,
                             adapt_delta = 0.995,
                             max_treedepth = 20)

log.a.t.chm=tv.chm_k2$draws(variable='log_a_t',format='draws_matrix')
m.log.a.t.chm=apply(log.a.t.chm,2,median)
l90.log.a.t.chm=apply(log.a.t.chm,2,quantile,0.025)
u90.log.a.t.chm=apply(log.a.t.chm,2,quantile,0.975)
log.a.jt.chm=tv.chm_k2$draws(variable='log_a_jt',format='draws_matrix')

plot(m.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.chm[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chm))],2,median)
  lines(logaj~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=chm.info$begin[j]-min(chm.dat$broodyear)+1
  end=chm.info$end[j]-min(chm.dat$broodyear)+1
  lines(logaj[st:end]~seq(chm.info$begin[j],chm.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=5,col='darkred')

log.a.r.chm=tv.chm_k2$draws(variable='log_a_rt',format='draws_matrix')
log.a.r.bs=log.a.r.chm[,grepl(',1]',colnames(log.a.r.chm))]
log.a.r.goa=log.a.r.chm[,grepl(',2]',colnames(log.a.r.chm))]
log.a.r.nc=log.a.r.chm[,grepl(',3]',colnames(log.a.r.chm))]
log.a.r.sc=log.a.r.chm[,grepl(',4]',colnames(log.a.r.chm))]

plot(m.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l',main='Chum - Ocean regions')
lines(apply(log.a.r.bs,2,median)~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),col='goldenrod',lwd=2)
lines(apply(log.a.r.goa,2,median)~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),col='darkgreen',lwd=2)
lines(apply(log.a.r.nc,2,median)~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),col='navy',lwd=2)
lines(apply(log.a.r.sc,2,median)~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),col='orange',lwd=2)


smax.chm=chm.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat.chm=make_design_matrix(chm.dat$spawners,grp=chm.dat$stock.id)

stk.year=expand.grid(levels(factor(chm.dat$stock.id)),seq(min(chm.dat$broodyear),max(chm.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
chm.dat$stk.year=match(paste(chm.dat$stock.id,chm.dat$broodyear,sep='.'),stk.year[,3])

dl.rw=list(N=nrow(chm.dat),
           J=length(unique(chm.dat$stock.id)),
           L=max(chm.dat$broodyear)-min(chm.dat$broodyear)+1,
           R=length(unique(chm.info$ocean.basin)),
           J_i=as.numeric(factor(chm.dat$stock.id)),
           J_ii=chm.dat$stk.year,
           J_or=as.numeric(factor(chm.info$ocean.basin)),
           R_S=chm.dat$logRS,
           S=S.mat.chm,
           pSmax_mean=0.5*smax.chm$max.S,
           pSmax_sig=2*smax.chm$max.S)

tv.chm_rw=     m.tv2$sample(data=dl.rw,
                            chains = 6, 
                            iter_warmup = 500,
                            iter_sampling =500,
                            refresh = 250,
                            adapt_delta = 0.999,
                            max_treedepth = 20)

tv.chm_rw$save_object('./outputs/model fits/ns_fit_chum.RDS')

log.a.t.chm=tv.chm_rw$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.chm=apply(log.a.t.chm,2,median)
l90.log.a.t.chm=apply(log.a.t.chm,2,quantile,0.025)
u90.log.a.t.chm=apply(log.a.t.chm,2,quantile,0.975)
log.a.jt.chm=tv.chm_rw$draws(variable='log_a_t',format='draws_matrix')

plot(m.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity log(alpha)',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.chm[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chm))],2,median)
  lines(logaj~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=chm.info$begin[j]-min(chm.dat$broodyear)+1
  end=chm.info$end[j]-min(chm.dat$broodyear)+1
  lines(logaj[st:end]~seq(chm.info$begin[j],chm.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=5,col='darkred')


log.a.r.chm=tv.chm_rw$draws(variable='mu_log_a_r',format='draws_matrix')
log.a.r.bs=log.a.r.chm[,grepl(',1]',colnames(log.a.r.chm))]
log.a.r.goa=log.a.r.chm[,grepl(',2]',colnames(log.a.r.chm))]
log.a.r.nc=log.a.r.chm[,grepl(',3]',colnames(log.a.r.chm))]
log.a.r.sc=log.a.r.chm[,grepl(',4]',colnames(log.a.r.chm))]

plot(m.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity log(alpha)',xlab='Brood Cohort Year',type='n',bty='l',main='Chum - Ocean regions')
lines(apply(log.a.r.bs,2,median)~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),col='goldenrod',lwd=2)
lines(apply(log.a.r.goa,2,median)~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),col='darkgreen',lwd=2)
lines(apply(log.a.r.nc,2,median)~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),col='navy',lwd=2)
lines(apply(log.a.r.sc,2,median)~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),col='orange',lwd=2)


loglik.gam.chm=apply(tv.chm_k2$draws(variable='log_lik',format='draws_matrix'),1,sum)
loglik.rw.chm=apply(tv.chm_rw$draws(variable='log_lik',format='draws_matrix'),1,sum)
summary(loglik.gam.chm)
summary(loglik.rw.chm)


smax1=apply(tv.chm_k2$draws(variables='Smax',format='draws_matrix'),2,median)
smax2=apply(tv.chm_rw$draws(variables='Smax',format='draws_matrix'),2,median)
plot(smax1~smax2)
abline(c(0,1))
summary(smax1/smax2)

smsy1=apply(tv.chm_k2$draws(variables='Smsy_jt',format='draws_matrix'),2,median)
smsy2=apply(tv.chm_rw$draws(variables='Smsy',format='draws_matrix'),2,median)
plot(smsy1~smsy2)
abline(c(0,1))
summary(smsy1/smsy2)

logat1=apply(tv.chm_k2$draws(variables='log_a_jt',format='draws_matrix'),2,median)
logat2=apply(tv.chm_rw$draws(variables='log_a_t',format='draws_matrix'),2,median)
plot(logat1~logat2)
abline(c(0,1))
summary(logat1/logat2)


tv.chm_rw$save_object('./outputs/model fits/tv_fit_chum.RDS')

chm.info=subset(stk.info,species=='Chum')
chm.dat=stk.dat[stk.dat$stock.id%in%chm.info$stock.id,]

N.s.chm=rag_n(chm.dat$stock.id)

smax.chm=chm.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat.chm=make_design_matrix(chm.dat$spawners,grp=chm.dat$stock.id)

stk.year=expand.grid(levels(factor(chm.dat$stock.id)),seq(min(chm.dat$broodyear),max(chm.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
chm.dat$stk.year=match(paste(chm.dat$stock.id,chm.dat$broodyear,sep='.'),stk.year[,3])

dl2=list(N=nrow(chm.dat),
         J=length(unique(chm.dat$stock.id)),
         L=max(chm.dat$broodyear)-min(chm.dat$broodyear)+1,
         J_i=as.numeric(factor(chm.dat$stock.id)),
         J_ii=chm.dat$stk.year,
         start_y=N.s.chm[,1],
         end_y=N.s.chm[,2],
         R_S=chm.dat$logRS,
         S=S.mat.chm,
         pSmax_mean=0.5*smax.chm$max.S,
         pSmax_sig=2*smax.chm$max.S)

tv.chm=     m.tv2$sample(data=dl2,
                         chains = 6, 
                         iter_warmup = 500,
                         iter_sampling =1000,
                         refresh = 250,
                         adapt_delta = 0.999,
                         max_treedepth = 20)

log.a.t.chm=tv.chm$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.chm=apply(exp(log.a.t.chm),2,median)
l90.log.a.t.chm=apply(exp(log.a.t.chm),2,quantile,0.025)
u90.log.a.t.chm=apply(exp(log.a.t.chm),2,quantile,0.975)
log.a.jt.chm=tv.chm$draws(variable='log_a_t',format='draws_matrix')

plot(m.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),ylim=c(0,20),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.chm[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chm))],2,median)
  lines(exp(logaj)~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=chm.info$begin[j]-min(chm.dat$broodyear)+1
  end=chm.info$end[j]-min(chm.dat$broodyear)+1
  lines(exp(logaj[st:end])~seq(chm.info$begin[j],chm.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=5,col='darkred')



#Coho####
cho.info=subset(stk.info,species=='Coho')
cho.dat=stk.dat[stk.dat$stock.id%in%cho.info$stock.id,]

smax.cho=cho.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

bspline_year=splines::bs(cho.dat$broodyear,knots=seq(min(cho.dat$broodyear),max(cho.dat$broodyear),by=2),degree=3,intercept=F)
predspline_year=splines::bs(seq(min(cho.dat$broodyear),max(cho.dat$broodyear)),knots=seq(min(cho.dat$broodyear),max(cho.dat$broodyear),by=2),degree=3,intercept=F)

sr.cho=dplyr::distinct(cho.info,region,.keep_all=T)
sr.cho=sr.cho[order(sr.cho$region),] #dataframe of ordered subregions

dl2=list(N=nrow(cho.dat),
         J=length(unique(cho.dat$stock.id)),
         L=max(cho.dat$broodyear)-min(cho.dat$broodyear)+1,
         J_i=as.numeric(factor(cho.dat$stock.id)),
         sr_or=as.numeric(factor(sr.cho$ocean.basin)),
         j_sr=as.numeric(factor(cho.info$region)),
         R=length(unique(cho.info$ocean.basin)),
         SR=length(unique(cho.info$region)),
         R_S=cho.dat$logRS,
         S=cho.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.cho$max.S,
         pSmax_sig=2*smax.cho$max.S)

tv.cho_k2=     m.tv33$sample(data=dl2,
                             chains = 5, 
                             init=0,
                             iter_warmup = 200,
                             iter_sampling =600,
                             refresh = 250,
                             adapt_delta = 0.999,
                             max_treedepth = 30)

log.a.t.cho=tv.cho_k2$draws(variable='log_a_t',format='draws_matrix')
m.log.a.t.cho=apply(exp(log.a.t.cho),2,median)
l90.log.a.t.cho=apply(exp(log.a.t.cho),2,quantile,0.025)
u90.log.a.t.cho=apply(exp(log.a.t.cho),2,quantile,0.975)
log.a.jt.cho=tv.cho_k2$draws(variable='log_a_jt',format='draws_matrix')

plot(m.log.a.t.cho~seq(min(cho.dat$broodyear),max(cho.dat$broodyear)),ylim=c(0,20),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.cho[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.cho))],2,median)
  lines(exp(logaj)~seq(min(cho.dat$broodyear),max(cho.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=cho.info$begin[j]-min(cho.dat$broodyear)+1
  end=cho.info$end[j]-min(cho.dat$broodyear)+1
  lines(exp(logaj[st:end])~seq(cho.info$begin[j],cho.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.cho~seq(min(cho.dat$broodyear),max(cho.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.cho~seq(min(cho.dat$broodyear),max(cho.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.cho~seq(min(cho.dat$broodyear),max(cho.dat$broodyear)),lty=5,col='darkred')

cho.info=subset(stk.info,species=='Coho')
cho.dat=stk.dat[stk.dat$stock.id%in%cho.info$stock.id,]

N.s.cho=rag_n(cho.dat$stock.id)

smax.cho=cho.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat.cho=make_design_matrix(cho.dat$spawners,grp=cho.dat$stock.id)

stk.year=expand.grid(levels(factor(cho.dat$stock.id)),seq(min(cho.dat$broodyear),max(cho.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
cho.dat$stk.year=match(paste(cho.dat$stock.id,cho.dat$broodyear,sep='.'),stk.year[,3])

dl.rw=list(N=nrow(cho.dat),
           J=length(unique(cho.dat$stock.id)),
           L=max(cho.dat$broodyear)-min(cho.dat$broodyear)+1,
           R=length(unique(cho.info$ocean.basin)),
           J_i=as.numeric(factor(cho.dat$stock.id)),
           J_ii=cho.dat$stk.year,
           J_or=as.numeric(factor(cho.info$ocean.basin)),
           R_S=cho.dat$logRS,
           S=S.mat.cho,
           pSmax_mean=0.5*smax.cho$max.S,
           pSmax_sig=2*smax.cho$max.S)

tv.cho_rw=     m.tv2$sample(data=dl.rw,
                            chains = 6, 
                            iter_warmup = 500,
                            iter_sampling =500,
                            refresh = 250,
                            adapt_delta = 0.999,
                            max_treedepth = 20)

tv.cho_rw$save_object('./outputs/model fits/ns_fit_coho.RDS')

log.a.t.cho=tv.cho$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.cho=apply(exp(log.a.t.cho),2,median)
l90.log.a.t.cho=apply(exp(log.a.t.cho),2,quantile,0.025)
u90.log.a.t.cho=apply(exp(log.a.t.cho),2,quantile,0.975)
log.a.jt.cho=tv.cho$draws(variable='log_a_t',format='draws_matrix')

plot(m.log.a.t.cho~seq(min(cho.dat$broodyear),max(cho.dat$broodyear)),ylim=c(0,20),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.cho[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.cho))],2,median)
  lines(exp(logaj)~seq(min(cho.dat$broodyear),max(cho.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=cho.info$begin[j]-min(cho.dat$broodyear)+1
  end=cho.info$end[j]-min(cho.dat$broodyear)+1
  lines(exp(logaj[st:end])~seq(cho.info$begin[j],cho.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.cho~seq(min(cho.dat$broodyear),max(cho.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.cho~seq(min(cho.dat$broodyear),max(cho.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.cho~seq(min(cho.dat$broodyear),max(cho.dat$broodyear)),lty=5,col='darkred')

ll1=apply(tv.cho_k2$draws(variables='log_lik',format='draws_matrix'),1,sum)
ll2=apply(tv.cho$draws(variables='log_lik',format='draws_matrix'),1,sum)

smax1=apply(tv.cho_k2$draws(variables='Smax',format='draws_matrix'),2,median)
smax2=apply(tv.cho$draws(variables='Smax',format='draws_matrix'),2,median)
plot(smax1~smax2)
abline(c(0,1))
summary(smax1/smax2)

smsy1=apply(tv.cho_k2$draws(variables='Smsy_jt',format='draws_matrix'),2,median)
smsy2=apply(tv.cho$draws(variables='Smsy',format='draws_matrix'),2,median)
plot(smsy1~smsy2)
abline(c(0,1))
summary(smsy1/smsy2)

logat1=apply(tv.cho_k2$draws(variables='log_a_jt',format='draws_matrix'),2,median)
logat2=apply(tv.cho$draws(variables='log_a_t',format='draws_matrix'),2,median)
plot(logat1~logat2)
abline(c(0,1))
summary(logat1/logat2)



tv.cho$save_object('./outputs/model fits/tv_fit_coho.RDS')

#Pink-Even####
pke.info=subset(stk.info,species=='Pink-Even')
pke.dat=stk.dat[stk.dat$stock.id%in%pke.info$stock.id,]

smax.pke=pke.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

bspline_year=splines::bs(pke.dat$broodyear,knots=seq(min(pke.dat$broodyear),max(pke.dat$broodyear),by=2),degree=3,intercept=F)
predspline_year=splines::bs(seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),knots=seq(min(pke.dat$broodyear),max(pke.dat$broodyear),by=2),degree=3,intercept=F)

sr.pke=dplyr::distinct(pke.info,region,.keep_all=T)
sr.pke=sr.pke[order(sr.pke$region),] #dataframe of ordered subregions

dl2=list(N=nrow(pke.dat),
         J=length(unique(pke.dat$stock.id)),
         L=max(pke.dat$broodyear)-min(pke.dat$broodyear)+1,
         J_i=as.numeric(factor(pke.dat$stock.id)),
         sr_or=as.numeric(factor(sr.pke$ocean.basin)),
         j_or=as.numeric(factor(pke.info$ocean.basin)),
         j_sr=as.numeric(factor(pke.info$region)),
         R=length(unique(pke.info$ocean.basin)),
         SR=length(unique(pke.info$region)),
         R_S=pke.dat$logRS,
         S=pke.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.pke$max.S,
         pSmax_sig=2*smax.pke$max.S)

tv.pke_k2=     m.tv33$sample(data=dl2,
                             chains = 5, 
                             iter_warmup = 500,
                             iter_sampling =500,
                             refresh = 250,
                             adapt_delta = 0.995,
                             max_treedepth = 20)

log.a.t.pke=tv.pke_k2$draws(variable='log_a_t',format='draws_matrix')
m.log.a.t.pke=apply(log.a.t.pke,2,median)
l90.log.a.t.pke=apply(log.a.t.pke,2,quantile,0.025)
u90.log.a.t.pke=apply(log.a.t.pke,2,quantile,0.975)
log.a.jt.pke=tv.pke_k2$draws(variable='log_a_jt',format='draws_matrix')

plot(m.log.a.t.pke~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.pke[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.pke))],2,median)
  lines(logaj~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=pke.info$begin[j]-min(pke.dat$broodyear)+1
  end=pke.info$end[j]-min(pke.dat$broodyear)+1
  lines(logaj[st:end]~seq(pke.info$begin[j],pke.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.pke~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.pke~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.pke~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),lty=5,col='darkred')

log.a.r.pke=tv.pke_k2$draws(variable='log_a_rt',format='draws_matrix')
log.a.r.bs=log.a.r.pke[,grepl(',1]',colnames(log.a.r.pke))]
log.a.r.goa=log.a.r.pke[,grepl(',2]',colnames(log.a.r.pke))]
log.a.r.nc=log.a.r.pke[,grepl(',3]',colnames(log.a.r.pke))]
log.a.r.sc=log.a.r.pke[,grepl(',4]',colnames(log.a.r.pke))]

plot(m.log.a.t.pke~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
lines(apply(log.a.r.bs,2,median)~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),col='goldenrod',lwd=2)
lines(apply(log.a.r.goa,2,median)~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),col='darkgreen',lwd=2)
lines(apply(log.a.r.nc,2,median)~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),col='navy',lwd=2)
lines(apply(log.a.r.sc,2,median)~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),col='orange',lwd=2)


tv.pke$save_object('./outputs/model fits/tv_fit_pinkeven.RDS')

pke.info=subset(stk.info,species=='Pink-Even')
pke.dat=stk.dat[stk.dat$stock.id%in%pke.info$stock.id,]

N.s.pke=rag_n(pke.dat$stock.id)

smax.pke=pke.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat.pke=make_design_matrix(pke.dat$spawners,grp=pke.dat$stock.id)

stk.year=expand.grid(levels(factor(pke.dat$stock.id)),seq(min(pke.dat$broodyear),max(pke.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
pke.dat$stk.year=match(paste(pke.dat$stock.id,pke.dat$broodyear,sep='.'),stk.year[,3])

dl.rw=list(N=nrow(pke.dat),
           J=length(unique(pke.dat$stock.id)),
           L=max(pke.dat$broodyear)-min(pke.dat$broodyear)+1,
           R=length(unique(pke.info$ocean.basin)),
           J_i=as.numeric(factor(pke.dat$stock.id)),
           J_ii=pke.dat$stk.year,
           J_or=as.numeric(factor(pke.info$ocean.basin)),
           R_S=pke.dat$logRS,
           S=S.mat.pke,
           pSmax_mean=0.5*smax.pke$max.S,
           pSmax_sig=2*smax.pke$max.S)

tv.pke_rw=     m.tv2$sample(data=dl.rw,
                            chains = 6, 
                            iter_warmup = 500,
                            iter_sampling =500,
                            refresh = 250,
                            adapt_delta = 0.999,
                            max_treedepth = 20)

tv.pke_rw$save_object('./outputs/model fits/ns_fit_pink_even.RDS')


log.a.t.pke=tv.pke$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.pke=apply(exp(log.a.t.pke),2,median)
l90.log.a.t.pke=apply(exp(log.a.t.pke),2,quantile,0.025)
u90.log.a.t.pke=apply(exp(log.a.t.pke),2,quantile,0.975)
log.a.jt.pke=tv.pke$draws(variable='log_a_t',format='draws_matrix')

plot(m.log.a.t.pke~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),ylim=c(0,20),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.pke[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.pke))],2,median)
  lines(exp(logaj)~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=pke.info$begin[j]-min(pke.dat$broodyear)+1
  end=pke.info$end[j]-min(pke.dat$broodyear)+1
  lines(exp(logaj[st:end])~seq(pke.info$begin[j],pke.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.pke~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.pke~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.pke~seq(min(pke.dat$broodyear),max(pke.dat$broodyear)),lty=5,col='darkred')

ll1=apply(tv.pke_k2$draws(variables='log_lik',format='draws_matrix'),1,sum)
ll2=apply(tv.pke$draws(variables='log_lik',format='draws_matrix'),1,sum)

#MARSS model test
file.tv=file.path(cmdstanr::cmdstan_path(),'sr models', "hier_tvalpha_ricker.stan")
m.tv2=cmdstanr::cmdstan_model(file.tv)

file.tv=file.path(cmdstanr::cmdstan_path(),'sr models', "hier_tvalpha_ricker_simple.stan")
m.tv22=cmdstanr::cmdstan_model(file.tv)

file.s1=file.path(cmdstanr::cmdstan_path(),'sr models', "static_hier_mvn.stan")
m.s1=cmdstanr::cmdstan_model(file.s1)

pke.info=subset(stk.info,species=='Pink-Even')
pke.dat=stk.dat[stk.dat$stock.id%in%pke.info$stock.id,]

smax.pke=pke.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S_mat=time_design_matrix_pinks(s=pke.dat$spawners,t=pke.dat$broodyear,j=pke.dat$stock.id)
R_S_mat=time_design_matrix_pinks(s=pke.dat$logRS,t=pke.dat$broodyear,j=pke.dat$stock.id)
miss_mat=miss_matrix(R_S_mat)

sr.pke=dplyr::distinct(pke.info,region,.keep_all=T)
sr.pke=sr.pke[order(sr.pke$region),] #dataframe of ordered subregions

dl2=list(N=nrow(pke.dat),
         J=length(unique(pke.dat$stock.id)),
         L=length(seq(min(pke.dat$broodyear),max(pke.dat$broodyear),by=2)),
         SR=length(unique(pke.info$region)),
         R=length(unique(pke.info$ocean.basin)),
         R_S=R_S_mat,
         S=S_mat,
         miss_matrix=miss_mat,
         sr_or=as.numeric(factor(sr.pke$ocean.basin)),
         j_sr=as.numeric(factor(pke.info$region)),
         pSmax_mean=0.5*smax.pke$max.S,
         pSmax_sig=2*smax.pke$max.S)

tv.pke_k2=     m.tv2$sample(data=dl2,
                             chains = 3,
                             iter_warmup = 200,
                             iter_sampling =200,
                             refresh = 10,
                             adapt_delta = 0.99,
                             max_treedepth = 30)
smax.pko=pko.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat.pko=make_design_matrix(pko.dat$spawners,grp=pko.dat$stock.id)



#Pink-Odd####
pko.info=subset(stk.info,species=='Pink-Odd')
pko.dat=stk.dat[stk.dat$stock.id%in%pko.info$stock.id,]

smax.pko=pko.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

bspline_year=splines::bs(pko.dat$broodyear,knots=seq(min(pko.dat$broodyear),max(pko.dat$broodyear),by=2),degree=3,intercept=F)
predspline_year=splines::bs(seq(min(pko.dat$broodyear),max(pko.dat$broodyear)),knots=seq(min(pko.dat$broodyear),max(pko.dat$broodyear),by=2),degree=3,intercept=F)

sr.pko=dplyr::distinct(pko.info,region,.keep_all=T)
sr.pko=sr.pko[order(sr.pko$region),] #dataframe of ordered subregions

dl2=list(N=nrow(pko.dat),
         J=length(unique(pko.dat$stock.id)),
         L=max(pko.dat$broodyear)-min(pko.dat$broodyear)+1,
         J_i=as.numeric(factor(pko.dat$stock.id)),
         sr_or=as.numeric(factor(sr.pko$ocean.basin)),
         j_or=as.numeric(factor(pko.info$ocean.basin)),
         j_sr=as.numeric(factor(pko.info$region)),
         R=length(unique(pko.info$ocean.basin)),
         SR=length(unique(pko.info$region)),
         R_S=pko.dat$logRS,
         S=pko.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.pko$max.S,
         pSmax_sig=2*smax.pko$max.S)

tv.pko_k2=     m.tv33$sample(data=dl2,
                             chains = 5, 
                             iter_warmup = 500,
                             iter_sampling =500,
                             refresh = 250,
                             adapt_delta = 0.995,
                             max_treedepth = 20)

log.a.t.pko=tv.pko_k2$draws(variable='log_a_t',format='draws_matrix')
m.log.a.t.pko=apply(exp(log.a.t.pko),2,median)
l90.log.a.t.pko=apply(exp(log.a.t.pko),2,quantile,0.025)
u90.log.a.t.pko=apply(exp(log.a.t.pko),2,quantile,0.975)
log.a.jt.pko=tv.pko_k2$draws(variable='log_a_jt',format='draws_matrix')

plot(m.log.a.t.pko~seq(min(pko.dat$broodyear),max(pko.dat$broodyear)),ylim=c(0,20),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.pko[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.pko))],2,median)
  lines(exp(logaj)~seq(min(pko.dat$broodyear),max(pko.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=pko.info$begin[j]-min(pko.dat$broodyear)+1
  end=pko.info$end[j]-min(pko.dat$broodyear)+1
  lines(exp(logaj[st:end])~seq(pko.info$begin[j],pko.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.pko~seq(min(pko.dat$broodyear),max(pko.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.pko~seq(min(pko.dat$broodyear),max(pko.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.pko~seq(min(pko.dat$broodyear),max(pko.dat$broodyear)),lty=5,col='darkred')

log.a.r.pko=tv.pko_k2$draws(variable='log_a_rt',format='draws_matrix')
log.a.r.bs=log.a.r.pko[,grepl(',1]',colnames(log.a.r.pko))]
log.a.r.goa=log.a.r.pko[,grepl(',2]',colnames(log.a.r.pko))]
log.a.r.nc=log.a.r.pko[,grepl(',3]',colnames(log.a.r.pko))]
log.a.r.sc=log.a.r.pko[,grepl(',4]',colnames(log.a.r.pko))]

plot(m.log.a.t.pko~seq(min(pko.dat$broodyear),max(pko.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
lines(apply(log.a.r.bs,2,median)~seq(min(pko.dat$broodyear),max(pko.dat$broodyear)),col='goldenrod',lwd=2)
lines(apply(log.a.r.goa,2,median)~seq(min(pko.dat$broodyear),max(pko.dat$broodyear)),col='darkgreen',lwd=2)
lines(apply(log.a.r.nc,2,median)~seq(min(pko.dat$broodyear),max(pko.dat$broodyear)),col='navy',lwd=2)
lines(apply(log.a.r.sc,2,median)~seq(min(pko.dat$broodyear),max(pko.dat$broodyear)),col='orange',lwd=2)

smax.pko=pko.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat.pko=make_design_matrix(pko.dat$spawners,grp=pko.dat$stock.id)

stk.year=expand.grid(levels(factor(pko.dat$stock.id)),seq(min(pko.dat$broodyear),max(pko.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
pko.dat$stk.year=match(paste(pko.dat$stock.id,pko.dat$broodyear,sep='.'),stk.year[,3])

dl.rw=list(N=nrow(pko.dat),
           J=length(unique(pko.dat$stock.id)),
           L=max(pko.dat$broodyear)-min(pko.dat$broodyear)+1,
           R=length(unique(pko.info$ocean.basin)),
           J_i=as.numeric(factor(pko.dat$stock.id)),
           J_ii=pko.dat$stk.year,
           J_or=as.numeric(factor(pko.info$ocean.basin)),
           R_S=pko.dat$logRS,
           S=S.mat.pko,
           pSmax_mean=0.5*smax.pko$max.S,
           pSmax_sig=2*smax.pko$max.S)

tv.pko_rw=     m.tv2$sample(data=dl.rw,
                            chains = 6, 
                            iter_warmup = 500,
                            iter_sampling =500,
                            refresh = 250,
                            adapt_delta = 0.999,
                            max_treedepth = 20)

tv.pko_rw$save_object('./outputs/model fits/ns_fit_pink_odd.RDS')

#Sockeye####
scy.info=subset(stk.info,species=='Sockeye')
scy.dat=stk.dat[stk.dat$stock.id%in%scy.info$stock.id,]

smax.scy=scy.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

bspline_year=splines::bs(scy.dat$broodyear,knots=seq(min(scy.dat$broodyear),max(scy.dat$broodyear),by=2),degree=3,intercept=F)
predspline_year=splines::bs(seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),knots=seq(min(scy.dat$broodyear),max(scy.dat$broodyear),by=2),degree=3,intercept=F)

sr.scy=dplyr::distinct(scy.info,region,.keep_all=T)
sr.scy=sr.scy[order(sr.scy$region),] #dataframe of ordered subregions

dl2=list(N=nrow(scy.dat),
         J=length(unique(scy.dat$stock.id)),
         L=max(scy.dat$broodyear)-min(scy.dat$broodyear)+1,
         J_i=as.numeric(factor(scy.dat$stock.id)),
         sr_or=as.numeric(factor(sr.scy$ocean.basin)),
         j_or=as.numeric(factor(scy.info$ocean.basin)),
         j_sr=as.numeric(factor(scy.info$region)),
         R=length(unique(scy.info$ocean.basin)),
         SR=length(unique(scy.info$region)),
         R_S=scy.dat$logRS,
         S=scy.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.scy$max.S,
         pSmax_sig=2*smax.scy$max.S)

tv.scy_k2=     m.tv33$sample(data=dl2,
                             chains = 5, 
                             iter_warmup = 500,
                             iter_sampling =500,
                             refresh = 250,
                             adapt_delta = 0.995,
                             max_treedepth = 20)


log.a.t.scy=tv.scy_k2$draws(variable='log_a_t',format='draws_matrix')
m.log.a.t.scy=apply(log.a.t.scy,2,median)
l90.log.a.t.scy=apply(log.a.t.scy,2,quantile,0.025)
u90.log.a.t.scy=apply(log.a.t.scy,2,quantile,0.975)
log.a.jt.scy=tv.scy_k2$draws(variable='log_a_jt',format='draws_matrix')

plot(m.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),ylim=c(-1,4),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l',main='Sockeye')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.scy[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.scy))],2,median)
  lines(logaj~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=scy.info$begin[j]-min(scy.dat$broodyear)+1
  end=scy.info$end[j]-min(scy.dat$broodyear)+1
  lines(logaj[st:end]~seq(scy.info$begin[j],scy.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
for(l in 1:dl.rw$L){
  lines(c(m.log.a.t.scy[c(l,l+1)])~seq((min(scy.dat$broodyear)+(l-1)),(min(scy.dat$broodyear)+l)),col=adjustcolor('darkred',alpha.f=test$SampleCount[l]/dl.rw$J),lwd=4)
  lines(c(l90.log.a.t.scy[c(l,l+1)])~seq((min(scy.dat$broodyear)+(l-1)),(min(scy.dat$broodyear)+l)),col=adjustcolor('darkred',alpha.f=test$SampleCount[l]/dl.rw$J),lwd=3,lty=4)
  lines(c(u90.log.a.t.scy[c(l,l+1)])~seq((min(scy.dat$broodyear)+(l-1)),(min(scy.dat$broodyear)+l)),col=adjustcolor('darkred',alpha.f=test$SampleCount[l]/dl.rw$J),lwd=3,lty=4)
}
lines(m.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),lty=5,col='darkred')

log.a.r.scy=tv.scy_k2$draws(variable='log_a_rt',format='draws_matrix')
log.a.r.bs=log.a.r.scy[,grepl(',1]',colnames(log.a.r.scy))]
log.a.r.goa=log.a.r.scy[,grepl(',2]',colnames(log.a.r.scy))]
log.a.r.nc=log.a.r.scy[,grepl(',3]',colnames(log.a.r.scy))]
log.a.r.sc=log.a.r.scy[,grepl(',4]',colnames(log.a.r.scy))]

plot(m.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l',main='Sockeye')
lines(apply(log.a.r.bs,2,median)~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),col='goldenrod',lwd=2)
lines(apply(log.a.r.goa,2,median)~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),col='darkgreen',lwd=2)
lines(apply(log.a.r.nc,2,median)~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),col='navy',lwd=2)
lines(apply(log.a.r.sc,2,median)~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),col='orange',lwd=2)

smax.scy=scy.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat.scy=make_design_matrix(scy.dat$spawners,grp=scy.dat$stock.id)

stk.year=expand.grid(levels(factor(scy.dat$stock.id)),seq(min(scy.dat$broodyear),max(scy.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
scy.dat$stk.year=match(paste(scy.dat$stock.id,scy.dat$broodyear,sep='.'),stk.year[,3])

dl.rw=list(N=nrow(scy.dat),
           J=length(unique(scy.dat$stock.id)),
           L=max(scy.dat$broodyear)-min(scy.dat$broodyear)+1,
           R=length(unique(scy.info$ocean.basin)),
           J_i=as.numeric(factor(scy.dat$stock.id)),
           J_ii=scy.dat$stk.year,
           J_or=as.numeric(factor(scy.info$ocean.basin)),
           R_S=scy.dat$logRS,
           S=S.mat.scy,
           pSmax_mean=0.5*smax.scy$max.S,
           pSmax_sig=2*smax.scy$max.S)

tv.scy_rw=     m.tv2$sample(data=dl.rw,
                            chains = 6, 
                            iter_warmup = 500,
                            iter_sampling =500,
                            refresh = 250,
                            adapt_delta = 0.999,
                            max_treedepth = 20)

tv.scy_rw$save_object('./outputs/model fits/ns_fit_sockeye.RDS')


log.a.t.scy=tv.scy_rw$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.scy=apply(log.a.t.scy,2,median)
l90.log.a.t.scy=apply(log.a.t.scy,2,quantile,0.025)
u90.log.a.t.scy=apply(log.a.t.scy,2,quantile,0.975)
log.a.jt.scy=tv.scy_rw$draws(variable='log_a_t',format='draws_matrix')

plot(m.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),ylim=c(-1,4),lwd=3,col='darkred',ylab='Productivity log(alpha)',xlab='Brood Cohort Year',type='n',bty='l',main='Sockeye')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.scy[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.scy))],2,median)
  lines(logaj~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=scy.info$begin[j]-min(scy.dat$broodyear)+1
  end=scy.info$end[j]-min(scy.dat$broodyear)+1
  lines(logaj[st:end]~seq(scy.info$begin[j],scy.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),lty=5,col='darkred')


log.a.r.scy=tv.scy_rw$draws(variable='mu_log_a_r',format='draws_matrix')
log.a.r.bs=log.a.r.scy[,grepl(',1]',colnames(log.a.r.scy))]
log.a.r.goa=log.a.r.scy[,grepl(',2]',colnames(log.a.r.scy))]
log.a.r.nc=log.a.r.scy[,grepl(',3]',colnames(log.a.r.scy))]
log.a.r.sc=log.a.r.scy[,grepl(',4]',colnames(log.a.r.scy))]

plot(m.log.a.t.scy~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),ylim=c(-1,3),lwd=3,col='darkred',ylab='Productivity log(alpha)',xlab='Brood Cohort Year',type='n',bty='l',main='Sockeye')
lines(apply(log.a.r.bs,2,median)~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),col='goldenrod',lwd=2)
lines(apply(log.a.r.goa,2,median)~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),col='darkgreen',lwd=2)
lines(apply(log.a.r.nc,2,median)~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),col='navy',lwd=2)
lines(apply(log.a.r.sc,2,median)~seq(min(scy.dat$broodyear),max(scy.dat$broodyear)),col='orange',lwd=2)

loglik.gam.scy=apply(tv.scy_k2$draws(variable='log_lik',format='draws_matrix'),1,sum)
loglik.rw.scy=apply(tv.scy_rw$draws(variable='log_lik',format='draws_matrix'),1,sum)
summary(loglik.gam.scy)
summary(loglik.rw.scy)

smax1=apply(tv.scy_k2$draws(variables='Smax',format='draws_matrix'),2,median)
smax2=apply(tv.scy_rw$draws(variables='Smax',format='draws_matrix'),2,median)
plot(smax1~smax2)
abline(c(0,1))
summary(smax1/smax2)

smsy1=apply(tv.scy_k2$draws(variables='Smsy_jt',format='draws_matrix'),2,median)
smsy2=apply(tv.scy_rw$draws(variables='Smsy',format='draws_matrix'),2,median)
plot(smsy1~smsy2)
abline(c(0,1))
summary(smsy1/smsy2)

logat1=apply(tv.scy_k2$draws(variables='log_a_jt',format='draws_matrix'),2,median)
logat2=apply(tv.scy_rw$draws(variables='log_a_t',format='draws_matrix'),2,median)
plot(logat1~logat2)
abline(c(0,1))
summary(logat1/logat2)

tv.scy$save_object('./outputs/model fits/tv_fit_sockeye.RDS')



get_samples_by_year <- function(data, year_col, series_id, value) {
  # Group by year and count non-NA observations per year
  library(dplyr)
  samples_by_year <- data %>%
    dplyr::group_by(!!sym(year_col)) %>%
    dplyr::summarize(SampleCount = sum(!is.na(.data[[value]])), .groups = "drop") %>%
    dplyr::rename(Year = !!sym(year_col))
  return(samples_by_year)
}
#Gam hierarchical time-varying productivity fits

#Chinook####
chk.info=subset(stk.info,species=='Chinook')
chk.dat=stk.dat[stk.dat$stock.id%in%chk.info$stock.id,]

N.s.chk=rag_n(chk.dat$stock.id)

smax.chk=chk.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

bspline_year=splines::bs(chk.dat$broodyear,knots=seq(min(chk.dat$broodyear),max(chk.dat$broodyear),by=3),degree=3,intercept=F)
predspline_year=splines::bs(seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),knots=seq(min(chk.dat$broodyear),max(chk.dat$broodyear),by=3),degree=3,intercept=F)

S.mat.chk=make_design_matrix(chk.dat$spawners,grp=chk.dat$stock.id)

dl2=list(N=nrow(chk.dat),
         J=length(unique(chk.dat$stock.id)),
         L=max(chk.dat$broodyear)-min(chk.dat$broodyear)+1,
         J_i=as.numeric(factor(chk.dat$stock.id)),
         R_S=chk.dat$logRS,
         S=chk.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.chk$max.S,
         pSmax_sig=2*smax.chk$max.S,
         start_y=N.s.chk[,1],
         end_y=N.s.chk[,2],
         R=length(unique(chk.info$ocean.basin)),
         R_j=as.numeric(factor(chk.info$ocean.basin)))

bspline_year=splines::bs(chk.dat$broodyear,knots=seq(min(chk.dat$broodyear),max(chk.dat$broodyear),by=1),degree=3,intercept=F)
predspline_year=splines::bs(seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),knots=seq(min(chk.dat$broodyear),max(chk.dat$broodyear),by=1),degree=3,intercept=F)

dl3=list(N=nrow(chk.dat),
         J=length(unique(chk.dat$stock.id)),
         L=max(chk.dat$broodyear)-min(chk.dat$broodyear)+1,
         J_i=as.numeric(factor(chk.dat$stock.id)),
         R_S=chk.dat$logRS,
         S=chk.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.chk$max.S,
         pSmax_sig=2*smax.chk$max.S,
         start_y=N.s.chk[,1],
         end_y=N.s.chk[,2],
         R=length(unique(chk.info$ocean.basin)),
         R_j=as.numeric(factor(chk.info$ocean.basin)))

tv.chk=     m.tv33$sample(data=dl2,
                         chains = 5, 
                         iter_warmup = 200,
                         iter_sampling =300,
                         refresh = 200,
                         adapt_delta = 0.99,
                         max_treedepth = 30)

tv.chk2=     m.tv33$sample(data=dl3,

                          chains = 5, 
                          iter_warmup = 200,
                          iter_sampling =300,
                          refresh = 200,
                          adapt_delta = 0.99,
                          max_treedepth = 20)

tv.chk3=readRDS('./outputs/model fits/tv_fit_chinook.RDS')


ll1=apply(tv.chk$draws(variable='log_lik',format='draws_matrix'),1,sum)
ll2=apply(tv.chk2$draws(variable='log_lik',format='draws_matrix'),1,sum)
ll3=apply(tv.chk3$draws(variable='log_lik',format='draws_matrix'),1,sum)
summary(ll1)
summary(ll2)
summary(ll3)

llr1=tv.chk$draws(variable='log_lik',format='draws_matrix')
llr3=tv.chk3$draws(variable='log_lik',format='draws_matrix')

loo::loo_model_weights(list(llr1,llr3),method='stacking')

#hier gam - k=1, df = 1
log.a.t.chk=tv.chk$draws(variable='log_a_t',format='draws_matrix')
m.log.a.t.chk=apply(exp(log.a.t.chk),2,median)
l90.log.a.t.chk=apply(exp(log.a.t.chk),2,quantile,0.025)
u90.log.a.t.chk=apply(exp(log.a.t.chk),2,quantile,0.975)
log.a.jt.chk=tv.chk$draws(variable='log_a_jt',format='draws_matrix')

plot(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),ylim=c(0,20),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.chk[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chk))],2,median)
  lines(exp(logaj)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=chk.info$begin[j]-min(chk.dat$broodyear)+1
  end=chk.info$end[j]-min(chk.dat$broodyear)+1
  lines(exp(logaj[st:end])~seq(chk.info$begin[j],chk.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,col='darkred')

#individual RW

log.a.t.chk=tv.chk3$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.chk=apply(exp(log.a.t.chk),2,median)
l90.log.a.t.chk=apply(exp(log.a.t.chk),2,quantile,0.025)
u90.log.a.t.chk=apply(exp(log.a.t.chk),2,quantile,0.975)
log.a.jt.chk=tv.chk3$draws(variable='log_a_t',format='draws_matrix')

plot(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),ylim=c(0,20),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.chk[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chk))],2,median)
  lines(exp(logaj)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=chk.info$begin[j]-min(chk.dat$broodyear)+1
  end=chk.info$end[j]-min(chk.dat$broodyear)+1
  lines(exp(logaj[st:end])~seq(chk.info$begin[j],chk.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,col='darkred')

#hier gam - k=3, df = 3
log.a.t.chk=tv.chk2$draws(variable='log_a_t',format='draws_matrix')
m.log.a.t.chk=apply(exp(log.a.t.chk),2,median)
l90.log.a.t.chk=apply(exp(log.a.t.chk),2,quantile,0.025)
u90.log.a.t.chk=apply(exp(log.a.t.chk),2,quantile,0.975)
log.a.jt.chk=tv.chk2$draws(variable='log_a_jt',format='draws_matrix')

plot(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),ylim=c(0,20),lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n',bty='l')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.chk[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chk))],2,median)
  lines(exp(logaj)~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=chk.info$begin[j]-min(chk.dat$broodyear)+1
  end=chk.info$end[j]-min(chk.dat$broodyear)+1
  lines(exp(logaj[st:end])~seq(chk.info$begin[j],chk.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lwd=3,col='darkred')
lines(l90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.chk~seq(min(chk.dat$broodyear),max(chk.dat$broodyear)),lty=5,col='darkred')


#Chum
chm.info=subset(stk.info,species=='Chum')
chm.dat=stk.dat[stk.dat$stock.id%in%chm.info$stock.id,]

N.s.chm=rag_n(chm.dat$stock.id)

smax.chm=chm.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

bspline_year=splines::bs(chm.dat$broodyear,knots=seq(min(chm.dat$broodyear),max(chm.dat$broodyear),by=1),degree=2,intercept=F)
predspline_year=splines::bs(seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),knots=seq(min(chm.dat$broodyear),max(chm.dat$broodyear),by=1),degree=2,intercept=F)

S.mat.chm=make_design_matrix(chm.dat$spawners,grp=chm.dat$stock.id)

dl2=list(N=nrow(chm.dat),
         J=length(unique(chm.dat$stock.id)),
         L=max(chm.dat$broodyear)-min(chm.dat$broodyear)+1,
         J_i=as.numeric(factor(chm.dat$stock.id)),
         R_S=chm.dat$logRS,
         S=chm.dat$spawners,
         Nb=ncol(bspline_year),
         B_year=bspline_year,
         B_pred_year=predspline_year,
         pSmax_mean=0.5*smax.chm$max.S,
         pSmax_sig=2*smax.chm$max.S,
         start_y=N.s.chm[,1],
         end_y=N.s.chm[,2],
         R=length(unique(chm.info$ocean.basin)),
         R_j=as.numeric(factor(chm.info$ocean.basin)))

tv.chm=     m.tv3$sample(data=dl2,
                         init=0,
                         chains = 6, 
                         iter_warmup = 200,
                         iter_sampling =300,
                         refresh = 200,
                         adapt_delta = 0.99,
                         max_treedepth = 20)

tv.chm3=     m.tv33$sample(data=dl2,
                         chains = 6, 
                         iter_warmup = 200,
                         iter_sampling =300,
                         refresh = 200,
                         adapt_delta = 0.99,
                         max_treedepth = 20)

ll.k4.df2=apply(tv.chm3$draws(variable='log_lik',format='draws_matrix'),1,sum)

log.a.t.chm=tv.chm3$draws(variable='log_a_t',format='draws_matrix')
m.log.a.t.chm=apply(exp(log.a.t.chm),2,median)
l90.log.a.t.chm=apply(exp(log.a.t.chm),2,quantile,0.025)
u90.log.a.t.chm=apply(exp(log.a.t.chm),2,quantile,0.975)
log.a.jt.chm=tv.chm3$draws(variable='log_a_jt',format='draws_matrix')


plot(m.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),ylim=c(min(l90.log.a.t.chm),max(u90.log.a.t.chm)),type='l',lwd=3,col='darkred',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year')
lines(l90.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=5,col='darkred')
lines(u90.log.a.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=5,col='darkred')
for(j in 1:dl2$J){
  logaj=apply(log.a.jt.chm[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chm))],2,median)
  lines(exp(logaj)~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=chm.info$begin[j]-min(chm.dat$broodyear)+1
  end=chm.info$end[j]-min(chm.dat$broodyear)+1
  lines(exp(logaj[st:end])~seq(chm.info$begin[j],chm.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}

umsy.t.chm=tv.chm3$draws(variable='Umsy',format='draws_matrix')
m.umsy.t.chm=apply(umsy.t.chm,2,median)
l90.umsy.t.chm=apply(umsy.t.chm,2,quantile,0.025)
u90.umsy.t.chm=apply(umsy.t.chm,2,quantile,0.975)
umsy.jt.chm=tv.chm3$draws(variable='Umsy_jt',format='draws_matrix')


plot(m.umsy.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),ylim=c(0,1),type='l',ylab='Productivity (alpha) - Recruits/Spawner',xlab='Brood Cohort Year',type='n')
for(j in 1:dl2$J){
  umsyj=apply(umsy.jt.chm[,grepl(paste(',',j,']',sep=''),colnames(umsy.jt.chm))],2,median)
  lines(umsyj~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=3,col=adjustcolor('darkgray',alpha.f=0.4))
  st=chm.info$begin[j]-min(chm.dat$broodyear)+1
  end=chm.info$end[j]-min(chm.dat$broodyear)+1
  lines(umsyj[st:end]~seq(chm.info$begin[j],chm.info$end[j]),col=adjustcolor('darkgray',alpha.f=0.8))
}
lines(m.umsy.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lwd=3,col='darkred')
lines(l90.umsy.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=5,col='darkred')
lines(u90.umsy.t.chm~seq(min(chm.dat$broodyear),max(chm.dat$broodyear)),lty=5,col='darkred')


##test single stock
file.gam4=file.path(cmdstanr::cmdstan_path(),'sr models', "ss_m3_gam.stan")
m.tv4=cmdstanr::cmdstan_model(file.gam4)

file.m3=file.path(cmdstanr::cmdstan_path(),'sr models', "m3f_ip.stan")
m.m3=cmdstanr::cmdstan_model(file.m3)

dat=subset(chk.dat,stock.id==260)

bspline_year=splines::bs(dat$broodyear,knots=seq(min(dat$broodyear),max(dat$broodyear),by=4),degree=2,intercept=F)
predspline_year=splines::bs(seq(min(dat$broodyear),max(dat$broodyear)),knots=seq(min(dat$broodyear),max(dat$broodyear),by=4),degree=2,intercept=F)


dl=list(N=nrow(dat),
        L=max(dat$broodyear)-min(dat$broodyear)+1,
        R_S=dat$logRS,
        S=dat$spawners,
        ii=dat$broodyear-min(dat$broodyear)+1,
        Nb=ncol(bspline_year),
        B_year=bspline_year,
        B_pred_year=predspline_year,
        pSmax_mean=0.5*max(dat$spawners),
        pSmax_sig=2*max(dat$spawners))

test=    m.tv4$sample(data=dl,
                            chains = 3, 
                            iter_warmup = 200,
                            iter_sampling =600,
                            refresh = 200,
                            adapt_delta = 0.999,
                            max_treedepth = 20)

test2=    m.m3$sample(data=dl,
                      chains = 3, 
                      iter_warmup = 200,
                      iter_sampling =600,
                      refresh = 200,
                      adapt_delta = 0.999,
                      max_treedepth = 20)

pars_gam=test$draws(variables=c('log_a_t','Smax','sigma'),format='draws_matrix')
pars_rw=test2$draws(variables=c('log_a','Smax','sigma'),format='draws_matrix')

log_a_gam_m=apply(pars_gam[,1:dl$L],2,median)
log_a_gam_l=apply(pars_gam[,1:dl$L],2,quantile,0.05)
log_a_gam_u=apply(pars_gam[,1:dl$L],2,quantile,0.95)

log_a_rw_m=apply(pars_rw[,1:dl$L],2,median)
log_a_rw_l=apply(pars_rw[,1:dl$L],2,quantile,0.05)
log_a_rw_u=apply(pars_rw[,1:dl$L],2,quantile,0.95)

ll_gam=test$draws(variables=c('log_lik'),format='draws_matrix')
ll_rw=test2$draws(variables=c('log_lik'),format='draws_matrix')

plot(log_a_gam_m~dat$broodyear,ylim=c(min(c(log_a_rw_l,log_a_gam_l)),max(max(c(log_a_rw_u,log_a_gam_u)))),type='l',lwd=3,col='darkred',main='Nicola Chinook',ylab='log(alpha)',xlab='brood year')
lines(log_a_rw_m~dat$broodyear,lwd=3,col='navy')
lines(log_a_rw_l~dat$broodyear,col='navy',lty=5)
lines(log_a_rw_u~dat$broodyear,col='navy',lty=5)
lines(log_a_gam_l~dat$broodyear,col='darkred',lty=5)
lines(log_a_gam_u~dat$broodyear,col='darkred',lty=5)
text(x=par('usr')[2]-5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.025,paste('GAM, df =',attr(bspline_year,'degree'),',k =',attr(bspline_year,'knots')[2]-attr(bspline_year,'knots')[1]),col='darkred')
text(x=par('usr')[2]-5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.05,'RW',col='navy')

text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.025,paste('Smax =',round(median(pars_gam[,dl$L+1]))),col='darkred')
text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.05,paste('Smax =',round(median(pars_rw[,dl$L+1]))),col='navy')

text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.075,paste('sigma =',round(median(pars_gam[,dl$L+2]),digits=3)),col='darkred')
text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.1,paste('sigma =',round(median(pars_rw[,dl$L+2]),digits=3)),col='navy')

text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.125,paste('LL =',round(median(apply(ll_gam,1,sum)),digits=1)),col='darkred')
text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.15,paste('LL =',round(median(apply(ll_rw,1,sum)),digits=1)),col='navy')



dat=subset(chk.dat,stock.id==261)

bspline_year=splines::bs(dat$broodyear,knots=seq(min(dat$broodyear),max(dat$broodyear),by=4),degree=2,intercept=F)
predspline_year=splines::bs(seq(min(dat$broodyear),max(dat$broodyear)),knots=seq(min(dat$broodyear),max(dat$broodyear),by=4),degree=2,intercept=F)


dl=list(N=nrow(dat),
        L=max(dat$broodyear)-min(dat$broodyear)+1,
        R_S=dat$logRS,
        S=dat$spawners,
        ii=dat$broodyear-min(dat$broodyear)+1,
        Nb=ncol(bspline_year),
        B_year=bspline_year,
        B_pred_year=predspline_year,
        pSmax_mean=0.5*max(dat$spawners),
        pSmax_sig=2*max(dat$spawners))

test=    m.tv4$sample(data=dl,
                      chains = 3, 
                      iter_warmup = 200,
                      iter_sampling =600,
                      refresh = 200,
                      adapt_delta = 0.999,
                      max_treedepth = 20)

test2=    m.m3$sample(data=dl,
                      chains = 3, 
                      iter_warmup = 200,
                      iter_sampling =600,
                      refresh = 200,
                      adapt_delta = 0.999,
                      max_treedepth = 20)

pars_gam=test$draws(variables=c('log_a_t','Smax','sigma'),format='draws_matrix')
pars_rw=test2$draws(variables=c('log_a','Smax','sigma'),format='draws_matrix')

log_a_gam_m=apply(pars_gam[,1:dl$L],2,median)
log_a_gam_l=apply(pars_gam[,1:dl$L],2,quantile,0.05)
log_a_gam_u=apply(pars_gam[,1:dl$L],2,quantile,0.95)

log_a_rw_m=apply(pars_rw[,1:dl$L],2,median)
log_a_rw_l=apply(pars_rw[,1:dl$L],2,quantile,0.05)
log_a_rw_u=apply(pars_rw[,1:dl$L],2,quantile,0.95)

ll_gam=test$draws(variables=c('log_lik'),format='draws_matrix')
ll_rw=test2$draws(variables=c('log_lik'),format='draws_matrix')

plot(log_a_gam_m~dat$broodyear,ylim=c(min(c(log_a_rw_l,log_a_gam_l)),max(max(c(log_a_rw_u,log_a_gam_u)))),type='l',lwd=3,col='darkred',main='Middle Skeena Chinook',ylab='log(alpha)',xlab='brood year')
lines(log_a_rw_m~dat$broodyear,lwd=3,col='navy')
lines(log_a_rw_l~dat$broodyear,col='navy',lty=5)
lines(log_a_rw_u~dat$broodyear,col='navy',lty=5)
lines(log_a_gam_l~dat$broodyear,col='darkred',lty=5)
lines(log_a_gam_u~dat$broodyear,col='darkred',lty=5)
text(x=par('usr')[2]-5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.025,paste('GAM, df =',attr(bspline_year,'degree'),',k =',attr(bspline_year,'knots')[2]-attr(bspline_year,'knots')[1]),col='darkred')
text(x=par('usr')[2]-5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.05,'RW',col='navy')

text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.025,paste('Smax =',round(median(pars_gam[,dl$L+1]))),col='darkred')
text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.05,paste('Smax =',round(median(pars_rw[,dl$L+1]))),col='navy')

text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.075,paste('sigma =',round(median(pars_gam[,dl$L+2]),digits=3)),col='darkred')
text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.1,paste('sigma =',round(median(pars_rw[,dl$L+2]),digits=3)),col='navy')

text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.125,paste('LL =',round(median(apply(ll_gam,1,sum)),digits=1)),col='darkred')
text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.15,paste('LL =',round(median(apply(ll_rw,1,sum)),digits=1)),col='navy')

dat=subset(chk.dat,stock.id==257)

bspline_year=splines::bs(dat$broodyear,knots=seq(min(dat$broodyear),max(dat$broodyear),by=3),degree=3,intercept=F)
predspline_year=splines::bs(seq(min(dat$broodyear),max(dat$broodyear)),knots=seq(min(dat$broodyear),max(dat$broodyear),by=3),degree=3,intercept=F)


dl=list(N=nrow(dat),
        L=max(dat$broodyear)-min(dat$broodyear)+1,
        R_S=dat$logRS,
        S=dat$spawners,
        ii=dat$broodyear-min(dat$broodyear)+1,
        Nb=ncol(bspline_year),
        B_year=bspline_year,
        B_pred_year=predspline_year,
        pSmax_mean=0.5*max(dat$spawners),
        pSmax_sig=2*max(dat$spawners))

test=    m.tv4$sample(data=dl,
                      chains = 3, 
                      iter_warmup = 200,
                      iter_sampling =600,
                      refresh = 200,
                      adapt_delta = 0.999,
                      max_treedepth = 20)

test2=    m.m3$sample(data=dl,
                      chains = 3, 
                      iter_warmup = 200,
                      iter_sampling =600,
                      refresh = 200,
                      adapt_delta = 0.999,
                      max_treedepth = 20)

pars_gam=test$draws(variables=c('log_a_t','Smax','sigma'),format='draws_matrix')
pars_rw=test2$draws(variables=c('log_a','Smax','sigma'),format='draws_matrix')

log_a_gam_m=apply(pars_gam[,1:dl$L],2,median)
log_a_gam_l=apply(pars_gam[,1:dl$L],2,quantile,0.05)
log_a_gam_u=apply(pars_gam[,1:dl$L],2,quantile,0.95)

log_a_rw_m=apply(pars_rw[,1:dl$L],2,median)
log_a_rw_l=apply(pars_rw[,1:dl$L],2,quantile,0.05)
log_a_rw_u=apply(pars_rw[,1:dl$L],2,quantile,0.95)

ll_gam=test$draws(variables=c('log_lik'),format='draws_matrix')
ll_rw=test2$draws(variables=c('log_lik'),format='draws_matrix')

plot(log_a_gam_m~dat$broodyear,ylim=c(min(c(log_a_rw_l,log_a_gam_l)),max(max(c(log_a_rw_u,log_a_gam_u)))),type='l',lwd=3,col='darkred',main='Cowichan Chinook',ylab='log(alpha)',xlab='brood year')
lines(log_a_rw_m~dat$broodyear,lwd=3,col='navy')
lines(log_a_rw_l~dat$broodyear,col='navy',lty=5)
lines(log_a_rw_u~dat$broodyear,col='navy',lty=5)
lines(log_a_gam_l~dat$broodyear,col='darkred',lty=5)
lines(log_a_gam_u~dat$broodyear,col='darkred',lty=5)
text(x=par('usr')[2]-5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.025,paste('GAM, df =',attr(bspline_year,'degree'),',k =',attr(bspline_year,'knots')[2]-attr(bspline_year,'knots')[1]),col='darkred')
text(x=par('usr')[2]-5,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.05,'RW',col='navy')

text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.025,paste('Smax =',round(median(pars_gam[,dl$L+1]))),col='darkred')
text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.05,paste('Smax =',round(median(pars_rw[,dl$L+1]))),col='navy')

text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.075,paste('sigma =',round(median(pars_gam[,dl$L+2]),digits=3)),col='darkred')
text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.1,paste('sigma =',round(median(pars_rw[,dl$L+2]),digits=3)),col='navy')

text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.125,paste('LL =',round(median(apply(ll_gam,1,sum)),digits=1)),col='darkred')
text(x=par('usr')[2]-2,y=par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.15,paste('LL =',round(median(apply(ll_rw,1,sum)),digits=1)),col='navy')
