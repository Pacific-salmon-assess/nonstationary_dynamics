library(here);library(dplyr)
source(here('code','util_func.R'))

#models####
cmdstanr::set_cmdstan_path(path='C:/Users/greenbergda/Documents/.cmdstan/cmdstan-2.29.2')
library(cmdstanr)
options(mc.cores=6)

#stationary model
file=file.path(cmdstanr::cmdstan_path(),'sr models', "hier_static_ricker.stan")
m.static=cmdstanr::cmdstan_model(file)

#non-stationary
file.tv=file.path(cmdstanr::cmdstan_path(),'sr models', "ind_tvalpha_ricker.stan")
m.tv=cmdstanr::cmdstan_model(file.tv)

#gam non-stationary
file.gam=file.path(cmdstanr::cmdstan_path(),'sr models', "gam_tvalpha_hier.stan")
m.gam=cmdstanr::cmdstan_model(file.gam)

#data ####
stk.dat<- read.csv(here('data','salmon_productivity_compilation2026-08-10.csv'))
stk.info=read.csv(here('data','stock_info2026-08-10.csv'))

#Stationary model set####
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
                            adapt_delta = 0.999,
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
                            adapt_delta = 0.999,
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
                             adapt_delta = 0.999,
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

#Non-stationary model sets####
##Chinook####
chk.info=subset(stk.info,species=='Chinook')
chk.dat=stk.dat[stk.dat$stock.id%in%chk.info$stock.id,]

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

tv.chk_rw=     m.tv$sample(data=dl.rw,
                           chains = 6, 
                           init=0,
                           iter_warmup = 500,
                           iter_sampling =500,
                           refresh = 250,
                           adapt_delta = 0.999,
                           max_treedepth = 20)

tv.chk_rw$save_object('./outputs/model fits/ns_fit_chinook.RDS')

##Chum####

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

tv.chm_rw=     m.tv$sample(data=dl.rw,
                           chains = 6, 
                           iter_warmup = 500,
                           iter_sampling =500,
                           refresh = 250,
                           adapt_delta = 0.999,
                           max_treedepth = 20)

tv.chm_rw$save_object('./outputs/model fits/ns_fit_chum.RDS')

##Coho####
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

tv.cho_rw=     m.tv$sample(data=dl.rw,
                           chains = 6, 
                           iter_warmup = 500,
                           iter_sampling =500,
                           refresh = 250,
                           adapt_delta = 0.999,
                           max_treedepth = 20)

tv.cho_rw$save_object('./outputs/model fits/ns_fit_coho.RDS')

##Sockeye####
scy.info=subset(stk.info,species=='Sockeye')
scy.dat=stk.dat[stk.dat$stock.id%in%scy.info$stock.id,]

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

tv.scy_rw=     m.tv$sample(data=dl.rw,
                           chains = 6, 
                           iter_warmup = 500,
                           iter_sampling =500,
                           refresh = 250,
                           adapt_delta = 0.999,
                           max_treedepth = 20)

tv.scy_rw$save_object('./outputs/model fits/ns_fit_sockeye.RDS')


#GAM alternative test with Chum####
chm.info=subset(stk.info,species=='Chum')
chm.dat=stk.dat[stk.dat$stock.id%in%chm.info$stock.id,]

smax.chm=chm.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

#cubic spline, knot biannually
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

tv.chm_k2=     m.gam$sample(data=dl2,
                            chains = 5, 
                            iter_warmup = 500,
                            iter_sampling =500,
                            refresh = 250,
                            adapt_delta = 0.995,
                            max_treedepth = 20)

tv.chm_k2$save_object('./outputs/model fits/gam_fit_chum_k2.RDS')

