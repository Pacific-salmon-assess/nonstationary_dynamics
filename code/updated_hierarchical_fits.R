#chinook hierarchical model
library(here);library(dplyr)
source(here('code','util_func.R'))

#models####
cmdstanr::set_cmdstan_path(path='C:/Users/greenbergda/Documents/.cmdstan/cmdstan-2.29.2')
library(cmdstanr)
options(mc.cores=6)

file=file.path(cmdstanr::cmdstan_path(),'sr models', "hier_static_ricker.stan")
m.static=cmdstanr::cmdstan_model(file)

file.tv=file.path(cmdstanr::cmdstan_path(),'sr models', "hier_tvalpha_ricker.stan")
m.tv=cmdstanr::cmdstan_model(file.tv)

file.tv=file.path(cmdstanr::cmdstan_path(),'sr models', "ind_tvalpha_ricker.stan")
m.tv2=cmdstanr::cmdstan_model(file.tv)

#stock data ####
stk.dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation2024-10-17.csv'))
stk.info=read.csv(here('data','filtered datasets','stock_info2024-10-17.csv'))


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
                  init=0,
                  iter_warmup = 500,
                  iter_sampling =1000,
                  refresh = 100,
                  adapt_delta = 0.999,
                  max_treedepth = 20)

tv.wc$save_object('./outputs/model fits/multivariate_fit_WC.RDS')
tv.wc=readRDS('./outputs/model fits/multivariate_fit_WC.RDS')


S.mat=make_design_matrix(ck.dat$spawners,grp=ck.dat$stock.id)

stk.year=expand.grid(levels(factor(ck.dat$stock.id)),seq(min(ck.dat$broodyear),max(ck.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,3]),]
ck.dat$stk.year=match(paste(ck.dat$stock.id,ck.dat$broodyear,sep='.'),stk.year[,3])

dl2=list(N=nrow(ck.dat),
         J=length(unique(ck.dat$stock.id)),
         L=max(ck.dat$broodyear)-min(ck.dat$broodyear)+1,
         J_i=as.numeric(factor(ck.dat$stock.id)),
         J_ii=ck.dat$stk.year,
         start_y=N.s.ck[,1],
         end_y=N.s.ck[,2],
         R_S=ck.dat$logRS,
         S=S.mat,
         pSmax_mean=0.5*smax.ck$max.S,
         pSmax_sig=2*smax.ck$max.S)

tv.ck=     m.tv$sample(data=dl2,
                            chains = 6, 
                            init=0,
                            iter_warmup = 500,
                            iter_sampling =1000,
                            refresh = 100,
                            adapt_delta = 0.999,
                            max_treedepth = 20)

mu_log_a=tv.ck$draws('mu_log_a',format='draws_matrix')
log_a_j=tv.ck$draws('log_a_t',format='draws_matrix')
plot(c(-1,4)~c(min(ck.dat$broodyear),max(ck.dat$broodyear)),bty='l',type='n',ylab='log(alpha)',xlab='Brood year')
for(i in 1:dl2$J){
 log_a_t=log_a_j[,grepl(paste(',',i,']',sep=''),colnames(log_a_j))]  
 lines(apply(log_a_t,2,median)~seq(min(ck.dat$broodyear),max(ck.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
 }

lines(apply(mu_log_a,2,median)~seq(min(ck.dat$broodyear),max(ck.dat$broodyear)),lwd=3)
lines(apply(mu_log_a,2,quantile,0.025)~seq(min(ck.dat$broodyear),max(ck.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_log_a,2,quantile,0.975)~seq(min(ck.dat$broodyear),max(ck.dat$broodyear)),lwd=2,lty=3)

umsy_j=tv.ck$draws('Umsy',format='draws_matrix')
plot(c(0,1)~c(min(ck.dat$broodyear),max(ck.dat$broodyear)),bty='l',type='n',ylab='Umsy',xlab='Brood year')
for(i in 1:dl2$J){
  umsy_t=umsy_j[,grepl(paste(',',i,']',sep=''),colnames(umsy_j))]  
  lines(apply(umsy_t,2,median)~seq(min(ck.dat$broodyear),max(ck.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
}
colnames(umsy_j)=paste('Umsy.',expand.grid(seq(1:dl2$L),seq(1:dl2$J))[,1],',',expand.grid(seq(1:dl2$L),seq(1:dl2$J))[,2],'.',sep='')


mu_umsy=matrix(ncol=dl2$L,nrow=nrow(umsy_j))
for(l in 1:dl2$L){
  mu_umsy[,l]=apply(umsy_j[,grepl(paste('Umsy.',l,',',sep=''),colnames(umsy_j))],1,mean)
}
lines(apply(mu_umsy,2,median)~seq(min(ck.dat$broodyear),max(ck.dat$broodyear)),lwd=3)
lines(apply(mu_umsy,2,quantile,0.025)~seq(min(ck.dat$broodyear),max(ck.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_umsy,2,quantile,0.975)~seq(min(ck.dat$broodyear),max(ck.dat$broodyear)),lwd=2,lty=3)

mu_umsy2=apply(mu_log_a,1:2,umsyCalc)

umsy_j=tv.ck$draws('Umsy',format='draws_matrix')
plot(c(0,1)~c(min(ck.dat$broodyear),max(ck.dat$broodyear)),bty='l',type='n',ylab='Umsy',xlab='Brood year')
for(i in 1:dl2$J){
  umsy_t=umsy_j[,grepl(paste(',',i,']',sep=''),colnames(umsy_j))]  
  lines(apply(umsy_t,2,median)~seq(min(ck.dat$broodyear),max(ck.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
}
colnames(umsy_j)=paste('Umsy.',expand.grid(seq(1:dl2$L),seq(1:dl2$J))[,1],',',expand.grid(seq(1:dl2$L),seq(1:dl2$J))[,2],'.',sep='')

lines(apply(mu_umsy2,2,median)~seq(min(ck.dat$broodyear),max(ck.dat$broodyear)),lwd=3)
lines(apply(mu_umsy2,2,quantile,0.025)~seq(min(ck.dat$broodyear),max(ck.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_umsy2,2,quantile,0.975)~seq(min(ck.dat$broodyear),max(ck.dat$broodyear)),lwd=2,lty=3)






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
                            init=0,
                            iter_warmup = 500,
                            iter_sampling =1000,
                            refresh = 100,
                            adapt_delta = 0.999,
                            max_treedepth = 20)

##coho####
coh.info=subset(stk.info,species=='Coho')
coh.dat=stk.dat[stk.dat$stock.id%in%coh.info$stock.id,]

N.s.coh=rag_n(coh.dat$stock.id)

smax.coh=coh.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

sr.coh=dplyr::distinct(coh.info,region,.keep_all=T)
sr.coh=sr.coh[order(sr.coh$region),] #dataframe of ordered subregions


dl1=list(N=nrow(coh.dat),
         J=length(unique(coh.dat$stock.id)),
         SR=length(unique(coh.info$region)),
         R=length(unique(coh.info$ocean.basin)),
         L=max(coh.dat$broodyear)-min(coh.dat$broodyear)+1,
         ii=coh.dat$broodyear-min(coh.dat$broodyear)+1,
         ij=as.numeric(factor(coh.dat$stock.id)),
         sr_or=as.numeric(factor(sr.coh$ocean.basin)),
         j_sr=as.numeric(factor(coh.info$region)),
         start_y=N.s.coh[,1],
         end_y=N.s.coh[,2],
         R_S=coh.dat$logRS,
         S=coh.dat$spawners,
         pSmax_mean=0.5*smax.coh$max.S,
         pSmax_sig=2*smax.coh$max.S)

stc.coh=     m.static$sample(data=dl1,
                            chains = 6, 
                            init=0,
                            iter_warmup = 500,
                            iter_sampling =1000,
                            refresh = 100,
                            adapt_delta = 0.999,
                            max_treedepth = 20)


##Pinks####
pk.info=subset(stk.info,species=='Pink-Even'|species=='Pink-Odd')
pk.dat=stk.dat[stk.dat$stock.id%in%pk.info$stock.id,]

#pink model
#sub-divide even-odd-two models?

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
                             init=0,
                             iter_warmup = 500,
                             iter_sampling =1000,
                             refresh = 100,
                             adapt_delta = 0.999,
                             max_treedepth = 20)


#Time-varying model set####

## West coast####
###Chinook####
wc.chk.info=subset(stk.info,ocean.basin=='WC'&species=='Chinook')
wc.chk.dat=stk.dat[stk.dat$stock.id%in%wc.chk.info$stock.id,]

N.s.wc=rag_n(wc.chk.dat$stock.id)

smax.wc=wc.chk.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat=make_design_matrix(wc.chk.dat$spawners,grp=wc.chk.dat$stock.id)

stk.year=expand.grid(levels(factor(wc.chk.dat$stock.id)),seq(min(wc.chk.dat$broodyear),max(wc.chk.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
wc.chk.dat$stk.year=match(paste(wc.chk.dat$stock.id,wc.chk.dat$broodyear,sep='.'),stk.year[,3])

dl2=list(N=nrow(wc.chk.dat),
         J=length(unique(wc.chk.dat$stock.id)),
         L=max(wc.chk.dat$broodyear)-min(wc.chk.dat$broodyear)+1,
         J_i=as.numeric(factor(wc.chk.dat$stock.id)),
         J_ii=wc.chk.dat$stk.year,
         start_y=N.s.wc[,1],
         end_y=N.s.wc[,2],
         R_S=wc.chk.dat$logRS,
         S=S.mat,
         pSmax_mean=0.5*smax.wc$max.S,
         pSmax_sig=2*smax.wc$max.S)

tv.wc.chk=     m.tv$sample(data=dl2,
                           chains = 6, 
                           iter_warmup = 300,
                           iter_sampling =500,
                           refresh = 100,
                           adapt_delta = 0.95,
                           max_treedepth = 20)


corr_heatmap(tv.wc.chk,info=wc.chk.info,init=T,title='West Coast Chinook')
corr_heatmap(tv.wc.chk,info=wc.chk.info,init=F,k=5,title='West Coast Chinook')




### Chum####
wc.chm.info=subset(stk.info,ocean.basin=='WC'&species=='Chum')
wc.chm.dat=stk.dat[stk.dat$stock.id%in%wc.chm.info$stock.id,]

N.s.wc=rag_n(wc.chm.dat$stock.id)

smax.wc=wc.chm.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat=make_design_matrix(wc.chm.dat$spawners,grp=wc.chm.dat$stock.id)

stk.year=expand.grid(levels(factor(wc.chm.dat$stock.id)),seq(min(wc.chm.dat$broodyear),max(wc.chm.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
wc.chm.dat$stk.year=match(paste(wc.chm.dat$stock.id,wc.chm.dat$broodyear,sep='.'),stk.year[,3])

dl2=list(N=nrow(wc.chm.dat),
         J=length(unique(wc.chm.dat$stock.id)),
         L=max(wc.chm.dat$broodyear)-min(wc.chm.dat$broodyear)+1,
         J_i=as.numeric(factor(wc.chm.dat$stock.id)),
         J_ii=wc.chm.dat$stk.year,
         start_y=N.s.wc[,1],
         end_y=N.s.wc[,2],
         R_S=wc.chm.dat$logRS,
         S=S.mat,
         pSmax_mean=0.5*smax.wc$max.S,
         pSmax_sig=2*smax.wc$max.S)

tv.wc.chm=     m.tv$sample(data=dl2,
                            chains = 6, 
                            init=0,
                            iter_warmup = 300,
                            iter_sampling =500,
                            refresh = 100,
                            adapt_delta = 0.95,
                            max_treedepth = 20)

corr_heatmap(tv.wc.chm,info=wc.chm.info,init=T,title='West Coast Chum')
corr_heatmap(tv.wc.chm,info=wc.chm.info,init=F,k=4,title='West Coast Chinook')

tv.wc.chm2=     m.tv2$sample(data=dl2,
                           chains = 6, 
                           init=0,
                           iter_warmup = 300,
                           iter_sampling =500,
                           refresh = 100,
                           adapt_delta = 0.95,
                           max_treedepth = 20)




###Coho####
wc.coh.info=subset(stk.info,ocean.basin=='WC'&species=='Coho')
wc.coh.dat=stk.dat[stk.dat$stock.id%in%wc.coh.info$stock.id,]

N.s.wc=rag_n(wc.coh.dat$stock.id)

smax.wc=wc.coh.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat=make_design_matrix(wc.coh.dat$spawners,grp=wc.coh.dat$stock.id)

stk.year=expand.grid(levels(factor(wc.coh.dat$stock.id)),seq(min(wc.coh.dat$broodyear),max(wc.coh.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
wc.coh.dat$stk.year=match(paste(wc.coh.dat$stock.id,wc.coh.dat$broodyear,sep='.'),stk.year[,3])

dl2=list(N=nrow(wc.coh.dat),
         J=length(unique(wc.coh.dat$stock.id)),
         L=max(wc.coh.dat$broodyear)-min(wc.coh.dat$broodyear)+1,
         J_i=as.numeric(factor(wc.coh.dat$stock.id)),
         J_ii=wc.coh.dat$stk.year,
         start_y=N.s.wc[,1],
         end_y=N.s.wc[,2],
         R_S=wc.coh.dat$logRS,
         S=S.mat,
         pSmax_mean=0.5*smax.wc$max.S,
         pSmax_sig=2*smax.wc$max.S)

tv.wc.coh=     m.tv$sample(data=dl2,
                           chains = 4, 
                           iter_warmup = 300,
                           iter_sampling =500,
                           refresh = 100,
                           adapt_delta = 0.99,
                           max_treedepth = 20)

corr_heatmap(tv.wc.coh,info=wc.coh.info,init=T,title='West Coast Pink (odd)')


###Pink####

####even####
wc.pke.info=subset(stk.info,ocean.basin=='WC'&species=='Pink-Even')
wc.pke.dat=stk.dat[stk.dat$stock.id%in%wc.pke.info$stock.id,]

N.s.wc=rag_n(wc.pke.dat$stock.id)

smax.wc=wc.pke.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat=make_design_matrix(wc.pke.dat$spawners,grp=wc.pke.dat$stock.id)

stk.year=expand.grid(levels(factor(wc.pke.dat$stock.id)),seq(min(wc.pke.dat$broodyear),max(wc.pke.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
wc.pke.dat$stk.year=match(paste(wc.pke.dat$stock.id,wc.pke.dat$broodyear,sep='.'),stk.year[,3])

dl2=list(N=nrow(wc.pke.dat),
         J=length(unique(wc.pke.dat$stock.id)),
         L=max(wc.pke.dat$broodyear)-min(wc.pke.dat$broodyear)+1,
         J_i=as.numeric(factor(wc.pke.dat$stock.id)),
         J_ii=wc.pke.dat$stk.year,
         start_y=N.s.wc[,1],
         end_y=N.s.wc[,2],
         R_S=wc.pke.dat$logRS,
         S=S.mat,
         pSmax_mean=0.5*smax.wc$max.S,
         pSmax_sig=2*smax.wc$max.S)

tv.wc.pke=     m.tv$sample(data=dl2,
                           chains = 6, 
                           iter_warmup = 300,
                           iter_sampling =500,
                           refresh = 100,
                           adapt_delta = 0.95,
                           max_treedepth = 20)


corr_heatmap(tv.wc.pke,info=wc.pke.info,init=T,title='West Coast Pink')




####odd####
wc.pko.info=subset(stk.info,ocean.basin=='WC'&species=='Pink-Odd')
wc.pko.dat=stk.dat[stk.dat$stock.id%in%wc.pko.info$stock.id,]

N.s.wc=rag_n(wc.pko.dat$stock.id)

smax.wc=wc.pko.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat=make_design_matrix(wc.pko.dat$spawners,grp=wc.pko.dat$stock.id)

stk.year=expand.grid(levels(factor(wc.pko.dat$stock.id)),seq(min(wc.pko.dat$broodyear),max(wc.pko.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
wc.pko.dat$stk.year=match(paste(wc.pko.dat$stock.id,wc.pko.dat$broodyear,sep='.'),stk.year[,3])

dl2=list(N=nrow(wc.pko.dat),
         J=length(unique(wc.pko.dat$stock.id)),
         L=max(wc.pko.dat$broodyear)-min(wc.pko.dat$broodyear)+1,
         J_i=as.numeric(factor(wc.pko.dat$stock.id)),
         J_ii=wc.pko.dat$stk.year,
         start_y=N.s.wc[,1],
         end_y=N.s.wc[,2],
         R_S=wc.pko.dat$logRS,
         S=S.mat,
         pSmax_mean=0.5*smax.wc$max.S,
         pSmax_sig=2*smax.wc$max.S)

tv.wc.pko=     m.tv$sample(data=dl2,
                           chains = 4, 
                           iter_warmup = 300,
                           iter_sampling =500,
                           refresh = 100,
                           adapt_delta = 0.99,
                           max_treedepth = 20)

corr_heatmap(tv.wc.pko,info=wc.pko.info,init=T,title='West Coast Pink (odd)')

pnk.r=which(wc.info$species=='Pink-Odd')

loga_wc_pnk=log_a_j[,c((pnk.r[1]-1)*dl2$L+1):c(dl2$L*pnk.r[1])]

for(r in 2:length(pnk.r)){
  loga_wc_pnk=cbind(loga_wc_pnk,log_a_j[,c((pnk.r[r]-1)*dl2$L+1):c(dl2$L*pnk.r[r])])
}

plot(c(-1,4)~c(min(wc.dat$broodyear),max(wc.dat$broodyear)),bty='l',type='n',ylab='log(alpha)',xlab='Brood year',main='West Coast (Pink - Odd)')
for(i in 1:length(pnk.r)){
  log_a_t=log_a_j[,grepl(paste(',',pnk.r[i],']',sep=''),colnames(log_a_j))]  
  lines(apply(log_a_t,2,median)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
}

colnames(loga_wc_pnk)=paste('log_a.',expand.grid(seq(1,dl2$L),seq(1,length(pnk.r)))[,1],',',expand.grid(seq(1,dl2$L),seq(1,length(pnk.r)))[,2],'.',sep='')
mu_log_a_wc_pnk=matrix(nrow=nrow(log_a_t),ncol=dl2$L)
for(l in 1:dl2$L){
  log_a_wc_t=loga_wc_pnk[,grepl(paste('log_a.',l,',',sep=''),colnames(loga_wc_pnk))]  
  mu_log_a_wc_pnk[,l]=apply(log_a_wc_t,1,mean)
}

lines(apply(mu_log_a_wc_pnk,2,median)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),lwd=3)
lines(apply(mu_log_a_wc_pnk,2,quantile,0.025)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_log_a_wc_pnk,2,quantile,0.975)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),lwd=2,lty=3)

plot(c(0,1)~c(min(wc.dat$broodyear),max(wc.dat$broodyear)),bty='l',type='n',ylab='Umsy',xlab='Brood year',main='West Coast (Pink - Odd)')
#abline(h=seq(0,1,by=0.2),col=adjustcolor('gray',alpha.f=0.5),lty=5)
for(i in 1:length(pnk.r)){
  umsy_t=umsy_j[,grepl(paste(',',pnk.r[i],']',sep=''),colnames(umsy_j))]  
  lines(apply(umsy_t,2,median)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
}

mu_umsy_wc_pnk=apply(mu_log_a_wc_pnk,1:2,umsyCalc)

lines(apply(mu_umsy_wc_pnk,2,median)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),lwd=3)
lines(apply(mu_umsy_wc_pnk,2,quantile,0.025)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_umsy_wc_pnk,2,quantile,0.975)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),lwd=2,lty=3)


###Sockeye####
wc.soc.info=subset(stk.info,ocean.basin=='WC'&species=='Sockeye')
wc.soc.dat=stk.dat[stk.dat$stock.id%in%wc.soc.info$stock.id,]

N.s.wc=rag_n(wc.soc.dat$stock.id)

smax.wc=wc.soc.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat=make_design_matrix(wc.soc.dat$spawners,grp=wc.soc.dat$stock.id)

stk.year=expand.grid(levels(factor(wc.soc.dat$stock.id)),seq(min(wc.soc.dat$broodyear),max(wc.soc.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,1],stk.year[,2]),]
wc.soc.dat$stk.year=match(paste(wc.soc.dat$stock.id,wc.soc.dat$broodyear,sep='.'),stk.year[,3])

dl2=list(N=nrow(wc.soc.dat),
         J=length(unique(wc.soc.dat$stock.id)),
         L=max(wc.soc.dat$broodyear)-min(wc.soc.dat$broodyear)+1,
         J_i=as.numeric(factor(wc.soc.dat$stock.id)),
         J_ii=wc.soc.dat$stk.year,
         start_y=N.s.wc[,1],
         end_y=N.s.wc[,2],
         R_S=wc.soc.dat$logRS,
         S=S.mat,
         pSmax_mean=0.5*smax.wc$max.S,
         pSmax_sig=2*smax.wc$max.S)

tv.wc.soc=     m.tv$sample(data=dl2,
                           chains = 6,
                           init=0,
                           iter_warmup = 500,
                           iter_sampling =1000,
                           refresh = 100,
                           adapt_delta = 0.99,
                           max_treedepth = 20)

corr_heatmap(tv.wc.soc,info=wc.soc.info,init=T,title='West Coast Sockeye')


soc.r=which(wc.info$species=='Sockeye')

loga_wc_soc=log_a_j[,1:c(dl2$L*max(soc.r))]
plot(c(-1,4)~c(min(wc.dat$broodyear),max(wc.dat$broodyear)),bty='l',type='n',ylab='log(alpha)',xlab='Brood year',main='West Coast (Sockeye)')
for(i in 1:length(soc.r)){
  log_a_t=loga_wc_soc[,grepl(paste(',',soc.r[i],']',sep=''),colnames(loga_wc_soc))]  
  lines(apply(log_a_t,2,median)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
}

colnames(loga_wc_soc)=paste('log_a.',expand.grid(seq(1,dl2$L),seq(1,length(soc.r)))[,1],',',expand.grid(seq(1,dl2$L),seq(1,length(soc.r)))[,2],'.',sep='')
mu_log_a_wc_soc=matrix(nrow=nrow(log_a_t),ncol=dl2$L)
for(l in 1:dl2$L){
  log_a_wc_t=loga_wc_soc[,grepl(paste('log_a.',l,',',sep=''),colnames(loga_wc_soc))]  
  mu_log_a_wc_soc[,l]=apply(log_a_wc_t,1,mean)
}

lines(apply(mu_log_a_wc_soc,2,median)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),lwd=3)
lines(apply(mu_log_a_wc_soc,2,quantile,0.025)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_log_a_wc_soc,2,quantile,0.975)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),lwd=2,lty=3)

umsy_wc_soc=umsy_j[,soc.r[1]:c(dl2$L*max(soc.r))]

plot(c(0,1)~c(min(wc.dat$broodyear),max(wc.dat$broodyear)),bty='l',type='n',ylab='Umsy',xlab='Brood year',main='West Coast (Sockeye)')
#abline(h=seq(0,1,by=0.2),col=adjustcolor('gray',alpha.f=0.5),lty=5)
for(i in 1:length(soc.r)){
  umsy_t=umsy_wc_soc[,grepl(paste(',',i,']',sep=''),colnames(umsy_wc_soc))]  
  lines(apply(umsy_t,2,median)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
}

mu_umsy_wc_soc=apply(mu_log_a_wc_soc,1:2,umsyCalc)

lines(apply(mu_umsy_wc_soc,2,median)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),lwd=3)
lines(apply(mu_umsy_wc_soc,2,quantile,0.025)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_umsy_wc_soc,2,quantile,0.975)~seq(min(wc.dat$broodyear),max(wc.dat$broodyear)),lwd=2,lty=3)




#Gulf of Alaska####
goa.info=subset(stk.info,ocean.basin=='GOA')
goa.dat=stk.dat[stk.dat$stock.id%in%goa.info$stock.id,]

N.s.goa=rag_n(goa.dat$stock.id)

smax.goa=goa.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat=make_design_matrix(goa.dat$spawners,grp=goa.dat$stock.id)

stk.year=expand.grid(levels(factor(goa.dat$stock.id)),seq(min(goa.dat$broodyear),max(goa.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,3]),]
goa.dat$stk.year=match(paste(goa.dat$stock.id,goa.dat$broodyear,sep='.'),stk.year[,3])

dl2=list(N=nrow(goa.dat),
         J=length(unique(goa.dat$stock.id)),
         L=max(goa.dat$broodyear)-min(goa.dat$broodyear)+1,
         J_i=as.numeric(factor(goa.dat$stock.id)),
         J_ii=goa.dat$stk.year,
         start_y=N.s.goa[,1],
         end_y=N.s.goa[,2],
         R_S=goa.dat$logRS,
         S=S.mat,
         pSmax_mean=0.5*smax.goa$max.S,
         pSmax_sig=2*smax.goa$max.S)

tv.goa=     m.tv$sample(data=dl2,
                       chains = 6, 
                       init=0,
                       iter_warmup = 500,
                       iter_sampling =1000,
                       refresh = 100,
                       adapt_delta = 0.999,
                       max_treedepth = 20)

tv.goa$save_object('./outputs/model fits/multivariate_fit_GOA.RDS')
tv.goa=readRDS('./outputs/model fits/multivariate_fit_GOA.RDS')


mu_log_a=tv.goa$draws('mu_log_a',format='draws_matrix')
log_a_j=tv.goa$draws('log_a_t',format='draws_matrix')
plot(c(-1,4)~c(min(goa.dat$broodyear),max(goa.dat$broodyear)),bty='l',type='n',ylab='log(alpha)',xlab='Brood year',main='Gulf of Alaska (all species)')
for(i in 1:dl2$J){
  log_a_t=log_a_j[,grepl(paste(',',i,']',sep=''),colnames(log_a_j))]  
  lines(apply(log_a_t,2,median)~seq(min(goa.dat$broodyear),max(goa.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
}

lines(apply(mu_log_a,2,median)~seq(min(goa.dat$broodyear),max(goa.dat$broodyear)),lwd=3)
lines(apply(mu_log_a,2,quantile,0.025)~seq(min(goa.dat$broodyear),max(goa.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_log_a,2,quantile,0.975)~seq(min(goa.dat$broodyear),max(goa.dat$broodyear)),lwd=2,lty=3)


mu_umsy2=apply(mu_log_a,1:2,umsyCalc)

umsy_j=tv.goa$draws('Umsy',format='draws_matrix')
plot(c(0,1)~c(min(goa.dat$broodyear),max(goa.dat$broodyear)),bty='l',type='n',ylab='Umsy',xlab='Brood year',main='Gulf of Alaska (all species)')
for(i in 1:dl2$J){
  umsy_t=umsy_j[,grepl(paste(',',i,']',sep=''),colnames(umsy_j))]  
  lines(apply(umsy_t,2,median)~seq(min(goa.dat$broodyear),max(goa.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
}

lines(apply(mu_umsy2,2,median)~seq(min(goa.dat$broodyear),max(goa.dat$broodyear)),lwd=3)
lines(apply(mu_umsy2,2,quantile,0.025)~seq(min(goa.dat$broodyear),max(goa.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_umsy2,2,quantile,0.975)~seq(min(goa.dat$broodyear),max(goa.dat$broodyear)),lwd=2,lty=3)

#Bering Sea####
bs.info=subset(stk.info,ocean.basin=='BS')
bs.dat=stk.dat[stk.dat$stock.id%in%bs.info$stock.id,]

N.s.bs=rag_n(bs.dat$stock.id)

smax.bs=bs.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

S.mat=make_design_matrix(bs.dat$spawners,grp=bs.dat$stock.id)

stk.year=expand.grid(levels(factor(bs.dat$stock.id)),seq(min(bs.dat$broodyear),max(bs.dat$broodyear)))
stk.year[,3]=paste(stk.year[,1],stk.year[,2],sep='.')
stk.year=stk.year[order(stk.year[,3]),]
bs.dat$stk.year=match(paste(bs.dat$stock.id,bs.dat$broodyear,sep='.'),stk.year[,3])

dl2=list(N=nrow(bs.dat),
         J=length(unique(bs.dat$stock.id)),
         L=max(bs.dat$broodyear)-min(bs.dat$broodyear)+1,
         J_i=as.numeric(factor(bs.dat$stock.id)),
         J_ii=bs.dat$stk.year,
         start_y=N.s.bs[,1],
         end_y=N.s.bs[,2],
         R_S=bs.dat$logRS,
         S=S.mat,
         pSmax_mean=0.5*smax.bs$max.S,
         pSmax_sig=2*smax.bs$max.S)

tv.bs=     m.tv$sample(data=dl2,
                        chains = 6, 
                        init=0,
                        iter_warmup = 500,
                        iter_sampling =1000,
                        refresh = 100,
                        adapt_delta = 0.999,
                        max_treedepth = 20)

tv.bs$save_object('./outputs/model fits/multivariate_fit_bs.RDS')
tv.bs=readRDS('./outputs/model fits/multivariate_fit_bs.RDS')

mu_log_a=tv.bs$draws('mu_log_a',format='draws_matrix')
log_a_j=tv.bs$draws('log_a_t',format='draws_matrix')
plot(c(-1,4)~c(min(bs.dat$broodyear),max(bs.dat$broodyear)),bty='l',type='n',ylab='log(alpha)',xlab='Brood year',main='Bering Sea (all species)')
for(i in 1:dl2$J){
  log_a_t=log_a_j[,grepl(paste(',',i,']',sep=''),colnames(log_a_j))]  
  lines(apply(log_a_t,2,median)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
}

lines(apply(mu_log_a,2,median)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),lwd=3)
lines(apply(mu_log_a,2,quantile,0.025)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_log_a,2,quantile,0.975)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),lwd=2,lty=3)

mu_umsy2=apply(mu_log_a,1:2,umsyCalc)

umsy_j=tv.bs$draws('Umsy',format='draws_matrix')
plot(c(0,1)~c(min(bs.dat$broodyear),max(bs.dat$broodyear)),bty='l',type='n',ylab='Umsy',xlab='Brood year',main='Bering Sea (all species)')
for(i in 1:dl2$J){
  umsy_t=umsy_j[,grepl(paste(',',i,']',sep=''),colnames(umsy_j))]  
  lines(apply(umsy_t,2,median)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
}

lines(apply(mu_umsy2,2,median)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),lwd=3)
lines(apply(mu_umsy2,2,quantile,0.025)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_umsy2,2,quantile,0.975)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),lwd=2,lty=3)

#By species
#sockeye
soc.r=which(bs.info$species=='Sockeye')


loga_bs_soc=log_a_j[,1:c(dl2$L*max(soc.r))]
plot(c(-1,4)~c(min(bs.dat$broodyear),max(bs.dat$broodyear)),bty='l',type='n',ylab='log(alpha)',xlab='Brood year',main='Bering Sea (Bristol Bay Sockeye)')
for(i in 1:length(soc.r)){
  log_a_t=loga_bs_soc[,grepl(paste(',',i,']',sep=''),colnames(loga_bs_soc))]  
  lines(apply(log_a_t,2,median)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
}

colnames(loga_bs_soc)=paste('log_a.',expand.grid(seq(1,dl2$L),seq(1,length(soc.r)))[,1],',',expand.grid(seq(1,dl2$L),seq(1,length(soc.r)))[,2],'.',sep='')
mu_log_a_bs_soc=matrix(nrow=nrow(log_a_t),ncol=dl2$L)
for(l in 1:dl2$L){
  log_a_bs_t=loga_bs_soc[,grepl(paste('log_a.',l,',',sep=''),colnames(loga_bs_soc))]  
  mu_log_a_bs_soc[,l]=apply(log_a_bs_t,1,mean)
}

lines(apply(mu_log_a_bs_soc,2,median)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),lwd=3)
lines(apply(mu_log_a_bs_soc,2,quantile,0.025)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_log_a_bs_soc,2,quantile,0.975)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),lwd=2,lty=3)

umsy_bs_soc=umsy_j[,1:c(dl2$L*max(soc.r))]

plot(c(0,1)~c(min(bs.dat$broodyear),max(bs.dat$broodyear)),bty='l',type='n',ylab='Umsy',xlab='Brood year',main='Bering Sea (Bristol Bay Sockeye)')
for(i in 1:length(soc.r)){
  umsy_t=umsy_bs_soc[,grepl(paste(',',i,']',sep=''),colnames(umsy_bs_soc))]  
  lines(apply(umsy_t,2,median)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),col=adjustcolor('darkgray',alpha.f=0.5))
}

mu_umsy_bs_soc=apply(mu_log_a_bs_soc,1:2,umsyCalc)

lines(apply(mu_umsy_bs_soc,2,median)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),lwd=3)
lines(apply(mu_umsy_bs_soc,2,quantile,0.025)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),lwd=2,lty=3)
lines(apply(mu_umsy_bs_soc,2,quantile,0.975)~seq(min(bs.dat$broodyear),max(bs.dat$broodyear)),lwd=2,lty=3)
