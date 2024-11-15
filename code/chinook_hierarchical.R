#chinook hierarchical model
library(here);library(dplyr)

#autocorrelated Ricker model
cmdstanr::set_cmdstan_path(path='C:/Users/greenbergda/Documents/.cmdstan/cmdstan-2.29.2')
library(cmdstanr)
options(mc.cores=6)

file=file.path(cmdstanr::cmdstan_path(),'sr models', "hier_static_ricker.stan")
m=cmdstanr::cmdstan_model(file)

#data load
stock_dat<- read.csv(here('data','salmon_productivity_compilation2024-09-24.csv'))
ck.info=read.csv(here('data','chinook_overview_cf.csv'))

ck.dat=subset(stock_dat,stock.id%in%ck.info$stock.id)

sr=dplyr::distinct(ck.info,region,.keep_all=T)
sr=sr[order(sr$region),] #dataframe of ordered subregions

rag_n=function(x){
  rle=rle(x) #get run lengths of strings
  start_n=c(1,1+cumsum(rle$lengths)) #starting points
  start_n=start_n[-length(start_n)] #drop last value
  end_n=c(start_n[-1]-1) #end points
  end_n[length(start_n)]=length(x) #final obs.
  
  return(cbind(start_n,end_n))
}

#ragged start and end points for each SR series
N_s=rag_n(ck.dat$stock.id)

or_n

smaxs=ck.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

dl=list(N=nrow(ck.dat),
        J=length(unique(ck.dat$stock.id)),
        SR=length(unique(ck.info$region)),
        R=length(unique(ck.info$ocean.basin)),
        L=max(ck.dat$broodyear)-min(ck.dat$broodyear)+1,
        ii=ck.dat$broodyear-min(ck.dat$broodyear)+1,
        ij=as.numeric(factor(ck.dat$stock.id)),
        sr_or=as.numeric(factor(sr$ocean.basin)),
        j_sr=as.numeric(factor(ck.info$region)),
        start_y=N_s[,1],
        end_y=N_s[,2],
        R_S=ck.dat$logRS,
        S=ck.dat$spawners,
        pSmax_mean=0.5*smaxs$max.S,
        pSmax_sig=2*smaxs$max.S)

ck1=     m$sample(data=dl,
                         chains = 6, 
                         init=0,
                         iter_warmup = 500,
                         iter_sampling =1000,
                         refresh = 100,
                         adapt_delta = 0.999,
                         max_treedepth = 20)

#Fraser....
levels(factor(ck.info$region))
#=8

log_a_sr=ck1$draws(variables='log_a_sr',format='df')

log_a_fr=log_a_sr$`log_a_sr[8]`
d=density(log_a_fr,bw=0.02)
summary(log_a_fr)

plot(d$y~d$x,type='l',xlab='log(alpha)',ylab='density',bty='l')

abline(v=quantile(log_a_fr,seq(0,1,by=0.05)))



#post hoc subyearling vs yearling
syl=rownames(ck.info[ck.info$life.history=='subyearling',])
yl=rownames(ck.info[ck.info$life.history=='yearling',])

log_a=ck1$draws(variables='log_a',format='df')

log_a_syl=as.data.frame(log_a[,as.numeric(syl)])

log_a_yl=as.data.frame(log_a[,as.numeric(yl)])

summary(unlist(log_a_syl))
summary(unlist(log_a_yl))

dsyl=density(unlist(log_a_syl),bw=0.02)
dyl=density(unlist(log_a_yl),bw=0.02)

par(mfrow=c(1,2))
plot(dsyl$y~dsyl$x,type='l',xlab='log(alpha)',ylab='density',bty='l',main='subyearling stocks')

plot(dyl$y~dyl$x,type='l',xlab='log(alpha)',ylab='density',bty='l',main='yearling stocks')

sldat=data.frame(stock=ck.info$stock.name[ck.info$life.history=='subyearling'],prod=apply(log_a_syl,2,median))
rownames(sldat)=NULL
yldat=data.frame(stock=ck.info$stock.name[ck.info$life.history=='yearling'],prod=apply(log_a_yl,2,median))
rownames(yldat)=NULL

sylprod=unlist(log_a_syl)
ylprod=unlist(log_a_yl)

diff=sylprod-sample(ylprod,size=length(sylprod))
for(i in 1:19){
  d=sylprod-sample(ylprod,size=length(sylprod))
  diff=c(diff,d)
}


#correlation between smax & log(alpha) - within stock
b=ck1$draws(variables='b',format='df')

cor_ab=NA
for(j in 1:nrow(ck.info)){
  cor_ab[j]=cor(b[,j],log_a[,j])
}

hist(cor_ab,main='Correlation between log(alpha) and beta',xlab='correlation',breaks=30)



###Chum
stock_dat<- read.csv(here('data','salmon_productivity_compilation2024-09-04.csv'))
s.info=read.csv(here('data','stock_info2024-09-04.csv'))
cm.info=subset(s.info,species=='Chum')

cm.dat=subset(stock_dat,stock.id%in%ck.info$stock.id)

sr=dplyr::distinct(ck.info,region,.keep_all=T)
sr=sr[order(sr$region),] #dataframe of ordered subregions

rag_n=function(x){
  rle=rle(x) #get run lengths of strings
  start_n=c(1,1+cumsum(rle$lengths)) #starting points
  start_n=start_n[-length(start_n)] #drop last value
  end_n=c(start_n[-1]-1) #end points
  end_n[length(start_n)]=length(x) #final obs.
  
  return(cbind(start_n,end_n))
}

#ragged start and end points for each SR series
N_s=rag_n(ck.dat$stock.id)

or_n

smaxs=ck.dat%>%group_by(stock.id)%>%summarize(max.S=max(spawners))

dl=list(N=nrow(ck.dat),
        J=length(unique(ck.dat$stock.id)),
        SR=length(unique(ck.info$region)),
        R=length(unique(ck.info$ocean.basin)),
        L=max(ck.dat$broodyear)-min(ck.dat$broodyear)+1,
        ii=ck.dat$broodyear-min(ck.dat$broodyear)+1,
        ij=as.numeric(factor(ck.dat$stock.id)),
        sr_or=as.numeric(factor(sr$ocean.basin)),
        j_sr=as.numeric(factor(ck.info$region)),
        start_y=N_s[,1],
        end_y=N_s[,2],
        R_S=ck.dat$logRS,
        S=ck.dat$spawners,
        pSmax_mean=0.5*smaxs$max.S,
        pSmax_sig=2*smaxs$max.S)

ck1=     m$sample(data=dl,
                  chains = 6, 
                  init=0,
                  iter_warmup = 500,
                  iter_sampling =1000,
                  refresh = 100,
                  adapt_delta = 0.999,
                  max_treedepth = 20)

