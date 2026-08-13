library(here);library(cmdstanr)
#results in text metrics
stock.dat<- read.csv(here('data','salmon_productivity_compilation2026-08-10.csv'))
stock.info<- read.csv(here('data','stock_info2026-08-10.csv'))
stock.info$state=factor(stock.info$state,levels=c('WA','BC','AK'))
stock.info$species2=gsub("-.*", "", stock.info$species) # collapse even/odd pinks together

stk.info=stock.info
chk.info=subset(stock.info,species=='Chinook')
scy.info=subset(stock.info,species=='Sockeye')
chm.info=subset(stock.info,species=='Chum')
cho.info=subset(stock.info,species=='Coho')
pke.info=subset(stock.info,species=='Pink-Even')
pko.info=subset(stock.info,species=='Pink-Odd')

#nonstationary####
tv.chk=readRDS(here('outputs','model fits','ns_fit_chinook.RDS'))
tv.cho=readRDS(here('outputs','model fits','ns_fit_coho.RDS'))
tv.chm=readRDS(here('outputs','model fits','ns_fit_chum.RDS'))
tv.pke=readRDS(here('outputs','model fits','ns_fit_pink_even.RDS'))
tv.pko=readRDS(here('outputs','model fits','ns_fit_pink_odd.RDS'))
tv.scy=readRDS(here('outputs','model fits','ns_fit_sockeye.RDS'))

#early productivity vs. late productivity - coastwide

#chinook
log.a.t.chk=tv.chk$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.chk=apply(exp(log.a.t.chk),2,median)
names(m.log.a.t.chk)=paste(seq(1:length(m.log.a.t.chk)),seq(min(chk.info$begin),max(chk.info$end)))
#prod - 1954:1986
median(apply(exp(log.a.t.chk[,1:31]),1,mean))
HDInterval::hdi(apply(exp(log.a.t.chk[,1:31]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(log.a.t.chk[,57:67]),1,mean))
HDInterval::hdi(apply(exp(log.a.t.chk[,57:67]),1,mean),credMass=0.9)
#Percent decline
median(1-apply(exp(log.a.t.chk[,57:67]),1,mean)/apply(exp(log.a.t.chk[,1:31]),1,mean))
HDInterval::hdi(1-apply(exp(log.a.t.chk[,57:67]),1,mean)/apply(exp(log.a.t.chk[,1:31]),1,mean),credMass=0.9)


#Sockeye
log.a.t.scy=tv.scy$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.scy=apply(exp(log.a.t.scy),2,median)
names(m.log.a.t.scy)=paste(seq(1:length(m.log.a.t.scy)),seq(min(scy.info$begin),max(scy.info$end)))
#prod - 1948:1978
median(apply(exp(log.a.t.scy[,1:31]),1,mean))
HDInterval::hdi(apply(exp(log.a.t.scy[,1:31]),1,mean),credMass=0.9)
#prod - 2010:2019
median(apply(exp(log.a.t.scy[,63:73]),1,mean))
HDInterval::hdi(apply(exp(log.a.t.scy[,63:73]),1,mean),credMass=0.9)
#Percent decline
median(1-apply(exp(log.a.t.scy[,63:73]),1,mean)/apply(exp(log.a.t.scy[,1:31]),1,mean))
HDInterval::hdi(1-apply(exp(log.a.t.scy[,63:73]),1,mean)/apply(exp(log.a.t.scy[,1:31]),1,mean),credMass=0.9)

#Chum
log.a.t.chm=tv.chm$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.chm=apply(exp(log.a.t.chm),2,median)
names(m.log.a.t.chm)=paste(seq(1:length(m.log.a.t.chm)),seq(min(chm.info$begin),max(chm.info$end)))
#prod - 1948:1977
median(apply(exp(log.a.t.chm[,1:32]),1,mean))
HDInterval::hdi(apply(exp(log.a.t.chm[,1:32]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(log.a.t.chm[,58:68]),1,mean))
HDInterval::hdi(apply(exp(log.a.t.chm[,58:68]),1,mean),credMass=0.9)
#Percent decline
median(1-apply(exp(log.a.t.chm[,58:68]),1,mean)/apply(exp(log.a.t.chm[,1:31]),1,mean))
HDInterval::hdi(1-apply(exp(log.a.t.chm[,58:68]),1,mean)/apply(exp(log.a.t.chm[,1:31]),1,mean),credMass=0.9)

#Coho
log.a.t.cho=tv.cho$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.cho=apply(exp(log.a.t.cho),2,median)
names(m.log.a.t.cho)=paste(seq(1:length(m.log.a.t.cho)),seq(min(cho.info$begin),max(cho.info$end)))
#prod - 1948:1977
median(apply(exp(log.a.t.cho[,1:37]),1,mean))
HDInterval::hdi(apply(exp(log.a.t.cho[,1:37]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(log.a.t.cho[,38:46]),1,mean))
HDInterval::hdi(apply(exp(log.a.t.cho[,38:46]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(log.a.t.cho[,38:46]),1,mean)/apply(exp(log.a.t.cho[,1:37]),1,mean))-1
HDInterval::hdi(apply(exp(log.a.t.cho[,38:46]),1,mean)/apply(exp(log.a.t.cho[,1:37]),1,mean)-1,credMass=0.9)

#Pink-e
log.a.t.pke=tv.pke$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.pke=apply(exp(log.a.t.pke),2,median)
names(m.log.a.t.pke)=paste(seq(1:length(m.log.a.t.pke)),seq(min(pke.info$begin),max(pke.info$end)))
#prod - 1948:1977
median(apply(exp(log.a.t.pke[,1:31]),1,mean))
HDInterval::hdi(apply(exp(log.a.t.pke[,1:31]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(log.a.t.pke[,51:71]),1,mean))
HDInterval::hdi(apply(exp(log.a.t.pke[,61:71]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(log.a.t.pke[,51:71]),1,mean)/apply(exp(log.a.t.pke[,1:31]),1,mean))-1
HDInterval::hdi(apply(exp(log.a.t.pke[,51:71]),1,mean)/apply(exp(log.a.t.pke[,1:31]),1,mean)-1,credMass=0.9)

#early productivity vs. late productivity - by ocean region

#chinook
log.a.jt.chk=tv.chk$draws(variable='log_a_t',format='draws_matrix')
chk.sc=which(chk.info$ocean.basin=='SC')
chk.nc=which(chk.info$ocean.basin=='NC')
chk.goa=which(chk.info$ocean.basin=='GoA')
chk.bs=which(chk.info$ocean.basin=='BS')

logat=log.a.jt.chk
colnames(logat)=gsub("]", ".", colnames(logat))
colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)

logat.sc=logat[,grepl(paste(paste(',',chk.sc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.nc=logat[,grepl(paste(paste(',',chk.nc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.goa=logat[,grepl(paste(paste(',',chk.goa,'.',sep=''),collapse = "|"),colnames(logat))]
logat.bs=logat[,grepl(paste(paste(',',chk.bs,'.',sep=''),collapse = "|"),colnames(logat))]

L=max(chk.info$end)-min(chk.info$begin)+1
mu_loga.chk.sc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.chk.nc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.chk.goa=matrix(nrow=nrow(logat),ncol=L)
mu_loga.chk.bs=matrix(nrow=nrow(logat),ncol=L)
for(t in 1:L){
  logas.sc=logat.sc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.sc))]
  mu_loga.chk.sc[,t]=apply(logas.sc,1,mean)
  
  logas.nc=logat.nc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.nc))]
  mu_loga.chk.nc[,t]=apply(logas.nc,1,mean)
  
  logas.goa=logat.goa[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.goa))]
  mu_loga.chk.goa[,t]=apply(logas.goa,1,mean)
  
  logas.bs=logat.bs[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.bs))]
  mu_loga.chk.bs[,t]=apply(logas.bs,1,mean)
}

#south coast

m.log.chk.sc=apply(exp(mu_loga.chk.sc),2,median)
names(m.log.chk.sc)=paste(seq(1:length(m.log.chk.sc)),seq(min(chk.info$begin),max(chk.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.chk.sc[,1:31]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chk.sc[,1:31]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.chk.sc[,57:67]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chk.sc[,57:67]),1,mean),credMass=0.9)
#Percent decline
median(1-apply(exp(mu_loga.chk.sc[,57:67]),1,mean)/apply(exp(mu_loga.chk.sc[,1:31]),1,mean))
HDInterval::hdi(1-apply(exp(mu_loga.chk.sc[,57:67]),1,mean)/apply(exp(mu_loga.chk.sc[,1:31]),1,mean),credMass=0.9)

#north coast
m.log.chk.nc=apply(exp(mu_loga.chk.nc),2,median)
names(m.log.chk.nc)=paste(seq(1:length(m.log.chk.nc)),seq(min(chk.info$begin),max(chk.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.chk.nc[,1:31]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chk.nc[,1:31]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.chk.nc[,57:67]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chk.nc[,57:67]),1,mean),credMass=0.9)
#Percent decline
median(1-apply(exp(mu_loga.chk.nc[,57:67]),1,mean)/apply(exp(mu_loga.chk.nc[,1:31]),1,mean))
HDInterval::hdi(1-apply(exp(mu_loga.chk.nc[,57:67]),1,mean)/apply(exp(mu_loga.chk.nc[,1:31]),1,mean),credMass=0.9)

#gulf of alaska
m.log.chk.goa=apply(exp(mu_loga.chk.goa),2,median)
names(m.log.chk.goa)=paste(seq(1:length(m.log.chk.goa)),seq(min(chk.info$begin),max(chk.info$end)))
#prod - 1954:1991
median(apply(exp(mu_loga.chk.goa[,1:47]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chk.goa[,1:47]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.chk.goa[,57:67]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chk.goa[,57:67]),1,mean),credMass=0.9)
#Percent decline
median(1-apply(exp(mu_loga.chk.goa[,57:67]),1,mean)/apply(exp(mu_loga.chk.goa[,1:47]),1,mean))
HDInterval::hdi(1-apply(exp(mu_loga.chk.goa[,57:67]),1,mean)/apply(exp(mu_loga.chk.goa[,1:47]),1,mean),credMass=0.9)

#Bering sea
m.log.chk.bs=apply(exp(mu_loga.chk.bs),2,median)
names(m.log.chk.bs)=paste(seq(1:length(m.log.chk.bs)),seq(min(chk.info$begin),max(chk.info$end)))
#prod - 1954:1991
median(apply(exp(mu_loga.chk.bs[,1:47]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chk.bs[,1:47]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.chk.bs[,57:67]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chk.bs[,57:67]),1,mean),credMass=0.9)
#Percent decline
median(1-apply(exp(mu_loga.chk.bs[,57:67]),1,mean)/apply(exp(mu_loga.chk.bs[,1:47]),1,mean))
HDInterval::hdi(1-apply(exp(mu_loga.chk.bs[,57:67]),1,mean)/apply(exp(mu_loga.chk.bs[,1:47]),1,mean),credMass=0.9)

#sockeye
log.a.jt.scy=tv.scy$draws(variable='log_a_t',format='draws_matrix')
scy.sc=which(scy.info$ocean.basin=='SC')
scy.nc=which(scy.info$ocean.basin=='NC')
scy.goa=which(scy.info$ocean.basin=='GoA')
scy.bs=which(scy.info$ocean.basin=='BS')

logat=log.a.jt.scy
colnames(logat)=gsub("]", ".", colnames(logat))
colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)

logat.sc=logat[,grepl(paste(paste(',',scy.sc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.nc=logat[,grepl(paste(paste(',',scy.nc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.goa=logat[,grepl(paste(paste(',',scy.goa,'.',sep=''),collapse = "|"),colnames(logat))]
logat.bs=logat[,grepl(paste(paste(',',scy.bs,'.',sep=''),collapse = "|"),colnames(logat))]

L=max(scy.info$end)-min(scy.info$begin)+1
mu_loga.scy.sc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.scy.nc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.scy.goa=matrix(nrow=nrow(logat),ncol=L)
mu_loga.scy.bs=matrix(nrow=nrow(logat),ncol=L)
for(t in 1:L){
  logas.sc=logat.sc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.sc))]
  mu_loga.scy.sc[,t]=apply(logas.sc,1,mean)
  
  logas.nc=logat.nc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.nc))]
  mu_loga.scy.nc[,t]=apply(logas.nc,1,mean)
  
  logas.goa=logat.goa[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.goa))]
  mu_loga.scy.goa[,t]=apply(logas.goa,1,mean)
  
  logas.bs=logat.bs[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.bs))]
  mu_loga.scy.bs[,t]=apply(logas.bs,1,mean)
}

#south coast

m.log.scy.sc=apply(exp(mu_loga.scy.sc),2,median)
names(m.log.scy.sc)=paste(seq(1:length(m.log.scy.sc)),seq(min(scy.info$begin),max(scy.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.scy.sc[,1:31]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.scy.sc[,1:31]),1,mean),credMass=0.9)
#prod - 2010:2019
median(apply(exp(mu_loga.scy.sc[,63:72]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.scy.sc[,63:72]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.scy.sc[,63:72]),1,mean)/apply(exp(mu_loga.scy.sc[,1:31]),1,mean))-1
HDInterval::hdi(apply(exp(mu_loga.scy.sc[,63:72]),1,mean)/apply(exp(mu_loga.scy.sc[,1:31]),1,mean)-1,credMass=0.9)

#north coast
m.log.scy.nc=apply(exp(mu_loga.scy.nc),2,median)
names(m.log.scy.nc)=paste(seq(1:length(m.log.scy.nc)),seq(min(scy.info$begin),max(scy.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.scy.nc[,1:43]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.scy.nc[,1:43]),1,mean),credMass=0.9)
#prod - 2010:2019
median(apply(exp(mu_loga.scy.nc[,63:72]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.scy.nc[,63:72]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.scy.nc[,63:72]),1,mean)/apply(exp(mu_loga.scy.nc[,1:43]),1,mean)-1)
HDInterval::hdi(apply(exp(mu_loga.scy.nc[,63:72]),1,mean)/apply(exp(mu_loga.scy.nc[,1:43]),1,mean)-1,credMass=0.9)

#gulf of alaska
m.log.scy.goa=apply(exp(mu_loga.scy.goa),2,median)
names(m.log.scy.goa)=paste(seq(1:length(m.log.scy.goa)),seq(min(scy.info$begin),max(scy.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.scy.goa[,1:43]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.scy.goa[,1:43]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.scy.goa[,57:67]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.scy.goa[,57:67]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.scy.goa[,63:72]),1,mean)/apply(exp(mu_loga.scy.goa[,1:43]),1,mean))-1
HDInterval::hdi(apply(exp(mu_loga.scy.goa[,63:72]),1,mean)/apply(exp(mu_loga.scy.goa[,1:43]),1,mean)-1,credMass=0.9)

#Bering sea
m.log.scy.bs=apply(exp(mu_loga.scy.bs),2,median)
names(m.log.scy.bs)=paste(seq(1:length(m.log.scy.bs)),seq(min(scy.info$begin),max(scy.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.scy.bs[,1:23]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.scy.bs[,1:23]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.scy.bs[,63:72]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.scy.bs[,63:72]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.scy.bs[,63:72]),1,mean)/apply(exp(mu_loga.scy.bs[,1:23]),1,mean))-1
HDInterval::hdi(apply(exp(mu_loga.scy.bs[,63:72]),1,mean)/apply(exp(mu_loga.scy.bs[,1:23]),1,mean)-1,credMass=0.9)

#chum
log.a.jt.chm=tv.chm$draws(variable='log_a_t',format='draws_matrix')
chm.sc=which(chm.info$ocean.basin=='SC')
chm.nc=which(chm.info$ocean.basin=='NC')
chm.goa=which(chm.info$ocean.basin=='GoA')
chm.bs=which(chm.info$ocean.basin=='BS')

logat=log.a.jt.chm
colnames(logat)=gsub("]", ".", colnames(logat))
colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)

logat.sc=logat[,grepl(paste(paste(',',chm.sc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.nc=logat[,grepl(paste(paste(',',chm.nc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.goa=logat[,grepl(paste(paste(',',chm.goa,'.',sep=''),collapse = "|"),colnames(logat))]
logat.bs=logat[,grepl(paste(paste(',',chm.bs,'.',sep=''),collapse = "|"),colnames(logat))]

L=max(chm.info$end)-min(chm.info$begin)+1
mu_loga.chm.sc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.chm.nc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.chm.goa=matrix(nrow=nrow(logat),ncol=L)
mu_loga.chm.bs=matrix(nrow=nrow(logat),ncol=L)
for(t in 1:L){
  logas.sc=logat.sc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.sc))]
  mu_loga.chm.sc[,t]=apply(logas.sc,1,mean)
  
  logas.nc=logat.nc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.nc))]
  mu_loga.chm.nc[,t]=apply(logas.nc,1,mean)
  
  logas.goa=logat.goa[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.goa))]
  mu_loga.chm.goa[,t]=apply(logas.goa,1,mean)
  
  logas.bs=logat.bs[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.bs))]
  mu_loga.chm.bs[,t]=apply(logas.bs,1,mean)
}

#south coast

m.log.chm.sc=apply(exp(mu_loga.chm.sc),2,median)
names(m.log.chm.sc)=paste(seq(1:length(m.log.chm.sc)),seq(min(chm.info$begin),max(chm.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.chm.sc[,1:32]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chm.sc[,1:32]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.chm.sc[,58:68]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chm.sc[,58:68]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.chm.sc[,58:68]),1,mean)/apply(exp(mu_loga.chm.sc[,1:31]),1,mean))-1
HDInterval::hdi(apply(exp(mu_loga.chm.sc[,58:68]),1,mean)/apply(exp(mu_loga.chm.sc[,1:32]),1,mean)-1,credMass=0.9)

#north coast
m.log.chm.nc=apply(exp(mu_loga.chm.nc),2,median)
names(m.log.chm.nc)=paste(seq(1:length(m.log.chm.nc)),seq(min(chm.info$begin),max(chm.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.chm.nc[,1:32]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chm.nc[,1:32]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.chm.nc[,58:68]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chm.nc[,58:68]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.chm.nc[,58:68]),1,mean)/apply(exp(mu_loga.chm.nc[,1:32]),1,mean)-1)
HDInterval::hdi(apply(exp(mu_loga.chm.nc[,58:68]),1,mean)/apply(exp(mu_loga.chm.nc[,1:32]),1,mean)-1,credMass=0.9)

#gulf of alaska
m.log.chm.goa=apply(exp(mu_loga.chm.goa),2,median)
names(m.log.chm.goa)=paste(seq(1:length(m.log.chm.goa)),seq(min(chm.info$begin),max(chm.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.chm.goa[,1:32]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chm.goa[,1:32]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.chm.goa[,58:68]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chm.goa[,58:68]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.chm.goa[,58:68]),1,mean)/apply(exp(mu_loga.chm.goa[,1:32]),1,mean))-1
HDInterval::hdi(apply(exp(mu_loga.chm.goa[,58:68]),1,mean)/apply(exp(mu_loga.chm.goa[,1:32]),1,mean)-1,credMass=0.9)

#Bering sea
m.log.chm.bs=apply(exp(mu_loga.chm.bs),2,median)
names(m.log.chm.bs)=paste(seq(1:length(m.log.chm.bs)),seq(min(chm.info$begin),max(chm.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.chm.bs[,1:32]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chm.bs[,1:32]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.chm.bs[,58:68]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.chm.bs[,58:68]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.chm.bs[,58:68]),1,mean)/apply(exp(mu_loga.chm.bs[,1:32]),1,mean))-1
HDInterval::hdi(apply(exp(mu_loga.chm.bs[,58:68]),1,mean)/apply(exp(mu_loga.chm.bs[,1:32]),1,mean)-1,credMass=0.9)

#coho
log.a.jt.cho=tv.cho$draws(variable='log_a_t',format='draws_matrix')
cho.sc=which(cho.info$ocean.basin=='SC')
cho.nc=which(cho.info$ocean.basin=='NC')
cho.goa=which(cho.info$ocean.basin=='GoA')
cho.bs=which(cho.info$ocean.basin=='BS')

logat=log.a.jt.cho
colnames(logat)=gsub("]", ".", colnames(logat))
colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)

logat.sc=logat[,grepl(paste(paste(',',cho.sc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.nc=logat[,grepl(paste(paste(',',cho.nc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.goa=logat[,grepl(paste(paste(',',cho.goa,'.',sep=''),collapse = "|"),colnames(logat))]
logat.bs=logat[,grepl(paste(paste(',',cho.bs,'.',sep=''),collapse = "|"),colnames(logat))]

L=max(cho.info$end)-min(cho.info$begin)+1
mu_loga.cho.sc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.cho.nc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.cho.goa=matrix(nrow=nrow(logat),ncol=L)
mu_loga.cho.bs=matrix(nrow=nrow(logat),ncol=L)
for(t in 1:L){
  logas.sc=logat.sc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.sc))]
  mu_loga.cho.sc[,t]=apply(logas.sc,1,mean)
  
  logas.nc=logat.nc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.nc))]
  mu_loga.cho.nc[,t]=apply(logas.nc,1,mean)
  
  logas.goa=logat.goa[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.goa))]
  mu_loga.cho.goa[,t]=apply(logas.goa,1,mean)
  
  logas.bs=logat.bs[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.bs))]
  mu_loga.cho.bs[,t]=apply(logas.bs,1,mean)
}

#south coast

m.log.cho.sc=apply(exp(mu_loga.cho.sc),2,median)
names(m.log.cho.sc)=paste(seq(1:length(m.log.cho.sc)),seq(min(cho.info$begin),max(cho.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.cho.sc[,1:20]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.cho.sc[,1:20]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.cho.sc[,37:46]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.cho.sc[,36:46]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.cho.sc[,37:46]),1,mean)/apply(exp(mu_loga.cho.sc[,1:20]),1,mean))-1
HDInterval::hdi(apply(exp(mu_loga.cho.sc[,36:46]),1,mean)/apply(exp(mu_loga.cho.sc[,1:20]),1,mean)-1,credMass=0.9)

#north coast
m.log.cho.nc=apply(exp(mu_loga.cho.nc),2,median)
names(m.log.cho.nc)=paste(seq(1:length(m.log.cho.nc)),seq(min(cho.info$begin),max(cho.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.cho.nc[,1:18]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.cho.nc[,1:18]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.cho.nc[,37:46]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.cho.nc[,37:46]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.cho.nc[,37:46]),1,mean)/apply(exp(mu_loga.cho.nc[,1:18]),1,mean)-1)
HDInterval::hdi(apply(exp(mu_loga.cho.nc[,37:46]),1,mean)/apply(exp(mu_loga.cho.nc[,1:18]),1,mean)-1,credMass=0.9)

#gulf of alaska
m.log.cho.goa=apply(exp(mu_loga.cho.goa),2,median)
names(m.log.cho.goa)=paste(seq(1:length(m.log.cho.goa)),seq(min(cho.info$begin),max(cho.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.cho.goa[,1:31]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.cho.goa[,1:31]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.cho.goa[,58:68]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.cho.goa[,58:68]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.cho.goa[,58:68]),1,mean)/apply(exp(mu_loga.cho.goa[,1:31]),1,mean))-1
HDInterval::hdi(apply(exp(mu_loga.cho.goa[,58:68]),1,mean)/apply(exp(mu_loga.cho.goa[,1:31]),1,mean)-1,credMass=0.9)


#Pink -even
log.a.jt.pke=tv.pke$draws(variable='log_a_t',format='draws_matrix')
pke.sc=which(pke.info$ocean.basin=='SC')
pke.nc=which(pke.info$ocean.basin=='NC')
pke.goa=which(pke.info$ocean.basin=='GoA')
pke.bs=which(pke.info$ocean.basin=='BS')

logat=log.a.jt.pke
colnames(logat)=gsub("]", ".", colnames(logat))
colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)

logat.sc=logat[,grepl(paste(paste(',',pke.sc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.nc=logat[,grepl(paste(paste(',',pke.nc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.goa=logat[,grepl(paste(paste(',',pke.goa,'.',sep=''),collapse = "|"),colnames(logat))]
logat.bs=logat[,grepl(paste(paste(',',pke.bs,'.',sep=''),collapse = "|"),colnames(logat))]

pke.dat.sc=pke.dat[pke.dat$stock%in%subset(pke.info,ocean.basin=='SC')$stock.name,]
pke.dat.nc=pke.dat[pke.dat$stock%in%subset(pke.info,ocean.basin=='NC')$stock.name,]
pke.dat.goa=pke.dat[pke.dat$stock%in%subset(pke.info,ocean.basin=='GoA')$stock.name,]
pke.dat.bs=pke.dat[pke.dat$stock%in%subset(pke.info,ocean.basin=='BS')$stock.name,]

L=max(pke.info$end)-min(pke.info$begin)+1
mu_loga.pke.sc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.pke.nc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.pke.goa=matrix(nrow=nrow(logat),ncol=L)
mu_loga.pke.bs=matrix(nrow=nrow(logat),ncol=L)
for(t in 1:L){
  logas.sc=logat.sc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.sc))]
  mu_loga.pke.sc[,t]=apply(logas.sc,1,mean)
  
  logas.nc=logat.nc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.nc))]
  mu_loga.pke.nc[,t]=apply(logas.nc,1,mean)
  
  logas.goa=logat.goa[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.goa))]
  mu_loga.pke.goa[,t]=apply(logas.goa,1,mean)
  
  logas.bs=logat.bs[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.bs))]
  mu_loga.pke.bs[,t]=apply(logas.bs,1,mean)
}

#south coast

m.log.pke.sc=apply(exp(mu_loga.pke.sc),2,median)
names(m.log.pke.sc)=paste(seq(1:length(m.log.pke.sc)),seq(min(pke.info$begin),max(pke.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.pke.sc[,1:17]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.pke.sc[,1:17]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.pke.sc[,57:71]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.pke.sc[,57:71]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.pke.sc[,57:71]),1,mean)/apply(exp(mu_loga.pke.sc[,1:17]),1,mean))-1
HDInterval::hdi(apply(exp(mu_loga.pke.sc[,57:71]),1,mean)/apply(exp(mu_loga.pke.sc[,1:17]),1,mean)-1,credMass=0.9)

#north coast
m.log.pke.nc=apply(exp(mu_loga.pke.nc),2,median)
names(m.log.pke.nc)=paste(seq(1:length(m.log.pke.nc)),seq(min(pke.info$begin),max(pke.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.pke.nc[,1:19]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.pke.nc[,1:19]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.pke.nc[,57:71]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.pke.nc[,57:71]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.pke.nc[,57:71]),1,mean)/apply(exp(mu_loga.pke.nc[,1:19]),1,mean)-1)
HDInterval::hdi(apply(exp(mu_loga.pke.nc[,57:71]),1,mean)/apply(exp(mu_loga.pke.nc[,1:19]),1,mean)-1,credMass=0.9)

#gulf of alaska
m.log.pke.goa=apply(exp(mu_loga.pke.goa),2,median)
names(m.log.pke.goa)=paste(seq(1:length(m.log.pke.goa)),seq(min(pke.info$begin),max(pke.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.pke.goa[,1:31]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.pke.goa[,1:31]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.pke.goa[,57:71]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.pke.goa[,57:71]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.pke.goa[,57:71]),1,mean)/apply(exp(mu_loga.pke.goa[,1:31]),1,mean))-1
HDInterval::hdi(apply(exp(mu_loga.pke.goa[,57:71]),1,mean)/apply(exp(mu_loga.pke.goa[,1:31]),1,mean)-1,credMass=0.9)

#Bering sea
m.log.pke.bs=apply(exp(mu_loga.pke.bs),2,median)
names(m.log.pke.bs)=paste(seq(1:length(m.log.pke.bs)),seq(min(pke.info$begin),max(pke.info$end)))
#prod - 1954:1986
median(apply(exp(mu_loga.pke.bs[,1:27]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.pke.bs[,1:27]),1,mean),credMass=0.9)
#prod - 2010:2020
median(apply(exp(mu_loga.pke.bs[,57:67]),1,mean))
HDInterval::hdi(apply(exp(mu_loga.pke.bs[,57:67]),1,mean),credMass=0.9)
#Percent decline
median(apply(exp(mu_loga.pke.bs[,57:67]),1,mean)/apply(exp(mu_loga.pke.bs[,1:31]),1,mean))-1
HDInterval::hdi(apply(exp(mu_loga.pke.bs[,57:67]),1,mean)/apply(exp(mu_loga.pke.bs[,1:31]),1,mean)-1,credMass=0.9)

#harvest rates####
loga.t.chk=tv.chk$draws(variables='log_a_t',format='draws_matrix')
umsy.st.chk=stc.chk$draws(variables='Umsy',format='draws_matrix')

umsy_comp_chk=data.frame(stock=chk.info$stock.name,st.umsy=NA,st.umsy.l90=NA,st.umsy.u90=NA,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,p.change=NA)
for(i in 1:nrow(chk.info)){
  logas=loga.t.chk[,grepl(paste(',',i,']',sep=''),colnames(loga.t.chk))]
  umsy.st=umsy.st.chk[,i]
  umsy.st=ifelse(umsy.st<0,0.005,umsy.st)
  
  o.logas=apply(logas[,c(chk.info$begin[i]-min(chk.info$begin)+1):c(chk.info$begin[i]-min(chk.info$begin)+5)],1,mean)
  r.logas=apply(logas[,c(chk.info$end[i]-min(chk.info$begin)+1-5):c(chk.info$end[i]-min(chk.info$begin))],1,mean)
  r.umsy=sapply(r.logas,umsyCalc)
  r.umsy=ifelse(r.umsy<0,0.005,r.umsy)
  
  umsy_comp_chk[i,2]=median(umsy.st)
  umsy_comp_chk[i,3]=HDInterval::hdi(umsy.st,credMass=0.9)[1]
  umsy_comp_chk[i,4]=HDInterval::hdi(umsy.st,credMass=0.9)[2]
  umsy_comp_chk[i,5]=median(r.umsy)
  umsy_comp_chk[i,6]=HDInterval::hdi(r.umsy,credMass=0.9)[1]
  umsy_comp_chk[i,7]=HDInterval::hdi(r.umsy,credMass=0.9)[2]
  umsy_comp_chk[i,8]=median((r.umsy/umsy.st)-1)
}

loga.t.scy=tv.scy$draws(variables='log_a_t',format='draws_matrix')
umsy.st.scy=stc.scy$draws(variables='Umsy',format='draws_matrix')

umsy_comp_scy=data.frame(stock=scy.info$stock.name,st.umsy=NA,st.umsy.l90=NA,st.umsy.u90=NA,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,p.change=NA)
for(i in 1:nrow(scy.info)){
  logas=loga.t.scy[,grepl(paste(',',i,']',sep=''),colnames(loga.t.scy))]
  umsy.st=umsy.st.scy[,i]
  umsy.st=ifelse(umsy.st<0,0.005,umsy.st)
  
  o.logas=apply(logas[,c(scy.info$begin[i]-min(scy.info$begin)+1):c(scy.info$begin[i]-min(scy.info$begin)+5)],1,mean)
  r.logas=apply(logas[,c(scy.info$end[i]-min(scy.info$begin)+1-5):c(scy.info$end[i]-min(scy.info$begin))],1,mean)
  r.umsy=sapply(r.logas,umsyCalc)
  r.umsy=ifelse(r.umsy<0,0.005,r.umsy)
  
  umsy_comp_scy[i,2]=median(umsy.st)
  umsy_comp_scy[i,3]=HDInterval::hdi(umsy.st,credMass=0.9)[1]
  umsy_comp_scy[i,4]=HDInterval::hdi(umsy.st,credMass=0.9)[2]
  umsy_comp_scy[i,5]=median(r.umsy)
  umsy_comp_scy[i,6]=HDInterval::hdi(r.umsy,credMass=0.9)[1]
  umsy_comp_scy[i,7]=HDInterval::hdi(r.umsy,credMass=0.9)[2]
  umsy_comp_scy[i,8]=median((r.umsy/umsy.st)-1)
}

loga.t.chm=tv.chm$draws(variables='log_a_t',format='draws_matrix')
umsy.st.chm=stc.chm$draws(variables='Umsy',format='draws_matrix')

umsy_comp_chm=data.frame(stock=chm.info$stock.name,st.umsy=NA,st.umsy.l90=NA,st.umsy.u90=NA,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,p.change=NA)
for(i in 1:nrow(chm.info)){
  logas=loga.t.chm[,grepl(paste(',',i,']',sep=''),colnames(loga.t.chm))]
  umsy.st=umsy.st.chm[,i]
  umsy.st=ifelse(umsy.st<0,0.005,umsy.st)
  
  o.logas=apply(logas[,c(chm.info$begin[i]-min(chm.info$begin)+1):c(chm.info$begin[i]-min(chm.info$begin)+5)],1,mean)
  r.logas=apply(logas[,c(chm.info$end[i]-min(chm.info$begin)+1-5):c(chm.info$end[i]-min(chm.info$begin))],1,mean)
  r.umsy=sapply(r.logas,umsyCalc)
  r.umsy=ifelse(r.umsy<0,0.005,r.umsy)
  
  umsy_comp_chm[i,2]=median(umsy.st)
  umsy_comp_chm[i,3]=HDInterval::hdi(umsy.st,credMass=0.9)[1]
  umsy_comp_chm[i,4]=HDInterval::hdi(umsy.st,credMass=0.9)[2]
  umsy_comp_chm[i,5]=median(r.umsy)
  umsy_comp_chm[i,6]=HDInterval::hdi(r.umsy,credMass=0.9)[1]
  umsy_comp_chm[i,7]=HDInterval::hdi(r.umsy,credMass=0.9)[2]
  umsy_comp_chm[i,8]=median((r.umsy/umsy.st)-1)
}

loga.t.cho=tv.cho$draws(variables='log_a_t',format='draws_matrix')
umsy.st.cho=stc.cho$draws(variables='Umsy',format='draws_matrix')

umsy_comp_cho=data.frame(stock=cho.info$stock.name,st.umsy=NA,st.umsy.l90=NA,st.umsy.u90=NA,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,p.change=NA)
for(i in 1:nrow(cho.info)){
  logas=loga.t.cho[,grepl(paste(',',i,']',sep=''),colnames(loga.t.cho))]
  umsy.st=umsy.st.cho[,i]
  umsy.st=ifelse(umsy.st<0,0.005,umsy.st)
  
  o.logas=apply(logas[,c(cho.info$begin[i]-min(cho.info$begin)+1):c(cho.info$begin[i]-min(cho.info$begin)+5)],1,mean)
  r.logas=apply(logas[,c(cho.info$end[i]-min(cho.info$begin)+1-5):c(cho.info$end[i]-min(cho.info$begin))],1,mean)
  r.umsy=sapply(r.logas,umsyCalc)
  r.umsy=ifelse(r.umsy<0,0.005,r.umsy)
  
  umsy_comp_cho[i,2]=median(umsy.st)
  umsy_comp_cho[i,3]=HDInterval::hdi(umsy.st,credMass=0.9)[1]
  umsy_comp_cho[i,4]=HDInterval::hdi(umsy.st,credMass=0.9)[2]
  umsy_comp_cho[i,5]=median(r.umsy)
  umsy_comp_cho[i,6]=HDInterval::hdi(r.umsy,credMass=0.9)[1]
  umsy_comp_cho[i,7]=HDInterval::hdi(r.umsy,credMass=0.9)[2]
  umsy_comp_cho[i,8]=median((r.umsy/umsy.st)-1)
}

loga.t.pke=tv.pke$draws(variables='log_a_t',format='draws_matrix')
umsy.st.pke=stc.pke$draws(variables='Umsy',format='draws_matrix')

umsy_comp_pke=data.frame(stock=pke.info$stock.name,st.umsy=NA,st.umsy.l90=NA,st.umsy.u90=NA,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,p.change=NA)
for(i in 1:nrow(pke.info)){
  logas=loga.t.pke[,grepl(paste(',',i,']',sep=''),colnames(loga.t.pke))]
  umsy.st=umsy.st.pke[,i]
  umsy.st=ifelse(umsy.st<0,0.005,umsy.st)
  
  o.logas=apply(logas[,c(pke.info$begin[i]-min(pke.info$begin)+1):c(pke.info$begin[i]-min(pke.info$begin)+5)],1,mean)
  r.logas=apply(logas[,c(pke.info$end[i]-min(pke.info$begin)+1-5):c(pke.info$end[i]-min(pke.info$begin))],1,mean)
  r.umsy=sapply(r.logas,umsyCalc)
  r.umsy=ifelse(r.umsy<0,0.005,r.umsy)
  
  umsy_comp_pke[i,2]=median(umsy.st)
  umsy_comp_pke[i,3]=HDInterval::hdi(umsy.st,credMass=0.9)[1]
  umsy_comp_pke[i,4]=HDInterval::hdi(umsy.st,credMass=0.9)[2]
  umsy_comp_pke[i,5]=median(r.umsy)
  umsy_comp_pke[i,6]=HDInterval::hdi(r.umsy,credMass=0.9)[1]
  umsy_comp_pke[i,7]=HDInterval::hdi(r.umsy,credMass=0.9)[2]
  umsy_comp_pke[i,8]=median((r.umsy/umsy.st)-1)
}

loga.t.pko=tv.pko$draws(variables='log_a_t',format='draws_matrix')
umsy.st.pko=stc.pko$draws(variables='Umsy',format='draws_matrix')

umsy_comp_pko=data.frame(stock=pko.info$stock.name,st.umsy=NA,st.umsy.l90=NA,st.umsy.u90=NA,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,p.change=NA)
for(i in 1:nrow(pko.info)){
  logas=loga.t.pko[,grepl(paste(',',i,']',sep=''),colnames(loga.t.pko))]
  umsy.st=umsy.st.pko[,i]
  umsy.st=ifelse(umsy.st<0,0.005,umsy.st)
  
  o.logas=apply(logas[,c(pko.info$begin[i]-min(pko.info$begin)+1):c(pko.info$begin[i]-min(pko.info$begin)+5)],1,mean)
  r.logas=apply(logas[,c(pko.info$end[i]-min(pko.info$begin)+1-5):c(pko.info$end[i]-min(pko.info$begin))],1,mean)
  r.umsy=sapply(r.logas,umsyCalc)
  r.umsy=ifelse(r.umsy<0,0.005,r.umsy)
  
  umsy_comp_pko[i,2]=median(umsy.st)
  umsy_comp_pko[i,3]=HDInterval::hdi(umsy.st,credMass=0.9)[1]
  umsy_comp_pko[i,4]=HDInterval::hdi(umsy.st,credMass=0.9)[2]
  umsy_comp_pko[i,5]=median(r.umsy)
  umsy_comp_pko[i,6]=HDInterval::hdi(r.umsy,credMass=0.9)[1]
  umsy_comp_pko[i,7]=HDInterval::hdi(r.umsy,credMass=0.9)[2]
  umsy_comp_pko[i,8]=median((r.umsy/umsy.st)-1)
}

nrow(subset(umsy_comp_chk,p.change< -0.25))/nrow(umsy_comp_chk)
nrow(subset(umsy_comp_scy,p.change< -0.25))/nrow(umsy_comp_scy)
nrow(subset(umsy_comp_chm,p.change< -0.25))/nrow(umsy_comp_chm)
nrow(subset(umsy_comp_cho,p.change< -0.25))/nrow(umsy_comp_cho)
nrow(subset(umsy_comp_pke,p.change< -0.25))/nrow(umsy_comp_pke)
nrow(subset(umsy_comp_pko,p.change< -0.25))/nrow(umsy_comp_pko)

nrow(subset(umsy_comp_chk,r.umsy<0.05))/nrow(umsy_comp_chk)
nrow(subset(umsy_comp_scy,r.umsy<0.05))/nrow(umsy_comp_scy)
nrow(subset(umsy_comp_chm,r.umsy<0.05))/nrow(umsy_comp_chm)
nrow(subset(umsy_comp_cho,r.umsy<0.05))/nrow(umsy_comp_cho)
nrow(subset(umsy_comp_pke,r.umsy<0.05))/nrow(umsy_comp_pke)
nrow(subset(umsy_comp_pko,r.umsy<0.05))/nrow(umsy_comp_pko)


#stationary - table 1####
stc.chk=readRDS(here('outputs','model fits','static_fit_chinook.RDS'))
stc.cho=readRDS(here('outputs','model fits','static_fit_coho.RDS'))
stc.chm=readRDS(here('outputs','model fits','static_fit_chum.RDS'))
stc.pke=readRDS(here('outputs','model fits','static_fit_pinke.RDS'))
stc.pko=readRDS(here('outputs','model fits','static_fit_pinko.RDS'))
stc.scy=readRDS(here('outputs','model fits','static_fit_sockeye.RDS'))

#productivity
#coastwide
glob.alpha.chk=stc.chk$draws(variable='log_a0',format='draws_matrix')
glob.alpha.chm=stc.chm$draws(variable='log_a0',format='draws_matrix')
glob.alpha.cho=stc.cho$draws(variable='log_a0',format='draws_matrix')
glob.alpha.pke=stc.pke$draws(variable='log_a0',format='draws_matrix')
glob.alpha.pko=stc.pko$draws(variable='log_a0',format='draws_matrix')
glob.alpha.scy=stc.scy$draws(variable='log_a0',format='draws_matrix')

#chk
mean(exp(glob.alpha.chk))
sd(exp(glob.alpha.chk))
mean(glob.alpha.chk)
sd(glob.alpha.chk)
#chm
mean(exp(glob.alpha.chm))
sd(exp(glob.alpha.chm))
mean(glob.alpha.chm)
sd(glob.alpha.chm)

#cho
mean(exp(glob.alpha.cho))
sd(exp(glob.alpha.cho))
mean(glob.alpha.cho)
sd(glob.alpha.cho)

#pke
mean(exp(glob.alpha.pke))
sd(exp(glob.alpha.pke))
mean(glob.alpha.pke)
sd(glob.alpha.pke)

#pko
mean(exp(glob.alpha.pko))
sd(exp(glob.alpha.pko))
mean(glob.alpha.pko)
sd(glob.alpha.pko)

#scy
mean(exp(glob.alpha.scy))
sd(exp(glob.alpha.scy))
mean(glob.alpha.scy)
sd(glob.alpha.scy)

#by ocean region
lme.alpha.chk1=exp(stc.chk$draws(variable='log_a_or',format='draws_matrix'))
lme.alpha.chk2=stc.chk$draws(variable='log_a_or',format='draws_matrix')
lme.alpha.chm1=exp(stc.chm$draws(variable='log_a_or',format='draws_matrix'))
lme.alpha.chm2=stc.chm$draws(variable='log_a_or',format='draws_matrix')
lme.alpha.cho1=exp(stc.cho$draws(variable='log_a_or',format='draws_matrix'))
lme.alpha.cho2=stc.cho$draws(variable='log_a_or',format='draws_matrix')
lme.alpha.pke1=exp(stc.pke$draws(variable='log_a_or',format='draws_matrix'))
lme.alpha.pke2=stc.pke$draws(variable='log_a_or',format='draws_matrix')
lme.alpha.pko1=exp(stc.pko$draws(variable='log_a_or',format='draws_matrix'))
lme.alpha.pko2=stc.pko$draws(variable='log_a_or',format='draws_matrix')
lme.alpha.scy1=exp(stc.scy$draws(variable='log_a_or',format='draws_matrix'))
lme.alpha.scy2=stc.scy$draws(variable='log_a_or',format='draws_matrix')

lme.n = stk.info %>% group_by(species,ocean.basin) %>% summarize(n=n())

data.frame(species=c(rep('Chinook',4),rep('Chum',4),rep('Coho',3),rep('Pink (even)',4),rep('Pink (odd)',4),rep('Sockeye',4)),ocean.region=c(rep(c('Bering Sea','Gulf of Alaska','North Coast','South Coast'),2),c('Gulf of Alaska','North Coast','South Coast'),rep(c('Bering Sea','Gulf of Alaska','North Coast','South Coast'),3)),
           n=lme.n[,3],
           a.mean=c(round(apply(lme.alpha.chk1,2,mean),2),round(apply(lme.alpha.chm1,2,mean),2),round(apply(lme.alpha.cho1,2,mean),2),round(apply(lme.alpha.pke1,2,mean),2),round(apply(lme.alpha.pko1,2,mean),2),round(apply(lme.alpha.scy1,2,mean),2)),
           a.sd=c(round(apply(lme.alpha.chk1,2,sd),2),round(apply(lme.alpha.chm1,2,sd),2),round(apply(lme.alpha.cho1,2,sd),2),round(apply(lme.alpha.pke1,2,sd),2),round(apply(lme.alpha.pko1,2,sd),2),round(apply(lme.alpha.scy1,2,sd),2)),
           loga.mean=c(round(apply(lme.alpha.chk2,2,mean),2),round(apply(lme.alpha.chm2,2,mean),2),round(apply(lme.alpha.cho2,2,mean),2),round(apply(lme.alpha.pke2,2,mean),2),round(apply(lme.alpha.pko2,2,mean),2),round(apply(lme.alpha.scy2,2,mean),2)),
           loga.sd=c(round(apply(lme.alpha.chk2,2,sd),2),round(apply(lme.alpha.chm2,2,sd),2),round(apply(lme.alpha.cho2,2,sd),2),round(apply(lme.alpha.pke2,2,sd),2),round(apply(lme.alpha.pko2,2,sd),2),round(apply(lme.alpha.scy2,2,sd),2)))


#sigma
glob.sigma.chk=stc.chk$draws(variable='sigma',format='draws_matrix')
glob.sigma.chm=stc.chm$draws(variable='sigma',format='draws_matrix')
glob.sigma.cho=stc.cho$draws(variable='sigma',format='draws_matrix')
glob.sigma.pke=stc.pke$draws(variable='sigma',format='draws_matrix')
glob.sigma.pko=stc.pko$draws(variable='sigma',format='draws_matrix')
glob.sigma.scy=stc.scy$draws(variable='sigma',format='draws_matrix')

#chk
mean(glob.sigma.chk)
sd(glob.sigma.chk)
#chm
mean(glob.sigma.chm)
sd(glob.sigma.chm)
#cho
mean(glob.sigma.cho)
sd(glob.sigma.cho)
#pke
mean(glob.sigma.pke)
sd(glob.sigma.pke)
#pko
mean(glob.sigma.pko)
sd(glob.sigma.pko)
#scy
mean(glob.sigma.scy)
sd(glob.sigma.scy)

#
chk.sc=which(chk.info$ocean.basin=='SC')
chm.sc=which(chm.info$ocean.basin=='SC')
cho.sc=which(cho.info$ocean.basin=='SC')
pke.sc=which(pke.info$ocean.basin=='SC')
pko.sc=which(pko.info$ocean.basin=='SC')
scy.sc=which(scy.info$ocean.basin=='SC')
#nbc-seak
chk.nc=which(chk.info$ocean.basin=='NC')
chm.nc=which(chm.info$ocean.basin=='NC')
cho.nc=which(cho.info$ocean.basin=='NC')
pke.nc=which(pke.info$ocean.basin=='NC')
pko.nc=which(pko.info$ocean.basin=='NC')
scy.nc=which(scy.info$ocean.basin=='NC')
#gulf
chk.goa=which(chk.info$ocean.basin=='GoA')
chm.goa=which(chm.info$ocean.basin=='GoA')
cho.goa=which(cho.info$ocean.basin=='GoA')
pke.goa=which(pke.info$ocean.basin=='GoA')
pko.goa=which(pko.info$ocean.basin=='GoA')
scy.goa=which(scy.info$ocean.basin=='GoA')
#bering
chk.bs=which(chk.info$ocean.basin=='BS')
chm.bs=which(chm.info$ocean.basin=='BS')
cho.bs=which(cho.info$ocean.basin=='BS')
pke.bs=which(pke.info$ocean.basin=='BS')
pko.bs=which(pko.info$ocean.basin=='BS')
scy.bs=which(scy.info$ocean.basin=='BS')


sigma.chk.sc=glob.sigma.chk[,chk.sc];sigma.scy.sc=glob.sigma.scy[,scy.sc];sigma.chm.sc=glob.sigma.chm[,chm.sc];sigma.cho.sc=glob.sigma.cho[,cho.sc];sigma.pke.sc=glob.sigma.pke[,pke.sc];sigma.pko.sc=glob.sigma.pko[pko.sc]
sigma.chk.nc=glob.sigma.chk[,chk.nc];sigma.scy.nc=glob.sigma.scy[,scy.nc];sigma.chm.nc=glob.sigma.chm[,chm.nc];sigma.cho.nc=glob.sigma.cho[,cho.nc];sigma.pke.nc=glob.sigma.pke[,pke.nc];sigma.pko.nc=glob.sigma.pko[pko.nc]
sigma.chk.goa=glob.sigma.chk[,chk.goa];sigma.scy.goa=glob.sigma.scy[,scy.goa];sigma.chm.goa=glob.sigma.chm[,chm.goa];sigma.cho.goa=glob.sigma.cho[,cho.goa];sigma.pke.goa=glob.sigma.pke[,pke.goa];sigma.pko.goa=glob.sigma.pko[pko.goa]
sigma.chk.bs=glob.sigma.chk[,chk.bs];sigma.scy.bs=glob.sigma.scy[,scy.bs];sigma.chm.bs=glob.sigma.chm[,chm.bs];sigma.cho.bs=glob.sigma.cho[,cho.bs];sigma.pke.bs=glob.sigma.pke[,pke.bs];sigma.pko.bs=glob.sigma.pko[pko.bs]

data.frame(species=c(rep('Chinook',4),rep('Chum',4),rep('Coho',3),rep('Pink (even)',4),rep('Pink (odd)',4),rep('Sockeye',4)),ocean.region=c(rep(c('Bering Sea','Gulf of Alaska','North Coast','South Coast'),2),c('Gulf of Alaska','North Coast','South Coast'),rep(c('Bering Sea','Gulf of Alaska','North Coast','South Coast'),3)),
             n=lme.n[,3],
              mean=c(round(mean(sigma.chk.bs),2),round(mean(sigma.chk.goa),2),round(mean(sigma.chk.nc),2),round(mean(sigma.chk.sc),2),
                     round(mean(sigma.chm.bs),2),round(mean(sigma.chm.goa),2),round(mean(sigma.chm.nc),2),round(mean(sigma.chm.sc),2),
                     round(mean(sigma.cho.goa),2),round(mean(sigma.cho.nc),2),round(mean(sigma.cho.sc),2),
                     round(mean(sigma.pke.bs),2),round(mean(sigma.pke.goa),2),round(mean(sigma.pke.nc),2),round(mean(sigma.pke.sc),2),
                     round(mean(sigma.pko.bs),2),round(mean(sigma.pko.goa),2),round(mean(sigma.pko.nc),2),round(mean(sigma.pko.sc),2),
                     round(mean(sigma.scy.bs),2),round(mean(sigma.scy.goa),2),round(mean(sigma.scy.nc),2),round(mean(sigma.scy.sc),2)
                     ),
           
              sd=c(round(sd(sigma.chk.bs),2),round(sd(sigma.chk.goa),2),round(sd(sigma.chk.nc),2),round(sd(sigma.chk.sc),2),
                         round(sd(sigma.chm.bs),2),round(sd(sigma.chm.goa),2),round(sd(sigma.chm.nc),2),round(sd(sigma.chm.sc),2),
                         round(sd(sigma.cho.goa),2),round(sd(sigma.cho.nc),2),round(sd(sigma.cho.sc),2),
                         round(sd(sigma.pke.bs),2),round(sd(sigma.pke.goa),2),round(sd(sigma.pke.nc),2),round(sd(sigma.pke.sc),2),
                         round(sd(sigma.pko.bs),2),round(sd(sigma.pko.goa),2),round(sd(sigma.pko.nc),2),round(sd(sigma.pko.sc),2),
                         round(sd(sigma.scy.bs),2),round(sd(sigma.scy.goa),2),round(sd(sigma.scy.nc),2),round(sd(sigma.scy.sc),2)
              )
)


glob.rho.chk=stc.chk$draws(variable='rho',format='draws_matrix')
glob.rho.chm=stc.chm$draws(variable='rho',format='draws_matrix')
glob.rho.cho=stc.cho$draws(variable='rho',format='draws_matrix')
glob.rho.pke=stc.pke$draws(variable='rho',format='draws_matrix')
glob.rho.pko=stc.pko$draws(variable='rho',format='draws_matrix')
glob.rho.scy=stc.scy$draws(variable='rho',format='draws_matrix')

#chk
mean(glob.rho.chk)
sd(glob.rho.chk)
#chm
mean(glob.rho.chm)
sd(glob.rho.chm)
#cho
mean(glob.rho.cho)
sd(glob.rho.cho)
#pke
mean(glob.rho.pke)
sd(glob.rho.pke)
#pko
mean(glob.rho.pko)
sd(glob.rho.pko)
#scy
mean(glob.rho.scy)
sd(glob.rho.scy)

rho.chk.sc=glob.rho.chk[,chk.sc];rho.scy.sc=glob.rho.scy[,scy.sc];rho.chm.sc=glob.rho.chm[,chm.sc];rho.cho.sc=glob.rho.cho[,cho.sc];rho.pke.sc=glob.rho.pke[,pke.sc];rho.pko.sc=glob.rho.pko[pko.sc]
rho.chk.nc=glob.rho.chk[,chk.nc];rho.scy.nc=glob.rho.scy[,scy.nc];rho.chm.nc=glob.rho.chm[,chm.nc];rho.cho.nc=glob.rho.cho[,cho.nc];rho.pke.nc=glob.rho.pke[,pke.nc];rho.pko.nc=glob.rho.pko[pko.nc]
rho.chk.goa=glob.rho.chk[,chk.goa];rho.scy.goa=glob.rho.scy[,scy.goa];rho.chm.goa=glob.rho.chm[,chm.goa];rho.cho.goa=glob.rho.cho[,cho.goa];rho.pke.goa=glob.rho.pke[,pke.goa];rho.pko.goa=glob.rho.pko[pko.goa]
rho.chk.bs=glob.rho.chk[,chk.bs];rho.scy.bs=glob.rho.scy[,scy.bs];rho.chm.bs=glob.rho.chm[,chm.bs];rho.cho.bs=glob.rho.cho[,cho.bs];rho.pke.bs=glob.rho.pke[,pke.bs];rho.pko.bs=glob.rho.pko[pko.bs]

data.frame(species=c(rep('Chinook',4),rep('Chum',4),rep('Coho',3),rep('Pink (even)',4),rep('Pink (odd)',4),rep('Sockeye',4)),ocean.region=c(rep(c('Bering Sea','Gulf of Alaska','North Coast','South Coast'),2),c('Gulf of Alaska','North Coast','South Coast'),rep(c('Bering Sea','Gulf of Alaska','North Coast','South Coast'),3)),
           n=lme.n[,3],
           mean=c(round(mean(rho.chk.bs),2),round(mean(rho.chk.goa),2),round(mean(rho.chk.nc),2),round(mean(rho.chk.sc),2),
                  round(mean(rho.chm.bs),2),round(mean(rho.chm.goa),2),round(mean(rho.chm.nc),2),round(mean(rho.chm.sc),2),
                  round(mean(rho.cho.goa),2),round(mean(rho.cho.nc),2),round(mean(rho.cho.sc),2),
                  round(mean(rho.pke.bs),2),round(mean(rho.pke.goa),2),round(mean(rho.pke.nc),2),round(mean(rho.pke.sc),2),
                  round(mean(rho.pko.bs),2),round(mean(rho.pko.goa),2),round(mean(rho.pko.nc),2),round(mean(rho.pko.sc),2),
                  round(mean(rho.scy.bs),2),round(mean(rho.scy.goa),2),round(mean(rho.scy.nc),2),round(mean(rho.scy.sc),2)
           ),
           
           sd=c(round(sd(rho.chk.bs),2),round(sd(rho.chk.goa),2),round(sd(rho.chk.nc),2),round(sd(rho.chk.sc),2),
                round(sd(rho.chm.bs),2),round(sd(rho.chm.goa),2),round(sd(rho.chm.nc),2),round(sd(rho.chm.sc),2),
                round(sd(rho.cho.goa),2),round(sd(rho.cho.nc),2),round(sd(rho.cho.sc),2),
                round(sd(rho.pke.bs),2),round(sd(rho.pke.goa),2),round(sd(rho.pke.nc),2),round(sd(rho.pke.sc),2),
                round(sd(rho.pko.bs),2),round(sd(rho.pko.goa),2),round(sd(rho.pko.nc),2),round(sd(rho.pko.sc),2),
                round(sd(rho.scy.bs),2),round(sd(rho.scy.goa),2),round(sd(rho.scy.nc),2),round(sd(rho.scy.sc),2)
           )
)


