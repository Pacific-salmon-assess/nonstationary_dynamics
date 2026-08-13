#main text figures


umsyCalc <- function(a) {
  #gives the min Ricker log-likelihood
  umsy <- (1 - gsl::lambert_W0(exp(1 - a)))
  
  return(umsy)
}

#data####
library(here);library(ggplot2);library(dplyr);library(cowplot);library(ggspatial);library(sf);library(rnaturalearth);library(cmdstanr)
options(dplyr.summarise.inform = FALSE)
stock.dat<- read.csv(here('data','salmon_productivity_compilation2026-08-05.csv'))
stock.info<- read.csv(here('data','stock_info2026-08-05.csv'))
stock.info$state=factor(stock.info$state,levels=c('WA','BC','AK'))
stock.info$ocean.basin=factor(stock.info$ocean.basin,levels=c('SC','NC','GoA','BS'))
stock.info$species2=gsub("-.*", "", stock.info$species) # collapse even/odd pinks together

stk.info=stock.info
chk.info=subset(stock.info,species=='Chinook')
scy.info=subset(stock.info,species=='Sockeye')
chm.info=subset(stock.info,species=='Chum')
cho.info=subset(stock.info,species=='Coho')
pke.info=subset(stock.info,species=='Pink-Even')
pko.info=subset(stock.info,species=='Pink-Odd')

#by ocean region
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


pke.sc=which(pke.info$ocean.basin=='SC')
pke.sc2=subset(pke.info,ocean.basin=='SC')
pke.nc=which(pke.info$ocean.basin=='NC')
pke.nc2=subset(pke.info,ocean.basin=='NC')
pke.goa=which(pke.info$ocean.basin=='GoA')
pke.goa2=subset(pke.info,ocean.basin=='GoA')
pke.bs=which(pke.info$ocean.basin=='BS')
pke.bs2=subset(pke.info,ocean.basin=='BS')

pke.dat=stock.dat[stock.dat$stock.id%in%pke.info$stock.id,]

pke.dat.sc=pke.dat[pke.dat$stock%in%subset(pke.info,ocean.basin=='SC')$stock.name,]
pke.dat.nc=pke.dat[pke.dat$stock%in%subset(pke.info,ocean.basin=='NC')$stock.name,]
pke.dat.goa=pke.dat[pke.dat$stock%in%subset(pke.info,ocean.basin=='GoA')$stock.name,]
pke.dat.bs=pke.dat[pke.dat$stock%in%subset(pke.info,ocean.basin=='BS')$stock.name,]

pko.dat=stock.dat[stock.dat$stock.id%in%pko.info$stock.id,]

pko.dat.sc=pko.dat[pko.dat$stock%in%subset(pko.info,ocean.basin=='SC')$stock.name,]
pko.dat.nc=pko.dat[pko.dat$stock%in%subset(pko.info,ocean.basin=='NC')$stock.name,]
pko.dat.goa=pko.dat[pko.dat$stock%in%subset(pko.info,ocean.basin=='GoA')$stock.name,]
pko.dat.bs=pko.dat[pko.dat$stock%in%subset(pko.info,ocean.basin=='BS')$stock.name,]

pko.sc=which(pko.info$ocean.basin=='SC')
pko.sc2=subset(pko.info,ocean.basin=='SC')
pko.nc=which(pko.info$ocean.basin=='NC')
pko.nc2=subset(pko.info,ocean.basin=='NC')
pko.goa=which(pko.info$ocean.basin=='GoA')
pko.goa2=subset(pko.info,ocean.basin=='GoA')
pko.bs=which(pko.info$ocean.basin=='BS')
pko.bs2=subset(pko.info,ocean.basin=='BS')


#static fits
stc.chk=readRDS(here('outputs','model fits','static_fit_chinook.RDS'))
stc.cho=readRDS(here('outputs','model fits','static_fit_coho.RDS'))
stc.chm=readRDS(here('outputs','model fits','static_fit_chum.RDS'))
stc.pke=readRDS(here('outputs','model fits','static_fit_pinke.RDS'))
stc.pko=readRDS(here('outputs','model fits','static_fit_pinko.RDS'))
stc.scy=readRDS(here('outputs','model fits','static_fit_sockeye.RDS'))

#nonstationary fits
tv.chk=readRDS(here('outputs','model fits','ns_fit_chinook.RDS'))
tv.cho=readRDS(here('outputs','model fits','ns_fit_coho.RDS'))
tv.chm=readRDS(here('outputs','model fits','ns_fit_chum.RDS'))
tv.pke=readRDS(here('outputs','model fits','ns_fit_pink_even.RDS'))
tv.pko=readRDS(here('outputs','model fits','ns_fit_pink_odd.RDS'))
tv.scy=readRDS(here('outputs','model fits','ns_fit_sockeye.RDS'))


#visual specs####
sp_cols=cbind(c(colorBlindness::SteppedSequential5Steps[16],colorBlindness::SteppedSequential5Steps[11],colorBlindness::SteppedSequential5Steps[19],colorBlindness::SteppedSequential5Steps[8],colorBlindness::SteppedSequential5Steps[1]),c('Chinook','Chum','Coho','Pink','Sockeye'))

basin_cols=cbind(c('#ffbe19',
                   '#4ca94b',
                   '#00787a',
                   '#003f5c')
                 ,c('SC','NC','GoA','BS'))


mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)


#Fig. 1 Map of time-series ####
#map of S-R time-series
world <- ne_countries(scale = "medium", returnclass = "sf")

all_lat_lon = stock.info %>% group_by(lat,lon,species2) %>% summarize(n=n(),.groups = 'keep')
all_lat_lon$lon=ifelse(all_lat_lon$lon>0,-all_lat_lon$lon,all_lat_lon$lon)
all_lat_lon$ocean.basin=stock.info$ocean.basin[match(all_lat_lon$lat,stock.info$lat)]
ll_p=subset(all_lat_lon,species2=='Pink')
ll_ch=subset(all_lat_lon,species2=='Chum')
ll_coh=subset(all_lat_lon,species2=='Coho')
ll_soc=subset(all_lat_lon,species2=='Sockeye')
ll_chi=subset(all_lat_lon,species2=='Chinook')

sp_sum<- stock.info %>% group_by(species2,ocean.basin) %>% summarize(n=n(),.groups='keep')
sp_sum$ocean.basin2=dplyr::recode_factor(sp_sum$ocean.basin,
                                         BS='BS',
                                         GOA='GoA',
                                         SEAK='NC',
                                         WC='SC')

inset=ggplot(sp_sum, aes(fill=species2, y=n, x=ocean.basin2))+
  scale_fill_manual(values=sp_cols[match(levels(factor(sp_sum$species2)),sp_cols[,2])],name='Species')+ 
  geom_bar(position="stack", stat="identity") + theme_minimal() +
  xlab('')+ylab('')+
  
  theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background  = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      strip.background = element_rect(fill = "transparent", colour = NA),
      panel.grid = element_blank(),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15),legend.position = "none",axis.text.x = element_text(
          color = c(basin_cols[4],basin_cols[3],basin_cols[2],basin_cols[1])))


main=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("Longitude") + ylab("Latitude") +
  geom_point(data = ll_p,shape=21, mapping = aes(x = lon, y = lat), fill  = alpha(sp_cols[4], 0.4), color =  alpha(sp_cols[4],0.7), size = 2.2*log(ll_p$n+1),stroke=1.1) +
  geom_point(data = ll_ch,shape=21, mapping = aes(x = lon, y = lat), fill  = alpha(sp_cols[2], 0.4), color =  alpha(sp_cols[2],0.7), size = 2.2*log(ll_ch$n+1),stroke=1.1) +
  geom_point(data =  ll_soc,shape=21, mapping = aes(x = lon, y = lat), fill  = alpha(sp_cols[5], 0.4), color =  alpha(sp_cols[5],0.7), size = 2.2*log(ll_soc$n+1),stroke=1.1) +
  geom_point(data =  ll_coh,shape=21, mapping = aes(x = lon, y = lat), fill  = alpha(sp_cols[3], 0.4), color =  alpha(sp_cols[3],0.7), size = 2.2*log(ll_coh$n+1),stroke=1.1) +
  geom_point(data = ll_chi,shape=21, mapping = aes(x = lon, y = lat),fill  = alpha(sp_cols[1], 0.4), color =  alpha(sp_cols[1],0.7), size = 2.2*log(ll_chi$n+1),stroke=1.1) +
  annotate(geom = 'text', x = -160, y = 53.5, label = 'Chinook', color = sp_cols[1], size = 4.5) + 
  annotate(geom = 'text', x = -160, y = 52, label = 'Chum', color = sp_cols[2], size = 4.5) + 
  annotate(geom = 'text', x = -160, y = 50.5, label = 'Coho', color = sp_cols[3], size = 4.5) + 
  annotate(geom = 'text', x = -160, y = 49, label = 'Pink', color = sp_cols[4], size = 4.5) + 
  annotate(geom = 'text', x = -160, y = 47.5, label = 'Sockeye', color = sp_cols[5], size = 4.5) + 
  annotate(geom = 'point', x = -154, y = 52, size = 2.2*log(1+1)) +
  annotate(geom = 'point', x = -154, y = 50.5, size = 2.2*log(3+1)) + 
  annotate(geom = 'point', x = -154, y = 49, size = 2.2*log(10+1)) + 
  annotate(geom = 'point', x = -154, y = 47.5, size = 2.2*log(20+1)) +
  annotate(geom = 'text', x = -164, y = 57.5,label='Bering Sea', size = 4.5,color=basin_cols[4,1])+
  annotate(geom = 'text', x = -150, y = 55.8,label='Gulf of Alaska', size = 4.5,color=basin_cols[3,1])+
  annotate(geom = 'text', x = -138.5, y = 55,label='North Coast', size = 4.5,color=basin_cols[2,1])+
  annotate(geom = 'text', x = -131.5, y = 49,label='South Coast', size = 4.5,color=basin_cols[1,1])+
  annotate(geom = 'text', x = -151.5, y = 53,label='n', size = 3.5) +
  annotate(geom = 'text', x = -151.5, y = 52,label='1', size = 3.5) +
  annotate(geom = 'text', x = -151.5, y = 50.5,label='3', size = 3.5) + 
  annotate(geom = 'text', x = -151.5, y = 49,label='10', size = 3.5) + 
  annotate(geom = 'text', x = -151.5, y = 47.5,label='20', size = 3.5) + 
  annotate(geom='segment',x=-144,xend=-140.2,y=56,yend=59.5,colour=basin_cols[2],linetype=3,linewidth = 0.75)+
  annotate(geom='segment',x=-144.2,xend=-140.4,y=56.2,yend=59.7,colour=basin_cols[3],linetype=3,linewidth = 0.75)+
  annotate(geom='segment',x=-136,xend=-130.5,y=51.8,yend=51.8,colour=basin_cols[2],linetype=3,linewidth = 0.75)+
  annotate(geom='segment',x=-130.5,xend=-130.5,y=51.8,yend=53.5,colour=basin_cols[2],linetype=3,linewidth = 0.75)+
  annotate(geom='segment',x=-136,xend=-130.25,y=51.6,yend=51.6,colour=basin_cols[1],linetype=3,linewidth = 0.75)+
  annotate(geom='segment',x=-130.25,xend=-130.25,y=51.6,yend=53.5,colour=basin_cols[1],linetype=3,linewidth = 0.75)+
  coord_sf(xlim = c(-170, -120), ylim = c(47, 71), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

inset_grob <- ggplotGrob(inset)

main + 
  annotation_custom(
    grob = inset_grob,
    xmin = -141, xmax = -118, ymin = 58.5, ymax = 69
  )
#pdf version
pdf('./outputs/main figs/figure_1.pdf',width=8,height=6)
main + 
  annotation_custom(
    grob = inset_grob,
    xmin = -141, xmax = -118, ymin = 58.5, ymax = 69
  )
dev.off()


#Fig. 2 Coastwide productivity trajectories ####

#Global productivity by species
pdf('./outputs/main figs/figure_2.pdf',width=8,height=9)
par(mfrow=c(3,2))
#chinook
par(mar=c(3.1,4.1,2.1,2.1))
log.a.t.chk=exp(tv.chk$draws(variable='mu_log_a',format='draws_matrix'))
m.log.a.t.chk=apply(log.a.t.chk,2,median)
names(m.log.a.t.chk)=seq(min(chk.info$begin),max(chk.info$end))
l95.log.a.t.chk=HDInterval::hdi(log.a.t.chk,credMass=0.95)[1,]
u95.log.a.t.chk=HDInterval::hdi(log.a.t.chk,credMass=0.95)[2,]
log.a.jt.chk=tv.chk$draws(variable='log_a_t',format='draws_matrix')

L=max(chk.info$end)-min(chk.info$begin)+1

plot(m.log.a.t.chk~seq(min(chk.info$begin),max(chk.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)
polygon(
  c(seq(min(chk.info$begin),max(chk.info$end)), rev(seq(min(chk.info$begin),max(chk.info$end)))),           # x coordinates forward and back
  c(u95.log.a.t.chk, rev(l95.log.a.t.chk)),   # y upper then lower reversed
  col = adjustcolor(sp_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
for(j in 1:nrow(chk.info)){
  logaj=apply(log.a.jt.chk[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chk))],2,median)
  st=chk.info$begin[j]-min(chk.info$begin)+1
  end=chk.info$end[j]-min(chk.info$begin)+1
  lines(exp(logaj[1:c(st)])~seq(min(chk.info$begin),names(m.log.a.t.chk[c(st)])),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[st:end])~seq(chk.info$begin[j],chk.info$end[j]),col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[c(end):length(logaj)])~seq(c(chk.info$end[j]),max(chk.info$end)),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  }
for(i in 1:nrow(chk.info)){
  ts=seq(chk.info$begin[i],chk.info$end[i])
    ms=m.log.a.t.chk[names(m.log.a.t.chk)%in%ts]
    lines(ms~ts,col=adjustcolor(sp_cols[1],alpha.f=0.025),lwd=3)
}
for(i in 1:5){lines(rep(7.6,2)~c(1950,1953),col=adjustcolor(sp_cols[1],alpha.f=0.025),lwd=2)}
for(i in 1:10){lines(rep(7.2,2)~c(1950,1953),col=adjustcolor(sp_cols[1],alpha.f=0.025),lwd=2)}
for(i in 1:20){lines(rep(6.8,2)~c(1950,1953),col=adjustcolor(sp_cols[1],alpha.f=0.025),lwd=2)}
for(i in 1:50){lines(rep(6.4,2)~c(1950,1953),col=adjustcolor(sp_cols[1],alpha.f=0.025),lwd=2)}

text(x=1956,y=8,'n',cex=1)
text(x=1956,y=7.6,'5',cex=0.9)
text(x=1956,y=7.2,'10',cex=0.9)
text(x=1956,y=6.8,'20',cex=0.9)
text(x=1956,y=6.4,'50',cex=0.9)
text(x=2022,y=8,paste('N = ',nrow(chk.info)),cex=1)

mtext(expression(bold('A')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Chinook')),side=3,adj=0.5,line=0.5,cex=0.75)

par(mar=c(3.1,2.1,2.1,2.1))
log.a.t.scy=tv.scy$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.scy=apply(exp(log.a.t.scy),2,median)
names(m.log.a.t.scy)=seq(min(scy.info$begin),max(scy.info$end))
l95.log.a.t.scy=HDInterval::hdi(exp(log.a.t.scy),credMass=0.95)[1,]
u95.log.a.t.scy=HDInterval::hdi(exp(log.a.t.scy),credMass=0.95)[2,]
log.a.jt.scy=tv.scy$draws(variable='log_a_t',format='draws_matrix')

L=max(scy.info$end)-min(scy.info$begin)+1

plot(m.log.a.t.scy~seq(min(scy.info$begin),max(scy.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)
polygon(
  c(seq(min(scy.info$begin),max(scy.info$end)), rev(seq(min(scy.info$begin),max(scy.info$end)))),           # x coordinates forward and back
  c(u95.log.a.t.scy, rev(l95.log.a.t.scy)),   # y upper then lower reversed
  col = adjustcolor(sp_cols[5],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
for(j in 1:nrow(scy.info)){
  logaj=apply(log.a.jt.scy[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.scy))],2,median)
  st=scy.info$begin[j]-min(scy.info$begin)+1
  end=scy.info$end[j]-min(scy.info$begin)+1
  lines(exp(logaj[1:c(st)])~seq(min(scy.info$begin),names(m.log.a.t.scy[c(st)])),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[st:end])~seq(scy.info$begin[j],scy.info$end[j]),col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[c(end):length(logaj)])~seq(c(scy.info$end[j]),max(scy.info$end)),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
}
for(i in 1:nrow(scy.info)){
  ts=seq(scy.info$begin[i],scy.info$end[i])
  ms=m.log.a.t.scy[names(m.log.a.t.scy)%in%ts]
  lines(ms~ts,col=adjustcolor(sp_cols[5],alpha.f=0.025),lwd=3)
}
text(x=2022,y=8,paste('N = ',nrow(scy.info)),cex=1)
mtext(expression(bold('B')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Sockeye')),side=3,adj=0.5,line=0.5,cex=0.75)

#Chum
par(mar=c(3.1,4.1,2.1,2.1))
log.a.t.chm=tv.chm$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.chm=apply(exp(log.a.t.chm),2,median)
names(m.log.a.t.chm)=seq(min(chm.info$begin),max(chm.info$end))
l95.log.a.t.chm=HDInterval::hdi(exp(log.a.t.chm),credMass=0.95)[1,]
u95.log.a.t.chm=HDInterval::hdi(exp(log.a.t.chm),credMass=0.95)[2,]
log.a.jt.chm=tv.chm$draws(variable='log_a_t',format='draws_matrix')

L=max(chm.info$end)-min(chm.info$begin)+1

plot(m.log.a.t.chm~seq(min(chm.info$begin),max(chm.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)

polygon(
  c(seq(min(chm.info$begin),max(chm.info$end)), rev(seq(min(chm.info$begin),max(chm.info$end)))),           # x coordinates forward and back
  c(u95.log.a.t.chm, rev(l95.log.a.t.chm)),   # y upper then lower reversed
  col = adjustcolor(sp_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
for(j in 1:nrow(chm.info)){
  logaj=apply(log.a.jt.chm[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chm))],2,median)
  st=chm.info$begin[j]-min(chm.info$begin)+1
  end=chm.info$end[j]-min(chm.info$begin)+1
  lines(exp(logaj[1:c(st)])~seq(min(chm.info$begin),names(m.log.a.t.chm[c(st)])),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[st:end])~seq(chm.info$begin[j],chm.info$end[j]),col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[c(end):length(logaj)])~seq(c(chm.info$end[j]),max(chm.info$end)),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
}
for(i in 1:nrow(chm.info)){
  ts=seq(chm.info$begin[i],chm.info$end[i])
  ms=m.log.a.t.chm[names(m.log.a.t.chm)%in%ts]
  lines(ms~ts,col=adjustcolor(sp_cols[2],alpha.f=0.025),lwd=3)
}
text(x=2022,y=8,paste('N = ',nrow(chm.info)),cex=1)
mtext(expression(bold('Productivity (max. recruits per spawner)')),side=2,adj=0.5,line=2)
mtext(expression(bold('C')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Chum')),side=3,adj=0.5,line=0.5,cex=0.75)

#Coho
par(mar=c(3.1,2.1,2.1,2.1))
log.a.t.cho=tv.cho$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.cho=apply(exp(log.a.t.cho),2,median)
names(m.log.a.t.cho)=seq(min(cho.info$begin),max(cho.info$end))
l95.log.a.t.cho=HDInterval::hdi(exp(log.a.t.cho),credMass=0.95)[1,]
u95.log.a.t.cho=HDInterval::hdi(exp(log.a.t.cho),credMass=0.95)[2,]
log.a.jt.cho=tv.cho$draws(variable='log_a_t',format='draws_matrix')

L=max(cho.info$end)-min(cho.info$begin)+1

plot(m.log.a.t.cho~seq(min(cho.info$begin),max(cho.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)
polygon(
  c(seq(min(cho.info$begin),max(cho.info$end)), rev(seq(min(cho.info$begin),max(cho.info$end)))),           # x coordinates forward and back
  c(u95.log.a.t.cho, rev(l95.log.a.t.cho)),   # y upper then lower reversed
  col = adjustcolor(sp_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
for(j in 1:nrow(cho.info)){
  logaj=apply(log.a.jt.cho[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.cho))],2,median)
  st=cho.info$begin[j]-min(cho.info$begin)+1
  end=cho.info$end[j]-min(cho.info$begin)+1
  lines(exp(logaj[1:c(st)])~seq(min(cho.info$begin),names(m.log.a.t.cho[c(st)])),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[st:end])~seq(cho.info$begin[j],cho.info$end[j]),col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[c(end):length(logaj)])~seq(c(cho.info$end[j]),max(cho.info$end)),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
}
for(i in 1:nrow(cho.info)){
  ts=seq(cho.info$begin[i],cho.info$end[i])
  ms=m.log.a.t.cho[names(m.log.a.t.cho)%in%ts]
  lines(ms~ts,col=adjustcolor(sp_cols[3],alpha.f=0.025),lwd=3)
}
text(x=2022,y=8,paste('N = ',nrow(cho.info)),cex=1)
mtext(expression(bold('D')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Coho')),side=3,adj=0.5,line=0.5,cex=0.75)

#Pink - even
par(mar=c(5.1,4.1,2.1,2.1))
log.a.t.pke=tv.pke$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.pke=apply(exp(log.a.t.pke),2,median)
names(m.log.a.t.pke)=seq(min(pke.info$begin),max(pke.info$end))

l95.log.a.t.pke=HDInterval::hdi(exp(log.a.t.pke),credMass=0.95)[1,]
u95.log.a.t.pke=HDInterval::hdi(exp(log.a.t.pke),credMass=0.95)[2,]
log.a.jt.pke=tv.pke$draws(variable='log_a_t',format='draws_matrix')

L=max(pke.info$end)-min(pke.info$begin)+1

plot(m.log.a.t.pke~seq(min(pke.info$begin),max(pke.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)

polygon(
  c(seq(min(pke.info$begin),max(pke.info$end)), rev(seq(min(pke.info$begin),max(pke.info$end)))),           # x coordinates forward and back
  c(u95.log.a.t.pke, rev(l95.log.a.t.pke)),   # y upper then lower reversed
  col = adjustcolor(sp_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
for(j in 1:nrow(pke.info)){
  logaj=apply(log.a.jt.pke[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.pke))],2,median)
  st=pke.info$begin[j]-min(pke.info$begin)+1
  end=pke.info$end[j]-min(pke.info$begin)+1
  lines(exp(logaj[1:c(st)])~seq(min(pke.info$begin),names(m.log.a.t.pke[c(st)])),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[st:end])~seq(pke.info$begin[j],pke.info$end[j]),col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[c(end):length(logaj)])~seq(c(pke.info$end[j]),max(pke.info$end)),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
}
for(i in 1:nrow(pke.info)){
  ts=seq(pke.info$begin[i],pke.info$end[i])
  ms=m.log.a.t.pke[names(m.log.a.t.pke)%in%ts]
  lines(ms~ts,col=adjustcolor(sp_cols[4],alpha.f=0.025),lwd=3)
}
text(x=2022,y=8,paste('N = ',nrow(pke.info)),cex=1)
mtext(expression(bold('E')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Pink - Even')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('Brood cohort year')),side=1,at=2035,line=3)

#Pink - odd
par(mar=c(5.1,2.1,2.1,2.1))
log.a.t.pko=tv.pko$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.pko=apply(exp(log.a.t.pko),2,median)
names(m.log.a.t.pko)=seq(min(pko.info$begin),max(pko.info$end))

l95.log.a.t.pko=HDInterval::hdi(exp(log.a.t.pko),credMass=0.95)[1,]
u95.log.a.t.pko=HDInterval::hdi(exp(log.a.t.pko),credMass=0.95)[2,]
log.a.jt.pko=tv.pko$draws(variable='log_a_t',format='draws_matrix')

L=max(pko.info$end)-min(pko.info$begin)+1

plot(m.log.a.t.pko~seq(min(pko.info$begin),max(pko.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)
polygon(
  c(seq(min(pko.info$begin),max(pko.info$end)), rev(seq(min(pko.info$begin),max(pko.info$end)))),           # x coordinates forward and back
  c(u95.log.a.t.pko, rev(l95.log.a.t.pko)),   # y upper then lower reversed
  col = adjustcolor(sp_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
for(j in 1:nrow(pko.info)){
  logaj=apply(log.a.jt.pko[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.pko))],2,median)
  st=pko.info$begin[j]-min(pko.info$begin)+1
  end=pko.info$end[j]-min(pko.info$begin)+1
  lines(exp(logaj[1:c(st)])~seq(min(pko.info$begin),names(m.log.a.t.pko[c(st)])),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[st:end])~seq(pko.info$begin[j],pko.info$end[j]),col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[c(end):length(logaj)])~seq(c(pko.info$end[j]),max(pko.info$end)),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
}
for(i in 1:nrow(pko.info)){
  ts=seq(pko.info$begin[i],pko.info$end[i])
  ms=m.log.a.t.pko[names(m.log.a.t.pko)%in%ts]
  lines(ms~ts,col=adjustcolor(sp_cols[4],alpha.f=0.025),lwd=3)
}
text(x=2022,y=8,paste('N = ',nrow(pko.info)),cex=1)
mtext(expression(bold('F')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Pink - Odd')),side=3,adj=0.5,line=0.5,cex=0.75)
dev.off()



#Fig. 3 Productivity trajectories by region #### 
pdf('./outputs/main figs/figure_3.pdf',width=8,height=9)
par(mfrow=c(3,2))
par(mar=c(3.1,4.1,2.1,2.1))
#Mean productivity trends by ocean basin
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

plot(apply(exp(mu_loga.chk.bs),2,median)~seq(min(chk.info$begin),max(chk.info$end)),ylim=c(0,8),ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)

l95.chk.bs=HDInterval::hdi(exp(mu_loga.chk.bs),credMass=0.9)[1,]
u95.chk.bs=HDInterval::hdi(exp(mu_loga.chk.bs),credMass=0.9)[2,]
l95.chk.goa=HDInterval::hdi(exp(mu_loga.chk.goa),credMass=0.9)[1,]
u95.chk.goa=HDInterval::hdi(exp(mu_loga.chk.goa),credMass=0.9)[2,]
l95.chk.nc=HDInterval::hdi(exp(mu_loga.chk.nc),credMass=0.9)[1,]
u95.chk.nc=HDInterval::hdi(exp(mu_loga.chk.nc),credMass=0.9)[2,]
l95.chk.sc=HDInterval::hdi(exp(mu_loga.chk.sc),credMass=0.9)[1,]
u95.chk.sc=HDInterval::hdi(exp(mu_loga.chk.sc),credMass=0.9)[2,]
polygon(
  c(seq(min(chk.info$begin),max(chk.info$end)), rev(seq(min(chk.info$begin),max(chk.info$end)))),           # x coordinates forward and back
  c(u95.chk.bs, rev(l95.chk.bs)),   # y upper then lower reversed
  col = adjustcolor(basin_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(chk.info$begin),max(chk.info$end)), rev(seq(min(chk.info$begin),max(chk.info$end)))),           # x coordinates forward and back
  c(u95.chk.goa, rev(l95.chk.goa)), # y upper then lower reversed
  col = adjustcolor(basin_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(chk.info$begin),max(chk.info$end)), rev(seq(min(chk.info$begin),max(chk.info$end)))),           # x coordinates forward and back
  c(u95.chk.nc, rev(l95.chk.nc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(chk.info$begin),max(chk.info$end)), rev(seq(min(chk.info$begin),max(chk.info$end)))),           # x coordinates forward and back
  c(u95.chk.sc, rev(l95.chk.sc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
mu_loga.chk.sc.rs=apply(exp(mu_loga.chk.sc),2,median)
names(mu_loga.chk.sc.rs)=seq(min(chk.info$begin),max(chk.info$end))
mu_loga.chk.nc.rs=apply(exp(mu_loga.chk.nc),2,median)
names(mu_loga.chk.nc.rs)=seq(min(chk.info$begin),max(chk.info$end))
mu_loga.chk.goa.rs=apply(exp(mu_loga.chk.goa),2,median)
names(mu_loga.chk.goa.rs)=seq(min(chk.info$begin),max(chk.info$end))
mu_loga.chk.bs.rs=apply(exp(mu_loga.chk.bs),2,median)
names(mu_loga.chk.bs.rs)=seq(min(chk.info$begin),max(chk.info$end))

for(i in 1:nrow(chk.info)){
  ts=seq(chk.info$begin[i],chk.info$end[i])
  if(chk.info$ocean.basin[i]=='SC'){
    ms=mu_loga.chk.sc.rs[names(mu_loga.chk.sc.rs)%in%ts]
    lines(ms~ts,col=adjustcolor(basin_cols[1],alpha.f=0.1),lwd=3)
  }else
    if(chk.info$ocean.basin[i]=='NC'){
      ms=mu_loga.chk.nc.rs[names(mu_loga.chk.nc.rs)%in%ts]
      lines(ms~ts,col=adjustcolor(basin_cols[2],alpha.f=0.1),lwd=3)
    }else
      if(chk.info$ocean.basin[i]=='GoA'){
        ms=mu_loga.chk.goa.rs[names(mu_loga.chk.goa.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[3],alpha.f=0.1),lwd=3)
      }else{
        ms=mu_loga.chk.bs.rs[names(mu_loga.chk.bs.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=3)
      }
}
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.25,'N',cex=1)
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.8,nrow(subset(chk.info,ocean.basin=='BS')),cex=1,col=basin_cols[4])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-1.4,nrow(subset(chk.info,ocean.basin=='GoA')),cex=1,col=basin_cols[3])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.0,nrow(subset(chk.info,ocean.basin=='NC')),cex=1,col=basin_cols[2])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.6,nrow(subset(chk.info,ocean.basin=='SC')),cex=1,col=basin_cols[1])

text(x=par('usr')[2]-6,y=par('usr')[4]-0.8,expression(bold('Bering Sea')),col=basin_cols[4],cex=1.25,adj=1)
text(x=par('usr')[2]-6,y=par('usr')[4]-1.4,expression(bold('Gulf of Alaska')),col=basin_cols[3],cex=1.25,adj=1)
text(x=par('usr')[2]-6,y=par('usr')[4]-2.0,expression(bold('North Coast')),col=basin_cols[2],cex=1.25,adj=1)
text(x=par('usr')[2]-6,y=par('usr')[4]-2.6,expression(bold('South Coast')),col=basin_cols[1],cex=1.25,adj=1)
lines(rep(7.6,2)~c(1950,1953),col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=2)
for(i in 1:5){lines(rep(7.2,2)~c(1950,1953),col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=2)}
for(i in 1:10){lines(rep(6.8,2)~c(1950,1953),col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=2)}
for(i in 1:20){lines(rep(6.4,2)~c(1950,1953),col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=2)}
text(x=1956,y=8,'n',cex=1)
text(x=1956,y=7.6,'1',cex=0.9)
text(x=1956,y=7.2,'5',cex=0.9)
text(x=1956,y=6.8,'10',cex=0.9)
text(x=1956,y=6.4,'20',cex=0.9)

mtext(expression(bold('A')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Chinook')),side=3,adj=0.5,line=0.5,cex=0.75)

par(mar=c(3.1,2.1,2.1,2.1))
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

plot(apply(exp(mu_loga.scy.bs),2,median)~seq(min(scy.info$begin),max(scy.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)

l95.scy.bs=HDInterval::hdi(exp(mu_loga.scy.bs),credMass=0.9)[1,]
u95.scy.bs=HDInterval::hdi(exp(mu_loga.scy.bs),credMass=0.9)[2,]
l95.scy.goa=HDInterval::hdi(exp(mu_loga.scy.goa),credMass=0.9)[1,]
u95.scy.goa=HDInterval::hdi(exp(mu_loga.scy.goa),credMass=0.9)[2,]
l95.scy.nc=HDInterval::hdi(exp(mu_loga.scy.nc),credMass=0.9)[1,]
u95.scy.nc=HDInterval::hdi(exp(mu_loga.scy.nc),credMass=0.9)[2,]
l95.scy.sc=HDInterval::hdi(exp(mu_loga.scy.sc),credMass=0.9)[1,]
u95.scy.sc=HDInterval::hdi(exp(mu_loga.scy.sc),credMass=0.9)[2,]
polygon(
  c(seq(min(scy.info$begin),max(scy.info$end)), rev(seq(min(scy.info$begin),max(scy.info$end)))),           # x coordinates forward and back
  c(u95.scy.bs, rev(l95.scy.bs)),   # y upper then lower reversed
  col = adjustcolor(basin_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(scy.info$begin),max(scy.info$end)), rev(seq(min(scy.info$begin),max(scy.info$end)))),           # x coordinates forward and back
  c(u95.scy.goa, rev(l95.scy.goa)), # y upper then lower reversed
  col = adjustcolor(basin_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(scy.info$begin),max(scy.info$end)), rev(seq(min(scy.info$begin),max(scy.info$end)))),           # x coordinates forward and back
  c(u95.scy.nc, rev(l95.scy.nc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(scy.info$begin),max(scy.info$end)), rev(seq(min(scy.info$begin),max(scy.info$end)))),           # x coordinates forward and back
  c(u95.scy.sc, rev(l95.scy.sc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
mu_loga.scy.sc.rs=apply(exp(mu_loga.scy.sc),2,median)
names(mu_loga.scy.sc.rs)=seq(min(scy.info$begin),max(scy.info$end))
mu_loga.scy.nc.rs=apply(exp(mu_loga.scy.nc),2,median)
names(mu_loga.scy.nc.rs)=seq(min(scy.info$begin),max(scy.info$end))
mu_loga.scy.goa.rs=apply(exp(mu_loga.scy.goa),2,median)
names(mu_loga.scy.goa.rs)=seq(min(scy.info$begin),max(scy.info$end))
mu_loga.scy.bs.rs=apply(exp(mu_loga.scy.bs),2,median)
names(mu_loga.scy.bs.rs)=seq(min(scy.info$begin),max(scy.info$end))

for(i in 1:nrow(scy.info)){
  ts=seq(scy.info$begin[i],scy.info$end[i])
  if(scy.info$ocean.basin[i]=='SC'){
    ms=mu_loga.scy.sc.rs[names(mu_loga.scy.sc.rs)%in%ts]
    lines(ms~ts,col=adjustcolor(basin_cols[1],alpha.f=0.1),lwd=3)
  }else
    if(scy.info$ocean.basin[i]=='NC'){
      ms=mu_loga.scy.nc.rs[names(mu_loga.scy.nc.rs)%in%ts]
      lines(ms~ts,col=adjustcolor(basin_cols[2],alpha.f=0.1),lwd=3)
    }else
      if(scy.info$ocean.basin[i]=='GoA'){
        ms=mu_loga.scy.goa.rs[names(mu_loga.scy.goa.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[3],alpha.f=0.1),lwd=3)
      }else{
        ms=mu_loga.scy.bs.rs[names(mu_loga.scy.bs.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=3)
      }
}
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.25,'N',cex=1)
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.8,nrow(subset(scy.info,ocean.basin=='BS')),cex=1,col=basin_cols[4])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-1.4,nrow(subset(scy.info,ocean.basin=='GoA')),cex=1,col=basin_cols[3])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.0,nrow(subset(scy.info,ocean.basin=='NC')),cex=1,col=basin_cols[2])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.6,nrow(subset(scy.info,ocean.basin=='SC')),cex=1,col=basin_cols[1])

mtext(expression(bold('B')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Sockeye')),side=3,adj=0.5,line=0.5,cex=0.75)

par(mar=c(3.1,4.1,2.1,2.1))
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

plot(apply(exp(mu_loga.chm.bs),2,median)~seq(min(chm.info$begin),max(chm.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)

l95.chm.bs=HDInterval::hdi(exp(mu_loga.chm.bs),credMass=0.9)[1,]
u95.chm.bs=HDInterval::hdi(exp(mu_loga.chm.bs),credMass=0.9)[2,]
l95.chm.goa=HDInterval::hdi(exp(mu_loga.chm.goa),credMass=0.9)[1,]
u95.chm.goa=HDInterval::hdi(exp(mu_loga.chm.goa),credMass=0.9)[2,]
l95.chm.nc=HDInterval::hdi(exp(mu_loga.chm.nc),credMass=0.9)[1,]
u95.chm.nc=HDInterval::hdi(exp(mu_loga.chm.nc),credMass=0.9)[2,]
l95.chm.sc=HDInterval::hdi(exp(mu_loga.chm.sc),credMass=0.9)[1,]
u95.chm.sc=HDInterval::hdi(exp(mu_loga.chm.sc),credMass=0.9)[2,]
polygon(
  c(seq(min(chm.info$begin),max(chm.info$end)), rev(seq(min(chm.info$begin),max(chm.info$end)))),           # x coordinates forward and back
  c(u95.chm.bs, rev(l95.chm.bs)),   # y upper then lower reversed
  col = adjustcolor(basin_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(chm.info$begin),max(chm.info$end)), rev(seq(min(chm.info$begin),max(chm.info$end)))),           # x coordinates forward and back
  c(u95.chm.goa, rev(l95.chm.goa)), # y upper then lower reversed
  col = adjustcolor(basin_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(chm.info$begin),max(chm.info$end)), rev(seq(min(chm.info$begin),max(chm.info$end)))),           # x coordinates forward and back
  c(u95.chm.nc, rev(l95.chm.nc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(chm.info$begin),max(chm.info$end)), rev(seq(min(chm.info$begin),max(chm.info$end)))),           # x coordinates forward and back
  c(u95.chm.sc, rev(l95.chm.sc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
mu_loga.chm.sc.rs=apply(exp(mu_loga.chm.sc),2,median)
names(mu_loga.chm.sc.rs)=seq(min(chm.info$begin),max(chm.info$end))
mu_loga.chm.nc.rs=apply(exp(mu_loga.chm.nc),2,median)
names(mu_loga.chm.nc.rs)=seq(min(chm.info$begin),max(chm.info$end))
mu_loga.chm.goa.rs=apply(exp(mu_loga.chm.goa),2,median)
names(mu_loga.chm.goa.rs)=seq(min(chm.info$begin),max(chm.info$end))
mu_loga.chm.bs.rs=apply(exp(mu_loga.chm.bs),2,median)
names(mu_loga.chm.bs.rs)=seq(min(chm.info$begin),max(chm.info$end))

for(i in 1:nrow(chm.info)){
  ts=seq(chm.info$begin[i],chm.info$end[i])
  if(chm.info$ocean.basin[i]=='SC'){
    ms=mu_loga.chm.sc.rs[names(mu_loga.chm.sc.rs)%in%ts]
    lines(ms~ts,col=adjustcolor(basin_cols[1],alpha.f=0.1),lwd=3)
  }else
    if(chm.info$ocean.basin[i]=='NC'){
      ms=mu_loga.chm.nc.rs[names(mu_loga.chm.nc.rs)%in%ts]
      lines(ms~ts,col=adjustcolor(basin_cols[2],alpha.f=0.1),lwd=3)
    }else
      if(chm.info$ocean.basin[i]=='GoA'){
        ms=mu_loga.chm.goa.rs[names(mu_loga.chm.goa.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[3],alpha.f=0.1),lwd=3)
      }else{
        ms=mu_loga.chm.bs.rs[names(mu_loga.chm.bs.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=3)
      }
}
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.25,'N',cex=1)
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.8,nrow(subset(chm.info,ocean.basin=='BS')),cex=1,col=basin_cols[4])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-1.4,nrow(subset(chm.info,ocean.basin=='GoA')),cex=1,col=basin_cols[3])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.0,nrow(subset(chm.info,ocean.basin=='NC')),cex=1,col=basin_cols[2])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.6,nrow(subset(chm.info,ocean.basin=='SC')),cex=1,col=basin_cols[1])

mtext(expression(bold('Mean productivity (max. recruits per spawner)')),side=2,adj=0.5,line=2)
mtext(expression(bold('C')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Chum')),side=3,adj=0.5,line=0.5,cex=0.75)

par(mar=c(3.1,2.1,2.1,2.1))
log.a.jt.cho=tv.cho$draws(variable='log_a_t',format='draws_matrix')
cho.sc=which(cho.info$ocean.basin=='SC')
cho.nc=which(cho.info$ocean.basin=='NC')
cho.goa=which(cho.info$ocean.basin=='GoA')

logat=log.a.jt.cho
colnames(logat)=gsub("]", ".", colnames(logat))
colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)

logat.sc=logat[,grepl(paste(paste(',',cho.sc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.nc=logat[,grepl(paste(paste(',',cho.nc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.goa=logat[,grepl(paste(paste(',',cho.goa,'.',sep=''),collapse = "|"),colnames(logat))]

L=max(cho.info$end)-min(cho.info$begin)+1
mu_loga.cho.sc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.cho.nc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.cho.goa=matrix(nrow=nrow(logat),ncol=L)
for(t in 1:L){
  logas.sc=logat.sc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.sc))]
  mu_loga.cho.sc[,t]=apply(logas.sc,1,mean)
  
  logas.nc=logat.nc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.nc))]
  mu_loga.cho.nc[,t]=apply(logas.nc,1,mean)
  
  logas.goa=logat.goa[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.goa))]
  mu_loga.cho.goa[,t]=apply(logas.goa,1,mean)
  
}

plot(apply(exp(mu_loga.cho.sc),2,median)~seq(min(cho.info$begin),max(cho.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)

l95.cho.goa=HDInterval::hdi(exp(mu_loga.cho.goa),credMass=0.9)[1,]
u95.cho.goa=HDInterval::hdi(exp(mu_loga.cho.goa),credMass=0.9)[2,]
l95.cho.nc=HDInterval::hdi(exp(mu_loga.cho.nc),credMass=0.9)[1,]
u95.cho.nc=HDInterval::hdi(exp(mu_loga.cho.nc),credMass=0.9)[2,]
l95.cho.sc=HDInterval::hdi(exp(mu_loga.cho.sc),credMass=0.9)[1,]
u95.cho.sc=HDInterval::hdi(exp(mu_loga.cho.sc),credMass=0.9)[2,]
polygon(
  c(seq(min(cho.info$begin),max(cho.info$end)), rev(seq(min(cho.info$begin),max(cho.info$end)))),           # x coordinates forward and back
  c(u95.cho.goa, rev(l95.cho.goa)), # y upper then lower reversed
  col = adjustcolor(basin_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(cho.info$begin),max(cho.info$end)), rev(seq(min(cho.info$begin),max(cho.info$end)))),           # x coordinates forward and back
  c(u95.cho.nc, rev(l95.cho.nc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(cho.info$begin),max(cho.info$end)), rev(seq(min(cho.info$begin),max(cho.info$end)))),           # x coordinates forward and back
  c(u95.cho.sc, rev(l95.cho.sc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
mu_loga.cho.sc.rs=apply(exp(mu_loga.cho.sc),2,median)
names(mu_loga.cho.sc.rs)=seq(min(cho.info$begin),max(cho.info$end))
mu_loga.cho.nc.rs=apply(exp(mu_loga.cho.nc),2,median)
names(mu_loga.cho.nc.rs)=seq(min(cho.info$begin),max(cho.info$end))
mu_loga.cho.goa.rs=apply(exp(mu_loga.cho.goa),2,median)
names(mu_loga.cho.goa.rs)=seq(min(cho.info$begin),max(cho.info$end))

for(i in 1:nrow(cho.info)){
  ts=seq(cho.info$begin[i],cho.info$end[i])
  if(cho.info$ocean.basin[i]=='SC'){
    ms=mu_loga.cho.sc.rs[names(mu_loga.cho.sc.rs)%in%ts]
    lines(ms~ts,col=adjustcolor(basin_cols[1],alpha.f=0.1),lwd=3)
  }else
    if(cho.info$ocean.basin[i]=='NC'){
      ms=mu_loga.cho.nc.rs[names(mu_loga.cho.nc.rs)%in%ts]
      lines(ms~ts,col=adjustcolor(basin_cols[2],alpha.f=0.1),lwd=3)
    }else
      if(cho.info$ocean.basin[i]=='GoA'){
        ms=mu_loga.cho.goa.rs[names(mu_loga.cho.goa.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[3],alpha.f=0.1),lwd=3)
      }
}
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.25,'N',cex=1)
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.8,nrow(subset(cho.info,ocean.basin=='GoA')),cex=1,col=basin_cols[3])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-1.4,nrow(subset(cho.info,ocean.basin=='NC')),cex=1,col=basin_cols[2])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.0,nrow(subset(cho.info,ocean.basin=='SC')),cex=1,col=basin_cols[1])

mtext(expression(bold('D')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Coho')),side=3,adj=0.5,line=0.5,cex=0.75)

par(mar=c(5.1,4.1,2.1,2.1))
log.a.jt.pke=tv.pke$draws(variable='log_a_t',format='draws_matrix')

logat=log.a.jt.pke
colnames(logat)=gsub("]", ".", colnames(logat))
colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)

logat.sc=logat[,grepl(paste(paste(',',pke.sc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.nc=logat[,grepl(paste(paste(',',pke.nc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.goa=logat[,grepl(paste(paste(',',pke.goa,'.',sep=''),collapse = "|"),colnames(logat))]
logat.bs=logat[,grepl(paste(paste(',',pke.bs,'.',sep=''),collapse = "|"),colnames(logat))]

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

plot(apply(exp(mu_loga.pke.bs),2,median)~seq(min(pke.info$begin),max(pke.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)
l95.pke.bs=HDInterval::hdi(exp(mu_loga.pke.bs),credMass=0.9)[1,]
u95.pke.bs=HDInterval::hdi(exp(mu_loga.pke.bs),credMass=0.9)[2,]
l95.pke.goa=HDInterval::hdi(exp(mu_loga.pke.goa),credMass=0.9)[1,]
u95.pke.goa=HDInterval::hdi(exp(mu_loga.pke.goa),credMass=0.9)[2,]
l95.pke.nc=HDInterval::hdi(exp(mu_loga.pke.nc),credMass=0.9)[1,]
u95.pke.nc=HDInterval::hdi(exp(mu_loga.pke.nc),credMass=0.9)[2,]
l95.pke.sc=HDInterval::hdi(exp(mu_loga.pke.sc),credMass=0.9)[1,]
u95.pke.sc=HDInterval::hdi(exp(mu_loga.pke.sc),credMass=0.9)[2,]
polygon(
  c(seq(min(pke.info$begin),max(pke.info$end)), rev(seq(min(pke.info$begin),max(pke.info$end)))),           # x coordinates forward and back
  c(u95.pke.bs, rev(l95.pke.bs)),   # y upper then lower reversed
  col = adjustcolor(basin_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(pke.info$begin),max(pke.info$end)), rev(seq(min(pke.info$begin),max(pke.info$end)))),           # x coordinates forward and back
  c(u95.pke.goa, rev(l95.pke.goa)), # y upper then lower reversed
  col = adjustcolor(basin_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(pke.info$begin),max(pke.info$end)), rev(seq(min(pke.info$begin),max(pke.info$end)))),           # x coordinates forward and back
  c(u95.pke.nc, rev(l95.pke.nc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(pke.info$begin),max(pke.info$end)), rev(seq(min(pke.info$begin),max(pke.info$end)))),           # x coordinates forward and back
  c(u95.pke.sc, rev(l95.pke.sc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
mu_loga.pke.sc.rs=apply(exp(mu_loga.pke.sc),2,median)
names(mu_loga.pke.sc.rs)=seq(min(pke.info$begin),max(pke.info$end))
mu_loga.pke.nc.rs=apply(exp(mu_loga.pke.nc),2,median)
names(mu_loga.pke.nc.rs)=seq(min(pke.info$begin),max(pke.info$end))
mu_loga.pke.goa.rs=apply(exp(mu_loga.pke.goa),2,median)
names(mu_loga.pke.goa.rs)=seq(min(pke.info$begin),max(pke.info$end))
mu_loga.pke.bs.rs=apply(exp(mu_loga.pke.bs),2,median)
names(mu_loga.pke.bs.rs)=seq(min(pke.info$begin),max(pke.info$end))

for(i in 1:nrow(pke.info)){
  ts=seq(pke.info$begin[i],pke.info$end[i])
  if(pke.info$ocean.basin[i]=='SC'){
    ms=mu_loga.pke.sc.rs[names(mu_loga.pke.sc.rs)%in%ts]
    lines(ms~ts,col=adjustcolor(basin_cols[1],alpha.f=0.1),lwd=3)
  }else
    if(pke.info$ocean.basin[i]=='NC'){
      ms=mu_loga.pke.nc.rs[names(mu_loga.pke.nc.rs)%in%ts]
      lines(ms~ts,col=adjustcolor(basin_cols[2],alpha.f=0.1),lwd=3)
    }else
      if(pke.info$ocean.basin[i]=='GoA'){
        ms=mu_loga.pke.goa.rs[names(mu_loga.pke.goa.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[3],alpha.f=0.1),lwd=3)
      }else{
        ms=mu_loga.pke.bs.rs[names(mu_loga.pke.bs.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=3)
      }
}
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.25,'N',cex=1)
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.8,nrow(subset(pke.info,ocean.basin=='BS')),cex=1,col=basin_cols[4])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-1.4,nrow(subset(pke.info,ocean.basin=='GoA')),cex=1,col=basin_cols[3])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.0,nrow(subset(pke.info,ocean.basin=='NC')),cex=1,col=basin_cols[2])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.6,nrow(subset(pke.info,ocean.basin=='SC')),cex=1,col=basin_cols[1])

mtext(expression(bold('E')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Pink - Even')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('Brood cohort year')),side=1,at=2035,line=3)

par(mar=c(5.1,2.1,2.1,2.1))
log.a.jt.pko=tv.pko$draws(variable='log_a_t',format='draws_matrix')

logat=log.a.jt.pko
colnames(logat)=gsub("]", ".", colnames(logat))
colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)

logat.sc=logat[,grepl(paste(paste(',',pko.sc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.nc=logat[,grepl(paste(paste(',',pko.nc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.goa=logat[,grepl(paste(paste(',',pko.goa,'.',sep=''),collapse = "|"),colnames(logat))]
logat.bs=logat[,grepl(paste(paste(',',pko.bs,'.',sep=''),collapse = "|"),colnames(logat))]


L=max(pko.info$end)-min(pko.info$begin)+1
mu_loga.pko.sc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.pko.nc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.pko.goa=matrix(nrow=nrow(logat),ncol=L)
mu_loga.pko.bs=matrix(nrow=nrow(logat),ncol=L)
for(t in 1:L){
  logas.sc=logat.sc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.sc))]
  mu_loga.pko.sc[,t]=apply(logas.sc,1,mean)
  
  logas.nc=logat.nc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.nc))]
  mu_loga.pko.nc[,t]=apply(logas.nc,1,mean)
  
  logas.goa=logat.goa[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.goa))]
  mu_loga.pko.goa[,t]=apply(logas.goa,1,mean)
  
  logas.bs=logat.bs[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.bs))]
  mu_loga.pko.bs[,t]=apply(logas.bs,1,mean)
}

plot(apply(exp(mu_loga.pko.bs),2,median)~seq(min(pko.info$begin),max(pko.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)
l95.pko.bs=HDInterval::hdi(exp(mu_loga.pko.bs),credMass=0.9)[1,]
u95.pko.bs=HDInterval::hdi(exp(mu_loga.pko.bs),credMass=0.9)[2,]
l95.pko.goa=HDInterval::hdi(exp(mu_loga.pko.goa),credMass=0.9)[1,]
u95.pko.goa=HDInterval::hdi(exp(mu_loga.pko.goa),credMass=0.9)[2,]
l95.pko.nc=HDInterval::hdi(exp(mu_loga.pko.nc),credMass=0.9)[1,]
u95.pko.nc=HDInterval::hdi(exp(mu_loga.pko.nc),credMass=0.9)[2,]
l95.pko.sc=HDInterval::hdi(exp(mu_loga.pko.sc),credMass=0.9)[1,]
u95.pko.sc=HDInterval::hdi(exp(mu_loga.pko.sc),credMass=0.9)[2,]
polygon(
  c(seq(min(pko.info$begin),max(pko.info$end)), rev(seq(min(pko.info$begin),max(pko.info$end)))),           # x coordinates forward and back
  c(u95.pko.bs, rev(l95.pko.bs)),   # y upper then lower reversed
  col = adjustcolor(basin_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(pko.info$begin),max(pko.info$end)), rev(seq(min(pko.info$begin),max(pko.info$end)))),           # x coordinates forward and back
  c(u95.pko.goa, rev(l95.pko.goa)), # y upper then lower reversed
  col = adjustcolor(basin_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(pko.info$begin),max(pko.info$end)), rev(seq(min(pko.info$begin),max(pko.info$end)))),           # x coordinates forward and back
  c(u95.pko.nc, rev(l95.pko.nc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(pko.info$begin),max(pko.info$end)), rev(seq(min(pko.info$begin),max(pko.info$end)))),           # x coordinates forward and back
  c(u95.pko.sc, rev(l95.pko.sc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
mu_loga.pko.sc.rs=apply(exp(mu_loga.pko.sc),2,median)
names(mu_loga.pko.sc.rs)=seq(min(pko.info$begin),max(pko.info$end))
mu_loga.pko.nc.rs=apply(exp(mu_loga.pko.nc),2,median)
names(mu_loga.pko.nc.rs)=seq(min(pko.info$begin),max(pko.info$end))
mu_loga.pko.goa.rs=apply(exp(mu_loga.pko.goa),2,median)
names(mu_loga.pko.goa.rs)=seq(min(pko.info$begin),max(pko.info$end))
mu_loga.pko.bs.rs=apply(exp(mu_loga.pko.bs),2,median)
names(mu_loga.pko.bs.rs)=seq(min(pko.info$begin),max(pko.info$end))

for(i in 1:nrow(pko.info)){
  ts=seq(pko.info$begin[i],pko.info$end[i])
  if(pko.info$ocean.basin[i]=='SC'){
    ms=mu_loga.pko.sc.rs[names(mu_loga.pko.sc.rs)%in%ts]
    lines(ms~ts,col=adjustcolor(basin_cols[1],alpha.f=0.1),lwd=3)
  }else
    if(pko.info$ocean.basin[i]=='NC'){
      ms=mu_loga.pko.nc.rs[names(mu_loga.pko.nc.rs)%in%ts]
      lines(ms~ts,col=adjustcolor(basin_cols[2],alpha.f=0.1),lwd=3)
    }else
      if(pko.info$ocean.basin[i]=='GoA'){
        ms=mu_loga.pko.goa.rs[names(mu_loga.pko.goa.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[3],alpha.f=0.1),lwd=3)
      }else{
        ms=mu_loga.pko.bs.rs[names(mu_loga.pko.bs.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=3)
      }
}
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.25,'N',cex=1)
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.8,nrow(subset(pko.info,ocean.basin=='BS')),cex=1,col=basin_cols[4])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-1.4,nrow(subset(pko.info,ocean.basin=='GoA')),cex=1,col=basin_cols[3])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.0,nrow(subset(pko.info,ocean.basin=='NC')),cex=1,col=basin_cols[2])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.6,nrow(subset(pko.info,ocean.basin=='SC')),cex=1,col=basin_cols[1])

mtext(expression(bold('F')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Pink - Odd')),side=3,adj=0.5,line=0.5,cex=0.75)
dev.off()

#Fig. 4 Recent vs stationary Umsy ####

#Chinook
loga.t.chk=tv.chk$draws(variables='log_a_t',format='draws_matrix')
umsy.st.chk=stc.chk$draws(variables='Umsy',format='draws_matrix')

trend_chk_l=matrix(nrow=nrow(umsy.st.chk),ncol=ncol(umsy.st.chk))
umsy_comp_chk=data.frame(stock=chk.info$stock.name,st.umsy=NA,st.umsy.l90=NA,st.umsy.u90=NA,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,p.change=NA)
for(i in 1:nrow(chk.info)){
  logas=loga.t.chk[,grepl(paste(',',i,']',sep=''),colnames(loga.t.chk))]
  umsy.st=umsy.st.chk[,i]
  umsy.st=ifelse(umsy.st<0,0.005,umsy.st)
  
  o.logas=apply(logas[,c(chk.info$begin[i]-min(chk.info$begin)+1):c(chk.info$begin[i]-min(chk.info$begin)+5)],1,mean)
  r.logas=apply(logas[,c(chk.info$end[i]-min(chk.info$begin)+1-5):c(chk.info$end[i]-min(chk.info$begin))],1,mean)
  r.umsy=sapply(r.logas,umsyCalc)
  r.umsy=ifelse(r.umsy<0,0.005,r.umsy)
   
  trend_chk_l[,i]=(r.logas-o.logas)/(chk.info$end[i]-chk.info$begin[i]+1-5)
  umsy_comp_chk[i,2]=median(umsy.st)
  umsy_comp_chk[i,3]=HDInterval::hdi(umsy.st,credMass=0.9)[1]
  umsy_comp_chk[i,4]=HDInterval::hdi(umsy.st,credMass=0.9)[2]
  umsy_comp_chk[i,5]=median(r.umsy)
  umsy_comp_chk[i,6]=HDInterval::hdi(r.umsy,credMass=0.9)[1]
  umsy_comp_chk[i,7]=HDInterval::hdi(r.umsy,credMass=0.9)[2]
  umsy_comp_chk[i,8]=median((r.umsy/umsy.st)-1)
  }

trend_chk_l.sc=trend_chk_l[,chk.sc]
trend_chk_l.nc=trend_chk_l[,chk.nc]
trend_chk_l.goa=trend_chk_l[,chk.goa]
trend_chk_l.bs=trend_chk_l[,chk.bs]


#Sockeye
loga.t.scy=tv.scy$draws(variables='log_a_t',format='draws_matrix')
umsy.st.scy=stc.scy$draws(variables='Umsy',format='draws_matrix')

trend_scy_l=matrix(nrow=nrow(umsy.st.scy),ncol=ncol(umsy.st.scy))
umsy_comp_scy=data.frame(stock=scy.info$stock.name,st.umsy=NA,st.umsy.l90=NA,st.umsy.u90=NA,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,p.change=NA)
for(i in 1:nrow(scy.info)){
  logas=loga.t.scy[,grepl(paste(',',i,']',sep=''),colnames(loga.t.scy))]
  umsy.st=umsy.st.scy[,i]
  umsy.st=ifelse(umsy.st<0,0.005,umsy.st)
  
  o.logas=apply(logas[,c(scy.info$begin[i]-min(scy.info$begin)+1):c(scy.info$begin[i]-min(scy.info$begin)+5)],1,mean)
  r.logas=apply(logas[,c(scy.info$end[i]-min(scy.info$begin)+1-5):c(scy.info$end[i]-min(scy.info$begin))],1,mean)
  r.umsy=sapply(r.logas,umsyCalc)
  r.umsy=ifelse(r.umsy<0,0.005,r.umsy)
  
  trend_scy_l[,i]=(r.logas-o.logas)/(scy.info$end[i]-scy.info$begin[i]+1-5)
  umsy_comp_scy[i,2]=median(umsy.st)
  umsy_comp_scy[i,3]=HDInterval::hdi(umsy.st,credMass=0.9)[1]
  umsy_comp_scy[i,4]=HDInterval::hdi(umsy.st,credMass=0.9)[2]
  umsy_comp_scy[i,5]=median(r.umsy)
  umsy_comp_scy[i,6]=HDInterval::hdi(r.umsy,credMass=0.9)[1]
  umsy_comp_scy[i,7]=HDInterval::hdi(r.umsy,credMass=0.9)[2]
  umsy_comp_scy[i,8]=median((r.umsy/umsy.st)-1)
}

trend_scy_l.sc=trend_scy_l[,scy.sc]
trend_scy_l.nc=trend_scy_l[,scy.nc]
trend_scy_l.goa=trend_scy_l[,scy.goa]
trend_scy_l.bs=trend_scy_l[,scy.bs]

#Chum
loga.t.chm=tv.chm$draws(variables='log_a_t',format='draws_matrix')
umsy.st.chm=stc.chm$draws(variables='Umsy',format='draws_matrix')

trend_chm_l=matrix(nrow=nrow(umsy.st.chm),ncol=ncol(umsy.st.chm))
umsy_comp_chm=data.frame(stock=chm.info$stock.name,st.umsy=NA,st.umsy.l90=NA,st.umsy.u90=NA,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,p.change=NA)
for(i in 1:nrow(chm.info)){
  logas=loga.t.chm[,grepl(paste(',',i,']',sep=''),colnames(loga.t.chm))]
  umsy.st=umsy.st.chm[,i]
  umsy.st=ifelse(umsy.st<0,0.005,umsy.st)
  
  o.logas=apply(logas[,c(chm.info$begin[i]-min(chm.info$begin)+1):c(chm.info$begin[i]-min(chm.info$begin)+5)],1,mean)
  r.logas=apply(logas[,c(chm.info$end[i]-min(chm.info$begin)+1-5):c(chm.info$end[i]-min(chm.info$begin))],1,mean)
  r.umsy=sapply(r.logas,umsyCalc)
  r.umsy=ifelse(r.umsy<0,0.005,r.umsy)
  
  trend_chm_l[,i]=(r.logas-o.logas)/(chm.info$end[i]-chm.info$begin[i]+1-5)
  umsy_comp_chm[i,2]=median(umsy.st)
  umsy_comp_chm[i,3]=HDInterval::hdi(umsy.st,credMass=0.9)[1]
  umsy_comp_chm[i,4]=HDInterval::hdi(umsy.st,credMass=0.9)[2]
  umsy_comp_chm[i,5]=median(r.umsy)
  umsy_comp_chm[i,6]=HDInterval::hdi(r.umsy,credMass=0.9)[1]
  umsy_comp_chm[i,7]=HDInterval::hdi(r.umsy,credMass=0.9)[2]
  umsy_comp_chm[i,8]=median((r.umsy/umsy.st)-1)
}

trend_chm_l.sc=trend_chm_l[,chm.sc]
trend_chm_l.nc=trend_chm_l[,chm.nc]
trend_chm_l.goa=trend_chm_l[,chm.goa]
trend_chm_l.bs=trend_chm_l[,chm.bs]

#Coho
loga.t.cho=tv.cho$draws(variables='log_a_t',format='draws_matrix')
umsy.st.cho=stc.cho$draws(variables='Umsy',format='draws_matrix')

trend_cho_l=matrix(nrow=nrow(umsy.st.cho),ncol=ncol(umsy.st.cho))
umsy_comp_cho=data.frame(stock=cho.info$stock.name,st.umsy=NA,st.umsy.l90=NA,st.umsy.u90=NA,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,p.change=NA)
for(i in 1:nrow(cho.info)){
  logas=loga.t.cho[,grepl(paste(',',i,']',sep=''),colnames(loga.t.cho))]
  umsy.st=umsy.st.cho[,i]
  umsy.st=ifelse(umsy.st<0,0.005,umsy.st)
  
  o.logas=apply(logas[,c(cho.info$begin[i]-min(cho.info$begin)+1):c(cho.info$begin[i]-min(cho.info$begin)+5)],1,mean)
  r.logas=apply(logas[,c(cho.info$end[i]-min(cho.info$begin)+1-5):c(cho.info$end[i]-min(cho.info$begin))],1,mean)
  r.umsy=sapply(r.logas,umsyCalc)
  r.umsy=ifelse(r.umsy<0,0.005,r.umsy)
  
  trend_cho_l[,i]=(r.logas-o.logas)/(cho.info$end[i]-cho.info$begin[i]+1-5)
  umsy_comp_cho[i,2]=median(umsy.st)
  umsy_comp_cho[i,3]=HDInterval::hdi(umsy.st,credMass=0.9)[1]
  umsy_comp_cho[i,4]=HDInterval::hdi(umsy.st,credMass=0.9)[2]
  umsy_comp_cho[i,5]=median(r.umsy)
  umsy_comp_cho[i,6]=HDInterval::hdi(r.umsy,credMass=0.9)[1]
  umsy_comp_cho[i,7]=HDInterval::hdi(r.umsy,credMass=0.9)[2]
  umsy_comp_cho[i,8]=median((r.umsy/umsy.st)-1)
}

trend_cho_l.sc=trend_cho_l[,cho.sc]
trend_cho_l.nc=trend_cho_l[,cho.nc]
trend_cho_l.goa=trend_cho_l[,cho.goa]
trend_cho_l.bs=trend_cho_l[,cho.bs]

#Pink - even
loga.t.pke=tv.pke$draws(variables='log_a_t',format='draws_matrix')
umsy.st.pke=stc.pke$draws(variables='Umsy',format='draws_matrix')

trend_pke_l=matrix(nrow=nrow(umsy.st.pke),ncol=ncol(umsy.st.pke))
umsy_comp_pke=data.frame(stock=pke.info$stock.name,st.umsy=NA,st.umsy.l90=NA,st.umsy.u90=NA,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,p.change=NA)
for(i in 1:nrow(pke.info)){
  logas=loga.t.pke[,grepl(paste(',',i,']',sep=''),colnames(loga.t.pke))]
  umsy.st=umsy.st.pke[,i]
  umsy.st=ifelse(umsy.st<0,0.005,umsy.st)
  
  o.logas=apply(logas[,c(pke.info$begin[i]-min(pke.info$begin)+1):c(pke.info$begin[i]-min(pke.info$begin)+5)],1,mean)
  r.logas=apply(logas[,c(pke.info$end[i]-min(pke.info$begin)+1-5):c(pke.info$end[i]-min(pke.info$begin))],1,mean)
  r.umsy=sapply(r.logas,umsyCalc)
  r.umsy=ifelse(r.umsy<0,0.005,r.umsy)
  
  trend_pke_l[,i]=(r.logas-o.logas)/(pke.info$end[i]-pke.info$begin[i]+1-5)
  umsy_comp_pke[i,2]=median(umsy.st)
  umsy_comp_pke[i,3]=HDInterval::hdi(umsy.st,credMass=0.9)[1]
  umsy_comp_pke[i,4]=HDInterval::hdi(umsy.st,credMass=0.9)[2]
  umsy_comp_pke[i,5]=median(r.umsy)
  umsy_comp_pke[i,6]=HDInterval::hdi(r.umsy,credMass=0.9)[1]
  umsy_comp_pke[i,7]=HDInterval::hdi(r.umsy,credMass=0.9)[2]
  umsy_comp_pke[i,8]=median((r.umsy/umsy.st)-1)
}

trend_pke_l.sc=trend_pke_l[,pke.sc]
trend_pke_l.nc=trend_pke_l[,pke.nc]
trend_pke_l.goa=trend_pke_l[,pke.goa]
trend_pke_l.bs=trend_pke_l[,pke.bs]

#Pink - odd
loga.t.pko=tv.pko$draws(variables='log_a_t',format='draws_matrix')
umsy.st.pko=stc.pko$draws(variables='Umsy',format='draws_matrix')

trend_pko_l=matrix(nrow=nrow(umsy.st.pko),ncol=ncol(umsy.st.pko))
umsy_comp_pko=data.frame(stock=pko.info$stock.name,st.umsy=NA,st.umsy.l90=NA,st.umsy.u90=NA,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,p.change=NA)
for(i in 1:nrow(pko.info)){
  logas=loga.t.pko[,grepl(paste(',',i,']',sep=''),colnames(loga.t.pko))]
  umsy.st=umsy.st.pko[,i]
  umsy.st=ifelse(umsy.st<0,0.005,umsy.st)
  
  o.logas=apply(logas[,c(pko.info$begin[i]-min(pko.info$begin)+1):c(pko.info$begin[i]-min(pko.info$begin)+5)],1,mean)
  r.logas=apply(logas[,c(pko.info$end[i]-min(pko.info$begin)+1-5):c(pko.info$end[i]-min(pko.info$begin))],1,mean)
  r.umsy=sapply(r.logas,umsyCalc)
  r.umsy=ifelse(r.umsy<0,0.005,r.umsy)
  
  trend_pko_l[,i]=(r.logas-o.logas)/(pko.info$end[i]-pko.info$begin[i]+1-5)
  umsy_comp_pko[i,2]=median(umsy.st)
  umsy_comp_pko[i,3]=HDInterval::hdi(umsy.st,credMass=0.9)[1]
  umsy_comp_pko[i,4]=HDInterval::hdi(umsy.st,credMass=0.9)[2]
  umsy_comp_pko[i,5]=median(r.umsy)
  umsy_comp_pko[i,6]=HDInterval::hdi(r.umsy,credMass=0.9)[1]
  umsy_comp_pko[i,7]=HDInterval::hdi(r.umsy,credMass=0.9)[2]
  umsy_comp_pko[i,8]=median((r.umsy/umsy.st)-1)
}

trend_pko_l.sc=trend_pko_l[,pko.sc]
trend_pko_l.nc=trend_pko_l[,pko.nc]
trend_pko_l.goa=trend_pko_l[,pko.goa]
trend_pko_l.bs=trend_pko_l[,pko.bs]

#Fig. 4 Recent vs stationary Umsy ####
pdf('./outputs/main figs/figure_4.pdf',width=8,height=9)
par(mfrow=c(3,2),mar=c(5.1,4.1,2.1,2.1))

umsy.diff.chk=umsy_comp_chk$p.change*100

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts.chk=cut(umsy.diff.chk, breaks = c(seq(-100,100,length.out=11)),right=F)
p_cols2.chk=p_cols1[as.numeric(cuts.chk)]


plot(r.umsy~st.umsy,data=umsy_comp_chk,ylab='',xlab='',bty='l',type='n',xlim=c(0,1),ylim=c(0,1),font.axis=2,las=1)
abline(0,1)
for(i in 1:nrow(umsy_comp_chk)){
  lines(rep(r.umsy[i],2)~c(st.umsy.l90[i],st.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_chk,lwd=0.6)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(st.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_chk,lwd=0.6)
  lines(c(st.umsy[i],r.umsy[i])~rep(st.umsy[i],2),col=adjustcolor(p_cols2.chk[i],alpha.f = 0.8),data=umsy_comp_chk,lwd=1)
  
}
points(r.umsy~st.umsy,data=umsy_comp_chk,bg=adjustcolor(p_cols2.chk,alpha.f=0.8),col=adjustcolor('black',alpha.f=0.6),pch=21,cex=1.4,lwd=0.5)
mtext(expression(bold('Chinook')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('A')),side=3,adj=0,line=0.5,cex=0.75)


u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v1 <-  c(v[1]+(v[2]-v[1])*0.05, v[1]+(v[2]-v[1])*0.3, v[3]+(v[4]-v[3])*0.75, v[4]-(v[4]-v[3])*0.01)

umsy.diff.scy=umsy_comp_scy$p.change*100

cuts.scy=cut(umsy.diff.scy, breaks = c(seq(-100,100,length.out=11)),right=F)
p_cols2.scy=p_cols1[as.numeric(cuts.scy)]

plot(r.umsy~st.umsy,data=umsy_comp_scy,ylab='',xlab='',bty='l',type='n',xlim=c(0,1),ylim=c(0,1),font.axis=2,las=1)
abline(0,1)
for(i in 1:nrow(umsy_comp_scy)){
  lines(rep(r.umsy[i],2)~c(st.umsy.l90[i],st.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_scy,lwd=0.6)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(st.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_scy,lwd=0.6)
  lines(c(st.umsy[i],r.umsy[i])~rep(st.umsy[i],2),col=adjustcolor(p_cols2.scy[i],alpha.f = 0.8),data=umsy_comp_scy,lwd=1)
  
}
points(r.umsy~st.umsy,data=umsy_comp_scy,bg=adjustcolor(p_cols2.scy,alpha.f=0.8),col=adjustcolor('black',alpha.f=0.6),pch=21,cex=1.4,lwd=0.5)
mtext(expression(bold('Sockeye')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('B')),side=3,adj=0,line=0.5,cex=0.75)

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v2 <-  c(v[1]+(v[2]-v[1])*0.05, v[1]+(v[2]-v[1])*0.3, v[3]+(v[4]-v[3])*0.75, v[4]-(v[4]-v[3])*0.01)

umsy.diff.chm=umsy_comp_chm$p.change*100

cuts.chm=cut(umsy.diff.chm, breaks = c(seq(-100,100,length.out=11)),right=F)
p_cols2.chm=p_cols1[as.numeric(cuts.chm)]

plot(r.umsy~st.umsy,data=umsy_comp_chm,ylab='',xlab='',bty='l',type='n',xlim=c(0,1),ylim=c(0,1),font=2,las=1)
abline(0,1)
for(i in 1:nrow(umsy_comp_chm)){
  lines(rep(r.umsy[i],2)~c(st.umsy.l90[i],st.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_chm,lwd=0.6)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(st.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_chm,lwd=0.6)
  lines(c(st.umsy[i],r.umsy[i])~rep(st.umsy[i],2),col=adjustcolor(p_cols2.chm[i],alpha.f = 0.8),data=umsy_comp_chm,lwd=1)
  
}
points(r.umsy~st.umsy,data=umsy_comp_chm,bg=adjustcolor(p_cols2.chm,alpha.f=0.8),col=adjustcolor('black',alpha.f=0.6),pch=21,cex=1.4,lwd=0.5)
mtext(expression(bold('Chum')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('C')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold(paste('Last observed 5-years mean ',italic(U['MSY,j,t']),sep=''))),side=2,adj=0.5,line=2.25,cex=0.9)

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v3 <-  c(v[1]+(v[2]-v[1])*0.05, v[1]+(v[2]-v[1])*0.3, v[3]+(v[4]-v[3])*0.75, v[4]-(v[4]-v[3])*0.01)

umsy.diff.cho=umsy_comp_cho$p.change*100

cuts.cho=cut(umsy.diff.cho, breaks = c(seq(-100,100,length.out=11)),right=F)
p_cols2.cho=p_cols1[as.numeric(cuts.cho)]

plot(r.umsy~st.umsy,data=umsy_comp_cho,ylab='',xlab='',bty='l',type='n',xlim=c(0,1),ylim=c(0,1),font=2,las=1)
abline(0,1)
for(i in 1:nrow(umsy_comp_cho)){
  lines(rep(r.umsy[i],2)~c(st.umsy.l90[i],st.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_cho,lwd=0.6)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(st.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_cho,lwd=0.6)
  lines(c(st.umsy[i],r.umsy[i])~rep(st.umsy[i],2),col=adjustcolor(p_cols2.cho[i],alpha.f = 0.8),data=umsy_comp_cho,lwd=1)
  
}
points(r.umsy~st.umsy,data=umsy_comp_cho,bg=adjustcolor(p_cols2.cho,alpha.f=0.8),col=adjustcolor('black',alpha.f=0.6),pch=21,cex=1.4,lwd=0.5)
mtext(expression(bold('Coho')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('D')),side=3,adj=0,line=0.5,cex=0.75)

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v4 <-  c(v[1]+(v[2]-v[1])*0.05, v[1]+(v[2]-v[1])*0.3, v[3]+(v[4]-v[3])*0.75, v[4]-(v[4]-v[3])*0.01)

umsy.diff.pke=umsy_comp_pke$p.change*100

cuts.pke=cut(umsy.diff.pke, breaks = c(seq(-100,100,length.out=11)),right=F)
p_cols2.pke=p_cols1[as.numeric(cuts.pke)]

plot(r.umsy~st.umsy,data=umsy_comp_pke,ylab='',xlab='',bty='l',type='n',xlim=c(0,1),ylim=c(0,1),font=2,las=1)
abline(0,1)
for(i in 1:nrow(umsy_comp_pke)){
  lines(rep(r.umsy[i],2)~c(st.umsy.l90[i],st.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_pke,lwd=0.6)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(st.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_pke,lwd=0.6)
  lines(c(st.umsy[i],r.umsy[i])~rep(st.umsy[i],2),col=adjustcolor(p_cols2.pke[i],alpha.f = 0.8),data=umsy_comp_pke,lwd=1)
  
}
points(r.umsy~st.umsy,data=umsy_comp_pke,bg=adjustcolor(p_cols2.pke,alpha.f=0.8),col=adjustcolor('black',alpha.f=0.6),pch=21,cex=1.4,lwd=0.5)
mtext(expression(bold('Pink - Even')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('E')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold(paste('Stationary ',italic(U['MSY,j']),sep=''))),side=1,at=1.25,line=3,cex=0.9)

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v5 <-  c(v[1]+(v[2]-v[1])*0.05, v[1]+(v[2]-v[1])*0.3, v[3]+(v[4]-v[3])*0.75, v[4]-(v[4]-v[3])*0.01)

umsy.diff.pko=umsy_comp_pko$p.change*100

cuts.pko=cut(umsy.diff.pko, breaks = c(seq(-100,100,length.out=11)),right=F)
p_cols2.pko=p_cols1[as.numeric(cuts.pko)]

plot(r.umsy~st.umsy,data=umsy_comp_pko,ylab='',xlab='',bty='l',type='n',xlim=c(0,1),ylim=c(0,1),font=2,las=1)
abline(0,1)
for(i in 1:nrow(umsy_comp_pko)){
  lines(rep(r.umsy[i],2)~c(st.umsy.l90[i],st.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_pko,lwd=0.8)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(st.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_pko,lwd=0.8)
  lines(c(st.umsy[i],r.umsy[i])~rep(st.umsy[i],2),col=adjustcolor(p_cols2.pko[i],alpha.f = 0.8),data=umsy_comp_pko,lwd=2)
  
}
points(r.umsy~st.umsy,data=umsy_comp_pko,bg=adjustcolor(p_cols2.pko,alpha.f=0.8),col=adjustcolor('black',alpha.f=0.6),pch=21,cex=1.4,lwd=0.5)
mtext(expression(bold('Pink - Odd')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('F')),side=3,adj=0,line=0.5,cex=0.75)

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v6 <-  c(v[1]+(v[2]-v[1])*0.05, v[1]+(v[2]-v[1])*0.3, v[3]+(v[4]-v[3])*0.75, v[4]-(v[4]-v[3])*0.01)


par(fig=v1, new=TRUE, mar=c(0,0,0,0),xpd=T)

breaks=seq(-100,100,by=10)
break_cols=p_cols1[cut(breaks,c(seq(-100,100,length.out=11)),right=F)]
break_cols=adjustcolor(break_cols,alpha.f=0.8)

hist(umsy.diff.chk,breaks=breaks,xlab='',main='',border='white',col=break_cols,yaxt='n',xaxt='n',ylab='',xlim=c(-100,100),cex.axis=0.5,cex.lab=0.5)
axis(side=1,at=c(seq(-100,100,by=25)),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
axis(side=2,at=seq(0,10,by=5),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
mtext(c('0','5','10'),at=c(0,5,10),cex=0.35,side=2,line=0.3,las=1)
mtext(c('-100','-50','0','+50','+100'),at=c(-100,-50,0,50,100),cex=0.35,side=1,line=0.05)
mtext('% change',side=1,cex=0.45,line=0.8)

par(fig=v2, new=TRUE, mar=c(0,0,0,0),xpd=T)

hist(umsy.diff.scy,breaks=breaks,xlab='',main='',border='white',col=break_cols,yaxt='n',xaxt='n',ylab='',xlim=c(-100,100),cex.axis=0.5,cex.lab=0.5)
axis(side=1,at=c(seq(-100,100,by=25)),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
axis(side=2,at=seq(0,10,by=5),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
mtext(c('0','5','10'),at=c(0,5,10),cex=0.35,side=2,line=0.3,las=1)
mtext(c('-100','-50','0','+50','+100'),at=c(-100,-50,0,50,100),cex=0.35,side=1,line=0.05)
mtext('% change',side=1,cex=0.45,line=0.8)

par(fig=v3, new=TRUE, mar=c(0,0,0,0),xpd=T)

hist(umsy.diff.chm,breaks=breaks,xlab='',main='',border='white',col=break_cols,yaxt='n',xaxt='n',ylab='',xlim=c(-100,100),cex.axis=0.5,cex.lab=0.5)
axis(side=1,at=c(seq(-100,100,by=25)),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
axis(side=2,at=seq(0,5,by=5),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
mtext(c('0','5'),at=c(0,5),cex=0.35,side=2,line=0.3,las=1)
mtext(c('-100','-50','0','+50','+100'),at=c(-100,-50,0,50,100),cex=0.35,side=1,line=0.05)
mtext('% change',side=1,cex=0.45,line=0.8)

par(fig=v4, new=TRUE, mar=c(0,0,0,0),xpd=T)

hist(umsy.diff.cho,breaks=breaks,xlab='',main='',border='white',col=break_cols,yaxt='n',xaxt='n',ylab='',xlim=c(-100,100),cex.axis=0.5,cex.lab=0.5)
axis(side=1,at=c(seq(-100,100,by=25)),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
axis(side=2,at=seq(0,10,by=5),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
mtext(c('0','5','10'),at=c(0,5,10),cex=0.35,side=2,line=0.3,las=1)
mtext(c('-100','-50','0','+50','+100'),at=c(-100,-50,0,50,100),cex=0.35,side=1,line=0.05)
mtext('% change',side=1,cex=0.45,line=0.8)

par(fig=v5, new=TRUE, mar=c(0,0,0,0),xpd=T)

hist(umsy.diff.pke,breaks=breaks,xlab='',main='',border='white',col=break_cols,yaxt='n',xaxt='n',ylab='',xlim=c(-100,100),cex.axis=0.5,cex.lab=0.5)
axis(side=1,at=c(seq(-100,100,by=25)),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
axis(side=2,at=seq(0,5,by=5),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
mtext(c('0','5'),at=c(0,5),cex=0.35,side=2,line=0.3,las=1)
mtext(c('-100','-50','0','+50','+100'),at=c(-100,-50,0,50,100),cex=0.35,side=1,line=0.05)
mtext('% change',side=1,cex=0.45,line=0.8)

par(fig=v6, new=TRUE, mar=c(0,0,0,0),xpd=T)

hist(umsy.diff.pko,breaks=breaks,xlab='',main='',border='white',col=break_cols,yaxt='n',xaxt='n',ylab='',xlim=c(-100,100),cex.axis=0.5,cex.lab=0.5)
axis(side=1,at=c(seq(-100,100,by=25)),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
axis(side=2,at=seq(0,5,by=5),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
mtext(c('0','5'),at=c(0,5),cex=0.35,side=2,line=0.3,las=1)
mtext(c('-100','-50','0','+50','+100'),at=c(-100,-50,0,50,100),cex=0.35,side=1,line=0.05)
mtext('% change',side=1,cex=0.45,line=0.8)

dev.off()


##


# Jpeg versions ####
##fig. 1####
jpeg('./outputs/main figs/figure_1.jpg',width=8,height=6,res=1200,units='in')
main + 
  annotation_custom(
    grob = inset_grob,
    xmin = -141, xmax = -118, ymin = 58.5, ymax = 69
  )
dev.off()

##fig. 2####
jpeg('./outputs/main figs/figure_2.jpg',width=8,height=9,res=1200,units='in')
par(mfrow=c(3,2))
#chinook
par(mar=c(3.1,4.1,2.1,2.1))
log.a.t.chk=exp(tv.chk$draws(variable='mu_log_a',format='draws_matrix'))
m.log.a.t.chk=apply(log.a.t.chk,2,median)
names(m.log.a.t.chk)=seq(min(chk.info$begin),max(chk.info$end))
l95.log.a.t.chk=HDInterval::hdi(log.a.t.chk,credMass=0.95)[1,]
u95.log.a.t.chk=HDInterval::hdi(log.a.t.chk,credMass=0.95)[2,]
log.a.jt.chk=tv.chk$draws(variable='log_a_t',format='draws_matrix')

L=max(chk.info$end)-min(chk.info$begin)+1

plot(m.log.a.t.chk~seq(min(chk.info$begin),max(chk.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)
polygon(
  c(seq(min(chk.info$begin),max(chk.info$end)), rev(seq(min(chk.info$begin),max(chk.info$end)))),           # x coordinates forward and back
  c(u95.log.a.t.chk, rev(l95.log.a.t.chk)),   # y upper then lower reversed
  col = adjustcolor(sp_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
for(j in 1:nrow(chk.info)){
  logaj=apply(log.a.jt.chk[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chk))],2,median)
  st=chk.info$begin[j]-min(chk.info$begin)+1
  end=chk.info$end[j]-min(chk.info$begin)+1
  lines(exp(logaj[1:c(st)])~seq(min(chk.info$begin),names(m.log.a.t.chk[c(st)])),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[st:end])~seq(chk.info$begin[j],chk.info$end[j]),col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[c(end):length(logaj)])~seq(c(chk.info$end[j]),max(chk.info$end)),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
}
for(i in 1:nrow(chk.info)){
  ts=seq(chk.info$begin[i],chk.info$end[i])
  ms=m.log.a.t.chk[names(m.log.a.t.chk)%in%ts]
  lines(ms~ts,col=adjustcolor(sp_cols[1],alpha.f=0.025),lwd=3)
}
for(i in 1:5){lines(rep(7.6,2)~c(1950,1953),col=adjustcolor(sp_cols[1],alpha.f=0.025),lwd=2)}
for(i in 1:10){lines(rep(7.2,2)~c(1950,1953),col=adjustcolor(sp_cols[1],alpha.f=0.025),lwd=2)}
for(i in 1:20){lines(rep(6.8,2)~c(1950,1953),col=adjustcolor(sp_cols[1],alpha.f=0.025),lwd=2)}
for(i in 1:50){lines(rep(6.4,2)~c(1950,1953),col=adjustcolor(sp_cols[1],alpha.f=0.025),lwd=2)}

text(x=1956,y=8,'n',cex=1)
text(x=1956,y=7.6,'5',cex=0.9)
text(x=1956,y=7.2,'10',cex=0.9)
text(x=1956,y=6.8,'20',cex=0.9)
text(x=1956,y=6.4,'50',cex=0.9)
text(x=2022,y=8,paste('N = ',nrow(chk.info)),cex=1)

mtext(expression(bold('A')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Chinook')),side=3,adj=0.5,line=0.5,cex=0.75)

par(mar=c(3.1,2.1,2.1,2.1))
log.a.t.scy=tv.scy$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.scy=apply(exp(log.a.t.scy),2,median)
names(m.log.a.t.scy)=seq(min(scy.info$begin),max(scy.info$end))
l95.log.a.t.scy=HDInterval::hdi(exp(log.a.t.scy),credMass=0.95)[1,]
u95.log.a.t.scy=HDInterval::hdi(exp(log.a.t.scy),credMass=0.95)[2,]
log.a.jt.scy=tv.scy$draws(variable='log_a_t',format='draws_matrix')

L=max(scy.info$end)-min(scy.info$begin)+1

plot(m.log.a.t.scy~seq(min(scy.info$begin),max(scy.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)
polygon(
  c(seq(min(scy.info$begin),max(scy.info$end)), rev(seq(min(scy.info$begin),max(scy.info$end)))),           # x coordinates forward and back
  c(u95.log.a.t.scy, rev(l95.log.a.t.scy)),   # y upper then lower reversed
  col = adjustcolor(sp_cols[5],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
for(j in 1:nrow(scy.info)){
  logaj=apply(log.a.jt.scy[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.scy))],2,median)
  st=scy.info$begin[j]-min(scy.info$begin)+1
  end=scy.info$end[j]-min(scy.info$begin)+1
  lines(exp(logaj[1:c(st)])~seq(min(scy.info$begin),names(m.log.a.t.scy[c(st)])),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[st:end])~seq(scy.info$begin[j],scy.info$end[j]),col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[c(end):length(logaj)])~seq(c(scy.info$end[j]),max(scy.info$end)),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
}
for(i in 1:nrow(scy.info)){
  ts=seq(scy.info$begin[i],scy.info$end[i])
  ms=m.log.a.t.scy[names(m.log.a.t.scy)%in%ts]
  lines(ms~ts,col=adjustcolor(sp_cols[5],alpha.f=0.025),lwd=3)
}
text(x=2022,y=8,paste('N = ',nrow(scy.info)),cex=1)
mtext(expression(bold('B')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Sockeye')),side=3,adj=0.5,line=0.5,cex=0.75)

#Chum
par(mar=c(3.1,4.1,2.1,2.1))
log.a.t.chm=tv.chm$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.chm=apply(exp(log.a.t.chm),2,median)
names(m.log.a.t.chm)=seq(min(chm.info$begin),max(chm.info$end))
l95.log.a.t.chm=HDInterval::hdi(exp(log.a.t.chm),credMass=0.95)[1,]
u95.log.a.t.chm=HDInterval::hdi(exp(log.a.t.chm),credMass=0.95)[2,]
log.a.jt.chm=tv.chm$draws(variable='log_a_t',format='draws_matrix')

L=max(chm.info$end)-min(chm.info$begin)+1

plot(m.log.a.t.chm~seq(min(chm.info$begin),max(chm.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)

polygon(
  c(seq(min(chm.info$begin),max(chm.info$end)), rev(seq(min(chm.info$begin),max(chm.info$end)))),           # x coordinates forward and back
  c(u95.log.a.t.chm, rev(l95.log.a.t.chm)),   # y upper then lower reversed
  col = adjustcolor(sp_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
for(j in 1:nrow(chm.info)){
  logaj=apply(log.a.jt.chm[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.chm))],2,median)
  st=chm.info$begin[j]-min(chm.info$begin)+1
  end=chm.info$end[j]-min(chm.info$begin)+1
  lines(exp(logaj[1:c(st)])~seq(min(chm.info$begin),names(m.log.a.t.chm[c(st)])),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[st:end])~seq(chm.info$begin[j],chm.info$end[j]),col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[c(end):length(logaj)])~seq(c(chm.info$end[j]),max(chm.info$end)),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
}
for(i in 1:nrow(chm.info)){
  ts=seq(chm.info$begin[i],chm.info$end[i])
  ms=m.log.a.t.chm[names(m.log.a.t.chm)%in%ts]
  lines(ms~ts,col=adjustcolor(sp_cols[2],alpha.f=0.025),lwd=3)
}
text(x=2022,y=8,paste('N = ',nrow(chm.info)),cex=1)
mtext(expression(bold('Productivity (max. recruits per spawner)')),side=2,adj=0.5,line=2)
mtext(expression(bold('C')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Chum')),side=3,adj=0.5,line=0.5,cex=0.75)

#Coho
par(mar=c(3.1,2.1,2.1,2.1))
log.a.t.cho=tv.cho$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.cho=apply(exp(log.a.t.cho),2,median)
names(m.log.a.t.cho)=seq(min(cho.info$begin),max(cho.info$end))
l95.log.a.t.cho=HDInterval::hdi(exp(log.a.t.cho),credMass=0.95)[1,]
u95.log.a.t.cho=HDInterval::hdi(exp(log.a.t.cho),credMass=0.95)[2,]
log.a.jt.cho=tv.cho$draws(variable='log_a_t',format='draws_matrix')

L=max(cho.info$end)-min(cho.info$begin)+1

plot(m.log.a.t.cho~seq(min(cho.info$begin),max(cho.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)
polygon(
  c(seq(min(cho.info$begin),max(cho.info$end)), rev(seq(min(cho.info$begin),max(cho.info$end)))),           # x coordinates forward and back
  c(u95.log.a.t.cho, rev(l95.log.a.t.cho)),   # y upper then lower reversed
  col = adjustcolor(sp_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
for(j in 1:nrow(cho.info)){
  logaj=apply(log.a.jt.cho[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.cho))],2,median)
  st=cho.info$begin[j]-min(cho.info$begin)+1
  end=cho.info$end[j]-min(cho.info$begin)+1
  lines(exp(logaj[1:c(st)])~seq(min(cho.info$begin),names(m.log.a.t.cho[c(st)])),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[st:end])~seq(cho.info$begin[j],cho.info$end[j]),col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[c(end):length(logaj)])~seq(c(cho.info$end[j]),max(cho.info$end)),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
}
for(i in 1:nrow(cho.info)){
  ts=seq(cho.info$begin[i],cho.info$end[i])
  ms=m.log.a.t.cho[names(m.log.a.t.cho)%in%ts]
  lines(ms~ts,col=adjustcolor(sp_cols[3],alpha.f=0.025),lwd=3)
}
text(x=2022,y=8,paste('N = ',nrow(cho.info)),cex=1)
mtext(expression(bold('D')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Coho')),side=3,adj=0.5,line=0.5,cex=0.75)

#Pink - even
par(mar=c(5.1,4.1,2.1,2.1))
log.a.t.pke=tv.pke$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.pke=apply(exp(log.a.t.pke),2,median)
names(m.log.a.t.pke)=seq(min(pke.info$begin),max(pke.info$end))

l95.log.a.t.pke=HDInterval::hdi(exp(log.a.t.pke),credMass=0.95)[1,]
u95.log.a.t.pke=HDInterval::hdi(exp(log.a.t.pke),credMass=0.95)[2,]
log.a.jt.pke=tv.pke$draws(variable='log_a_t',format='draws_matrix')

L=max(pke.info$end)-min(pke.info$begin)+1

plot(m.log.a.t.pke~seq(min(pke.info$begin),max(pke.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)

polygon(
  c(seq(min(pke.info$begin),max(pke.info$end)), rev(seq(min(pke.info$begin),max(pke.info$end)))),           # x coordinates forward and back
  c(u95.log.a.t.pke, rev(l95.log.a.t.pke)),   # y upper then lower reversed
  col = adjustcolor(sp_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
for(j in 1:nrow(pke.info)){
  logaj=apply(log.a.jt.pke[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.pke))],2,median)
  st=pke.info$begin[j]-min(pke.info$begin)+1
  end=pke.info$end[j]-min(pke.info$begin)+1
  lines(exp(logaj[1:c(st)])~seq(min(pke.info$begin),names(m.log.a.t.pke[c(st)])),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[st:end])~seq(pke.info$begin[j],pke.info$end[j]),col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[c(end):length(logaj)])~seq(c(pke.info$end[j]),max(pke.info$end)),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
}
for(i in 1:nrow(pke.info)){
  ts=seq(pke.info$begin[i],pke.info$end[i])
  ms=m.log.a.t.pke[names(m.log.a.t.pke)%in%ts]
  lines(ms~ts,col=adjustcolor(sp_cols[4],alpha.f=0.025),lwd=3)
}
text(x=2022,y=8,paste('N = ',nrow(pke.info)),cex=1)
mtext(expression(bold('E')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Pink - Even')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('Brood cohort year')),side=1,at=2035,line=3)

#Pink - odd
par(mar=c(5.1,2.1,2.1,2.1))
log.a.t.pko=tv.pko$draws(variable='mu_log_a',format='draws_matrix')
m.log.a.t.pko=apply(exp(log.a.t.pko),2,median)
names(m.log.a.t.pko)=seq(min(pko.info$begin),max(pko.info$end))

l95.log.a.t.pko=HDInterval::hdi(exp(log.a.t.pko),credMass=0.95)[1,]
u95.log.a.t.pko=HDInterval::hdi(exp(log.a.t.pko),credMass=0.95)[2,]
log.a.jt.pko=tv.pko$draws(variable='log_a_t',format='draws_matrix')

L=max(pko.info$end)-min(pko.info$begin)+1

plot(m.log.a.t.pko~seq(min(pko.info$begin),max(pko.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)
polygon(
  c(seq(min(pko.info$begin),max(pko.info$end)), rev(seq(min(pko.info$begin),max(pko.info$end)))),           # x coordinates forward and back
  c(u95.log.a.t.pko, rev(l95.log.a.t.pko)),   # y upper then lower reversed
  col = adjustcolor(sp_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
for(j in 1:nrow(pko.info)){
  logaj=apply(log.a.jt.pko[,grepl(paste(',',j,']',sep=''),colnames(log.a.jt.pko))],2,median)
  st=pko.info$begin[j]-min(pko.info$begin)+1
  end=pko.info$end[j]-min(pko.info$begin)+1
  lines(exp(logaj[1:c(st)])~seq(min(pko.info$begin),names(m.log.a.t.pko[c(st)])),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[st:end])~seq(pko.info$begin[j],pko.info$end[j]),col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
  lines(exp(logaj[c(end):length(logaj)])~seq(c(pko.info$end[j]),max(pko.info$end)),lty=3,col=adjustcolor('black',alpha.f=0.1),lwd=0.5)
}
for(i in 1:nrow(pko.info)){
  ts=seq(pko.info$begin[i],pko.info$end[i])
  ms=m.log.a.t.pko[names(m.log.a.t.pko)%in%ts]
  lines(ms~ts,col=adjustcolor(sp_cols[4],alpha.f=0.025),lwd=3)
}
text(x=2022,y=8,paste('N = ',nrow(pko.info)),cex=1)
mtext(expression(bold('F')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Pink - Odd')),side=3,adj=0.5,line=0.5,cex=0.75)
dev.off()

##fig. 3####
jpeg('./outputs/main figs/figure_3.jpg',width=8,height=9,res=1200,units='in')
par(mfrow=c(3,2))
par(mar=c(3.1,4.1,2.1,2.1))
#Mean productivity trends by ocean basin
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

plot(apply(exp(mu_loga.chk.bs),2,median)~seq(min(chk.info$begin),max(chk.info$end)),ylim=c(0,8),ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)

l95.chk.bs=HDInterval::hdi(exp(mu_loga.chk.bs),credMass=0.9)[1,]
u95.chk.bs=HDInterval::hdi(exp(mu_loga.chk.bs),credMass=0.9)[2,]
l95.chk.goa=HDInterval::hdi(exp(mu_loga.chk.goa),credMass=0.9)[1,]
u95.chk.goa=HDInterval::hdi(exp(mu_loga.chk.goa),credMass=0.9)[2,]
l95.chk.nc=HDInterval::hdi(exp(mu_loga.chk.nc),credMass=0.9)[1,]
u95.chk.nc=HDInterval::hdi(exp(mu_loga.chk.nc),credMass=0.9)[2,]
l95.chk.sc=HDInterval::hdi(exp(mu_loga.chk.sc),credMass=0.9)[1,]
u95.chk.sc=HDInterval::hdi(exp(mu_loga.chk.sc),credMass=0.9)[2,]
polygon(
  c(seq(min(chk.info$begin),max(chk.info$end)), rev(seq(min(chk.info$begin),max(chk.info$end)))),           # x coordinates forward and back
  c(u95.chk.bs, rev(l95.chk.bs)),   # y upper then lower reversed
  col = adjustcolor(basin_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(chk.info$begin),max(chk.info$end)), rev(seq(min(chk.info$begin),max(chk.info$end)))),           # x coordinates forward and back
  c(u95.chk.goa, rev(l95.chk.goa)), # y upper then lower reversed
  col = adjustcolor(basin_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(chk.info$begin),max(chk.info$end)), rev(seq(min(chk.info$begin),max(chk.info$end)))),           # x coordinates forward and back
  c(u95.chk.nc, rev(l95.chk.nc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(chk.info$begin),max(chk.info$end)), rev(seq(min(chk.info$begin),max(chk.info$end)))),           # x coordinates forward and back
  c(u95.chk.sc, rev(l95.chk.sc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
mu_loga.chk.sc.rs=apply(exp(mu_loga.chk.sc),2,median)
names(mu_loga.chk.sc.rs)=seq(min(chk.info$begin),max(chk.info$end))
mu_loga.chk.nc.rs=apply(exp(mu_loga.chk.nc),2,median)
names(mu_loga.chk.nc.rs)=seq(min(chk.info$begin),max(chk.info$end))
mu_loga.chk.goa.rs=apply(exp(mu_loga.chk.goa),2,median)
names(mu_loga.chk.goa.rs)=seq(min(chk.info$begin),max(chk.info$end))
mu_loga.chk.bs.rs=apply(exp(mu_loga.chk.bs),2,median)
names(mu_loga.chk.bs.rs)=seq(min(chk.info$begin),max(chk.info$end))

for(i in 1:nrow(chk.info)){
  ts=seq(chk.info$begin[i],chk.info$end[i])
  if(chk.info$ocean.basin[i]=='SC'){
    ms=mu_loga.chk.sc.rs[names(mu_loga.chk.sc.rs)%in%ts]
    lines(ms~ts,col=adjustcolor(basin_cols[1],alpha.f=0.1),lwd=3)
  }else
    if(chk.info$ocean.basin[i]=='NC'){
      ms=mu_loga.chk.nc.rs[names(mu_loga.chk.nc.rs)%in%ts]
      lines(ms~ts,col=adjustcolor(basin_cols[2],alpha.f=0.1),lwd=3)
    }else
      if(chk.info$ocean.basin[i]=='GoA'){
        ms=mu_loga.chk.goa.rs[names(mu_loga.chk.goa.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[3],alpha.f=0.1),lwd=3)
      }else{
        ms=mu_loga.chk.bs.rs[names(mu_loga.chk.bs.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=3)
      }
}
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.25,'N',cex=1)
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.8,nrow(subset(chk.info,ocean.basin=='BS')),cex=1,col=basin_cols[4])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-1.4,nrow(subset(chk.info,ocean.basin=='GoA')),cex=1,col=basin_cols[3])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.0,nrow(subset(chk.info,ocean.basin=='NC')),cex=1,col=basin_cols[2])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.6,nrow(subset(chk.info,ocean.basin=='SC')),cex=1,col=basin_cols[1])

text(x=par('usr')[2]-6,y=par('usr')[4]-0.8,expression(bold('Bering Sea')),col=basin_cols[4],cex=1.25,adj=1)
text(x=par('usr')[2]-6,y=par('usr')[4]-1.4,expression(bold('Gulf of Alaska')),col=basin_cols[3],cex=1.25,adj=1)
text(x=par('usr')[2]-6,y=par('usr')[4]-2.0,expression(bold('North Coast')),col=basin_cols[2],cex=1.25,adj=1)
text(x=par('usr')[2]-6,y=par('usr')[4]-2.6,expression(bold('South Coast')),col=basin_cols[1],cex=1.25,adj=1)
lines(rep(7.6,2)~c(1950,1953),col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=2)
for(i in 1:5){lines(rep(7.2,2)~c(1950,1953),col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=2)}
for(i in 1:10){lines(rep(6.8,2)~c(1950,1953),col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=2)}
for(i in 1:20){lines(rep(6.4,2)~c(1950,1953),col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=2)}
text(x=1956,y=8,'n',cex=1)
text(x=1956,y=7.6,'1',cex=0.9)
text(x=1956,y=7.2,'5',cex=0.9)
text(x=1956,y=6.8,'10',cex=0.9)
text(x=1956,y=6.4,'20',cex=0.9)

mtext(expression(bold('A')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Chinook')),side=3,adj=0.5,line=0.5,cex=0.75)

par(mar=c(3.1,2.1,2.1,2.1))
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

plot(apply(exp(mu_loga.scy.bs),2,median)~seq(min(scy.info$begin),max(scy.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)

l95.scy.bs=HDInterval::hdi(exp(mu_loga.scy.bs),credMass=0.9)[1,]
u95.scy.bs=HDInterval::hdi(exp(mu_loga.scy.bs),credMass=0.9)[2,]
l95.scy.goa=HDInterval::hdi(exp(mu_loga.scy.goa),credMass=0.9)[1,]
u95.scy.goa=HDInterval::hdi(exp(mu_loga.scy.goa),credMass=0.9)[2,]
l95.scy.nc=HDInterval::hdi(exp(mu_loga.scy.nc),credMass=0.9)[1,]
u95.scy.nc=HDInterval::hdi(exp(mu_loga.scy.nc),credMass=0.9)[2,]
l95.scy.sc=HDInterval::hdi(exp(mu_loga.scy.sc),credMass=0.9)[1,]
u95.scy.sc=HDInterval::hdi(exp(mu_loga.scy.sc),credMass=0.9)[2,]
polygon(
  c(seq(min(scy.info$begin),max(scy.info$end)), rev(seq(min(scy.info$begin),max(scy.info$end)))),           # x coordinates forward and back
  c(u95.scy.bs, rev(l95.scy.bs)),   # y upper then lower reversed
  col = adjustcolor(basin_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(scy.info$begin),max(scy.info$end)), rev(seq(min(scy.info$begin),max(scy.info$end)))),           # x coordinates forward and back
  c(u95.scy.goa, rev(l95.scy.goa)), # y upper then lower reversed
  col = adjustcolor(basin_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(scy.info$begin),max(scy.info$end)), rev(seq(min(scy.info$begin),max(scy.info$end)))),           # x coordinates forward and back
  c(u95.scy.nc, rev(l95.scy.nc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(scy.info$begin),max(scy.info$end)), rev(seq(min(scy.info$begin),max(scy.info$end)))),           # x coordinates forward and back
  c(u95.scy.sc, rev(l95.scy.sc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
mu_loga.scy.sc.rs=apply(exp(mu_loga.scy.sc),2,median)
names(mu_loga.scy.sc.rs)=seq(min(scy.info$begin),max(scy.info$end))
mu_loga.scy.nc.rs=apply(exp(mu_loga.scy.nc),2,median)
names(mu_loga.scy.nc.rs)=seq(min(scy.info$begin),max(scy.info$end))
mu_loga.scy.goa.rs=apply(exp(mu_loga.scy.goa),2,median)
names(mu_loga.scy.goa.rs)=seq(min(scy.info$begin),max(scy.info$end))
mu_loga.scy.bs.rs=apply(exp(mu_loga.scy.bs),2,median)
names(mu_loga.scy.bs.rs)=seq(min(scy.info$begin),max(scy.info$end))

for(i in 1:nrow(scy.info)){
  ts=seq(scy.info$begin[i],scy.info$end[i])
  if(scy.info$ocean.basin[i]=='SC'){
    ms=mu_loga.scy.sc.rs[names(mu_loga.scy.sc.rs)%in%ts]
    lines(ms~ts,col=adjustcolor(basin_cols[1],alpha.f=0.1),lwd=3)
  }else
    if(scy.info$ocean.basin[i]=='NC'){
      ms=mu_loga.scy.nc.rs[names(mu_loga.scy.nc.rs)%in%ts]
      lines(ms~ts,col=adjustcolor(basin_cols[2],alpha.f=0.1),lwd=3)
    }else
      if(scy.info$ocean.basin[i]=='GoA'){
        ms=mu_loga.scy.goa.rs[names(mu_loga.scy.goa.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[3],alpha.f=0.1),lwd=3)
      }else{
        ms=mu_loga.scy.bs.rs[names(mu_loga.scy.bs.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=3)
      }
}
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.25,'N',cex=1)
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.8,nrow(subset(scy.info,ocean.basin=='BS')),cex=1,col=basin_cols[4])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-1.4,nrow(subset(scy.info,ocean.basin=='GoA')),cex=1,col=basin_cols[3])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.0,nrow(subset(scy.info,ocean.basin=='NC')),cex=1,col=basin_cols[2])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.6,nrow(subset(scy.info,ocean.basin=='SC')),cex=1,col=basin_cols[1])

mtext(expression(bold('B')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Sockeye')),side=3,adj=0.5,line=0.5,cex=0.75)

par(mar=c(3.1,4.1,2.1,2.1))
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

plot(apply(exp(mu_loga.chm.bs),2,median)~seq(min(chm.info$begin),max(chm.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)

l95.chm.bs=HDInterval::hdi(exp(mu_loga.chm.bs),credMass=0.9)[1,]
u95.chm.bs=HDInterval::hdi(exp(mu_loga.chm.bs),credMass=0.9)[2,]
l95.chm.goa=HDInterval::hdi(exp(mu_loga.chm.goa),credMass=0.9)[1,]
u95.chm.goa=HDInterval::hdi(exp(mu_loga.chm.goa),credMass=0.9)[2,]
l95.chm.nc=HDInterval::hdi(exp(mu_loga.chm.nc),credMass=0.9)[1,]
u95.chm.nc=HDInterval::hdi(exp(mu_loga.chm.nc),credMass=0.9)[2,]
l95.chm.sc=HDInterval::hdi(exp(mu_loga.chm.sc),credMass=0.9)[1,]
u95.chm.sc=HDInterval::hdi(exp(mu_loga.chm.sc),credMass=0.9)[2,]
polygon(
  c(seq(min(chm.info$begin),max(chm.info$end)), rev(seq(min(chm.info$begin),max(chm.info$end)))),           # x coordinates forward and back
  c(u95.chm.bs, rev(l95.chm.bs)),   # y upper then lower reversed
  col = adjustcolor(basin_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(chm.info$begin),max(chm.info$end)), rev(seq(min(chm.info$begin),max(chm.info$end)))),           # x coordinates forward and back
  c(u95.chm.goa, rev(l95.chm.goa)), # y upper then lower reversed
  col = adjustcolor(basin_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(chm.info$begin),max(chm.info$end)), rev(seq(min(chm.info$begin),max(chm.info$end)))),           # x coordinates forward and back
  c(u95.chm.nc, rev(l95.chm.nc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(chm.info$begin),max(chm.info$end)), rev(seq(min(chm.info$begin),max(chm.info$end)))),           # x coordinates forward and back
  c(u95.chm.sc, rev(l95.chm.sc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
mu_loga.chm.sc.rs=apply(exp(mu_loga.chm.sc),2,median)
names(mu_loga.chm.sc.rs)=seq(min(chm.info$begin),max(chm.info$end))
mu_loga.chm.nc.rs=apply(exp(mu_loga.chm.nc),2,median)
names(mu_loga.chm.nc.rs)=seq(min(chm.info$begin),max(chm.info$end))
mu_loga.chm.goa.rs=apply(exp(mu_loga.chm.goa),2,median)
names(mu_loga.chm.goa.rs)=seq(min(chm.info$begin),max(chm.info$end))
mu_loga.chm.bs.rs=apply(exp(mu_loga.chm.bs),2,median)
names(mu_loga.chm.bs.rs)=seq(min(chm.info$begin),max(chm.info$end))

for(i in 1:nrow(chm.info)){
  ts=seq(chm.info$begin[i],chm.info$end[i])
  if(chm.info$ocean.basin[i]=='SC'){
    ms=mu_loga.chm.sc.rs[names(mu_loga.chm.sc.rs)%in%ts]
    lines(ms~ts,col=adjustcolor(basin_cols[1],alpha.f=0.1),lwd=3)
  }else
    if(chm.info$ocean.basin[i]=='NC'){
      ms=mu_loga.chm.nc.rs[names(mu_loga.chm.nc.rs)%in%ts]
      lines(ms~ts,col=adjustcolor(basin_cols[2],alpha.f=0.1),lwd=3)
    }else
      if(chm.info$ocean.basin[i]=='GoA'){
        ms=mu_loga.chm.goa.rs[names(mu_loga.chm.goa.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[3],alpha.f=0.1),lwd=3)
      }else{
        ms=mu_loga.chm.bs.rs[names(mu_loga.chm.bs.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=3)
      }
}
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.25,'N',cex=1)
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.8,nrow(subset(chm.info,ocean.basin=='BS')),cex=1,col=basin_cols[4])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-1.4,nrow(subset(chm.info,ocean.basin=='GoA')),cex=1,col=basin_cols[3])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.0,nrow(subset(chm.info,ocean.basin=='NC')),cex=1,col=basin_cols[2])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.6,nrow(subset(chm.info,ocean.basin=='SC')),cex=1,col=basin_cols[1])

mtext(expression(bold('Mean productivity (max. recruits per spawner)')),side=2,adj=0.5,line=2)
mtext(expression(bold('C')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Chum')),side=3,adj=0.5,line=0.5,cex=0.75)

par(mar=c(3.1,2.1,2.1,2.1))
log.a.jt.cho=tv.cho$draws(variable='log_a_t',format='draws_matrix')
cho.sc=which(cho.info$ocean.basin=='SC')
cho.nc=which(cho.info$ocean.basin=='NC')
cho.goa=which(cho.info$ocean.basin=='GoA')

logat=log.a.jt.cho
colnames(logat)=gsub("]", ".", colnames(logat))
colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)

logat.sc=logat[,grepl(paste(paste(',',cho.sc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.nc=logat[,grepl(paste(paste(',',cho.nc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.goa=logat[,grepl(paste(paste(',',cho.goa,'.',sep=''),collapse = "|"),colnames(logat))]

L=max(cho.info$end)-min(cho.info$begin)+1
mu_loga.cho.sc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.cho.nc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.cho.goa=matrix(nrow=nrow(logat),ncol=L)
for(t in 1:L){
  logas.sc=logat.sc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.sc))]
  mu_loga.cho.sc[,t]=apply(logas.sc,1,mean)
  
  logas.nc=logat.nc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.nc))]
  mu_loga.cho.nc[,t]=apply(logas.nc,1,mean)
  
  logas.goa=logat.goa[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.goa))]
  mu_loga.cho.goa[,t]=apply(logas.goa,1,mean)
  
}

plot(apply(exp(mu_loga.cho.sc),2,median)~seq(min(cho.info$begin),max(cho.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)

l95.cho.goa=HDInterval::hdi(exp(mu_loga.cho.goa),credMass=0.9)[1,]
u95.cho.goa=HDInterval::hdi(exp(mu_loga.cho.goa),credMass=0.9)[2,]
l95.cho.nc=HDInterval::hdi(exp(mu_loga.cho.nc),credMass=0.9)[1,]
u95.cho.nc=HDInterval::hdi(exp(mu_loga.cho.nc),credMass=0.9)[2,]
l95.cho.sc=HDInterval::hdi(exp(mu_loga.cho.sc),credMass=0.9)[1,]
u95.cho.sc=HDInterval::hdi(exp(mu_loga.cho.sc),credMass=0.9)[2,]
polygon(
  c(seq(min(cho.info$begin),max(cho.info$end)), rev(seq(min(cho.info$begin),max(cho.info$end)))),           # x coordinates forward and back
  c(u95.cho.goa, rev(l95.cho.goa)), # y upper then lower reversed
  col = adjustcolor(basin_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(cho.info$begin),max(cho.info$end)), rev(seq(min(cho.info$begin),max(cho.info$end)))),           # x coordinates forward and back
  c(u95.cho.nc, rev(l95.cho.nc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(cho.info$begin),max(cho.info$end)), rev(seq(min(cho.info$begin),max(cho.info$end)))),           # x coordinates forward and back
  c(u95.cho.sc, rev(l95.cho.sc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
mu_loga.cho.sc.rs=apply(exp(mu_loga.cho.sc),2,median)
names(mu_loga.cho.sc.rs)=seq(min(cho.info$begin),max(cho.info$end))
mu_loga.cho.nc.rs=apply(exp(mu_loga.cho.nc),2,median)
names(mu_loga.cho.nc.rs)=seq(min(cho.info$begin),max(cho.info$end))
mu_loga.cho.goa.rs=apply(exp(mu_loga.cho.goa),2,median)
names(mu_loga.cho.goa.rs)=seq(min(cho.info$begin),max(cho.info$end))

for(i in 1:nrow(cho.info)){
  ts=seq(cho.info$begin[i],cho.info$end[i])
  if(cho.info$ocean.basin[i]=='SC'){
    ms=mu_loga.cho.sc.rs[names(mu_loga.cho.sc.rs)%in%ts]
    lines(ms~ts,col=adjustcolor(basin_cols[1],alpha.f=0.1),lwd=3)
  }else
    if(cho.info$ocean.basin[i]=='NC'){
      ms=mu_loga.cho.nc.rs[names(mu_loga.cho.nc.rs)%in%ts]
      lines(ms~ts,col=adjustcolor(basin_cols[2],alpha.f=0.1),lwd=3)
    }else
      if(cho.info$ocean.basin[i]=='GoA'){
        ms=mu_loga.cho.goa.rs[names(mu_loga.cho.goa.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[3],alpha.f=0.1),lwd=3)
      }
}
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.25,'N',cex=1)
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.8,nrow(subset(cho.info,ocean.basin=='GoA')),cex=1,col=basin_cols[3])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-1.4,nrow(subset(cho.info,ocean.basin=='NC')),cex=1,col=basin_cols[2])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.0,nrow(subset(cho.info,ocean.basin=='SC')),cex=1,col=basin_cols[1])

mtext(expression(bold('D')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Coho')),side=3,adj=0.5,line=0.5,cex=0.75)

par(mar=c(5.1,4.1,2.1,2.1))
log.a.jt.pke=tv.pke$draws(variable='log_a_t',format='draws_matrix')

logat=log.a.jt.pke
colnames(logat)=gsub("]", ".", colnames(logat))
colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)

logat.sc=logat[,grepl(paste(paste(',',pke.sc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.nc=logat[,grepl(paste(paste(',',pke.nc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.goa=logat[,grepl(paste(paste(',',pke.goa,'.',sep=''),collapse = "|"),colnames(logat))]
logat.bs=logat[,grepl(paste(paste(',',pke.bs,'.',sep=''),collapse = "|"),colnames(logat))]

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

plot(apply(exp(mu_loga.pke.bs),2,median)~seq(min(pke.info$begin),max(pke.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)
l95.pke.bs=HDInterval::hdi(exp(mu_loga.pke.bs),credMass=0.9)[1,]
u95.pke.bs=HDInterval::hdi(exp(mu_loga.pke.bs),credMass=0.9)[2,]
l95.pke.goa=HDInterval::hdi(exp(mu_loga.pke.goa),credMass=0.9)[1,]
u95.pke.goa=HDInterval::hdi(exp(mu_loga.pke.goa),credMass=0.9)[2,]
l95.pke.nc=HDInterval::hdi(exp(mu_loga.pke.nc),credMass=0.9)[1,]
u95.pke.nc=HDInterval::hdi(exp(mu_loga.pke.nc),credMass=0.9)[2,]
l95.pke.sc=HDInterval::hdi(exp(mu_loga.pke.sc),credMass=0.9)[1,]
u95.pke.sc=HDInterval::hdi(exp(mu_loga.pke.sc),credMass=0.9)[2,]
polygon(
  c(seq(min(pke.info$begin),max(pke.info$end)), rev(seq(min(pke.info$begin),max(pke.info$end)))),           # x coordinates forward and back
  c(u95.pke.bs, rev(l95.pke.bs)),   # y upper then lower reversed
  col = adjustcolor(basin_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(pke.info$begin),max(pke.info$end)), rev(seq(min(pke.info$begin),max(pke.info$end)))),           # x coordinates forward and back
  c(u95.pke.goa, rev(l95.pke.goa)), # y upper then lower reversed
  col = adjustcolor(basin_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(pke.info$begin),max(pke.info$end)), rev(seq(min(pke.info$begin),max(pke.info$end)))),           # x coordinates forward and back
  c(u95.pke.nc, rev(l95.pke.nc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(pke.info$begin),max(pke.info$end)), rev(seq(min(pke.info$begin),max(pke.info$end)))),           # x coordinates forward and back
  c(u95.pke.sc, rev(l95.pke.sc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
mu_loga.pke.sc.rs=apply(exp(mu_loga.pke.sc),2,median)
names(mu_loga.pke.sc.rs)=seq(min(pke.info$begin),max(pke.info$end))
mu_loga.pke.nc.rs=apply(exp(mu_loga.pke.nc),2,median)
names(mu_loga.pke.nc.rs)=seq(min(pke.info$begin),max(pke.info$end))
mu_loga.pke.goa.rs=apply(exp(mu_loga.pke.goa),2,median)
names(mu_loga.pke.goa.rs)=seq(min(pke.info$begin),max(pke.info$end))
mu_loga.pke.bs.rs=apply(exp(mu_loga.pke.bs),2,median)
names(mu_loga.pke.bs.rs)=seq(min(pke.info$begin),max(pke.info$end))

for(i in 1:nrow(pke.info)){
  ts=seq(pke.info$begin[i],pke.info$end[i])
  if(pke.info$ocean.basin[i]=='SC'){
    ms=mu_loga.pke.sc.rs[names(mu_loga.pke.sc.rs)%in%ts]
    lines(ms~ts,col=adjustcolor(basin_cols[1],alpha.f=0.1),lwd=3)
  }else
    if(pke.info$ocean.basin[i]=='NC'){
      ms=mu_loga.pke.nc.rs[names(mu_loga.pke.nc.rs)%in%ts]
      lines(ms~ts,col=adjustcolor(basin_cols[2],alpha.f=0.1),lwd=3)
    }else
      if(pke.info$ocean.basin[i]=='GoA'){
        ms=mu_loga.pke.goa.rs[names(mu_loga.pke.goa.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[3],alpha.f=0.1),lwd=3)
      }else{
        ms=mu_loga.pke.bs.rs[names(mu_loga.pke.bs.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=3)
      }
}
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.25,'N',cex=1)
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.8,nrow(subset(pke.info,ocean.basin=='BS')),cex=1,col=basin_cols[4])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-1.4,nrow(subset(pke.info,ocean.basin=='GoA')),cex=1,col=basin_cols[3])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.0,nrow(subset(pke.info,ocean.basin=='NC')),cex=1,col=basin_cols[2])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.6,nrow(subset(pke.info,ocean.basin=='SC')),cex=1,col=basin_cols[1])

mtext(expression(bold('E')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Pink - Even')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('Brood cohort year')),side=1,at=2035,line=3)

par(mar=c(5.1,2.1,2.1,2.1))
log.a.jt.pko=tv.pko$draws(variable='log_a_t',format='draws_matrix')

logat=log.a.jt.pko
colnames(logat)=gsub("]", ".", colnames(logat))
colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)

logat.sc=logat[,grepl(paste(paste(',',pko.sc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.nc=logat[,grepl(paste(paste(',',pko.nc,'.',sep=''),collapse = "|"),colnames(logat))]
logat.goa=logat[,grepl(paste(paste(',',pko.goa,'.',sep=''),collapse = "|"),colnames(logat))]
logat.bs=logat[,grepl(paste(paste(',',pko.bs,'.',sep=''),collapse = "|"),colnames(logat))]


L=max(pko.info$end)-min(pko.info$begin)+1
mu_loga.pko.sc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.pko.nc=matrix(nrow=nrow(logat),ncol=L)
mu_loga.pko.goa=matrix(nrow=nrow(logat),ncol=L)
mu_loga.pko.bs=matrix(nrow=nrow(logat),ncol=L)
for(t in 1:L){
  logas.sc=logat.sc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.sc))]
  mu_loga.pko.sc[,t]=apply(logas.sc,1,mean)
  
  logas.nc=logat.nc[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.nc))]
  mu_loga.pko.nc[,t]=apply(logas.nc,1,mean)
  
  logas.goa=logat.goa[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.goa))]
  mu_loga.pko.goa[,t]=apply(logas.goa,1,mean)
  
  logas.bs=logat.bs[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat.bs))]
  mu_loga.pko.bs[,t]=apply(logas.bs,1,mean)
}

plot(apply(exp(mu_loga.pko.bs),2,median)~seq(min(pko.info$begin),max(pko.info$end)),ylim=c(0,8),lwd=3,col='navy',ylab='',xlab='',type='n',main='',xlim=c(1950,2025),yaxt='n',xaxt='n')
axis(side=1,at=seq(1950,2025,by=5),font=2,cex=0.6,labels=ifelse(seq(1950,2025,by=5) %% 10 == 0, seq(1950,2025,by=5), ""))
axis(side=2,at=seq(0,8,by=1),font=2,las=1,cex=0.6)
l95.pko.bs=HDInterval::hdi(exp(mu_loga.pko.bs),credMass=0.9)[1,]
u95.pko.bs=HDInterval::hdi(exp(mu_loga.pko.bs),credMass=0.9)[2,]
l95.pko.goa=HDInterval::hdi(exp(mu_loga.pko.goa),credMass=0.9)[1,]
u95.pko.goa=HDInterval::hdi(exp(mu_loga.pko.goa),credMass=0.9)[2,]
l95.pko.nc=HDInterval::hdi(exp(mu_loga.pko.nc),credMass=0.9)[1,]
u95.pko.nc=HDInterval::hdi(exp(mu_loga.pko.nc),credMass=0.9)[2,]
l95.pko.sc=HDInterval::hdi(exp(mu_loga.pko.sc),credMass=0.9)[1,]
u95.pko.sc=HDInterval::hdi(exp(mu_loga.pko.sc),credMass=0.9)[2,]
polygon(
  c(seq(min(pko.info$begin),max(pko.info$end)), rev(seq(min(pko.info$begin),max(pko.info$end)))),           # x coordinates forward and back
  c(u95.pko.bs, rev(l95.pko.bs)),   # y upper then lower reversed
  col = adjustcolor(basin_cols[4],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(pko.info$begin),max(pko.info$end)), rev(seq(min(pko.info$begin),max(pko.info$end)))),           # x coordinates forward and back
  c(u95.pko.goa, rev(l95.pko.goa)), # y upper then lower reversed
  col = adjustcolor(basin_cols[3],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(pko.info$begin),max(pko.info$end)), rev(seq(min(pko.info$begin),max(pko.info$end)))),           # x coordinates forward and back
  c(u95.pko.nc, rev(l95.pko.nc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[2],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
polygon(
  c(seq(min(pko.info$begin),max(pko.info$end)), rev(seq(min(pko.info$begin),max(pko.info$end)))),           # x coordinates forward and back
  c(u95.pko.sc, rev(l95.pko.sc)), # y upper then lower reversed
  col = adjustcolor(basin_cols[1],alpha.f=0.15),  # transparent fill
  border = NA             # no border line
)
mu_loga.pko.sc.rs=apply(exp(mu_loga.pko.sc),2,median)
names(mu_loga.pko.sc.rs)=seq(min(pko.info$begin),max(pko.info$end))
mu_loga.pko.nc.rs=apply(exp(mu_loga.pko.nc),2,median)
names(mu_loga.pko.nc.rs)=seq(min(pko.info$begin),max(pko.info$end))
mu_loga.pko.goa.rs=apply(exp(mu_loga.pko.goa),2,median)
names(mu_loga.pko.goa.rs)=seq(min(pko.info$begin),max(pko.info$end))
mu_loga.pko.bs.rs=apply(exp(mu_loga.pko.bs),2,median)
names(mu_loga.pko.bs.rs)=seq(min(pko.info$begin),max(pko.info$end))

for(i in 1:nrow(pko.info)){
  ts=seq(pko.info$begin[i],pko.info$end[i])
  if(pko.info$ocean.basin[i]=='SC'){
    ms=mu_loga.pko.sc.rs[names(mu_loga.pko.sc.rs)%in%ts]
    lines(ms~ts,col=adjustcolor(basin_cols[1],alpha.f=0.1),lwd=3)
  }else
    if(pko.info$ocean.basin[i]=='NC'){
      ms=mu_loga.pko.nc.rs[names(mu_loga.pko.nc.rs)%in%ts]
      lines(ms~ts,col=adjustcolor(basin_cols[2],alpha.f=0.1),lwd=3)
    }else
      if(pko.info$ocean.basin[i]=='GoA'){
        ms=mu_loga.pko.goa.rs[names(mu_loga.pko.goa.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[3],alpha.f=0.1),lwd=3)
      }else{
        ms=mu_loga.pko.bs.rs[names(mu_loga.pko.bs.rs)%in%ts]
        lines(ms~ts,col=adjustcolor(basin_cols[4],alpha.f=0.1),lwd=3)
      }
}
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.25,'N',cex=1)
text(x=par('usr')[2]-2.5,y=par('usr')[4]-0.8,nrow(subset(pko.info,ocean.basin=='BS')),cex=1,col=basin_cols[4])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-1.4,nrow(subset(pko.info,ocean.basin=='GoA')),cex=1,col=basin_cols[3])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.0,nrow(subset(pko.info,ocean.basin=='NC')),cex=1,col=basin_cols[2])
text(x=par('usr')[2]-2.5,y=par('usr')[4]-2.6,nrow(subset(pko.info,ocean.basin=='SC')),cex=1,col=basin_cols[1])

mtext(expression(bold('F')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold('Pink - Odd')),side=3,adj=0.5,line=0.5,cex=0.75)
dev.off()

##fig. 4####
jpeg('./outputs/main figs/figure_4.jpg',width=8,height=9,res=1200,units='in')
par(mfrow=c(3,2),mar=c(5.1,4.1,2.1,2.1))

umsy.diff.chk=umsy_comp_chk$p.change*100

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts.chk=cut(umsy.diff.chk, breaks = c(seq(-100,100,length.out=11)),right=F)
p_cols2.chk=p_cols1[as.numeric(cuts.chk)]


plot(r.umsy~st.umsy,data=umsy_comp_chk,ylab='',xlab='',bty='l',type='n',xlim=c(0,1),ylim=c(0,1),font.axis=2,las=1)
abline(0,1)
for(i in 1:nrow(umsy_comp_chk)){
  lines(rep(r.umsy[i],2)~c(st.umsy.l90[i],st.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_chk,lwd=0.6)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(st.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_chk,lwd=0.6)
  lines(c(st.umsy[i],r.umsy[i])~rep(st.umsy[i],2),col=adjustcolor(p_cols2.chk[i],alpha.f = 0.8),data=umsy_comp_chk,lwd=1)
  
}
points(r.umsy~st.umsy,data=umsy_comp_chk,bg=adjustcolor(p_cols2.chk,alpha.f=0.8),col=adjustcolor('black',alpha.f=0.6),pch=21,cex=1.4,lwd=0.5)
mtext(expression(bold('Chinook')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('A')),side=3,adj=0,line=0.5,cex=0.75)


u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v1 <-  c(v[1]+(v[2]-v[1])*0.05, v[1]+(v[2]-v[1])*0.3, v[3]+(v[4]-v[3])*0.75, v[4]-(v[4]-v[3])*0.01)

umsy.diff.scy=umsy_comp_scy$p.change*100

cuts.scy=cut(umsy.diff.scy, breaks = c(seq(-100,100,length.out=11)),right=F)
p_cols2.scy=p_cols1[as.numeric(cuts.scy)]

plot(r.umsy~st.umsy,data=umsy_comp_scy,ylab='',xlab='',bty='l',type='n',xlim=c(0,1),ylim=c(0,1),font.axis=2,las=1)
abline(0,1)
for(i in 1:nrow(umsy_comp_scy)){
  lines(rep(r.umsy[i],2)~c(st.umsy.l90[i],st.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_scy,lwd=0.6)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(st.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_scy,lwd=0.6)
  lines(c(st.umsy[i],r.umsy[i])~rep(st.umsy[i],2),col=adjustcolor(p_cols2.scy[i],alpha.f = 0.8),data=umsy_comp_scy,lwd=1)
  
}
points(r.umsy~st.umsy,data=umsy_comp_scy,bg=adjustcolor(p_cols2.scy,alpha.f=0.8),col=adjustcolor('black',alpha.f=0.6),pch=21,cex=1.4,lwd=0.5)
mtext(expression(bold('Sockeye')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('B')),side=3,adj=0,line=0.5,cex=0.75)

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v2 <-  c(v[1]+(v[2]-v[1])*0.05, v[1]+(v[2]-v[1])*0.3, v[3]+(v[4]-v[3])*0.75, v[4]-(v[4]-v[3])*0.01)

umsy.diff.chm=umsy_comp_chm$p.change*100

cuts.chm=cut(umsy.diff.chm, breaks = c(seq(-100,100,length.out=11)),right=F)
p_cols2.chm=p_cols1[as.numeric(cuts.chm)]

plot(r.umsy~st.umsy,data=umsy_comp_chm,ylab='',xlab='',bty='l',type='n',xlim=c(0,1),ylim=c(0,1),font=2,las=1)
abline(0,1)
for(i in 1:nrow(umsy_comp_chm)){
  lines(rep(r.umsy[i],2)~c(st.umsy.l90[i],st.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_chm,lwd=0.6)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(st.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_chm,lwd=0.6)
  lines(c(st.umsy[i],r.umsy[i])~rep(st.umsy[i],2),col=adjustcolor(p_cols2.chm[i],alpha.f = 0.8),data=umsy_comp_chm,lwd=1)
  
}
points(r.umsy~st.umsy,data=umsy_comp_chm,bg=adjustcolor(p_cols2.chm,alpha.f=0.8),col=adjustcolor('black',alpha.f=0.6),pch=21,cex=1.4,lwd=0.5)
mtext(expression(bold('Chum')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('C')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold(paste('Last observed 5-years mean ',italic(U['MSY,j,t']),sep=''))),side=2,adj=0.5,line=2.25,cex=0.9)

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v3 <-  c(v[1]+(v[2]-v[1])*0.05, v[1]+(v[2]-v[1])*0.3, v[3]+(v[4]-v[3])*0.75, v[4]-(v[4]-v[3])*0.01)

umsy.diff.cho=umsy_comp_cho$p.change*100

cuts.cho=cut(umsy.diff.cho, breaks = c(seq(-100,100,length.out=11)),right=F)
p_cols2.cho=p_cols1[as.numeric(cuts.cho)]

plot(r.umsy~st.umsy,data=umsy_comp_cho,ylab='',xlab='',bty='l',type='n',xlim=c(0,1),ylim=c(0,1),font=2,las=1)
abline(0,1)
for(i in 1:nrow(umsy_comp_cho)){
  lines(rep(r.umsy[i],2)~c(st.umsy.l90[i],st.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_cho,lwd=0.6)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(st.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_cho,lwd=0.6)
  lines(c(st.umsy[i],r.umsy[i])~rep(st.umsy[i],2),col=adjustcolor(p_cols2.cho[i],alpha.f = 0.8),data=umsy_comp_cho,lwd=1)
  
}
points(r.umsy~st.umsy,data=umsy_comp_cho,bg=adjustcolor(p_cols2.cho,alpha.f=0.8),col=adjustcolor('black',alpha.f=0.6),pch=21,cex=1.4,lwd=0.5)
mtext(expression(bold('Coho')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('D')),side=3,adj=0,line=0.5,cex=0.75)

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v4 <-  c(v[1]+(v[2]-v[1])*0.05, v[1]+(v[2]-v[1])*0.3, v[3]+(v[4]-v[3])*0.75, v[4]-(v[4]-v[3])*0.01)

umsy.diff.pke=umsy_comp_pke$p.change*100

cuts.pke=cut(umsy.diff.pke, breaks = c(seq(-100,100,length.out=11)),right=F)
p_cols2.pke=p_cols1[as.numeric(cuts.pke)]

plot(r.umsy~st.umsy,data=umsy_comp_pke,ylab='',xlab='',bty='l',type='n',xlim=c(0,1),ylim=c(0,1),font=2,las=1)
abline(0,1)
for(i in 1:nrow(umsy_comp_pke)){
  lines(rep(r.umsy[i],2)~c(st.umsy.l90[i],st.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_pke,lwd=0.6)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(st.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_pke,lwd=0.6)
  lines(c(st.umsy[i],r.umsy[i])~rep(st.umsy[i],2),col=adjustcolor(p_cols2.pke[i],alpha.f = 0.8),data=umsy_comp_pke,lwd=1)
  
}
points(r.umsy~st.umsy,data=umsy_comp_pke,bg=adjustcolor(p_cols2.pke,alpha.f=0.8),col=adjustcolor('black',alpha.f=0.6),pch=21,cex=1.4,lwd=0.5)
mtext(expression(bold('Pink - Even')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('E')),side=3,adj=0,line=0.5,cex=0.75)
mtext(expression(bold(paste('Stationary ',italic(U['MSY,j']),sep=''))),side=1,at=1.25,line=3,cex=0.9)

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v5 <-  c(v[1]+(v[2]-v[1])*0.05, v[1]+(v[2]-v[1])*0.3, v[3]+(v[4]-v[3])*0.75, v[4]-(v[4]-v[3])*0.01)

umsy.diff.pko=umsy_comp_pko$p.change*100

cuts.pko=cut(umsy.diff.pko, breaks = c(seq(-100,100,length.out=11)),right=F)
p_cols2.pko=p_cols1[as.numeric(cuts.pko)]

plot(r.umsy~st.umsy,data=umsy_comp_pko,ylab='',xlab='',bty='l',type='n',xlim=c(0,1),ylim=c(0,1),font=2,las=1)
abline(0,1)
for(i in 1:nrow(umsy_comp_pko)){
  lines(rep(r.umsy[i],2)~c(st.umsy.l90[i],st.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_pko,lwd=0.8)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(st.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp_pko,lwd=0.8)
  lines(c(st.umsy[i],r.umsy[i])~rep(st.umsy[i],2),col=adjustcolor(p_cols2.pko[i],alpha.f = 0.8),data=umsy_comp_pko,lwd=2)
  
}
points(r.umsy~st.umsy,data=umsy_comp_pko,bg=adjustcolor(p_cols2.pko,alpha.f=0.8),col=adjustcolor('black',alpha.f=0.6),pch=21,cex=1.4,lwd=0.5)
mtext(expression(bold('Pink - Odd')),side=3,adj=0.5,line=0.5,cex=0.75)
mtext(expression(bold('F')),side=3,adj=0,line=0.5,cex=0.75)

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
v6 <-  c(v[1]+(v[2]-v[1])*0.05, v[1]+(v[2]-v[1])*0.3, v[3]+(v[4]-v[3])*0.75, v[4]-(v[4]-v[3])*0.01)


par(fig=v1, new=TRUE, mar=c(0,0,0,0),xpd=T)

breaks=seq(-100,100,by=10)
break_cols=p_cols1[cut(breaks,c(seq(-100,100,length.out=11)),right=F)]
break_cols=adjustcolor(break_cols,alpha.f=0.8)

hist(umsy.diff.chk,breaks=breaks,xlab='',main='',border='white',col=break_cols,yaxt='n',xaxt='n',ylab='',xlim=c(-100,100),cex.axis=0.5,cex.lab=0.5)
axis(side=1,at=c(seq(-100,100,by=25)),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
axis(side=2,at=seq(0,10,by=5),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
mtext(c('0','5','10'),at=c(0,5,10),cex=0.35,side=2,line=0.3,las=1)
mtext(c('-100','-50','0','+50','+100'),at=c(-100,-50,0,50,100),cex=0.35,side=1,line=0.05)
mtext('% change',side=1,cex=0.45,line=0.8)

par(fig=v2, new=TRUE, mar=c(0,0,0,0),xpd=T)

hist(umsy.diff.scy,breaks=breaks,xlab='',main='',border='white',col=break_cols,yaxt='n',xaxt='n',ylab='',xlim=c(-100,100),cex.axis=0.5,cex.lab=0.5)
axis(side=1,at=c(seq(-100,100,by=25)),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
axis(side=2,at=seq(0,10,by=5),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
mtext(c('0','5','10'),at=c(0,5,10),cex=0.35,side=2,line=0.3,las=1)
mtext(c('-100','-50','0','+50','+100'),at=c(-100,-50,0,50,100),cex=0.35,side=1,line=0.05)
mtext('% change',side=1,cex=0.45,line=0.8)

par(fig=v3, new=TRUE, mar=c(0,0,0,0),xpd=T)

hist(umsy.diff.chm,breaks=breaks,xlab='',main='',border='white',col=break_cols,yaxt='n',xaxt='n',ylab='',xlim=c(-100,100),cex.axis=0.5,cex.lab=0.5)
axis(side=1,at=c(seq(-100,100,by=25)),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
axis(side=2,at=seq(0,5,by=5),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
mtext(c('0','5'),at=c(0,5),cex=0.35,side=2,line=0.3,las=1)
mtext(c('-100','-50','0','+50','+100'),at=c(-100,-50,0,50,100),cex=0.35,side=1,line=0.05)
mtext('% change',side=1,cex=0.45,line=0.8)

par(fig=v4, new=TRUE, mar=c(0,0,0,0),xpd=T)

hist(umsy.diff.cho,breaks=breaks,xlab='',main='',border='white',col=break_cols,yaxt='n',xaxt='n',ylab='',xlim=c(-100,100),cex.axis=0.5,cex.lab=0.5)
axis(side=1,at=c(seq(-100,100,by=25)),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
axis(side=2,at=seq(0,10,by=5),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
mtext(c('0','5','10'),at=c(0,5,10),cex=0.35,side=2,line=0.3,las=1)
mtext(c('-100','-50','0','+50','+100'),at=c(-100,-50,0,50,100),cex=0.35,side=1,line=0.05)
mtext('% change',side=1,cex=0.45,line=0.8)

par(fig=v5, new=TRUE, mar=c(0,0,0,0),xpd=T)

hist(umsy.diff.pke,breaks=breaks,xlab='',main='',border='white',col=break_cols,yaxt='n',xaxt='n',ylab='',xlim=c(-100,100),cex.axis=0.5,cex.lab=0.5)
axis(side=1,at=c(seq(-100,100,by=25)),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
axis(side=2,at=seq(0,5,by=5),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
mtext(c('0','5'),at=c(0,5),cex=0.35,side=2,line=0.3,las=1)
mtext(c('-100','-50','0','+50','+100'),at=c(-100,-50,0,50,100),cex=0.35,side=1,line=0.05)
mtext('% change',side=1,cex=0.45,line=0.8)

par(fig=v6, new=TRUE, mar=c(0,0,0,0),xpd=T)

hist(umsy.diff.pko,breaks=breaks,xlab='',main='',border='white',col=break_cols,yaxt='n',xaxt='n',ylab='',xlim=c(-100,100),cex.axis=0.5,cex.lab=0.5)
axis(side=1,at=c(seq(-100,100,by=25)),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
axis(side=2,at=seq(0,5,by=5),labels=F,cex.lab=0.5,cex.axis=0.5,tcl=-0.2)
mtext(c('0','5'),at=c(0,5),cex=0.35,side=2,line=0.3,las=1)
mtext(c('-100','-50','0','+50','+100'),at=c(-100,-50,0,50,100),cex=0.35,side=1,line=0.05)
mtext('% change',side=1,cex=0.45,line=0.8)

dev.off()
