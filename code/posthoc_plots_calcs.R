#post-hoc visualizations
library(here);library(ggplot2);library(cowplot);library(ggspatial);library(sf)
library(dplyr)
#Long-term static model results####
#static models
stc.chk=readRDS('./outputs/model fits/hier_fit_chinook.RDS')
stc.cho=readRDS('./outputs/model fits/hier_fit_coho.RDS')
stc.chm=readRDS('./outputs/model fits/hier_fit_chum.RDS')
stc.pke=readRDS('./outputs/model fits/hier_fit_pinke.RDS')
stc.pko=readRDS('./outputs/model fits/hier_fit_pinko.RDS')
stc.scy=readRDS('./outputs/model fits/hier_fit_sockeye.RDS')

#visual assets
sp_cols=cbind(c('#6f9ba7',
                '#316024',
                '#295f8b',
                '#a2450c',
                '#811b0a'),c('Chinook','Chum','Coho','Pink','Sockeye'))

basin_cols=cbind(c('#ffc454',
                   '#64a95d',
                   '#007779',
                   '#003f5c')
                 ,c('WC','SEAK','GOA','BS'))

#metadata
stk.info=read.csv(here('data','filtered datasets','stock_info2024-10-17.csv'))
chk.info=subset(stk.info,species=='Chinook')
scy.info=subset(stk.info,species=='Sockeye')
chm.info=subset(stk.info,species=='Chum')
cho.info=subset(stk.info,species=='Coho')
pke.info=subset(stk.info,species=='Pink-Even')
pko.info=subset(stk.info,species=='Pink-Odd')

#broken down by LME
#wc
chk.wc=which(chk.info$ocean.basin=='WC')
chm.wc=which(chm.info$ocean.basin=='WC')
cho.wc=which(cho.info$ocean.basin=='WC')
pke.wc=which(pke.info$ocean.basin=='WC')
pko.wc=which(pko.info$ocean.basin=='WC')
scy.wc=which(scy.info$ocean.basin=='WC')
#nbc-seak
chk.seak=which(chk.info$ocean.basin=='SEAK')
chm.seak=which(chm.info$ocean.basin=='SEAK')
cho.seak=which(cho.info$ocean.basin=='SEAK')
pke.seak=which(pke.info$ocean.basin=='SEAK')
pko.seak=which(pko.info$ocean.basin=='SEAK')
scy.seak=which(scy.info$ocean.basin=='SEAK')
#gulf
chk.goa=which(chk.info$ocean.basin=='GOA')
chm.goa=which(chm.info$ocean.basin=='GOA')
cho.goa=which(cho.info$ocean.basin=='GOA')
pke.goa=which(pke.info$ocean.basin=='GOA')
pko.goa=which(pko.info$ocean.basin=='GOA')
scy.goa=which(scy.info$ocean.basin=='GOA')
#bering
chk.bs=which(chk.info$ocean.basin=='BS')
chm.bs=which(chm.info$ocean.basin=='BS')
cho.bs=which(cho.info$ocean.basin=='BS')
pke.bs=which(pke.info$ocean.basin=='BS')
pko.bs=which(pko.info$ocean.basin=='BS')
scy.bs=which(scy.info$ocean.basin=='BS')




#Productivity####
##mean productivity by species####
glob.alpha.chk=stc.chk$draws(variable='log_a0',format='draws_matrix')
glob.alpha.chm=stc.chm$draws(variable='log_a0',format='draws_matrix')
glob.alpha.cho=stc.cho$draws(variable='log_a0',format='draws_matrix')
glob.alpha.pke=stc.pke$draws(variable='log_a0',format='draws_matrix')
glob.alpha.pko=stc.pko$draws(variable='log_a0',format='draws_matrix')
glob.alpha.scy=stc.scy$draws(variable='log_a0',format='draws_matrix')

par(mfrow=c(3,2))
par(mar=c(4,5,1,1))
min=min(c(glob.alpha.chk,glob.alpha.chm,glob.alpha.cho,glob.alpha.pke,glob.alpha.pko,glob.alpha.scy))
max=max(c(glob.alpha.chk,glob.alpha.chm,glob.alpha.cho,glob.alpha.pke,glob.alpha.pko,glob.alpha.scy))

hist(as.vector(glob.alpha.chk),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),main='Chinook',breaks=50,freq=T,xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(glob.alpha.chk),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(as.vector(glob.alpha.scy),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),main='Sockeye',breaks=50,freq=T,xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(glob.alpha.scy),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(as.vector(glob.alpha.chm),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),main='Chum',breaks=50,freq=T,xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(glob.alpha.chm),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(as.vector(glob.alpha.cho),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[3,1],alpha=0.4),main='Coho',breaks=50,freq=T,xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(glob.alpha.cho),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(as.vector(glob.alpha.pke),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),main='Pink - Even',breaks=50,freq=T,xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(glob.alpha.pke),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(as.vector(glob.alpha.pko),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),main='Pink - Odd',breaks=50,freq=T,xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(glob.alpha.pko),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)


##productivity breakdown - LME####
lme.alpha.chk=stc.chk$draws(variable='log_a_or',format='draws_matrix')
lme.alpha.chm=stc.chm$draws(variable='log_a_or',format='draws_matrix')
lme.alpha.cho=stc.cho$draws(variable='log_a_or',format='draws_matrix')
lme.alpha.pke=stc.pke$draws(variable='log_a_or',format='draws_matrix')
lme.alpha.pko=stc.pko$draws(variable='log_a_or',format='draws_matrix')
lme.alpha.scy=stc.scy$draws(variable='log_a_or',format='draws_matrix')

par(mfrow=c(3,2))
par(mar=c(4,5,3,1))

g4=density(lme.alpha.chk[,4],bw=0.03)
g3=density(lme.alpha.chk[,3],bw=0.03)
g2=density(lme.alpha.chk[,2],bw=0.03)
g1=density(lme.alpha.chk[,1],bw=0.03)

ymax=max(c(g4$y,g3$y,g2$y,g1$y))

plot(g4$y~g4$x,type='l',bty='l',col=basin_cols[4,1],yaxt='n',ylab='',lwd=3,xlab=expression(paste('log(',alpha,')')),main='Chinook',ylim=c(0,ymax),xlim=c(0,3))
lines(g3$y~g3$x,lwd=3,col=basin_cols[3,1])
lines(g2$y~g2$x,lwd=3,col=basin_cols[2,1])
lines(g1$y~g1$x,lwd=3,col=basin_cols[1,1])
abline(v=median(glob.alpha.chk[,4]),col=adjustcolor(basin_cols[4,1],alpha=0.95),lwd=2)
abline(v=median(glob.alpha.chk[,3]),col=adjustcolor(basin_cols[3,1],alpha=0.95),lwd=2)
abline(v=median(glob.alpha.chk[,2]),col=adjustcolor(basin_cols[2,1],alpha=0.95),lwd=2)
abline(v=median(glob.alpha.chk[,1]),col=adjustcolor(basin_cols[1,1],alpha=0.95),lwd=2)
ticksat <- seq(0,3,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,3,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
mtext('Posterior density',side=2,line=0.5)

g4=density(lme.alpha.scy[,4],bw=0.03)
g3=density(lme.alpha.scy[,3],bw=0.03)
g2=density(lme.alpha.scy[,2],bw=0.03)
g1=density(lme.alpha.scy[,1],bw=0.03)

ymax=max(c(g4$y,g3$y,g2$y,g1$y))

plot(g4$y~g4$x,type='l',bty='l',col=basin_cols[4,1],yaxt='n',ylab='',lwd=3,xlab=expression(paste('log(',alpha,')')),main='Sockeye',ylim=c(0,ymax),xlim=c(0,3))
lines(g3$y~g3$x,lwd=3,col=basin_cols[3,1])
lines(g2$y~g2$x,lwd=3,col=basin_cols[2,1])
lines(g1$y~g1$x,lwd=3,col=basin_cols[1,1])
abline(v=median(lme.alpha.scy[,4]),col=adjustcolor(basin_cols[4,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.scy[,3]),col=adjustcolor(basin_cols[3,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.scy[,2]),col=adjustcolor(basin_cols[2,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.scy[,1]),col=adjustcolor(basin_cols[1,1],alpha=0.95),lwd=2)
ticksat <- seq(0,3,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,3,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
mtext('Posterior density',side=2,line=0.5)

g4=density(lme.alpha.chm[,4],bw=0.03)
g3=density(lme.alpha.chm[,3],bw=0.03)
g2=density(lme.alpha.chm[,2],bw=0.03)
g1=density(lme.alpha.chm[,1],bw=0.03)

ymax=max(c(g4$y,g3$y,g2$y,g1$y))

plot(g4$y~g4$x,type='l',bty='l',col=basin_cols[4,1],yaxt='n',ylab='',lwd=3,xlab=expression(paste('log(',alpha,')')),main='Chum',ylim=c(0,ymax),xlim=c(0,3))
lines(g3$y~g3$x,lwd=3,col=basin_cols[3,1])
lines(g2$y~g2$x,lwd=3,col=basin_cols[2,1])
lines(g1$y~g1$x,lwd=3,col=basin_cols[1,1])
abline(v=median(lme.alpha.chm[,4]),col=adjustcolor(basin_cols[4,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.chm[,2]),col=adjustcolor(basin_cols[2,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.chm[,1]),col=adjustcolor(basin_cols[1,1],alpha=0.95),lwd=2)
ticksat <- seq(0,3,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,3,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
mtext('Posterior density',side=2,line=0.5)
text(x=2.6,y=par('usr')[4]-0.25,'Bering Sea',col=basin_cols[1,1],cex=2)
text(x=2.6,y=par('usr')[4]-0.75,'Gulf of Alaska',col=basin_cols[2,1],cex=2)
text(x=2.6,y=par('usr')[4]-1.25,'SE Alaska',col=basin_cols[3,1],cex=2)
text(x=2.6,y=par('usr')[4]-1.75,'West Coast',col=basin_cols[4,1],cex=2)


g3=density(lme.alpha.cho[,3],bw=0.03)
g2=density(lme.alpha.cho[,2],bw=0.03)
g1=density(lme.alpha.cho[,1],bw=0.03)

ymax=max(c(g3$y,g2$y,g1$y))

plot(g3$y~g3$x,type='l',bty='l',col=basin_cols[4,1],yaxt='n',ylab='',lwd=3,xlab=expression(paste('log(',alpha,')')),main='Coho',ylim=c(0,ymax),xlim=c(0,3))
lines(g2$y~g2$x,lwd=3,col=basin_cols[3,1])
lines(g1$y~g1$x,lwd=3,col=basin_cols[2,1])
abline(v=median(lme.alpha.cho[,3]),col=adjustcolor(basin_cols[4,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.cho[,2]),col=adjustcolor(basin_cols[3,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.cho[,1]),col=adjustcolor(basin_cols[2,1],alpha=0.95),lwd=2)
ticksat <- seq(0,3,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,3,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
mtext('Posterior density',side=2,line=0.5)

g4=density(lme.alpha.pke[,4],bw=0.03)
g3=density(lme.alpha.pke[,3],bw=0.03)
g2=density(lme.alpha.pke[,2],bw=0.03)
g1=density(lme.alpha.pke[,1],bw=0.03)

ymax=max(c(g4$y,g3$y,g2$y,g1$y))

plot(g4$y~g4$x,type='l',bty='l',col=basin_cols[4,1],yaxt='n',ylab='',lwd=3,xlab=expression(paste('log(',alpha,')')),main='Pink (even)',ylim=c(0,ymax),xlim=c(0,3))
lines(g3$y~g3$x,lwd=3,col=basin_cols[3,1])
lines(g2$y~g2$x,lwd=3,col=basin_cols[2,1])
lines(g1$y~g1$x,lwd=3,col=basin_cols[1,1])
abline(v=median(lme.alpha.pke[,4]),col=adjustcolor(basin_cols[4,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.pke[,3]),col=adjustcolor(basin_cols[3,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.pke[,2]),col=adjustcolor(basin_cols[2,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.pke[,1]),col=adjustcolor(basin_cols[1,1],alpha=0.95),lwd=2)
ticksat <- seq(0,3,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,3,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
mtext('Posterior density',side=2,line=0.5)

g4=density(lme.alpha.pko[,4],bw=0.03)
g3=density(lme.alpha.pko[,3],bw=0.03)
g2=density(lme.alpha.pko[,2],bw=0.03)
g1=density(lme.alpha.pko[,1],bw=0.03)

ymax=max(c(g4$y,g3$y,g2$y,g1$y))

plot(g4$y~g4$x,type='l',bty='l',col=basin_cols[4,1],yaxt='n',ylab='',lwd=3,xlab=expression(paste('log(',alpha,')')),main='Pink (odd)',ylim=c(0,ymax),xlim=c(0,3))
lines(g3$y~g3$x,lwd=3,col=basin_cols[3,1])
lines(g2$y~g2$x,lwd=3,col=basin_cols[2,1])
lines(g1$y~g1$x,lwd=3,col=basin_cols[1,1])
abline(v=median(lme.alpha.pko[,4]),col=adjustcolor(basin_cols[4,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.pko[,3]),col=adjustcolor(basin_cols[3,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.pko[,2]),col=adjustcolor(basin_cols[2,1],alpha=0.95),lwd=2)
abline(v=median(lme.alpha.pko[,1]),col=adjustcolor(basin_cols[1,1],alpha=0.95),lwd=2)
ticksat <- seq(0,3,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,3,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
mtext('Posterior density',side=2,line=0.5)




##productivity breakdown - subregions####

sr.alpha.chk=stc.chk$draws(variable='log_a_sr',format='draws_matrix')
sr.alpha.chm=stc.chm$draws(variable='log_a_sr',format='draws_matrix')
sr.alpha.cho=stc.cho$draws(variable='log_a_sr',format='draws_matrix')
sr.alpha.pke=stc.pke$draws(variable='log_a_sr',format='draws_matrix')
sr.alpha.pko=stc.pko$draws(variable='log_a_sr',format='draws_matrix')
sr.alpha.scy=stc.scy$draws(variable='log_a_sr',format='draws_matrix')

###chinook####
srs.chk=levels(factor(chk.info$region))
length(srs.chk)

min=0
max=3

par(mfrow=c(5,3))
par(mar=c(4,5,3,1))

hist(sr.alpha.chk[,1],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[1],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,1]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chk[,2],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[2],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,2]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chk[,3],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[3],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,3]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chk[,4],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[4],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,4]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chk[,5],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[5],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,5]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chk[,6],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[6],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,6]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(sr.alpha.chk[,7],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[7],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,7]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chk[,8],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[8],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,8]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chk[,9],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[9],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,9]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chk[,10],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[10],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,10]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chk[,11],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[11],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,11]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chk[,12],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[12],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,12]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chk[,13],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main=srs.chk[13],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chk[,13]),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

###sockeye####
srs.scy=levels(factor(scy.info$region))
length(srs.scy)

min=0
max=3

par(mfrow=c(5,3))
par(mar=c(4,5,3,1))

hist(sr.alpha.scy[,1],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[1],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,1]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.scy[,2],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[2],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,2]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.scy[,3],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[3],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,3]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.scy[,4],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[4],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,4]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.scy[,5],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[5],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,5]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.scy[,6],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[6],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,6]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(sr.alpha.scy[,7],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[7],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,7]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.scy[,8],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[8],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,8]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.scy[,9],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[9],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,9]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.scy[,10],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[10],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,10]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.scy[,11],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[11],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,11]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.scy[,12],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[12],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,12]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.scy[,13],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[13],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,13]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.scy[,14],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main=srs.scy[14],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.scy[,14]),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

###chum####
srs.chm=levels(factor(chm.info$region))
length(srs.chm)

min=0
max=3

par(mfrow=c(5,3))
par(mar=c(4,5,3,1))

hist(sr.alpha.chm[,1],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[1],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,1]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,2],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[2],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,2]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,3],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[3],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,3]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,4],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[4],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,4]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,5],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[5],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,5]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,6],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[6],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,6]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(sr.alpha.chm[,7],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[7],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,7]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,8],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[8],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,8]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,9],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[9],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,9]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,10],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[10],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,10]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,11],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[11],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,11]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,12],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[12],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,12]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,13],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[13],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,13]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,14],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[14],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,14]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.chm[,15],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main=srs.chm[15],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.chm[,15]),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

###coho####
srs.cho=levels(factor(cho.info$region))
length(srs.cho)

min=0
max=3

par(mfrow=c(5,3))
par(mar=c(4,5,3,1))

hist(sr.alpha.cho[,1],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[3,1],alpha=0.4),breaks=50,freq=T,main=srs.cho[1],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.cho[,1]),col=adjustcolor(sp_cols[3,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.cho[,2],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[3,1],alpha=0.4),breaks=50,freq=T,main=srs.cho[2],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.cho[,2]),col=adjustcolor(sp_cols[3,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.cho[,3],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[3,1],alpha=0.4),breaks=50,freq=T,main=srs.cho[3],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.cho[,3]),col=adjustcolor(sp_cols[3,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.cho[,4],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[3,1],alpha=0.4),breaks=50,freq=T,main=srs.cho[4],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.cho[,4]),col=adjustcolor(sp_cols[3,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.cho[,5],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[3,1],alpha=0.4),breaks=50,freq=T,main=srs.cho[5],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.cho[,5]),col=adjustcolor(sp_cols[3,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.cho[,6],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[3,1],alpha=0.4),breaks=50,freq=T,main=srs.cho[6],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.cho[,6]),col=adjustcolor(sp_cols[3,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(sr.alpha.cho[,7],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[3,1],alpha=0.4),breaks=50,freq=T,main=srs.cho[7],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.cho[,7]),col=adjustcolor(sp_cols[3,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

###Pink - Even####
srs.pke=levels(factor(pke.info$region))
length(srs.pke)

min=0
max=3

par(mfrow=c(5,3))
par(mar=c(4,5,3,1))

hist(sr.alpha.pke[,1],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[1],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,1]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pke[,2],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[2],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,2]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pke[,3],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[3],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,3]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pke[,4],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[4],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,4]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pke[,5],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[5],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,5]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pke[,6],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[6],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,6]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(sr.alpha.pke[,7],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[7],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,7]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pke[,8],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[8],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,8]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pke[,9],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[9],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,9]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pke[,10],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[10],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,10]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pke[,11],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[11],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,11]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pke[,12],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[12],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,12]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pke[,13],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pke[13],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pke[,13]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

###Pink - Odd####
srs.pko=levels(factor(pko.info$region))
length(srs.pko)

min=0
max=3

par(mfrow=c(5,3))
par(mar=c(4,5,3,1))

hist(sr.alpha.pko[,1],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[1],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,1]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pko[,2],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[2],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,2]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pko[,3],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[3],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,3]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pko[,4],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[4],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,4]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pko[,5],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[5],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,5]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pko[,6],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[6],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,6]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(sr.alpha.pko[,7],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[7],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,7]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pko[,8],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[8],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,8]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pko[,9],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[9],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,9]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pko[,10],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[10],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,10]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pko[,11],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[11],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,11]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pko[,12],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[12],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,12]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(sr.alpha.pko[,13],border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main=srs.pko[13],xlab=expression(paste('log(',alpha,')')),yaxt='n',ylab='')
abline(v=median(sr.alpha.pko[,13]),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)


#Capacity####
#hyperdistributions by species
smax.chk=stc.chk$draws(variable='Smax',format='draws_matrix')
smax.chm=stc.chm$draws(variable='Smax',format='draws_matrix')
smax.cho=stc.cho$draws(variable='Smax',format='draws_matrix')
smax.pke=stc.pke$draws(variable='Smax',format='draws_matrix')
smax.pko=stc.pko$draws(variable='Smax',format='draws_matrix')
smax.scy=stc.scy$draws(variable='Smax',format='draws_matrix')

par(mfrow=c(3,2))
par(mar=c(3,5,3,1))
hist(as.vector(log10(smax.chk)),border=adjustcolor('white',alpha=0.4),xlim=c(0,8),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main='Chinook',xlab='',yaxt='n',ylab='',xaxt='n')
pow <- 0:7
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(0,8,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(1),expression(10),expression(100),expression('1K'),expression('10K'),expression('100K'),expression('1M'),expression('10M'),expression('100M')))

hist(as.vector(log10(smax.scy)),border=adjustcolor('white',alpha=0.4),xlim=c(0,8),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main='Sockeye',xlab='',yaxt='n',ylab='',xaxt='n')
pow <- 0:7
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(0,8,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(1),expression(10),expression(100),expression('1K'),expression('10K'),expression('100K'),expression('1M'),expression('10M'),expression('100M')))

hist(as.vector(log10(smax.chm)),border=adjustcolor('white',alpha=0.4),xlim=c(0,max(log10(smax.chm))),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main='Chum',xlab='',yaxt='n',ylab='',xaxt='n')
pow <- 0:7
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(0,8,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(1),expression(10),expression(100),expression('1K'),expression('10K'),expression('100K'),expression('1M'),expression('10M'),expression('100M')))

hist(as.vector(log10(smax.cho)),border=adjustcolor('white',alpha=0.4),xlim=c(0,max(log10(smax.cho))),col=adjustcolor(sp_cols[3,1],alpha=0.4),breaks=50,freq=T,main='Coho',xlab='',yaxt='n',ylab='',xaxt='n')
pow <- 0:6
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(0,8,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(1),expression(10),expression(100),expression('1K'),expression('10K'),expression('100K'),expression('1M'),expression('10M'),expression('100M')))
par(mar=c(5,5,3,1))
hist(as.vector(log10(smax.pke)),border=adjustcolor('white',alpha=0.4),xlim=c(0,max(log10(smax.pke))),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main='Pink-Even',xlab='Capacity (Smax)',yaxt='n',ylab='',xaxt='n',cex.lab=1.5)
pow <- 0:8
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(0,8,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(1),expression(10),expression(100),expression('1K'),expression('10K'),expression('100K'),expression('1M'),expression('10M'),expression('100M')))

hist(as.vector(log10(smax.pko)),border=adjustcolor('white',alpha=0.4),xlim=c(0,max(log10(smax.pko))),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main='Pink-Odd',xlab='Capacity (Smax)',yaxt='n',ylab='',xaxt='n',cex.lab=1.5)
pow <- 0:8
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(0,8,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(1),expression(10),expression(100),expression('1K'),expression('10K'),expression('100K'),expression('1M'),expression('10M'),expression('100M')))

#Sigma####
glob.sigma.chk=stc.chk$draws(variable='sigma',format='draws_matrix')
glob.sigma.chm=stc.chm$draws(variable='sigma',format='draws_matrix')
glob.sigma.cho=stc.cho$draws(variable='sigma',format='draws_matrix')
glob.sigma.pke=stc.pke$draws(variable='sigma',format='draws_matrix')
glob.sigma.pko=stc.pko$draws(variable='sigma',format='draws_matrix')
glob.sigma.scy=stc.scy$draws(variable='sigma',format='draws_matrix')

par(mfrow=c(3,2))
par(mar=c(4,5,1,1))
min=0
max=4

hist(as.vector(glob.sigma.chk),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main='',xlab=expression(sigma),yaxt='n',ylab='')
abline(v=median(glob.sigma.chk),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(as.vector(glob.sigma.scy),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main='',xlab=expression(sigma),yaxt='n',ylab='')
abline(v=median(glob.sigma.scy),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(as.vector(glob.sigma.chm),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main='',xlab=expression(sigma),yaxt='n',ylab='')
abline(v=median(glob.sigma.chm),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(as.vector(glob.sigma.cho),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[3,1],alpha=0.4),breaks=50,freq=T,main='',xlab=expression(sigma),yaxt='n',ylab='')
abline(v=median(glob.sigma.cho),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(as.vector(glob.sigma.pke),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main='',xlab=expression(sigma),yaxt='n',ylab='')
abline(v=median(glob.sigma.pke),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(as.vector(glob.sigma.pko),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main='',xlab=expression(sigma),yaxt='n',ylab='')
abline(v=median(glob.sigma.pko),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)


#Autocorrelation (rho)####
glob.rho.chk=stc.chk$draws(variable='rho',format='draws_matrix')
glob.rho.chm=stc.chm$draws(variable='rho',format='draws_matrix')
glob.rho.cho=stc.cho$draws(variable='rho',format='draws_matrix')
glob.rho.pke=stc.pke$draws(variable='rho',format='draws_matrix')
glob.rho.pko=stc.pko$draws(variable='rho',format='draws_matrix')
glob.rho.scy=stc.scy$draws(variable='rho',format='draws_matrix')

par(mfrow=c(3,2))
par(mar=c(4,5,1,1))
min=-1
max=1

hist(as.vector(glob.rho.chk),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main='Chinook',xlab=expression(rho),yaxt='n',ylab='')
abline(v=median(glob.rho.chk),col=adjustcolor(sp_cols[1,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(as.vector(glob.rho.scy),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main='Sockeye',xlab=expression(rho),yaxt='n',ylab='')
abline(v=median(glob.rho.scy),col=adjustcolor(sp_cols[5,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)
hist(as.vector(glob.rho.chm),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main='Chum',xlab=expression(rho),yaxt='n',ylab='')
abline(v=median(glob.rho.chm),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(as.vector(glob.rho.cho),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[3,1],alpha=0.4),breaks=50,freq=T,main='Coho',xlab=expression(rho),yaxt='n',ylab='')
abline(v=median(glob.rho.cho),col=adjustcolor(sp_cols[2,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(as.vector(glob.rho.pke),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main='Pink-Even',xlab=expression(rho),yaxt='n',ylab='')
abline(v=median(glob.rho.pke),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

hist(as.vector(glob.rho.pko),border=adjustcolor('white',alpha=0.4),xlim=c(min,max),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main='Pink-Odd',xlab=expression(rho),yaxt='n',ylab='')
abline(v=median(glob.rho.pko),col=adjustcolor(sp_cols[4,1],alpha=0.8))
ticksat <- seq(0,4,by=0.1)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
ticksat <- seq(0,4,by=0.5)
axis(1, ticksat, col="black", labels=NA,
     tcl=-0.5, lwd=0, lwd.ticks=1)

#Smsy####
#hyperdistributions by species
smsy.chk=stc.chk$draws(variable='Smsy',format='draws_matrix')
smsy.chm=stc.chm$draws(variable='Smsy',format='draws_matrix')
smsy.cho=stc.cho$draws(variable='Smsy',format='draws_matrix')
smsy.pke=stc.pke$draws(variable='Smsy',format='draws_matrix')
smsy.pko=stc.pko$draws(variable='Smsy',format='draws_matrix')
smsy.scy=stc.scy$draws(variable='Smsy',format='draws_matrix')

par(mfrow=c(3,2))
par(mar=c(4,5,3,1))

hist(as.vector(log10(smsy.chk)),border=adjustcolor('white',alpha=0.4),xlim=c(0,8),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main='Chinook',xlab='Smsy',yaxt='n',ylab='',xaxt='n')
pow <- 0:7
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(0,8,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(1),expression(10),expression(100),expression('1K'),expression('10K'),expression('100K'),expression('1M'),expression('10M'),expression('100M')))

hist(as.vector(log10(smsy.scy)),border=adjustcolor('white',alpha=0.4),xlim=c(0,8),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main='Sockeye',xlab='Smsy',yaxt='n',ylab='',xaxt='n')
pow <- 0:7
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(0,8,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(1),expression(10),expression(100),expression('1K'),expression('10K'),expression('100K'),expression('1M'),expression('10M'),expression('100M')))

hist(as.vector(log10(smsy.chm)),border=adjustcolor('white',alpha=0.4),xlim=c(0,8),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main='Chum',xlab='Smsy',yaxt='n',ylab='',xaxt='n')
pow <- 0:7
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(0,8,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(1),expression(10),expression(100),expression('1K'),expression('10K'),expression('100K'),expression('1M'),expression('10M'),expression('100M')))

hist(as.vector(log10(smsy.cho)),border=adjustcolor('white',alpha=0.4),xlim=c(0,8),col=adjustcolor(sp_cols[3,1],alpha=0.4),breaks=50,freq=T,main='Coho',xlab='Smsy',yaxt='n',ylab='',xaxt='n')
pow <- 0:7
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(0,8,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(1),expression(10),expression(100),expression('1K'),expression('10K'),expression('100K'),expression('1M'),expression('10M'),expression('100M')))

hist(as.vector(log10(smsy.pke)),border=adjustcolor('white',alpha=0.4),xlim=c(0,8),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main='Pink - Even',xlab='Smsy',yaxt='n',ylab='',xaxt='n')
pow <- 0:7
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(0,8,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(1),expression(10),expression(100),expression('1K'),expression('10K'),expression('100K'),expression('1M'),expression('10M'),expression('100M')))

hist(as.vector(log10(smsy.pko)),border=adjustcolor('white',alpha=0.4),xlim=c(0,8),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main='Pink - Odd',xlab='Smsy',yaxt='n',ylab='',xaxt='n')
pow <- 0:7
ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
axis(1, log10(ticksat), col="black", labels=NA,
     tcl=-0.2, lwd=0, lwd.ticks=1)
axis(1, col="black", at=seq(0,8,by=1),   tcl=-0.45, cex.axis=1.2,
     labels=c(expression(1),expression(10),expression(100),expression('1K'),expression('10K'),expression('100K'),expression('1M'),expression('10M'),expression('100M')))

#hyperdistributions by species & LME

#Umsy####
umsy.chk=stc.chk$draws(variable='Umsy',format='draws_matrix')
umsy.chm=stc.chm$draws(variable='Umsy',format='draws_matrix')
umsy.cho=stc.cho$draws(variable='Umsy',format='draws_matrix')
umsy.pke=stc.pke$draws(variable='Umsy',format='draws_matrix')
umsy.pko=stc.pko$draws(variable='Umsy',format='draws_matrix')
umsy.scy=stc.scy$draws(variable='Umsy',format='draws_matrix')

par(mfrow=c(3,2))
par(mar=c(4,5,3,1))

hist(as.vector(umsy.chk),border=adjustcolor('white',alpha=0.4),xlim=c(0,1),col=adjustcolor(sp_cols[1,1],alpha=0.4),breaks=50,freq=T,main='Chinook',xlab='',yaxt='n',ylab='',cex.axis=1.25)

hist(as.vector(umsy.scy),border=adjustcolor('white',alpha=0.4),xlim=c(0,1),col=adjustcolor(sp_cols[5,1],alpha=0.4),breaks=50,freq=T,main='Sockeye',xlab='',yaxt='n',ylab='',cex.axis=1.25)

hist(as.vector(umsy.chm),border=adjustcolor('white',alpha=0.4),xlim=c(0,1),col=adjustcolor(sp_cols[2,1],alpha=0.4),breaks=50,freq=T,main='Chum',xlab='',yaxt='n',ylab='',cex.axis=1.25)

hist(as.vector(umsy.cho),border=adjustcolor('white',alpha=0.4),xlim=c(0,1),col=adjustcolor(sp_cols[3,1],alpha=0.4),breaks=50,freq=T,main='Coho',xlab='',yaxt='n',ylab='',cex.axis=1.25)

hist(as.vector(umsy.pke),border=adjustcolor('white',alpha=0.4),xlim=c(0,1),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main='Pink - Even',xlab='Umsy',yaxt='n',ylab='',cex.axis=1.25,cex.lab=1.5)
hist(as.vector(umsy.pko),border=adjustcolor('white',alpha=0.4),xlim=c(0,1),col=adjustcolor(sp_cols[4,1],alpha=0.4),breaks=50,freq=T,main='Pink - Odd',xlab='Umsy',yaxt='n',ylab='',cex.axis=1.25,cex.lab=1.5)

#hyperdistributions by species & LME


#comparative patterns####
par(mfrow=c(3,2))
par(mar=c(4,5,3,1))

##productivity vs. capacity####
alphas.chk=stc.chk$draws(variables='log_a',format='draws_matrix')
alphas.scy=stc.scy$draws(variables='log_a',format='draws_matrix')
alphas.chm=stc.chm$draws(variables='log_a',format='draws_matrix')
alphas.cho=stc.cho$draws(variables='log_a',format='draws_matrix')
alphas.pke=stc.pke$draws(variables='log_a',format='draws_matrix')
alphas.pko=stc.pko$draws(variables='log_a',format='draws_matrix')

corr.across.chk=NA
for(i in 1:nrow(alphas.chk)){
  corr.across.chk[i]=cor(as.vector(alphas.chk[i,]),as.vector(log(smax.chk[i,])))
}
corr.within.chk=NA
for(i in 1:ncol(alphas.chk)){
  corr.within.chk[i]=cor(as.vector(alphas.chk[,i]),as.vector(log(smax.chk[,i])))
}

corr.across.scy=NA
for(i in 1:nrow(alphas.scy)){
  corr.across.scy[i]=cor(as.vector(alphas.scy[i,]),as.vector(log(smax.scy[i,])))
}
corr.within.scy=NA
for(i in 1:ncol(alphas.scy)){
  corr.within.scy[i]=cor(as.vector(alphas.scy[,i]),as.vector(log(smax.scy[,i])))
}

corr.across.chm=NA
for(i in 1:nrow(alphas.chm)){
  corr.across.chm[i]=cor(as.vector(alphas.chm[i,]),as.vector(log(smax.chm[i,])))
}
corr.within.chm=NA
for(i in 1:ncol(alphas.chm)){
  corr.within.chm[i]=cor(as.vector(alphas.chm[,i]),as.vector(log(smax.chm[,i])))
}

corr.across.cho=NA
for(i in 1:nrow(alphas.cho)){
  corr.across.cho[i]=cor(as.vector(alphas.cho[i,]),as.vector(log(smax.cho[i,])))
}
corr.within.cho=NA
for(i in 1:ncol(alphas.cho)){
  corr.within.cho[i]=cor(as.vector(alphas.cho[,i]),as.vector(log(smax.cho[,i])))
}

corr.across.pke=NA
for(i in 1:nrow(alphas.pke)){
  corr.across.pke[i]=cor(as.vector(alphas.pke[i,]),as.vector(log(smax.pke[i,])))
}
corr.within.pke=NA
for(i in 1:ncol(alphas.pke)){
  corr.within.pke[i]=cor(as.vector(alphas.pke[,i]),as.vector(log(smax.pke[,i])))
}

corr.across.pko=NA
for(i in 1:nrow(alphas.pko)){
  corr.across.pko[i]=cor(as.vector(alphas.pko[i,]),as.vector(log(smax.pko[i,])))
}
corr.within.pko=NA
for(i in 1:ncol(alphas.pko)){
  corr.within.pko[i]=cor(as.vector(alphas.pko[,i]),as.vector(log(smax.pko[,i])))
}

par(mfrow=c(3,2))
par(mar=c(3,5,3,1))
#Correlations between alpha/beta across stocks
alphas.chk=stc.chk$draws(variables='log_a',format='draws_matrix')
plot(apply(alphas.chk,2,median)~apply(log10(smax.chk),2,median),bty='l',ylab=expression(paste('Productivity - log(',alpha,')')),xlab='',pch=21,col=sp_cols[1,1],cex=1.25,cex.lab=1.25,main='Chinook')
text(x=par('usr')[2]*0.9,y=par('usr')[4]*0.9,paste('r = ',round(median(corr.across.chk),2)),cex=1.2)

alphas.scy=stc.scy$draws(variables='log_a',format='draws_matrix')
plot(apply(alphas.scy,2,median)~apply(log10(smax.scy),2,median),bty='l',ylab='',xlab='',pch=21,col=sp_cols[5,1],cex=1.2,main='Sockeye')
text(x=par('usr')[2]*0.9,y=par('usr')[4]*0.9,paste('r = ',round(median(corr.across.scy),2)),cex=1.2)

alphas.chm=stc.chm$draws(variables='log_a',format='draws_matrix')
plot(apply(alphas.chm,2,median)~apply(log10(smax.chm),2,median),bty='l',ylab=expression(paste('Productivity - log(',alpha,')')),xlab='',pch=21,col=sp_cols[2,1],cex=1.25,cex.lab=1.25,main='Chum')
text(x=par('usr')[2]*0.9,y=par('usr')[4]*0.9,paste('r = ',round(median(corr.across.chm),2)),cex=1.2)

alphas.cho=stc.cho$draws(variables='log_a',format='draws_matrix')
plot(apply(alphas.cho,2,median)~apply(log10(smax.cho),2,median),bty='l',ylab='',xlab='',pch=21,col=sp_cols[3,1],cex=1.25,cex.lab=1.25,main='Coho')
text(x=par('usr')[2]*0.9,y=par('usr')[4]*0.9,paste('r = ',round(median(corr.across.cho),2)),cex=1.2)
par(mar=c(5,5,3,1))
alphas.pke=stc.pke$draws(variables='log_a',format='draws_matrix')
plot(apply(alphas.pke,2,median)~apply(log10(smax.pke),2,median),bty='l',ylab=expression(paste('Productivity - log(',alpha,')')),xlab='Capacity (Smax)',pch=21,col=sp_cols[4,1],cex=1.25,cex.lab=1.25,main='Pink-Even')
text(x=par('usr')[2]*0.9,y=par('usr')[4]*0.9,paste('r = ',round(median(corr.across.pke),2)),cex=1.2)

alphas.pko=stc.pko$draws(variables='log_a',format='draws_matrix')
plot(apply(alphas.pko,2,median)~apply(log10(smax.pko),2,median),bty='l',ylab='',xlab='Capacity (Smax)',pch=21,col=sp_cols[4,1],cex=1.25,cex.lab=1.25,main='Pink-Odd')
text(x=par('usr')[2]*0.9,y=par('usr')[4]*0.9,paste('r = ',round(median(corr.across.pko),2)),cex=1.2)

#among stock correlations
par(mfrow=c(3,2))
par(mar=c(4,5,3,1))
#productivity vs. capacity
hist(corr.across.chk,breaks=30,xlab='Correlation estimate',yaxt='n',ylab='',main='Chinook',border=adjustcolor('white',alpha=0.4),col=adjustcolor(sp_cols[1,1],alpha=0.4))
abline(v=median(corr.across.chk),col=adjustcolor(sp_cols[1,1],alpha=0.8))
hist(corr.across.scy,breaks=30,xlab='Correlation estimate',yaxt='n',ylab='',main='Sockeye',border=adjustcolor('white',alpha=0.4),col=adjustcolor(sp_cols[5,1],alpha=0.4))
abline(v=median(corr.across.scy),col=adjustcolor(sp_cols[5,1],alpha=0.8))

hist(corr.across.chm,breaks=30,xlab='Correlation estimate',yaxt='n',ylab='',main='Chum',border=adjustcolor('white',alpha=0.4),col=adjustcolor(sp_cols[2,1],alpha=0.4))
abline(v=median(corr.across.chm),col=adjustcolor(sp_cols[2,1],alpha=0.8))

hist(corr.across.cho,breaks=30,xlab='Correlation estimate',yaxt='n',ylab='',main='Coho',border=adjustcolor('white',alpha=0.4),col=adjustcolor(sp_cols[3,1],alpha=0.4))
abline(v=median(corr.across.cho),col=adjustcolor(sp_cols[3,1],alpha=0.8))

hist(corr.across.pke,breaks=30,xlab='Correlation estimate',yaxt='n',ylab='',main='Pink - even',border=adjustcolor('white',alpha=0.4),col=adjustcolor(sp_cols[4,1],alpha=0.4))
abline(v=median(corr.across.pke),col=adjustcolor(sp_cols[4,1],alpha=0.8))

hist(corr.across.pko,breaks=30,xlab='Correlation estimate',yaxt='n',ylab='',main='Pink - odd',border=adjustcolor('white',alpha=0.4),col=adjustcolor(sp_cols[4,1],alpha=0.4))
abline(v=median(corr.across.pko),col=adjustcolor(sp_cols[4,1],alpha=0.8))

#within-stock posterior correlations
par(mfrow=c(3,2))
par(mar=c(4,5,3,1))
#productivity vs. capacity
hist(corr.within.chk,breaks=30,xlab='Correlation estimate',yaxt='n',ylab='',main='Chinook',border=adjustcolor('white',alpha=0.4),col=adjustcolor(sp_cols[1,1],alpha=0.4))
abline(v=median(corr.within.chk),col=adjustcolor(sp_cols[1,1],alpha=0.8))
hist(corr.within.scy,breaks=30,xlab='Correlation estimate',yaxt='n',ylab='',main='Sockeye',border=adjustcolor('white',alpha=0.4),col=adjustcolor(sp_cols[5,1],alpha=0.4))
abline(v=median(corr.within.scy),col=adjustcolor(sp_cols[5,1],alpha=0.8))

hist(corr.within.chm,breaks=30,xlab='Correlation estimate',yaxt='n',ylab='',main='Chum',border=adjustcolor('white',alpha=0.4),col=adjustcolor(sp_cols[2,1],alpha=0.4))
abline(v=median(corr.within.chm),col=adjustcolor(sp_cols[2,1],alpha=0.8))

hist(corr.within.cho,breaks=30,xlab='Correlation estimate',yaxt='n',ylab='',main='Coho',border=adjustcolor('white',alpha=0.4),col=adjustcolor(sp_cols[3,1],alpha=0.4))
abline(v=median(corr.within.cho),col=adjustcolor(sp_cols[3,1],alpha=0.8))

hist(corr.within.pke,breaks=30,xlab='Correlation estimate',yaxt='n',ylab='',main='Pink - even',border=adjustcolor('white',alpha=0.4),col=adjustcolor(sp_cols[4,1],alpha=0.4))
abline(v=median(corr.within.pke),col=adjustcolor(sp_cols[4,1],alpha=0.8))

hist(corr.within.pko,breaks=30,xlab='Correlation estimate',yaxt='n',ylab='',main='Pink - odd',border=adjustcolor('white',alpha=0.4),col=adjustcolor(sp_cols[4,1],alpha=0.4))
abline(v=median(corr.within.pko),col=adjustcolor(sp_cols[4,1],alpha=0.8))


#Time-varying fits####
tv.chk=readRDS('./outputs/model fits/tv_fit_chinook.RDS')
tv.cho=readRDS('./outputs/model fits/tv_fit_coho.RDS')
tv.chm=readRDS('./outputs/model fits/tv_fit_chum.RDS')
tv.pke=readRDS('./outputs/model fits/tv_fit_pinke.RDS')
tv.pko=readRDS('./outputs/model fits/tv_fit_pinko.RDS')
tv.scy=readRDS('./outputs/model fits/tv_fit_sockeye.RDS')

#Productivity trajectories####

## Global estimates ####
mu_loga.chk=tv.chk$draws(variables='mu_log_a',format='draws_matrix')
mu_loga.cho=tv.cho$draws(variables='mu_log_a',format='draws_matrix')
mu_loga.chm=tv.chm$draws(variables='mu_log_a',format='draws_matrix')
mu_loga.pke=tv.pke$draws(variables='mu_log_a',format='draws_matrix')
mu_loga.pko=tv.pko$draws(variables='mu_log_a',format='draws_matrix')
mu_loga.scy=tv.scy$draws(variables='mu_log_a',format='draws_matrix')

par(mar=c(5,5,1,1))
ylim=c(2,7)
plot(c(ylim[1],ylim[2])~c(min(stk.info$begin),max(stk.info$end)),bty='l',type='n',ylab=expression(paste("Mean Productivity, max. R/S")),xlab='Brood year',cex.axis=1.3,cex.lab=1.3,ylim=c(2,7))
lines(apply(exp(mu_loga.chk),2,median)~seq(min(chk.info$begin),max(chk.info$end)),lwd=3,col=adjustcolor(sp_cols[1],alpha.f=0.9))
polygon(c(seq(min(chk.info$begin),max(chk.info$end)), rev(seq(min(chk.info$begin),max(chk.info$end)))), c(apply(exp(mu_loga.chk),2,quantile,0.05), rev(apply(exp(mu_loga.chk),2,quantile,0.95))), col = adjustcolor(sp_cols[1],alpha.f=0.1), border = NA)
lines(apply(exp(mu_loga.chm),2,median)~seq(min(chm.info$begin),max(chm.info$end)),lwd=3,col=adjustcolor(sp_cols[2],alpha.f=0.9))
polygon(c(seq(min(chm.info$begin),max(chm.info$end)), rev(seq(min(chm.info$begin),max(chm.info$end)))), c(apply(exp(mu_loga.chm),2,quantile,0.05), rev(apply(exp(mu_loga.chm),2,quantile,0.95))), col = adjustcolor(sp_cols[2],alpha.f=0.1), border = NA)
lines(apply(exp(mu_loga.cho),2,median)~seq(min(cho.info$begin),max(cho.info$end)),lwd=3,col=adjustcolor(sp_cols[3],alpha.f=0.9))
polygon(c(seq(min(cho.info$begin),max(cho.info$end)), rev(seq(min(cho.info$begin),max(cho.info$end)))), c(apply(exp(mu_loga.cho),2,quantile,0.05), rev(apply(exp(mu_loga.cho),2,quantile,0.95))), col = adjustcolor(sp_cols[3],alpha.f=0.1), border = NA)
lines(apply(exp(mu_loga.pke),2,median)~seq(min(pke.info$begin),max(pke.info$end)),lwd=3,col=adjustcolor(sp_cols[4],alpha.f=0.9))
polygon(c(seq(min(pke.info$begin),max(pke.info$end)), rev(seq(min(pke.info$begin),max(pke.info$end)))), c(apply(exp(mu_loga.pke),2,quantile,0.05), rev(apply(exp(mu_loga.pke),2,quantile,0.95))), col = adjustcolor(sp_cols[4],alpha.f=0.1), border = NA)
lines(apply(exp(mu_loga.pko),2,median)~seq(min(pko.info$begin),max(pko.info$end)),lwd=3,col=adjustcolor(sp_cols[4],alpha.f=0.9))
polygon(c(seq(min(pko.info$begin),max(pko.info$end)), rev(seq(min(pko.info$begin),max(pko.info$end)))), c(apply(exp(mu_loga.pko),2,quantile,0.05), rev(apply(exp(mu_loga.pko),2,quantile,0.95))), col = adjustcolor(sp_cols[4],alpha.f=0.1), border = NA)
lines(apply(exp(mu_loga.scy),2,median)~seq(min(scy.info$begin),max(scy.info$end)),lwd=3,col=adjustcolor(sp_cols[5],alpha.f=0.9))
polygon(c(seq(min(scy.info$begin),max(scy.info$end)), rev(seq(min(scy.info$begin),max(scy.info$end)))), c(apply(exp(mu_loga.scy),2,quantile,0.05), rev(apply(exp(mu_loga.scy),2,quantile,0.95))), col = adjustcolor(sp_cols[5],alpha.f=0.1), border = NA)


## Regional breakdown ####
### Chinook####
loga.chk=tv.chk$draws(variables='log_a_t',format='draws_matrix')
loga.chk.wc=loga.chk[,chk.wc]
loga.chk.seak=loga.chk[,chk.seak]
loga.chk.goa=loga.chk[,chk.goa]
loga.chk.bs=loga.chk[,chk.bs]

L=max(chk.info$end)-min(chk.info$begin)+1

#mean log(a) - WEst coast
mu_loga.chk.wc=matrix(nrow=nrow(loga.chk),ncol=L)
subset.rows=paste(',',chk.wc,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
wc.loga.chk=loga.chk[,grepl(subset.rows2,colnames(loga.chk))] #all log-a estimates for WC stocks
logat=wc.loga.chk
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.chk.wc[,t]=apply(logas,1,mean)
}    

#mean log(a) - SE Alaska/NBC
mu_loga.chk.seak=matrix(nrow=nrow(loga.chk),ncol=L)
subset.rows=paste(',',chk.seak,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
seak.loga.chk=loga.chk[,grepl(subset.rows2,colnames(loga.chk))] #all log-a estimates for WC stocks
logat=seak.loga.chk
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.chk.seak[,t]=apply(logas,1,mean)
}    

#mean log(a) - Gulf of Alaska
mu_loga.chk.goa=matrix(nrow=nrow(loga.chk),ncol=L)
subset.rows=paste(',',chk.goa,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
goa.loga.chk=loga.chk[,grepl(subset.rows2,colnames(loga.chk))] #all log-a estimates for WC stocks
logat=goa.loga.chk
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.chk.goa[,t]=apply(logas,1,mean)
}    

#mean log(a) - Bering sea
mu_loga.chk.bs=matrix(nrow=nrow(loga.chk),ncol=L)
subset.rows=paste(',',chk.bs,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
bs.loga.chk=loga.chk[,grepl(subset.rows2,colnames(loga.chk))] #all log-a estimates for WC stocks
logat=bs.loga.chk
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.chk.bs[,t]=apply(logas,1,mean)
}    

df=data.frame(by=seq(min(chk.info$begin),max(chk.info$end)),
              y.bs=apply(exp(mu_loga.chk.bs),2,median),y.bsl90=apply(exp(mu_loga.chk.bs),2,quantile,0.05),y.bsu90=apply(exp(mu_loga.chk.bs),2,quantile,0.95),
              y.goa=apply(exp(mu_loga.chk.goa),2,median),y.goal90=apply(exp(mu_loga.chk.goa),2,quantile,0.05),y.goau90=apply(exp(mu_loga.chk.goa),2,quantile,0.95),
              y.seak=apply(exp(mu_loga.chk.seak),2,median),y.seakl90=apply(exp(mu_loga.chk.seak),2,quantile,0.05),y.seaku90=apply(exp(mu_loga.chk.seak),2,quantile,0.95),
              y.wc=apply(exp(mu_loga.chk.wc),2,median),y.wcl90=apply(exp(mu_loga.chk.wc),2,quantile,0.05),y.wcu90=apply(exp(mu_loga.chk.wc),2,quantile,0.95)
              )
mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_blank(),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

prod=ggplot(df,aes(x=by,y=y.bs))+
  ylim(1,10)+
  geom_line(aes(y=y.bs),color=basin_cols[4],size=2)+
  geom_ribbon(aes(ymin = y.bsl90, ymax = y.bsu90), fill = basin_cols[4], alpha = 0.1) +
  geom_line(aes(y=y.goa),color=basin_cols[3],size=2)+
  geom_ribbon(aes(ymin = y.goal90, ymax = y.goau90), fill = basin_cols[3], alpha = 0.1) +
  geom_line(aes(y=y.seak),color=basin_cols[2],size=2)+
  geom_ribbon(aes(ymin = y.seakl90, ymax = y.seaku90), fill = basin_cols[2], alpha = 0.1) +
  geom_line(aes(y=y.wc),color=basin_cols[1],size=2)+
  geom_ribbon(aes(ymin = y.wcl90, ymax = y.wcu90), fill = basin_cols[1], alpha = 0.1) +
  xlab('Brood cohort year')+
  ylab(expression(paste('Mean productivity, max. R/S')))+
  ggtitle('Chinook')+
  mytheme

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

all_lat_lon.chk = chk.info %>% group_by(lat,lon) %>% summarize(n=n(),.groups = 'keep')
all_lat_lon.chk$lon=ifelse(all_lat_lon.chk$lon>0,-all_lat_lon.chk$lon,all_lat_lon.chk$lon) #need to fix the lat lons
all_lat_lon.chk$ocean.basin=chk.info$ocean.basin[match(all_lat_lon.chk$lat,chk.info$lat)]

ll_chk1=subset(all_lat_lon.chk,ocean.basin=='WC')
ll_chk2=subset(all_lat_lon.chk,ocean.basin=='SEAK')
ll_chk3=subset(all_lat_lon.chk,ocean.basin=='GOA')
ll_chk4=subset(all_lat_lon.chk,ocean.basin=='BS')


map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  ll_chk4, mapping = aes(x = lon, y = lat), color = basin_cols[4], size = 2.5*log(ll_chk4$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk3, mapping = aes(x = lon, y = lat), color =basin_cols[3], size = 2.5*log(ll_chk3$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk2, mapping = aes(x = lon, y = lat), color = basin_cols[2], size = 2.5*log(ll_chk2$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk1, mapping = aes(x = lon, y = lat), color =basin_cols[1], size = 2.5*log(ll_chk1$n+2), alpha = 0.7) +
  annotate(geom = 'point', x = -165, y = 50, size = 2.5*log(1+2)) +
  annotate(geom = 'point', x = -165, y = 49, size = 2.5*log(3+2)) + 
  annotate(geom = 'point', x = -165, y = 48, size = 2.5*log(10+2)) + 
  annotate(geom = 'point', x = -165, y = 47, size = 2.5*log(20+2)) + 
  annotate(geom = 'text', x = -165, y = 51,label='n', size = 6) +
  annotate(geom = 'text', x = -162, y = 50,label='1', size = 5) +
  annotate(geom = 'text', x = -162, y = 49,label='3', size = 5) + 
  annotate(geom = 'text', x = -162, y = 48,label='10', size = 5) + 
  annotate(geom = 'text', x = -162, y = 47,label='20', size = 5) + 
  annotate(geom = 'text', x = -134, y = 50.5, label = 'West Coast', color = basin_cols[1], size = 7) + 
  annotate(geom = 'text', x = -134, y = 49.5, label = paste('n =', length(chk.wc)), color = basin_cols[1], size = 5) + 
  annotate(geom = 'text', x = -138, y = 54.5, label = 'SE Alaska', color = basin_cols[2], size = 7) + 
  annotate(geom = 'text', x = -138, y = 53.5, label = paste('n =', length(chk.seak)), color = basin_cols[2], size = 5) + 
  annotate(geom = 'text', x = -145, y = 57.5, label = 'Gulf of Alaska', color = basin_cols[3], size = 7) + 
  annotate(geom = 'text', x = -145, y = 56.5, label = paste('n =', length(chk.goa)), color = basin_cols[3], size = 5) + 
  annotate(geom = 'text', x = -164, y = 58, label = 'Bering Sea', color = basin_cols[4], size = 7) + 
  annotate(geom = 'text', x = -164, y = 57, label = paste('n =', length(chk.bs)), color = basin_cols[4], size = 5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

cowplot::plot_grid(prod,map,labels=c('A.','B.'),rel_widths = c(0.5,0.5))

### Sockeye####
loga.scy=tv.scy$draws(variables='log_a_t',format='draws_matrix')
loga.scy.wc=loga.scy[,scy.wc]
loga.scy.seak=loga.scy[,scy.seak]
loga.scy.goa=loga.scy[,scy.goa]
loga.scy.bs=loga.scy[,scy.bs]

L=max(scy.info$end)-min(scy.info$begin)+1

#mean log(a) - WEst coast
mu_loga.scy.wc=matrix(nrow=nrow(loga.scy),ncol=L)
subset.rows=paste(',',scy.wc,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
wc.loga.scy=loga.scy[,grepl(subset.rows2,colnames(loga.scy))] #all log-a estimates for WC stocks
logat=wc.loga.scy
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.scy.wc[,t]=apply(logas,1,mean)
}    

#mean log(a) - SE Alaska/NBC
mu_loga.scy.seak=matrix(nrow=nrow(loga.scy),ncol=L)
subset.rows=paste(',',scy.seak,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
seak.loga.scy=loga.scy[,grepl(subset.rows2,colnames(loga.scy))] #all log-a estimates for WC stocks
logat=seak.loga.scy
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.scy.seak[,t]=apply(logas,1,mean)
}    

#mean log(a) - Gulf of Alaska
mu_loga.scy.goa=matrix(nrow=nrow(loga.scy),ncol=L)
subset.rows=paste(',',scy.goa,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
goa.loga.scy=loga.scy[,grepl(subset.rows2,colnames(loga.scy))] #all log-a estimates for WC stocks
logat=goa.loga.scy
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.scy.goa[,t]=apply(logas,1,mean)
}    

#mean log(a) - Bering sea
mu_loga.scy.bs=matrix(nrow=nrow(loga.scy),ncol=L)
subset.rows=paste(',',scy.bs,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
bs.loga.scy=loga.scy[,grepl(subset.rows2,colnames(loga.scy))] #all log-a estimates for WC stocks
logat=bs.loga.scy
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.scy.bs[,t]=apply(logas,1,mean)
}    

df=data.frame(by=seq(min(scy.info$begin),max(scy.info$end)),
              y.bs=apply(exp(mu_loga.scy.bs),2,median),y.bsl90=apply(exp(mu_loga.scy.bs),2,quantile,0.05),y.bsu90=apply(exp(mu_loga.scy.bs),2,quantile,0.95),
              y.goa=apply(exp(mu_loga.scy.goa),2,median),y.goal90=apply(exp(mu_loga.scy.goa),2,quantile,0.05),y.goau90=apply(exp(mu_loga.scy.goa),2,quantile,0.95),
              y.seak=apply(exp(mu_loga.scy.seak),2,median),y.seakl90=apply(exp(mu_loga.scy.seak),2,quantile,0.05),y.seaku90=apply(exp(mu_loga.scy.seak),2,quantile,0.95),
              y.wc=apply(exp(mu_loga.scy.wc),2,median),y.wcl90=apply(exp(mu_loga.scy.wc),2,quantile,0.05),y.wcu90=apply(exp(mu_loga.scy.wc),2,quantile,0.95)
)
mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_blank(),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

prod=ggplot(df,aes(x=by,y=y.bs))+
  ylim(1,10)+
  geom_line(aes(y=y.bs),color=basin_cols[4],size=2)+
  geom_ribbon(aes(ymin = y.bsl90, ymax = y.bsu90), fill = basin_cols[4], alpha = 0.1) +
  geom_line(aes(y=y.goa),color=basin_cols[3],size=2)+
  geom_ribbon(aes(ymin = y.goal90, ymax = y.goau90), fill = basin_cols[3], alpha = 0.1) +
  geom_line(aes(y=y.seak),color=basin_cols[2],size=2)+
  geom_ribbon(aes(ymin = y.seakl90, ymax = y.seaku90), fill = basin_cols[2], alpha = 0.1) +
  geom_line(aes(y=y.wc),color=basin_cols[1],size=2)+
  geom_ribbon(aes(ymin = y.wcl90, ymax = y.wcu90), fill = basin_cols[1], alpha = 0.1) +
  xlab('Brood cohort year')+
  ylab(expression(paste('Mean productivity, max. R/S')))+
  ggtitle('Sockeye')+
  mytheme

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

all_lat_lon.scy = scy.info %>% group_by(lat,lon) %>% summarize(n=n(),.groups = 'keep')
all_lat_lon.scy$lon=ifelse(all_lat_lon.scy$lon>0,-all_lat_lon.scy$lon,all_lat_lon.scy$lon) #need to fix the lat lons
all_lat_lon.scy$ocean.basin=scy.info$ocean.basin[match(all_lat_lon.scy$lat,scy.info$lat)]

ll_chk1=subset(all_lat_lon.scy,ocean.basin=='WC')
ll_chk2=subset(all_lat_lon.scy,ocean.basin=='SEAK')
ll_chk3=subset(all_lat_lon.scy,ocean.basin=='GOA')
ll_chk4=subset(all_lat_lon.scy,ocean.basin=='BS')


map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  ll_chk4, mapping = aes(x = lon, y = lat), color = basin_cols[4], size = 2.5*log(ll_chk4$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk3, mapping = aes(x = lon, y = lat), color =basin_cols[3], size = 2.5*log(ll_chk3$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk2, mapping = aes(x = lon, y = lat), color = basin_cols[2], size = 2.5*log(ll_chk2$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk1, mapping = aes(x = lon, y = lat), color =basin_cols[1], size = 2.5*log(ll_chk1$n+2), alpha = 0.7) +
  annotate(geom = 'point', x = -165, y = 50, size = 2.5*log(1+2)) +
  annotate(geom = 'point', x = -165, y = 49, size = 2.5*log(3+2)) + 
  annotate(geom = 'point', x = -165, y = 48, size = 2.5*log(10+2)) + 
  annotate(geom = 'point', x = -165, y = 47, size = 2.5*log(20+2)) + 
  annotate(geom = 'text', x = -165, y = 51,label='n', size = 6) +
  annotate(geom = 'text', x = -162, y = 50,label='1', size = 5) +
  annotate(geom = 'text', x = -162, y = 49,label='3', size = 5) + 
  annotate(geom = 'text', x = -162, y = 48,label='10', size = 5) + 
  annotate(geom = 'text', x = -162, y = 47,label='20', size = 5) + 
  annotate(geom = 'text', x = -134, y = 50.5, label = 'West Coast', color = basin_cols[1], size = 7) + 
  annotate(geom = 'text', x = -134, y = 49.5, label = paste('n =', length(scy.wc)), color = basin_cols[1], size = 5) + 
  annotate(geom = 'text', x = -138, y = 54.5, label = 'SE Alaska', color = basin_cols[2], size = 7) + 
  annotate(geom = 'text', x = -138, y = 53.5, label = paste('n =', length(scy.seak)), color = basin_cols[2], size = 5) + 
  annotate(geom = 'text', x = -145, y = 57.5, label = 'Gulf of Alaska', color = basin_cols[3], size = 7) + 
  annotate(geom = 'text', x = -145, y = 56.5, label = paste('n =', length(scy.goa)), color = basin_cols[3], size = 5) + 
  annotate(geom = 'text', x = -164, y = 58, label = 'Bering Sea', color = basin_cols[4], size = 7) + 
  annotate(geom = 'text', x = -164, y = 57, label = paste('n =', length(scy.bs)), color = basin_cols[4], size = 5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

cowplot::plot_grid(prod,map,labels=c('A.','B.'),rel_widths = c(0.5,0.5))

### Chum####
loga.chm=tv.chm$draws(variables='log_a_t',format='draws_matrix')
loga.chm.wc=loga.chm[,chm.wc]
loga.chm.seak=loga.chm[,chm.seak]
loga.chm.goa=loga.chm[,chm.goa]
loga.chm.bs=loga.chm[,chm.bs]

L=max(chm.info$end)-min(chm.info$begin)+1

#mean log(a) - WEst coast
mu_loga.chm.wc=matrix(nrow=nrow(loga.chm),ncol=L)
subset.rows=paste(',',chm.wc,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
wc.loga.chm=loga.chm[,grepl(subset.rows2,colnames(loga.chm))] #all log-a estimates for WC stocks
logat=wc.loga.chm
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.chm.wc[,t]=apply(logas,1,mean)
}    

#mean log(a) - SE Alaska/NBC
mu_loga.chm.seak=matrix(nrow=nrow(loga.chm),ncol=L)
subset.rows=paste(',',chm.seak,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
seak.loga.chm=loga.chm[,grepl(subset.rows2,colnames(loga.chm))] #all log-a estimates for WC stocks
logat=seak.loga.chm
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.chm.seak[,t]=apply(logas,1,mean)
}    

#mean log(a) - Gulf of Alaska
mu_loga.chm.goa=matrix(nrow=nrow(loga.chm),ncol=L)
subset.rows=paste(',',chm.goa,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
goa.loga.chm=loga.chm[,grepl(subset.rows2,colnames(loga.chm))] #all log-a estimates for WC stocks
logat=goa.loga.chm
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.chm.goa[,t]=apply(logas,1,mean)
}    

#mean log(a) - Bering sea
mu_loga.chm.bs=matrix(nrow=nrow(loga.chm),ncol=L)
subset.rows=paste(',',chm.bs,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
bs.loga.chm=loga.chm[,grepl(subset.rows2,colnames(loga.chm))] #all log-a estimates for WC stocks
logat=bs.loga.chm
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.chm.bs[,t]=apply(logas,1,mean)
}    

df=data.frame(by=seq(min(chm.info$begin),max(chm.info$end)),
              y.bs=apply(exp(mu_loga.chm.bs),2,median),y.bsl90=apply(exp(mu_loga.chm.bs),2,quantile,0.05),y.bsu90=apply(exp(mu_loga.chm.bs),2,quantile,0.95),
              y.goa=apply(exp(mu_loga.chm.goa),2,median),y.goal90=apply(exp(mu_loga.chm.goa),2,quantile,0.05),y.goau90=apply(exp(mu_loga.chm.goa),2,quantile,0.95),
              y.seak=apply(exp(mu_loga.chm.seak),2,median),y.seakl90=apply(exp(mu_loga.chm.seak),2,quantile,0.05),y.seaku90=apply(exp(mu_loga.chm.seak),2,quantile,0.95),
              y.wc=apply(exp(mu_loga.chm.wc),2,median),y.wcl90=apply(exp(mu_loga.chm.wc),2,quantile,0.05),y.wcu90=apply(exp(mu_loga.chm.wc),2,quantile,0.95)
)
mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_blank(),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

prod=ggplot(df,aes(x=by,y=y.bs))+
  ylim(1,10)+
  geom_line(aes(y=y.bs),color=basin_cols[4],size=2)+
  geom_ribbon(aes(ymin = y.bsl90, ymax = y.bsu90), fill = basin_cols[4], alpha = 0.1) +
  geom_line(aes(y=y.goa),color=basin_cols[3],size=2)+
  geom_ribbon(aes(ymin = y.goal90, ymax = y.goau90), fill = basin_cols[3], alpha = 0.1) +
  geom_line(aes(y=y.seak),color=basin_cols[2],size=2)+
  geom_ribbon(aes(ymin = y.seakl90, ymax = y.seaku90), fill = basin_cols[2], alpha = 0.1) +
  geom_line(aes(y=y.wc),color=basin_cols[1],size=2)+
  geom_ribbon(aes(ymin = y.wcl90, ymax = y.wcu90), fill = basin_cols[1], alpha = 0.1) +
  xlab('Brood cohort year')+
  ylab(expression(paste('Mean productivity, max. R/S')))+
  ggtitle('Chum')+
  mytheme

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

all_lat_lon.chm = chm.info %>% group_by(lat,lon) %>% summarize(n=n(),.groups = 'keep')
all_lat_lon.chm$lon=ifelse(all_lat_lon.chm$lon>0,-all_lat_lon.chm$lon,all_lat_lon.chm$lon) #need to fix the lat lons
all_lat_lon.chm$ocean.basin=chm.info$ocean.basin[match(all_lat_lon.chm$lat,chm.info$lat)]

ll_chk1=subset(all_lat_lon.chm,ocean.basin=='WC')
ll_chk2=subset(all_lat_lon.chm,ocean.basin=='SEAK')
ll_chk3=subset(all_lat_lon.chm,ocean.basin=='GOA')
ll_chk4=subset(all_lat_lon.chm,ocean.basin=='BS')


map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  ll_chk4, mapping = aes(x = lon, y = lat), color = basin_cols[4], size = 2.5*log(ll_chk4$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk3, mapping = aes(x = lon, y = lat), color =basin_cols[3], size = 2.5*log(ll_chk3$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk2, mapping = aes(x = lon, y = lat), color = basin_cols[2], size = 2.5*log(ll_chk2$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk1, mapping = aes(x = lon, y = lat), color =basin_cols[1], size = 2.5*log(ll_chk1$n+2), alpha = 0.7) +
  annotate(geom = 'point', x = -165, y = 50, size = 2.5*log(1+2)) +
  annotate(geom = 'point', x = -165, y = 49, size = 2.5*log(3+2)) + 
  annotate(geom = 'point', x = -165, y = 48, size = 2.5*log(10+2)) + 
  annotate(geom = 'point', x = -165, y = 47, size = 2.5*log(20+2)) + 
  annotate(geom = 'text', x = -165, y = 51,label='n', size = 6) +
  annotate(geom = 'text', x = -162, y = 50,label='1', size = 5) +
  annotate(geom = 'text', x = -162, y = 49,label='3', size = 5) + 
  annotate(geom = 'text', x = -162, y = 48,label='10', size = 5) + 
  annotate(geom = 'text', x = -162, y = 47,label='20', size = 5) + 
  annotate(geom = 'text', x = -134, y = 50.5, label = 'West Coast', color = basin_cols[1], size = 7) + 
  annotate(geom = 'text', x = -134, y = 49.5, label = paste('n =', length(chm.wc)), color = basin_cols[1], size = 5) + 
  annotate(geom = 'text', x = -138, y = 54.5, label = 'SE Alaska', color = basin_cols[2], size = 7) + 
  annotate(geom = 'text', x = -138, y = 53.5, label = paste('n =', length(chm.seak)), color = basin_cols[2], size = 5) + 
  annotate(geom = 'text', x = -145, y = 57.5, label = 'Gulf of Alaska', color = basin_cols[3], size = 7) + 
  annotate(geom = 'text', x = -145, y = 56.5, label = paste('n =', length(chm.goa)), color = basin_cols[3], size = 5) + 
  annotate(geom = 'text', x = -164, y = 58, label = 'Bering Sea', color = basin_cols[4], size = 7) + 
  annotate(geom = 'text', x = -164, y = 57, label = paste('n =', length(chm.bs)), color = basin_cols[4], size = 5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

cowplot::plot_grid(prod,map,labels=c('A.','B.'),rel_widths = c(0.5,0.5))

### Coho####
loga.cho=tv.cho$draws(variables='log_a_t',format='draws_matrix')
loga.cho.wc=loga.cho[,cho.wc]
loga.cho.seak=loga.cho[,cho.seak]
loga.cho.goa=loga.cho[,cho.goa]
loga.cho.bs=loga.cho[,cho.bs]

L=max(cho.info$end)-min(cho.info$begin)+1

#mean log(a) - WEst coast
mu_loga.cho.wc=matrix(nrow=nrow(loga.cho),ncol=L)
subset.rows=paste(',',cho.wc,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
wc.loga.cho=loga.cho[,grepl(subset.rows2,colnames(loga.cho))] #all log-a estimates for WC stocks
logat=wc.loga.cho
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.cho.wc[,t]=apply(logas,1,mean)
}    

#mean log(a) - SE Alaska/NBC
mu_loga.cho.seak=matrix(nrow=nrow(loga.cho),ncol=L)
subset.rows=paste(',',cho.seak,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
seak.loga.cho=loga.cho[,grepl(subset.rows2,colnames(loga.cho))] #all log-a estimates for WC stocks
logat=seak.loga.cho
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.cho.seak[,t]=apply(logas,1,mean)
}    

#mean log(a) - Gulf of Alaska
mu_loga.cho.goa=matrix(nrow=nrow(loga.cho),ncol=L)
subset.rows=paste(',',cho.goa,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
goa.loga.cho=loga.cho[,grepl(subset.rows2,colnames(loga.cho))] #all log-a estimates for WC stocks
logat=goa.loga.cho
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.cho.goa[,t]=apply(logas,1,mean)
}    

df=data.frame(by=seq(min(cho.info$begin),max(cho.info$end)),
              y.goa=apply(exp(mu_loga.cho.goa),2,median),y.goal90=apply(exp(mu_loga.cho.goa),2,quantile,0.05),y.goau90=apply(exp(mu_loga.cho.goa),2,quantile,0.95),
              y.seak=apply(exp(mu_loga.cho.seak),2,median),y.seakl90=apply(exp(mu_loga.cho.seak),2,quantile,0.05),y.seaku90=apply(exp(mu_loga.cho.seak),2,quantile,0.95),
              y.wc=apply(exp(mu_loga.cho.wc),2,median),y.wcl90=apply(exp(mu_loga.cho.wc),2,quantile,0.05),y.wcu90=apply(exp(mu_loga.cho.wc),2,quantile,0.95)
)
mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_blank(),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

prod=ggplot(df,aes(x=by,y=y.goa))+
  ylim(1,10)+
  geom_line(aes(y=y.goa),color=basin_cols[3],size=2)+
  geom_ribbon(aes(ymin = y.goal90, ymax = y.goau90), fill = basin_cols[3], alpha = 0.1) +
  geom_line(aes(y=y.seak),color=basin_cols[2],size=2)+
  geom_ribbon(aes(ymin = y.seakl90, ymax = y.seaku90), fill = basin_cols[2], alpha = 0.1) +
  geom_line(aes(y=y.wc),color=basin_cols[1],size=2)+
  geom_ribbon(aes(ymin = y.wcl90, ymax = y.wcu90), fill = basin_cols[1], alpha = 0.1) +
  xlab('Brood cohort year')+
  ylab(expression(paste('Mean productivity, max. R/S')))+
  ggtitle('Coho')+
  mytheme

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

all_lat_lon.cho = cho.info %>% group_by(lat,lon) %>% summarize(n=n(),.groups = 'keep')
all_lat_lon.cho$lon=ifelse(all_lat_lon.cho$lon>0,-all_lat_lon.cho$lon,all_lat_lon.cho$lon) #need to fix the lat lons
all_lat_lon.cho$ocean.basin=cho.info$ocean.basin[match(all_lat_lon.cho$lat,cho.info$lat)]

ll_chk1=subset(all_lat_lon.cho,ocean.basin=='WC')
ll_chk2=subset(all_lat_lon.cho,ocean.basin=='SEAK')
ll_chk3=subset(all_lat_lon.cho,ocean.basin=='GOA')


map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  ll_chk3, mapping = aes(x = lon, y = lat), color =basin_cols[3], size = 2.5*log(ll_chk3$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk2, mapping = aes(x = lon, y = lat), color = basin_cols[2], size = 2.5*log(ll_chk2$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk1, mapping = aes(x = lon, y = lat), color =basin_cols[1], size = 2.5*log(ll_chk1$n+2), alpha = 0.7) +
  annotate(geom = 'point', x = -165, y = 50, size = 2.5*log(1+2)) +
  annotate(geom = 'point', x = -165, y = 49, size = 2.5*log(3+2)) + 
  annotate(geom = 'point', x = -165, y = 48, size = 2.5*log(10+2)) + 
  annotate(geom = 'point', x = -165, y = 47, size = 2.5*log(20+2)) + 
  annotate(geom = 'text', x = -165, y = 51,label='n', size = 6) +
  annotate(geom = 'text', x = -162, y = 50,label='1', size = 5) +
  annotate(geom = 'text', x = -162, y = 49,label='3', size = 5) + 
  annotate(geom = 'text', x = -162, y = 48,label='10', size = 5) + 
  annotate(geom = 'text', x = -162, y = 47,label='20', size = 5) + 
  annotate(geom = 'text', x = -134, y = 50.5, label = 'West Coast', color = basin_cols[1], size = 7) + 
  annotate(geom = 'text', x = -134, y = 49.5, label = paste('n =', length(cho.wc)), color = basin_cols[1], size = 5) + 
  annotate(geom = 'text', x = -138, y = 54.5, label = 'SE Alaska', color = basin_cols[2], size = 7) + 
  annotate(geom = 'text', x = -138, y = 53.5, label = paste('n =', length(cho.seak)), color = basin_cols[2], size = 5) + 
  annotate(geom = 'text', x = -145, y = 57.5, label = 'Gulf of Alaska', color = basin_cols[3], size = 7) + 
  annotate(geom = 'text', x = -145, y = 56.5, label = paste('n =', length(cho.goa)), color = basin_cols[3], size = 5) + 
  annotate(geom = 'text', x = -164, y = 58, label = 'Bering Sea', color = basin_cols[4], size = 7) + 
  annotate(geom = 'text', x = -164, y = 57, label = paste('n =', length(cho.bs)), color = basin_cols[4], size = 5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

cowplot::plot_grid(prod,map,labels=c('A.','B.'),rel_widths = c(0.5,0.5))

### Pink-Even####
loga.pke=tv.pke$draws(variables='log_a_t',format='draws_matrix')
loga.pke.wc=loga.pke[,pke.wc]
loga.pke.seak=loga.pke[,pke.seak]
loga.pke.goa=loga.pke[,pke.goa]
loga.pke.bs=loga.pke[,pke.bs]

L=max(pke.info$end)-min(pke.info$begin)+1

#mean log(a) - WEst coast
mu_loga.pke.wc=matrix(nrow=nrow(loga.pke),ncol=L)
subset.rows=paste(',',pke.wc,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
wc.loga.pke=loga.pke[,grepl(subset.rows2,colnames(loga.pke))] #all log-a estimates for WC stocks
logat=wc.loga.pke
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.pke.wc[,t]=apply(logas,1,mean)
}    

#mean log(a) - SE Alaska/NBC
mu_loga.pke.seak=matrix(nrow=nrow(loga.pke),ncol=L)
subset.rows=paste(',',pke.seak,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
seak.loga.pke=loga.pke[,grepl(subset.rows2,colnames(loga.pke))] #all log-a estimates for WC stocks
logat=seak.loga.pke
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.pke.seak[,t]=apply(logas,1,mean)
}    

#mean log(a) - Gulf of Alaska
mu_loga.pke.goa=matrix(nrow=nrow(loga.pke),ncol=L)
subset.rows=paste(',',pke.goa,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
goa.loga.pke=loga.pke[,grepl(subset.rows2,colnames(loga.pke))] #all log-a estimates for WC stocks
logat=goa.loga.pke
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.pke.goa[,t]=apply(logas,1,mean)
}    

#mean log(a) - Bering sea
mu_loga.pke.bs=matrix(nrow=nrow(loga.pke),ncol=L)
subset.rows=paste(',',pke.bs,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
bs.loga.pke=loga.pke[,grepl(subset.rows2,colnames(loga.pke))] #all log-a estimates for WC stocks
logat=bs.loga.pke
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.pke.bs[,t]=apply(logas,1,mean)
}    

df=data.frame(by=seq(min(pke.info$begin),max(pke.info$end)),
              y.bs=apply(exp(mu_loga.pke.bs),2,median),y.bsl90=apply(exp(mu_loga.pke.bs),2,quantile,0.05),y.bsu90=apply(exp(mu_loga.pke.bs),2,quantile,0.95),
              y.goa=apply(exp(mu_loga.pke.goa),2,median),y.goal90=apply(exp(mu_loga.pke.goa),2,quantile,0.05),y.goau90=apply(exp(mu_loga.pke.goa),2,quantile,0.95),
              y.seak=apply(exp(mu_loga.pke.seak),2,median),y.seakl90=apply(exp(mu_loga.pke.seak),2,quantile,0.05),y.seaku90=apply(exp(mu_loga.pke.seak),2,quantile,0.95),
              y.wc=apply(exp(mu_loga.pke.wc),2,median),y.wcl90=apply(exp(mu_loga.pke.wc),2,quantile,0.05),y.wcu90=apply(exp(mu_loga.pke.wc),2,quantile,0.95)
)
mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_blank(),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

prod=ggplot(df,aes(x=by,y=y.bs))+
  ylim(1,10)+
  geom_line(aes(y=y.bs),color=basin_cols[4],size=2)+
  geom_ribbon(aes(ymin = y.bsl90, ymax = y.bsu90), fill = basin_cols[4], alpha = 0.1) +
  geom_line(aes(y=y.goa),color=basin_cols[3],size=2)+
  geom_ribbon(aes(ymin = y.goal90, ymax = y.goau90), fill = basin_cols[3], alpha = 0.1) +
  geom_line(aes(y=y.seak),color=basin_cols[2],size=2)+
  geom_ribbon(aes(ymin = y.seakl90, ymax = y.seaku90), fill = basin_cols[2], alpha = 0.1) +
  geom_line(aes(y=y.wc),color=basin_cols[1],size=2)+
  geom_ribbon(aes(ymin = y.wcl90, ymax = y.wcu90), fill = basin_cols[1], alpha = 0.1) +
  xlab('Brood cohort year')+
  ylab(expression(paste('Mean productivity, max. R/S')))+
  ggtitle('Pink - even year lines')+
  mytheme

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

all_lat_lon.pke = pke.info %>% group_by(lat,lon) %>% summarize(n=n(),.groups = 'keep')
all_lat_lon.pke$lon=ifelse(all_lat_lon.pke$lon>0,-all_lat_lon.pke$lon,all_lat_lon.pke$lon) #need to fix the lat lons
all_lat_lon.pke$ocean.basin=pke.info$ocean.basin[match(all_lat_lon.pke$lat,pke.info$lat)]

ll_chk1=subset(all_lat_lon.pke,ocean.basin=='WC')
ll_chk2=subset(all_lat_lon.pke,ocean.basin=='SEAK')
ll_chk3=subset(all_lat_lon.pke,ocean.basin=='GOA')
ll_chk4=subset(all_lat_lon.pke,ocean.basin=='BS')


map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  ll_chk4, mapping = aes(x = lon, y = lat), color = basin_cols[4], size = 2.5*log(ll_chk4$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk3, mapping = aes(x = lon, y = lat), color =basin_cols[3], size = 2.5*log(ll_chk3$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk2, mapping = aes(x = lon, y = lat), color = basin_cols[2], size = 2.5*log(ll_chk2$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk1, mapping = aes(x = lon, y = lat), color =basin_cols[1], size = 2.5*log(ll_chk1$n+2), alpha = 0.7) +
  annotate(geom = 'point', x = -165, y = 50, size = 2.5*log(1+2)) +
  annotate(geom = 'point', x = -165, y = 49, size = 2.5*log(3+2)) + 
  annotate(geom = 'point', x = -165, y = 48, size = 2.5*log(10+2)) + 
  annotate(geom = 'point', x = -165, y = 47, size = 2.5*log(20+2)) + 
  annotate(geom = 'text', x = -165, y = 51,label='n', size = 6) +
  annotate(geom = 'text', x = -162, y = 50,label='1', size = 5) +
  annotate(geom = 'text', x = -162, y = 49,label='3', size = 5) + 
  annotate(geom = 'text', x = -162, y = 48,label='10', size = 5) + 
  annotate(geom = 'text', x = -162, y = 47,label='20', size = 5) + 
  annotate(geom = 'text', x = -134, y = 50.5, label = 'West Coast', color = basin_cols[1], size = 7) + 
  annotate(geom = 'text', x = -134, y = 49.5, label = paste('n =', length(pke.wc)), color = basin_cols[1], size = 5) + 
  annotate(geom = 'text', x = -138, y = 54.5, label = 'SE Alaska', color = basin_cols[2], size = 7) + 
  annotate(geom = 'text', x = -138, y = 53.5, label = paste('n =', length(pke.seak)), color = basin_cols[2], size = 5) + 
  annotate(geom = 'text', x = -145, y = 57.5, label = 'Gulf of Alaska', color = basin_cols[3], size = 7) + 
  annotate(geom = 'text', x = -145, y = 56.5, label = paste('n =', length(pke.goa)), color = basin_cols[3], size = 5) + 
  annotate(geom = 'text', x = -164, y = 58, label = 'Bering Sea', color = basin_cols[4], size = 7) + 
  annotate(geom = 'text', x = -164, y = 57, label = paste('n =', length(pke.bs)), color = basin_cols[4], size = 5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

cowplot::plot_grid(prod,map,labels=c('A.','B.'),rel_widths = c(0.5,0.5))

### Pink-Odd####
loga.pko=tv.pko$draws(variables='log_a_t',format='draws_matrix')
loga.pko.wc=loga.pko[,pko.wc]
loga.pko.seak=loga.pko[,pko.seak]
loga.pko.goa=loga.pko[,pko.goa]
loga.pko.bs=loga.pko[,pko.bs]

L=max(pko.info$end)-min(pko.info$begin)+1

#mean log(a) - WEst coast
mu_loga.pko.wc=matrix(nrow=nrow(loga.pko),ncol=L)
subset.rows=paste(',',pko.wc,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
wc.loga.pko=loga.pko[,grepl(subset.rows2,colnames(loga.pko))] #all log-a estimates for WC stocks
logat=wc.loga.pko
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.pko.wc[,t]=apply(logas,1,mean)
}    

#mean log(a) - SE Alaska/NBC
mu_loga.pko.seak=matrix(nrow=nrow(loga.pko),ncol=L)
subset.rows=paste(',',pko.seak,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
seak.loga.pko=loga.pko[,grepl(subset.rows2,colnames(loga.pko))] #all log-a estimates for WC stocks
logat=seak.loga.pko
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.pko.seak[,t]=apply(logas,1,mean)
}    

#mean log(a) - Gulf of Alaska
mu_loga.pko.goa=matrix(nrow=nrow(loga.pko),ncol=L)
subset.rows=paste(',',pko.goa,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
goa.loga.pko=loga.pko[,grepl(subset.rows2,colnames(loga.pko))] #all log-a estimates for WC stocks
logat=goa.loga.pko
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.pko.goa[,t]=apply(logas,1,mean)
}    

#mean log(a) - Bering sea
mu_loga.pko.bs=matrix(nrow=nrow(loga.pko),ncol=L)
subset.rows=paste(',',pko.bs,']',sep='') #column name construction
subset.rows2=paste(subset.rows, collapse = "|") #conditional for grepl
bs.loga.pko=loga.pko[,grepl(subset.rows2,colnames(loga.pko))] #all log-a estimates for WC stocks
logat=bs.loga.pko
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.pko.bs[,t]=apply(logas,1,mean)
}    

df=data.frame(by=seq(min(pko.info$begin),max(pko.info$end)),
              y.bs=apply(exp(mu_loga.pko.bs),2,median),y.bsl90=apply(exp(mu_loga.pko.bs),2,quantile,0.05),y.bsu90=apply(exp(mu_loga.pko.bs),2,quantile,0.95),
              y.goa=apply(exp(mu_loga.pko.goa),2,median),y.goal90=apply(exp(mu_loga.pko.goa),2,quantile,0.05),y.goau90=apply(exp(mu_loga.pko.goa),2,quantile,0.95),
              y.seak=apply(exp(mu_loga.pko.seak),2,median),y.seakl90=apply(exp(mu_loga.pko.seak),2,quantile,0.05),y.seaku90=apply(exp(mu_loga.pko.seak),2,quantile,0.95),
              y.wc=apply(exp(mu_loga.pko.wc),2,median),y.wcl90=apply(exp(mu_loga.pko.wc),2,quantile,0.05),y.wcu90=apply(exp(mu_loga.pko.wc),2,quantile,0.95)
)
mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_blank(),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

prod=ggplot(df,aes(x=by,y=y.bs))+
  ylim(1,10)+
  geom_ribbon(aes(ymin = y.bsl90, ymax = y.bsu90), fill = basin_cols[4], alpha = 0.1) +
  geom_ribbon(aes(ymin = y.goal90, ymax = y.goau90), fill = basin_cols[3], alpha = 0.1) +
  geom_ribbon(aes(ymin = y.seakl90, ymax = y.seaku90), fill = basin_cols[2], alpha = 0.1) +
  geom_ribbon(aes(ymin = y.wcl90, ymax = y.wcu90), fill = basin_cols[1], alpha = 0.1) +
  geom_line(aes(y=y.bs),color=basin_cols[4],size=2)+
  geom_line(aes(y=y.goa),color=basin_cols[3],size=2)+
  geom_line(aes(y=y.seak),color=basin_cols[2],size=2)+
  geom_line(aes(y=y.wc),color=basin_cols[1],size=2)+
  xlab('Brood cohort year')+
  ylab(expression(paste('Mean productivity, max. R/S')))+
  ggtitle('Pink - Odd')+
  mytheme

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

all_lat_lon.pko = pko.info %>% group_by(lat,lon) %>% summarize(n=n(),.groups = 'keep')
all_lat_lon.pko$lon=ifelse(all_lat_lon.pko$lon>0,-all_lat_lon.pko$lon,all_lat_lon.pko$lon) #need to fix the lat lons
all_lat_lon.pko$ocean.basin=pko.info$ocean.basin[match(all_lat_lon.pko$lat,pko.info$lat)]

ll_chk1=subset(all_lat_lon.pko,ocean.basin=='WC')
ll_chk2=subset(all_lat_lon.pko,ocean.basin=='SEAK')
ll_chk3=subset(all_lat_lon.pko,ocean.basin=='GOA')
ll_chk4=subset(all_lat_lon.pko,ocean.basin=='BS')


map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  ll_chk4, mapping = aes(x = lon, y = lat), color = basin_cols[4], size = 2.5*log(ll_chk4$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk3, mapping = aes(x = lon, y = lat), color =basin_cols[3], size = 2.5*log(ll_chk3$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk2, mapping = aes(x = lon, y = lat), color = basin_cols[2], size = 2.5*log(ll_chk2$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk1, mapping = aes(x = lon, y = lat), color =basin_cols[1], size = 2.5*log(ll_chk1$n+2), alpha = 0.7) +
  annotate(geom = 'point', x = -165, y = 50, size = 2.5*log(1+2)) +
  annotate(geom = 'point', x = -165, y = 49, size = 2.5*log(3+2)) + 
  annotate(geom = 'point', x = -165, y = 48, size = 2.5*log(10+2)) + 
  annotate(geom = 'point', x = -165, y = 47, size = 2.5*log(20+2)) + 
  annotate(geom = 'text', x = -165, y = 51,label='n', size = 6) +
  annotate(geom = 'text', x = -162, y = 50,label='1', size = 5) +
  annotate(geom = 'text', x = -162, y = 49,label='3', size = 5) + 
  annotate(geom = 'text', x = -162, y = 48,label='10', size = 5) + 
  annotate(geom = 'text', x = -162, y = 47,label='20', size = 5) + 
  annotate(geom = 'text', x = -134, y = 50.5, label = 'West Coast', color = basin_cols[1], size = 7) + 
  annotate(geom = 'text', x = -134, y = 49.5, label = paste('n =', length(pko.wc)), color = basin_cols[1], size = 5) + 
  annotate(geom = 'text', x = -138, y = 54.5, label = 'SE Alaska', color = basin_cols[2], size = 7) + 
  annotate(geom = 'text', x = -138, y = 53.5, label = paste('n =', length(pko.seak)), color = basin_cols[2], size = 5) + 
  annotate(geom = 'text', x = -145, y = 57.5, label = 'Gulf of Alaska', color = basin_cols[3], size = 7) + 
  annotate(geom = 'text', x = -145, y = 56.5, label = paste('n =', length(pko.goa)), color = basin_cols[3], size = 5) + 
  annotate(geom = 'text', x = -164, y = 58, label = 'Bering Sea', color = basin_cols[4], size = 7) + 
  annotate(geom = 'text', x = -164, y = 57, label = paste('n =', length(pko.bs)), color = basin_cols[4], size = 5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

cowplot::plot_grid(prod,map,labels=c('A.','B.'),rel_widths = c(0.5,0.5))



#Productivity changes####

##Chinook####
loga.chk=tv.chk$draws(variables='log_a_t',format='draws_matrix')
st.loga.chk=stc.chk$draws(variables='log_a',format='draws_matrix')

alpha_comp=data.frame(stock=chk.info$stock.name,LME=chk.info$ocean.basin,begin=chk.info$begin,end=chk.info$end,r.alpha=NA,r.alpha.l90=NA,r.alpha.u90=NA,o.alpha=NA,o.alpha.l90=NA,o.alpha.u90=NA,diff.alpha=NA,diff.alpha.l90=NA,diff.alpha.u90=NA,pct.diff.alpha=NA,pct.diff.alpha.l90=NA,pct.diff.alpha.u90=NA,trend.pct.alpha=NA,trend.pct.alpha.l90=NA,trend.pct.alpha.u90=NA)
for(i in 1:nrow(chk.info)){
  logas=loga.chk[,grepl(paste(',',i,']',sep=''),colnames(loga.chk))]
  logas.st=exp(st.loga.chk[,i])
  
  o.logas=apply(exp(logas[,c(alpha_comp$begin[i]-min(alpha_comp$begin)+1):c(alpha_comp$begin[i]-min(alpha_comp$begin)+5)]),1,mean)
  r.logas=apply(exp(logas[,c(alpha_comp$end[i]-min(alpha_comp$begin)+1-5):c(alpha_comp$end[i]-min(alpha_comp$begin))]),1,mean)
  
  alpha_comp[i,5]=median(r.logas)
  alpha_comp[i,6]=quantile(r.logas,0.05)
  alpha_comp[i,7]=quantile(r.logas,0.95)
  alpha_comp[i,8]=median(logas.st)
  alpha_comp[i,9]=quantile(logas.st,0.05)
  alpha_comp[i,10]=quantile(logas.st,0.95)
  alpha_comp[i,11]=median(r.logas-logas.st)
  alpha_comp[i,12]=quantile(r.logas-logas.st,0.05)
  alpha_comp[i,13]=quantile(r.logas-logas.st,0.95)
  alpha_comp[i,14]=median((r.logas-logas.st)/logas.st)
  alpha_comp[i,15]=quantile((r.logas-logas.st)/logas.st,0.05)
  alpha_comp[i,16]=quantile((r.logas-logas.st)/logas.st,0.95)
  alpha_comp[i,17]=median((r.logas-logas.st)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas)
  alpha_comp[i,18]=quantile((r.logas-o.logas)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas,0.05)
  alpha_comp[i,19]=quantile((r.logas-o.logas)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas,0.95)
  
}
#transform % differences above a certain cut-off as a transformed bin category - for broken axis
t.alpha.diff=ifelse(alpha_comp$pct.diff.alpha*100>125,135,alpha_comp$pct.diff.alpha*100)

#colour palette specifications
p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts=cut(t.alpha.diff, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)
p_cols2=p_cols1[as.numeric(cuts)]
breaks=seq(-100,140,by=10)

par(mfrow=c(1,3),mar=c(5,5,3,1))
plot(r.alpha~o.alpha,data=alpha_comp,ylab=expression(paste('Recent productivity (last 5 cohorts mean  ',alpha['j,t'],')',sep=' ')),xlab=expression(paste('Long-term average productivity ',alpha[j],'',sep=' ')),bty='l',type='n',xlim=c(0,10),ylim=c(0,10))
abline(0,1)
for(i in 1:nrow(alpha_comp)){
  lines(rep(r.alpha[i],2)~c(o.alpha.l90[i],o.alpha.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=alpha_comp,lwd=0.8)
  lines(c(r.alpha.l90[i],r.alpha.u90[i])~rep(o.alpha[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=alpha_comp,lwd=0.8)
  lines(c(o.alpha[i],r.alpha[i])~rep(o.alpha[i],2),col=p_cols2[i],data=alpha_comp,lwd=2)
  
   }
points(r.alpha~o.alpha,data=alpha_comp,bg=adjustcolor(p_cols2,alpha.f=0.4),col=adjustcolor(p_cols2,alpha.f=1),pch=21,cex=1.6)

break_cols=p_cols1[cut(breaks,c(seq(-100,0,length.out=6),seq(20,100,length.out=5),141),right=F)]

hist(t.alpha.diff,breaks=breaks,xlab='% change from average',main='',border='white',col=break_cols,xaxt='n',xlim=c(-100,140))
axis(side=1,at=c(seq(-100,100,by=25),135),labels=c('-100','','-50','','0','','+50','','+100','>+125'))
plotrix::axis.break(1,breakpos=125,style='zigzag')
mtext(side=3,'Chinook',line=0.75,cex=1.25)

t.alpha.trend=ifelse(alpha_comp$trend.pct.alpha*100>5,6.5,alpha_comp$trend.pct.alpha*100)
breaks=seq(-5,7,by=0.5)
break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)]

hist(t.alpha.trend,breaks=breaks,xlab='Trend over series (annual % change)',main='',border='white',col=break_cols,xaxt='n',xlim=c(-5,6.5))
axis(side=1,at=c(seq(-5,5,by=2.5),6.5),labels=c('-5','-2.5','0','+2.5','+5','+>5'))
plotrix::axis.break(1,breakpos=6,style='zigzag')

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

alpha_comp$lat=chk.info$lat
alpha_comp$lon=chk.info$lon
c=cut(seq(-5,7),c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)
levels(c)[11]='5+'
leg=cbind(levels(c),p_cols1)

map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  alpha_comp, mapping = aes(x = lon, y = lat), color = p_cols2, size = 4, alpha = 0.7) +
  annotate(geom = 'text', x = -144, y = 57, label = 'Trend over series (%)', size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 46, size = 3.5, color = leg[1,2]) +
  annotate(geom = 'text', x = -143, y = 46, label = leg[1,1], color = leg[1,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 47, size = 3.5, color = leg[2,2]) +
  annotate(geom = 'text', x = -143, y = 47, label = leg[2,1], color = leg[2,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 48, size = 3.5, color = leg[3,2]) +
  annotate(geom = 'text', x = -143, y = 48, label = leg[3,1], color = leg[3,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 49, size = 3.5, color = leg[4,2]) +
  annotate(geom = 'text', x = -143, y = 49, label = leg[4,1], color = leg[4,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 50, size = 3.5, color = leg[5,2]) +
  annotate(geom = 'text', x = -143, y = 50, label = leg[5,1], color = leg[5,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 51, size = 3.5, color = leg[6,2]) +
  annotate(geom = 'text', x = -143, y = 51, label = leg[6,1], color = leg[6,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 52, size = 3.5, color = leg[7,2]) +
  annotate(geom = 'text', x = -143, y = 52, label = leg[7,1], color = leg[7,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 53, size = 3.5, color = leg[8,2]) +
  annotate(geom = 'text', x = -143, y = 53, label = leg[8,1], color = leg[8,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 54, size = 3.5, color = leg[9,2]) +
  annotate(geom = 'text', x = -143, y = 54, label = leg[9,1], color = leg[9,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 55, size = 3.5, color = leg[10,2]) +
  annotate(geom = 'text', x = -143, y = 55, label = leg[10,1], color = leg[10,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 56, size = 3.5, color = leg[11,2]) +
  annotate(geom = 'text', x = -143, y = 56, label = leg[11,1], color = leg[11,2], size = 3.5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

map

##Sockeye####
loga.scy=tv.scy$draws(variables='log_a_t',format='draws_matrix')
st.loga.scy=stc.scy$draws(variables='log_a',format='draws_matrix')

alpha_comp=data.frame(stock=scy.info$stock.name,LME=scy.info$ocean.basin,begin=scy.info$begin,end=scy.info$end,r.alpha=NA,r.alpha.l90=NA,r.alpha.u90=NA,o.alpha=NA,o.alpha.l90=NA,o.alpha.u90=NA,diff.alpha=NA,diff.alpha.l90=NA,diff.alpha.u90=NA,pct.diff.alpha=NA,pct.diff.alpha.l90=NA,pct.diff.alpha.u90=NA,trend.pct.alpha=NA,trend.pct.alpha.l90=NA,trend.pct.alpha.u90=NA)
for(i in 1:nrow(scy.info)){
  logas=loga.scy[,grepl(paste(',',i,']',sep=''),colnames(loga.scy))]
  logas.st=exp(st.loga.scy[,i])
  
  o.logas=apply(exp(logas[,c(alpha_comp$begin[i]-min(alpha_comp$begin)+1):c(alpha_comp$begin[i]-min(alpha_comp$begin)+5)]),1,mean)
  r.logas=apply(exp(logas[,c(alpha_comp$end[i]-min(alpha_comp$begin)+1-5):c(alpha_comp$end[i]-min(alpha_comp$begin))]),1,mean)
  
  alpha_comp[i,5]=median(r.logas)
  alpha_comp[i,6]=quantile(r.logas,0.05)
  alpha_comp[i,7]=quantile(r.logas,0.95)
  alpha_comp[i,8]=median(logas.st)
  alpha_comp[i,9]=quantile(logas.st,0.05)
  alpha_comp[i,10]=quantile(logas.st,0.95)
  alpha_comp[i,11]=median(r.logas-logas.st)
  alpha_comp[i,12]=quantile(r.logas-logas.st,0.05)
  alpha_comp[i,13]=quantile(r.logas-logas.st,0.95)
  alpha_comp[i,14]=median((r.logas-logas.st)/logas.st)
  alpha_comp[i,15]=quantile((r.logas-logas.st)/logas.st,0.05)
  alpha_comp[i,16]=quantile((r.logas-logas.st)/logas.st,0.95)
  alpha_comp[i,17]=median((r.logas-logas.st)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas)
  alpha_comp[i,18]=quantile((r.logas-o.logas)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas,0.05)
  alpha_comp[i,19]=quantile((r.logas-o.logas)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas,0.95)
  
}
t.alpha.diff=ifelse(alpha_comp$pct.diff.alpha*100>125,135,alpha_comp$pct.diff.alpha*100)

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts=cut(t.alpha.diff, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)
p_cols2=p_cols1[as.numeric(cuts)]
breaks=seq(-100,140,by=10)

par(mfrow=c(1,3),mar=c(5,5,3,1))
plot(r.alpha~o.alpha,data=alpha_comp,ylab=expression(paste('Recent productivity (last 5 cohorts mean  ',alpha['j,t'],')',sep=' ')),xlab=expression(paste('Long-term average productivity ',alpha[j],'',sep=' ')),bty='l',type='n',xlim=c(0,30),ylim=c(0,30))
abline(0,1)
for(i in 1:nrow(alpha_comp)){
  lines(rep(r.alpha[i],2)~c(o.alpha.l90[i],o.alpha.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=alpha_comp,lwd=0.8)
  lines(c(r.alpha.l90[i],r.alpha.u90[i])~rep(o.alpha[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=alpha_comp,lwd=0.8)
  lines(c(o.alpha[i],r.alpha[i])~rep(o.alpha[i],2),col=p_cols2[i],data=alpha_comp,lwd=2)
  
}
points(r.alpha~o.alpha,data=alpha_comp,bg=adjustcolor(p_cols2,alpha.f=0.4),col=adjustcolor(p_cols2,alpha.f=1),pch=21,cex=1.6)

break_cols=p_cols1[cut(breaks,c(seq(-100,0,length.out=6),seq(20,100,length.out=5),141),right=F)]

hist(t.alpha.diff,breaks=breaks,xlab='% change from average',main='',border='white',col=break_cols,xaxt='n',xlim=c(-100,140))
axis(side=1,at=c(seq(-100,100,by=25),135),labels=c('-100','','-50','','0','','+50','','+100','>+125'))
plotrix::axis.break(1,breakpos=125,style='zigzag')
mtext(side=3,'Sockeye',line=0.75,cex=1.25)


t.alpha.trend=ifelse(alpha_comp$trend.pct.alpha*100>5,6.5,alpha_comp$trend.pct.alpha*100)
breaks=seq(-5,7,by=0.5)
break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)]

hist(t.alpha.trend,breaks=breaks,xlab='Trend over series (annual % change)',main='',border='white',col=break_cols,xaxt='n',xlim=c(-5,6.5))
axis(side=1,at=c(seq(-5,5,by=2.5),6.5),labels=c('-5','-2.5','0','+2.5','+5','+>5'))
plotrix::axis.break(1,breakpos=6,style='zigzag')

alpha_comp$lat=scy.info$lat
alpha_comp$lon=scy.info$lon
c=cut(seq(-5,7),c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)
levels(c)[11]='5+'
leg=cbind(levels(c),p_cols1)

map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  alpha_comp, mapping = aes(x = lon, y = lat), color = p_cols2, size = 4, alpha = 0.7) +
  annotate(geom = 'text', x = -144, y = 57, label = 'Trend over series (%)', size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 46, size = 3.5, color = leg[1,2]) +
  annotate(geom = 'text', x = -143, y = 46, label = leg[1,1], color = leg[1,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 47, size = 3.5, color = leg[2,2]) +
  annotate(geom = 'text', x = -143, y = 47, label = leg[2,1], color = leg[2,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 48, size = 3.5, color = leg[3,2]) +
  annotate(geom = 'text', x = -143, y = 48, label = leg[3,1], color = leg[3,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 49, size = 3.5, color = leg[4,2]) +
  annotate(geom = 'text', x = -143, y = 49, label = leg[4,1], color = leg[4,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 50, size = 3.5, color = leg[5,2]) +
  annotate(geom = 'text', x = -143, y = 50, label = leg[5,1], color = leg[5,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 51, size = 3.5, color = leg[6,2]) +
  annotate(geom = 'text', x = -143, y = 51, label = leg[6,1], color = leg[6,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 52, size = 3.5, color = leg[7,2]) +
  annotate(geom = 'text', x = -143, y = 52, label = leg[7,1], color = leg[7,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 53, size = 3.5, color = leg[8,2]) +
  annotate(geom = 'text', x = -143, y = 53, label = leg[8,1], color = leg[8,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 54, size = 3.5, color = leg[9,2]) +
  annotate(geom = 'text', x = -143, y = 54, label = leg[9,1], color = leg[9,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 55, size = 3.5, color = leg[10,2]) +
  annotate(geom = 'text', x = -143, y = 55, label = leg[10,1], color = leg[10,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 56, size = 3.5, color = leg[11,2]) +
  annotate(geom = 'text', x = -143, y = 56, label = leg[11,1], color = leg[11,2], size = 3.5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

map

##Chum####
loga.chm=tv.chm$draws(variables='log_a_t',format='draws_matrix')
st.loga.chm=stc.chm$draws(variables='log_a',format='draws_matrix')

alpha_comp=data.frame(stock=chm.info$stock.name,LME=chm.info$ocean.basin,begin=chm.info$begin,end=chm.info$end,r.alpha=NA,r.alpha.l90=NA,r.alpha.u90=NA,o.alpha=NA,o.alpha.l90=NA,o.alpha.u90=NA,diff.alpha=NA,diff.alpha.l90=NA,diff.alpha.u90=NA,pct.diff.alpha=NA,pct.diff.alpha.l90=NA,pct.diff.alpha.u90=NA,trend.pct.alpha=NA,trend.pct.alpha.l90=NA,trend.pct.alpha.u90=NA)
for(i in 1:nrow(chm.info)){
  logas=loga.chm[,grepl(paste(',',i,']',sep=''),colnames(loga.chm))]
  logas.st=exp(st.loga.chm[,i])
  
  o.logas=apply(exp(logas[,c(alpha_comp$begin[i]-min(alpha_comp$begin)+1):c(alpha_comp$begin[i]-min(alpha_comp$begin)+5)]),1,mean)
  r.logas=apply(exp(logas[,c(alpha_comp$end[i]-min(alpha_comp$begin)+1-5):c(alpha_comp$end[i]-min(alpha_comp$begin))]),1,mean)
  
  alpha_comp[i,5]=median(r.logas)
  alpha_comp[i,6]=quantile(r.logas,0.05)
  alpha_comp[i,7]=quantile(r.logas,0.95)
  alpha_comp[i,8]=median(logas.st)
  alpha_comp[i,9]=quantile(logas.st,0.05)
  alpha_comp[i,10]=quantile(logas.st,0.95)
  alpha_comp[i,11]=median(r.logas-logas.st)
  alpha_comp[i,12]=quantile(r.logas-logas.st,0.05)
  alpha_comp[i,13]=quantile(r.logas-logas.st,0.95)
  alpha_comp[i,14]=median((r.logas-logas.st)/logas.st)
  alpha_comp[i,15]=quantile((r.logas-logas.st)/logas.st,0.05)
  alpha_comp[i,16]=quantile((r.logas-logas.st)/logas.st,0.95)
  alpha_comp[i,17]=median((r.logas-logas.st)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas)
  alpha_comp[i,18]=quantile((r.logas-o.logas)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas,0.05)
  alpha_comp[i,19]=quantile((r.logas-o.logas)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas,0.95)
  
}
t.alpha.diff=ifelse(alpha_comp$pct.diff.alpha*100>125,135,alpha_comp$pct.diff.alpha*100)

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts=cut(t.alpha.diff, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)
p_cols2=p_cols1[as.numeric(cuts)]
breaks=seq(-100,140,by=10)

par(mfrow=c(1,3),mar=c(5,5,3,1))
plot(r.alpha~o.alpha,data=alpha_comp,ylab=expression(paste('Recent productivity (last 5 cohorts mean  ',alpha['j,t'],')',sep=' ')),xlab=expression(paste('Long-term average productivity ',alpha[j],'',sep=' ')),bty='l',type='n',xlim=c(0,30),ylim=c(0,30))
abline(0,1)
for(i in 1:nrow(alpha_comp)){
  lines(rep(r.alpha[i],2)~c(o.alpha.l90[i],o.alpha.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=alpha_comp,lwd=0.8)
  lines(c(r.alpha.l90[i],r.alpha.u90[i])~rep(o.alpha[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=alpha_comp,lwd=0.8)
  lines(c(o.alpha[i],r.alpha[i])~rep(o.alpha[i],2),col=p_cols2[i],data=alpha_comp,lwd=2)
  
}
points(r.alpha~o.alpha,data=alpha_comp,bg=adjustcolor(p_cols2,alpha.f=0.4),col=adjustcolor(p_cols2,alpha.f=1),pch=21,cex=1.6)

break_cols=p_cols1[cut(breaks, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),141),right=F)]

hist(t.alpha.diff,breaks=breaks,xlab='% change from average',main='',border='white',col=break_cols,xaxt='n',xlim=c(-100,140),right = F)
axis(side=1,at=c(seq(-100,100,by=25),135),labels=c('-100','','-50','','0','','+50','','+100','>+125'))
plotrix::axis.break(1,breakpos=125,style='zigzag')
mtext(side=3,'Chum',line=0.75,cex=1.25)

t.alpha.trend=ifelse(alpha_comp$trend.pct.alpha*100>5,6.5,alpha_comp$trend.pct.alpha*100)
breaks=seq(-5,7,by=0.5)
break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)]

hist(t.alpha.trend,breaks=breaks,xlab='Trend over series (annual % change)',main='',border='white',col=break_cols,xaxt='n',xlim=c(-5,7))
axis(side=1,at=c(seq(-5,5,by=2.5),6.5),labels=c('-5','-2.5','0','+2.5','+5','+>5'))
plotrix::axis.break(1,breakpos=6,style='zigzag')

alpha_comp$lat=chm.info$lat
alpha_comp$lon=chm.info$lon
c=cut(seq(-5,7),c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)
levels(c)[11]='5+'
leg=cbind(levels(c),p_cols1)

map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  alpha_comp, mapping = aes(x = lon, y = lat), color = p_cols2, size = 4, alpha = 0.7) +
  annotate(geom = 'text', x = -144, y = 57, label = 'Trend over series (%)', size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 46, size = 3.5, color = leg[1,2]) +
  annotate(geom = 'text', x = -143, y = 46, label = leg[1,1], color = leg[1,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 47, size = 3.5, color = leg[2,2]) +
  annotate(geom = 'text', x = -143, y = 47, label = leg[2,1], color = leg[2,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 48, size = 3.5, color = leg[3,2]) +
  annotate(geom = 'text', x = -143, y = 48, label = leg[3,1], color = leg[3,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 49, size = 3.5, color = leg[4,2]) +
  annotate(geom = 'text', x = -143, y = 49, label = leg[4,1], color = leg[4,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 50, size = 3.5, color = leg[5,2]) +
  annotate(geom = 'text', x = -143, y = 50, label = leg[5,1], color = leg[5,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 51, size = 3.5, color = leg[6,2]) +
  annotate(geom = 'text', x = -143, y = 51, label = leg[6,1], color = leg[6,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 52, size = 3.5, color = leg[7,2]) +
  annotate(geom = 'text', x = -143, y = 52, label = leg[7,1], color = leg[7,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 53, size = 3.5, color = leg[8,2]) +
  annotate(geom = 'text', x = -143, y = 53, label = leg[8,1], color = leg[8,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 54, size = 3.5, color = leg[9,2]) +
  annotate(geom = 'text', x = -143, y = 54, label = leg[9,1], color = leg[9,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 55, size = 3.5, color = leg[10,2]) +
  annotate(geom = 'text', x = -143, y = 55, label = leg[10,1], color = leg[10,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 56, size = 3.5, color = leg[11,2]) +
  annotate(geom = 'text', x = -143, y = 56, label = leg[11,1], color = leg[11,2], size = 3.5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

map

##Coho####
loga.cho=tv.cho$draws(variables='log_a_t',format='draws_matrix')
st.loga.cho=stc.cho$draws(variables='log_a',format='draws_matrix')

alpha_comp=data.frame(stock=cho.info$stock.name,LME=cho.info$ocean.basin,begin=cho.info$begin,end=cho.info$end,r.alpha=NA,r.alpha.l90=NA,r.alpha.u90=NA,o.alpha=NA,o.alpha.l90=NA,o.alpha.u90=NA,diff.alpha=NA,diff.alpha.l90=NA,diff.alpha.u90=NA,pct.diff.alpha=NA,pct.diff.alpha.l90=NA,pct.diff.alpha.u90=NA,trend.pct.alpha=NA,trend.pct.alpha.l90=NA,trend.pct.alpha.u90=NA)
for(i in 1:nrow(cho.info)){
  logas=loga.cho[,grepl(paste(',',i,']',sep=''),colnames(loga.cho))]
  logas.st=exp(st.loga.cho[,i])
  
  o.logas=apply(exp(logas[,c(alpha_comp$begin[i]-min(alpha_comp$begin)+1):c(alpha_comp$begin[i]-min(alpha_comp$begin)+5)]),1,mean)
  r.logas=apply(exp(logas[,c(alpha_comp$end[i]-min(alpha_comp$begin)+1-5):c(alpha_comp$end[i]-min(alpha_comp$begin))]),1,mean)
  
  alpha_comp[i,5]=median(r.logas)
  alpha_comp[i,6]=quantile(r.logas,0.05)
  alpha_comp[i,7]=quantile(r.logas,0.95)
  alpha_comp[i,8]=median(logas.st)
  alpha_comp[i,9]=quantile(logas.st,0.05)
  alpha_comp[i,10]=quantile(logas.st,0.95)
  alpha_comp[i,11]=median(r.logas-logas.st)
  alpha_comp[i,12]=quantile(r.logas-logas.st,0.05)
  alpha_comp[i,13]=quantile(r.logas-logas.st,0.95)
  alpha_comp[i,14]=median((r.logas-logas.st)/logas.st)
  alpha_comp[i,15]=quantile((r.logas-logas.st)/logas.st,0.05)
  alpha_comp[i,16]=quantile((r.logas-logas.st)/logas.st,0.95)
  alpha_comp[i,17]=median((r.logas-logas.st)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas)
  alpha_comp[i,18]=quantile((r.logas-o.logas)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas,0.05)
  alpha_comp[i,19]=quantile((r.logas-o.logas)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas,0.95)
  
}
t.alpha.diff=ifelse(alpha_comp$pct.diff.alpha*100>125,135,alpha_comp$pct.diff.alpha*100)

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts=cut(t.alpha.diff, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)
p_cols2=p_cols1[as.numeric(cuts)]
breaks=seq(-100,140,by=10)

par(mfrow=c(1,3),mar=c(5,5,3,1))
plot(r.alpha~o.alpha,data=alpha_comp,ylab=expression(paste('Recent productivity (last 5 cohorts mean  ',alpha['j,t'],')',sep=' ')),xlab=expression(paste('Long-term average productivity ',alpha[j],'',sep=' ')),bty='l',type='n',xlim=c(0,45),ylim=c(0,45))
abline(0,1)
for(i in 1:nrow(alpha_comp)){
  lines(rep(r.alpha[i],2)~c(o.alpha.l90[i],o.alpha.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=alpha_comp,lwd=0.8)
  lines(c(r.alpha.l90[i],r.alpha.u90[i])~rep(o.alpha[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=alpha_comp,lwd=0.8)
  lines(c(o.alpha[i],r.alpha[i])~rep(o.alpha[i],2),col=p_cols2[i],data=alpha_comp,lwd=2)
  
}
points(r.alpha~o.alpha,data=alpha_comp,bg=adjustcolor(p_cols2,alpha.f=0.4),col=adjustcolor(p_cols2,alpha.f=1),pch=21,cex=1.6)

break_cols=p_cols1[cut(breaks,c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)]

hist(t.alpha.diff,breaks=breaks,xlab='% change from average',main='',border='white',col=break_cols,xaxt='n',xlim=c(-100,140))
axis(side=1,at=c(seq(-100,100,by=25),135),labels=c('-100','','-50','','0','','+50','','+100','>+125'))
plotrix::axis.break(1,breakpos=125,style='zigzag')
mtext(side=3,'Coho',line=0.75,cex=1.25)

t.alpha.trend=ifelse(alpha_comp$trend.pct.alpha*100>5,6.5,alpha_comp$trend.pct.alpha*100)
breaks=seq(-5,7,by=0.5)
break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)]

hist(t.alpha.trend,breaks=breaks,xlab='Trend over series (annual % change)',main='',border='white',col=break_cols,xaxt='n',xlim=c(-5,6.5))
axis(side=1,at=c(seq(-5,5,by=2.5),6.5),labels=c('-5','-2.5','0','+2.5','+5','+>5'))
plotrix::axis.break(1,breakpos=6,style='zigzag')

alpha_comp$lat=cho.info$lat
alpha_comp$lon=cho.info$lon
c=cut(seq(-5,7),c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)
levels(c)[11]='5+'
leg=cbind(levels(c),p_cols1)

map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  alpha_comp, mapping = aes(x = lon, y = lat), color = p_cols2, size = 4, alpha = 0.7) +
  annotate(geom = 'text', x = -144, y = 57, label = 'Trend over series (%)', size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 46, size = 3.5, color = leg[1,2]) +
  annotate(geom = 'text', x = -143, y = 46, label = leg[1,1], color = leg[1,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 47, size = 3.5, color = leg[2,2]) +
  annotate(geom = 'text', x = -143, y = 47, label = leg[2,1], color = leg[2,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 48, size = 3.5, color = leg[3,2]) +
  annotate(geom = 'text', x = -143, y = 48, label = leg[3,1], color = leg[3,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 49, size = 3.5, color = leg[4,2]) +
  annotate(geom = 'text', x = -143, y = 49, label = leg[4,1], color = leg[4,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 50, size = 3.5, color = leg[5,2]) +
  annotate(geom = 'text', x = -143, y = 50, label = leg[5,1], color = leg[5,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 51, size = 3.5, color = leg[6,2]) +
  annotate(geom = 'text', x = -143, y = 51, label = leg[6,1], color = leg[6,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 52, size = 3.5, color = leg[7,2]) +
  annotate(geom = 'text', x = -143, y = 52, label = leg[7,1], color = leg[7,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 53, size = 3.5, color = leg[8,2]) +
  annotate(geom = 'text', x = -143, y = 53, label = leg[8,1], color = leg[8,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 54, size = 3.5, color = leg[9,2]) +
  annotate(geom = 'text', x = -143, y = 54, label = leg[9,1], color = leg[9,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 55, size = 3.5, color = leg[10,2]) +
  annotate(geom = 'text', x = -143, y = 55, label = leg[10,1], color = leg[10,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 56, size = 3.5, color = leg[11,2]) +
  annotate(geom = 'text', x = -143, y = 56, label = leg[11,1], color = leg[11,2], size = 3.5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

map

##Pink - Even####
loga.pke=tv.pke$draws(variables='log_a_t',format='draws_matrix')
st.loga.pke=stc.pke$draws(variables='log_a',format='draws_matrix')

alpha_comp=data.frame(stock=pke.info$stock.name,LME=pke.info$ocean.basin,begin=pke.info$begin,end=pke.info$end,r.alpha=NA,r.alpha.l90=NA,r.alpha.u90=NA,o.alpha=NA,o.alpha.l90=NA,o.alpha.u90=NA,diff.alpha=NA,diff.alpha.l90=NA,diff.alpha.u90=NA,pct.diff.alpha=NA,pct.diff.alpha.l90=NA,pct.diff.alpha.u90=NA,trend.pct.alpha=NA,trend.pct.alpha.l90=NA,trend.pct.alpha.u90=NA)
for(i in 1:nrow(pke.info)){
  logas=loga.pke[,grepl(paste(',',i,']',sep=''),colnames(loga.pke))]
  logas.st=exp(st.loga.pke[,i])
  
  o.logas=apply(exp(logas[,c(alpha_comp$begin[i]-min(alpha_comp$begin)+1):c(alpha_comp$begin[i]-min(alpha_comp$begin)+5)]),1,mean)
  r.logas=apply(exp(logas[,c(alpha_comp$end[i]-min(alpha_comp$begin)+1-5):c(alpha_comp$end[i]-min(alpha_comp$begin))]),1,mean)
  
  alpha_comp[i,5]=median(r.logas)
  alpha_comp[i,6]=quantile(r.logas,0.05)
  alpha_comp[i,7]=quantile(r.logas,0.95)
  alpha_comp[i,8]=median(logas.st)
  alpha_comp[i,9]=quantile(logas.st,0.05)
  alpha_comp[i,10]=quantile(logas.st,0.95)
  alpha_comp[i,11]=median(r.logas-logas.st)
  alpha_comp[i,12]=quantile(r.logas-logas.st,0.05)
  alpha_comp[i,13]=quantile(r.logas-logas.st,0.95)
  alpha_comp[i,14]=median((r.logas-logas.st)/logas.st)
  alpha_comp[i,15]=quantile((r.logas-logas.st)/logas.st,0.05)
  alpha_comp[i,16]=quantile((r.logas-logas.st)/logas.st,0.95)
  alpha_comp[i,17]=median((r.logas-logas.st)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas)
  alpha_comp[i,18]=quantile((r.logas-o.logas)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas,0.05)
  alpha_comp[i,19]=quantile((r.logas-o.logas)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas,0.95)
  
}
t.alpha.diff=ifelse(alpha_comp$pct.diff.alpha*100>125,135,alpha_comp$pct.diff.alpha*100)

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts=cut(t.alpha.diff, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)
p_cols2=p_cols1[as.numeric(cuts)]
breaks=seq(-100,140,by=10)

par(mfrow=c(1,3),mar=c(5,5,3,1))
plot(r.alpha~o.alpha,data=alpha_comp,ylab=expression(paste('Recent productivity (last 5 cohorts mean  ',alpha['j,t'],')',sep=' ')),xlab=expression(paste('Long-term average productivity ',alpha[j],'',sep=' ')),bty='l',type='n',xlim=c(0,20),ylim=c(0,20))
abline(0,1)
for(i in 1:nrow(alpha_comp)){
  lines(rep(r.alpha[i],2)~c(o.alpha.l90[i],o.alpha.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=alpha_comp,lwd=0.8)
  lines(c(r.alpha.l90[i],r.alpha.u90[i])~rep(o.alpha[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=alpha_comp,lwd=0.8)
  lines(c(o.alpha[i],r.alpha[i])~rep(o.alpha[i],2),col=p_cols2[i],data=alpha_comp,lwd=2)
  
}
points(r.alpha~o.alpha,data=alpha_comp,bg=adjustcolor(p_cols2,alpha.f=0.4),col=adjustcolor(p_cols2,alpha.f=1),pch=21,cex=1.6)

break_cols=p_cols1[cut(breaks,c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)]

hist(t.alpha.diff,breaks=breaks,xlab='% change from average',main='',border='white',col=break_cols,xaxt='n',xlim=c(-100,140))
axis(side=1,at=c(seq(-100,100,by=25),135),labels=c('-100','','-50','','0','','+50','','+100','>+125'))
plotrix::axis.break(1,breakpos=125,style='zigzag')
mtext(side=3,'Pink - even year lines',line=0.75,cex=1.25)

t.alpha.trend=ifelse(alpha_comp$trend.pct.alpha*100>5,6.5,alpha_comp$trend.pct.alpha*100)
breaks=seq(-5,7,by=0.5)
break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)]

hist(t.alpha.trend,breaks=breaks,xlab='Trend over series (annual % change)',main='',border='white',col=break_cols,xaxt='n',xlim=c(-5,6.5))
axis(side=1,at=c(seq(-5,5,by=2.5),6.5),labels=c('-5','-2.5','0','+2.5','+5','+>5'))
plotrix::axis.break(1,breakpos=6,style='zigzag')

alpha_comp$lat=pke.info$lat
alpha_comp$lon=pke.info$lon
c=cut(seq(-5,7),c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)
levels(c)[11]='5+'
leg=cbind(levels(c),p_cols1)

map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  alpha_comp, mapping = aes(x = lon, y = lat), color = p_cols2, size = 4, alpha = 0.7) +
  annotate(geom = 'text', x = -144, y = 57, label = 'Trend over series (%)', size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 46, size = 3.5, color = leg[1,2]) +
  annotate(geom = 'text', x = -143, y = 46, label = leg[1,1], color = leg[1,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 47, size = 3.5, color = leg[2,2]) +
  annotate(geom = 'text', x = -143, y = 47, label = leg[2,1], color = leg[2,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 48, size = 3.5, color = leg[3,2]) +
  annotate(geom = 'text', x = -143, y = 48, label = leg[3,1], color = leg[3,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 49, size = 3.5, color = leg[4,2]) +
  annotate(geom = 'text', x = -143, y = 49, label = leg[4,1], color = leg[4,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 50, size = 3.5, color = leg[5,2]) +
  annotate(geom = 'text', x = -143, y = 50, label = leg[5,1], color = leg[5,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 51, size = 3.5, color = leg[6,2]) +
  annotate(geom = 'text', x = -143, y = 51, label = leg[6,1], color = leg[6,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 52, size = 3.5, color = leg[7,2]) +
  annotate(geom = 'text', x = -143, y = 52, label = leg[7,1], color = leg[7,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 53, size = 3.5, color = leg[8,2]) +
  annotate(geom = 'text', x = -143, y = 53, label = leg[8,1], color = leg[8,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 54, size = 3.5, color = leg[9,2]) +
  annotate(geom = 'text', x = -143, y = 54, label = leg[9,1], color = leg[9,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 55, size = 3.5, color = leg[10,2]) +
  annotate(geom = 'text', x = -143, y = 55, label = leg[10,1], color = leg[10,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 56, size = 3.5, color = leg[11,2]) +
  annotate(geom = 'text', x = -143, y = 56, label = leg[11,1], color = leg[11,2], size = 3.5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

map


##Pink - Odd####
loga.pko=tv.pko$draws(variables='log_a_t',format='draws_matrix')
st.loga.pko=stc.pko$draws(variables='log_a',format='draws_matrix')

alpha_comp=data.frame(stock=pko.info$stock.name,LME=pko.info$ocean.basin,begin=pko.info$begin,end=pko.info$end,r.alpha=NA,r.alpha.l90=NA,r.alpha.u90=NA,o.alpha=NA,o.alpha.l90=NA,o.alpha.u90=NA,diff.alpha=NA,diff.alpha.l90=NA,diff.alpha.u90=NA,pct.diff.alpha=NA,pct.diff.alpha.l90=NA,pct.diff.alpha.u90=NA,trend.pct.alpha=NA,trend.pct.alpha.l90=NA,trend.pct.alpha.u90=NA)
for(i in 1:nrow(pko.info)){
  logas=loga.pko[,grepl(paste(',',i,']',sep=''),colnames(loga.pko))]
  logas.st=exp(st.loga.pko[,i])
  
  o.logas=apply(exp(logas[,c(alpha_comp$begin[i]-min(alpha_comp$begin)+1):c(alpha_comp$begin[i]-min(alpha_comp$begin)+5)]),1,mean)
  r.logas=apply(exp(logas[,c(alpha_comp$end[i]-min(alpha_comp$begin)+1-5):c(alpha_comp$end[i]-min(alpha_comp$begin))]),1,mean)
  
  alpha_comp[i,5]=median(r.logas)
  alpha_comp[i,6]=quantile(r.logas,0.05)
  alpha_comp[i,7]=quantile(r.logas,0.95)
  alpha_comp[i,8]=median(logas.st)
  alpha_comp[i,9]=quantile(logas.st,0.05)
  alpha_comp[i,10]=quantile(logas.st,0.95)
  alpha_comp[i,11]=median(r.logas-logas.st)
  alpha_comp[i,12]=quantile(r.logas-logas.st,0.05)
  alpha_comp[i,13]=quantile(r.logas-logas.st,0.95)
  alpha_comp[i,14]=median((r.logas-logas.st)/logas.st)
  alpha_comp[i,15]=quantile((r.logas-logas.st)/logas.st,0.05)
  alpha_comp[i,16]=quantile((r.logas-logas.st)/logas.st,0.95)
  alpha_comp[i,17]=median((r.logas-logas.st)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas)
  alpha_comp[i,18]=quantile((r.logas-o.logas)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas,0.05)
  alpha_comp[i,19]=quantile((r.logas-o.logas)/(alpha_comp$end[i]-alpha_comp$begin[i]+1)/o.logas,0.95)
  
}
t.alpha.diff=ifelse(alpha_comp$pct.diff.alpha*100>125,135,alpha_comp$pct.diff.alpha*100)

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts=cut(t.alpha.diff, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)
p_cols2=p_cols1[as.numeric(cuts)]
breaks=seq(-100,140,by=10)

par(mfrow=c(1,3),mar=c(5,5,3,1))
plot(r.alpha~o.alpha,data=alpha_comp,ylab=expression(paste('Recent productivity (last 5 cohorts mean  ',alpha['j,t'],')',sep=' ')),xlab=expression(paste('Long-term average productivity ',alpha[j],'',sep=' ')),bty='l',type='n',xlim=c(0,20),ylim=c(0,20))
abline(0,1)
for(i in 1:nrow(alpha_comp)){
  lines(rep(r.alpha[i],2)~c(o.alpha.l90[i],o.alpha.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=alpha_comp,lwd=0.8)
  lines(c(r.alpha.l90[i],r.alpha.u90[i])~rep(o.alpha[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=alpha_comp,lwd=0.8)
  lines(c(o.alpha[i],r.alpha[i])~rep(o.alpha[i],2),col=p_cols2[i],data=alpha_comp,lwd=2)
  
}
points(r.alpha~o.alpha,data=alpha_comp,bg=adjustcolor(p_cols2,alpha.f=0.4),col=adjustcolor(p_cols2,alpha.f=1),pch=21,cex=1.6)

break_cols=p_cols1[cut(breaks,c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)]

hist(t.alpha.diff,breaks=breaks,xlab='% change from average',main='',border='white',col=break_cols,xaxt='n',xlim=c(-100,140))
axis(side=1,at=c(seq(-100,100,by=25),135),labels=c('-100','','-50','','0','','+50','','+100','>+125'))
plotrix::axis.break(1,breakpos=125,style='zigzag')
mtext(side=3,'Pink - odd year lines',line=0.75,cex=1.25)

t.alpha.trend=ifelse(alpha_comp$trend.pct.alpha*100>5,6.5,alpha_comp$trend.pct.alpha*100)
breaks=seq(-5,7,by=0.5)
break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)]

hist(t.alpha.trend,breaks=breaks,xlab='Trend over series (annual % change)',main='',border='white',col=break_cols,xaxt='n',xlim=c(-5,6.5))
axis(side=1,at=c(seq(-5,5,by=2.5),6.5),labels=c('-5','-2.5','0','+2.5','+5','+>5'))
plotrix::axis.break(1,breakpos=6,style='zigzag')

alpha_comp$lat=pko.info$lat
alpha_comp$lon=pko.info$lon
c=cut(seq(-5,7),c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)
levels(c)[11]='5+'
leg=cbind(levels(c),p_cols1)

map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  alpha_comp, mapping = aes(x = lon, y = lat), color = p_cols2, size = 4, alpha = 0.7) +
  annotate(geom = 'text', x = -144, y = 57, label = 'Trend over series (%)', size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 46, size = 3.5, color = leg[1,2]) +
  annotate(geom = 'text', x = -143, y = 46, label = leg[1,1], color = leg[1,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 47, size = 3.5, color = leg[2,2]) +
  annotate(geom = 'text', x = -143, y = 47, label = leg[2,1], color = leg[2,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 48, size = 3.5, color = leg[3,2]) +
  annotate(geom = 'text', x = -143, y = 48, label = leg[3,1], color = leg[3,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 49, size = 3.5, color = leg[4,2]) +
  annotate(geom = 'text', x = -143, y = 49, label = leg[4,1], color = leg[4,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 50, size = 3.5, color = leg[5,2]) +
  annotate(geom = 'text', x = -143, y = 50, label = leg[5,1], color = leg[5,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 51, size = 3.5, color = leg[6,2]) +
  annotate(geom = 'text', x = -143, y = 51, label = leg[6,1], color = leg[6,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 52, size = 3.5, color = leg[7,2]) +
  annotate(geom = 'text', x = -143, y = 52, label = leg[7,1], color = leg[7,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 53, size = 3.5, color = leg[8,2]) +
  annotate(geom = 'text', x = -143, y = 53, label = leg[8,1], color = leg[8,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 54, size = 3.5, color = leg[9,2]) +
  annotate(geom = 'text', x = -143, y = 54, label = leg[9,1], color = leg[9,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 55, size = 3.5, color = leg[10,2]) +
  annotate(geom = 'text', x = -143, y = 55, label = leg[10,1], color = leg[10,2], size = 3.5) + 
  annotate(geom = 'point', x = -145, y = 56, size = 3.5, color = leg[11,2]) +
  annotate(geom = 'text', x = -143, y = 56, label = leg[11,1], color = leg[11,2], size = 3.5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.background = element_rect(fill = 'white'),legend.title = element_blank())

map
#Umsy changes####

##Chinook####
umsy.chk=tv.chk$draws(variables='Umsy',format='draws_matrix')
st.umsy.chk=stc.chk$draws(variables='Umsy',format='draws_matrix')

#set a lower limit just above 0 to prevent negative Umsy values
umsy.chk=apply(umsy.chk,1:2,function(x) {ifelse(x<0,0.005,x)})
st.umsy.chk=apply(st.umsy.chk,1:2,function(x) {ifelse(x<0,0.005,x)})

umsy_comp=data.frame(stock=chk.info$stock.name,LME=chk.info$ocean.basin,begin=chk.info$begin,end=chk.info$end,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,o.umsy=NA,o.umsy.l90=NA,o.umsy.u90=NA,diff.umsy=NA,diff.umsy.l90=NA,diff.umsy.u90=NA,pct.diff.umsy=NA,pct.diff.umsy.l90=NA,pct.diff.umsy.u90=NA,trend.pct.umsy=NA,trend.pct.umsy.l90=NA,trend.pct.umsy.u90=NA)
for(i in 1:nrow(chk.info)){
  umsys=umsy.chk[,grepl(paste(',',i,']',sep=''),colnames(umsy.chk))]
  umsy.st=st.umsy.chk[,i]
  
  o.umsys=apply(umsys[,c(umsy_comp$begin[i]-min(umsy_comp$begin)+1):c(umsy_comp$begin[i]-min(umsy_comp$begin)+5)],1,mean)
  r.umsys=apply(umsys[,c(umsy_comp$end[i]-min(umsy_comp$begin)+1-5):c(umsy_comp$end[i]-min(umsy_comp$begin))],1,mean)
 # o.umsys=sapply(o.umsys,function(x) {ifelse(x<0,0.001,x)})#set lower limit of 0
#  r.umsys=sapply(r.umsys,function(x) {ifelse(x<0,0.001,x)})#set lower limit of 0
  
  umsy_comp[i,5]=median(r.umsys)
  umsy_comp[i,6]=quantile(r.umsys,0.05)
  umsy_comp[i,7]=quantile(r.umsys,0.95)
  umsy_comp[i,8]=median(umsy.st)
  umsy_comp[i,9]=quantile(umsy.st,0.05)
  umsy_comp[i,10]=quantile(umsy.st,0.95)
  umsy_comp[i,11]=median(r.umsys-umsy.st)
  umsy_comp[i,12]=quantile(r.umsys-umsy.st,0.05)
  umsy_comp[i,13]=quantile(r.umsys-umsy.st,0.95)
  umsy_comp[i,14]=median((r.umsys-umsy.st)/umsy.st)
  umsy_comp[i,15]=quantile((r.umsys-umsy.st)/umsy.st,0.05)
  umsy_comp[i,16]=quantile((r.umsys-umsy.st)/umsy.st,0.95)
  umsy_comp[i,17]=median((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys)
  umsy_comp[i,18]=quantile((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys,0.05)
  umsy_comp[i,19]=quantile((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys,0.95)
}
t.umsy.diff=ifelse(umsy_comp$pct.diff.umsy*100>125,135,umsy_comp$pct.diff.umsy*100)

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts=cut(t.umsy.diff, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)
p_cols2=p_cols1[as.numeric(cuts)]
breaks=seq(-100,140,by=10)

par(mfrow=c(1,3),mar=c(5,5,3,1))
plot(r.umsy~o.umsy,data=umsy_comp,ylab=expression(paste('Recent productivity, last 5 cohorts mean ',U['MSY,j,t'],sep='')),xlab=expression(paste('Long-term average ',U['MSY,j'],sep='')),bty='l',type='n',xlim=c(0,1),ylim=c(0,1))
abline(0,1)
for(i in 1:nrow(umsy_comp)){
  lines(rep(r.umsy[i],2)~c(o.umsy.l90[i],o.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp,lwd=0.8)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(o.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp,lwd=0.8)
  lines(c(o.umsy[i],r.umsy[i])~rep(o.umsy[i],2),col=p_cols2[i],data=umsy_comp,lwd=2)
  
}
points(r.umsy~o.umsy,data=umsy_comp,bg=adjustcolor(p_cols2,alpha.f=0.4),col=adjustcolor(p_cols2,alpha.f=1),pch=21,cex=1.6)

break_cols=p_cols1[cut(breaks,c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)]

hist(t.umsy.diff,breaks=breaks,xlab='% change from average',main='',border='white',col=break_cols,xaxt='n',xlim=c(-100,140))
axis(side=1,at=c(seq(-100,100,by=25),135),labels=c('-100','','-50','','0','','+50','','+100','>+125'))
plotrix::axis.break(1,breakpos=125,style='zigzag')
mtext(side=3,'Chinook',line=0.75,cex=1.25)

t.umsy.trend=ifelse(umsy_comp$trend.pct.umsy*100>5,6.5,umsy_comp$trend.pct.umsy*100)
breaks=seq(-5,7,by=0.5)
break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)]

hist(t.umsy.trend,breaks=breaks,xlab='Trend over series (annual % change)',main='',border='white',col=break_cols,xaxt='n',xlim=c(-5,6.5))
axis(side=1,at=c(seq(-5,5,by=2.5),6.5),labels=c('-5','-2.5','0','+2.5','+5','+>5'))
plotrix::axis.break(1,breakpos=6,style='zigzag')

##Sockeye####
umsy.scy=tv.scy$draws(variables='Umsy',format='draws_matrix')
st.umsy.scy=stc.scy$draws(variables='Umsy',format='draws_matrix')

#set a lower limit just above 0 to prevent negative Umsy values
umsy.scy=apply(umsy.scy,1:2,function(x) {ifelse(x<0,0.005,x)})
st.umsy.scy=apply(st.umsy.scy,1:2,function(x) {ifelse(x<0,0.005,x)})

umsy_comp=data.frame(stock=scy.info$stock.name,LME=scy.info$ocean.basin,begin=scy.info$begin,end=scy.info$end,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,o.umsy=NA,o.umsy.l90=NA,o.umsy.u90=NA,diff.umsy=NA,diff.umsy.l90=NA,diff.umsy.u90=NA,pct.diff.umsy=NA,pct.diff.umsy.l90=NA,pct.diff.umsy.u90=NA,trend.pct.umsy=NA,trend.pct.umsy.l90=NA,trend.pct.umsy.u90=NA)
for(i in 1:nrow(scy.info)){
  umsys=umsy.scy[,grepl(paste(',',i,']',sep=''),colnames(umsy.scy))]
  umsy.st=st.umsy.scy[,i]
  
  o.umsys=apply(umsys[,c(umsy_comp$begin[i]-min(umsy_comp$begin)+1):c(umsy_comp$begin[i]-min(umsy_comp$begin)+5)],1,mean)
  r.umsys=apply(umsys[,c(umsy_comp$end[i]-min(umsy_comp$begin)+1-5):c(umsy_comp$end[i]-min(umsy_comp$begin))],1,mean)
  # o.umsys=sapply(o.umsys,function(x) {ifelse(x<0,0.001,x)})#set lower limit of 0
  #  r.umsys=sapply(r.umsys,function(x) {ifelse(x<0,0.001,x)})#set lower limit of 0
  
  umsy_comp[i,5]=median(r.umsys)
  umsy_comp[i,6]=quantile(r.umsys,0.05)
  umsy_comp[i,7]=quantile(r.umsys,0.95)
  umsy_comp[i,8]=median(umsy.st)
  umsy_comp[i,9]=quantile(umsy.st,0.05)
  umsy_comp[i,10]=quantile(umsy.st,0.95)
  umsy_comp[i,11]=median(r.umsys-umsy.st)
  umsy_comp[i,12]=quantile(r.umsys-umsy.st,0.05)
  umsy_comp[i,13]=quantile(r.umsys-umsy.st,0.95)
  umsy_comp[i,14]=median((r.umsys-umsy.st)/umsy.st)
  umsy_comp[i,15]=quantile((r.umsys-umsy.st)/umsy.st,0.05)
  umsy_comp[i,16]=quantile((r.umsys-umsy.st)/umsy.st,0.95)
  umsy_comp[i,17]=median((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys)
  umsy_comp[i,18]=quantile((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys,0.05)
  umsy_comp[i,19]=quantile((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys,0.95)
}
t.umsy.diff=ifelse(umsy_comp$pct.diff.umsy*100>125,135,umsy_comp$pct.diff.umsy*100)

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts=cut(t.umsy.diff, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)
p_cols2=p_cols1[as.numeric(cuts)]
breaks=seq(-100,140,by=10)

par(mfrow=c(1,3),mar=c(5,5,3,1))
plot(r.umsy~o.umsy,data=umsy_comp,ylab=expression(paste('Recent productivity, last 5 cohorts mean ',U['MSY,j,t'],sep='')),xlab=expression(paste('Long-term average ',U['MSY,j'],sep='')),bty='l',type='n',xlim=c(0,1),ylim=c(0,1))
abline(0,1)
for(i in 1:nrow(umsy_comp)){
  lines(rep(r.umsy[i],2)~c(o.umsy.l90[i],o.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp,lwd=0.8)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(o.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp,lwd=0.8)
  lines(c(o.umsy[i],r.umsy[i])~rep(o.umsy[i],2),col=p_cols2[i],data=umsy_comp,lwd=2)
  
}
points(r.umsy~o.umsy,data=umsy_comp,bg=adjustcolor(p_cols2,alpha.f=0.4),col=adjustcolor(p_cols2,alpha.f=1),pch=21,cex=1.6)

break_cols=p_cols1[cut(breaks,c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)]

hist(t.umsy.diff,breaks=breaks,xlab='% change from average',main='',border='white',col=break_cols,xaxt='n',xlim=c(-100,140))
axis(side=1,at=c(seq(-100,100,by=25),135),labels=c('-100','','-50','','0','','+50','','+100','>+125'))
plotrix::axis.break(1,breakpos=125,style='zigzag')
mtext(side=3,'Sockeye',line=0.75,cex=1.25)

t.umsy.trend=ifelse(umsy_comp$trend.pct.umsy*100>5,6.5,umsy_comp$trend.pct.umsy*100)
breaks=seq(-5,7,by=0.5)
break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)]

hist(t.umsy.trend,breaks=breaks,xlab='Trend over series (annual % change)',main='',border='white',col=break_cols,xaxt='n',xlim=c(-5,6.5))
axis(side=1,at=c(seq(-5,5,by=2.5),6.5),labels=c('-5','-2.5','0','+2.5','+5','+>5'))
plotrix::axis.break(1,breakpos=6,style='zigzag')


##Chum####
umsy.chm=tv.chm$draws(variables='Umsy',format='draws_matrix')
st.umsy.chm=stc.chm$draws(variables='Umsy',format='draws_matrix')

#set a lower limit just above 0 to prevent negative Umsy values
umsy.chm=apply(umsy.chm,1:2,function(x) {ifelse(x<0,0.005,x)})
st.umsy.chm=apply(st.umsy.chm,1:2,function(x) {ifelse(x<0,0.005,x)})

umsy_comp=data.frame(stock=chm.info$stock.name,LME=chm.info$ocean.basin,begin=chm.info$begin,end=chm.info$end,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,o.umsy=NA,o.umsy.l90=NA,o.umsy.u90=NA,diff.umsy=NA,diff.umsy.l90=NA,diff.umsy.u90=NA,pct.diff.umsy=NA,pct.diff.umsy.l90=NA,pct.diff.umsy.u90=NA,trend.pct.umsy=NA,trend.pct.umsy.l90=NA,trend.pct.umsy.u90=NA)
for(i in 1:nrow(chm.info)){
  umsys=umsy.chm[,grepl(paste(',',i,']',sep=''),colnames(umsy.chm))]
  umsy.st=st.umsy.chm[,i]
  
  o.umsys=apply(umsys[,c(umsy_comp$begin[i]-min(umsy_comp$begin)+1):c(umsy_comp$begin[i]-min(umsy_comp$begin)+5)],1,mean)
  r.umsys=apply(umsys[,c(umsy_comp$end[i]-min(umsy_comp$begin)+1-5):c(umsy_comp$end[i]-min(umsy_comp$begin))],1,mean)
  # o.umsys=sapply(o.umsys,function(x) {ifelse(x<0,0.001,x)})#set lower limit of 0
  #  r.umsys=sapply(r.umsys,function(x) {ifelse(x<0,0.001,x)})#set lower limit of 0
  
  umsy_comp[i,5]=median(r.umsys)
  umsy_comp[i,6]=quantile(r.umsys,0.05)
  umsy_comp[i,7]=quantile(r.umsys,0.95)
  umsy_comp[i,8]=median(umsy.st)
  umsy_comp[i,9]=quantile(umsy.st,0.05)
  umsy_comp[i,10]=quantile(umsy.st,0.95)
  umsy_comp[i,11]=median(r.umsys-umsy.st)
  umsy_comp[i,12]=quantile(r.umsys-umsy.st,0.05)
  umsy_comp[i,13]=quantile(r.umsys-umsy.st,0.95)
  umsy_comp[i,14]=median((r.umsys-umsy.st)/umsy.st)
  umsy_comp[i,15]=quantile((r.umsys-umsy.st)/umsy.st,0.05)
  umsy_comp[i,16]=quantile((r.umsys-umsy.st)/umsy.st,0.95)
  umsy_comp[i,17]=median((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys)
  umsy_comp[i,18]=quantile((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys,0.05)
  umsy_comp[i,19]=quantile((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys,0.95)
}
t.umsy.diff=ifelse(umsy_comp$pct.diff.umsy*100>125,135,umsy_comp$pct.diff.umsy*100)

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts=cut(t.umsy.diff, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)
p_cols2=p_cols1[as.numeric(cuts)]
breaks=seq(-100,140,by=10)

par(mfrow=c(1,3),mar=c(5,5,3,1))
plot(r.umsy~o.umsy,data=umsy_comp,ylab=expression(paste('Recent productivity, last 5 cohorts mean ',U['MSY,j,t'],sep='')),xlab=expression(paste('Long-term average ',U['MSY,j'],sep='')),bty='l',type='n',xlim=c(0,1),ylim=c(0,1))
abline(0,1)
for(i in 1:nrow(umsy_comp)){
  lines(rep(r.umsy[i],2)~c(o.umsy.l90[i],o.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp,lwd=0.8)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(o.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp,lwd=0.8)
  lines(c(o.umsy[i],r.umsy[i])~rep(o.umsy[i],2),col=p_cols2[i],data=umsy_comp,lwd=2)
  
}
points(r.umsy~o.umsy,data=umsy_comp,bg=adjustcolor(p_cols2,alpha.f=0.4),col=adjustcolor(p_cols2,alpha.f=1),pch=21,cex=1.6)

break_cols=p_cols1[cut(breaks,c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)]

hist(t.umsy.diff,breaks=breaks,xlab='% change from average',main='',border='white',col=break_cols,xaxt='n',xlim=c(-100,140))
axis(side=1,at=c(seq(-100,100,by=25),135),labels=c('-100','','-50','','0','','+50','','+100','>+125'))
plotrix::axis.break(1,breakpos=125,style='zigzag')
mtext(side=3,'Chum',line=0.75,cex=1.25)

t.umsy.trend=ifelse(umsy_comp$trend.pct.umsy*100>5,6.5,umsy_comp$trend.pct.umsy*100)
breaks=seq(-5,7,by=0.5)
break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)]

hist(t.umsy.trend,breaks=breaks,xlab='Trend over series (annual % change)',main='',border='white',col=break_cols,xaxt='n',xlim=c(-5,6.5))
axis(side=1,at=c(seq(-5,5,by=2.5),6.5),labels=c('-5','-2.5','0','+2.5','+5','+>5'))
plotrix::axis.break(1,breakpos=6,style='zigzag')


##Coho####
umsy.cho=tv.cho$draws(variables='Umsy',format='draws_matrix')
st.umsy.cho=stc.cho$draws(variables='Umsy',format='draws_matrix')

#set a lower limit just above 0 to prevent negative Umsy values
umsy.cho=apply(umsy.cho,1:2,function(x) {ifelse(x<0,0.005,x)})
st.umsy.cho=apply(st.umsy.cho,1:2,function(x) {ifelse(x<0,0.005,x)})

umsy_comp=data.frame(stock=cho.info$stock.name,LME=cho.info$ocean.basin,begin=cho.info$begin,end=cho.info$end,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,o.umsy=NA,o.umsy.l90=NA,o.umsy.u90=NA,diff.umsy=NA,diff.umsy.l90=NA,diff.umsy.u90=NA,pct.diff.umsy=NA,pct.diff.umsy.l90=NA,pct.diff.umsy.u90=NA,trend.pct.umsy=NA,trend.pct.umsy.l90=NA,trend.pct.umsy.u90=NA)
for(i in 1:nrow(cho.info)){
  umsys=umsy.cho[,grepl(paste(',',i,']',sep=''),colnames(umsy.cho))]
  umsy.st=st.umsy.cho[,i]
  
  o.umsys=apply(umsys[,c(umsy_comp$begin[i]-min(umsy_comp$begin)+1):c(umsy_comp$begin[i]-min(umsy_comp$begin)+5)],1,mean)
  r.umsys=apply(umsys[,c(umsy_comp$end[i]-min(umsy_comp$begin)+1-5):c(umsy_comp$end[i]-min(umsy_comp$begin))],1,mean)
  # o.umsys=sapply(o.umsys,function(x) {ifelse(x<0,0.001,x)})#set lower limit of 0
  #  r.umsys=sapply(r.umsys,function(x) {ifelse(x<0,0.001,x)})#set lower limit of 0
  
  umsy_comp[i,5]=median(r.umsys)
  umsy_comp[i,6]=quantile(r.umsys,0.05)
  umsy_comp[i,7]=quantile(r.umsys,0.95)
  umsy_comp[i,8]=median(umsy.st)
  umsy_comp[i,9]=quantile(umsy.st,0.05)
  umsy_comp[i,10]=quantile(umsy.st,0.95)
  umsy_comp[i,11]=median(r.umsys-umsy.st)
  umsy_comp[i,12]=quantile(r.umsys-umsy.st,0.05)
  umsy_comp[i,13]=quantile(r.umsys-umsy.st,0.95)
  umsy_comp[i,14]=median((r.umsys-umsy.st)/umsy.st)
  umsy_comp[i,15]=quantile((r.umsys-umsy.st)/umsy.st,0.05)
  umsy_comp[i,16]=quantile((r.umsys-umsy.st)/umsy.st,0.95)
  umsy_comp[i,17]=median((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys)
  umsy_comp[i,18]=quantile((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys,0.05)
  umsy_comp[i,19]=quantile((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys,0.95)
}
t.umsy.diff=ifelse(umsy_comp$pct.diff.umsy*100>125,135,umsy_comp$pct.diff.umsy*100)

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts=cut(t.umsy.diff, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)
p_cols2=p_cols1[as.numeric(cuts)]
breaks=seq(-100,140,by=10)

par(mfrow=c(1,3),mar=c(5,5,3,1))
plot(r.umsy~o.umsy,data=umsy_comp,ylab=expression(paste('Recent productivity, last 5 cohorts mean ',U['MSY,j,t'],sep='')),xlab=expression(paste('Long-term average ',U['MSY,j'],sep='')),bty='l',type='n',xlim=c(0,1),ylim=c(0,1))
abline(0,1)
for(i in 1:nrow(umsy_comp)){
  lines(rep(r.umsy[i],2)~c(o.umsy.l90[i],o.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp,lwd=0.8)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(o.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp,lwd=0.8)
  lines(c(o.umsy[i],r.umsy[i])~rep(o.umsy[i],2),col=p_cols2[i],data=umsy_comp,lwd=2)
  
}
points(r.umsy~o.umsy,data=umsy_comp,bg=adjustcolor(p_cols2,alpha.f=0.4),col=adjustcolor(p_cols2,alpha.f=1),pch=21,cex=1.6)

break_cols=p_cols1[cut(breaks,c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)]

hist(t.umsy.diff,breaks=breaks,xlab='% change from average',main='',border='white',col=break_cols,xaxt='n',xlim=c(-100,140))
axis(side=1,at=c(seq(-100,100,by=25),135),labels=c('-100','','-50','','0','','+50','','+100','>+125'))
plotrix::axis.break(1,breakpos=125,style='zigzag')
mtext(side=3,'Coho',line=0.75,cex=1.25)

t.umsy.trend=ifelse(umsy_comp$trend.pct.umsy*100>5,6.5,umsy_comp$trend.pct.umsy*100)
breaks=seq(-5,7,by=0.5)
break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)]

hist(t.umsy.trend,breaks=breaks,xlab='Trend over series (annual % change)',main='',border='white',col=break_cols,xaxt='n',xlim=c(-5,6.5))
axis(side=1,at=c(seq(-5,5,by=2.5),6.5),labels=c('-5','-2.5','0','+2.5','+5','+>5'))
plotrix::axis.break(1,breakpos=6,style='zigzag')

##Pink - Even####
umsy.pke=tv.pke$draws(variables='Umsy',format='draws_matrix')
st.umsy.pke=stc.pke$draws(variables='Umsy',format='draws_matrix')

#set a lower limit just above 0 to prevent negative Umsy values
umsy.pke=apply(umsy.pke,1:2,function(x) {ifelse(x<0,0.005,x)})
st.umsy.pke=apply(st.umsy.pke,1:2,function(x) {ifelse(x<0,0.005,x)})

umsy_comp=data.frame(stock=pke.info$stock.name,LME=pke.info$ocean.basin,begin=pke.info$begin,end=pke.info$end,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,o.umsy=NA,o.umsy.l90=NA,o.umsy.u90=NA,diff.umsy=NA,diff.umsy.l90=NA,diff.umsy.u90=NA,pct.diff.umsy=NA,pct.diff.umsy.l90=NA,pct.diff.umsy.u90=NA,trend.pct.umsy=NA,trend.pct.umsy.l90=NA,trend.pct.umsy.u90=NA)
for(i in 1:nrow(pke.info)){
  umsys=umsy.pke[,grepl(paste(',',i,']',sep=''),colnames(umsy.pke))]
  umsy.st=st.umsy.pke[,i]
  
  o.umsys=apply(umsys[,c(umsy_comp$begin[i]-min(umsy_comp$begin)+1):c(umsy_comp$begin[i]-min(umsy_comp$begin)+5)],1,mean)
  r.umsys=apply(umsys[,c(umsy_comp$end[i]-min(umsy_comp$begin)+1-5):c(umsy_comp$end[i]-min(umsy_comp$begin))],1,mean)
  # o.umsys=sapply(o.umsys,function(x) {ifelse(x<0,0.001,x)})#set lower limit of 0
  #  r.umsys=sapply(r.umsys,function(x) {ifelse(x<0,0.001,x)})#set lower limit of 0
  
  umsy_comp[i,5]=median(r.umsys)
  umsy_comp[i,6]=quantile(r.umsys,0.05)
  umsy_comp[i,7]=quantile(r.umsys,0.95)
  umsy_comp[i,8]=median(umsy.st)
  umsy_comp[i,9]=quantile(umsy.st,0.05)
  umsy_comp[i,10]=quantile(umsy.st,0.95)
  umsy_comp[i,11]=median(r.umsys-umsy.st)
  umsy_comp[i,12]=quantile(r.umsys-umsy.st,0.05)
  umsy_comp[i,13]=quantile(r.umsys-umsy.st,0.95)
  umsy_comp[i,14]=median((r.umsys-umsy.st)/umsy.st)
  umsy_comp[i,15]=quantile((r.umsys-umsy.st)/umsy.st,0.05)
  umsy_comp[i,16]=quantile((r.umsys-umsy.st)/umsy.st,0.95)
  umsy_comp[i,17]=median((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys)
  umsy_comp[i,18]=quantile((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys,0.05)
  umsy_comp[i,19]=quantile((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys,0.95)
}
t.umsy.diff=ifelse(umsy_comp$pct.diff.umsy*100>125,135,umsy_comp$pct.diff.umsy*100)

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts=cut(t.umsy.diff, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)
p_cols2=p_cols1[as.numeric(cuts)]
breaks=seq(-100,140,by=10)

par(mfrow=c(1,3),mar=c(5,5,3,1))
plot(r.umsy~o.umsy,data=umsy_comp,ylab=expression(paste('Recent productivity, last 5 cohorts mean ',U['MSY,j,t'],sep='')),xlab=expression(paste('Long-term average ',U['MSY,j'],sep='')),bty='l',type='n',xlim=c(0,1),ylim=c(0,1))
abline(0,1)
for(i in 1:nrow(umsy_comp)){
  lines(rep(r.umsy[i],2)~c(o.umsy.l90[i],o.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp,lwd=0.8)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(o.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp,lwd=0.8)
  lines(c(o.umsy[i],r.umsy[i])~rep(o.umsy[i],2),col=p_cols2[i],data=umsy_comp,lwd=2)
  
}
points(r.umsy~o.umsy,data=umsy_comp,bg=adjustcolor(p_cols2,alpha.f=0.4),col=adjustcolor(p_cols2,alpha.f=1),pch=21,cex=1.6)

break_cols=p_cols1[cut(breaks,c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)]

hist(t.umsy.diff,breaks=breaks,xlab='% change from average',main='',border='white',col=break_cols,xaxt='n',xlim=c(-100,140))
axis(side=1,at=c(seq(-100,100,by=25),135),labels=c('-100','','-50','','0','','+50','','+100','>+125'))
plotrix::axis.break(1,breakpos=125,style='zigzag')
mtext(side=3,'Pink - even year lines',line=0.75,cex=1.25)

t.umsy.trend=ifelse(umsy_comp$trend.pct.umsy*100>5,6.5,umsy_comp$trend.pct.umsy*100)
breaks=seq(-5,7,by=0.5)
break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)]

hist(t.umsy.trend,breaks=breaks,xlab='Trend over series (annual % change)',main='',border='white',col=break_cols,xaxt='n',xlim=c(-5,6.5))
axis(side=1,at=c(seq(-5,5,by=2.5),6.5),labels=c('-5','-2.5','0','+2.5','+5','+>5'))
plotrix::axis.break(1,breakpos=6,style='zigzag')

##Pink - Odd####
umsy.pko=tv.pko$draws(variables='Umsy',format='draws_matrix')
st.umsy.pko=stc.pko$draws(variables='Umsy',format='draws_matrix')

#set a lower limit just above 0 to prevent negative Umsy values
umsy.pko=apply(umsy.pko,1:2,function(x) {ifelse(x<0,0.005,x)})
st.umsy.pko=apply(st.umsy.pko,1:2,function(x) {ifelse(x<0,0.005,x)})

umsy_comp=data.frame(stock=pko.info$stock.name,LME=pko.info$ocean.basin,begin=pko.info$begin,end=pko.info$end,r.umsy=NA,r.umsy.l90=NA,r.umsy.u90=NA,o.umsy=NA,o.umsy.l90=NA,o.umsy.u90=NA,diff.umsy=NA,diff.umsy.l90=NA,diff.umsy.u90=NA,pct.diff.umsy=NA,pct.diff.umsy.l90=NA,pct.diff.umsy.u90=NA,trend.pct.umsy=NA,trend.pct.umsy.l90=NA,trend.pct.umsy.u90=NA)
for(i in 1:nrow(pko.info)){
  umsys=umsy.pko[,grepl(paste(',',i,']',sep=''),colnames(umsy.pko))]
  umsy.st=st.umsy.pko[,i]
  
  o.umsys=apply(umsys[,c(umsy_comp$begin[i]-min(umsy_comp$begin)+1):c(umsy_comp$begin[i]-min(umsy_comp$begin)+5)],1,mean)
  r.umsys=apply(umsys[,c(umsy_comp$end[i]-min(umsy_comp$begin)+1-5):c(umsy_comp$end[i]-min(umsy_comp$begin))],1,mean)
  # o.umsys=sapply(o.umsys,function(x) {ifelse(x<0,0.001,x)})#set lower limit of 0
  #  r.umsys=sapply(r.umsys,function(x) {ifelse(x<0,0.001,x)})#set lower limit of 0
  
  umsy_comp[i,5]=median(r.umsys)
  umsy_comp[i,6]=quantile(r.umsys,0.05)
  umsy_comp[i,7]=quantile(r.umsys,0.95)
  umsy_comp[i,8]=median(umsy.st)
  umsy_comp[i,9]=quantile(umsy.st,0.05)
  umsy_comp[i,10]=quantile(umsy.st,0.95)
  umsy_comp[i,11]=median(r.umsys-umsy.st)
  umsy_comp[i,12]=quantile(r.umsys-umsy.st,0.05)
  umsy_comp[i,13]=quantile(r.umsys-umsy.st,0.95)
  umsy_comp[i,14]=median((r.umsys-umsy.st)/umsy.st)
  umsy_comp[i,15]=quantile((r.umsys-umsy.st)/umsy.st,0.05)
  umsy_comp[i,16]=quantile((r.umsys-umsy.st)/umsy.st,0.95)
  umsy_comp[i,17]=median((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys)
  umsy_comp[i,18]=quantile((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys,0.05)
  umsy_comp[i,19]=quantile((r.umsys-o.umsys)/(umsy_comp$end[i]-umsy_comp$begin[i]+1)/o.umsys,0.95)
}
t.umsy.diff=ifelse(umsy_comp$pct.diff.umsy*100>125,135,umsy_comp$pct.diff.umsy*100)

p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu');p_cols1[11]='#31006b'
cuts=cut(t.umsy.diff, breaks = c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)
p_cols2=p_cols1[as.numeric(cuts)]
breaks=seq(-100,140,by=10)

par(mfrow=c(1,3),mar=c(5,5,3,1))
plot(r.umsy~o.umsy,data=umsy_comp,ylab=expression(paste('Recent productivity, last 5 cohorts mean ',U['MSY,j,t'],sep='')),xlab=expression(paste('Long-term average ',U['MSY,j'],sep='')),bty='l',type='n',xlim=c(0,1),ylim=c(0,1))
abline(0,1)
for(i in 1:nrow(umsy_comp)){
  lines(rep(r.umsy[i],2)~c(o.umsy.l90[i],o.umsy.u90[i]),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp,lwd=0.8)
  lines(c(r.umsy.l90[i],r.umsy.u90[i])~rep(o.umsy[i],2),col=adjustcolor('darkgray',alpha.f=0.4),data=umsy_comp,lwd=0.8)
  lines(c(o.umsy[i],r.umsy[i])~rep(o.umsy[i],2),col=p_cols2[i],data=umsy_comp,lwd=2)
  
}
points(r.umsy~o.umsy,data=umsy_comp,bg=adjustcolor(p_cols2,alpha.f=0.4),col=adjustcolor(p_cols2,alpha.f=1),pch=21,cex=1.6)

break_cols=p_cols1[cut(breaks,c(seq(-100,0,length.out=6),seq(20,100,length.out=5),140),right=F)]

hist(t.umsy.diff,breaks=breaks,xlab='% change from average',main='',border='white',col=break_cols,xaxt='n',xlim=c(-100,140))
axis(side=1,at=c(seq(-100,100,by=25),135),labels=c('-100','','-50','','0','','+50','','+100','>+125'))
plotrix::axis.break(1,breakpos=125,style='zigzag')
mtext(side=3,'Pink - odd year lines',line=0.75,cex=1.25)

t.umsy.trend=ifelse(umsy_comp$trend.pct.umsy*100>5,6.5,umsy_comp$trend.pct.umsy*100)
breaks=seq(-5,7,by=0.5)
break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-5,0,length.out=6),seq(1,5,length.out=5),7),right=F)]

hist(t.umsy.trend,breaks=breaks,xlab='Trend over series (annual % change)',main='',border='white',col=break_cols,xaxt='n',xlim=c(-5,6.5))
axis(side=1,at=c(seq(-5,5,by=2.5),6.5),labels=c('-5','-2.5','0','+2.5','+5','+>5'))
plotrix::axis.break(1,breakpos=6,style='zigzag')

#Covariance patterns####


calc_dist<- function(lat1, lon1, lat2, lon2) {
  # Coordinates are in degrees, so we use the distVincentySphere function
  distance <- geosphere::distVincentySphere(c(lon1, lat1), c(lon2, lat2))
  return(distance)  # Distance in meters
}

chk.info$l


p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu')
p_change=exp(lt_prod_chi$r.alpha)/exp(lt_prod_chi$alpha)-1
cuts=cut(p_change, breaks = c(seq(-1,1,length.out=10),3.5))
p_cols2=p_cols1[as.numeric(cuts)]
breaks=hist(p_change,breaks=30)$breaks

par(mfrow=c(1,2))
plot(exp(r.alpha)~exp(alpha),data=lt_prod,ylab='Recent (last 10y) productivity (max. recruits/spawner)',xlab='Long-term average productivity (max. recruits/spawner)',bty='l',type='n',xlim=c(0,28),ylim=c(0,28))
abline(0,1)
for(i in 1:nrow(stock_info_filtered)){
  lines(c(exp(alpha[i]),exp(r.alpha[i]))~rep(exp(alpha[i]),2),col=p_cols2[i],data=lt_prod_chi)
}
#points(exp(r.alpha)~exp(alpha),data=lt_prod,bg=adjustcolor('darkgray',alpha.f=0.2),col=adjustcolor('black',alpha.f=0.8),pch=21,cex=1)
points(exp(r.alpha)~exp(alpha),data=lt_prod_chi,bg=adjustcolor(sp_cols[1],alpha.f=0.2),col=adjustcolor(sp_cols[1],alpha.f=0.8),pch=21,cex=1)
#points(exp(r.alpha)~exp(alpha),data=lt_prod_chu,bg=adjustcolor(sp_cols[2],alpha.f=0.2),col=adjustcolor(sp_cols[2],alpha.f=0.8),pch=21,cex=1.2)
#points(exp(r.alpha)~exp(alpha),data=lt_prod_coh,bg=adjustcolor(sp_cols[3],alpha.f=0.2),col=adjustcolor(sp_cols[3],alpha.f=0.8),pch=21,cex=1.2)
#points(exp(r.alpha)~exp(alpha),data=lt_prod_pi,bg=adjustcolor(sp_cols[4],alpha.f=0.2),col=adjustcolor(sp_cols[4],alpha.f=0.8),pch=21,cex=1.2)
#points(exp(r.alpha)~exp(alpha),data=lt_prod_soc,bg=adjustcolor(sp_cols[5],alpha.f=0.2),col=adjustcolor(sp_cols[5],alpha.f=0.8),pch=21,cex=1.2)

break_cols=p_cols1[cut(as.numeric(breaks),breaks = c(seq(-1,1,length.out=10),3.5))]

hist(p_change,breaks=30,xlab='Proportional change relative to long-term average',main='',col=break_cols,xaxt='n')
axis(side=1,at=c(seq(-1,3,by=0.25)),labels=c('-1','','','','0','','','','1','','','','2','','','','3'))

summary(p_change)

#Individual Spawner-recruit model fits####
##Chinook####
###West Coast####

###SE Alaska####
###Gulf of Alaska####
###Bering Sea####

