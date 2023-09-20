#multipop S-R curves
#stock characteristics
library(here);library(dplyr);library(rstan);library(ggplot2)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_may2023.csv'))
stock_info<- read.csv(here('data','filtered datasets','all_stocks_info_may2023.csv'))

library(samEst)

#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')
source(here('code','util_func.R'))

###Load in data####
#Remove stocks with less than 15 years of recruitment data
stock_info_filtered=subset(stock_info,n.years>=16) #242 stocks
stock_info_filtered$stock.name=gsub('/','_',stock_info_filtered$stock.name)
stock_info_filtered$stock.name=gsub('&','and',stock_info_filtered$stock.name)

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock.id)) #242
stock_info_filtered$stock.id2=seq(1:nrow(stock_info_filtered))
if(any(stock_dat2$spawners==0)){stock_dat2$spawners=stock_dat2$spawners+1}
if(any(stock_dat2$recruits==0)){stock_dat2$recruits=stock_dat2$recruits+1}

stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)
stock_dat2=stock_dat2[complete.cases(stock_dat2$logR_S),]



##Individual fits via TMB####
stock_info_filtered$state2=stock_info_filtered$state
stock_info_filtered$state2=ifelse(stock_info_filtered$state=='OR'|stock_info_filtered$state=='WA','OR/WA',stock_info_filtered$state2)

chi_stocks=subset(stock_info_filtered,species=='Chinook')
chi_dat=subset(stock_dat2,stock.id %in% chi_stocks$stock.id)
length(unique(chi_dat$stock.id))

coh_stocks=subset(stock_info_filtered,species=='Coho')
coh_dat=subset(stock_dat2,stock.id %in% coh_stocks$stock.id)
length(unique(coh_dat$stock.id))

soc_stocks=subset(stock_info_filtered,species=='Sockeye')
soc_dat=subset(stock_dat2,stock.id %in% soc_stocks$stock.id)
length(unique(soc_dat$stock.id))

pink_stocks1=subset(stock_info_filtered,species=='Pink-Even')
pink_dat1=subset(stock_dat2,stock.id %in% pink_stocks$stock.id)
length(unique(pink_dat1$stock.id))

pink_stocks2=subset(stock_info_filtered,species=='Pink-Even')
pink_dat2=subset(stock_dat2,stock.id %in% pink_stocks$stock.id)
length(unique(pink_dat1$stock.id))

#Colours####
sp_cols=cbind(c('#6f9ba7',
                '#316024',
                '#295f8b',
                '#a98b1a',
                '#975a0d',
                '#811b0a'),c('Chinook','Chum','Coho','Pink-Even','Pink-Odd','Sockeye'))

mod_col=cbind(c("azure4", "azure3", "cyan4","cadetblue3","deepskyblue4",'seagreen','olivedrab','darkgreen'),seq=1:8)

mc_col=cbind(c("azure3","deepskyblue4",'darkgreen'),c('static','dynamic','regime'))


mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

#Chinook####
chi_stocks=subset(stock_info_filtered,species=='Chinook')
chi_dat=subset(stock_dat2,stock.id %in% chi_stocks$stock.id)
length(unique(chi_dat$stock.id))

alpha_mat=matrix(nrow=nrow(chi_stocks),ncol=(max(chi_stocks$end)-min(chi_stocks$begin))+1)
colnames(alpha_mat)=seq(min(chi_stocks$begin),max(chi_stocks$end))


for(u in 1:nrow(chi_stocks)){
  s<- subset(chi_dat,stock.id==chi_stocks$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  
  TMBtva <- ricker_rw_TMB(data=df,tv.par='a')
  
  m=match(s$broodyear,colnames(alpha_mat))
  alpha_mat[u,m]=TMBtva$alpha
}

pdf('chin_all.pdf',height=8,width=10)
par(mar=c(4,5,1,1),oma=c(0,0,0,0))
glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
plot(exp(glob_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(chi_stocks)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkgray',alpha.f=0.3),lwd=2)
}
lines(exp(glob_a)~as.numeric(colnames(alpha_mat)),lwd=4,col=sp_cols[1])
dev.off()

state=chi_stocks$state2
AK_mat=alpha_mat[state=='AK',]
AK_a=apply(AK_mat,2,mean,na.rm=TRUE)
BC_mat=alpha_mat[state=='BC',]
BC_a=apply(BC_mat,2,mean,na.rm=TRUE)
WA_mat=alpha_mat[state=='OR/WA',]
WA_a=apply(WA_mat,2,mean,na.rm=TRUE)

pdf('chin_bystate.pdf',height=8,width=10)
par(mar=c(4,5,1,1),oma=c(0,0,0,0))
plot(exp(WA_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(AK_mat)){
  alph=AK_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('steelblue',alpha.f=0.3),lwd=2)
}
for(u in 1:nrow(BC_mat)){
  alph=BC_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('goldenrod',alpha.f=0.3),lwd=2)
}
for(u in 1:nrow(WA_mat)){
  alph=WA_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('red4',alpha.f=0.3),lwd=2)
}
lines(exp(WA_a)~as.numeric(colnames(alpha_mat)),lwd=4,col='red4')
lines(exp(BC_a)~as.numeric(colnames(alpha_mat)),lwd=4,col='goldenrod')
lines(exp(AK_a)~as.numeric(colnames(alpha_mat)),lwd=4,col='steelblue')
dev.off()

##Chum####
chu_stocks=subset(stock_info_filtered,species=='Chum')
chu_dat=subset(stock_dat2,stock.id %in% chu_stocks$stock.id)
length(unique(chu_dat$stock.id))

alpha_mat=matrix(nrow=nrow(chu_stocks),ncol=(max(chu_stocks$end)-min(chu_stocks$begin))+1)
colnames(alpha_mat)=seq(min(chu_stocks$begin),max(chu_stocks$end))

for(u in 1:nrow(chu_stocks)){
  s<- subset(chu_dat,stock.id==chu_stocks$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  
  TMBtva <- ricker_rw_TMB(data=df,tv.par='a')
  
  m=match(s$broodyear,colnames(alpha_mat))
  alpha_mat[u,m]=TMBtva$alpha
}

par(mar=c(4,4.5,1,1),oma=c(0.3,0.3,0,0))
glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
plot(exp(glob_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(chi_stocks)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkgray',alpha.f=0.3),lwd=2)
}
lines(exp(glob_a)~as.numeric(colnames(alpha_mat)),lwd=4,col=sp_cols[2])

pdf('chum_all.pdf',height=8,width=10)
par(mar=c(4,5,2,2),oma=c(0,0,0,0))
glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
plot(exp(glob_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(chu_stocks)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkgray',alpha.f=0.3),lwd=2)
}
lines(exp(glob_a)~as.numeric(colnames(alpha_mat)),lwd=4,col=sp_cols[2])
dev.off()

state=chu_stocks$state2
AK_mat=alpha_mat[state=='AK',]
AK_a=apply(AK_mat,2,mean,na.rm=TRUE)
BC_mat=alpha_mat[state=='BC',]
BC_a=apply(BC_mat,2,mean,na.rm=TRUE)
WA_mat=alpha_mat[state=='OR/WA',]
WA_a=apply(WA_mat,2,mean,na.rm=TRUE)

pdf('chum_bystate.pdf',height=8,width=10)
par(mar=c(4,5,2,2),oma=c(0,0,0,0))
plot(exp(WA_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(AK_mat)){
  alph=AK_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('steelblue',alpha.f=0.3),lwd=2)
}
for(u in 1:nrow(BC_mat)){
  alph=BC_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('goldenrod',alpha.f=0.3),lwd=2)
}
for(u in 1:nrow(WA_mat)){
  alph=WA_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('red4',alpha.f=0.3),lwd=2)
}
lines(exp(WA_a)~as.numeric(colnames(alpha_mat)),lwd=4,col='red4')
lines(exp(BC_a)~as.numeric(colnames(alpha_mat)),lwd=4,col='goldenrod')
lines(exp(AK_a)~as.numeric(colnames(alpha_mat)),lwd=4,col='steelblue')
dev.off()



#Coho
coh_stocks=subset(stock_info_filtered,species=='Coho')
coh_dat=subset(stock_dat2,stock.id %in% coh_stocks$stock.id)
length(unique(coh_dat$stock.id))

alpha_mat=matrix(nrow=nrow(coh_stocks),ncol=(max(coh_stocks$end)-min(coh_stocks$begin))+1)
colnames(alpha_mat)=seq(min(coh_stocks$begin),max(coh_stocks$end))

for(u in 1:nrow(coh_stocks)){
  s<- subset(coh_dat,stock.id==coh_stocks$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  
  TMBtva <- ricker_rw_TMB(data=df,tv.par='a')
  
  m=match(s$broodyear,colnames(alpha_mat))
  alpha_mat[u,m]=TMBtva$alpha
}

par(mar=c(4,4.5,1,1),oma=c(0.3,0.3,0,0))
glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
plot(exp(glob_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(coh_stocks)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkgray',alpha.f=0.3),lwd=2)
}
lines(exp(glob_a)~as.numeric(colnames(alpha_mat)),lwd=4,col=sp_cols[3])


#Pink####
pin_stocks_even=subset(stock_info_filtered,species=='Pink-Even')
pin_stocks_odd=subset(stock_info_filtered,species=='Pink-Odd')
pin_dat_even=subset(stock_dat2,stock.id %in% pin_stocks_even$stock.id)
pin_dat_odd=subset(stock_dat2,stock.id %in% pin_stocks_odd$stock.id)

length(unique(pin_dat_even$stock.id))
length(unique(pin_dat_odd$stock.id))

alpha_mat=matrix(nrow=nrow(pin_stocks_even),ncol=(max(pin_stocks_even$end)-min(pin_stocks_even$begin))+1)
colnames(alpha_mat)=seq(min(pin_stocks_even$begin),max(pin_stocks_even$end))
for(u in 1:nrow(pin_stocks_even)){
  s<- subset(pin_dat_even,stock.id==pin_stocks_even$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  
  TMBtva <- ricker_rw_TMB(data=df,tv.par='a')
  
  m=match(s$broodyear,colnames(alpha_mat))
  alpha_mat[u,m]=TMBtva$alpha
}

pdf('pink_all.pdf',height=8,width=10)
par(mar=c(4,5,2,2),oma=c(0,0,0,0))
glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
glob_a=glob_a[complete.cases(glob_a)]
plot(exp(glob_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(pin_stocks_even)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkgray',alpha.f=0.3),lwd=2)
}
lines(exp(glob_a)~as.numeric(names(glob_a)),lwd=4,col=sp_cols[4])
dev.off()

alpha_mat=matrix(nrow=nrow(pin_stocks_odd),ncol=(max(pin_stocks_odd$end)-min(pin_stocks_odd$begin))+1)
colnames(alpha_mat)=seq(min(pin_stocks_odd$begin),max(pin_stocks_odd$end))
for(u in 1:nrow(pin_stocks_odd)){
  s<- subset(pin_dat_odd,stock.id==pin_stocks_odd$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  
  TMBtva <- ricker_rw_TMB(data=df,tv.par='a')
  
  m=match(s$broodyear,colnames(alpha_mat))
  alpha_mat[u,m]=TMBtva$alpha
}

par(mar=c(4,5,2,2),oma=c(0,0,0,0))
glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
glob_a=glob_a[complete.cases(glob_a)]
plot(exp(glob_a)~as.numeric(names(glob_a)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(pin_stocks_odd)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkgray',alpha.f=0.3),lwd=2)
}
lines(exp(glob_a)~as.numeric(names(glob_a)),lwd=4,col=sp_cols[5])
dev.off()



state=pin_stocks$state2
AK_mat=alpha_mat[state=='AK',]
AK_a=apply(AK_mat,2,mean,na.rm=TRUE)
BC_mat=alpha_mat[state=='BC',]
BC_a=apply(BC_mat,2,mean,na.rm=TRUE)
WA_mat=alpha_mat[state=='OR/WA',]
WA_a=apply(WA_mat,2,mean,na.rm=TRUE)

pdf('pink_bystate.pdf',height=8,width=10)
par(mar=c(4,5,2,2),oma=c(0,0,0,0))
plot(exp(BC_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(AK_mat)){
  alph=AK_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('steelblue',alpha.f=0.3),lwd=2)
}
for(u in 1:nrow(BC_mat)){
  alph=BC_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('goldenrod',alpha.f=0.3),lwd=2)
}
for(u in 1:nrow(WA_mat)){
  alph=WA_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('red4',alpha.f=0.3),lwd=2)
}
WA_a=WA_a[complete.cases(WA_a)]
lines(exp(WA_a)~as.numeric(names(WA_a)),lwd=4,col='red4')
lines(exp(BC_a)~as.numeric(colnames(alpha_mat)),lwd=4,col='goldenrod')
lines(exp(AK_a)~as.numeric(colnames(alpha_mat)),lwd=4,col='steelblue')
dev.off()

par(mar=c(4,4.5,1,1),oma=c(0.3,0.3,0,0))
glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
plot(exp(glob_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(pin_stocks)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkgray',alpha.f=0.3),lwd=2)
}
lines(exp(glob_a)~as.numeric(colnames(alpha_mat)),lwd=4,col=sp_cols[4])

###Odd vs Even####
even_mat=alpha_mat[pin_stocks$odd==0,]
even_a=apply(even_mat,2,mean,na.rm=TRUE)
odd_mat=alpha_mat[pin_stocks$odd==1,]
odd_a=apply(odd_mat,2,mean,na.rm=TRUE)

odd_a=apply(alpha_mat_odd,2,mean,na.rm=TRUE)
even_a=apply(alpha_mat_even,2,mean,na.rm=TRUE)


pdf('pink_even_odd.pdf',height=8,width=10)
par(mar=c(4,5,2,2),oma=c(0,0,0,0))

plot(exp(even_a)~as.numeric(colnames(alpha_mat_even)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(alpha_mat_odd)){
  alph=alpha_mat_odd[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkorange4',alpha.f=0.3),lwd=2)
}
for(u in 1:nrow(alpha_mat_even)){
  alph=alpha_mat_even[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkblue',alpha.f=0.3),lwd=2)
}

odd_a=odd_a[complete.cases(odd_a)]
even_a=odd_a[complete.cases(even_a)]

lines(exp(odd_a)~as.numeric(names(odd_a)),lwd=4,col='darkorange4')
lines(exp(even_a)~as.numeric(names(even_a)),lwd=4,col='darkblue')
dev.off()

#Sockeye####
soc_stocks=subset(stock_info_filtered,species=='Sockeye')
soc_dat=subset(stock_dat2,stock.id %in% soc_stocks$stock.id)
length(unique(soc_dat$stock.id))

alpha_mat=matrix(nrow=nrow(soc_stocks),ncol=(max(soc_stocks$end)-min(soc_stocks$begin))+1)
colnames(alpha_mat)=seq(min(soc_stocks$begin),max(soc_stocks$end))

for(u in 1:nrow(soc_stocks)){
  s<- subset(soc_dat,stock.id==soc_stocks$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  
  TMBtva <- ricker_rw_TMB(data=df,tv.par='a')
  
  m=match(s$broodyear,colnames(alpha_mat))
  alpha_mat[u,m]=TMBtva$alpha
}


pdf('sock_all.pdf',height=8,width=10)
par(mar=c(4,5,2,2),oma=c(0,0,0,0))
glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
plot(exp(glob_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(pin_stocks)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkgray',alpha.f=0.3),lwd=2)
}
lines(exp(glob_a)~as.numeric(colnames(alpha_mat)),lwd=4,col=sp_cols[5])
dev.off()

state=soc_stocks$state2
AK_mat=alpha_mat[state=='AK',]
AK_a=apply(AK_mat,2,mean,na.rm=TRUE)
BC_mat=alpha_mat[state=='BC',]
BC_a=apply(BC_mat,2,mean,na.rm=TRUE)

pdf('soc_bystate.pdf',height=8,width=10)
par(mar=c(4,5,2,2),oma=c(0,0,0,0))
plot(exp(BC_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(AK_mat)){
  alph=AK_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('steelblue',alpha.f=0.3),lwd=2)
}
for(u in 1:nrow(BC_mat)){
  alph=BC_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('goldenrod',alpha.f=0.3),lwd=2)
}
alph=alpha_mat[state=='OR/WA',];alph=alph[is.na(alph)==F]
lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('red4',alpha.f=0.3),lwd=2)

lines(exp(BC_a)~as.numeric(colnames(alpha_mat)),lwd=4,col='goldenrod')
lines(exp(AK_a)~as.numeric(colnames(alpha_mat)),lwd=4,col='steelblue')
dev.off()


par(mar=c(4,4.5,1,1),oma=c(0.3,0.3,0,0))
glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
plot(exp(glob_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(soc_stocks)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkgray',alpha.f=0.3),lwd=2)
}
lines(exp(glob_a)~as.numeric(colnames(alpha_mat)),lwd=4,col=sp_cols[5])


soc_stocks=subset(stock_info_filtered,species=='Sockeye'&state2=='AK')
soc_dat=subset(stock_dat2,stock.id %in% soc_stocks$stock.id)
length(unique(soc_dat$stock.id))

alpha_mat=matrix(nrow=nrow(soc_stocks),ncol=(max(soc_stocks$end)-min(soc_stocks$begin))+1)
colnames(alpha_mat)=seq(min(soc_stocks$begin),max(soc_stocks$end))

for(u in 1:nrow(soc_stocks)){
  s<- subset(soc_dat,stock.id==soc_stocks$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  
  TMBtva <- ricker_rw_TMB(data=df,tv.par='a')
  
  m=match(s$broodyear,colnames(alpha_mat))
  alpha_mat[u,m]=TMBtva$alpha
}

par(mar=c(4,4.5,1,1),oma=c(0.3,0.3,0,0))
glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
plot(exp(glob_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(soc_stocks)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkgray',alpha.f=0.3),lwd=2)
}
lines(exp(glob_a)~as.numeric(colnames(alpha_mat)),lwd=4,col=sp_cols[5])

soc_stocks=subset(stock_info_filtered,species=='Sockeye'&state2=='BC')
soc_dat=subset(stock_dat2,stock.id %in% soc_stocks$stock.id)
length(unique(soc_dat$stock.id))

alpha_mat=matrix(nrow=nrow(soc_stocks),ncol=(max(soc_stocks$end)-min(soc_stocks$begin))+1)
colnames(alpha_mat)=seq(min(soc_stocks$begin),max(soc_stocks$end))

for(u in 1:nrow(soc_stocks)){
  s<- subset(soc_dat,stock.id==soc_stocks$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  
  TMBtva <- ricker_rw_TMB(data=df,tv.par='a')
  
  m=match(s$broodyear,colnames(alpha_mat))
  alpha_mat[u,m]=TMBtva$alpha
}

par(mar=c(4,4.5,1,1),oma=c(0.3,0.3,0,0))
glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
plot(exp(glob_a)~as.numeric(colnames(alpha_mat)),type='n',ylim=c(0,20),xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)',cex.lab=1.5,cex.axis=1.5)
for(u in 1:nrow(soc_stocks)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(exp(alph)~as.numeric(names(alph)),col=adjustcolor('darkgray',alpha.f=0.3),lwd=2)
}
lines(exp(glob_a)~as.numeric(colnames(alpha_mat)),lwd=4,col=sp_cols[5])

soc_stocks=subset(stock_info_filtered,species=='Sockeye'&state2=='BC')
soc_dat=subset(stock_dat2,stock.id %in% soc_stocks$stock.id)
length(unique(soc_dat$stock.id))

alpha_mat=matrix(nrow=nrow(soc_stocks),ncol=(max(soc_stocks$end)-min(soc_stocks$begin))+1)
colnames(alpha_mat)=seq(min(soc_stocks$begin),max(soc_stocks$end))

for(u in 1:nrow(soc_stocks)){
  s<- subset(soc_dat,stock.id==soc_stocks$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  
  TMBtva <- ricker_rw_TMB(data=df,tv.par='a')
  
  m=match(s$broodyear,colnames(alpha_mat))
  alpha_mat[u,m]=TMBtva$alpha
}

glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
plot(glob_a~as.numeric(colnames(alpha_mat)),type='b',ylim=c(-2,4))
for(u in 1:nrow(soc_stocks)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(alph~as.numeric(names(alph)),col=adjustcolor('black',alpha.f=0.2))
}

soc_stocks=subset(stock_info_filtered,species=='Sockeye'&state2=='OR/WA')
soc_dat=subset(stock_dat2,stock.id %in% soc_stocks$stock.id)
length(unique(soc_dat$stock.id))

alpha_mat=matrix(nrow=nrow(soc_stocks),ncol=(max(soc_stocks$end)-min(soc_stocks$begin))+1)
colnames(alpha_mat)=seq(min(soc_stocks$begin),max(soc_stocks$end))

for(u in 1:nrow(soc_stocks)){
  s<- subset(soc_dat,stock.id==soc_stocks$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  
  TMBtva <- ricker_rw_TMB(data=df,tv.par='a')
  
  m=match(s$broodyear,colnames(alpha_mat))
  alpha_mat[u,m]=TMBtva$alpha
}

glob_a=apply(alpha_mat,2,mean,na.rm=TRUE)
plot(glob_a~as.numeric(colnames(alpha_mat)),type='b',ylim=c(-2,4))
for(u in 1:nrow(soc_stocks)){
  alph=alpha_mat[u,];alph=alph[is.na(alph)==F]
  lines(alph~as.numeric(names(alph)),col=adjustcolor('black',alpha.f=0.2))
}






#Chinook
chi_stocks=subset(stock_info_filtered,species=='Chinook')
chi_dat=subset(stock_dat2,stock.id %in% chi_stocks$stock.id)
length(unique(chi_dat$stock.id))

stock_year=expand.grid(unique(chi_dat$stock),unique(chi_dat$broodyear))
stock_year=stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")

chi_dat$stock_yr=stock_year[match(paste(chi_dat$stock,chi_dat$broodyear,sep='_'),stock_year[,3]),3]

X_s=make_design_matrix(chi_dat$spawners,grp=chi_dat$stock)

test_run3 = rstan::sampling(test2, 
                            data = list(N=nrow(chi_dat),
                                        L=max(chi_dat$broodyear)-min(chi_dat$broodyear)+1,
                                        J_i=as.numeric(factor(chi_dat$stock)),
                                        J_ii=as.numeric(factor(chi_dat$stock_yr)),
                                        J=length(unique(chi_dat$stock)),
                                        R_S=chi_dat$logR_S,
                                        S=X_s),
                            control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 100, chains = 4, iter = 300,init=1)
d=extract(test_run3)

summary(test_run3,pars=c('sigma'))

plot(apply(d$log_a,2,median),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:16){
  lines(apply(d$log_a_t[,,i],2,median),col=adjustcolor('black',alpha.f=0.2))
}



##Chinook
#Fit individually
chi_stocks=subset(stock_info_filtered,species=='Chinook')
chi_dat=subset(stock_dat2,stock.id %in% chi_stocks$stock.id)
length(unique(chi_dat$stock.id))

m3f=samEst::sr_mod(type='rw',par='a',lfo=F)

prod_mat=matrix(ncol=max(chi_stocks$end)-min(chi_stocks$begin)+1,nrow=nrow(chi_stocks))
colnames(prod_mat)=seq(min(chi_stocks$begin),max(chi_stocks$end))
rownames(prod_mat)=chi_stocks$stock.name
for(i in 1:nrow(chi_stocks)){
  s=subset(chi_dat,stock.id==chi_stocks$stock.id[i])
  
  f3 = rstan::sampling(m3f, 
                       data = list(N=nrow(s),
                                   L=max(s$broodyear)-min(s$broodyear)+1,
                                   ii=s$broodyear-min(s$broodyear)+1,
                                   R_S=s$logR_S,
                                   S=s$spawners),
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 600)
  
  d_3=rstan::extract(f3)
  med_log_a=apply(d_3$log_a,2,median)
  names(med_log_a)=seq(chi_stocks$begin[i],chi_stocks$end[i])
  
  prod_mat[i,c(match(names(med_log_a),colnames(prod_mat)))]=med_log_a
}

prod_avg=apply(prod_mat,2,mean,na.rm=TRUE)

plot(prod_avg~colnames(prod_mat),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:16){
  st_pr=prod_mat[i,];st_pr=st_pr[complete.cases(st_pr)]
  lines(st_pr~names(st_pr),col=adjustcolor('black',alpha.f=0.2))
  points(st_pr~names(st_pr),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}



#By species
mc="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int J_i[N];//index of stocks
  int J_ii[N];//index of stock-year combination
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  }
parameters{
  real log_a0; //inital productivity (global)
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<upper = 0>[J] log_b; // rate capacity - fixed in this

 //variance components  
 vector<lower = 0>[J] sigma;

//MVN parameters  
  cholesky_factor_corr[J] Lcorr;
  vector[L-1] z_dev; //average deviation in stock productivity among stocks
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  real<lower = 0> sigma_a; //variance in average productivity among years
  vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years
  
}
transformed parameters{
  vector[J] b; //capacity rate
  //global state process
  vector[L] log_a; //global (ie. average) productivity among stocks over time
  //stock state process
  matrix[L,J] log_a_s; //stock deviation from global log_a over time
  matrix[L-1,J] a_dev; //stock-level deviations in year-to-year productivity
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  
  b=exp(log_b);
  
  //initial productivities
  log_a[1] = log_a0; //average productivity among stocks at t = 1
  log_a_s[1,] = to_row_vector(log_a0_s); //stock deviation in productivity at t =1
  log_a_t[1,] = log_a[1]+to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  
  for(t in 1:L-1){
   a_dev[t,] = (diag_pre_multiply(sigma_a_s, Lcorr) * to_vector(z_dev_s[t,]))';
  }
  
  for(t in 2:L){
    log_a[t] = log_a[t-1] + z_dev[t-1]*sigma_a; //global prod. random walk
    log_a_s[t,] = log_a_s[t-1,] + a_dev[t-1,]; //stock-specific deviation in random walk
    log_a_t[t,] = log_a[t] + log_a_s[t,]; //final estimate of stock prod. (global + stock)
  }
  
}  
model{
  //priors
  log_a0 ~ gamma(3,1.5); //initial globalproductivity - wide prior
  log_a0_s ~ normal(0,1); //stock-level deviations around the global productivity
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  z_dev ~ std_normal(); //standardized stock-level deviances in prod
  to_vector(z_dev_s) ~ std_normal(); //global deviation in prod.
  
  Lcorr ~ lkj_corr_cholesky(2.0); // prior for correlation of process deviances
 
  //variance terms
  target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior  
  target += normal_lpdf(sigma_a_s | 0, 0.5) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior
  target += normal_lpdf(sigma_a| 0, 0.5) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior
 

  R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities {
  vector[N] log_lik;

  for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
  corr_matrix[J] Cor_t = multiply_lower_tri_self_transpose(Lcorr);
}


"


#Colours
sp_cols=cbind(c('#6f9ba7',
                '#316024',
                '#295f8b',
                '#a2450c',
                '#811b0a'),c('Chinook','Chum','Coho','Pink','Sockeye'))

mod_col=cbind(c("azure4", "azure3", "cyan4","cadetblue3","deepskyblue4",'seagreen','olivedrab','darkgreen'),seq=1:8)

mc_col=cbind(c("azure3","deepskyblue4",'darkgreen'),c('static','dynamic','regime'))


mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)


#Chinook####

chi_stocks=subset(stock_info_filtered,species=='Chinook')
chi_dat=subset(stock_dat2,stock.id %in% chi_stocks$stock.id)
length(unique(chi_dat$stock.id))

test=rstan::stan_model(model_code=multi_avg)


stock_year=expand.grid(unique(chi_dat$stock),unique(chi_dat$broodyear))
stock_year=stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")

chi_dat$stock_yr=match(paste(chi_dat$stock,chi_dat$broodyear,sep='_'),stock_year[,3])

X_s=make_design_matrix(chi_dat$spawners,grp=chi_dat$stock)

chi_indfit_avg = rstan::sampling(test, 
                            data = list(N=nrow(chi_dat),
                                        init=4,
                                        L=max(chi_dat$broodyear)-min(chi_dat$broodyear)+1,
                                        J_i=as.numeric(factor(chi_dat$stock)),
                                        J_ii=chi_dat$stock_yr,
                                        J=length(unique(chi_dat$stock)),
                                        R_S=chi_dat$logR_S,
                                        S=X_s),
                            control = list(adapt_delta = 0.99,max_treedepth=20),init=2, warmup = 100, chains = 1, iter = 300)

coh_stocks=subset(stock_info_filtered,species=='Coho')
coh_dat=subset(stock_dat2,stock.id %in% coh_stocks$stock.id)
length(unique(coh_dat$stock.id))

test=rstan::stan_model(model_code=multi_avg)


stock_year=expand.grid(unique(coh_dat$stock),unique(coh_dat$broodyear))
stock_year=stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")

coh_dat$stock_yr=match(paste(coh_dat$stock,coh_dat$broodyear,sep='_'),stock_year[,3])

X_s=make_design_matrix(coh_dat$spawners,grp=coh_dat$stock)

coh_indfit_avg = rstan::sampling(test, 
                                 data = list(N=nrow(coh_dat),
                                             init=0,
                                             L=max(coh_dat$broodyear)-min(coh_dat$broodyear)+1,
                                             J_i=as.numeric(factor(coh_dat$stock)),
                                             J_ii=coh_dat$stock_yr,
                                             J=length(unique(coh_dat$stock)),
                                             R_S=coh_dat$logR_S,
                                             S=X_s),
                                 control = list(adapt_delta = 0.99,max_treedepth=20),init=2, warmup = 100, chains = 4, iter = 300)

d=extract(coh_indfit_avg)

plot(apply(d$log_a_g,2,median)~seq(min(coh_stocks$begin),max(coh_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2,xlab='year')
for(i in 1:nrow(coh_stocks)){
  s=coh_stocks$begin[i]-min(coh_stocks$begin)+1
  e=coh_stocks$end[i]-min(coh_stocks$begin)+1
  length=max(coh_stocks$end)-min(coh_stocks$begin)+1
  
  lines(apply(d$log_a_t[,,i][,s:e],2,median)~seq(coh_stocks$begin[i],coh_stocks$end[i]),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i][,s:e],2,median)~seq(coh_stocks$begin[i],coh_stocks$end[i]),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}


soc_stocks=subset(stock_info_filtered,species=='Sockeye')
soc_dat=subset(stock_dat2,stock.id %in% soc_stocks$stock.id)
length(unique(soc_dat$stock.id))

stock_year=expand.grid(unique(soc_dat$stock),unique(soc_dat$broodyear))
stock_year=stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")

soc_dat$stock_yr=match(paste(soc_dat$stock,soc_dat$broodyear,sep='_'),stock_year[,3])

X_s=make_design_matrix(soc_dat$spawners,grp=soc_dat$stock)

soc_indfit_avg = rstan::sampling(test, 
                                 data = list(N=nrow(soc_dat),
                                             init=2,
                                             L=max(soc_dat$broodyear)-min(soc_dat$broodyear)+1,
                                             J_i=as.numeric(factor(soc_dat$stock)),
                                             J_ii=soc_dat$stock_yr,
                                             J=length(unique(soc_dat$stock)),
                                             R_S=soc_dat$logR_S,
                                             S=X_s),
                                 control = list(adapt_delta = 0.99,max_treedepth=20),init=2, warmup = 100, chains = 4, iter = 300)

d=extract(soc_indfit_avg)

plot(apply(d$log_a_g,2,median)~seq(min(soc_stocks$begin),max(soc_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2,xlab='year')
for(i in 1:nrow(soc_stocks)){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(soc_stocks$begin),max(soc_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(soc_stocks$begin),max(soc_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}

pink_stocks=subset(stock_info_filtered,species=='Pink')
pink_dat=subset(stock_dat2,stock.id %in% pink_stocks$stock.id)
length(unique(pink_dat$stock.id))

stock_year=expand.grid(unique(pink_dat$stock),unique(pink_dat$broodyear))
stock_year=stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")

pink_dat$stock_yr=match(paste(pink_dat$stock,pink_dat$broodyear,sep='_'),stock_year[,3])

X_s=make_design_matrix(pink_dat$spawners,grp=pink_dat$stock)

pink_indfit_avg = rstan::sampling(test, 
                                 data = list(N=nrow(pink_dat),
                                             init=2,
                                             L=max(pink_dat$broodyear)-min(pink_dat$broodyear)+1,
                                             J_i=as.numeric(factor(pink_dat$stock)),
                                             J_ii=pink_dat$stock_yr,
                                             J=length(unique(pink_dat$stock)),
                                             R_S=pink_dat$logR_S,
                                             S=X_s),
                                 control = list(adapt_delta = 0.99,max_treedepth=20),init=0, warmup = 100, chains = 4, iter = 300)

d=extract(pink_indfit_avg)

plot(apply(d$log_a_g,2,median)~seq(min(pink_stocks$begin),max(pink_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2,xlab='year')
for(i in 1:nrow(pink_stocks)){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(pink_stocks$begin),max(pink_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(pink_stocks$begin),max(pink_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}





#Presy plots####
sd_chin=subset(stock_dat2,stock.id==202)

plot(recruits~spawners,d)

mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

l1=lm(sd_chin$logR_S~sd_chin$spawners)
ndata=data.frame(spawners=seq(0,max(sd_chin$spawners)))
ndata$recruits=exp(l1$coefficients[1]+l1$coefficients[2]*ndata$spawners)*newdata$spawners

pdf('s_r_1.pdf',width=8,height=6)
par(mar=c(4,5,2,2),oma=c(0,0,0,0))
plot(sd_chin$recruits~sd_chin$spawners,xlab='Spawners',ylab='Recruits',cex.axis=1.5,cex.lab=1.5,pch=21,bg=adjustcolor('black',alpha.f = 0.5),cex=1.7)
dev.off()

pdf('s_r_2.pdf',width=8,height=6)
par(mar=c(4,5,2,2),oma=c(0,0,0,0))
plot(sd_chin$recruits~sd_chin$spawners,xlab='Spawners',ylab='Recruits',cex.axis=1.5,cex.lab=1.5,pch=21,bg=adjustcolor('black',alpha.f = 0.5),cex=1.7)
lines(ndata$recruits~ndata$spawners,lwd=4)
dev.off()


pdf('resids.pdf',width=8,height=6)
par(mar=c(4,5,2,2),oma=c(0,0,0,0))
plot(residuals(l1)~df$by,xlab='Year',ylab='Productivity Residuals',cex.axis=1.5,cex.lab=1.5,pch=21,bg=adjustcolor('black',alpha.f = 0.5),cex=1.7)
abline(h=0)
dev.off()

pdf('resids.pdf',width=8,height=6)
par(mar=c(4,5,2,2),oma=c(0,0,0,0))
resids=residuals(l1)
resids=resids[sample(nrow(df))]



df=data.frame(S=sd_chin$spawners,R=sd_chin$recruits,by=sd_chin$broodyear,logRS=sd_chin$logR_S)

TMBtva <- ricker_rw_TMB(data=df,tv.par='a')
TMBtvb <- ricker_rw_TMB(data=df,tv.par='b')
TMBtvab <- ricker_rw_TMB(data=df,tv.par='both')

pdf('sr_tva.pdf',width=14,height=6)
samEst::sr_plot(df,mod=TMBtva,type='rw',par='a',form='tmb',title='')
dev.off()

pdf('sr_tvb.pdf',width=14,height=6)
samEst::sr_plot(df,mod=TMBtvb,type='rw',par='b',form='tmb',title='')
dev.off()

pdf('sr_tvab.pdf',width=14,height=6)
samEst::sr_plot(df,mod=TMBtvb,type='rw',par='both',form='tmb',title='')
dev.off()


p=ggplot(ndata) +
  geom_line(aes(x=spawners,y=recruits),linewidth=2) +
  mytheme + 
  theme(legend.position="right") +
  labs(col = "year") +
  geom_point(data=sd_chin,aes(x=spawners,y=recruits,col=as.factor(broodyear),size=7),alpha=.5)

p=ggplot(sd_chin) +
  mytheme + 
  geom_point(data=sd_chin,aes(x=spawners,y=recruits,size=7),alpha=.5)

#Skeena/Nass stocks####
skeena_stocks=subset(stock_info_filtered,lat==54.2237)
skeena_dat=subset(stock_dat2,stock.id %in% skeena_stocks$stock.id)
length(unique(skeena_dat$stock.id))

stock_year=expand.grid(unique(skeena_dat$stock),unique(skeena_dat$broodyear))
stock_year=stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")

#Find the position of each stock-year combination in all possible stock-year combinations (ie. where in the total vector should we expect this estimate)
skeena_dat$stock_yr=match(paste(skeena_dat$stock,skeena_dat$broodyear,sep='_'),stock_year[,3])

X_s=make_design_matrix(skeena_dat$spawners,grp=skeena_dat$stock)

mc="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int J_i[N];//index of stocks
  int J_ii[N];//index of stock-year combination
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  }
parameters{
  real log_a0; //inital productivity (global)
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<upper = 0>[J] log_b; // rate capacity - fixed in this

 //variance components  
 vector<lower = 0>[J] sigma;

//MVN parameters  
  cholesky_factor_corr[J] Lcorr;
  vector[L-1] z_dev; //average deviation in stock productivity among stocks
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  real<lower = 0> sigma_a; //variance in average productivity among years
  vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years
  
}
transformed parameters{
  vector[J] b; //capacity rate
  //global state process
  vector[L] log_a; //global (ie. average) productivity among stocks over time
  
  //stock state process
  matrix[L,J] log_a_s; //stock deviation from global log_a over time
  matrix[L-1,J] a_dev; //stock-level deviations in year-to-year productivity
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  
  b=exp(log_b);
  
  //initial productivities
  log_a[1] = log_a0; //average productivity among stocks at t = 1
  log_a_s[1,] = to_row_vector(log_a0_s); //stock deviation in productivity at t =1
  log_a_t[1,] = log_a[1]+to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  
  for(t in 1:L-1){
   a_dev[t,] = (diag_pre_multiply(sigma_a_s, Lcorr) * to_vector(z_dev_s[t,]))';
  }
  
  for(t in 2:L){
    log_a[t] = log_a[t-1] + z_dev[t-1]*sigma_a; //global prod. random walk
    log_a_s[t,] = log_a_s[t-1,] + a_dev[t-1,]; //stock-specific deviation in random walk
    log_a_t[t,] = log_a[t] + log_a_s[t,]; //final estimate of stock prod. (global + stock)
  }
  
}  
model{
  //priors
  log_a0 ~ gamma(3,1.5); //initial globalproductivity - wide prior
  log_a0_s ~ normal(0,1); //stock-level deviations around the global productivity
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  z_dev ~ std_normal(); //standardized stock-level deviances in prod
  to_vector(z_dev_s) ~ std_normal(); //global deviation in prod.
  
  Lcorr ~ lkj_corr_cholesky(2.0); // prior for correlation of process deviances
 
  //variance terms
  target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior  
  target += normal_lpdf(sigma_a_s | 0, 0.5) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior
  target += normal_lpdf(sigma_a| 0, 0.5) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior
 

  R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities {
  vector[N] log_lik;

  for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
  corr_matrix[J] Cor_t = multiply_lower_tri_self_transpose(Lcorr);
}


"

mc2="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int J_i[N];//index of stocks
  int J_ii[N];//index of stock-year combination
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  }
parameters{
  real log_a0; //inital productivity (global)
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<upper = 0>[J] log_b; // rate capacity - fixed in this

 //variance components  
 vector<lower = 0>[J] sigma;

//MVN parameters  
  vector[L-1] z_dev; //average deviation in stock productivity among stocks
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  real<lower = 0> sigma_a; //variance in average productivity among years
  vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years
  
}
transformed parameters{
  vector[J] b; //capacity rate
  //global state process
  vector[L] log_a; //global (ie. average) productivity among stocks over time
  //stock state process
  matrix[L,J] log_a_s; //stock deviation from global log_a over time
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  
  b=exp(log_b);
  
  //initial productivities
  log_a[1] = log_a0; //average productivity among stocks at t = 1
  
  log_a_s[1,] = to_row_vector(log_a0_s); //stock deviation in productivity at t =1
  log_a_t[1,] = log_a[1]+to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  for(t in 2:L){
    log_a[t] = log_a[t-1] + z_dev[t-1]*sigma_a; //global prod. random walk
    log_a_s[t,] = log_a_s[t-1,] + z_dev_s[t-1,].*to_row_vector(sigma_a_s); //stock-specific deviation in random walk
    log_a_t[t,] = log_a[t] + log_a_s[t,]; //final estimate of stock prod. (global + stock)
    }
  
}  
model{
  //priors
  log_a0 ~ gamma(3,1); //initial globalproductivity - wide prior
  log_a0_s ~ normal(0,1); //stock-level deviations around the initial global productivity
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  z_dev ~ std_normal(); //standardized stock-level deviances in prod
  to_vector(z_dev_s) ~ std_normal(); //global deviation in prod.
  
  //variance terms
  target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior  
  target += normal_lpdf(sigma_a_s | 0, 0.5) - normal_lcdf(0 | 0, 0.5); //remove density below zero to accomodate half-normal prior
  target += normal_lpdf(sigma_a| 0, 0.5) - normal_lcdf(0 | 0, 0.5); //remove density below zero to accomodate half-normal prior
 
  R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities{
vector[N] log_lik;

for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
}
"

mc3="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int J_i[N];//index of stocks
  int J_ii[N];//index of stock-year combination
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  }
parameters{
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<upper = 0>[J] log_b; // rate capacity - fixed in this

 //variance components  
 vector<lower = 0>[J] sigma;

//MVN parameters  
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years
  
}
transformed parameters{
  vector[J] b; //capacity rate
  //global state process
  //stock state process
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  
  b=exp(log_b);
  
  //initial productivities
  
  log_a_t[1,] = to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  for(t in 2:L){
    log_a_t[t,] = log_a_t[t-1,] + z_dev_s[t-1,].*to_row_vector(sigma_a_s); //stock-specific deviation in random walk
    }
}  
model{
  //priors
  log_a0_s ~ normal(3,2); //stock-level deviations around the initial global productivity
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  to_vector(z_dev_s) ~ std_normal(); //global deviation in prod.
  
  //variance terms
  target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior  
  target += normal_lpdf(sigma_a_s | 0, 0.5) - normal_lcdf(0 | 0, 0.5); //remove density below zero to accomodate half-normal prior

  R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities{
vector[N] log_lik;

for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
}
"

mc_avg="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int J_i[N];//index of stocks
  int J_ii[N];//index of stock-year combination
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  }
parameters{
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<upper = 0>[J] log_b; // rate capacity - fixed in this

 //variance components  
 vector<lower = 0>[J] sigma;

//MVN parameters  
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years
  
}
transformed parameters{
  vector[J] b; //capacity rate
  //global state process
  //stock state process
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  vector[L] log_a_g; //global (mean) productivity
  
  b=exp(log_b);
  
  //initial productivities
  
  log_a_t[1,] = to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  log_a_g[1] = mean(to_array_1d(log_a_t[1,]));
  for(t in 2:L){
    log_a_t[t,] = log_a_t[t-1,] + z_dev_s[t-1,].*to_row_vector(sigma_a_s); //stock-specific deviation in random walk
    log_a_g[t]= mean(to_array_1d(log_a_t[t,]));
  }
    
    
}  
model{
  //priors
  log_a0_s ~ normal(3,2); //stock-level deviations around the initial global productivity
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  to_vector(z_dev_s) ~ std_normal(); //global deviation in prod.
  
  //variance terms
  target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior  
  target += normal_lpdf(sigma_a_s | 0, 0.5) - normal_lcdf(0 | 0, 0.5); //remove density below zero to accomodate half-normal prior

  R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities{
vector[N] log_lik;

for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
}
"

multi_avg="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int J_i[N];//index of stocks
  int J_ii[N];//index of stock-year combination
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  }
parameters{
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<upper = 0>[J] log_b; // rate capacity - fixed in this

 //variance components  
 vector<lower = 0>[J] sigma;

 
//MVN parameters  
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  cholesky_factor_corr[J] Lcorr;
  vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years
  
}
transformed parameters{
  vector[J] b; //capacity rate
  //global state process
  //stock state process
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  matrix[L-1,J] a_dev; //stock-level deviations in year-to-year productivity
  vector[L] log_a_g; //global (mean) productivity
  
  b=exp(log_b);
  
  //initial productivities
  
  log_a_t[1,] = to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  log_a_g[1] = mean(to_array_1d(log_a_t[1,]));
  
    
  for(t in 1:L-1){
   a_dev[t,] = (diag_pre_multiply(sigma_a_s, Lcorr) * to_vector(z_dev_s[t,]))';
  }
  
  
  for(t in 2:L){
    log_a_t[t,] = log_a_t[t-1,] + a_dev[t-1,]; //stock-specific deviation in random walk
    log_a_g[t]= mean(to_array_1d(log_a_t[t,]));
  }
    
    
}  
model{
  //priors
  log_a0_s ~ normal(3,2); //stock-level deviations around the initial global productivity
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  to_vector(z_dev_s) ~ std_normal(); //global deviation in prod.
  
  Lcorr ~ lkj_corr_cholesky(2.0); // prior for correlation of process deviances
 
  //variance terms
  target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior  
  target += normal_lpdf(sigma_a_s | 0, 0.5) - normal_lcdf(0 | 0, 0.5); //remove density below zero to accomodate half-normal prior

  R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities{
vector[N] log_lik;

for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
corr_matrix[J] Cor_t = multiply_lower_tri_self_transpose(Lcorr);
}
"

#Multi-var with global avg
mv_avg_m=rstan::stan_model(model_code=multi_avg)

test_mv = rstan::sampling(mv_avg_m, 
                          data = list(N=nrow(skeena_dat),
                                      L=max(skeena_dat$broodyear)-min(skeena_dat$broodyear)+1,
                                      J_i=as.numeric(factor(skeena_dat$stock)),
                                      J_ii=skeena_dat$stock_yr,
                                      J=length(unique(skeena_dat$stock)),
                                      R_S=skeena_dat$logR_S,
                                      S=X_s),
                          control = list(adapt_delta = 0.99,max_treedepth=20), warmup = 100, chains = 4, iter = 300)
d=rstan::extract(test_mv)

plot(apply(d$log_a_g,2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2,xlab='year')
for(i in 1:13){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}

sum(apply(d$log_lik,2,log_mean_exp))

loo1=loo(test_mv)

#Multivar with evolving global trend
mv_ag_m=rstan::stan_model(model_code=mc)

test_mv2 = rstan::sampling(mv_ag_m, 
                           data = list(N=nrow(skeena_dat),
                                       L=max(skeena_dat$broodyear)-min(skeena_dat$broodyear)+1,
                                       J_i=as.numeric(factor(skeena_dat$stock)),
                                       J_ii=skeena_dat$stock_yr,
                                       J=length(unique(skeena_dat$stock)),
                                       R_S=skeena_dat$logR_S,
                                       S=X_s),
                           control = list(adapt_delta = 0.99,max_treedepth=20), warmup = 100, chains = 4, iter = 300)
d=extract(test_mv2)

plot(apply(d$log_a,2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2,xlab='year')
for(i in 1:13){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}

loo2=loo(test_mv2)
sum(apply(d$log_lik,2,log_mean_exp))

#Independent fit w/ average
test3=rstan::stan_model(model_code=mc_avg)

test_run3 = rstan::sampling(test3, 
                            data = list(N=nrow(skeena_dat),
                                        L=max(skeena_dat$broodyear)-min(skeena_dat$broodyear)+1,
                                        J_i=as.numeric(factor(skeena_dat$stock)),
                                        J_ii=skeena_dat$stock_yr,
                                        J=length(unique(skeena_dat$stock)),
                                        R_S=skeena_dat$logR_S,
                                        S=X_s),
                            control = list(adapt_delta = 0.99,max_treedepth=20), warmup = 100, chains = 4, iter = 300)
d=extract(test_run3)

plot(apply(d$log_a_g,2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:13){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}


loo3=loo(test_run3)

sum(apply(d$log_lik,2,log_mean_exp))

test4=rstan::stan_model(model_code=mc4)

test_run4 = rstan::sampling(test4, 
                            data = list(N=nrow(sk_nass_dat),
                                        L=max(sk_nass_dat$broodyear)-min(sk_nass_dat$broodyear)+1,
                                        J_i=as.numeric(factor(sk_nass_dat$stock)),
                                        J_ii=sk_nass_dat$stock_yr,
                                        J=length(unique(sk_nass_dat$stock)),
                                        R_S=sk_nass_dat$logR_S,
                                        S=X_s),
                            control = list(adapt_delta = 0.99,max_treedepth=20), warmup = 100, chains = 4, iter = 300)
d=extract(test_run4)

plot(apply(d$log_a_t,2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:16){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}

##Compare
loo_mlt=loo(test)


#Fit individually

m3f=samEst::sr_mod(type='rw',par='a',lfo=F)

prod_mat=matrix(ncol=max(sk_nass_stocks$end)-min(sk_nass_stocks$begin)+1,nrow=nrow(sk_nass_stocks))
colnames(prod_mat)=seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end))
rownames(prod_mat)=sk_nass_stocks$stock.name
for(i in 1:16){
  s=subset(sk_nass_dat,stock.id==sk_nass_stocks$stock.id[i])
  
  f3 = rstan::sampling(m3f, 
                       data = list(N=nrow(s),
                                   L=max(s$broodyear)-min(s$broodyear)+1,
                                   ii=s$broodyear-min(s$broodyear)+1,
                                   R_S=s$logR_S,
                                   S=s$spawners),
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 600)
  
  d_3=rstan::extract(f3)
  med_log_a=apply(d_3$log_a,2,median)
  names(med_log_a)=seq(sk_nass_stocks$begin[i],sk_nass_stocks$end[i])
  
  prod_mat[i,c(match(names(med_log_a),colnames(prod_mat)))]=med_log_a
}

prod_avg=apply(prod_mat,2,mean,na.rm=TRUE)

plot(prod_avg~colnames(prod_mat),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:16){
  st_pr=prod_mat[i,];st_pr=st_pr[complete.cases(st_pr)]
  lines(st_pr~names(st_pr),col=adjustcolor('black',alpha.f=0.2))
  points(st_pr~names(st_pr),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}



test2=rstan::stan_model(model_code=mc2)

test_run2 = rstan::sampling(test2, 
                            data = list(N=nrow(sk_nass_dat),
                                        L=max(sk_nass_dat$broodyear)-min(sk_nass_dat$broodyear)+1,
                                        J_i=as.numeric(factor(sk_nass_dat$stock)),
                                        J_ii=sk_nass_dat$stock_yr,
                                        J=length(unique(sk_nass_dat$stock)),
                                        R_S=sk_nass_dat$logR_S,
                                        S=X_s),
                            control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 100, chains = 4, iter = 300)
d=extract(test_run2)

plot(apply(d$log_a,2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:16){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
}

test3=rstan::stan_model(model_code=mc_avg)

test_run3 = rstan::sampling(test3, 
                            data = list(N=nrow(skeena_dat),
                                        L=max(skeena_dat$broodyear)-min(skeena_dat$broodyear)+1,
                                        J_i=as.numeric(factor(skeena_dat$stock)),
                                        J_ii=skeena_dat$stock_yr,
                                        J=length(unique(skeena_dat$stock)),
                                        R_S=skeena_dat$logR_S,
                                        S=X_s),
                            control = list(adapt_delta = 0.99,max_treedepth=20), warmup = 100, chains = 4, iter = 300)
d=extract(test_run3)

plot(apply(d$log_a_g,2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:13){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}

sum(apply(d$log_lik,2,log_mean_exp))

test4=rstan::stan_model(model_code=mc4)

test_run4 = rstan::sampling(test4, 
                            data = list(N=nrow(sk_nass_dat),
                                        L=max(sk_nass_dat$broodyear)-min(sk_nass_dat$broodyear)+1,
                                        J_i=as.numeric(factor(sk_nass_dat$stock)),
                                        J_ii=sk_nass_dat$stock_yr,
                                        J=length(unique(sk_nass_dat$stock)),
                                        R_S=sk_nass_dat$logR_S,
                                        S=X_s),
                            control = list(adapt_delta = 0.99,max_treedepth=20), warmup = 100, chains = 4, iter = 300)
d=extract(test_run4)

plot(apply(d$log_a_t,2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:16){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}


loo_ind=loo(test_run2)
loo_multi=loo(test_run)
loo_ind_2=loo(test_run3)

