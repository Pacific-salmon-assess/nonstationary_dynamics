library(here);library(samEst)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation2024-10-17.csv'))
stock_info<- read.csv(here('data','filtered datasets','stock_info2024-10-17.csv'))

head(stock_info)

stock_info2=subset(stock_info,n.years>=10)
stock_dat2=subset(stock_dat,stock.id %in% stock_info2$stock.id)

df=matrix(ncol=8,nrow=nrow(stock_info2))
df=as.data.frame(df)
names(df)=c('alpha','smax','smsy','umsy','sig','rho','maxS','scale_smax')
df$stock=stock_info2$stock.name
df$species=stock_info2$species

stock_dat2$species2=gsub("-.*","",stock_dat2$species)

for(j in 1:nrow(stock_info2)){
  s=subset(stock_dat2,stock==stock_info2$stock.name[j])
  d=data.frame(R=s$recruits,S=s$spawners,by=s$broodyear,logRS=s$logRS)
  
  ric_ac=samEst::ricker_TMB(d,AC=T,priors_flag = 0)
  
  df[j,1]=ric_ac$logalpha;df[j,2]=ric_ac$Smax;df[j,3]=ric_ac$Smsy;df[j,4]=ric_ac$umsy;df[j,5]=ric_ac$sigar;df[j,6]=ric_ac$rho;df[j,7]=max(d$S);df[j,8]=df[j,2]/df[j,7]

  pdf(here('outputs','plots',unique(s$species2),paste(unique(s$stock.id),'_',unique(s$stock),'.pdf',sep='')),height=6,width=8)
  plot(d$R~d$S,xlim=c(0,max(d$S)),ylim=c(0,max(d$R)),type='n',bty='l',xlab='spawners',ylab='recruits',main=gsub('_',' ',unique(s$stock)))
  abline(c(0,1),lty=5)
  lines(x=rep(ric_ac$Smax,2),y=c(par('usr')[1],exp(ric_ac$logalpha-ric_ac$beta*ric_ac$Smax)*ric_ac$Smax),col='darkred')
  lines(x=rep(ric_ac$Smsy,2),y=c(par('usr')[1],exp(ric_ac$logalpha-ric_ac$beta*ric_ac$Smsy)*ric_ac$Smsy),col='navy')
  x_n=seq(0,max(d$S))
  p_n=exp(ric_ac$logalpha-ric_ac$beta*x_n)*x_n
  lines(p_n~x_n,lwd=3,col=adjustcolor('black',alpha.f=0.8))
  by_q=round(quantile(d$by,seq(0,1,length.out=8)))
  cols=RColorBrewer::brewer.pal(9, 'YlGnBu')
  cols=cols[-1]
  col.p=cols[cut(d$by, 8, labels = FALSE)]
  lines(d$R~d$S,col=adjustcolor('black',alpha.f=0.2),lwd=0.5)
  points(d$R~d$S,pch=21,bg=col.p,cex=1.5)

  text(x=d$S-max(d$S)*0.01,y=d$R+max(d$R)*0.03,d$by,cex=0.7)
  text(y=par('usr')[4]*0.95,x=par('usr')[2]*0.01,paste('log(a):',round(ric_ac$logalpha,2),sep=' '),adj=0)
  text(y=par('usr')[4]*0.9,x=par('usr')[2]*0.01,paste('Smax:',round(ric_ac$Smax),sep=' '),adj=0,col='darkred')
  text(y=par('usr')[4]*0.85,x=par('usr')[2]*0.01,paste('Smsy:',round(ric_ac$Smsy),sep=' '),adj=0,col='navy')
  text(y=par('usr')[4]*0.80,x=par('usr')[2]*0.01,paste('Umsy:',round(ric_ac$umsy,2),sep=' '),adj=0)
  text(y=par('usr')[4]*0.75,x=par('usr')[2]*0.01,paste('sigma:',round(ric_ac$sigar,2),sep=' '),adj=0)
  text(y=par('usr')[4]*0.70,x=par('usr')[2]*0.01,paste('rho:',round(ric_ac$rho,2),sep=' '),adj=0)
  
  dev.off()

  
  }


hist(df2$scale_smax,breaks=50,xlab='Smax/max(Spawners)',main='')

plot(log(df$smax)~df$alpha)

#species-scaled curves
chin_i=subset(stock_info2,species=='Chinook')
x=c(0,1)
y=c(0,1)

plot(y~x,ylab='scaled recruitment',xlab='scaled spawners',type='n',bty='l')
for(j in 1:nrow(chin_i)){
  s=subset(stock_dat2,stock==chin_i$stock.name[j])
  d=data.frame(R=s$recruits,S=s$spawners,by=s$broodyear,logRS=s$logRS)
  
  ric_ac=samEst::ricker_TMB(d,AC=T,priors_flag = 0)
  
  sn=seq(0,max(d$S))
  sc=sn/max(sn)
  pr=exp(ric_ac$alpha-ric_ac$beta*sn)*sn
  prsc=pr/max(d$R)
  
  lines(prsc~sc)
}

