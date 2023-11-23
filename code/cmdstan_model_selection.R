library(here);library(dplyr);library(ggplot2)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation2023-11-20.csv'))
stock_info<- read.csv(here('data','filtered datasets','stock_info2023-11-20.csv'))

stock_info_filtered=subset(stock_info,n.years>=15) #242 stocks
stock_info_filtered$stock.name=gsub('/','_',stock_info_filtered$stock.name)
stock_info_filtered$stock.name=gsub('&','and',stock_info_filtered$stock.name)

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock.id)) #289
stock_info_filtered$stock.id2=seq(1:nrow(stock_info_filtered))

stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)

library(cmdstanr)
options(mc.cores=6)
#Define models (helps prevent crashing)
file1=file.path(cmdstanr::cmdstan_path(),'sr models', "m1f_ip.stan")
m1=cmdstanr::cmdstan_model(file1)
file2=file.path(cmdstanr::cmdstan_path(),'sr models', "m2f_ip.stan")
m2=cmdstanr::cmdstan_model(file2)
file3=file.path(cmdstanr::cmdstan_path(),'sr models', "m3f_ip.stan")
m3=cmdstanr::cmdstan_model(file3)
file3b=file.path(cmdstanr::cmdstan_path(),'sr models', "m3f_ip_b.stan")
m3b=cmdstanr::cmdstan_model(file3b)
file3s=file.path(cmdstanr::cmdstan_path(),'sr models', "m3f_ip_siga.stan")
m3s=cmdstanr::cmdstan_model(file3s)

file4=file.path(cmdstanr::cmdstan_path(),'sr models', "m4f_smax_ip.stan")
m4=cmdstanr::cmdstan_model(file4)
file5=file.path(cmdstanr::cmdstan_path(),'sr models', "m5f_ip.stan")
m5=cmdstanr::cmdstan_model(file5)
file6=file.path(cmdstanr::cmdstan_path(),'sr models', "m6f_ip.stan")
m6=cmdstanr::cmdstan_model(file6)
file7=file.path(cmdstanr::cmdstan_path(),'sr models', "m7f_ip.stan")
m7=cmdstanr::cmdstan_model(file7)
file8=file.path(cmdstanr::cmdstan_path(),'sr models', "m8f_ip.stan")
m8=cmdstanr::cmdstan_model(file8)

file1=file.path(cmdstanr::cmdstan_path(),'sr models', "m1loo_ip.stan")
m1=cmdstanr::cmdstan_model(file1)
file2=file.path(cmdstanr::cmdstan_path(),'sr models', "m2loo_ip.stan")
m2=cmdstanr::cmdstan_model(file2)
file3=file.path(cmdstanr::cmdstan_path(),'sr models', "m3loo_ip.stan")
m3=cmdstanr::cmdstan_model(file3)
file4=file.path(cmdstanr::cmdstan_path(),'sr models', "m4loo_smax_ip.stan")
m4=cmdstanr::cmdstan_model(file4)
file5=file.path(cmdstanr::cmdstan_path(),'sr models', "m5loo_smax_ip.stan")
m5=cmdstanr::cmdstan_model(file5)
file6=file.path(cmdstanr::cmdstan_path(),'sr models', "m6loo_ip.stan")
m6=cmdstanr::cmdstan_model(file6)
file7=file.path(cmdstanr::cmdstan_path(),'sr models', "m7loo_ip.stan")
m7=cmdstanr::cmdstan_model(file7)
file8=file.path(cmdstanr::cmdstan_path(),'sr models', "m8loo_ip.stan")
m8=cmdstanr::cmdstan_model(file8)


#Create directories for outputs
for(i in 1:nrow(stock_info_filtered)){
  dir.create(here('outputs','sr plots',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep='')))
  dir.create(here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep='')))
  dir.create(here('outputs','loo',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep='')))
  dir.create(here('outputs','model summary',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep='')))
}

loo_comp=list()#list for loo comparisons

f1_l=list();f2_l=list();f3_l=list();f4_l=list();f5_l=list();f6_l=list();f7_l=list();
mfit_sum=data.frame(stock=stock_info_filtered$stock.name,m1.rhat=NA,m1.neff=NA,m2.rhat=NA,m2.neff=NA,m3.rhat=NA,m3.neff=NA,m4.rhat=NA,m4.neff=NA,m5.rhat=NA,m5.neff=NA,m6.rhat=NA,m6.neff=NA,m7.rhat=NA,m7.neff=NA)
for(i in 231:289){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  s<- s[complete.cases(s$logR_S),]
  s<- subset(s,spawners!=0&recruits!=0)
  
  df=data.frame(S=s$spawners,R=s$recruits,by=s$broodyear,logRS=s$logR_S)
  dl=list(S=s$spawners,R=s$recruits,by=s$broodyear,T=nrow(s),N=nrow(s),R_S=s$logR_S,L=max(s$broodyear)-min(s$broodyear)+1,ii=s$broodyear-min(s$broodyear)+1,K=2,alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2),pSmax_mean=0.5*max(df$S),pSmax_sig=0.5*max(df$S),psig_b=max(df$S)*0.5,ll=ifelse(is.na(match(seq(min(df$by),max(df$by)),df$by))==T,0,1))
  

  
  f2 <- m2$sample(data=dl,
                  seed = 123,
                  chains = 6, 
                  iter_warmup = 200,
                  iter_sampling = 1200,
                  refresh=500,
                  adapt_delta = 0.95,
                  max_treedepth = 20)
  
  f3 <- m3$sample(data=dl,
                  seed = 123,
                  chains = 6, 
                  iter_warmup = 200,
                  iter_sampling = 1200,
                  refresh=500,
                  adapt_delta = 0.99,
                  max_treedepth = 20)
  
  
  
  mfit_sum[i,2]=max(f1$summary()[,grepl('rhat',colnames(f1$summary()))])
  mfit_sum[i,3]=min(f1$summary()[,grepl('ess_bulk',colnames(f1$summary()))])
  mfit_sum[i,4]=max(f2$summary()[,grepl('rhat',colnames(f2$summary()))])
  mfit_sum[i,5]=min(f2$summary()[,grepl('ess_bulk',colnames(f2$summary()))])
  mfit_sum[i,6]=max(f3$summary()[,grepl('rhat',colnames(f3$summary()))])
  mfit_sum[i,7]=min(f3$summary()[,grepl('ess_bulk',colnames(f3$summary()))])
  mfit_sum[i,8]=max(f4$summary()[,grepl('rhat',colnames(f4$summary()))])
  mfit_sum[i,9]=min(f4$summary()[,grepl('ess_bulk',colnames(f4$summary()))])
  mfit_sum[i,10]=max(na.omit(f5$summary()[,grepl('rhat',colnames(f5$summary()))]))
  mfit_sum[i,11]=min(na.omit(f5$summary()[,grepl('ess_bulk',colnames(f5$summary()))]))
  mfit_sum[i,12]=max(na.omit(f6$summary()[,grepl('rhat',colnames(f6$summary()))]))
  mfit_sum[i,13]=min(na.omit(f6$summary()[,grepl('ess_bulk',colnames(f6$summary()))]))
  mfit_sum[i,14]=max(na.omit(f7$summary()[,grepl('rhat',colnames(f7$summary()))]))
  mfit_sum[i,15]=min(na.omit(f7$summary()[,grepl('ess_bulk',colnames(f7$summary()))]))
 
  
 
  #store posterior for key variables
  f2_l[[i]]=f2$draws(format='draws_matrix',variables=c('log_a','b','Smax','Umsy','Smsy','sigma','sigma_AR','rho'))
  write.csv(f2_l[[i]],here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m2.csv'))
  f3_l[[i]]=f3$draws(format='draws_matrix',variables=c('log_a','b','Smax','Umsy','Smsy','sigma','sigma_a'))
  write.csv(f3_l[[i]],here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m3.csv'))
  
  #store model summary - neff, diagnostics, etc.
  write.csv(f2_l[[i]],here('outputs','model summary',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m2.csv'))
  write.csv(f3_l[[i]],here('outputs','model summary',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m3.csv'))

}
f3_b=list();
for(i in 1:289){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  s<- s[complete.cases(s$logR_S),]
  s<- subset(s,spawners!=0&recruits!=0)
  
  df=data.frame(S=s$spawners,R=s$recruits,by=s$broodyear,logRS=s$logR_S)
  dl=list(S=s$spawners,R=s$recruits,by=s$broodyear,T=nrow(s),N=nrow(s),R_S=s$logR_S,L=max(s$broodyear)-min(s$broodyear)+1,ii=s$broodyear-min(s$broodyear)+1,K=2,alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2),pSmax_mean=0.5*max(df$S),pSmax_sig=0.5*max(df$S),psig_b=max(df$S)*0.5,ll=ifelse(is.na(match(seq(min(df$by),max(df$by)),df$by))==T,0,1))
  
  f3b <- m3b$sample(data=dl,
                  seed = 123,
                  chains = 6, 
                  iter_warmup = 200,
                  iter_sampling = 1200,
                  refresh=500,
                  adapt_delta = 0.99,
                  max_treedepth = 20)
  
  f3_b[[i]]=f3b$draws(format='draws_matrix',variables=c('log_a','b','Smax','Umsy','Smsy','sigma','sigma_a'))
  write.csv(f3_b[[i]],here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m3b.csv'))
}
  
  
  f3 <- m3$sample(data=dl,
                  seed = 123,
                  chains = 6, 
                  iter_warmup = 200,
                  iter_sampling = 1200,
                  refresh=500,
                  adapt_delta = 0.99,
                  max_treedepth = 20)

df1 <- m1$sample(data=dl,
                seed = 123,
                chains = 6, 
                iter_warmup = 200,
                iter_sampling = 1200,
                refresh=500,
                adapt_delta = 0.95,
                max_treedepth = 20)

f4 <- m4$sample(data=dl,
                seed = 1234,
                chains = 6, 
                iter_warmup = 200,
                iter_sampling = 1200,
                refresh=500,
                adapt_delta = 0.95,
                max_treedepth = 20)

f5 <- m6$sample(data=dl,
                seed = 1234,
                chains = 6, 
                iter_warmup = 200,
                iter_sampling = 1200,
                refresh=500,
                adapt_delta = 0.95,
                max_treedepth = 20)

f6 <- m7$sample(data=dl,
                seed = 123,
                chains = 6, 
                iter_warmup = 200,
                iter_sampling = 1200,
                refresh=500,
                adapt_delta = 0.95,
                max_treedepth = 20)

f7<- m8$sample(data=dl,
               seed = 123,
               chains = 6, 
               iter_warmup = 200,
               iter_sampling = 1200,
               refresh=500,
               adapt_delta = 0.95,
               max_treedepth = 20)
##outputs###



pl1=list() #params list - m1
pl2=list() #params list - m2
pl3=list() #params list - m3
pl4=list() #params list - m4
pl5=list() #params list - m5
pl6=list() #params list - m6
pl7=list() #params list - m7


for(i in 1:nrow(stock_info_filtered)){
  pl1[[i]]=read.csv(here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m1.csv'))  
  pl2[[i]]=read.csv(here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m2.csv'))  
  pl3[[i]]=read.csv(here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m3.csv'))  
  pl4[[i]]=read.csv(here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m4.csv'))  
  pl5[[i]]=read.csv(here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m5.csv'))  
  pl6[[i]]=read.csv(here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m6.csv'))  
  pl7[[i]]=read.csv(here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m7.csv'))  
}

m_sigma=unlist(lapply(pl1,function(x) median(x[,intersect(names(x), 'sigma')])))
m_rho=unlist(lapply(pl2,function(x) median(x[,intersect(names(x), 'rho')])))
m_alpha=unlist(lapply(pl2,function(x) median(x[,intersect(names(x), 'log_a')])))

d_sigma_l=lapply(pl1,function(x) x[,intersect(names(x), 'sigma')])
d_sigma=do.call('cbind',d_sigma_l)
ncol(d_sigma)

par(mar=c(4,5,1,1))
hist(d_sigma[1,],col=adjustcolor('gray',alpha=0.02),xlim=c(min(d_sigma),max(d_sigma)),border=adjustcolor('white',alpha=0.02),breaks=30,freq=T,ylim=c(0,50),main='',xlab=expression(paste(sigma)))
for(t in 2:nrow(d_sigma)){
  par(new=T)
  hist(d_sigma[t,],col=adjustcolor('gray',alpha=0.02),xlim=c(min(d_sigma),max(d_sigma)),border=adjustcolor('white',alpha=0.02),breaks=30,xlab='',ylab='',yaxt='n',main='',xaxt='n',freq=F,ylim=c(0,2.5))
  
}
par(new=T)
hist(m_sigma,col=adjustcolor('black',alpha=0.4),xlim=c(min(d_sigma),max(d_sigma)),border=adjustcolor('white',alpha=0.2),breaks=30,freq=T,ylim=c(0,50),main='',xlab=expression(paste(sigma)))
summary(as.vector(d_sigma))


d_rho_l=lapply(pl2,function(x) x[,intersect(names(x), 'rho')])
d_rho=do.call('cbind',d_rho_l)
ncol(d_rho)


par(mar=c(4,5,1,1))
hist(d_rho[1,],col=adjustcolor('gray',alpha=0.02),xlim=c(min(d_rho),max(d_rho)),border=adjustcolor('white',alpha=0.02),breaks=30,freq=T,ylim=c(0,50),main='',xlab=expression(paste(rho)))
for(t in 2:nrow(d_rho)){
  par(new=T)
  hist(d_rho[t,],col=adjustcolor('gray',alpha=0.02),xlim=c(min(d_rho),max(d_rho)),border=adjustcolor('white',alpha=0.02),breaks=30,xlab='',ylab='',yaxt='n',main='',xaxt='n',freq=F,ylim=c(0,2.5))
  
}
par(new=T)
hist(m_rho,col=adjustcolor('black',alpha=0.4),xlim=c(min(d_rho),max(d_rho)),border=adjustcolor('white',alpha=0.2),breaks=30,freq=T,ylim=c(0,50),main='',xlab=expression(paste(rho)))
summary(as.vector(d_rho))

hist(as.vector(d_rho))

x_new=seq(0,max(df$S),length.out=200)


pdf(here('outputs','figures','est comp',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],'-model_fits.pdf',sep='')),width=8,height=20)

#static
par(mfrow=c(6,2))
by_q=round(quantile(df$by,seq(0,1,by=0.1)))
cols=rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))
col.p=cols[cut(df$by, 11, labels = FALSE)]
plot(R~S,data=df,bty='l',pch=21,cex=1.5,bg='black',main=paste(stock_info_filtered$stock.name[i]),xlim=c(0,max(df$S)),ylim=c(0,max(df$R)),type='n')
lines(R~S,data=df)
text(y=par('usr')[4]*0.9,x=par('usr')[2]*0.05,paste('alpha:',round(median(f1_l[[i]][,grepl('log_a',colnames(f1_l[[i]]))]),3),sep=' '))
text(y=par('usr')[4]*0.85,x=par('usr')[2]*0.05,paste('smax:',round(median(f1_l[[i]][,grepl('Smax',colnames(f1_l[[i]]))]),0),sep=' '))
text(y=par('usr')[4]*0.8,x=par('usr')[2]*0.05,paste('sigma:',round(median(f1_l[[i]][,grepl('sigma',colnames(f1_l[[i]]))]),2),sep=' '))
text(y=par('usr')[4]*0.75,x=par('usr')[2]*0.05,paste('rho:',round(median(f2_l[[i]][,grepl('rho',colnames(f2_l[[i]]))]),2),sep=' '))
mtext(side=3,at=(0.65)*max(df$S),'begins',padj=-1.4,cex=0.8)
mtext(side=3,at=(0.65)*max(df$S),by_q[1],col=cols[1],padj=.2,cex=0.8)
mtext(side=3,at=(0.70)*max(df$S),by_q[3],col=cols[3],padj=.2,cex=0.8)
mtext(side=3,at=(0.75)*max(df$S),by_q[5],col=cols[5],padj=.2,cex=0.8)
mtext(side=3,at=(0.80)*max(df$S),by_q[7],col=cols[7],padj=.2,cex=0.8)
mtext(side=3,at=(0.85)*max(df$S),by_q[9],col=cols[9],padj=.2,cex=0.8)
mtext(side=3,at=(0.9)*max(df$S),by_q[11],col=cols[11],padj=.2,cex=0.8)
mtext(side=3,at=(0.9)*max(df$S),'ends',padj=-1.4,cex=0.8)

p1=exp(median(f1_l[[i]][,grepl('log_a',colnames(f1_l[[i]]))])-median(f1$draws(format='draws_matrix',variables=c('b'))*x_new))*x_new
p2=exp(median(f2_l[[i]][,grepl('log_a',colnames(f2_l[[i]]))])-median(f2$draws(format='draws_matrix',variables=c('b'))*x_new))*x_new

lines(p1~x_new,lwd=3)
lines(p2~x_new,lwd=3,col='darkgray',lty=5)
points(R~S,data=df,pch=21,cex=1.5,bg=col.p)

pr1=median(f1_l[[i]][,grepl('log_a',colnames(f1_l[[i]]))])-median(f1$draws(format='draws_matrix',variables=c('b')))*df$S
plot(df$logRS-pr1~df$by,bty='l',pch=21,bg=col.p,cex=1.5,ylab='residuals',xlab='year')
abline(h=0)

#dynamic prod
#par(mfrow=c(1,2))
by_q=round(quantile(df$by,seq(0,1,by=0.1)))
cols=rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))
col.p=cols[cut(df$by, 11, labels = FALSE)]
plot(R~S,data=df,bty='l',pch=21,cex=1.5,bg='black',main='dynamic prod.',xlim=c(0,max(df$S)),ylim=c(0,max(df$R)),type='n')
lines(R~S,data=df)
for(t in 1:nrow(df)){
  p=exp(median(f3_l[[i]][,grepl('log_a',colnames(f3_l[[i]]))][,t])-median(f3$draws(format='draws_matrix',variables=c('b'))*x_new))*x_new
  lines(p~x_new,lwd=2,col=col.p[t])
}
points(R~S,data=df,pch=21,cex=1.5,bg=col.p)

plot(apply(f3_l[[i]][,grepl('log_a',colnames(f3_l[[i]]))],2,median)~seq(min(df$by),max(df$by)),data=df,bty='l',pch=21,cex=1.5,bg=col.p,ylab='alpha',xlab='year',ylim=c(quantile(as.vector(f3_l[[i]][,grepl('log_a',colnames(f3_l[[i]]))]),0.01),quantile(f3_l[[i]][,grepl('log_a',colnames(f3_l[[i]]))],0.99)))
lci=apply(f3_l[[i]][,grepl('log_a',colnames(f3_l[[i]]))],2,quantile,0.05)
uci=apply(f3_l[[i]][,grepl('log_a',colnames(f3_l[[i]]))],2,quantile,0.95)
x<- c(seq(min(df$by),max(df$by)), rev(seq(min(df$by),max(df$by))))
y<- c(lci, rev(uci))
polygon(x, y, col = adjustcolor('gray', alpha = 0.2), border=NA) # Add uncertainty polygon
lines(apply(f3_l[[i]][,grepl('log_a',colnames(f3_l[[i]]))],2,median)~seq(min(df$by),max(df$by)),lwd=1)


#dynamic smax
#  par(mfrow=c(1,2))
by_q=round(quantile(df$by,seq(0,1,by=0.1)))
cols=rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))
col.p=cols[cut(seq(min(df$by),max(df$by)), 11, labels = FALSE)]
plot(R~S,data=df,bty='l',pch=21,cex=1.5,bg='black',main='dynamic cap.',xlim=c(0,max(df$S)),ylim=c(0,max(df$R)),type='n')
lines(R~S,data=df)
for(t in 1:nrow(df)){
  p=exp(median(f4_l[[i]][,grepl('log_a',colnames(f4_l[[i]]))])-median(f4$draws(format='draws_matrix',variables=c('b'))[,t])*x_new)*x_new
  lines(p~x_new,lwd=2,col=col.p[t])
}
points(R~S,data=df,pch=21,cex=1.5,bg=col.p)


plot(apply(f4_l[[i]][,grepl('Smax',colnames(f4_l[[i]]))],2,median)~seq(min(df$by),max(df$by)),data=df,bty='l',pch=21,cex=1.5,bg=col.p,ylab='Smax',xlab='year',ylim=c(0,quantile(f4_l[[i]][,grepl('Smax',colnames(f4_l[[i]]))],0.95)))
lci=apply(f4_l[[i]][,grepl('Smax',colnames(f4_l[[i]]))],2,quantile,0.05)
uci=apply(f4_l[[i]][,grepl('Smax',colnames(f4_l[[i]]))],2,quantile,0.95)
x<- c(seq(min(df$by),max(df$by)), rev(seq(min(df$by),max(df$by))))
y<- c(lci, rev(uci))
polygon(x, y, col = adjustcolor('gray', alpha = 0.2), border=NA) # Add uncertainty polygon
lines(apply(f4_l[[i]][,grepl('Smax',colnames(f4_l[[i]]))],2,median)~seq(min(df$by),max(df$by)),lwd=1)


#prod regimes
#  par(mfrow=c(1,2))
plot(R~S,data=df,bty='l',pch=21,cex=1.5,type='n',main='regime prod.',xlim=c(0,max(df$S)),ylim=c(0,max(df$R)))
cols=RColorBrewer::brewer.pal(10, 'BrBG')
col.p=cols[cut(1-apply(f5_l[[i]][,grepl('gamma',colnames(f5_l[[i]]))][,1:nrow(df)],2,median), breaks = seq(0,1,by=0.1), labels = FALSE)]
lines(R~S,data=df)
p=exp(median(f5_l[[i]][,grepl('log_a',colnames(f5_l[[i]]))][,1])-median(f5$draws(format='draws_matrix',variables=c('b')))*x_new)*x_new
lines(p~x_new,lwd=2,col=cols[1])
p=exp(median(f5_l[[i]][,grepl('log_a',colnames(f5_l[[i]]))][,2])-median(f5$draws(format='draws_matrix',variables=c('b')))*x_new)*x_new
lines(p~x_new,lwd=2,col=cols[length(cols)])
points(R~S,data=df,pch=21,cex=1.5,bg=col.p)
mtext(side=3,at=(0.65)*max(df$S),'Low prod. Regime',padj=-1.4,cex=0.8)
mtext(side=3,at=(0.65)*max(df$S),by_q[1],col=cols[1],padj=.2,cex=0.8)
mtext(side=3,at=(0.70)*max(df$S),by_q[3],col=cols[3],padj=.2,cex=0.8)
mtext(side=3,at=(0.75)*max(df$S),by_q[5],col=cols[5],padj=.2,cex=0.8)
mtext(side=3,at=(0.75)*max(df$S),'p(R)',padj=-1.4,cex=0.8)
mtext(side=3,at=(0.80)*max(df$S),by_q[7],col=cols[7],padj=.2,cex=0.8)
mtext(side=3,at=(0.85)*max(df$S),by_q[9],col=cols[9],padj=.2,cex=0.8)
mtext(side=3,at=(0.9)*max(df$S),by_q[11],col=cols[11],padj=.2,cex=0.8)
mtext(side=3,at=(0.9)*max(df$S),'High prod. Regime',padj=-1.4,cex=0.8)

plot(1-apply(f5_l[[i]][,grepl('gamma',colnames(f5_l[[i]]))][,1:nrow(df)],2,median)~df$by,data=df,bty='l',pch=21,cex=1.5,bg='black',ylab='p(high prod. regime)',xlab='year',ylim=c(0,1),type='n')
lci=1-apply(f5_l[[i]][,grepl('gamma',colnames(f5_l[[i]]))][,1:nrow(df)],2,quantile,0.05)
uci=1-apply(f5_l[[i]][,grepl('gamma',colnames(f5_l[[i]]))][,1:nrow(df)],2,quantile,0.95)
x<- c(df$by, rev(df$by))
y<- c(lci, rev(uci))
polygon(x, y, col = adjustcolor('gray', alpha = 0.2), border=NA) # Add uncertainty polygon
lines(1-apply(f5_l[[i]][,grepl('gamma',colnames(f5_l[[i]]))][,1:nrow(df)],2,median)~df$by,lwd=1)
points(1-apply(f5_l[[i]][,grepl('gamma',colnames(f5_l[[i]]))][,1:nrow(df)],2,median)~df$by,pch=21,cex=1.5,bg=col.p)

#smax regimes
# par(mfrow=c(1,2))
plot(R~S,data=df,bty='l',pch=21,cex=1.5,type='n',main='regime cap.',xlim=c(0,max(df$S)),ylim=c(0,max(df$R)))
cols=RColorBrewer::brewer.pal(10,'PuOr')
col.p=rev(cols)[cut(1-apply(f6_l[[i]][,grepl('gamma',colnames(f6_l[[i]]))][,1:nrow(df)],2,median), breaks = seq(0,1,by=0.1), labels = FALSE)]
lines(R~S,data=df)
p=exp(median(f6_l[[i]][,grepl('log_a',colnames(f6_l[[i]]))])-median(f6$draws(format='draws_matrix',variables=c('b'))[,1])*x_new)*x_new
lines(p~x_new,lwd=2,col=cols[length(cols)])
p=exp(median(f6_l[[i]][,grepl('log_a',colnames(f6_l[[i]]))])-median(f6$draws(format='draws_matrix',variables=c('b'))[,2])*x_new)*x_new
lines(p~x_new,lwd=2,col=cols[1])
points(R~S,data=df,pch=21,cex=1.5,bg=col.p)

plot(apply(f6_l[[i]][,grepl('gamma',colnames(f6_l[[i]]))][,1:nrow(df)],2,median)~df$by,data=df,bty='l',pch=21,cex=1.5,bg='black',ylab='p(high cap. regime)',xlab='year',ylim=c(0,1),type='n')
lci=apply(f6_l[[i]][,grepl('gamma',colnames(f6_l[[i]]))][,1:nrow(df)],2,quantile,0.05)
uci=apply(f6_l[[i]][,grepl('gamma',colnames(f6_l[[i]]))][,1:nrow(df)],2,quantile,0.95)
x<- c(df$by, rev(df$by))
y<- c(lci, rev(uci))
polygon(x, y, col = adjustcolor('gray', alpha = 0.2), border=NA) # Add uncertainty polygon
lines(apply(f6_l[[i]][,grepl('gamma',colnames(f6_l[[i]]))][,1:nrow(df)],2,median)~df$by,lwd=1)
points(apply(f6_l[[i]][,grepl('gamma',colnames(f6_l[[i]]))][,1:nrow(df)],2,median)~df$by,pch=21,cex=1.5,bg=col.p)

#prod & smax regimes
# par(mfrow=c(1,2))
plot(R~S,data=df,bty='l',pch=21,cex=1.5,type='n',main='regime prod. cap.',xlim=c(0,max(df$S)),ylim=c(0,max(df$R)))
cols=RColorBrewer::brewer.pal(11,'PuOr')
col.p=cols[cut(1-apply(f7_l[[i]][,grepl('gamma',colnames(f7_l[[i]]))][,1:nrow(df)],2,median), breaks = seq(0,1,by=0.1), labels = FALSE)]
lines(R~S,data=df)
p=exp(median(f7_l[[i]][,grepl('log_a',colnames(f7_l[[i]]))][,1])-median(f7$draws(format='draws_matrix',variables=c('b'))[,1])*x_new)*x_new
lines(p~x_new,lwd=2,col=cols[1])
p=exp(median(f7_l[[i]][,grepl('log_a',colnames(f7_l[[i]]))][,2])-median(f7$draws(format='draws_matrix',variables=c('b'))[,2])*x_new)*x_new
lines(p~x_new,lwd=2,col=cols[length(cols)])
points(R~S,data=df,pch=21,cex=1.5,bg=col.p)

plot(1-apply(f7_l[[i]][,grepl('gamma',colnames(f7_l[[i]]))][,1:nrow(df)],2,median)~df$by,data=df,bty='l',pch=21,cex=1.5,bg='black',ylab='p(high prod. regime)',xlab='year',ylim=c(0,1),type='n')
lci=1-apply(f7_l[[i]][,grepl('gamma',colnames(f7_l[[i]]))][,1:nrow(df)],2,quantile,0.05)
uci=1-apply(f7_l[[i]][,grepl('gamma',colnames(f7_l[[i]]))][,1:nrow(df)],2,quantile,0.95)
x<- c(df$by, rev(df$by))
y<- c(lci, rev(uci))
polygon(x, y, col = adjustcolor('gray', alpha = 0.2), border=NA) # Add uncertainty polygon
lines(1-apply(f7_l[[i]][,grepl('gamma',colnames(f7_l[[i]]))][,1:nrow(df)],2,median)~df$by,lwd=1)
points(1-apply(f7_l[[i]][,grepl('gamma',colnames(f7_l[[i]]))][,1:nrow(df)],2,median)~df$by,pch=21,cex=1.5,bg=col.p)

dev.off()

###multi-stock insights####

normalize=function(x){
  x1=(na.omit(x)-min(na.omit(x)))/(max(na.omit(x))-min(na.omit(x)))
  return(as.numeric(x1))
}

normalize2=function(x){
  x1=(na.omit(x))/mean(na.omit(x))
  return(as.numeric(x1))
}


alpha_mat=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_info_filtered$end)-min(stock_info_filtered$begin)+1)
smax_mat=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_info_filtered$end)-min(stock_info_filtered$begin)+1)
hmm_alphas=matrix(nrow=nrow(stock_info_filtered),ncol=2)
alpha_probs_mat=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_info_filtered$end)-min(stock_info_filtered$begin)+1)
hmm_alpha_t=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_info_filtered$end)-min(stock_info_filtered$begin)+1)

hmm_smax=matrix(nrow=nrow(stock_info_filtered),ncol=2)
smax_probs_mat=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_info_filtered$end)-min(stock_info_filtered$begin)+1)
hmm_smax_t=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_info_filtered$end)-min(stock_info_filtered$begin)+1)

colnames(alpha_mat)=seq(min(stock_info_filtered$begin),max(stock_info_filtered$end))
colnames(smax_mat)=seq(min(stock_info_filtered$begin),max(stock_info_filtered$end))
scaled_alpha_mat2=alpha_mat
scaled_alpha_mat=alpha_mat
scaled_smax_mat2=smax_mat
scaled_smax_mat=smax_mat

for(u in 1:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[u])
  s<- s[complete.cases(s$logR_S),]
  
  prod=as.data.frame(d3_a[[u]])[grepl('log_a',colnames(d3_a[[u]]))] #extract productivity series
  smax=as.data.frame(d4_b[[u]])[grepl('S_max',colnames(d4_b[[u]]))] #extract productivity series
  
  hmm_alphas[u,]=apply(as.data.frame(d6_a[[u]])[grepl('log_a',colnames(d6_a[[u]]))],2,median) #alpha distributions for each regime
  alpha_probs=as.data.frame(d6_a[[u]])[grepl('gamma',colnames(d6_a[[u]]))][,1:nrow(s)] #prob. of low regime
  alpha_probs=1-alpha_probs
  hmm_alpha_t[u,]=hmm_alphas[u,1]*(1-alpha_probs)+hmm_alphas[u,2]*(alpha_probs)
  
  hmm_smax[u,]=apply(as.data.frame(d7_b[[u]])[grepl('S_max',colnames(d7_b[[u]]))],2,median) #smax distributions for each regime
  smax_probs=as.data.frame(d7_b[[u]])[grepl('gamma',colnames(d7_b[[u]]))][,1:nrow(s)] #prob. of low regime
  smax_probs=1-smax_probs
  hmm_smax_t[u,]=hmm_alphas[u,1]*(1-alpha_probs)+hmm_alphas[u,2]*(alpha_probs)
  
  t=seq(stock_info_filtered$begin[u],stock_info_filtered$end[u]) #create brood cohorts for start/end of time-series
  t2=s$broodyear #this is for the hmm probabilities, where you do not infill
  alpha_mat[u,match(t,colnames(alpha_mat))]=apply(exp(prod),2,median) #median of productivity (on normal scale - R/S, 0 - inf)
  smax_mat[u,match(t,colnames(smax_mat))]=apply(smax,2,median) #median of capacity
  
  scaled_alpha_mat2[u,match(t,colnames(alpha_mat))]=normalize2(alpha_mat[u,]) #normalize productivity (as % of median)
  scaled_alpha_mat[u,match(t,colnames(alpha_mat))]=normalize(alpha_mat[u,]) #normalize productivity (as % of maximum)
  scaled_smax_mat2[u,match(t,colnames(smax_mat))]=normalize2(smax_mat[u,]) #normalize capacity (as % of median)
  scaled_smax_mat[u,match(t,colnames(smax_mat))]=normalize(smax_mat[u,]) #normalize capacity (as % of maximum)
  alpha_probs_mat[u,match(t2,colnames(smax_mat))]=apply(alpha_probs,2,median)
  smax_probs_mat[u,match(t2,colnames(smax_mat))]=apply(smax_probs,2,median)
}

plot(c(0,1)~c(min(stock_info_filtered$begin),max(stock_info_filtered$end)),type='n',bty='l',ylab='Scaled Capacity')
for(l in 1:nrow(smax_mat)){
  smax=scaled_smax_mat[l,];smax=smax[is.na(smax)==F]
  lines(smax~as.numeric(names(smax)),col=adjustcolor('darkgray',alpha.f=0.2),lwd=2)
}


moving_average_mat <- function(x, lag = 2) {             # Create user-defined function
  m=NA
  se=NA
  for(t in c(lag+1):c(ncol(x)-lag)){
   x1=x[,c(t-lag):c(t+lag)] 
   m[t]=mean(as.vector(na.omit(x1)))
   se[t]=sd(as.vector(na.omit(x1)))/sqrt(length(as.vector(na.omit(x1))))
  }
  return(rbind(m,se))
}

#Productivity - normalized as % of maximum
plot(c(0,1)~c(min(stock_info_filtered$begin),max(stock_info_filtered$end)),type='n',bty='l',ylab='Scaled Productivity',xlab='Brood Cohort Year')
for(l in 1:nrow(alpha_mat)){
  alpha=scaled_alpha_mat[l,];alpha=alpha[is.na(alpha)==F]
  lines(alpha~as.numeric(names(alpha)),col=adjustcolor('darkgray',alpha.f=0.2),lwd=2)
}

glob_a=moving_average_mat(scaled_alpha_mat)

lines(glob_a[1,]~as.numeric(colnames(scaled_alpha_mat[,1:c(ncol(scaled_alpha_mat)-2)])),lwd=2)
x<- c(as.numeric(colnames(scaled_alpha_mat[,1:c(ncol(scaled_alpha_mat)-2)])), rev(as.numeric(colnames(scaled_alpha_mat[,1:c(ncol(scaled_alpha_mat)-2)]))))
y<- c(glob_a[1,]-glob_a[2,]*2, rev(glob_a[1,]+glob_a[2,]*2))
polygon(x, y, col = adjustcolor('darkgray', alpha = 0.5), border=NA) # Add uncertainty polygon

#Productivity - untransformed
plot(c(0,30)~c(min(stock_info_filtered$begin),max(stock_info_filtered$end)),type='n',bty='l',ylab='Productivity - R/S',xlab='Brood Cohort Year')
for(l in 1:nrow(alpha_mat)){
  alpha=alpha_mat[l,];alpha=alpha[is.na(alpha)==F]
  lines(alpha~as.numeric(names(alpha)),col=adjustcolor('darkgray',alpha.f=0.2),lwd=2)
}

glob_a=moving_average_mat(alpha_mat)

lines(glob_a[1,]~as.numeric(colnames(scaled_alpha_mat[,1:c(ncol(scaled_alpha_mat)-2)])),lwd=2)
x<- c(as.numeric(colnames(scaled_alpha_mat[,1:c(ncol(scaled_alpha_mat)-2)])), rev(as.numeric(colnames(scaled_alpha_mat[,1:c(ncol(scaled_alpha_mat)-2)]))))
y<- c(glob_a[1,]-glob_a[2,]*2, rev(glob_a[1,]+glob_a[2,]*2))
polygon(x, y, col = adjustcolor('darkgray', alpha = 0.5), border=NA) # Add uncertainty polygon


#Capacity - normalized as % of maximum
plot(c(0,1)~c(min(stock_info_filtered$begin),max(stock_info_filtered$end)),type='n',bty='l',ylab='Scaled Capacity',xlab='Brood Cohort Year')
for(l in 1:nrow(smax_mat)){
  smax=scaled_smax_mat[l,];smax=smax[is.na(smax)==F]
  lines(smax~as.numeric(names(smax)),col=adjustcolor('darkgray',alpha.f=0.2),lwd=2)
}

glob_a=moving_average_mat(scaled_smax_mat)

lines(glob_a[1,]~as.numeric(colnames(scaled_smax_mat[,1:c(ncol(scaled_smax_mat)-2)])),lwd=2)
x<- c(as.numeric(colnames(scaled_smax_mat[,1:c(ncol(scaled_smax_mat)-2)])), rev(as.numeric(colnames(scaled_smax_mat[,1:c(ncol(scaled_smax_mat)-2)]))))
y<- c(glob_a[1,]-glob_a[2,]*2, rev(glob_a[1,]+glob_a[2,]*2))
polygon(x, y, col = adjustcolor('darkgray', alpha.f = 0.5), border=NA) # Add uncertainty polygon

#Capacity - untransformed
plot(c(0,3000000)~c(min(stock_info_filtered$begin),max(stock_info_filtered$end)),type='n',bty='l',ylab='Capacity',xlab='Brood Cohort Year')
for(l in 1:nrow(smax_mat)){
  smax=smax_mat[l,];smax=smax[is.na(smax)==F]
  lines(smax~as.numeric(names(smax)),col=adjustcolor('darkgray',alpha.f=0.2),lwd=2)
}

glob_a=moving_average_mat(smax_mat)

lines(glob_a[1,]~as.numeric(colnames(scaled_smax_mat[,1:c(ncol(scaled_smax_mat)-2)])),lwd=2)
x<- c(as.numeric(colnames(scaled_smax_mat[,1:c(ncol(scaled_smax_mat)-2)])), rev(as.numeric(colnames(scaled_smax_mat[,1:c(ncol(scaled_smax_mat)-2)]))))
y<- c(glob_a[1,]-glob_a[2,]*2, rev(glob_a[1,]+glob_a[2,]*2))
polygon(x, y, col = adjustcolor('darkgray', alpha.f = 0.5), border=NA) # Add uncertainty polygon



#Prob. of high prod. regime
plot(c(0,1)~c(min(stock_info_filtered$begin),max(stock_info_filtered$end)),type='n',bty='l',ylab='p(High Prod. Regime)',xlab='Brood Cohort Year')
for(l in 1:nrow(alpha_probs_mat)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[l])
  s<- s[complete.cases(s$logR_S),]
  alpha=alpha_probs_mat[l,];alpha=alpha[is.na(alpha)==F]
  lines(alpha~s$broodyear,col=adjustcolor('darkgray',alpha.f=0.2),lwd=2)
}

glob_a=moving_average_mat(alpha_probs_mat)

lines(glob_a[1,]~as.numeric(colnames(scaled_alpha_mat[,1:c(ncol(scaled_alpha_mat)-2)])),lwd=2)
x<- c(as.numeric(colnames(scaled_alpha_mat[,1:c(ncol(scaled_alpha_mat)-2)])), rev(as.numeric(colnames(scaled_alpha_mat[,1:c(ncol(scaled_alpha_mat)-2)]))))
y<- c(glob_a[1,]-glob_a[2,]*2, rev(glob_a[1,]+glob_a[2,]*2))
polygon(x, y, col = adjustcolor('darkgray', alpha = 0.5), border=NA) # Add uncertainty polygon

#alpha estimate: prob regime x regime
plot(c(-1,4)~c(min(stock_info_filtered$begin),max(stock_info_filtered$end)),type='n',bty='l',ylab='p(High Prod. Regime)',xlab='Brood Cohort Year')

for(l in 1:nrow(alpha_probs_mat)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[l])
  s<- s[complete.cases(s$logR_S),]

  alpha=alpha_probs_mat[l,];alpha=alpha[is.na(alpha)==F]
  alpha_t=alpha*hmm_alphas[l,2]+(1-alpha)*hmm_alphas[l,1]
  
  
  lines(alpha_t~s$broodyear,col=adjustcolor('darkgray',alpha.f=0.2),lwd=2)
}

glob_a=moving_average_mat(alpha_t)

lines(glob_a[1,]~as.numeric(colnames(scaled_alpha_mat[,1:c(ncol(scaled_alpha_mat)-2)])),lwd=2)
x<- c(as.numeric(colnames(scaled_alpha_mat[,1:c(ncol(scaled_alpha_mat)-2)])), rev(as.numeric(colnames(scaled_alpha_mat[,1:c(ncol(scaled_alpha_mat)-2)]))))
y<- c(glob_a[1,]-glob_a[2,]*2, rev(glob_a[1,]+glob_a[2,]*2))
polygon(x, y, col = adjustcolor('darkgray', alpha = 0.5), border=NA) # Add uncertainty polygon


