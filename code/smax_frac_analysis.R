library(here);library(dplyr);library(ggplot2)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_may2023.csv'))
stock_info<- read.csv(here('data','filtered datasets','all_stocks_info_may2023.csv'))

stock_info_filtered=subset(stock_info,n.years>=16) #242 stocks
stock_info_filtered$stock.name=gsub('/','_',stock_info_filtered$stock.name)
stock_info_filtered$stock.name=gsub('&','and',stock_info_filtered$stock.name)

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock.id)) #242
stock_info_filtered$stock.id2=seq(1:nrow(stock_info_filtered))

if(any(stock_dat2$spawners==0)){stock_dat2$spawners=stock_dat2$spawners+1;stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)}
if(any(stock_dat2$recruits==0)){stock_dat2$recruits=stock_dat2$recruits+1;stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)}
stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)

library(samEst)

samEst::ricker_TMB(data=df)

smax_dat=data.frame(stock=stock_info_filtered$stock.name,maxS=NA,peakS=NA,Smax.s=NA,Smax.ac=NA,maxSfrac=NA,peakSfrac=NA,conv1=NA,conv2=NA)

for(i in 1:nrow(stock_info_filtered)){
  dir.create(here('outputs','figures','sr plots','TMB',paste(i,stock_info_filtered$stock.name[i],sep='-')))
}

for(i in 1:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  s<- s[complete.cases(s$logR_S),] 
  
  df=data.frame(S=s$spawners,R=s$recruits,by=s$broodyear,logRS=s$logR_S)
  
  f=samEst::ricker_TMB(data=df,AC=F)
  fac=samEst::ricker_TMB(data=df,AC=T)
 
   smax_dat$maxS[i]=max(df$S)
  smax_dat$peakS[i]=df$S[which.max(df$R)]
  smax_dat$Smax.s[i]=f$Smax
  smax_dat$Smax.ac[i]=fac$Smax
  smax_dat$maxSfrac[i]=f$Smax/max(df$S)
  smax_dat$peakSfrac[i]=f$Smax/df$S[which.max(df$R)]
  smax_dat$conv1[i]=f$conv_problem
  smax_dat$conv2[i]=fac$conv_problem
  
  x_new=seq(min(df$S),max(df$S),length.out=200)
  pred_df1=data.frame(pred=exp(f$alpha-f$beta*x_new)*x_new,x_new=x_new)
  pred_df2=data.frame(pred=exp(fac$alpha-fac$beta*x_new)*x_new,x_new=x_new)
  
  plot1=ggplot2::ggplot(df, aes(S, R)) +
    geom_line(data=pred_df1,aes(x=x_new,y=pred),linewidth=1.3)+
    geom_point(aes(colour = by),size=2.5) +
    scale_colour_viridis_c(name='Year')+
    ggtitle(paste(stock_info_filtered$stock.name[i],'- simple'))+
    xlab("Spawners") + 
    ylab("Recruits")+
    xlim(0, max(df$S))+
    ylim(0, max(df$R))+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
  
  plot2=ggplot2::ggplot(df, aes(S, R)) +
    geom_line(data=pred_df2,aes(x=x_new,y=pred),linewidth=1.3)+
    geom_point(aes(colour = by),size=2.5) +
    scale_colour_viridis_c(name='Year')+
    ggtitle(paste(stock_info_filtered$stock.name[i],'- autocorr'))+
    xlab("Spawners") + 
    ylab("Recruits")+
    xlim(0, max(df$S))+
    ylim(0, max(df$R))+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
  
  plot_comp=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                               plot2 + theme(legend.position="none"),
                               ncol=2,nrow=1,labels=c("A","B"))
  
  pdf(here('outputs','figures','sr plots','TMB',paste(paste(i,stock_info_filtered$stock.name[i],sep='-'),'.pdf',sep='')),width=12,height=5)
  print(plot_comp)
  dev.off()
  #dev.off(here('outputs','figures','sr plots','TMB',paste(paste(i,stock_info_filtered$stock.name[i],sep='-'),'.pdf',sep='')))
}
hist(smax_dat$maxSfrac,breaks=30)
subset(smax_dat,maxSfrac>10)

rest=subset(smax_dat,maxSfrac<10)

hist(rest$maxSfrac,breaks=30)
hist(rest$peakSfrac,breaks=30)
summary(rest$peakSfrac)
summary(rest$maxSfrac)


smax_dat2=subset(smax_dat,conv1==FALSE&&conv2==FALSE)
smax_dat$diff=(smax_dat$Smax.ac-smax_dat$Smax.s)/smax_dat$Smax.s
summary(smax_dat$diff)
subset(smax_dat,diff< -0.3)
hist(smax_dat$Smax.ac-smax_dat$Smax.s,breaks=30)
