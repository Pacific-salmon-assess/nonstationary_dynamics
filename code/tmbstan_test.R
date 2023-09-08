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


file3=file.path(cmdstanr::cmdstan_path(),'sr models', "m3f.stan")
m3=cmdstanr::cmdstan_model(file3)


file3cp=file.path(cmdstanr::cmdstan_path(),'sr models', "m3f_cp.stan")
m3cp=cmdstanr::cmdstan_model(file3cp)

file4sm=file.path(cmdstanr::cmdstan_path(),'sr models', "m4f_smax2.stan")
m4sm=cmdstanr::cmdstan_model(file4sm)

file4=file.path(cmdstanr::cmdstan_path(),'sr models', "m4f.stan")
m4=cmdstanr::cmdstan_model(file4)

file4ip=file.path(cmdstanr::cmdstan_path(),'sr models', "m4f_ip.stan")
m4ip=cmdstanr::cmdstan_model(file4ip)

mode=function(x){
  d=density(x,bw=0.02)
  return(d$x[which.max(d$y)])
}

pdf(here('outputs','figures','est comp','test.pdf'))
for(i in 1:30){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  s<- s[complete.cases(s$logR_S),] 
  
  df=data.frame(S=s$spawners,R=s$recruits,by=s$broodyear,logRS=s$logR_S)
  dl=list(S=s$spawners,R=s$recruits,by=s$broodyear,T=nrow(s),N=nrow(s),R_S=s$logR_S,L=max(s$broodyear)-min(s$broodyear)+1,ii=s$broodyear-min(s$broodyear)+1,K=2,alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2),pSmax_mean=0.5*max(df$S),pSmax_sig=0.25*max(df$S))
  
  
  f3 <- m3$sample(data=dl,
                  seed = 123,
                  chains = 6, 
                  iter_warmup = 200,
                  iter_sampling = 600,
                  refresh = 0,
                  adapt_delta = 0.99,
                  max_treedepth = 15)
  
 
  testkf=samEst::ricker_kf_TMB(data=df)
  testrw=samEst::ricker_rw_TMB(data=df,tv.par="a",priors_flag=0)
  testrw_st=samEst::ricker_rw_TMBstan(data=df,tv.par="a",priors_flag=0)
  
  f3_a=apply(f3$draws(variable='log_a',format='draws_matrix'),2,mode)
  f3_b=median(f3$draws(variable='b',format='draws_matrix'))

  tmbst_alpha=apply(as.matrix(testrw_st$alpha),1,mode)
  
  
  plot(testrw$alpha,type='l',col='navy',lwd=2,ylim=c(min(c(f3_a,f3_acp,testkf$alpha,testrw$alpha,tmbst_alpha)),max(c(f3_a,f3_acp,testkf$alpha,testrw$alpha,tmbst_alpha))*1.25),main=stock_info_filtered$stock.name[i])
  lines(apply(f3$draws(variable='log_a',format='draws_matrix'),2,mode),lwd=2,col='darkred')
  lines(testkf$alpha,lwd=2,col='cyan4')
  lines(tmbst_alpha,lwd=2,col='darkorange')
  text(nrow(s)-2,max(testrw$alpha)*1.25,'TMB rw',col='navy')
  text(nrow(s)-2,max(testrw$alpha)*1.2,'TMB kf',col='cyan4')
  text(nrow(s)-2,max(testrw$alpha)*1.15,'Stan rw',col='darkred')
  text(nrow(s)-2,max(testrw$alpha)*1.05,'TMB stan',col='darkorange')
  
}
dev.off()

#smax prior testing
pdf(here('outputs','figures','est comp','test_newsmax4.pdf'))

for(i in 13:20){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  s<- s[complete.cases(s$logR_S),] 
  
  df=data.frame(S=s$spawners,R=s$recruits,by=s$broodyear,logRS=s$logR_S)
  dl=list(S=s$spawners,R=s$recruits,by=s$broodyear,T=nrow(s),N=nrow(s),R_S=s$logR_S,L=max(s$broodyear)-min(s$broodyear)+1,ii=s$broodyear-min(s$broodyear)+1,K=2,alpha_dirichlet=matrix(c(2,1,1,2),ncol=2,nrow=2),pSmax_mean=0.5*max(df$S),pSmax_sig=0.25*max(df$S))
  
  f4sm <- m4sm$sample(data=dl,
                      seed = 123,
                      chains = 6, 
                      iter_warmup = 200,
                      iter_sampling = 600,
                      refresh = 0,
                      adapt_delta = 0.99,
                      max_treedepth = 15)
  
  f4 <- m4ip$sample(data=dl,
                      seed = 1234,
                      chains = 6, 
                      iter_warmup = 200,
                      iter_sampling = 600,
                      refresh = 0,
                      adapt_delta = 0.99,
                      max_treedepth = 15)
  
  sm2=apply(f4sm$draws(variable='Smax',format='draws_matrix'),2,median)
  sm=apply(f4$draws(variable='S_max',format='draws_matrix'),2,median)
  
  

  x_new=seq(0,max(df$S),length.out=200)
  by_q=round(quantile(df$by,seq(0,1,by=0.1)))
  f4_a=median(f4$draws(variable='log_a',format='draws_matrix'))
  f4_b=apply(f4$draws(variable='b',format='draws_matrix'),2,median)
  f4s_a=median(f4sm$draws(variable='log_a',format='draws_matrix'))
  f4s_b=apply(f4sm$draws(variable='b',format='draws_matrix'),2,median)
  
  
  pred_df_f4=data.frame(x_new)
  pred_df_f4s=data.frame(x_new)
  
  for(n in 1:length(by_q)){
    pred_df_f4[,1+n]=exp(f4_a-f4_b[n]*x_new)*x_new
    pred_df_f4s[,1+n]=exp(f4s_a-f4s_b[n]*x_new)*x_new
  }
  
  plot1=ggplot2::ggplot(df, aes(S, R)) +
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,2],colour = by_q[1]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,3],colour = by_q[2]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,4],colour = by_q[3]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,5],colour = by_q[4]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,6],colour = by_q[5]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,7],colour = by_q[6]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,2],colour = by_q[1]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,3],colour = by_q[2]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,4],colour = by_q[3]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,5],colour = by_q[4]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,6],colour = by_q[5]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,7],colour = by_q[6]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,8],colour = by_q[7]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,9],colour = by_q[8]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,10],colour = by_q[9]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,11],colour = by_q[10]),linewidth=1.3)+
    geom_line(data=pred_df_f4,aes(x=x_new,y=pred_df_f4[,12],colour = by_q[11]),linewidth=1.3)+
    geom_point(aes(colour = by),size=2.5) +
    scale_colour_viridis_c(name='Year')+
    xlab("Spawners") + 
    ylab("Recruits")+
    ggtitle(paste(stock_info_filtered$stock.name[i],'logb-wide',sep='-'))+
    xlim(0, max(df$S))+
    ylim(0, max(df$R))+
    theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
  
  plot2=ggplot2::ggplot(df, aes(S, R)) +
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,2],colour = by_q[1]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,3],colour = by_q[2]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,4],colour = by_q[3]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,5],colour = by_q[4]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,6],colour = by_q[5]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,7],colour = by_q[6]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,2],colour = by_q[1]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,3],colour = by_q[2]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,4],colour = by_q[3]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,5],colour = by_q[4]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,6],colour = by_q[5]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,7],colour = by_q[6]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,8],colour = by_q[7]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,9],colour = by_q[8]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,10],colour = by_q[9]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,11],colour = by_q[10]),linewidth=1.3)+
    geom_line(data=pred_df_f4s,aes(x=x_new,y=pred_df_f4s[,12],colour = by_q[11]),linewidth=1.3)+
    geom_point(aes(colour = by),size=2.5) +
    scale_colour_viridis_c(name='Year')+
    ggtitle('smax_ip')+
    xlab("Spawners") + 
    ylab("Recruits")+
    xlim(0, max(df$S))+
    ylim(0, max(df$R))+
    theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
  
  
  gamma_df=data.frame(gamma=sm,by=s$broodyear)
  plot3=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
    geom_line(aes(x=by,y=gamma),linewidth=1.3)+
    geom_point(aes(colour = by),size=4) +
    scale_colour_viridis_c(name='Smax')+
    theme_classic(14)+
    ylab("Smax")+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
  
  gamma_df2=data.frame(gamma=sm2,by=s$broodyear)
  plot4=ggplot2::ggplot(gamma_df2, aes(by,gamma)) +
    geom_line(aes(x=by,y=gamma),linewidth=1.3)+
    geom_point(aes(colour = by),size=4) +
    ylab("Smax")+
    scale_colour_viridis_c(name='Smax')+
    theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
  
  
  legend = cowplot::get_legend(plot1)
  
  plot_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                               plot2 + theme(legend.position="none"),
                            plot3 + theme(legend.position="none"),
                            plot4 + theme(legend.position="none"),
                               ncol=2,nrow=2)
 plot=cowplot::plot_grid(plot_a,legend,rel_widths = c(3,.3))
 
 print(plot)
 
}
dev.off()



plot(seq(0.5,nrow(s)+0.5,length.out=nrow(s))~seq(min(c()),max=min(c()),length.out=nrow(s)),bty='l',xlab='',ylab='',type='n',yaxt='n',ylim=c(1,12.8))
lines(c(0.5,12.5)~rep(0,2),lty=5)
text(y=12.5,x=-1.2,xpd=T,labels="n",adj=0)
for(i in 1:11){
  f_dist=f_r_dist[[f_dat2$x[i]]]
  f_dist$y=f_dist$y+i
  polygon(f_dist,col=adjustcolor('gray',alpha.f=0.7),border=NA)
  lines(rep(i+0.4,2)~c(f_dat2[i,3],f_dat2[i,4]),lwd=1)
  lines(rep(i+0.4,2)~c(f_dat2[i,5],f_dat2[i,6]),lwd=3)
  lines(rep(i+0.4,2)~c(f_dat2[i,7],f_dat2[i,8]),lwd=5)
  #  peak_r=f_dist$x[f_dist$y==max(f_dist$y)]
  #  lines(c(min(f_dist$y),max(f_dist$y))~rep(peak_r,2),lwd=1,col='navy')
  points(i+0.4~f_dat2[i,2],pch=21,bg=adjustcolor('black',alpha.f=0.6),col='white',cex=1.5)
  text(y=i+0.8,x=-1.8,xpd=T,labels=f_dat2$grp[i],adj=0)
  text(y=i+0.8,x=-1.2,xpd=T,labels=f_dat2$n[i],adj=0)
}

  