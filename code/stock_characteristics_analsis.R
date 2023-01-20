#stock characteristics
library(here);library(dplyr);library(rstan);library(ggplot2)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_aug2022.csv'))
stock_info<- read.csv(here('data','filtered datasets','all_stocks_info_aug2022.csv'))

#source(here('code','samEst code','stan_functions.R'))
#source(here('code','samEst code','lfo_stan_functions.R'))
#source(here('code','samEst code','lfo_TMB_functions.R'))
##source(here('code','samEst code','TMB_functions.R'))
#source(here('code','samEst code','util_functions.R'))
library(samEst)
options(mc.cores = parallel::detectCores())
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')

###Load in data####
#Remove stocks with less than 15 years of recruitment data
stock_info_filtered=subset(stock_info,n.years>=18) #242 stocks
stock_info_filtered$stock.name=gsub('/','_',stock_info_filtered$stock.name)
stock_info_filtered$stock.name=gsub('&','and',stock_info_filtered$stock.name)

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock.id)) #242
stock_info_filtered$stock.id2=seq(1:nrow(stock_info_filtered))

stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)


log_a=NA
Smax=NA
sd_R=NA
AR1=NA
log_a_AR1=NA
rw_log_a_0=NA
rw_log_a_f=NA
rw_log_a_trend=NA
rw_log_a_med=NA
rw_log_a_min=NA
rw_log_a_max=NA
rw_Smax_med=NA
rw_Smax_min=NA
rw_Smax_max=NA
rw_Smax_0=NA
rw_Smax_f=NA
rw_Smax_trend=NA
hmm_log_a=matrix(ncol=2,nrow=nrow(stock_info_filtered))
colnames(hmm_log_a)=c('hmm_log_a1','hmm_log_a2')
hmm_Smax=matrix(ncol=2,nrow=nrow(stock_info_filtered))
colnames(hmm_Smax)=c('hmm_Smax1','hmm_Smax2')

for(i in 1:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  s<- s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}

  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  TMBstatic <- ricker_TMB(data=df)
  log_a[i]=TMBstatic$alpha
  Smax[i]=TMBstatic$Smax
  sd_R[i]=TMBstatic$sig
  TMBac <- ricker_TMB(data=df, AC=TRUE)
  AR1[i]=TMBac$rho
  log_a_AR1[i]=TMBac$alpha
  TMBtva <- tryCatch(ricker_rw_TMB(data=df,tv.par='a'),error = function(e) {TMBtva=list(conv_problem=TRUE)})
  if(is.null(TMBtva$alpha)==FALSE){
    rw_log_a_0[i]=TMBtva$alpha[1]
    rw_log_a_f[i]=TMBtva$alpha[nrow(df)]
    rw_log_a_trend[i]=((TMBtva$alpha[nrow(df)]-TMBtva$alpha[1])/TMBtva$alpha[1])/nrow(df)
    rw_log_a_med[i]=median(TMBtva$alpha)
    rw_log_a_min[i]=min(TMBtva$alpha)
    rw_log_a_max[i]=max(TMBtva$alpha)
  }
  
  TMBtvb <- tryCatch(ricker_rw_TMB(data=df, tv.par='b'),error = function(e){TMBtvb=list(conv_problem=TRUE)})
  if(is.null(TMBtvb$Smax)==FALSE){
    
  rw_Smax_0[i]=TMBtvb$Smax[1]
  rw_Smax_f[i]=TMBtvb$Smax[nrow(df)]
  rw_Smax_trend[i]=((TMBtvb$Smax[nrow(df)]-TMBtvb$Smax[1])/TMBtvb$Smax[1])/nrow(df)
  rw_Smax_med[i]=median(TMBtvb$Smax)
  rw_Smax_min[i]=min(TMBtvb$Smax)
  rw_Smax_max[i]=max(TMBtvb$Smax)
  }
  
  TMBhmma <- tryCatch(ricker_hmm_TMB(data=df, tv.par='a'),error = function(e){TMBhmma=list(conv_problem=TRUE)})
  if(TMBhmma$conv_problem==TRUE){hmm_log_a[i,]=rep(NA,2)}else{
    hmm_log_a[i,]=TMBhmma$alpha
  }
  TMBhmmb <- tryCatch(ricker_hmm_TMB(data=df, tv.par='b'),error = function(e){TMBhmmb=list(conv_problem=TRUE)})
  if(TMBhmmb$conv_problem==TRUE){hmm_Smax[i,]=rep(NA,2)}else{
    hmm_Smax[i,]=TMBhmmb$Smax
    }
}

stock_char=cbind(stock_info_filtered,log_a,Smax,sd_R,AR1,log_a_AR1,rw_log_a_0,rw_log_a_f,rw_log_a_trend,rw_log_a_med,rw_log_a_min,rw_log_a_max,rw_Smax_0,rw_Smax_f,rw_Smax_trend,rw_Smax_med,rw_Smax_min,rw_Smax_max,hmm_log_a,hmm_Smax)

#histograms
par(mfrow=c(1,1))
hist(stock_char$log_a,breaks=30)

sp_cols=cbind(c("darkgray", "darkgreen", "darkblue","sienna","darkred"),c('Chinook','Chum','Coho','Pink','Sockeye'))

mod_col=cbind(c("azure4", "azure3", "cyan4","cadetblue3","deepskyblue4",'seagreen','olivedrab','darkgreen'),seq=1:8)

mc_col=cbind(c("azure3","deepskyblue4",'darkgreen'),c('static','dynamic','regime'))

p= subset(stock_char,species=='Pink')
chi= subset(stock_char,species=='Chinook')
co= subset(stock_char,species=='Coho')
so= subset(stock_char,species=='Sockeye')
chu= subset(stock_char,species=='Chum')

pl_chi=ggplot(chi, aes(x=log_a))+
  geom_histogram(color="white",fill='darkgray',position="identity",alpha=0.5) + theme_minimal() +
  xlab('log(alpha)')
theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
      strip.text = element_text(face="bold", size=14),
      axis.text=element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=15))

pl_27=ggplot(sp_27, aes(fill=species, y=n, x=factor(best_mod)))+
  scale_fill_manual(values=sp_cols[match(levels(factor(sp_27$species)),sp_cols[,2])])+
  ggtitle("27 - 35 years")+
  geom_bar(position="stack", stat="identity") + 
  theme_minimal() +
  xlim(factor(seq(1:8)))+
  xlab('')+
  ylab('No. populations')+ 
  theme(text=element_text(size=12),axis.title=element_text(size=12,face="bold"),axis.line = element_line(colour = "black", 
                                                                                                         size = 1, linetype = "solid"),plot.title = element_text(size=14))

pl_35=ggplot(sp_35, aes(fill=species, y=n, x=factor(best_mod)))+
  scale_fill_manual(values=sp_cols[match(levels(factor(sp_35$species)),sp_cols[,2])])+
  ggtitle("35 - 45 years")+
  geom_bar(position="stack", stat="identity") + 
  theme_minimal() +
  xlim(factor(seq(1:8)))+
  xlab('')+
  ylab('')+ 
  theme(text=element_text(size=12), axis.title=element_text(size=12,face="bold"),axis.line = element_line(colour = "black", 
                                                                                                          size = 1, linetype = "solid"),plot.title = element_text(size=14))

pl_45=ggplot(sp_45, aes(fill=species, y=n, x=factor(best_mod)))+
  scale_fill_manual(values=sp_cols[match(levels(factor(sp_45$species)),sp_cols[,2])])+
  ggtitle("46 - 50 years")+
  geom_bar(position="stack", stat="identity") + 
  theme_minimal() +
  xlim(factor(seq(1:8)))+
  xlab('')+
  ylab('No. populations')+ 
  theme(text=element_text(size=12), axis.title=element_text(size=12,face="bold"),axis.line = element_line(colour = "black", 
                                                                                                          size = 1, linetype = "solid"),plot.title = element_text(size=14))

pl_50=ggplot(sp_50, aes(fill=species, y=n, x=factor(best_mod)))+
  scale_fill_manual(values=sp_cols[match(levels(factor(sp_50$species)),sp_cols[,2])])+
  ggtitle("50+ years")+
  geom_bar(position="stack", stat="identity") + 
  theme_minimal() +
  xlim(factor(seq(1:8)))+
  xlab('')+
  ylab('')+ 
  theme(text=element_text(size=12), axis.title=element_text(size=12,face="bold"),axis.line = element_line(colour = "black", 
                                                                                                          size = 1, linetype = "solid"),plot.title = element_text(size=14))

legend= get_legend(pl_18)

plot_grid(pl_18+ theme(legend.position="none"),
          pl_27+ theme(legend.position="none"),
          pl_35+ theme(legend.position="none"),
          pl_45+ theme(legend.position="none"),
          pl_50+ theme(legend.position="none"),
          legend,
          ncol=3,nrow=2,labels=c("A","B","C","D","E"))

ggplot(stock_char, aes(fill=species, x=log_a))+
  scale_fill_manual(values=sp_cols[match(levels(factor(stock_char$species)),sp_cols[,2])])+ 
  geom_histogram(color="white", position="identity",alpha=0.5) + theme_minimal() +
  xlab('log(alpha)')
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=14),
        axis.text=element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=15))

hist(log10(stock_char$Smax),breaks=30)
hist(stock_char$sd_R,breaks=30)
hist(stock_char$AR1,breaks=30)
t=subset(stock_char,AR1< -0.5)

par(mfrow=c(2,3))

hist(chi$log_a,breaks=30)
hist(chu$log_a,breaks=30)
hist(co$log_a,breaks=30)
hist(p$log_a,breaks=30)
hist(so$log_a,breaks=30)

#histograms - rw 
hist(stock_char$rw_log_a_trend*100,breaks=30,xlim = c(-0.15*100,0.15*100))
summary(stock_char$rw_log_a_trend*100)
hist(stock_char$rw_Smax_trend*100,breaks=30)
summary(stock_char$rw_Smax_trend*100)
hist(exp(stock_char$rw_log_a_f)/exp(stock_char$rw_log_a_0),breaks=30,xlim = c(0,10))
abline(v=1,lty=5)
hist(stock_char$rw_Smax_f/stock_char$rw_Smax_0,breaks=30,xlim=c(0,10))
abline(v=1,lty=5)
hist(stock_char$rw_Smax_f/stock_char$rw_Smax_0,breaks=30,xlim=c(0,10))

#
par(mfrow=c(1,1))
hist(stock_char$rw_log_a_max,breaks=30)
hist(stock_char$rw_log_a_min,breaks=30)
hist(stock_char$rw_log_a_max/stock_char$rw_log_a_med,breaks=30)
hist(exp(stock_char$rw_log_a_min)/exp(stock_char$rw_log_a_med),breaks=30)

hist(stock_char$rw_Smax_max,breaks=30)
hist(stock_char$rw_log_a_min,breaks=30)
hist(stock_char$rw_log_a_max/stock_char$rw_log_a_med,breaks=30)
hist(exp(stock_char$rw_log_a_min)/exp(stock_char$rw_log_a_med),breaks=30)

#hmms
par(mfrow=c(1,2))
hist(stock_char$hmm_log_a1,breaks=30)
hist(stock_char$hmm_log_a2,breaks=30)
par(mfrow=c(1,1))
hist(stock_char$hmm_log_a2-stock_char$hmm_log_a1,breaks=30)
hist((stock_char$hmm_log_a2-stock_char$hmm_log_a1)/stock_char$hmm_log_a1,breaks=30)

par(mfrow=c(1,2))
hist(log10(stock_char$hmm_Smax1),breaks=30)
hist(log10(stock_char$hmm_Smax2),breaks=30)
par(mfrow=c(1,1))
#hist(stock_char$hmm_Smax2-stock_char$hmm_Smax1,breaks=30)
hist(stock_char$hmm_Smax2/stock_char$hmm_Smax1,breaks=30)



#Stan fits####
log_a=NA
Smax=NA
sd_R=NA
AR1=NA
rw_log_a_0=NA
rw_log_a_f=NA
rw_log_a_trend=NA
rw_log_a_med=NA
rw_log_a_min=NA
rw_log_a_max=NA
rw_Smax_med=NA
rw_Smax_min=NA
rw_Smax_max=NA
rw_Smax_0=NA
rw_Smax_f=NA
rw_Smax_trend=NA
hmm_log_a=matrix(ncol=2,nrow=nrow(stock_info_filtered))
colnames(hmm_log_a)=c('hmm_log_a1','hmm_log_a2')
hmm_Smax=matrix(ncol=2,nrow=nrow(stock_info_filtered))
colnames(hmm_Smax)=c('hmm_Smax1','hmm_Smax2')



#Define models (helps prevent crashing)
m1f=sr_mod2(type='static',ac = FALSE,par='n',lfo =F)
m2f=sr_mod2(type='static',ac = TRUE,par='n',lfo=F)
m3f=sr_mod2(type='rw',par='a',lfo=F)
m4f=sr_mod2(type='rw',par='b',lfo=F)
m5f=sr_mod2(type='rw',par='both',lfo=F)
m6f=sr_mod2(type='hmm',par='a',lfo=F)
m7f=sr_mod2(type='hmm',par='b',lfo=F)
m8f=sr_mod2(type='hmm',par='both',lfo=F)


for(i in 1:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
  s<- s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  
  f1 = rstan::sampling(m1f, 
                       data = list(N=nrow(s),
                                   L=max(s$broodyear)-min(s$broodyear)+1,
                                   ii=s$broodyear-min(s$broodyear)+1,
                                   R_S=s$logR_S,
                                   S=s$spawners),
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 600)
  
  d_1=rstan::extract(f1)
  log_a[i]=median(d_1$log_a)
  Smax[i]=median(d_1$S_max)
  sd_R[i]=median(d_1$sigma)
  
  #model 2 - static autocorrelated Ricker
  f2 = rstan::sampling(m2f, 
                       data = list(N=nrow(s),
                                   L=max(s$broodyear)-min(s$broodyear)+1,
                                   ii=s$broodyear-min(s$broodyear)+1,
                                   R_S=s$logR_S,
                                   S=s$spawners),
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
  
  d_2=rstan::extract(f2)
  AR1[i]=median(d_2$rho)
  
}

stock_char=cbind(stock_info_filtered,log_a,Smax,sd_R,AR1)
write.csv(stock_char,here('outputs','static_stock_characteristics.csv'))


##
#



#model 2 - static autocorrelated Ricker
f2 = rstan::sampling(m2f, 
                     data = list(N=nrow(s),
                                 L=max(s$broodyear)-min(s$broodyear)+1,
                                 ii=s$broodyear-min(s$broodyear)+1,
                                 R_S=s$logR_S,
                                 S=s$spawners),
                     control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)

d_2=rstan::extract(f2)
AR1[i]=median(d_2$rho)
#model 3 - dynamic productivity Ricker
f3 = rstan::sampling(m3f, 
                     data = list(N=nrow(s),
                                 L=max(s$broodyear)-min(s$broodyear)+1,
                                 ii=s$broodyear-min(s$broodyear)+1,
                                 R_S=s$logR_S,
                                 S=s$spawners),
                     control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
d_3=rstan::extract(f3)
rw_log_a_med[i]=median(apply(d_3$log_a,2,median))
rw_log_a_min[i]=min(apply(d_3$log_a,2,median))
rw_log_a_max[i]=max(apply(d_3$log_a,2,median))
rw_log_a_trend[i]=((apply(d_3$log_a,2,median)[nrow(s)]-apply(d_3$log_a,2,median)[1])/apply(d_3$log_a,2,median)[1])/nrow(s)
#model 4 - dynamic capacity Ricker
f4 = rstan::sampling(m4f, 
                     data = list(N=nrow(s),
                                 L=max(s$broodyear)-min(s$broodyear)+1,
                                 ii=s$broodyear-min(s$broodyear)+1,
                                 R_S=s$logR_S,
                                 S=s$spawners),
                     control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
d_4=rstan::extract(f4)
rw_Smax_med[i]=median(apply(d_4$S_max,2,median))
rw_Smax_min[i]=min(apply(d_4$S_max,2,median))
rw_Smax_max[i]=max(apply(d_4$S_max,2,median))
rw_Smax_trend[i]=((apply(d_4$S_max,2,median)[nrow(s)]-apply(d_4$S_max,2,median)[1])/apply(d_4$S_max,2,median)[1])/nrow(s)

#model 6 - productivity regime shift - 2 regimes
f6 = rstan::sampling(m6f, 
                     data = list(N=nrow(s),
                                 R_S=s$logR_S,
                                 S=s$spawners,
                                 K=2,
                                 alpha_dirichlet=rep(1,2)), #prior for state transition probabilities (this makes them equal)
                     control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
d_6=rstan::extract(f6)
hmm_log_a[i,]=apply(d_6$log_a,2,median)

#model 7 - capacity regime shift
f7 = rstan::sampling(m7f, 
                     data = list(N=nrow(s),
                                 R_S=s$logR_S,
                                 S=s$spawners,
                                 K=2,
                                 alpha_dirichlet=rep(1,2)), #prior for state transition probabilities (this makes them equal)
                     control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)

d_7=rstan::extract(f7)
hmm_Smax[i,]=apply(d_7$S_max,2,median)
