library(here);library(dplyr);library(ggplot2)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation2023-10-18.csv'))
stock_info<- read.csv(here('data','filtered datasets','stock_info2023-10-18.csv'))

#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)

library(samEst)
#(mc.cores = parallel::detectCores())
#1


###Load in data####
#Remove stocks with less than 15 years of recruitment data
stock_info_filtered=subset(stock_info,n.years>=15) #252 stocks
stock_info_filtered$stock.name=gsub('/','_',stock_info_filtered$stock.name)
stock_info_filtered$stock.name=gsub('&','and',stock_info_filtered$stock.name)

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock.id)) #284
stock_info_filtered$stock.id2=seq(1:nrow(stock_info_filtered))

if(any(stock_dat2$spawners==0)){stock_dat2$spawners=stock_dat2$spawners+1;stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)}
if(any(stock_dat2$recruits==0)){stock_dat2$recruits=stock_dat2$recruits+1;stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)}
stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)
stock_dat2=stock_dat2[complete.cases(stock_dat2$logR_S),]

#Create directories for outputs
for(i in 1:nrow(stock_info_filtered)){
  dir.create(here('outputs','sr plots','TMB',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep='')))
}

prod_mat=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_info_filtered$end)-min(stock_info_filtered$begin)+1)
colnames(prod_mat)=seq(min(stock_info_filtered$begin),max(stock_info_filtered$end))
umsy_mat=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_info_filtered$end)-min(stock_info_filtered$begin)+1)
colnames(umsy_mat)=seq(min(stock_info_filtered$begin),max(stock_info_filtered$end))
smsy_mat=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_info_filtered$end)-min(stock_info_filtered$begin)+1)
colnames(smsy_mat)=seq(min(stock_info_filtered$begin),max(stock_info_filtered$end))
static_df=data.frame(stock=stock_info_filtered$stock.name,log_a=NA,smax=NA,umsy=NA,smsy=NA,sigma=NA,rho=NA,log_a.ac=NA,smax.ac=NA,umsy.ac=NA,smsy.ac=NA,sigma.ac=NA)

LLdf=matrix(ncol=3,nrow=nrow(stock_info_filtered))
AICdf=matrix(ncol=3,nrow=nrow(stock_info_filtered))
BICdf=matrix(ncol=3,nrow=nrow(stock_info_filtered))

for(u in 1:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[u])
  s=s[complete.cases(s$logR_S),]
  if(any(s$spawners==0)){s$spawners=s$spawners+1;s$logR_S=log(s$recruits/s$spawners)}
  if(any(s$recruits==0)){s$recruits=s$recruits+1;s$logR_S=log(s$recruits/s$spawners)}
  s<- s[complete.cases(s$spawners),]
  
  df <- data.frame(by=s$broodyear,
                   S=s$spawners,
                   R=s$recruits,
                   logRS=s$logR_S)
  
  #lfo comparison
  lfostatic<-samEst::tmb_mod_lfo_cv(data=df,model='static', L=10)
  lfoac <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='staticAC', L=10),error = function(e) {lfoac=list(lastparam=rep(-999,length(lfoac$lastparam)))})
  lfoalpha <- tryCatch(samEst::tmb_mod_lfo_cv(data=df,model='rw_a', siglfo="obs", L=10),error = function(e) {lfoalpha=list(lastparam=rep(-999,length(lfoac$lastparam)), 
                                                                                                                           last3param=rep(-999,length(lfoac$lastparam)))})
  TMBstatic <- ricker_TMB(data=df)
  TMBac <- ricker_TMB(data=df, AC=TRUE)
  TMBtva <- tryCatch(samEst::ricker_rw_TMB(data=df,tv.par='a'),error = function(e) {TMBtva=list(conv_problem=TRUE)})
  
  prod_mat[u,match(s$broodyear,colnames(prod_mat))]=TMBtva$alpha
  umsy_mat[u,match(s$broodyear,colnames(prod_mat))]=TMBtva$umsy
  smsy_mat[u,match(s$broodyear,colnames(prod_mat))]=TMBtva$Smsy
  static_df[u,2]=TMBstatic$alpha
  static_df[u,3]=TMBstatic$Smax
  static_df[u,4]=TMBstatic$umsy
  static_df[u,5]=TMBstatic$Smsy
  static_df[u,6]=TMBstatic$sig
  static_df[u,7]=TMBac$rho
  static_df[u,8]=TMBac$alpha
  static_df[u,9]=TMBac$Smax
  static_df[u,10]=TMBac$umsy
  static_df[u,11]=TMBac$Smsy
  static_df[u,12]=TMBac$sigar
  
  LL<-rbind(lfostatic$lastparam,lfoac$lastparam,
              lfoalpha$lastparam)
  LL<- ifelse(is.infinite(LL)==T,-999,LL)
  LLdf[u,]=apply(LL,1,sum)

  AIC<-c(TMBstatic$AICc,TMBac$AICc,TMBtva$AICc)
  AICdf[u,]=AIC
  BIC<-c(TMBstatic$BIC,TMBac$BIC,TMBtva$BIC)
  BICdf[u,]=BIC
}

write.csv(LLdf,here('outputs','TMB','lfo_df.csv'))
write.csv(AICdf,here('outputs','TMB','aic_df.csv'))
write.csv(BICdf,here('outputs','TMB','bic_df.csv'))
write.csv(as.data.frame(prod_mat),here('outputs','TMB','log_alpha.csv'))
write.csv(as.data.frame(umsy_mat),here('outputs','TMB','umsy.csv'))
write.csv(as.data.frame(smsy_mat),here('outputs','TMB','smsy.csv'))
write.csv(as.data.frame(static_df),here('outputs','TMB','static_df.csv'))



#plots
sp_cols=cbind(c('#6f9ba7',
                '#316024',
                '#295f8b',
                '#a2450c',
                '#811b0a'),c('Chinook','Chum','Coho','Pink','Sockeye'))

mod_col=cbind(c("azure3","deepskyblue4",'darkgreen'),c(1,2,3))

#
#model selection####
library(cowplot)


colnames(AICdf)=c('a1','a2','a3')
mod_weights=t(apply(AICdf,1,samEst::model_weights,form='AIC'))
colnames(mod_weights)=c('w1','w2','w3')
names(mod_weights)
aic=cbind(stock_info_filtered,AICdf,mod_weights)
aic$best_mod=apply(aic[,18:20],1,which.min)
aic$lat=round(aic$lat,2)
aic$lon=round(aic$lon,2)

dAIC=  t(apply(AICdf,1,function(x) x-min(x)))
aic$best_mod2=aic$best_mod
aic$best_mod2=ifelse(aic$best_mod==2&dAIC[,1]<2,1,aic$best_mod2)
aic$best_mod2=ifelse(aic$best_mod==3&dAIC[,2]<2,1,aic$best_mod2)

table(aic$best_mod)
round(table(aic$best_mod)/nrow(aic),2)
aic$state2=aic$state
aic=aic %>% mutate(state2=ifelse(state=='OR','OR-WA',state2),
                   state2=ifelse(state=='WA','OR-WA',state2))
aic$species=gsub('-.*','',aic$species)

sp = aic %>% group_by(species,best_mod) %>% summarize(n=n(),.groups='keep')
sp1 = aic %>% group_by(species) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$species,sp1$species)]

ggplot(sp, aes(fill=factor(best_mod), y=prop, x=species))+
  scale_fill_manual(values=mod_col[match(levels(factor(sp$best_mod)),mod_col[,2])],name='Model')+ 
  ggtitle("Top-ranked model (TMB AIC)")+
  geom_bar(position="stack", stat="identity",color='white') +
  geom_text(label=sp$n, position = position_stack(vjust = 0.5),colour='white',size=2.2)+
  ylim(0,1)+
  theme_minimal() +
  xlab('')+
  ylab('Prop. populations')+ 
  theme(text=element_text(size=14),axis.text=element_text(size=14),axis.line = element_line(colour = "black", 
                                                                                            size = 1, linetype = "solid"),plot.title = element_text(size=14))+
  annotation_custom(grid::textGrob(sp1$n[1]), xmin = 1, xmax = 1, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[2]), xmin = 2, xmax = 2, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[3]), xmin = 3, xmax = 3, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[4]), xmin = 4, xmax = 4, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[6]), xmin = 6, xmax = 6, ymin = 1.03, ymax = 1.03)

sp_state = aic %>% group_by(state2,species,best_mod) %>% summarize(n=n(),.groups='keep')
sp_state1 = aic %>% group_by(state2,species) %>% summarize(n=n(),.groups='keep')
sp_state$prop = sp_state$n/sp_state1$n[match(paste(sp_state$state2,sp_state$species,sep='_'),paste(sp_state1$state2,sp_state1$species,sep='_'))]
sp_state2 = aic %>% group_by(state2,best_mod) %>% summarize(n=n(),.groups='keep')
sp_state3 = aic %>% group_by(state2) %>% summarize(n=n())
sp_state2$prop = sp_state2$n/sp_state3$n[match(paste(sp_state2$state2,sp_state2$species,sep='_'),paste(sp_state3$state2,sp_state3$species,sep='_'))]

sp_ak=subset(sp_state,state2=='AK')
sp1_ak=subset(sp_state1,state2=='AK')
sp_bc=subset(sp_state,state2=='BC')
sp1_bc=subset(sp_state1,state2=='BC')
sp_wa=subset(sp_state,state2=='OR-WA')
sp1_wa=subset(sp_state1,state2=='OR-WA')


wa_pl=ggplot(sp_wa, aes(fill=factor(best_mod), y=prop, x=species))+
  scale_fill_manual(values=mod_col[match(levels(factor(sp_wa$best_mod)),mod_col[,2])],name='Model')+ 
  ggtitle("Oregon/Washington")+
  geom_bar(position="stack", stat="identity") +
  geom_text(label=sp_wa$n, position = position_stack(vjust = 0.5),colour='white',size=2.2)+
  ylim(0,1.05)+
  theme_minimal() +
  xlab('')+
  ylab('')+ 
  theme(text=element_text(size=10),axis.text=element_text(size=10),axis.line = element_line(colour = "black", 
                                                                                            size = 1, linetype = "solid"),plot.title = element_text(size=12))+
  annotation_custom(grid::textGrob(sp1_wa$n[1],gp=grid::gpar(fontsize=8)), xmin = 1, xmax = 1, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sp1_wa$n[2],gp=grid::gpar(fontsize=8)), xmin = 2, xmax = 2, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sp1_wa$n[3],gp=grid::gpar(fontsize=8)), xmin = 3, xmax = 3, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sp1_wa$n[4],gp=grid::gpar(fontsize=8)), xmin = 4, xmax = 4, ymin = 1.05, ymax = 1.05)


bc_pl=ggplot(sp_bc, aes(fill=factor(best_mod), y=prop, x=species))+
  scale_fill_manual(values=mod_col[match(levels(factor(sp_bc$best_mod)),mod_col[,2])],name='Model')+ 
  ggtitle("British Columbia")+
  geom_bar(position="stack", stat="identity") +
  geom_text(label=sp_bc$n, position = position_stack(vjust = 0.5),colour='white',size=2.2)+
  ylim(0,1.05)+
  theme_minimal() +
  xlab('')+
  ylab('Prop. populations')+ 
  theme(text=element_text(size=10),axis.text=element_text(size=10),axis.line = element_line(colour = "black", 
                                                                                            size = 1, linetype = "solid"),plot.title = element_text(size=12))+
  annotation_custom(grid::textGrob(sp1_bc$n[1],gp=grid::gpar(fontsize=8)), xmin = 1, xmax = 1, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sp1_bc$n[2],gp=grid::gpar(fontsize=8)), xmin = 2, xmax = 2, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sp1_bc$n[3],gp=grid::gpar(fontsize=8)), xmin = 3, xmax = 3, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sp1_bc$n[4],gp=grid::gpar(fontsize=8)), xmin = 4, xmax = 4, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sp1_bc$n[5],gp=grid::gpar(fontsize=8)), xmin = 5, xmax = 5, ymin = 1.05, ymax = 1.05)

ak_pl=ggplot(sp_ak, aes(fill=factor(best_mod), y=prop, x=species))+
  scale_fill_manual(values=mod_col[match(levels(factor(sp_ak$best_mod)),mod_col[,2])],name='Model')+ 
  ggtitle("Alaska")+
  geom_bar(position="stack", stat="identity") +
  geom_text(label=sp_ak$n, position = position_stack(vjust = 0.5),colour='white',size=2.2)+
  ylim(0,1.05)+
  theme_minimal() +
  xlab('')+
  ylab('')+ 
  theme(text=element_text(size=10),axis.text=element_text(size=10),axis.line = element_line(colour = "black", 
                                                                                            size = 1, linetype = "solid"),plot.title = element_text(size=12))+
  annotation_custom(grid::textGrob(sp1_ak$n[1],gp=grid::gpar(fontsize=8)), xmin = 1, xmax = 1, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sp1_ak$n[2],gp=grid::gpar(fontsize=8)), xmin = 2, xmax = 2, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sp1_ak$n[3],gp=grid::gpar(fontsize=8)), xmin = 3, xmax = 3, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sp1_ak$n[4],gp=grid::gpar(fontsize=8)), xmin = 4, xmax = 4, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sp1_ak$n[5],gp=grid::gpar(fontsize=8)), xmin = 5, xmax = 5, ymin = 1.05, ymax = 1.05)

total_pl=ggplot(sp_state2, aes(fill=factor(best_mod), y=prop, x=state2))+
  scale_fill_manual(values=mod_col[match(levels(factor(sp_state2$best_mod)),mod_col[,2])],name='Model')+
  ggtitle("Total")+
  geom_bar(position="stack", stat="identity") +
  geom_text(label=sp_state2$n, position = position_stack(vjust = 0.5),colour='white',size=2.2)+
  ylim(0,1.05)+
  theme_minimal() +
  xlab('')+
  ylab('Prop. populations')+ 
  theme(text=element_text(size=10),axis.text=element_text(size=10),axis.line = element_line(colour = "black", 
                                                                                            size = 1, linetype = "solid"),plot.title = element_text(size=12))+
  annotation_custom(grid::textGrob(sum(sp1_wa$n),gp=grid::gpar(fontsize=8)), xmin = 1, xmax = 1, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sum(sp1_bc$n),gp=grid::gpar(fontsize=8)), xmin = 2, xmax = 2, ymin = 1.05, ymax = 1.05)+
  annotation_custom(grid::textGrob(sum(sp1_ak$n),gp=grid::gpar(fontsize=8)), xmin = 3, xmax = 3, ymin = 1.05, ymax = 1.05)

pg1=plot_grid(
  wa_pl+ theme(legend.position="none"),
  bc_pl+ theme(legend.position="none"),
  ak_pl+ theme(legend.position="none"),
  total_pl+ theme(legend.position="none"),
  ncol=2,nrow=2,labels=c("A","B","C","D"))

legend= get_legend(total_pl)

plot_grid(pg1,legend,rel_widths = c(3,.5))



#
move_avg_mat<- function(x, lag = 2) {             # Create user-defined function
  m=NA
  se=NA
  for(t in c(lag+1):c(ncol(x)-lag)){
    x1=x[,c(t-lag):c(t+lag)] 
    m[t]=mean(as.vector(na.omit(x1)))
    se[t]=sd(as.vector(na.omit(x1)))/sqrt(length(as.vector(na.omit(x1))))
  }
  return(data.frame(m,se))
}
#

#productivity####
m_prod=move_avg_mat(exp(prod_mat))
plot(m_prod$m~colnames(prod_mat)[1:72],type='n',bty='l',xlab='Brood Year',ylab='Max. Productivity (Recruits/Spawner)')
for(i in 1:nrow(prod_mat)){
  lines(prod_mat[i,]~colnames(prod_mat),col=adjustcolor('darkgray',alpha.f=0.05))
}
lines(m_prod$m~colnames(prod_mat)[1:72])
prod_matr=exp(prod_mat)

plot(rep(0,ncol(prod_matr))~colnames(prod_matr)[1:72],ylim=c(min(na.omit(as.vector(prod_matr))),25),ylab='Max. recruits/spawner',type='n',bty='l',xlab='Brood cohort year')
for(i in 1:nrow(prod_matr)){
  lines(prod_matr[i,]~colnames(prod_matr),col=adjustcolor('darkgray',alpha.f=0.3))
}
m.prod.5f=move_avg_mat(prod_matr);m.prod.5f$y=colnames(prod_matr)[1:c(ncol(prod_matr)-2)]
m.prod.5f=m.prod.5f[complete.cases(m.prod.5f[,1]),]
x<- c(seq(m.prod.5f$y[1],m.prod.5f$y[nrow(m.prod.5f)]), rev(seq(m.prod.5f$y[1],m.prod.5f$y[nrow(m.prod.5f)])))
y<- c(m.prod.5f[,1]-2*m.prod.5f[,2], rev(m.prod.5f[,1]+2*m.prod.5f[,2]))
polygon(x, y, col = adjustcolor('navy', alpha = 0.3), border=NA) # Add uncertainty polygon
lines(m.prod.5f[,1]~seq(min(x),max(x)),lwd=2,col='navy')


#10 year v longterm - u
pr10=NA;se=NA
for(i in 1:nrow(prod_matr)){
  pr=prod_matr[i,complete.cases(prod_matr[i,])]
  pr10[i]=mean(pr[c(length(pr)-9):length(pr)])
  se[i]=sd(pr[c(length(pr)-9):length(pr)])/sqrt(10)
}
plot(pr10~exp(static_df$log_a.ac),ylim=c(0,25),xlim=c(0,25),type='n',xlab='long term productivity (alpha)',ylab='mean productivity (last 10 years)')
lines(seq(0:25)~seq(0:25))
lines(seq(0:25)*0.75~seq(0:25),lty=5,col='goldenrod')
lines(seq(0:25)*1.25~seq(0:25),lty=5,col='goldenrod')
lines(seq(0:25)*0.5~seq(0:25),lty=5,col='darkorange')
lines(seq(0:25)*1.5~seq(0:25),lty=5,col='darkorange')
lines(seq(0:25)*0.25~seq(0:25),lty=5,col='darkred')
lines(seq(0:25)*1.75~seq(0:25),lty=5,col='darkred')
points(pr10[aic$species=='Chinook']~exp(static_df$log_a[aic$species=='Chinook']),cex=1.5,bg=adjustcolor(sp_cols[1],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Chum']~exp(static_df$log_a[aic$species=='Chum']),cex=1.5,bg=adjustcolor(sp_cols[2],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Coho']~exp(static_df$log_a[aic$species=='Coho']),cex=1.5,bg=adjustcolor(sp_cols[3],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Pink']~exp(static_df$log_a[aic$species=='Pink']),cex=1.5,bg=adjustcolor(sp_cols[4],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Sockeye']~exp(static_df$log_a[aic$species=='Sockeye']),cex=1.5,bg=adjustcolor(sp_cols[5],alpha.f = 0.6),pch=21)


m3=prod_matr[aic$best_mod2==3,]
pr10=NA;se=NA
for(i in 1:nrow(m3)){
  pr=m3[i,complete.cases(m3[i,])]
  pr10[i]=mean(pr[c(length(pr)-9):length(pr)])
  se[i]=sd(pr[c(length(pr)-9):length(pr)])/sqrt(10)
}
plot(pr10~exp(static_df$log_a.ac[aic$best_mod2==3]),ylim=c(0,25),xlim=c(0,25),type='n',xlab='long term productivity (alpha)',ylab='mean productivity (last 10 years)')
lines(seq(0:25)~seq(0:25))
lines(seq(0:25)*0.75~seq(0:25),lty=5,col='goldenrod')
lines(seq(0:25)*1.25~seq(0:25),lty=5,col='goldenrod')
lines(seq(0:25)*0.5~seq(0:25),lty=5,col='darkorange')
lines(seq(0:25)*1.5~seq(0:25),lty=5,col='darkorange')
lines(seq(0:25)*0.25~seq(0:25),lty=5,col='darkred')
lines(seq(0:25)*1.75~seq(0:25),lty=5,col='darkred')
points(pr10[aic$species=='Chinook']~exp(static_df$log_a[aic$species=='Chinook']),cex=1.5,bg=adjustcolor(sp_cols[1],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Chum']~exp(static_df$log_a[aic$species=='Chum']),cex=1.5,bg=adjustcolor(sp_cols[2],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Coho']~exp(static_df$log_a[aic$species=='Coho']),cex=1.5,bg=adjustcolor(sp_cols[3],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Pink']~exp(static_df$log_a[aic$species=='Pink']),cex=1.5,bg=adjustcolor(sp_cols[4],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Sockeye']~exp(static_df$log_a[aic$species=='Sockeye']),cex=1.5,bg=adjustcolor(sp_cols[5],alpha.f = 0.6),pch=21)



#umsy
pr10=NA;se=NA
for(i in 1:nrow(umsy_mat)){
  pr=umsy_mat[i,complete.cases(umsy_mat[i,])]
  pr10[i]=mean(pr[c(length(pr)-9):length(pr)])
  se[i]=sd(pr[c(length(pr)-9):length(pr)])/sqrt(10)
}
plot(pr10~static_df$umsy.ac,ylim=c(0,1),xlim=c(0,1),type='n',xlab='long term productivity (alpha)',ylab='mean productivity (last 10 years)')
lines(seq(0,1)~seq(0,1))
lines(seq(0,1)*0.75~seq(0,1),lty=5,col='goldenrod')
lines(seq(0,1)*1.25~seq(0,1),lty=5,col='goldenrod')
lines(seq(0,1)*0.5~seq(0,1),lty=5,col='darkorange')
lines(seq(0,1)*1.5~seq(0,1),lty=5,col='darkorange')
lines(seq(0,1)*0.25~seq(0,1),lty=5,col='darkred')
lines(seq(0,1)*1.75~seq(0,1),lty=5,col='darkred')
points(pr10[aic$species=='Chinook']~static_df$umsy.ac[aic$species=='Chinook'],cex=1.5,bg=adjustcolor(sp_cols[1],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Chum']~static_df$umsy.ac[aic$species=='Chum'],cex=1.5,bg=adjustcolor(sp_cols[2],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Coho']~static_df$umsy.ac[aic$species=='Coho'],cex=1.5,bg=adjustcolor(sp_cols[3],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Pink']~static_df$umsy.ac[aic$species=='Pink'],cex=1.5,bg=adjustcolor(sp_cols[4],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Sockeye']~static_df$umsy.ac[aic$species=='Sockeye'],cex=1.5,bg=adjustcolor(sp_cols[5],alpha.f = 0.6),pch=21)



m3=prod_matr[aic$best_mod2==3,]
pr10=NA;se=NA
for(i in 1:nrow(m3)){
  pr=m3[i,complete.cases(m3[i,])]
  pr10[i]=mean(pr[c(length(pr)-9):length(pr)])
  se[i]=sd(pr[c(length(pr)-9):length(pr)])/sqrt(10)
}
plot(pr10~exp(static_df$log_a.ac[aic$best_mod2==3]),ylim=c(0,25),xlim=c(0,25),type='n',xlab='long term productivity (alpha)',ylab='mean productivity (last 10 years)')
lines(seq(0:25)~seq(0:25))
lines(seq(0:25)*0.75~seq(0:25),lty=5,col='goldenrod')
lines(seq(0:25)*1.25~seq(0:25),lty=5,col='goldenrod')
lines(seq(0:25)*0.5~seq(0:25),lty=5,col='darkorange')
lines(seq(0:25)*1.5~seq(0:25),lty=5,col='darkorange')
lines(seq(0:25)*0.25~seq(0:25),lty=5,col='darkred')
lines(seq(0:25)*1.75~seq(0:25),lty=5,col='darkred')
points(pr10[aic$species=='Chinook']~exp(static_df$log_a[aic$species=='Chinook']),cex=1.5,bg=adjustcolor(sp_cols[1],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Chum']~exp(static_df$log_a[aic$species=='Chum']),cex=1.5,bg=adjustcolor(sp_cols[2],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Coho']~exp(static_df$log_a[aic$species=='Coho']),cex=1.5,bg=adjustcolor(sp_cols[3],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Pink']~exp(static_df$log_a[aic$species=='Pink']),cex=1.5,bg=adjustcolor(sp_cols[4],alpha.f = 0.6),pch=21)
points(pr10[aic$species=='Sockeye']~exp(static_df$log_a[aic$species=='Sockeye']),cex=1.5,bg=adjustcolor(sp_cols[5],alpha.f = 0.6),pch=21)



#
chi_alpha=prod_matr[aic$species=='Chinook',]
chu_alpha=prod_matr[aic$species=='Chum',]
coh_alpha=prod_matr[aic$species=='Coho',]
pi_alpha1=prod_matr[aic$species=='Pink',]
pi_alpha2=prod_matr[aic$species=='Sockeye',]
