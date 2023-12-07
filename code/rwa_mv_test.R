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

#test out mv rw model with coho
test1=subset(stock_dat2,species=='Coho')

#util function for spawners design matrix
make_design_matrix=function(x,grp){
  x2=matrix(nrow=length(x),ncol=length(unique(grp)))
  for(i in 1:length(unique(grp))){
    x2[,i]=ifelse(grp==levels(factor(grp))[i],1,0)*x
  }
  return(x2)
}

X_s=make_design_matrix(test1$spawners,grp=test1$stock)

S_mat=matrix(0,ncol=length(unique(test1$stock)),nrow=max(test1$broodyear)-min(test1$broodyear)+1)
colnames(S_mat)=unique(test1$stock);rownames(S_mat)=seq(min(test1$broodyear),max(test1$broodyear))
for(j in 1:length(unique(test1$stock))){
  sk=subset(test1,stock==unique(test1$stock)[j])
  S_mat[match(sk$broodyear,rownames(S_mat)),j]=sk$spawners
}
#stock-year index for model
stock_year=expand.grid(unique(test1$stock),unique(test1$broodyear))
stock_year= stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")

test1$stock_yr=match(paste(test1$stock,test1$broodyear,sep='_'),stock_year[,3])

library(cmdstanr)
options(mc.cores=6)
#Define models (helps prevent crashing)
file1=file.path(cmdstanr::cmdstan_path(),'sr models', "rwa_mv2.stan")
mv_mod=cmdstanr::cmdstan_model(file1)

file1.1=file.path(cmdstanr::cmdstan_path(),'sr models', "rwa_mv.stan")
mv_mod2=cmdstanr::cmdstan_model(file1.1)

smax_prior=test1%>%group_by(stock) %>%summarize(m=max(spawners))

dl = list(N=nrow(test1), #number of observations
          L=max(test1$broodyear)-min(test1$broodyear)+1, #total years
          J_i=as.numeric(factor(test1$stock)), #stock index
          J_ii=test1$stock_yr, #stock-year index
          #L_i=test1$broodyear-min(test1$broodyear)+1, #year index
          Y_i=as.numeric(factor(test1$broodyear)),
          #start_y=start_n,
          #end_y=end_n,
          #start_t=s_y$s,
          #end_t=s_y$e,
          J=length(unique(test1$stock)),
          R_S=test1$logR_S,
          S_X=X_s,
          S=X_s,
          pSmax_mean=smax_prior$m*0.5,
          pSmax_sig=smax_prior$m*0.5)

ftest2 <- mv_mod$sample(data=dl,
                        seed = 12345,
                        init=0,
                        chains = 6, 
                        iter_warmup = 200,
                        iter_sampling = 500,
                        refresh = 0,
                        adapt_delta = 0.99,
                        max_treedepth = 20)

ftest3 <- mv_mod2$sample(data=dl,
                        seed = 12345,
                        init=0,
                        chains = 6, 
                        iter_warmup = 200,
                        iter_sampling = 500,
                        refresh = 0,
                        adapt_delta = 0.99,
                        max_treedepth = 20)

par(mar = c(4, 4, 1, 1),mfrow=c(1,1))
plot(rep(0,ncol(coh_alpha))~colnames(coh_alpha),ylim=c(0,25),ylab='Max. recruits/spawner',type='n',bty='l',xlab='Brood cohort year',xaxt='n')
axis(side=1,at=seq(1950,2020,by=5))
for(i in 1:nrow(coh_alpha)){
  lines(exp(coh_alpha)[i,]~colnames(coh_alpha),col=adjustcolor('darkgray',alpha.f=0.3))
}
m.prod.5f=move_avg_mat(exp(coh_alpha));m.prod.5f$y=colnames(coh_alpha)[1:c(ncol(coh_alpha)-2)]
m.prod.5f=m.prod.5f[complete.cases(m.prod.5f[,1]),]
x<- c(seq(m.prod.5f$y[1],m.prod.5f$y[nrow(m.prod.5f)]), rev(seq(m.prod.5f$y[1],m.prod.5f$y[nrow(m.prod.5f)])))
y<- c(m.prod.5f[,1]-2*m.prod.5f[,2], rev(m.prod.5f[,1]+2*m.prod.5f[,2]))
polygon(x, y, col = adjustcolor(sp_cols[1], alpha = 0.3), border=NA) # Add uncertainty polygon
lines(m.prod.5f[,1]~seq(min(x),max(x)),lwd=2,col=sp_cols[3])


par(mar = c(4, 4, 1, 1),mfrow=c(1,1))
plot(rep(0,ncol(coh_alpha))~colnames(coh_alpha),ylim=c(0,50),ylab='Max. recruits/spawner',type='n',bty='l',xlab='Brood cohort year',xaxt='n')
axis(side=1,at=seq(1950,2020,by=5))

log_a_t=ftest3$draws(variables='log_a_t',format='draws_matrix')
log_a_t2=apply(log_a_t,2,median)
for(i in 1:nrow(coh_alpha)){
  l=log_a_t2[grepl(paste(',',i,']',sep=''),names(log_a_t2))]
  lines(exp(l)~seq(min(test1$broodyear),max(test1$broodyear)),col=adjustcolor('darkgray',alpha.f=0.3))
}
log_a_g=ftest3$draws(variables='log_a_g',format='draws_matrix')
x<- c(seq(min(test1$broodyear),max(test1$broodyear)), rev(seq(min(test1$broodyear),max(test1$broodyear))))
y<- c(apply(exp(log_a_g),2,quantile,0.025), rev(apply(exp(log_a_g),2,quantile,0.975)))
polygon(x, y, col = adjustcolor(sp_cols[1], alpha = 0.3), border=NA) # Add uncertainty polygon
lines(apply(exp(log_a_g),2,median)~seq(min(test1$broodyear),max(test1$broodyear)),lwd=2,col=sp_cols[3])

par(mar = c(4, 4, 1, 1),mfrow=c(1,1))
plot(rep(0,ncol(coh_alpha))~colnames(coh_alpha),ylim=c(0,35),ylab='Max. recruits/spawner',type='n',bty='l',xlab='Brood cohort year',xaxt='n')
axis(side=1,at=seq(1950,2020,by=5))

log_a_t=ftest2$draws(variables='log_a_t',format='draws_matrix')
log_a_t2=apply(log_a_t,2,median)
log_a_g=ftest2$draws(variables='log_a_g',format='draws_matrix')
log_a_g2=apply(log_a_g,2,median)
for(i in 1:nrow(coh_alpha)){
  l=log_a_g2+log_a_t2[grepl(paste(',',i,']',sep=''),names(log_a_t2))]
  lines(exp(l)~seq(min(test1$broodyear),max(test1$broodyear)),col=adjustcolor('darkgray',alpha.f=0.3))
}
x<- c(seq(min(test1$broodyear),max(test1$broodyear)), rev(seq(min(test1$broodyear),max(test1$broodyear))))
y<- c(apply(exp(log_a_g),2,quantile,0.025), rev(apply(exp(log_a_g),2,quantile,0.975)))
polygon(x, y, col = adjustcolor(sp_cols[1], alpha = 0.3), border=NA) # Add uncertainty polygon
lines(apply(exp(log_a_g),2,median)~seq(min(test1$broodyear),max(test1$broodyear)),lwd=2,col=sp_cols[3])
