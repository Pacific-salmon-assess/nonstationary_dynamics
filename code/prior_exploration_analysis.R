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

library(cmdstanr)
file1=file.path(cmdstanr::cmdstan_path(),'sr models', "m1f.stan")
m1=cmdstanr::cmdstan_model(file1)
file1_ip=file.path(cmdstanr::cmdstan_path(),'sr models', "m1f_ip.stan")
m1_ip=cmdstanr::cmdstan_model(file1_ip)
file2=file.path(cmdstanr::cmdstan_path(),'sr models', "m2f.stan")
m2=cmdstanr::cmdstan_model(file2)
file2=file.path(cmdstanr::cmdstan_path(),'sr models', "m2f_ip.stan")
m2=cmdstanr::cmdstan_model(file2_ip)
file3=file.path(cmdstanr::cmdstan_path(),'sr models', "m3f.stan")
m3=cmdstanr::cmdstan_model(file3)
file3_ip=file.path(cmdstanr::cmdstan_path(),'sr models', "m3f_ip.stan")
m3_ip=cmdstanr::cmdstan_model(file3_ip)
file4=file.path(cmdstanr::cmdstan_path(),'sr models', "m4f.stan")
m4=cmdstanr::cmdstan_model(file4)
file4_ip=file.path(cmdstanr::cmdstan_path(),'sr models', "m4f_ip.stan")
m4_ip=cmdstanr::cmdstan_model(file4_ip)

file42=file.path(cmdstanr::cmdstan_path(),'sr models', "m4f2.stan")
m42=cmdstanr::cmdstan_model(file4)
file4_ip2=file.path(cmdstanr::cmdstan_path(),'sr models', "m4f_ip2.stan")
m4_ip2=cmdstanr::cmdstan_model(file4_ip2)


file5=file.path(cmdstanr::cmdstan_path(),'sr models', "m5f.stan")
m5=cmdstanr::cmdstan_model(file5)
file5_ip=file.path(cmdstanr::cmdstan_path(),'sr models', "m5f_ip.stan")
m5_ip=cmdstanr::cmdstan_model(file5_ip)
file6=file.path(cmdstanr::cmdstan_path(),'sr models', "m6f.stan")
m6=cmdstanr::cmdstan_model(file6)
file6_ip=file.path(cmdstanr::cmdstan_path(),'sr models', "m6f_ip.stan")
m6_ip=cmdstanr::cmdstan_model(file6_ip)
file7=file.path(cmdstanr::cmdstan_path(),'sr models', "m7f.stan")
m7=cmdstanr::cmdstan_model(file7)
file7_ip=file.path(cmdstanr::cmdstan_path(),'sr models', "m7f_ip.stan")
m7_ip=cmdstanr::cmdstan_model(file7_ip)
file8=file.path(cmdstanr::cmdstan_path(),'sr models', "m8f.stan")
m8=cmdstanr::cmdstan_model(file8)
file8_ip=file.path(cmdstanr::cmdstan_path(),'sr models', "m8f_ip.stan")
m8_ip=cmdstanr::cmdstan_model(file8_ip)

post_f3=list()
post_f3_ip=list()
post_f4=list()
post_f4_ip=list()
post_f42=list()
post_f4_ip2=list()
for(i in 31:40){

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
                  adapt_delta = 0.95,
                  max_treedepth = 15)
  
  f3_ip <- m3_ip$sample(data=dl,
                  seed = 123,
                  chains = 6, 
                  iter_warmup = 200,
                  iter_sampling = 600,
                  refresh = 0,
                  adapt_delta = 0.95,
                  max_treedepth = 15)
  
  f4 <- m4$sample(data=dl,
                  seed = 123,
                  chains = 6, 
                  iter_warmup = 200,
                  iter_sampling = 600,
                  refresh = 0,
                  adapt_delta = 0.95,
                  max_treedepth = 15)
  
  f4_ip <- m4_ip$sample(data=dl,
                        seed = 123,
                        chains = 6, 
                        iter_warmup = 200,
                        iter_sampling = 600,
                        refresh = 0,
                        adapt_delta = 0.95,
                        max_treedepth = 15)
  
  f42 <- m42$sample(data=dl,
                  seed = 123,
                  chains = 6, 
                  iter_warmup = 200,
                  iter_sampling = 600,
                  refresh = 0,
                  adapt_delta = 0.95,
                  max_treedepth = 15)
  
  f4_ip2 <- m4_ip2$sample(data=dl,
                        seed = 123,
                        chains = 6, 
                        iter_warmup = 200,
                        iter_sampling = 600,
                        refresh = 0,
                        adapt_delta = 0.95,
                        max_treedepth = 15)
  
  post_f3[[i]]=f3$draws(variable=c('log_a','S_max'),format='draws_matrix')
  post_f3_ip[[i]]=f3_ip$draws(variable=c('log_a','S_max'),format='draws_matrix')
  post_f4[[i]]=f4$draws(variable=c('log_a','S_max'),format='draws_matrix')
  post_f4_ip[[i]]=f4_ip$draws(variable=c('log_a','S_max'),format='draws_matrix')
  post_f42[[i]]=f42$draws(variable=c('log_a','S_max'),format='draws_matrix')
  post_f4_ip2[[i]]=f4_ip2$draws(variable=c('log_a','S_max'),format='draws_matrix')
}

pdf(paste(here('outputs','figures','prior test'),'/prior_test2.pdf',sep=''))
for(i in 31:40){
  
  par(mfrow=c(2,2))
  hist(log10(post_f3[[i]][,grepl('S_max',colnames(post_f3[[i]]))]),breaks=30,xlab='log10(smax)',main='wide prior-m3')
  hist(log10(post_f3_ip[[i]][,grepl('S_max',colnames(post_f3_ip[[i]]))]),breaks=30,xlab='log10(smax)',main='inf. prior-m3')
  plot(apply(post_f3[[i]][,grepl('log_a',colnames(post_f3[[i]]))],2,median),type='l',lwd=2,ylab='log_a',main=stock_info_filtered$stock.name[i])
  lines(apply(post_f3_ip[[i]][,grepl('log_a',colnames(post_f3_ip[[i]]))],2,median),lwd=2,col='red')
  text(3,par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.75,'inf.',col='red')
  
  par(mfrow=c(2,2))
  hist(post_f4[[i]][,grepl('log_a',colnames(post_f4[[i]]))],breaks=30,xlab='log_a',main='wide prior-m3')
  hist(post_f4_ip[[i]][,grepl('log_a',colnames(post_f4_ip[[i]]))],breaks=30,xlab='log_a',main='inf. prior-m3')
  hist(post_f42[[i]][,grepl('log_a',colnames(post_f42[[i]]))],breaks=30,xlab='log_a',main='wide prior-bounded')
  hist(post_f4_ip2[[i]][,grepl('log_a',colnames(post_f4_ip2[[i]]))],breaks=30,xlab='log_a',main='inf. prior-m3-bounded')

  #
  par(mfrow=c(1,1))
  plot(apply(log10(post_f4[[i]][,grepl('S_max',colnames(post_f4[[i]]))]),2,median),type='l',lwd=2,ylab='log10(smax)')
  lines(apply(log10(post_f4_ip[[i]][,grepl('S_max',colnames(post_f4_ip[[i]]))]),2,median),lwd=2,col='red')
  text(3,par('usr')[4]-(par('usr')[4]-par('usr')[3])*0.75,'inf.',col='red')
}
dev.off()
