rm(list=ls())
library(here);library(dplyr);library(ggplot2)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_may2023.csv'))
stock_info<- read.csv(here('data','filtered datasets','all_stocks_info_may2023.csv'))

#source(here('code','samEst code','stan_functions.R'))
#source(here('code','samEst code','lfo_stan_functions.R'))
#source(here('code','samEst code','lfo_TMB_functions.R'))
##source(here('code','samEst code','TMB_functions.R'))
#source(here('code','samEst code','util_functions.R'))
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)

library(samEst)
options(mc.cores = parallel::detectCores())
#1

stock_info_filtered=subset(stock_info,n.years>=16) #252 stocks
stock_info_filtered$stock.name=gsub('/','_',stock_info_filtered$stock.name)
stock_info_filtered$stock.name=gsub('&','and',stock_info_filtered$stock.name)

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock.id)) #267
stock_info_filtered$stock.id2=seq(1:nrow(stock_info_filtered))

if(any(stock_dat2$spawners==0)){stock_dat2$spawners=stock_dat2$spawners+1;stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)}
if(any(stock_dat2$recruits==0)){stock_dat2$recruits=stock_dat2$recruits+1;stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)}
stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)
stock_dat2=stock_dat2[complete.cases(stock_dat2$logR_S),]

#AIC/BIC/LOO stan####
#Define models (helps prevent crashing)
m1f=samEst::sr_mod(type='static',ac = FALSE,par='n',lfo =F)
m2f=samEst::sr_mod(type='static',ac = TRUE,par='n',lfo=F)
m3f=samEst::sr_mod(type='rw',par='a',lfo=F)
m4f=samEst::sr_mod(type='rw',par='b',lfo=F)
m5f=samEst::sr_mod(type='rw',par='both',lfo=F)
m6f=samEst::sr_mod(type='hmm',par='a',lfo=F)
m7f=samEst::sr_mod(type='hmm',par='b',lfo=F)
m8f=samEst::sr_mod(type='hmm',par='both',lfo=F)

i=176 #cowichan chinook

s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
s<- s[complete.cases(s$spawners),]

df=data.frame(S=s$spawners,R=s$recruits,by=s$broodyear)

f3 = rstan::sampling(m3f, 
                     data = list(N=nrow(s),
                                 L=max(s$broodyear)-min(s$broodyear)+1,
                                 ii=s$broodyear-min(s$broodyear)+1,
                                 R_S=s$logR_S,
                                 S=s$spawners),
                     control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 4, iter = 700)



d3=rstan::extract(f3)

Smsy_df=data.frame(by=df$by,med=apply(d3$S_msy,2,median),l90=apply(d3$S_msy,2,quantile,0.1),u90=apply(d3$S_msy,2,quantile,0.9))

ggplot2::ggplot(Smsy_df, aes(by,med)) +
  geom_line(aes(x=by,y=med),linewidth=1.3)+
  geom_point(aes(colour = by),size=4) +
  scale_colour_viridis_c(name='Year')+
  geom_ribbon(aes(ymin =l90, ymax =u90), alpha = 0.2)+
  xlab("Year") + 
  ylab("Smsy")+
  theme_classic(14)+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))

print(Smsy_df[c(nrow(Smsy_df)-10):nrow(Smsy_df),])