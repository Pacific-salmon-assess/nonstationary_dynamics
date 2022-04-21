#Salmon spawner=recruit time-series; dan.greenberg@dfo-mpo.gc.ca
library(here);library(dplyr);library(ggplot2);library(ggthemes);library(ggspatial);library(sf);library(rnaturalearth);library(rnaturalearthdata);library(wesanderson)
library(nlme); library(rlist)
#plot functions
source(here('code','plot_functions.R'))

#read in data for each species
sockeye<- read.csv(here('data','sockeye','sockeye_data.csv')) 
chum<- read.csv(here('data','chum','chum_data.csv')) 
pink<- read.csv(here('data','pink','pink_data.csv')) 
coho<- read.csv(here('data','coho','coho_data_fwmar.csv')) 
chinook<- read.csv(here('data','chinook','chinook_data_totalage.csv'))
 
#Synonymize column names so I don't go crazy
colnames(sockeye)<- tolower(names(sockeye))
names(chum)<- gsub('.yr','year',names(chum)) #make the brood years equivalent
names(pink)<- gsub('.yr','year',names(pink))  #make the brood years equivalent
colnames(coho)<- tolower(names(coho))
colnames(chinook)<- tolower(names(chinook))

#Remove estimates that are flagged as unreliable ('useflag'/'use')
sockeye<- subset(sockeye,useflag==1)
chum<- subset(chum,use==1)
chinook<- subset(chinook,useflag==1)
coho<- subset(coho,useflag==1)
pink<- subset(pink,use==1)
sockeye<- subset(sockeye,useflag==1)

#Number of Stock time-series for each species (before filtering)
length(unique(sockeye$stock.id)) #51
length(unique(chum$stock.id)) #53
length(unique(pink$stock.id)) #46
length(unique(coho$stock.id)) #13
length(unique(chinook$stock)) #20

#stock info with coordinates
sockeye_info<- read.csv(here('data','sockeye','sockeye_info.csv'))
colnames(sockeye_info)<- tolower(names(sockeye_info)) 
sockeye_info<- subset(sockeye_info,stock.id %in% sockeye$stock.id) #Removes some extra stocks that do not have actual data
chum_info<- read.csv(here('data','chum','chum_info.csv'))
chum_info$lon<- -chum_info$lon #make western longitudes negative
pink_info<- read.csv(here('data','pink','pink_info.csv'))
pink_info$lon<- -pink_info$lon
coho_info<- read.csv(here('data','coho','coho_info.csv'))
colnames(coho_info)<- tolower(names(coho_info))
chinook_info<- read.csv(here('data','chinook','chinook_info.csv'))
colnames(chinook_info)<- tolower(names(chinook_info))

###time-series lengths, putting each series into a list, and estimating Spawner-Recruit fits and residual autocorrelation
chinook_list= list()
chum_list= list()
coho_list= list()
pink_list= list()
sockeye_list= list()
for(i in 1:length(unique(chinook_info$stock.id))){ 
  c<- chinook %>% subset(stock.id==unique(stock.id)[i])
  c$logR_S<- log(c$recruits/c$spawners) #log R/S to use later
  c_sr<- c[complete.cases(c$logR_S),] #only keep the years with full spawner-recruit data
  chinook_info$ts.length[i]=nrow(c_sr) #time-series length
  chinook_info$ts.start[i]=min(c_sr$broodyear) #start year
  chinook_info$ts.end[i]=max(c_sr$broodyear) #end year
  
  chinook_list[[i]]=c_sr
}

for(i in 1:length(unique(chum_info$stock.id))){
  c<- chum %>% subset(stock.id==unique(stock.id)[i])
  c$logR_S<- log(c$recruits/c$spawners)
  c_sr<- c[complete.cases(c$logR_S),]
  chum_info$ts.length[i]=nrow(c_sr)
  chum_info$ts.start[i]=min(c_sr$broodyear)
  chum_info$ts.end[i]=max(c_sr$broodyear)
  
  chum_list[[i]]=c_sr
}

for(i in 1:length(unique(coho_info$stock.id))){
  c<- coho %>% subset(stock.id==unique(stock.id)[i])
  c$logR_S<- log(c$recruits/c$spawners)
  c_sr<- c[complete.cases(c$logR_S),]
  
  coho_info$ts.length[i]=nrow(c_sr)
  coho_info$ts.start[i]=min(c_sr$broodyear)
  coho_info$ts.end[i]=max(c_sr$broodyear)
  
  coho_list[[i]]=c_sr
}

for(i in 1:length(unique(pink_info$stock.id))){
  c<- subset(pink,stock.id==unique(pink_info$stock.id)[i])
  c$odd_even<- ifelse(as.integer(c$broodyear)%%2, "odd", "even")
  c$logR_S<- log(c$recruits/c$spawners)
  c_sr<- c[complete.cases(c$logR_S),]
  
  
  if(length(unique(c$odd_even))==2){
    c_odd<- subset(c_sr,odd_even=='odd')
    c_even<- subset(c_sr,odd_even=='even')
    
    pink_info$ts.length[i]=nrow(c_odd)
    pink_info$ts.start[i]=min(c_odd$broodyear)
    pink_info$ts.end[i]=max(c_odd$broodyear)
    
    even_info<- pink_info[i,]
    even_info$ts.length=nrow(c_even)
    even_info$ts.start=min(c_even$broodyear)
    even_info$ts.end=max(c_even$broodyear)
    even_info$stock.id=max(pink_info$stock.id)+1
    even_info$stock<- paste(pink_info$stock[i],'even',sep='_')
    
    pink_info$stock[i]<- paste(pink_info$stock[i],'odd',sep='_')
    
    pink_info=rbind(pink_info,even_info)
    
    pink_list[[i]]=c_odd
    pink_list[[nrow(pink_info)]]=c_even
    
  }
  else{
    pink_info$ts.length[i]=nrow(c_sr)
    pink_info$ts.start[i]=min(c_sr$broodyear)
    pink_info$ts.end[i]=max(c_sr$broodyear)
    
    pink_list[[i]]=c_sr
  }
}

for(i in 1:length(unique(sockeye_info$stock.id))){
  c<- subset(sockeye,stock.id==unique(sockeye_info$stock.id)[i])  
  c$logR_S<- log(c$recruits/c$spawners)
  c_sr<- c[complete.cases(c$logR_S),]
  if(nrow(c_sr)==0){ 
  sockeye_info$ts.length[i]=NA
  sockeye_info$ts.start[i]=NA
  sockeye_info$ts.end[i]=NA
  }else{
    sockeye_info$ts.length[i]=nrow(c_sr)
    sockeye_info$ts.start[i]=min(c_sr$broodyear)
    sockeye_info$ts.end[i]=max(c_sr$broodyear)
    
    sockeye_list[[i]]=c_sr
  }
}

chinook_filtered<- subset(chinook_info, ts.length>=20) #9 timeseries
chum_filtered<- subset(chum_info, ts.length>=20) #50 timeseries
coho_filtered<- subset(coho_info, ts.length>=20) #9 timseries
pink_filtered<- subset(pink_info, ts.length>=20) #39 timeseries
sockeye_filtered<- subset(sockeye_info, ts.length>=20) #51 timeseries

chinook_filtered_list<- rlist::list.filter(chinook_list, length(logR_S)>=20)
chum_filtered_list<- rlist::list.filter(chum_list, length(logR_S)>=20)
coho_filtered_list<- rlist::list.filter(coho_list, length(logR_S)>=20)
pink_filtered_list<- rlist::list.filter(pink_list, length(logR_S)>=20)
sockeye_filtered_list<- rlist::list.filter(sockeye_list, length(logR_S)>=20)

for(i in 1:length(chinook_filtered_list)){
  c_sr=chinook_filtered_list[[i]]
  
  sr_fit_ac<- gls(logR_S~spawners,data=c_sr,correlation = corAR1(form = ~ 1),method='ML')
  sr_fit<- gls(logR_S~spawners,data=c_sr,method='ML')
  acf_sr<- acf(residuals(sr_fit))
  chinook_filtered$logA[i]<- sr_fit_ac$coefficients[1]
  chinook_filtered$b[i]<- -sr_fit_ac$coefficients[2]
  chinook_filtered$K[i]<- 1/-sr_fit_ac$coefficients[2]
  chinook_filtered$rmax[i]<- (exp(sr_fit_ac$coefficients[1])/-sr_fit_ac$coefficients[2])*exp(-1)
  chinook_filtered$sigma[i]<- sr_fit_ac$sigma
  chinook_filtered$ar1[i]<-intervals(sr_fit_ac)$corStruct[,2] #ML estimate of AR1 correlation
  chinook_filtered$ar1.l95[i]<-intervals(sr_fit_ac)$corStruct[,1] #l95% CI of ML estimate
  chinook_filtered$ar1.u95[i]<-intervals(sr_fit_ac)$corStruct[,3] #u95% CI of ML estimate
  chinook_filtered$p_S_belowK[i]<- nrow(subset(c,spawners<(1/-sr_fit_ac$coefficients[2])))/nrow(c)
  chinook_filtered$logA2[i]<- sr_fit$coefficients[1]
  chinook_filtered$b2[i]<- -sr_fit$coefficients[2]
  chinook_filtered$K2[i]<- 1/-sr_fit$coefficients[2]
  chinook_filtered$rmax2[i]<- (exp(sr_fit$coefficients[1])/-sr_fit$coefficients[2])*exp(-1)
  chinook_filtered$sigma2[i]<- sr_fit$sigma
  chinook_filtered$acf_1[i]<- acf_sr$acf[2]
  chinook_filtered$acf_2[i]<- acf_sr$acf[3]
  chinook_filtered$acf_3[i]<- acf_sr$acf[4]
  chinook_filtered$acf_4[i]<- acf_sr$acf[5]
  chinook_filtered$m.spawners[i]<- mean(c_sr$spawners)
  chinook_filtered$sd.spawners[i]<- sd(c_sr$spawners)
  #  plot_SR(c_sr,sr_fit_ac,path=here('outputs','figures','stock-recruit','chinook'))
}

for(i in 1:length(unique(chum_filtered$stock.id))){
  c_sr=chum_filtered_list[[i]]
  
  sr_fit_ac<- gls(logR_S~spawners,data=c_sr,correlation = corAR1(form = ~ 1),method='ML')
  sr_fit<- gls(logR_S~spawners,data=c_sr,method='ML')
  acf_sr<- acf(residuals(sr_fit))
  chum_filtered$logA[i]<- sr_fit_ac$coefficients[1]
  chum_filtered$b[i]<- -sr_fit_ac$coefficients[2]
  chum_filtered$K[i]<- 1/-sr_fit_ac$coefficients[2]
  chum_filtered$rmax[i]<- (exp(sr_fit_ac$coefficients[1])/-sr_fit_ac$coefficients[2])*exp(-1)
  chum_filtered$sigma[i]<- sr_fit_ac$sigma
  chum_filtered$ar1[i]<-intervals(sr_fit_ac)$corStruct[,2] #ML estimate of AR1 correlation
  chum_filtered$ar1.l95[i]<-intervals(sr_fit_ac)$corStruct[,1] #l95% CI of ML estimate
  chum_filtered$ar1.u95[i]<-intervals(sr_fit_ac)$corStruct[,3] #u95% CI of ML estimate
  chum_filtered$p_S_belowK[i]<- nrow(subset(c,spawners<(1/-sr_fit_ac$coefficients[2])))/nrow(c)
  chum_filtered$logA2[i]<- sr_fit$coefficients[1]
  chum_filtered$b2[i]<- -sr_fit$coefficients[2]
  chum_filtered$K2[i]<- 1/-sr_fit$coefficients[2]
  chum_filtered$rmax2[i]<- (exp(sr_fit$coefficients[1])/-sr_fit$coefficients[2])*exp(-1)
  chum_filtered$sigma2[i]<- sr_fit$sigma
  chum_filtered$acf_1[i]<- acf_sr$acf[2]
  chum_filtered$acf_2[i]<- acf_sr$acf[3]
  chum_filtered$acf_3[i]<- acf_sr$acf[4]
  chum_filtered$acf_4[i]<- acf_sr$acf[5]
  chum_filtered$m.spawners[i]<- mean(c_sr$spawners)
  chum_filtered$sd.spawners[i]<- sd(c_sr$spawners)
  
 # plot_SR(c_sr,sr_fit_ac,path=here('outputs','figures','stock-recruit','chum'))
}

for(i in 1:length(unique(coho_filtered$stock.id))){
  c_sr=coho_filtered_list[[i]]
  c_sr<- subset(c_sr,recruits>0)
  
  sr_fit_ac<- gls(logR_S~spawners,data=c_sr,correlation = corAR1(form = ~ 1),method='ML')
  sr_fit<- gls(logR_S~spawners,data=c_sr,method='ML')
  acf_sr<- acf(residuals(sr_fit))
  coho_filtered$logA[i]<- sr_fit_ac$coefficients[1]
  coho_filtered$b[i]<- -sr_fit_ac$coefficients[2]
  coho_filtered$K[i]<- 1/-sr_fit_ac$coefficients[2]
  coho_filtered$rmax[i]<- (exp(sr_fit_ac$coefficients[1])/-sr_fit_ac$coefficients[2])*exp(-1)
  coho_filtered$sigma[i]<- sr_fit_ac$sigma
  coho_filtered$ar1[i]<-intervals(sr_fit_ac)$corStruct[,2] #ML estimate of AR1 correlation
  coho_filtered$ar1.l95[i]<-intervals(sr_fit_ac)$corStruct[,1] #l95% CI of ML estimate
  coho_filtered$ar1.u95[i]<-intervals(sr_fit_ac)$corStruct[,3] #u95% CI of ML estimate
  coho_filtered$p_S_belowK[i]<- nrow(subset(c,spawners<(1/-sr_fit_ac$coefficients[2])))/nrow(c)
  coho_filtered$logA2[i]<- sr_fit$coefficients[1]
  coho_filtered$b2[i]<- -sr_fit$coefficients[2]
  coho_filtered$K2[i]<- 1/-sr_fit$coefficients[2]
  coho_filtered$rmax2[i]<- (exp(sr_fit$coefficients[1])/-sr_fit$coefficients[2])*exp(-1)
  coho_filtered$sigma2[i]<- sr_fit$sigma
  coho_filtered$acf_1[i]<- acf_sr$acf[2]
  coho_filtered$acf_2[i]<- acf_sr$acf[3]
  coho_filtered$acf_3[i]<- acf_sr$acf[4]
  coho_filtered$acf_4[i]<- acf_sr$acf[5]
  coho_filtered$m.spawners[i]<- mean(c_sr$spawners)
  coho_filtered$sd.spawners[i]<- sd(c_sr$spawners)
  
 # plot_SR(c_sr,sr_fit_ac,path=here('outputs','figures','stock-recruit','coho'))
}

for(i in 1:length(pink_filtered_list)){
  c_sr=pink_filtered_list[[i]]
  
  sr_fit_ac<- gls(logR_S~spawners,data=c_sr,correlation = corAR1(form = ~ 1),method='ML')
  sr_fit<- gls(logR_S~spawners,data=c_sr,method='ML')
  acf_sr<- acf(residuals(sr_fit))
  pink_filtered$logA[i]<- sr_fit_ac$coefficients[1]
  pink_filtered$b[i]<- -sr_fit_ac$coefficients[2]
  pink_filtered$K[i]<- 1/-sr_fit_ac$coefficients[2]
  pink_filtered$rmax[i]<- (exp(sr_fit_ac$coefficients[1])/-sr_fit_ac$coefficients[2])*exp(-1)
  pink_filtered$sigma[i]<- sr_fit_ac$sigma
  pink_filtered$ar1[i]<-intervals(sr_fit_ac)$corStruct[,2] #ML estimate of AR1 correlation
  pink_filtered$ar1.l95[i]<-intervals(sr_fit_ac)$corStruct[,1] #l95% CI of ML estimate
  pink_filtered$ar1.u95[i]<-intervals(sr_fit_ac)$corStruct[,3] #u95% CI of ML estimate
  pink_filtered$p_S_belowK[i]<- nrow(subset(c,spawners<(1/-sr_fit_ac$coefficients[2])))/nrow(c)
  pink_filtered$logA2[i]<- sr_fit$coefficients[1]
  pink_filtered$b2[i]<- -sr_fit$coefficients[2]
  pink_filtered$K2[i]<- 1/-sr_fit$coefficients[2]
  pink_filtered$rmax2[i]<- (exp(sr_fit$coefficients[1])/-sr_fit$coefficients[2])*exp(-1)
  pink_filtered$sigma2[i]<- sr_fit$sigma
  pink_filtered$acf_1[i]<- acf_sr$acf[2]
  pink_filtered$acf_2[i]<- acf_sr$acf[3]
  pink_filtered$acf_3[i]<- acf_sr$acf[4]
  pink_filtered$acf_4[i]<- acf_sr$acf[5]
  pink_filtered$m.spawners[i]<- mean(c_sr$spawners)
  pink_filtered$sd.spawners[i]<- sd(c_sr$spawners)
  
#  plot_SR(c_sr,sr_fit_ac,path=here('outputs','figures','stock-recruit','pink'))
}

for(i in 1:length(unique(sockeye_filtered$stock.id))){
  c<- subset(sockeye,stock.id==unique(sockeye_filtered$stock.id)[i])  %>% subset(useflag==1)
  c$logR_S<- log(c$recruits/c$spawners)
  c_sr<- c[complete.cases(c$logR_S),]
  if(nrow(c_sr)==0){next}
  
  sr_fit_ac<- gls(logR_S~spawners,data=c_sr,correlation = corAR1(form = ~ 1),method='ML')
  sr_fit<- gls(logR_S~spawners,data=c_sr,method='ML')
  acf_sr<- acf(residuals(sr_fit))
  sockeye_filtered$logA[i]<- sr_fit_ac$coefficients[1]
  sockeye_filtered$b[i]<- -sr_fit_ac$coefficients[2]
  sockeye_filtered$K[i]<- 1/-sr_fit_ac$coefficients[2]
  sockeye_filtered$rmax[i]<- (exp(sr_fit_ac$coefficients[1])/-sr_fit_ac$coefficients[2])*exp(-1)
  sockeye_filtered$sigma[i]<- sr_fit_ac$sigma
  sockeye_filtered$ar1[i]<-intervals(sr_fit_ac)$corStruct[,2] #ML estimate of AR1 correlation
  sockeye_filtered$ar1.l95[i]<-intervals(sr_fit_ac)$corStruct[,1] #l95% CI of ML estimate
  sockeye_filtered$ar1.u95[i]<-intervals(sr_fit_ac)$corStruct[,3] #u95% CI of ML estimate
  sockeye_filtered$p_S_belowK[i]<- nrow(subset(c,spawners<(1/-sr_fit_ac$coefficients[2])))/nrow(c)
  sockeye_filtered$logA2[i]<- sr_fit$coefficients[1]
  sockeye_filtered$b2[i]<- -sr_fit$coefficients[2]
  sockeye_filtered$K2[i]<- 1/-sr_fit$coefficients[2]
  sockeye_filtered$rmax2[i]<- (exp(sr_fit$coefficients[1])/-sr_fit$coefficients[2])*exp(-1)
  sockeye_filtered$sigma2[i]<- sr_fit$sigma
  sockeye_filtered$acf_1[i]<- acf_sr$acf[2]
  sockeye_filtered$acf_2[i]<- acf_sr$acf[3]
  sockeye_filtered$acf_3[i]<- acf_sr$acf[4]
  sockeye_filtered$acf_4[i]<- acf_sr$acf[5]
  sockeye_filtered$m.spawners[i]<- mean(c_sr$spawners)
  sockeye_filtered$sd.spawners[i]<- sd(c_sr$spawners)
  
#  plot_SR(c_sr,sr_fit_ac,path=here('outputs','figures','stock-recruit','sockeye'))
}

write.csv(chinook_filtered,file.path(here('data','filtered datasets'),'filtered_chinook_info.csv'))
write.csv(coho_filtered,file.path(here('data','filtered datasets'),'filtered_coho_info.csv'))
write.csv(chum_filtered,file.path(here('data','filtered datasets'),'filtered_chum_info.csv'))
write.csv(pink_filtered,file.path(here('data','filtered datasets'),'filtered_pink_info.csv'))
write.csv(sockeye_filtered,file.path(here('data','filtered datasets'),'filtered_sockeye_info.csv'))

#Time-series length
par(mfrow=c(3,2))
hist(chinook_filtered$ts.length, main = 'Chinook (n = 9)',col='darkgray',xlab='',ylab='Count')
hist(chum_filtered$ts.length, main = 'Chum (n = 50)',col='darkgreen',xlab='',ylab='')
hist(coho_filtered$ts.length, main = 'Coho (n = 9)',col='darkblue',xlab='',ylab='Count')
hist(pink_filtered$ts.length, main = 'Pink (n = 39)',col='darksalmon',xlab='Time-series Length (years)',ylab='')
hist(sockeye_filtered$ts.length, main = 'Sockeye (n = 51)',col='darkred',xlab='Time-series Length (years)',ylab='Count')

#autocorrelation in residual productivity
par(mfrow=c(3,2))
plot_acf(chinook_filtered)
plot_acf(chum_filtered)
plot_acf(coho_filtered)
plot_acf(pink_filtered)
plot_acf(sockeye_filtered)

summary(chinook_filtered$ar1) #mean = 0.33
summary(chum_filtered$ar1) #mean = 0.3
summary(coho_filtered$ar1) #mean = 0.17
summary(pink_filtered$ar1) #mean = -0.01
summary(sockeye_filtered$ar1) #mean = 0.36

summary(chinook_filtered$acf_2) #mean = 0.11
summary(chum_filtered$acf_2) #mean = 0.16
summary(coho_filtered$acf_2) #mean = 0.03
summary(pink_filtered$acf_2) #mean = -0.08
summary(sockeye_filtered$acf_2) #mean = 0.17

summary(chinook_filtered$acf_3) #mean = -0.04
summary(chum_filtered$acf_3) #mean = 0.03
summary(coho_filtered$acf_3) #mean = 0.05
summary(pink_filtered$acf_3) #mean = -0.02
summary(sockeye_filtered$acf_3) #mean = 0.05




#PSE data
pse_dat<- read.csv(here('data','PSF data','PSE_RS.csv'))
pse_spawners<- subset(pse_dat, parameter %in% c('Spawners'))
pse_recruits<- subset(pse_dat, parameter %in% c('Recruits'))
pse_dat2<- left_join(pse_spawners,pse_recruits,by=c('species','location','year'))
pse_dat2$stock<- paste(pse_dat2$species,pse_dat2$location,sep='_')
pse_dat2<- pse_dat2[,c(-1,-6)] #remove excess columns
colnames(pse_dat2)[4:5]<- c('spawners','recruits')

pse_dq<- read.csv(here('data','PSF data', 'PSE_data_quality.csv'))
summary(factor(pse_dq$parameter))
pse_cq<- subset(pse_dq, parameter %in% 'catch_quality')
pse_sc<- subset(pse_dq, parameter %in% 'survey_coverage')
pse_se<- subset(pse_dq, parameter %in% 'survey_execution')
pse_sq<- subset(pse_dq, parameter %in% 'survey_quality')
pse_tq<- subset(pse_dq, parameter %in% 'dq_score')

pse_dq2<- left_join(pse_tq,pse_cq %>% dplyr::select(species,location,datavalue),by=c('species','location'))
pse_dq2<- left_join(pse_dq2,pse_sc %>% dplyr::select(species,location,datavalue),by=c('species','location'))
pse_dq2<- left_join(pse_dq2,pse_se %>% dplyr::select(species,location,datavalue),by=c('species','location'))
pse_dq2<- left_join(pse_dq2,pse_sq %>% dplyr::select(species,location,datavalue),by=c('species','location'))
pse_dq2<- pse_dq2[,c(-1,-4)] #remove excess columns
colnames(pse_dq2)[c(3:7)]<- c('total_dq_score','catch_quality','survey_coverage','survey_execution','survey_quality')

pse_dq_sub1<- subset(pse_dq2,total_dq_score>=12)
pse_dq_sub2<- subset(pse_dq2,catch_quality>=3&survey_coverage>=3&survey_execution>=3&survey_quality>=3)

pse_info<- data.frame(index=NA,species=NA,location=NA,ts.length=NA,ts.start=NA,ts.end=NA)
for(i in 1:length(unique(pse_dat2$stock))){
  c<- subset(pse_dat2,stock==unique(pse_dat2$stock)[i])
  c$logR_S<- log(c$recruits/c$spawners)
  c_sr<- c[complete.cases(c$logR_S),]
  pse_info[i,1]=i
  pse_info[i,2]=unique(c$species)
  pse_info[i,3]=unique(c$location)
  pse_info[i,4]=nrow(c_sr)
  pse_info[i,5]=min(c_sr[,3])
  pse_info[i,6]=max(c_sr[,3])
}

pse_filtered<- subset(pse_info,ts.length>=20) #116 stocks
summary(factor(pse_filtered$species))
#11 chinook
#16 chum
#16 coho
#50 sockeye
#23 pink
pse_chinook<- subset(pse_filtered,species=='Chinook')
pse_chum<- subset(pse_filtered,species=='Chum')
pse_coho<- subset(pse_filtered,species=='Coho')
pse_pink<- subset(pse_filtered,species=='Pink (even)'|species=='Pink (odd)')
pse_sockeye<- subset(pse_filtered,species=='Lake sockeye'|species=='River sockeye')

par(mfrow=c(3,2))
hist(pse_chinook$ts.length, main = 'Chinook (n = 11)',col='darkgray',xlab='',ylab='Count')
hist(pse_chum$ts.length, main = 'Chum (n = 16)',col='darkgreen',xlab='',ylab='')
hist(pse_coho$ts.length, main = 'Coho (n = 16)',col='darkblue',xlab='',ylab='Count')
hist(pse_pink$ts.length, main = 'Pink (n = 23)',col='darksalmon',xlab='Time-series Length (years)',ylab='')
hist(pse_sockeye$ts.length, main = 'Sockeye (n = 50)',col='darkred',xlab='Time-series Length (years)',ylab='Count')

write.csv(pse_filtered,here('data','filtered datasets','PSE_stockrecruit_info.csv'))

chinook_final<- do.call(rbind, lapply(chinook_filtered_list, data.frame, stringsAsFactors=FALSE))
chum_final<- do.call(rbind, lapply(chum_filtered_list, data.frame, stringsAsFactors=FALSE))
coho_final<- do.call(rbind, lapply(coho_filtered_list, data.frame, stringsAsFactors=FALSE))
pink_final<- do.call(rbind, lapply(pink_filtered_list, data.frame, stringsAsFactors=FALSE))
sockeye_final<- do.call(rbind, lapply(sockeye_filtered_list, data.frame, stringsAsFactors=FALSE))
write.csv(chinook_final,here('data','filtered datasets','chinook_final.csv'))
write.csv(chum_final,here('data','filtered datasets','chum_final.csv'))
write.csv(coho_final,here('data','filtered datasets','coho_final.csv'))
write.csv(pink_final,here('data','filtered datasets','pink_final.csv'))
write.csv(sockeye_final,here('data','filtered datasets','sockeye_final.csv'))


###Curry Cunningham compilation
cc_dat<- read.csv(here('data','cunningham dataset','AK-WCoast-Salmon-SR.csv'))
cc_dat$sp_stock<- paste(cc_dat$species,cc_dat$stock,sep='_')
length(unique(cc_dat$sp_stock))
cc_sp<- cc_dat %>% group_by(species) %>% summarize(n_distinct(stock))

#Find overlap with other datasets
chinook_info$sp_stock<- paste(chinook_info$species,chinook_info$stock,sep='_')
chum_info$sp_stock<- paste(chum_info$species,chum_info$stock,sep='_')
coho_info$sp_stock<- paste(coho_info$species,coho_info$stock,sep='_')
pink_info$sp_stock<- paste(pink_info$species,pink_info$stock,sep='_')
sockeye_info$sp_stock<- paste(sockeye_info$species,sockeye_info$stock,sep='_')

cc_chin<- subset(cc_dat,species=='Chinook')
cc_stocks_chi<- unique(cc_chin$sp_stock) #75
cc_stocks_chi[cc_stocks_chi %in% chinook_info$sp_stock] #19 of 20

cc_chum<- subset(cc_dat,species=='Chum')
cc_stocks_chu<- unique(cc_chum$sp_stock) #53
cc_stocks_chu[cc_stocks_chu %in% chum_info$sp_stock] #53 of 53

cc_coho<- subset(cc_dat,species=='Coho')
cc_stocks_coh<- unique(cc_coho$sp_stock) #6
cc_stocks_coh[cc_stocks_coh %in% coho_info$sp_stock] #0 shared

cc_pink<- subset(cc_dat,species=='Pink')
cc_stocks_pin<- unique(cc_pink$sp_stock) #46
cc_stocks_pin[cc_stocks_pin %in% pink_info$sp_stock] #46 shared

cc_sock<- subset(cc_dat,species=='Sockeye')
cc_stocks_sock<- unique(cc_sock$sp_stock) #48
cc_stocks_sock[cc_stocks_sock %in% sockeye_info$sp_stock] #41 shared
setdiff(cc_stocks_sock,sockeye_info$sp_stock) #checked these - all just name changes and are in original dataset
#early stuart, late stuart, early upper station, late upper station in original dataset; alagnak suggested not to be used


#Check and process these new time-series
cc_chin_info<- data.frame(index=NA,species=NA,location=NA,ts.length=NA,ts.start=NA,ts.end=NA)
for(i in 1:length(cc_stocks_chi)){
  c<- subset(cc_chin,stock==unique(cc_chin$stock)[i])
  c$logR_S<- log(c$rec/c$spawn)
  c_sr<- c[complete.cases(c$logR_S),]
  cc_chin_info[i,1]=i
  cc_chin_info[i,2]=unique(c$species)
  cc_chin_info[i,3]=unique(c$stock)
  cc_chin_info[i,4]=nrow(c_sr)
  cc_chin_info[i,5]=min(c_sr$broodYr)
  cc_chin_info[i,6]=max(c_sr$broodYr)
}
cc_chin_filtered<- subset(cc_chin_info,ts.length>=20) #64 chinook stocks
cc_chin_filtered<- subset(cc_chin_filtered, location %notin% chinook_info$stock) #55 stocks
write.csv(cc_chin_filtered,here('data','filtered datasets','cc_chinook_filtered.csv'))

cc_coho_info<- data.frame(index=NA,species=NA,location=NA,ts.length=NA,ts.start=NA,ts.end=NA)
for(i in 1:length(cc_stocks_coh)){
  c<- subset(cc_coho,stock==unique(cc_coho$stock)[i])
  c$logR_S<- log(c$rec/c$spawn)
  c_sr<- c[complete.cases(c$logR_S),]
  cc_coho_info[i,1]=i
  cc_coho_info[i,2]=unique(c$species)
  cc_coho_info[i,3]=unique(c$stock)
  cc_coho_info[i,4]=nrow(c_sr)
  cc_coho_info[i,5]=min(c_sr$broodYr)
  cc_coho_info[i,6]=max(c_sr$broodYr)
}
cc_coho_filtered<- subset(cc_coho_info,ts.length>=20) #6 coho stocks

write.csv(cc_coho_filtered,here('data','filtered datasets','cc_coho_filtered.csv'))

#Athena Ogden datasets
ao_dat<- read.csv(here('data','AthenaOgden datasets','AthenaOgden_trimmed_RS_timeseries.csv'))
unique(ao_dat$CU.or.PFMA.or.aggregate)
pse_overlap= c(unique(ao_dat$CU.or.PFMA.or.aggregate)[c(1:7,12:13)])
ao_dat_f<- subset(ao_dat, CU.or.PFMA.or.aggregate %notin% pse_overlap)
unique(ao_dat_f$CU.or.PFMA.or.aggregate) #5

ao_dat_info<- data.frame(index=NA,species=NA,location=NA,ts.length=NA,ts.start=NA,ts.end=NA,m.spawners=NA,sd.spawners=NA)
for(i in 1:length(unique(ao_dat_f$CU.or.PFMA.or.aggregate))){
  c<- subset(ao_dat_f,CU.or.PFMA.or.aggregate==unique(ao_dat_f$CU.or.PFMA.or.aggregate)[i])
  if(is.na(c$Spawner)==T){
    c$logR_S<- log(c$NTotRec/c$EFS)
  }else{
    c$logR_S<- log(c$NTotRec/c$Spawner)
  }
  c_sr<- c[complete.cases(c$logR_S),]
  
  ao_dat_info[i,1]=i
  ao_dat_info[i,2]=unique(c$Species)
  ao_dat_info[i,3]=unique(c$CU.or.PFMA.or.aggregate)
  ao_dat_info[i,4]=nrow(c_sr)
  ao_dat_info[i,5]=min(c_sr$Year)
  ao_dat_info[i,6]=max(c_sr$Year)
  if(is.na(c_sr$Spawner[1])==T){
    ao_dat_info[i,7]=mean(na.omit(c_sr$EFS))
    ao_dat_info[i,8]=sd(na.omit(c_sr$EFS))
  }else{
    ao_dat_info[i,7]=mean(na.omit(c_sr$Spawner))
    ao_dat_info[i,8]=sd(na.omit(c_sr$Spawner))
  }
  }
write.csv(ao_dat_info,here('data','filtered datasets','ao_dat_filtered.csv'))


##Maps

#Combined coordinate data
#Lat lons for each species
ao_lat<- read.csv(here('data','filtered datasets','ao_dat_filtered.csv')) %>% select(species,lat,lon)
pse_lat<- read.csv(here('data','filtered datasets','PSE_stockrecruit_info_latlons.csv')) %>% select(species,lat,lon)
for(i in 1:nrow(pse_lat)){if(pse_lat$species[i]=='Lake sockeye'|pse_lat$species[i]=='River sockeye'){pse_lat$species[i]='Sockeye'}}
for(i in 1:nrow(pse_lat)){if(pse_lat$species[i]=='Pink (even)'|pse_lat$species[i]=='Pink (odd)'){pse_lat$species[i]='Pink'}}
ch_lat<- read.csv(here('data','filtered datasets','filtered_chinook_info.csv')) %>% select(species,lat,lon)
co_lat<- read.csv(here('data','filtered datasets','filtered_coho_info.csv')) %>% select(species,lat,lon)
cu_lat<- read.csv(here('data','filtered datasets','filtered_chum_info.csv')) %>% select(species,lat,lon)
pi_lat<- read.csv(here('data','filtered datasets','filtered_pink_info.csv')) %>% select(species,lat,lon)
so_lat<- read.csv(here('data','filtered datasets','filtered_sockeye_info.csv')) %>% select(species,lat,lon)

sockeye_lat_lon<- rbind(ao_lat,so_lat,subset(pse_lat,species=='Sockeye'))
chinook_lat_lon<- rbind(ch_lat,subset(pse_lat,species=='Chinook'))
chum_lat_lon<- rbind(cu_lat,subset(pse_lat,species=='Chum'))
coho_lat_lon<- rbind(co_lat,subset(pse_lat,species=='Coho'))
pink_lat_lon<- rbind(pi_lat,subset(pse_lat,species=='Pink'))

all_lat_lon<- rbind(sockeye_lat_lon,chinook_lat_lon,chum_lat_lon,coho_lat_lon,pink_lat_lon)
length(unique(all_lat_lon[,2])); length(unique(all_lat_lon[,3])) #124 unique sites

##
world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("Longitude") + ylab("Latitude") +
  geom_point(data = subset(all_lat_lon,species=='Pink'), mapping = aes(x = lon, y = lat), color = "darksalmon", size = 3, alpha = 0.7) +
  geom_point(data = subset(all_lat_lon,species=='Chum'), mapping = aes(x = lon, y = lat), color = "darkgreen", size = 3, alpha = 0.7) +
  geom_point(data =  subset(all_lat_lon,species=='Sockeye'), mapping = aes(x = lon, y = lat), color = "darkred", size = 3, alpha = 0.7) +
  geom_point(data =  subset(all_lat_lon,species=='Coho'), mapping = aes(x = lon, y = lat), color = "darkblue", size = 3, alpha = 0.7) +
  geom_point(data = subset(all_lat_lon,species=='Chinook'), mapping = aes(x = lon, y = lat), color = "darkgray", size = 3, alpha = 0.7) +
  #  geom_text(data= world_points,aes(x=X, y=Y, label=name), color = 'darkblue', fontface = 'bold', check_overlap = FALSE) + 
  #  annotate(geom = “text”, x = -90, y = 26, label = “Gulf of Mexico”, fontface = “italic”, color = “grey22”, size = 6) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.75, 'in'), pad_y = unit(0.5, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -120), ylim = c(47, 72), expand = T) +
theme(panel.grid.major = element_line(color = gray(.5), linetype = 'dashed', size = 0.5), panel.background = element_rect(fill = 'lightblue'),legend.title = element_blank())
  
 # +                               # Add text element to plot
 # annotate("text", x = -155, y = 50, label = "Sockeye",color='darkred')

  #xlab(“Longitude”) + 
  #  ylab(“Latitude”) + 
  + # omit plot title saying 'color'
    scale_fill_manual(values = c('darkgray', 'darkgreen', 'darkblue','salmon','darkred'),
                      labels = c('Chinook','Chum','Coho','Pink','Sockeye'))

map + theme(legend.position = "bottom")
##To do on this:
#Remove junky time-series (ie. filter data FIRST)
#make it easier to see points with multiple species?
#add legend
#add salmon pics?

library(hrbrthemes)
filtered_all_dat<- read.csv(here('data','filtered datasets','overview_combined_mar10.csv'))
for(i in 1:nrow(filtered_all_dat)){if(filtered_all_dat$species[i]=='Pink (even)'|filtered_all_dat$species[i]=='Pink (odd)'){filtered_all_dat$species[i]='Pink'}}
sp_sum<- filtered_all_dat %>% group_by(species,region) %>% summarize(n=n())
ggplot(sp_sum, aes(fill=species, y=n, x=region))+
scale_fill_manual(values=c("darkgray", "darkgreen", "darkblue","darksalmon","darkred"))+ 
  geom_bar(position="stack", stat="identity") + theme_minimal() +
xlab('')+ylab('No. time-series')
