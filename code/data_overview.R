rm(list=ls())
##Salmon spawner-recruit time-series; dan.greenberg@dfo-mpo.gc.ca
library(here);library(dplyr);library(ggplot2);library(ggthemes);library(ggspatial);library(sf);library(rnaturalearth);library(rnaturalearthdata);library(wesanderson)
library(nlme); library(rlist)
#plot functions
source(here('code','functions.R'))

#Read datasets####
#sockeye
sockeye<- read.csv(here('data','raw data','sockeye','sockeye_data.csv'));sockeye_info<- read.csv(here('data','raw data','sockeye','sockeye_info.csv'));sockeye_source<- read.csv(here('data','raw data','sockeye','sockeye_sources.csv'))
psc_fraser_sockeye<- read.csv(here('data','raw data','sockeye','PSC_Fraser_broodtables.csv'))
skeena_nass_sockeye<- read.csv(here('data','raw data','sockeye','Skeena_Nass_sockeye.csv'))
ogden_comp<-  read.csv(here('data','raw data','multispecies','Salmon_RS_Database.csv')); ogden_info<-read.csv(here('data','raw data','multispecies','Salmon_RS_Time_Series_Summary.csv'))
pse_comp<- read.csv(here('data','raw data','multispecies','PSE_RS.csv'));pse_dq=read.csv(here('data','raw data','multispecies','PSE_data_quality.csv'))
#chum
chum<- read.csv(here('data','raw data','chum','chum_data.csv'));chum_info<- read.csv(here('data','raw data','chum','chum_info.csv'));chum_source<- read.csv(here('data','raw data','chum','chum_sources.csv'))
#pink
pink<- read.csv(here('data','raw data','pink','pink_data.csv'));pink_info<- read.csv(here('data','raw data','pink','pink_info.csv'));pink_source<- read.csv(here('data','raw data','pink','pink_sources.csv'))
#coho
coho<- read.csv(here('data','raw data','coho','coho_data_fwmar.csv'));coho_info<- read.csv(here('data','raw data','coho','coho_info.csv'));coho_source<- read.csv(here('data','raw data','coho','coho_sources.csv'))
ifr_coho<- read.csv(here('data','raw data','coho','SR_IFC_BY_98-16_created_2021-07-19.csv'))
#chinook
chinook<- read.csv(here('data','raw data','chinook','chinook_data_totalage.csv'));chinook_info<- read.csv(here('data','raw data','chinook','chinook_info.csv'));chinook_source<- read.csv(here('data','raw data','chinook','chinook_sources.csv')) 
cc_comp<- read.csv(here('data','raw data','multispecies','AK-WCoast-Salmon-SR.csv'))
cow_chin<- read.csv(here('data','raw data','chinook','Cowichan_chinook_broodtable.csv'),na.strings = c('#N/A'))
harrison_chin<-read.csv(here('data','raw data','chinook','Harrison_chinook_broodtable.csv'))
shuswap_chin<-read.csv(here('data','raw data','chinook','Lower_Shuswap_chinook_broodtable.csv'))
nicola_chin<-read.csv(here('data','raw data','chinook','Nicola_chinook_broodtable.csv'))

#(partially) synonymize column names for sanity
colnames(sockeye)<- tolower(names(sockeye));colnames(sockeye_info)<- tolower(names(sockeye_info));colnames(sockeye_source)<- tolower(names(sockeye_source))
names(chum)<- gsub('.yr','year',names(chum)) #make the brood years equivalent
names(pink)<- gsub('.yr','year',names(pink))  #make the brood years equivalent
colnames(coho)<- tolower(names(coho));colnames(coho_info)<- tolower(names(coho_info));colnames(coho_source)<- tolower(names(coho_source));names(coho_source)[1]='source.id'
colnames(chinook)<- tolower(names(chinook));colnames(chinook_info)<- tolower(names(chinook_info));colnames(chinook_source)<- tolower(names(chinook_source))
colnames(ogden_comp)<- tolower(names(ogden_comp));colnames(ogden_info)<- tolower(names(ogden_info))
names(cc_comp)<- gsub('broodYr','broodyear',names(cc_comp));names(cc_comp)<- gsub('spawn','spawners',names(cc_comp));names(cc_comp)<- gsub('rec','recruits',names(cc_comp));colnames(cc_comp)[1]='stock.id'  #make the brood years equivalent
names(cow_chin)[1]='stock.id'
names(psc_fraser_sockeye)[1]='stock';names(psc_fraser_sockeye)[3]='broodyear'
names(skeena_nass_sockeye)=tolower(names(skeena_nass_sockeye));names(skeena_nass_sockeye)[5]='broodyear';names(skeena_nass_sockeye)[6]='spawners';names(skeena_nass_sockeye)[10]='recruits'
names(ifr_coho)[1]='stock';names(ifr_coho)[4]='spawners';names(ifr_coho)[5]='recruits'
names(harrison_chin)[1]='broodyear';names(harrison_chin)[4]='spawners';names(harrison_chin)[5]='recruits'
names(shuswap_chin)[1]='broodyear';names(shuswap_chin)[4]='spawners';names(shuswap_chin)[5]='recruits'
names(nicola_chin)[3]='broodyear';names(shuswap_chin)[4]='spawners';names(shuswap_chin)[5]='recruits'

#add ocean regions to info
chum_info$ocean.region=sockeye_info$ocean.region[match(chum_info$region,sockeye_info$region)]
chum_info$ocean.region[1:9]='WC'
pink_info$ocean.region=chum_info$ocean.region[match(pink_info$region,chum_info$region)]

#Sockeye####
stock_dat=data.frame(stock.id=NA,species=NA,stock.name=NA,lat=NA,lon=NA,ocean.basin=NA,state=NA,begin=NA,end=NA,n.years=NA,max.spawners=NA,max.recruits=NA,source=NA,url=NA,comments=NA)
#General stock characteristics, start and end of time-series, number of useable years (removing years with useflags),average spawner and recruit 

#From the top - sockeye compilation

#Merge in new PSC data for Fraser stocks to existing compilation
fr_sock<- subset(sockeye, stock %in% psc_fraser_sockeye$stock)

sockeye2<- subset(sockeye, stock %notin% psc_fraser_sockeye$stock) #Drop out older data for Fraser R stocks

#Process the remaining stock data
sockeye_list=list()
for(i in 1:length(unique(sockeye2$stock.id))){
  s=subset(sockeye,stock.id==unique(sockeye2$stock.id)[i])
  s_info<- subset(sockeye_info,stock.id==unique(sockeye2$stock.id)[i])
  s_use=subset(s,useflag==1) %>% subset(is.na(spawners)==F&is.na(recruits)==F)
  
  stock_dat[i,1]=unique(s$stock.id)
  stock_dat[i,2]=unique(s$species)
  stock_dat[i,3]=paste(unique(s$stock),unique(s$species),sep='-')
  stock_dat[i,4]=unique(s_info$lat)
  stock_dat[i,5]=unique(s_info$lon)
  stock_dat[i,6]=unique(s_info$ocean.region)
  stock_dat[i,7]=unique(s_info$jurisdiction)
  
  if(nrow(s_use)!=0){
    stock_dat[i,8]=min(s_use$broodyear)
    stock_dat[i,9]=max(s_use$broodyear)
    stock_dat[i,10]=length(s_use$broodyear)
    stock_dat[i,11]=max(s_use$spawners)
    stock_dat[i,12]=max(s_use$recruits)
  }else{
    stock_dat[i,8]=NA
    stock_dat[i,9]=NA
    stock_dat[i,10]=0
    stock_dat[i,11:12]=NA
  }
  
  stock_dat[i,13]=sockeye_source$source[match(s_info$source.id,sockeye_source$source.id)]
  stock_dat[i,15]=s_info$comment..we.will.update.this.later.
  
  sockeye_list[[i]]=s_use[,c('stock','species','broodyear','recruits','spawners')]
}

#Fraser sockeye stocks - PSC 2022 production dataset 
for(i in 1:length(unique(psc_fraser_sockeye$stock))){
  s=subset(psc_fraser_sockeye,production_stock_code==unique(psc_fraser_sockeye$production_stock_code)[i])
  s= s %>% mutate(recruits = rowSums(s[,9:21],na.rm = TRUE))
  s$species<- rep('Sockeye',nrow(s))
  #to determine whether to use effective female spawners or total spawners based on data availability
  if(length(na.omit(s$total_broodyr_EFS))<length(na.omit(s$total_broodyr_spawners))){
   names(s)[4]='spawners' 
  }
  names(s)[5]='spawners' #use effective female spawners for these stocks
  
  stock_dat_temp=data.frame(stock.id=NA,species=NA,stock.name=NA,lat=NA,lon=NA,ocean.basin=NA,state=NA,begin=NA,end=NA,n.years=NA,max.spawners=NA,max.recruits=NA,source=NA,url=NA,comments=NA)

  stock_dat_temp[,1]=NA
  stock_dat_temp[,2]='Sockeye'
  stock_dat_temp[,3]=paste(unique(s$stock),'Sockeye',sep='-')
  stock_dat_temp[,4]=sockeye_info$lat[2] #lat for mouth of Fraser
  stock_dat_temp[,5]=sockeye_info$lon[2] #lon for mouth of Fraser
  stock_dat_temp[,6]='WC' #West Coast
  stock_dat_temp[,7]='BC' #British Columbia
  
  stock_dat_temp[,8]=min(s$broodyear)
  stock_dat_temp[,9]=max(s$broodyear)
  stock_dat_temp[,10]=length(s$broodyear)
  stock_dat_temp[,11]=max(s$spawners)
  stock_dat_temp[,12]=max(s$recruits)
  stock_dat_temp[,13]='Eric Taylor, Pacific Salmon Commission, 2022'
  stock_dat_temp[,15]=NA
  
  stock_dat=rbind(stock_dat,stock_dat_temp)
  
  sockeye_list[[nrow(stock_dat)]]=s[,c('stock','species','broodyear','recruits','spawners')]
}

#Skeena-Nass sockeye
skeena_nass_sockeye_4k=subset(skeena_nass_sockeye,label=='Filter45') #Use data without infill, but filtered out productivity (R/S) estimates over 45
skeena_nass_sockeye_4k=skeena_nass_sockeye_4k[complete.cases(skeena_nass_sockeye_4k$stock),]
skeena_nass_sockeye_4k=subset(skeena_nass_sockeye_4k,is.na(stock)==F)
skeena_nass_sockeye_4k=subset(skeena_nass_sockeye_4k,stock!='Meziadin') #longer time-series in other dataset, so removing this
skeena_nass_sockeye_4k$stock=gsub('Bear','Bear-Skeena',skeena_nass_sockeye_4k$stock) #renaming Bear stock to avoid conflict with Bear R. Alaska
for(i in 1:length(unique(skeena_nass_sockeye_4k$stock))){
  s=subset(skeena_nass_sockeye_4k,stock==unique(skeena_nass_sockeye_4k$stock)[i])
  s= s[complete.cases(s$spawners),];s=s[complete.cases(s$recruits),] #keep years with both spawner & recruit estimates
  if(nrow(s)==0){next}
  s$species<- rep('Sockeye',nrow(s))
 
  stock_dat_temp=data.frame(stock.id=NA,species=NA,stock.name=NA,lat=NA,lon=NA,ocean.basin=NA,state=NA,begin=NA,end=NA,n.years=NA,max.spawners=NA,max.recruits=NA,source=NA,url=NA,comments=NA)
  
  stock_dat_temp[,1]=unique(s$stock)
  stock_dat_temp[,2]='Sockeye'
  stock_dat_temp[,3]=paste(unique(s$stock),'Sockeye',sep='-')
  if(s$stock[1]=='Lower Nass Sea & River Type'|s$stock[1]=='Meziadin'|s$stock[1]=='Kwinageese'|s$stock[1]=='Damdochax'|s$stock[1]=='Upper Nass River Type'){ #Nass stocks
    stock_dat_temp[,4]=54.9898 #lat need to do these (one for Skeena estuary and one for Nass)
    stock_dat_temp[,5]=-130.02 #lon 
    
  }else{
    stock_dat_temp[,4]=54.2237 #lat need to do these (one for Skeena estuary and one for Nass)
    stock_dat_temp[,5]=-129.831 #lon 
  }
  stock_dat_temp[,6]='WC'
  stock_dat_temp[,7]='BC'
  
  stock_dat_temp[,8]=min(s$broodyear)
  stock_dat_temp[,9]=max(s$broodyear)
  stock_dat_temp[,10]=length(s$broodyear)
  stock_dat_temp[,11]=max(s$spawners)
  stock_dat_temp[,12]=max(s$recruits)
  stock_dat_temp[,13]='Charmaine Carr-Harris, DFO, 2022'
  stock_dat_temp[,15]='Used filtered data - R/S >45 were removed'
  
  stock_dat=rbind(stock_dat,stock_dat_temp)
  
  sockeye_list[[nrow(stock_dat)]]=s[,c('stock','species','broodyear','recruits','spawners')]
}

sockeye_filtered<- do.call("rbind", sockeye_list)

#Chum####
chum_list=list()
for(i in 1:length(unique(chum$stock.id))){
  s=subset(chum,stock.id==unique(chum$stock.id)[i])
  s_info<- subset(chum_info,stock.id==unique(chum$stock.id)[i])
  s_use=subset(s,use==1) %>% subset(is.na(spawners)==F&is.na(recruits)==F)
  
  stock_dat_temp=data.frame(stock.id=NA,species=NA,stock.name=NA,lat=NA,lon=NA,ocean.basin=NA,state=NA,begin=NA,end=NA,n.years=NA,max.spawners=NA,max.recruits=NA,source=NA,url=NA,comments=NA)
  
  stock_dat_temp[,1]=unique(s$stock.id)
  stock_dat_temp[,2]=unique(s$species)
  stock_dat_temp[,3]=paste(unique(s$stock),unique(s$species),sep='-')
  stock_dat_temp[,4]=unique(s_info$lat)
  stock_dat_temp[,5]=unique(s_info$lon)
  stock_dat_temp[,6]=unique(s_info$ocean.region)
  stock_dat_temp[,7]=unique(s_info$jurisdiction)
  
  if(nrow(s_use)!=0){
    stock_dat_temp[,8]=min(s_use$broodyear)
    stock_dat_temp[,9]=max(s_use$broodyear)
    stock_dat_temp[,10]=length(s_use$broodyear)
    stock_dat_temp[,11]=max(s_use$spawners)
    stock_dat_temp[,12]=max(s_use$recruits)
    }else{
    stock_dat_temp[,8]=NA
    stock_dat_temp[,9]=NA
    stock_dat_temp[,10]=0
    stock_dat_temp[,11:12]=NA
  }
 
  s.id=as.numeric(strsplit(s_info$source.id,',')[[1]])
  if(length(s.id)==1){ stock_dat_temp[,13]=chum_source$source[match(s_info$source.id,sockeye_source$source.id)]
}
  if(length(s.id)==2){
    source<- subset(chum_source, source.id %in% s.id)
    stock_dat_temp[,13]=paste(source$source[1],source$source[2],sep='; ')
  }
  stock_dat_temp[,15]=s_info$comment
  
  stock_dat=rbind(stock_dat,stock_dat_temp)
  
  chum_list[[i]]=s_use[,c('stock','species','broodyear','recruits','spawners')]
}
chum_filtered<- do.call("rbind", chum_list)

#Pink####
pink_list=list()
for(i in 1:length(unique(pink$stock.id))){
  s=subset(pink,stock.id==unique(pink$stock.id)[i])
  s_info<- subset(pink_info,stock.id==unique(pink$stock.id)[i])
  s_use=subset(s,use==1) %>% subset(is.na(spawners)==F&is.na(recruits)==F)
  
  stock_dat_temp=data.frame(stock.id=NA,species=NA,stock.name=NA,lat=NA,lon=NA,ocean.basin=NA,state=NA,begin=NA,end=NA,n.years=NA,max.spawners=NA,max.recruits=NA,source=NA,url=NA,comments=NA)
  
  stock_dat_temp[,1]=unique(s$stock.id)
  stock_dat_temp[,2]=unique(s$species)
  stock_dat_temp[,3]=paste(unique(s$stock),unique(s$species),sep='-')
  stock_dat_temp[,4]=unique(s_info$lat)
  stock_dat_temp[,5]=unique(s_info$lon)
  stock_dat_temp[,6]=unique(s_info$ocean.region)
  stock_dat_temp[,7]=unique(s_info$jurisdiction)
  
  if(nrow(s_use)!=0){
    stock_dat_temp[,8]=min(s_use$broodyear)
    stock_dat_temp[,9]=max(s_use$broodyear)
    stock_dat_temp[,10]=length(s_use$broodyear)
    stock_dat_temp[,11]=max(s_use$spawners)
    stock_dat_temp[,12]=max(s_use$recruits)
  }else{
    stock_dat_temp[,8]=NA
    stock_dat_temp[,9]=NA
    stock_dat_temp[,10]=0
    stock_dat_temp[,11:12]=NA
  }
  
  s.id=as.numeric(strsplit(s_info$source.id,',')[[1]])
  if(length(s.id)==1){stock_dat_temp[,13]=pink_source$source[match(s_info$source.id,sockeye_source$source.id)]
  }
  if(length(s.id)==2){
    source<- subset(pink_source, source.id %in% s.id)
    stock_dat_temp[,13]=paste(source$source[1],source$source[2],sep='; ')
  }
  stock_dat_temp[,15]=s_info$comment
  
  stock_dat=rbind(stock_dat,stock_dat_temp)
  
  pink_list[[i]]=s_use[,c('stock','species','broodyear','recruits','spawners')]
}
pink_filtered<- do.call("rbind", pink_list)

#Chinook####
#Add in Cowichan chinook
chinook=rbind(chinook,cow_chin) #add in S-R data
chinook_info[21,1:7]=c(333,'Chinook','Cowichan','WC','Vancouver Island','Vancouver Island','BC');chinook_info$lat[21]=48.7581;chinook_info$lon[21]=-123.6242;chinook_info$source.id[21]=10 #add in metadata

chinook_info[22,1:7]=c(334,'Chinook','Harrison','WC','Fraser River','Vancouver Island','BC');chinook_info$lat[22]=sockeye_info$lat[2];chinook_info$lon[22]=sockeye_info$lon[2];chinook_info$source.id[22]=11 #add in metadata
chinook_info[23,1:7]=c(335,'Chinook','Lower Shuswap','WC','Fraser River','Fraser River','BC');chinook_info$lat[23]=sockeye_info$lat[2];chinook_info$lon[23]=sockeye_info$lon[2];chinook_info$source.id[23]=11 #add in metadata
chinook_info[24,1:7]=c(336,'Chinook','Nicola','WC','Fraser River','Fraser River','BC');chinook_info$lat[24]=sockeye_info$lat[2];chinook_info$lon[24]=sockeye_info$lon[2];chinook_info$source.id[24]=12 #add in metadata
chinook_source[10,1]=10;chinook_source$title[10]='Karalea Cantera, DFO, 2022' #add in source
chinook_source[11,1]=11;chinook_source$title[11]='Chucken Parken, DFO, 2022' #add in source
chinook_source[12,1]=12;chinook_source$title[12]='Luke Warkentin, DFO, 2022' #add in source

chinook_list=list()
for(i in 1:length(unique(chinook$stock.id))){
  s=subset(chinook,stock.id==unique(chinook$stock.id)[i])
  s_info<- subset(chinook_info,stock.id==unique(chinook$stock.id)[i])
  s_use=subset(s,useflag==1) %>% subset(is.na(spawners)==F&is.na(recruits)==F)
  
  stock_dat_temp=data.frame(stock.id=NA,species=NA,stock.name=NA,lat=NA,lon=NA,ocean.basin=NA,state=NA,begin=NA,end=NA,n.years=NA,max.spawners=NA,max.recruits=NA,source=NA,url=NA,comments=NA)
  
  stock_dat_temp[,1]=unique(s$stock.id)
  stock_dat_temp[,2]=unique(s$species)
  stock_dat_temp[,3]=paste(unique(s$stock),unique(s$species),sep='-')
  stock_dat_temp[,4]=unique(s_info$lat)
  stock_dat_temp[,5]=unique(s_info$lon)
  stock_dat_temp[,6]=unique(s_info$ocean.region)
  stock_dat_temp[,7]=unique(s_info$jurisdiction)
  
  if(nrow(s_use)!=0){
    stock_dat_temp[,8]=min(s_use$broodyear)
    stock_dat_temp[,9]=max(s_use$broodyear)
    stock_dat_temp[,10]=length(s_use$broodyear)
    stock_dat_temp[,11]=max(s_use$spawners)
    stock_dat_temp[,12]=max(s_use$recruits)
  }else{
    stock_dat_temp[,8]=NA
    stock_dat_temp[,9]=NA
    stock_dat_temp[,10]=0
    stock_dat_temp[,11:12]=NA
  }
  
  
  s.id=as.numeric(strsplit(as.character(s_info$source.id),',')[[1]])
  if(length(s.id)==1){
    stock_dat_temp[,13]=chinook_source$title[match(s_info$source.id,sockeye_source$source.id)]
    stock_dat_temp[,14]=chinook_source$url[match(s_info$source.id,sockeye_source$source.id)]
    
  }
  if(length(s.id)==2){
    source<- subset(chinook_source, source.id %in% s.id)
    stock_dat_temp[,13]=paste(source$title[1],source$title[2],sep='; ')
    stock_dat_temp[,14]=paste(source$url[1],source$url[2],sep='; ')
  }
  stock_dat_temp[,15]=NA #no comments
  
  stock_dat=rbind(stock_dat,stock_dat_temp)
  
  chinook_list[[i]]=s_use[,c('stock','species','broodyear','recruits','spawners')]
}

#Harrison chinook
harrison_chin3=subset(harrison_chin,is.na(recruits)==F&is.na(spawners)==F) #remove years without recruit data

row.n=nrow(stock_dat)+1
stock_dat[row.n,1]=NA #no pre-assigned stock id
stock_dat[row.n,2]=chinook_info[22,2]
stock_dat[row.n,3]=paste(chinook_info[22,3],chinook_info[22,2],sep='-')
stock_dat[row.n,4]=chinook_info[22,8]
stock_dat[row.n,5]=chinook_info[22,9]
stock_dat[row.n,6]='WC'
stock_dat[row.n,7]='BC'
stock_dat[row.n,8]=min(harrison_chin3$broodyear)
stock_dat[row.n,9]=max(harrison_chin3$broodyear)
stock_dat[row.n,10]=length(harrison_chin3$broodyear)
stock_dat[row.n,11]=max(harrison_chin3$spawners)
stock_dat[row.n,12]=max(harrison_chin3$recruits)
stock_dat[row.n,13]=chinook_source$title[chinook_info$source.id[22]]
stock_dat[row.n,14]='NA'
stock_dat[row.n,15]='NA'

harrison_chin3$stock=rep('Harrison',nrow(harrison_chin3))
harrison_chin3$species=rep('Chinook',nrow(harrison_chin3))

chinook_list[[22]]=harrison_chin3[,c('stock','species','broodyear','recruits','spawners')]

#Lower Shuswap chinook
shuswap_chin$stock=rep('Lower Shuswap',nrow(shuswap_chin))
shuswap_chin$species=rep('Chinook',nrow(shuswap_chin))
shuswap_chin=subset(shuswap_chin,is.na(recruits)==F&is.na(spawners)==F) #remove years without recruit data

row.n=nrow(stock_dat)+1
stock_dat[row.n,1]=NA #no pre-assigned stock id
stock_dat[row.n,2]=chinook_info[23,2]
stock_dat[row.n,3]=paste(chinook_info[23,3],chinook_info[23,2],sep='-')
stock_dat[row.n,4]=chinook_info[23,8]
stock_dat[row.n,5]=chinook_info[23,9]
stock_dat[row.n,6]='WC'
stock_dat[row.n,7]='BC'
stock_dat[row.n,8]=min(shuswap_chin$broodyear)
stock_dat[row.n,9]=max(shuswap_chin$broodyear)
stock_dat[row.n,10]=length(shuswap_chin$broodyear)
stock_dat[row.n,11]=max(shuswap_chin$spawners)
stock_dat[row.n,12]=max(shuswap_chin$recruits)
stock_dat[row.n,13]=chinook_source$title[chinook_info$source.id[23]]

shuswap_chin$stock=rep('Lower Shuswap',nrow(shuswap_chin))
shuswap_chin$species=rep('Chinook',nrow(shuswap_chin))
shuswap_chin=subset(shuswap_chin,is.na(recruits)==F&is.na(spawners)==F) #remove years without recruit data


chinook_list[[23]]=shuswap_chin[,c('stock','species','broodyear','recruits','spawners')]

#Nicola chinook
#use total spawners as per Warkentin et al. 2020
nicola_recruits=nicola_chin %>% group_by(broodyear) %>% summarize(recruits=sum(recruits),spawners=sum(total_spawners))
#follow Warkentin et al. 2020 and use cohorts from 1992 to 2013
nicola_recruits=subset(nicola_recruits,broodyear>1991)
nicola_recruits$stock=rep('Nicola',nrow(nicola_recruits))
nicola_recruits$species=rep('Chinook',nrow(nicola_recruits))

row.n=nrow(stock_dat)+1
stock_dat[row.n,1]=NA #no pre-assigned stock id
stock_dat[row.n,2]=chinook_info[24,2]
stock_dat[row.n,3]=paste(chinook_info[24,3],chinook_info[24,2],sep='-')
stock_dat[row.n,4]=chinook_info[24,8]
stock_dat[row.n,5]=chinook_info[24,9]
stock_dat[row.n,6]='WC'
stock_dat[row.n,7]='BC'
stock_dat[row.n,8]=min(nicola_recruits$broodyear)
stock_dat[row.n,9]=max(nicola_recruits$broodyear)
stock_dat[row.n,10]=length(nicola_recruits$broodyear)
stock_dat[row.n,11]=max(nicola_recruits$spawners)
stock_dat[row.n,12]=max(nicola_recruits$recruits)
stock_dat[row.n,13]=chinook_source$title[chinook_info$source.id[24]]

chinook_list[[24]]=nicola_recruits[,c('stock','species','broodyear','recruits','spawners')]

chinook_filtered<- do.call("rbind", chinook_list)

#Coho####

coho_list=list()
for(i in 1:length(unique(coho$stock.id))){
  s=subset(coho,stock.id==unique(coho$stock.id)[i])
  #Some stocks from have repeated time-series with multiple estimates of escapement based extrapolations of true escapement from peak abundance estimated from aerial surveys. These are eg. 25%, 50%, 100% of true escapement.
  #Just taking a single time-series for these - the particular estimate does not matter for our purposes in comparing stationary and time-varying stock-recruitment dynamics
  if(any(duplicated(s$broodyear))==T){
    s=distinct(s,broodyear,.keep_all = T)
  }
  s_info<- subset(coho_info,stock.id==unique(coho$stock.id)[i])
  s_use=subset(s,useflag==1) %>% subset(is.na(spawners)==F&is.na(recruits)==F)
  
  stock_dat_temp=data.frame(stock.id=NA,species=NA,stock.name=NA,lat=NA,lon=NA,ocean.basin=NA,state=NA,begin=NA,end=NA,n.years=NA,max.spawners=NA,max.recruits=NA,source=NA,url=NA,comments=NA)
  
  stock_dat_temp[,1]=unique(s$stock.id)
  stock_dat_temp[,2]=unique(s$species)
  stock_dat_temp[,3]=paste(unique(s$stock),unique(s$species),sep='-')
  stock_dat_temp[,4]=unique(s_info$lat)
  stock_dat_temp[,5]=unique(s_info$lon)
  stock_dat_temp[,6]=unique(s_info$ocean.region)
  stock_dat_temp[,7]=unique(s_info$jurisdiction)
  
  if(nrow(s_use)!=0){
    stock_dat_temp[,8]=min(s_use$broodyear)
    stock_dat_temp[,9]=max(s_use$broodyear)
    stock_dat_temp[,10]=length(s_use$broodyear)
    stock_dat_temp[,11]=max(s_use$spawners)
    stock_dat_temp[,12]=max(s_use$recruits)
  }else{
    stock_dat_temp[,8]=NA
    stock_dat_temp[,9]=NA
    stock_dat_temp[,10]=0
    stock_dat_temp[,11:12]=NA
  }
  
  s.id=as.numeric(strsplit(as.character(s_info$source.id),',')[[1]])
  if(length(s.id)==1){
    stock_dat_temp[,13]=coho_source$title[match(s_info$source.id,sockeye_source$source.id)]
    stock_dat_temp[,14]=coho_source$url[match(s_info$source.id,sockeye_source$source.id)]
  }
  if(length(s.id)==2){
    source<- subset(coho_source, source.id %in% s.id)
    stock_dat_temp[,13]=paste(source$source[1],source$source[2],sep='; ')
    stock_dat_temp[,14]=paste(source$url[1],source$url[2],sep='; ')
  }
  stock_dat_temp[,15]=NA #no comments
  
  stock_dat=rbind(stock_dat,stock_dat_temp)
  
  coho_list[[i]]=s_use[,c('stock','species','broodyear','recruits','spawners')]
}

#Interior Fraser coho brood tables
for(i in 1:length(unique(ifr_coho$stock))){
  s=subset(ifr_coho,stock==unique(ifr_coho$stock)[i])
  s$species=rep('Coho',nrow(s))
  
  stock_dat_temp=data.frame(stock.id=NA,species=NA,stock.name=NA,lat=NA,lon=NA,ocean.basin=NA,state=NA,begin=NA,end=NA,n.years=NA,max.spawners=NA,max.recruits=NA,source=NA,url=NA,comments=NA)
  
  stock_dat_temp[,1]=unique(s$stock)
  stock_dat_temp[,2]='Coho'
  stock_dat_temp[,3]=paste(unique(s$stock),'Coho',sep='-')
  stock_dat_temp[,4]=sockeye_info$lat[2] #lat for mouth of Fraser
  stock_dat_temp[,5]=sockeye_info$lon[2] #lon for mouth of Fraser
  stock_dat_temp[,6]='WC'
  stock_dat_temp[,7]='BC'
  stock_dat_temp[,8]=min(s$broodyear)
  stock_dat_temp[,9]=max(s$broodyear)
  stock_dat_temp[,10]=length(s$broodyear)
  stock_dat_temp[,11]=max(s$spawners)
  stock_dat_temp[,12]=max(s$recruits)
  stock_dat_temp[,13]='Michael Arbeider, DFO, 2022'
  stock_dat_temp[,14]=NA #no comments

  stock_dat=rbind(stock_dat,stock_dat_temp)
  
  coho_list[[length(unique(coho$stock.id))+i]]=s[,c('stock','species','broodyear','recruits','spawners')]
}

coho_filtered<- do.call("rbind", coho_list)

##Issue 1: many redundant stocks across datasets - first pass to filter these out
###Curry Cunningham compilation
cc_comp$stock.name<- paste(cc_comp$stock,cc_comp$species,sep='-')
cc_comp2<- subset(cc_comp, stock.name %notin% stock_dat$stock.name)

unique(cc_comp2$stock.name)
#Sockeye stocks in here are repeats from the other dataset(early Stuart/Upper Station, etc.)
cc_comp2<- subset(cc_comp2,species!='Sockeye')

#Approximate lat/lons for these stocks taken from google maps - at mouth of the river system
alsek_klukshu=cbind(59.1348,-138.6067)
columbia=cbind(46.2508,-124.0143) #for all stocks in Columbia river system (Columbia, snake river, williamette)
oregon_coast=cbind(44.6175, -124.1093) #ESU for chinook south of Columbia/north of Cape Blanco - location from Newport (approx. mid range)
#puget sound stock lat/lons:
puget_sound_lat_lons<- data.frame(river.system=NA,lat=NA,lon=NA)
puget_sound_lat_lons[1,]=cbind('Green River',47.5925,-122.3600) #Green River lat/lon
puget_sound_lat_lons[2,]=cbind('Cedar River',47.5039,-122.2169) #Cedar River lat/lon
puget_sound_lat_lons[3,]=cbind('Skagit River',48.3677,-122.5141) #Skagit river & sauk/suiattle tributaries lat/lon
puget_sound_lat_lons[4,]=cbind('Stillaguamish River',48.2465,-122.3957) #Stillaguimish river
puget_sound_lat_lons[5,]=cbind('Nisqually River',47.0996,-122.7029) #Nisqually river
puget_sound_lat_lons[6,]=cbind('Puyallup River',47.2472,-122.4289) #Puyallup river & white river tributary
puget_sound_lat_lons[7,]=cbind('Skokomish River',47.3436,-123.1212) #Skokomish river & white river tributary
puget_sound_lat_lons[8,]=cbind('Snohomisin River',48.0206,-122.2122) #Tributaries of Snohomisin (Snoqualmie & Skykomish)
puget_sound_lat_lons[9,]=cbind('Mid-Hood Canal',47.6393,-122.9299) #mid hood canal - diswallups/duckabush/hamma hamma watersheds -coordinate is at mouth of the Duckabush system (mid way)

#putting in lat/lon info for these stocks manually:
cc_info<- distinct(cc_comp2,stock.name,.keep_all=T)
cc_info$lat=NA;cc_info$lon=NA
cc_info[1,12:13]=alsek_klukshu #Alsek-Chinook
cc_info[2:28,12]=columbia[,1];cc_info[2:28,13]=columbia[,2] #Columbia river chinook (26 pop)
cc_info[29,12:13]=puget_sound_lat_lons[2,2:3] #Cedar river chinook
cc_info[30,12:13]=puget_sound_lat_lons[1,2:3] #Green river
cc_info[31,12:13]=puget_sound_lat_lons[3,2:3] #Sauk river
cc_info[32,12:13]=puget_sound_lat_lons[3,2:3] #skagit river
cc_info[33,12:13]=puget_sound_lat_lons[9,2:3] #mid-hood canal
cc_info[34,12:13]=puget_sound_lat_lons[5,2:3] #Nisqually river
cc_info[35,12:13]=puget_sound_lat_lons[4,2:3] #Stillaguimish river
cc_info[36,12:13]=puget_sound_lat_lons[6,2:3] #Puyallup river
cc_info[37,12:13]=puget_sound_lat_lons[7,2:3] #Skokomish river
cc_info[38,12:13]=puget_sound_lat_lons[8,2:3] #Skykomish River
cc_info[39,12:13]=puget_sound_lat_lons[8,2:3] # Snoqualmie River
cc_info[40,12:13]=puget_sound_lat_lons[4,2:3] #Stillaguimish river
cc_info[41,12:13]=puget_sound_lat_lons[3,2:3] #Suiattle River
cc_info[42,12:13]=puget_sound_lat_lons[3,2:3] #Upper Sauk River
cc_info[43,12:13]=puget_sound_lat_lons[3,2:3] #Skagit River
cc_info[44,12:13]=puget_sound_lat_lons[6,2:3] #White River
cc_info[45:58,12]=columbia[,1];cc_info[45:58,13]=columbia[,2] #Willamette (Columbia)
cc_info[59:60,12]=oregon_coast[,1];cc_info[59:60,13]=oregon_coast[,2] #Oregon coast coho
cc_info[61:62,12]=columbia[,1];cc_info[61:62,13]=columbia[,2] #Willamette (Columbia) coho

#add in these new stocks
cc_comp_list=list()
for(i in 1:length(unique(cc_comp2$stock.name))){
  s=subset(cc_comp2,stock.name==unique(cc_comp2$stock.name)[i])
  s_info=subset(cc_info,stock.name==unique(cc_comp2$stock.name)[i])
  #Some stocks from have repeated time-series with multiple estimates of escapement based extrapolations of true escapement from peak abundance estimated from aerial surveys. These are eg. 25%, 50%, 100% of true escapement.
  #Just taking a single time-series for these - the particular estimate does not matter for our purposes in comparing stationary and time-varying stock-recruitment dynamics
  if(any(duplicated(s$broodyear))==T){
    s=distinct(s,broodyear,.keep_all = T)
  }
  s_use=subset(s,is.na(spawners)==F&is.na(recruits)==F)
  
  stock_dat_temp=data.frame(stock.id=NA,species=NA,stock.name=NA,lat=NA,lon=NA,ocean.basin=NA,state=NA,begin=NA,end=NA,n.years=NA,max.spawners=NA,max.recruits=NA,source=NA,url=NA,comments=NA)
  
  stock_dat_temp[,1]=unique(s$stock.id)
  stock_dat_temp[,2]=unique(s$species)
  stock_dat_temp[,3]=paste(unique(s$stock),unique(s$species),sep='-')
  stock_dat_temp[,4]=unique(s_info$lat)
  stock_dat_temp[,5]=unique(s_info$lon)
  stock_dat_temp[,6]=unique(s_info$large.region)
  if(s_info$region=='Southeast'){stock_dat_temp[,7]='AK'}
  if(s_info$region=='Interior Columbia'|s_info$region=='Puget Sound'){stock_dat_temp[,7]='WA'}
  if(s_info$region=='Willamette/Lower Columbia'|s_info$region=='Oregon Coast'){stock_dat_temp[,7]='OR'}
  if(nrow(s_use)!=0){
    stock_dat_temp[,8]=min(s_use$broodyear)
    stock_dat_temp[,9]=max(s_use$broodyear)
    stock_dat_temp[,10]=length(s_use$broodyear)
    stock_dat_temp[,11]=max(s_use$spawners)
    stock_dat_temp[,12]=max(s_use$recruits)
  }else{
    stock_dat_temp[,8]=NA
    stock_dat_temp[,9]=NA
    stock_dat_temp[,10]=0
    stock_dat_temp[,11:12]=NA
  }
  stock_dat_temp[,13]='Curry Cunningham/Brian Burke, 2022' #source to be confirmed with Curry
  
  stock_dat=rbind(stock_dat,stock_dat_temp)
  
  cc_comp_list[[i]]=s_use[,c('stock','species','broodyear','recruits','spawners')]
}
cc_comp_filtered<- do.call("rbind", cc_comp_list)


#Print out data####
filtered_productivity_data=rbind(sockeye_filtered,chum_filtered,pink_filtered,coho_filtered,chinook_filtered,cc_comp_filtered)

#Stock overview
stock_dat=subset(stock_dat,n.years>0)
stock_dat$stock.id=seq(1:nrow(stock_dat))
filtered_productivity_data$stock.id=stock_dat$stock.id[match(paste(filtered_productivity_data$stock,filtered_productivity_data$species,sep='-'),stock_dat$stock.name)]

#Write datasets
write.csv(filtered_productivity_data,here('data','filtered datasets','salmon_productivity_compilation_feb2023.csv'))
write.csv(stock_dat,here('data','filtered datasets','all_stocks_info_feb2023.csv'))








bc=subset(stock_dat,state=='BC') %>% subset(n.years>17)
nrow(subset(bc,species=='Sockeye'))
nrow(subset(bc,species=='Chum'))
nrow(subset(bc,species=='Pink'))
nrow(subset(bc,species=='Coho'))
nrow(subset(bc,species=='Chinook'))

##Barplot of stokcs
stock_dat2<- subset(stock_dat,n.years>17)
sp_sum<- stock_dat2 %>% group_by(species,state) %>% summarize(n=n()) %>% mutate(state2=factor(state,levels=c('OR','WA','BC','AK')))
ggplot(sp_sum, aes(fill=species, y=n, x=state2))+
  scale_fill_manual(values=c("darkgray", "darkgreen", "darkblue","darksalmon","darkred"))+ 
  geom_bar(position="stack", stat="identity") + theme_minimal() +
  xlab('')+ylab('No. time-series')+ theme(text=element_text(size=20), #change font size of all text
                                          axis.text=element_text(size=20))



















####Explorations
#Ogden et al. 2015 - compilation
ogden_info$stock.name=paste(ogden_info$cu.or.pfma,ogden_info$species,sep='-')
#All pink/chum stocks captured in previous compilations, removing these
ogden_info2<- subset(ogden_info,species %in% c('Sockeye','River-Sockeye','Lake Sockeye'))
#Remove series with particularly low data quality (<2)
ogden_info3<- subset(ogden_info2,sqmean>2&cqmean>2&aqmean>2&rqmean>2)
ogden_info3$stock.name;sockeye_info$stock

#PSF data
#PSE data
pse_dat<- read.csv(here('data','raw data','multispecies','PSE_RS.csv'))
pse_spawners<- subset(pse_dat, parameter %in% c('Spawners'))
pse_recruits<- subset(pse_dat, parameter %in% c('Recruits'))
pse_dat2<- left_join(pse_spawners,pse_recruits,by=c('species','location','year'))
pse_dat2$stock<- paste(pse_dat2$species,pse_dat2$location,sep='_')
pse_dat2<- pse_dat2[,c(-1,-6)] #remove excess columns
colnames(pse_dat2)[4:5]<- c('spawners','recruits')

pse_dq<- read.csv(here('data','raw data','multispecies','PSE_data_quality.csv'))
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

#pse_dq_sub1<- subset(pse_dq2,catch_quality>=2&survey_coverage>=2&survey_execution>=2&survey_quality>=2)
pse_dq_sub<- subset(pse_dq2,catch_quality>=3&survey_coverage>=3&survey_execution>=3&survey_quality>=3)
pse_dq_sub$stock.name<- paste(pse_dq_sub$location,pse_dq_sub$species,sep='-')
#All these pink stocks are included in the pink compilation
pse_dg_sub<- subset(pse_dq_sub,species %notin% c('Pink (even)','Pink (odd)'))
pse_dat3<- subset(pse_dat2, stock %in% pse_dq_sub$stock)


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



#Summary overviews for each stock


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
sockeye_info<- read.csv(here('data','raw data','sockeye','sockeye_info.csv'))
colnames(sockeye_info)<- tolower(names(sockeye_info)) 
sockeye_info<- subset(sockeye_info,stock.id %in% sockeye$stock.id) #Removes some extra stocks that do not have actual data
chum_info<- read.csv(here('data','raw data','chum','chum_info.csv'))
chum_info$lon<- -chum_info$lon #make western longitudes negative
pink_info<- read.csv(here('data','raw data','pink','pink_info.csv'))
pink_info$lon<- -pink_info$lon
coho_info<- read.csv(here('data','raw data','coho','coho_info.csv'))
colnames(coho_info)<- tolower(names(coho_info))
chinook_info<- read.csv(here('data','raw data','chinook','chinook_info.csv'))
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

#summary(chinook_filtered$ar1) #mean = 0.33
#summary(chum_filtered$ar1) #mean = 0.3
#summary(coho_filtered$ar1) #mean = 0.17
#summary(pink_filtered$ar1) #mean = -0.01
#summary(sockeye_filtered$ar1) #mean = 0.36

#summary(chinook_filtered$acf_2) #mean = 0.11
#summary(chum_filtered$acf_2) #mean = 0.16
#summary(coho_filtered$acf_2) #mean = 0.03
#summary(pink_filtered$acf_2) #mean = -0.08
#summary(sockeye_filtered$acf_2) #mean = 0.17

#summary(chinook_filtered$acf_3) #mean = -0.04
#summary(chum_filtered$acf_3) #mean = 0.03
#summary(coho_filtered$acf_3) #mean = 0.05
#summary(pink_filtered$acf_3) #mean = -0.02
#summary(sockeye_filtered$acf_3) #mean = 0.05




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

#pse_dq_sub1<- subset(pse_dq2,catch_quality>=2&survey_coverage>=2&survey_execution>=2&survey_quality>=2)
pse_dq_sub<- subset(pse_dq2,catch_quality>=3&survey_coverage>=3&survey_execution>=3&survey_quality>=3)
pse_dq_sub$stock<- paste(pse_dq_sub$species,pse_dq_sub$location,sep='_')
pse_dat3<- subset(pse_dat2, stock %in% pse_dq_sub$stock)

pse_info<- data.frame(index=NA,species=NA,location=NA,ts.length=NA,ts.start=NA,ts.end=NA)
for(i in 1:length(unique(pse_dat3$stock))){
  c<- subset(pse_dat3,stock==unique(pse_dat3$stock)[i])
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
ll_p=subset(all_lat_lon,species=='Pink')
ll_ch=subset(all_lat_lon,species=='Chum')
ll_soc=subset(all_lat_lon,species=='Sockeye')
ll_coh=subset(all_lat_lon,species=='Coho')
ll_chi=subset(all_lat_lon,species=='Chinook')

ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("Longitude") + ylab("Latitude") +
  geom_point(data = ll_p, mapping = aes(x = lon, y = lat), color = "sienna", size = 2.2*log(ll_p$n+1), alpha = 0.7) +
  geom_point(data = ll_ch, mapping = aes(x = lon, y = lat), color = "darkgreen", size = 2.2*log(ll_ch$n+1), alpha = 0.7) +
  geom_point(data =  ll_soc, mapping = aes(x = lon, y = lat), color = "darkred", size = 2.2*log(ll_soc$n+1), alpha = 0.7) +
  geom_point(data =  ll_coh, mapping = aes(x = lon, y = lat), color = "darkblue", size = 2.2*log(ll_coh$n+1), alpha = 0.7) +
  geom_point(data = ll_chi, mapping = aes(x = lon, y = lat), color = "darkgray", size = 2.2*log(ll_chi$n+1), alpha = 0.7) +
    annotate(geom = 'text', x = -155, y = 54.5, label = 'Chinook', color = 'darkgray', size = 3) + 
  annotate(geom = 'text', x = -155, y = 53.5, label = 'Chum', color = 'darkgreen', size = 3) + 
  annotate(geom = 'text', x = -155, y = 52.5, label = 'Coho', color = 'darkblue', size = 3) + 
  annotate(geom = 'text', x = -155, y = 51.5, label = 'Pink', color = 'sienna', size = 3) + 
  annotate(geom = 'text', x = -155, y = 50.5, label = 'Sockeye', color = 'darkred', size = 3) + 
  annotate(geom = 'point', x = -145, y = 53.5, size = 2.2*log(1+1)) +
  annotate(geom = 'point', x = -145, y = 52.5, size = 2.2*log(3+1)) + 
  annotate(geom = 'point', x = -145, y = 51.5, size = 2.2*log(10+1)) + 
  annotate(geom = 'point', x = -145, y = 50.5, size = 2.2*log(20+1)) + 
  annotate(geom = 'text', x = -143, y = 54.5,label='n', size = 2.5) +
  annotate(geom = 'text', x = -143, y = 53.5,label='1', size = 2.5) +
  annotate(geom = 'text', x = -143, y = 52.5,label='3', size = 2.5) + 
  annotate(geom = 'text', x = -143, y = 51.5,label='10', size = 2.5) + 
  annotate(geom = 'text', x = -143, y = 50.5,label='20', size = 2.5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.75, 'in'), pad_y = unit(0.5, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -120), ylim = c(47, 72), expand = T) +
theme(panel.grid.major = element_line(color = gray(.5), linetype = 'dashed', size = 0.5), panel.background = element_rect(fill = 'lightblue'),legend.title = element_blank())

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
