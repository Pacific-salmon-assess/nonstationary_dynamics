#Salmon spawner=recruit time-series; dan.greenberg@dfo-mpo.gc.ca
library(here);library(dplyr);library(ggplot2);library(ggthemes);library(ggspatial);library(sf);library(rnaturalearth);library(rnaturalearthdata);library(wesanderson)
library(nlme)
#plot functions
source(here('code','plot_functions.R'))

#read in data for each species
sockeye<- read.csv(here('data','sockeye','sockeye_data.csv'))
chum<- read.csv(here('data','chum','chum_data.csv'))
pink<- read.csv(here('data','pink','pink_data.csv'))
coho_1<- read.csv(here('data','coho','coho_data_fwmar.csv'))
coho_2<- read.csv(here('data','coho','coho_data_totalage.csv'))
chinook_1<- read.csv(here('data','chinook','chinook_data_fwmar.csv'))
chinook_2<- read.csv(here('data','chinook','chinook_data_totalage.csv'))

#Synonymize column names so I don't go crazy
colnames(sockeye)<- tolower(names(sockeye))
names(chum)<- gsub('.yr','year',names(chum)) #make the brood years equivalent
names(pink)<- gsub('.yr','year',names(pink))  #make the brood years equivalent
colnames(coho_1)<- tolower(names(coho_1))
colnames(coho_2)<- tolower(names(coho_2))
colnames(chinook_1)<- tolower(names(chinook_1))
colnames(chinook_2)<- tolower(names(chinook_2))

#Number of Stock time-series for each species (before filtering)
length(unique(sockeye$stock.id)) #51
length(unique(chum$stock.id)) #53
length(unique(pink$stock.id)) #46
length(unique(coho_1$stock.id)) #13
length(unique(chinook_2$stock)) #20

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
chinook_list<- list()
for(i in 1:length(unique(chinook_info$stock.id))){
  c<- chinook_2 %>% subset(stock.id==unique(stock.id)[i]) %>% subset(useflag==1)
  c$logR_S<- log(c$recruits/c$spawners)
  c_sr<- c[complete.cases(c$logR_S),]
  chinook_info$ts.length[i]=nrow(c_sr)
  chinook_info$ts.start[i]=min(c_sr$broodyear)
  chinook_info$ts.end[i]=max(c_sr$broodyear)
  
  if(chinook_info$ts.length[i]>=20){ #only keep decently long (min. 20 year) timeseries
    chinook_list[[i]]=c
  }else{next}
  
  sr_fit_ac<- gls(logR_S~spawners,data=c_sr,correlation = corAR1(form = ~ 1),method='ML')
  sr_fit<- gls(logR_S~spawners,data=c_sr,method='ML')
  acf_sr<- acf(residuals(sr_fit))
  chinook_info$logA[i]<- sr_fit_ac$coefficients[1]
  chinook_info$b[i]<- -sr_fit_ac$coefficients[2]
  chinook_info$K[i]<- 1/-sr_fit_ac$coefficients[2]
  chinook_info$rmax[i]<- (exp(sr_fit_ac$coefficients[1])/-sr_fit_ac$coefficients[2])*exp(-1)
  chinook_info$sigma[i]<- sr_fit_ac$sigma
  chinook_info$ar1[i]<-intervals(sr_fit_ac)$corStruct[,2] #ML estimate of AR1 correlation
  chinook_info$ar1.l95[i]<-intervals(sr_fit_ac)$corStruct[,1] #l95% CI of ML estimate
  chinook_info$ar1.u95[i]<-intervals(sr_fit_ac)$corStruct[,3] #u95% CI of ML estimate
  chinook_info$p_S_belowK[i]<- nrow(subset(c,spawners<(1/-sr_fit_ac$coefficients[2])))/nrow(c)
  chinook_info$logA2[i]<- sr_fit$coefficients[1]
  chinook_info$b2[i]<- -sr_fit$coefficients[2]
  chinook_info$K2[i]<- 1/-sr_fit$coefficients[2]
  chinook_info$rmax2[i]<- (exp(sr_fit$coefficients[1])/-sr_fit$coefficients[2])*exp(-1)
  chinook_info$sigma2[i]<- sr_fit$sigma
  chinook_info$acf_1[i]<- acf_sr$acf[2]
  chinook_info$acf_2[i]<- acf_sr$acf[3]
  chinook_info$acf_3[i]<- acf_sr$acf[4]
  chinook_info$acf_4[i]<- acf_sr$acf[5]
  
  plot_SR(c_sr,sr_fit_ac,path=here('outputs','figures','stock-recruit','chinook'))
}

chum_list<- list()
for(i in 1:length(unique(chum_info$stock.id))){
  c<- chum %>% subset(stock.id==unique(stock.id)[i]) %>% subset(use==1)
  c$logR_S<- log(c$recruits/c$spawners)
  chum_list[[i]]=c
  c_sr<- c[complete.cases(c$logR_S),]
  chum_info$ts.length[i]=nrow(c_sr)
  chum_info$ts.start[i]=min(c_sr[,3])
  chum_info$ts.end[i]=max(c_sr[,3])
  
  sr_fit_ac<- gls(logR_S~spawners,data=c_sr,correlation = corAR1(form = ~ 1),method='ML')
  sr_fit<- gls(logR_S~spawners,data=c_sr,method='ML')
  acf_sr<- acf(residuals(sr_fit))
  chum_info$logA[i]<- sr_fit_ac$coefficients[1]
  chum_info$b[i]<- -sr_fit_ac$coefficients[2]
  chum_info$K[i]<- 1/-sr_fit_ac$coefficients[2]
  chum_info$rmax[i]<- (exp(sr_fit_ac$coefficients[1])/-sr_fit_ac$coefficients[2])*exp(-1)
  chum_info$sigma[i]<- sr_fit_ac$sigma
  chum_info$ar1[i]<-intervals(sr_fit_ac)$corStruct[,2] #ML estimate of AR1 correlation
  chum_info$ar1.l95[i]<-intervals(sr_fit_ac)$corStruct[,1] #l95% CI of ML estimate
  chum_info$ar1.u95[i]<-intervals(sr_fit_ac)$corStruct[,3] #u95% CI of ML estimate
  chum_info$p_S_belowK[i]<- nrow(subset(c,spawners<(1/-sr_fit_ac$coefficients[2])))/nrow(c)
  chum_info$logA2[i]<- sr_fit$coefficients[1]
  chum_info$b2[i]<- -sr_fit$coefficients[2]
  chum_info$K2[i]<- 1/-sr_fit$coefficients[2]
  chum_info$rmax2[i]<- (exp(sr_fit$coefficients[1])/-sr_fit$coefficients[2])*exp(-1)
  chum_info$sigma2[i]<- sr_fit$sigma
  chum_info$acf_1[i]<- acf_sr$acf[2]
  chum_info$acf_2[i]<- acf_sr$acf[3]
  chum_info$acf_3[i]<- acf_sr$acf[4]
  chum_info$acf_4[i]<- acf_sr$acf[5]
  
  plot_SR(c_sr,sr_fit_ac,path=here('outputs','figures','stock-recruit','chum'))
}


coho_list<- list()
for(i in 1:length(unique(coho_info$stock.id))){
  c<- coho_1 %>% subset(stock.id==unique(stock.id)[i]) %>% subset(useflag==1)
  c$logR_S<- log(c$recruits/c$spawners)
  coho_list[[i]]=c
  c_sr<- c[complete.cases(c$logR_S),]
  coho_info$ts.length[i]=nrow(c_sr)
  coho_info$ts.start[i]=min(c_sr[,3])
  coho_info$ts.end[i]=max(c_sr[,3])
  
  sr_fit_ac<- gls(logR_S~spawners,data=c_sr,correlation = corAR1(form = ~ 1),method='ML')
  sr_fit<- gls(logR_S~spawners,data=c_sr,method='ML')
  acf_sr<- acf(residuals(sr_fit))
  coho_info$logA[i]<- sr_fit_ac$coefficients[1]
  coho_info$b[i]<- -sr_fit_ac$coefficients[2]
  coho_info$K[i]<- 1/-sr_fit_ac$coefficients[2]
  coho_info$rmax[i]<- (exp(sr_fit_ac$coefficients[1])/-sr_fit_ac$coefficients[2])*exp(-1)
  coho_info$sigma[i]<- sr_fit_ac$sigma
  coho_info$ar1[i]<-intervals(sr_fit_ac)$corStruct[,2] #ML estimate of AR1 correlation
  coho_info$ar1.l95[i]<-intervals(sr_fit_ac)$corStruct[,1] #l95% CI of ML estimate
  coho_info$ar1.u95[i]<-intervals(sr_fit_ac)$corStruct[,3] #u95% CI of ML estimate
  coho_info$p_S_belowK[i]<- nrow(subset(c,spawners<(1/-sr_fit_ac$coefficients[2])))/nrow(c)
  coho_info$logA2[i]<- sr_fit$coefficients[1]
  coho_info$b2[i]<- -sr_fit$coefficients[2]
  coho_info$K2[i]<- 1/-sr_fit$coefficients[2]
  coho_info$rmax2[i]<- (exp(sr_fit$coefficients[1])/-sr_fit$coefficients[2])*exp(-1)
  coho_info$sigma2[i]<- sr_fit$sigma
  coho_info$acf_1[i]<- acf_sr$acf[2]
  coho_info$acf_2[i]<- acf_sr$acf[3]
  coho_info$acf_3[i]<- acf_sr$acf[4]
  coho_info$acf_4[i]<- acf_sr$acf[5]
  
  plot_SR(c_sr,sr_fit_ac,path=here('outputs','figures','stock-recruit','coho'))
}

pink_list<- list()
for(i in 1:length(unique(pink_info$stock.id))){
  c<- subset(pink,stock.id==unique(pink_info$stock.id)[i])  %>% subset(use==1)
  c$logR_S<- log(c$recruits/c$spawners)
  c_sr<- c[complete.cases(c$logR_S),]
  pink_list[[i]]=c
  pink_info$ts.length[i]=nrow(c_sr)
  pink_info$ts.start[i]=min(c_sr[,3])
  pink_info$ts.end[i]=max(c_sr[,3])
  
  sr_fit_ac<- gls(logR_S~spawners,data=c_sr,correlation = corAR1(form = ~ 1),method='ML')
  sr_fit<- gls(logR_S~spawners,data=c_sr,method='ML')
  acf_sr<- acf(residuals(sr_fit))
  pink_info$logA[i]<- sr_fit_ac$coefficients[1]
  pink_info$b[i]<- -sr_fit_ac$coefficients[2]
  pink_info$K[i]<- 1/-sr_fit_ac$coefficients[2]
  pink_info$rmax[i]<- (exp(sr_fit_ac$coefficients[1])/-sr_fit_ac$coefficients[2])*exp(-1)
  pink_info$sigma[i]<- sr_fit_ac$sigma
  pink_info$ar1[i]<-intervals(sr_fit_ac)$corStruct[,2] #ML estimate of AR1 correlation
  pink_info$ar1.l95[i]<-intervals(sr_fit_ac)$corStruct[,1] #l95% CI of ML estimate
  pink_info$ar1.u95[i]<-intervals(sr_fit_ac)$corStruct[,3] #u95% CI of ML estimate
  pink_info$p_S_belowK[i]<- nrow(subset(c,spawners<(1/-sr_fit_ac$coefficients[2])))/nrow(c)
  pink_info$logA2[i]<- sr_fit$coefficients[1]
  pink_info$b2[i]<- -sr_fit$coefficients[2]
  pink_info$K2[i]<- 1/-sr_fit$coefficients[2]
  pink_info$rmax2[i]<- (exp(sr_fit$coefficients[1])/-sr_fit$coefficients[2])*exp(-1)
  pink_info$sigma2[i]<- sr_fit$sigma
  pink_info$acf_1[i]<- acf_sr$acf[2]
  pink_info$acf_2[i]<- acf_sr$acf[3]
  pink_info$acf_3[i]<- acf_sr$acf[4]
  pink_info$acf_4[i]<- acf_sr$acf[5]
  
  plot_SR(c_sr,sr_fit_ac,path=here('outputs','figures','stock-recruit','pink'))
}

sockeye_list<- list()
for(i in 1:length(unique(sockeye_info$stock.id))){
  c<- subset(sockeye,stock.id==unique(sockeye_info$stock.id)[i])  %>% subset(useflag==1)
  if(nrow(c)==0){next}
  c$logR_S<- log(c$recruits/c$spawners)
  sockeye_list[[i]]=c
  c_sr<- c[complete.cases(c$logR_S),]
  sockeye_info$ts.length[i]=nrow(c_sr)
  sockeye_info$ts.start[i]=min(c_sr[,3])
  sockeye_info$ts.end[i]=max(c_sr[,3])
  
  sr_fit_ac<- gls(logR_S~spawners,data=c_sr,correlation = corAR1(form = ~ 1),method='ML')
  sr_fit<- gls(logR_S~spawners,data=c_sr,method='ML')
  acf_sr<- acf(residuals(sr_fit))
  sockeye_info$logA[i]<- sr_fit_ac$coefficients[1]
  sockeye_info$b[i]<- -sr_fit_ac$coefficients[2]
  sockeye_info$K[i]<- 1/-sr_fit_ac$coefficients[2]
  sockeye_info$rmax[i]<- (exp(sr_fit_ac$coefficients[1])/-sr_fit_ac$coefficients[2])*exp(-1)
  sockeye_info$sigma[i]<- sr_fit_ac$sigma
  sockeye_info$ar1[i]<-intervals(sr_fit_ac)$corStruct[,2] #ML estimate of AR1 correlation
  sockeye_info$ar1.l95[i]<-intervals(sr_fit_ac)$corStruct[,1] #l95% CI of ML estimate
  sockeye_info$ar1.u95[i]<-intervals(sr_fit_ac)$corStruct[,3] #u95% CI of ML estimate
  sockeye_info$p_S_belowK[i]<- nrow(subset(c,spawners<(1/-sr_fit_ac$coefficients[2])))/nrow(c)
  sockeye_info$logA2[i]<- sr_fit$coefficients[1]
  sockeye_info$b2[i]<- -sr_fit$coefficients[2]
  sockeye_info$K2[i]<- 1/-sr_fit$coefficients[2]
  sockeye_info$rmax2[i]<- (exp(sr_fit$coefficients[1])/-sr_fit$coefficients[2])*exp(-1)
  sockeye_info$sigma2[i]<- sr_fit$sigma
  sockeye_info$acf_1[i]<- acf_sr$acf[2]
  sockeye_info$acf_2[i]<- acf_sr$acf[3]
  sockeye_info$acf_3[i]<- acf_sr$acf[4]
  sockeye_info$acf_4[i]<- acf_sr$acf[5]
  
  plot_SR(c_sr,sr_fit_ac,path=here('outputs','figures','stock-recruit','sockeye'))
}

#Time-series overlap
par(mfrow=c(3,2))
hist(chinook_info$ts.length)
summary(chinook_info$ts.length)
nrow(subset(chinook_info,ts.length>20))

hist(chum_info$ts.length)
summary(coho_info$ts.length)
nrow(subset(chum_info,ts.length>20))

hist(coho_info$ts.length)
summary(coho_info$ts.length)
nrow(subset(coho_info,ts.length>20))

hist(pink_info$ts.length)
summary(pink_info$ts.length)
nrow(subset(pink_info,ts.length>20))

hist(sockeye_info$ts.length)
summary(sockeye_info$ts.length)
nrow(subset(sockeye_info,ts.length>20))

dev.off()
#autocorrelation in residual productivity

par(mfrow=c(3,2))
plot_acf(chinook_info)
plot_acf(chum_info)
plot_acf(coho_info)
plot_acf(pink_info)
plot_acf(sockeye_info)

summary(chinook_info$ar1)
summary(chum_info$ar1)
summary(coho_info$ar1)
summary(pink_info$ar1)
summary(sockeye_info$ar1)

#Combined coordinate data
#Lat lons for each species
sockeye_lat_lon<- cbind(sockeye_info$lat,sockeye_info$lon)
chinook_lat_lon<- cbind(chinook_info$lat,chinook_info$lon)
chum_lat_lon<- cbind(chum_info$lat,chum_info$lon)
coho_lat_lon<- cbind(coho_info$lat,coho_info$lon)
pink_lat_lon<- cbind(pink_info$lat,pink_info$lon)

all_lat_lon<- rbind(sockeye_lat_lon,chinook_lat_lon,chum_lat_lon,coho_lat_lon,pink_lat_lon)
length(unique(all_lat_lon[,1])); length(unique(all_lat_lon[,2])) #124 unique sites

##
world <- ne_countries(scale = "medium", returnclass = "sf")
           
ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("Longitude") + ylab("Latitude") +
  geom_point(data = pink_info, mapping = aes(x = lon, y = lat), color = "salmon", size = 3, alpha = 0.7) +
  geom_point(data = chum_info, mapping = aes(x = lon, y = lat), color = "darkgreen", size = 3, alpha = 0.7) +
  geom_point(data = sockeye_info, mapping = aes(x = lon, y = lat), color = "darkred", size = 3, alpha = 0.7) +
  geom_point(data = coho_info, mapping = aes(x = lon, y = lat), color = "darkblue", size = 3, alpha = 0.7) +
  geom_point(data = chinook_info, mapping = aes(x = lon, y = lat), color = "darkgray", size = 3, alpha = 0.7) +
#  geom_text(data= world_points,aes(x=X, y=Y, label=name), color = 'darkblue', fontface = 'bold', check_overlap = FALSE) + 
#  annotate(geom = “text”, x = -90, y = 26, label = “Gulf of Mexico”, fontface = “italic”, color = “grey22”, size = 6) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.75, 'in'), pad_y = unit(0.5, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -120), ylim = c(47, 72), expand = T) + 
  #xlab(“Longitude”) + 
#  ylab(“Latitude”) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = 'dashed', size = 0.5), panel.background = element_rect(fill = 'lightblue'),legend.title = element_blank())+
  scale_fill_manual(values = c('darkgray', 'darkgreen', 'darkblue','salmon','darkred'),
                    labels = c('Chinook','Chum','Coho','Pink','Sockeye'))

##To do on this:
#Remove junky time-series (ie. filter data FIRST)
#make it easier to see points with multiple species?
#add legend
#add salmon pics?






