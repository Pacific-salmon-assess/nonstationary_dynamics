library(here);library(ggplot2);library(dplyr);library(cowplot);library(ggspatial);library(sf);library(rnaturalearth);library(rnaturalearthdata)
options(dplyr.summarise.inform = FALSE)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation2023-10-18.csv'))
stock_dat$stock=gsub('/','_',stock_dat$stock)
stock_dat$stock=gsub('&','and',stock_dat$stock)
stock_info<- read.csv(here('data','filtered datasets','stock_info2023-10-18.csv'))
#Remove stocks with less than 15 years of recruitment data
stock_info_filtered=subset(stock_info,n.years>=16) #242 stocks
stock_info_filtered$stock.name=gsub('/','_',stock_info_filtered$stock.name)
stock_info_filtered$stock.name=gsub('&','and',stock_info_filtered$stock.name)

stock_info_filtered$state=factor(stock_info_filtered$state,levels=c('OR','WA','BC','AK'))
stock_info_filtered$ocean.basin=factor(stock_info_filtered$ocean.basin,levels=c('WC','GOA','BS'))

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock)) #267

stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)

stock_info_filtered$state2=as.character(stock_info_filtered$state)
stock_info_filtered=stock_info_filtered %>% mutate(state2=ifelse(state=='OR','OR-WA',state2),
                                                   state2=ifelse(state=='WA','OR-WA',state2))
rownames(stock_info_filtered)=seq(1:nrow(stock_info_filtered))
stock_info_filtered$species2=gsub("-.*", "", stock_info_filtered$species) # collapse even/odd pinks


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

sp_cols=cbind(c('#bdc9e1',
                '#016c59',
                '#1c9099',
                '#fd8d3c',
                '#b30000'),c('Chinook','Chum','Coho','Pink','Sockeye'))

mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

world <- ne_countries(scale = "medium", returnclass = "sf")

##adding watersheds? will try later
#hyb3=sf::st_read(here('data','geographic','hybas_lake_na_lev03_v1c.shp'))
#hyb3ar=sf::st_read(here('data','geographic','hybas_lake_ar_lev03_v1c.shp'))
#hyb4=sf::st_read(here('data','geographic','hybas_lake_na_lev04_v1c.shp'))
#hyb4ar=sf::st_read(here('data','geographic','hybas_lake_ar_lev04_v1c.shp'))
#hyb4_2=sf::st_join(hyb4,hyb4ar,join = st_intersects)
#hyb5=sf::st_read(here('data','geographic','hybas_na_lev04_v1c.shp'))
#hy_r=sf::st_read(here('data','geographic','HydroRIVERS_v10_na.shp'))
#test=hy_r %>% dplyr::select(ORD_FLOW) %>% sf::st_intersection(hyb3)
#hy_l=sf::st_read(here('data','geographic','HydroLAKES_polys_v10.shp'))

all_lat_lon = stock_info_filtered %>% group_by(lat,lon,species2) %>% summarize(n=n(),.groups = 'keep')
all_lat_lon$lon=ifelse(all_lat_lon$lon>0,-all_lat_lon$lon,all_lat_lon$lon)
ll_p=subset(all_lat_lon,species2=='Pink')
ll_ch=subset(all_lat_lon,species2=='Chum')
ll_soc=subset(all_lat_lon,species2=='Sockeye')
ll_coh=subset(all_lat_lon,species2=='Coho')
ll_chi=subset(all_lat_lon,species2=='Chinook')

ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("Longitude") + ylab("Latitude") +
#  +geom_sf(data=hyb3,aes(color=factor(HYBAS_ID)))+
  geom_point(data = ll_p, mapping = aes(x = lon, y = lat), color = sp_cols[4], size = 2.2*log(ll_p$n+1), alpha = 0.7) +
  geom_point(data = ll_ch, mapping = aes(x = lon, y = lat), color = sp_cols[2], size = 2.2*log(ll_ch$n+1), alpha = 0.7) +
  geom_point(data =  ll_soc, mapping = aes(x = lon, y = lat), color = sp_cols[5], size = 2.2*log(ll_soc$n+1), alpha = 0.7) +
  geom_point(data =  ll_coh, mapping = aes(x = lon, y = lat), color = sp_cols[3], size = 2.2*log(ll_coh$n+1), alpha = 0.7) +
  geom_point(data = ll_chi, mapping = aes(x = lon, y = lat), color = sp_cols[1], size = 2.2*log(ll_chi$n+1), alpha = 0.7) +
  annotate(geom = 'text', x = -155, y = 54.5, label = 'Chinook', color = sp_cols[1], size = 4.5) + 
  annotate(geom = 'text', x = -155, y = 53, label = 'Chum', color = sp_cols[2], size = 4.5) + 
  annotate(geom = 'text', x = -155, y = 51.5, label = 'Coho', color = sp_cols[3], size = 4.5) + 
  annotate(geom = 'text', x = -155, y = 50, label = 'Pink', color = sp_cols[4], size = 4.5) + 
  annotate(geom = 'text', x = -155, y = 48.5, label = 'Sockeye', color = sp_cols[5], size = 4.5) + 
  annotate(geom = 'point', x = -145, y = 54, size = 2.2*log(1+1)) +
  annotate(geom = 'point', x = -145, y = 53, size = 2.2*log(3+1)) + 
  annotate(geom = 'point', x = -145, y = 51.5, size = 2.2*log(10+1)) + 
  annotate(geom = 'point', x = -145, y = 50, size = 2.2*log(20+1)) + 
  annotate(geom = 'text', x = -142, y = 55,label='n', size = 3.5) +
  annotate(geom = 'text', x = -142, y = 54,label='1', size = 3.5) +
  annotate(geom = 'text', x = -142, y = 53,label='3', size = 3.5) + 
  annotate(geom = 'text', x = -142, y = 51.5,label='10', size = 3.5) + 
  annotate(geom = 'text', x = -142, y = 50,label='20', size = 3.5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -120), ylim = c(47, 72), expand = T) +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = 'dashed', size = 0.1), panel.background = element_rect(fill = 'lightblue1'),legend.title = element_blank())

sp_sum<- stock_info_filtered %>% group_by(species2,state) %>% summarize(n=n(),.groups='keep')
sp_sum$state=factor(sp_sum$state,levels=c('OR','WA','BC','AK'))

ggplot(sp_sum, aes(fill=species2, y=n, x=state))+
  scale_fill_manual(values=sp_cols[match(levels(factor(sp_sum$species2)),sp_cols[,2])])+ 
  geom_bar(position="stack", stat="identity") +
  xlab('')+ylab('No. time-series')+
  mytheme

sp_sum<- stock_info_filtered %>% group_by(species2,ocean.basin) %>% summarize(n=n(),.groups='keep')

ggplot(sp_sum, aes(fill=species2, y=n, x=ocean.basin))+
  scale_fill_manual(values=sp_cols[match(levels(factor(sp_sum$species2)),sp_cols[,2])])+ 
  geom_bar(position="stack", stat="identity") + theme_minimal() +
  xlab('')+ylab('No. time-series')+
  mytheme

sp_sum<- stock_info_filtered %>% group_by(n.years,species2) %>% summarize(n=n(),.groups='keep')

ggplot(sp_sum, aes(fill=species2, y=n, x=n.years))+
  scale_fill_manual(values=sp_cols[match(levels(factor(sp_sum$species2)),sp_cols[,2])])+ 
  geom_bar(position="stack", stat="identity") + theme_minimal() +
  xlab('Years')+ylab('No. time-series')+
  mytheme

pl2=list() #params list - m2
pl3=list() #params list - m3

for(i in 1:nrow(stock_info_filtered)){
  pl2[[i]]=read.csv(here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m2.csv'))  
  pl3[[i]]=read.csv(here('outputs','parameters',paste(sprintf("%03d", i),stock_info_filtered$stock.name[i],sep=''),'m3.csv'))  
}

prod_matr=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_dat2$broodyear)-min(stock_dat2$broodyear)+1)
colnames(prod_matr)=seq(min(stock_dat2$broodyear),max(stock_dat2$broodyear))
prod_matr2=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_dat2$broodyear)-min(stock_dat2$broodyear)+1)
colnames(prod_matr2)=seq(min(stock_dat2$broodyear),max(stock_dat2$broodyear))
umsy_matr=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_dat2$broodyear)-min(stock_dat2$broodyear)+1)
colnames(umsy_matr)=seq(min(stock_dat2$broodyear),max(stock_dat2$broodyear))
umsy_matr2=matrix(nrow=nrow(stock_info_filtered),ncol=max(stock_dat2$broodyear)-min(stock_dat2$broodyear)+1)
colnames(umsy_matr2)=seq(min(stock_dat2$broodyear),max(stock_dat2$broodyear))

for(i in 1:nrow(stock_info_filtered)){
  s<- subset(stock_dat2,stock==stock_info_filtered$stock.name[i])
  s<- s[complete.cases(s$logR_S),] 
  
  df=data.frame(S=s$spawners,R=s$recruits,by=s$broodyear,logRS=s$logR_S)
  
  prod=pl3[[i]][,grepl('log_a',colnames(pl3[[i]]))]
  prod_m=apply(prod,2,median)
  prod_scale=prod-median(unlist(prod))
  prod_m_scale=apply(prod_scale,2,median) #scaled as multiplier of median
  
  prod_matr[i,match(seq(min(df$by),max(df$by)),colnames(prod_matr))]=prod_m
  prod_matr2[i,match(seq(min(df$by),max(df$by)),colnames(prod_matr2))]=prod_m_scale
  
  umsy=pl3[[i]][,grepl('Umsy',colnames(pl3[[i]]))]
  umsy_m=apply(umsy,2,median)
  umsy_scale=umsy-median(unlist(umsy))
  umsy_m_scale=apply(umsy_scale,2,median) #scaled as multiplier of median
  
  umsy_matr[i,match(seq(min(df$by),max(df$by)),colnames(umsy_matr))]=umsy_m
  umsy_matr2[i,match(seq(min(df$by),max(df$by)),colnames(umsy_matr2))]=umsy_m_scale
}

#overall prod
plot(rep(0,ncol(prod_matr))~colnames(prod_matr),ylim=c(0,25),ylab='Max. recruits/spawner',type='n',bty='l',xlab='Brood cohort year',xaxt='n')
axis(side=1,at=seq(1950,2020,by=5))
for(i in 1:nrow(prod_matr)){
  lines(exp(prod_matr[i,])~colnames(prod_matr),col=adjustcolor('darkgray',alpha.f=0.3))
}
m.prod.5f=move_avg_mat(exp(prod_matr));m.prod.5f$y=colnames(prod_matr)[1:c(ncol(prod_matr)-2)]
m.prod.5f=m.prod.5f[complete.cases(m.prod.5f[,1]),]
x<- c(seq(m.prod.5f$y[1],m.prod.5f$y[nrow(m.prod.5f)]), rev(seq(m.prod.5f$y[1],m.prod.5f$y[nrow(m.prod.5f)])))
y<- c(m.prod.5f[,1]-2*m.prod.5f[,2], rev(m.prod.5f[,1]+2*m.prod.5f[,2]))
polygon(x, y, col = adjustcolor('navy', alpha = 0.3), border=NA) # Add uncertainty polygon
lines(m.prod.5f[,1]~seq(min(x),max(x)),lwd=2,col='navy')

#scaled
plot(rep(0,ncol(prod_matr2))~colnames(prod_matr2),ylim=c(-3,3),ylab='Deviations from median log(alpha)',type='n',bty='l',xlab='Brood cohort year',xaxt='n')
axis(side=1,at=seq(1950,2020,by=5))
for(i in 1:nrow(prod_matr2)){
  lines(prod_matr2[i,]~colnames(prod_matr2),col=adjustcolor('darkgray',alpha.f=0.3))
}
m.prod.5f=move_avg_mat(prod_matr2);m.prod.5f$y=colnames(prod_matr)[1:c(ncol(prod_matr)-2)]
m.prod.5f=m.prod.5f[complete.cases(m.prod.5f[,1]),]
x<- c(seq(m.prod.5f$y[1],m.prod.5f$y[nrow(m.prod.5f)]), rev(seq(m.prod.5f$y[1],m.prod.5f$y[nrow(m.prod.5f)])))
y<- c(m.prod.5f[,1]-2*m.prod.5f[,2], rev(m.prod.5f[,1]+2*m.prod.5f[,2]))
polygon(x, y, col = adjustcolor('navy', alpha = 0.3), border=NA) # Add uncertainty polygon
lines(m.prod.5f[,1]~seq(min(x),max(x)),lwd=2,col='navy')

#umsy
plot(rep(0,ncol(umsy_matr))~colnames(umsy_matr),ylim=c(0,1),ylab='Sustainable harvest rate (Umsy)',type='n',bty='l',xlab='Brood cohort year',xaxt='n')
axis(side=1,at=seq(1950,2020,by=5))
for(i in 1:nrow(umsy_matr)){
  lines(umsy_matr[i,]~colnames(umsy_matr),col=adjustcolor('darkgray',alpha.f=0.3))
}
m.umsy.5f=move_avg_mat(umsy_matr);m.umsy.5f$y=colnames(umsy_matr)[1:c(ncol(umsy_matr)-2)]
m.umsy.5f=m.umsy.5f[complete.cases(m.umsy.5f[,1]),]
x<- c(seq(m.umsy.5f$y[1],m.umsy.5f$y[nrow(m.umsy.5f)]), rev(seq(m.umsy.5f$y[1],m.umsy.5f$y[nrow(m.umsy.5f)])))
y<- c(m.umsy.5f[,1]-2*m.umsy.5f[,2], rev(m.umsy.5f[,1]+2*m.umsy.5f[,2]))
polygon(x, y, col = adjustcolor('navy', alpha = 0.3), border=NA) # Add uncertainty polygon
lines(m.umsy.5f[,1]~seq(min(x),max(x)),lwd=2,col='navy')

plot(rep(0,ncol(umsy_matr2))~colnames(umsy_matr2),ylim=c(-1,1),ylab='Deviations from median Umsy',type='n',bty='l',xlab='Brood cohort year',xaxt='n')
axis(side=1,at=seq(1950,2020,by=5))
for(i in 1:nrow(umsy_matr2)){
  lines(umsy_matr2[i,]~colnames(umsy_matr2),col=adjustcolor('darkgray',alpha.f=0.3))
}
m.umsy.5f=move_avg_mat(umsy_matr2);m.umsy.5f$y=colnames(umsy_matr)[1:c(ncol(umsy_matr)-2)]
m.umsy.5f=m.umsy.5f[complete.cases(m.umsy.5f[,1]),]
x<- c(seq(m.umsy.5f$y[1],m.umsy.5f$y[nrow(m.umsy.5f)]), rev(seq(m.umsy.5f$y[1],m.umsy.5f$y[nrow(m.umsy.5f)])))
y<- c(m.umsy.5f[,1]-2*m.umsy.5f[,2], rev(m.umsy.5f[,1]+2*m.umsy.5f[,2]))
polygon(x, y, col = adjustcolor('navy', alpha = 0.3), border=NA) # Add uncertainty polygon
lines(m.umsy.5f[,1]~seq(min(x),max(x)),lwd=2,col='navy')

#prod matrices - breakdown
chi_alpha=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook',])),]
chi_alpha_goa=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$ocean.basin=='GOA',])),]
chi_alpha_bs=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$ocean.basin=='BS',])),]
chi_alpha_wc=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$ocean.basin=='WC',])),]
chi_alpha_ak=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$state2=='AK',])),]
chi_alpha_bc=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$state2=='BC',])),]
chi_alpha_orwa=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$state2=='OR-WA',])),]

chu_alpha=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum',])),]
chu_alpha_goa=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$ocean.basin=='GOA',])),]
chu_alpha_bs=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$ocean.basin=='BS',])),]
chu_alpha_wc=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$ocean.basin=='WC',])),]
chu_alpha_ak=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$state2=='AK',])),]
chu_alpha_bc=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$state2=='BC',])),]
chu_alpha_orwa=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$state2=='OR-WA',])),]

coh_alpha=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho',])),]
coh_alpha_goa=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$ocean.basin=='GOA',])),]
coh_alpha_bs=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$ocean.basin=='BS',])),]
coh_alpha_wc=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$ocean.basin=='WC',])),]
coh_alpha_ak=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$state2=='AK',])),]
coh_alpha_bc=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$state2=='BC',])),]
coh_alpha_orwa=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$state2=='OR-WA',])),]

pi_alpha=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species2=='Pink',])),]
pi_alpha_goa=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$ocean.basin=='GOA',])),]
pi_alpha_bs=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$ocean.basin=='BS',])),]
pi_alpha_wc=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$ocean.basin=='WC',])),]
pi_alpha_ak=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$state2=='AK',])),]
pi_alpha_bc=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$state2=='BC',])),]
pi_alpha_orwa=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$state2=='OR-WA',])),]

soc_alpha=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species2=='Sockeye',])),]
soc_alpha_goa=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$ocean.basin=='GOA',])),]
soc_alpha_bs=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$ocean.basin=='BS',])),]
soc_alpha_wc=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$ocean.basin=='WC',])),]
soc_alpha_ak=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$state2=='AK',])),]
soc_alpha_bc=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$state2=='BC',])),]
soc_alpha_orwa=prod_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$state2=='OR-WA',])),]

#prod scaled matrices
chi_alpha2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook',])),]
chi_alpha_goa2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$ocean.basin=='GOA',])),]
chi_alpha_bs2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$ocean.basin=='BS',])),]
chi_alpha_wc2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$ocean.basin=='WC',])),]
chi_alpha_ak2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$state2=='AK',])),]
chi_alpha_bc2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$state2=='BC',])),]
chi_alpha_orwa2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$state2=='OR-WA',])),]

chu_alpha2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum',])),]
chu_alpha_goa2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$ocean.basin=='GOA',])),]
chu_alpha_bs2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$ocean.basin=='BS',])),]
chu_alpha_wc2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$ocean.basin=='WC',])),]
chu_alpha_ak2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$state2=='AK',])),]
chu_alpha_bc2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$state2=='BC',])),]
chu_alpha_orwa2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$state2=='OR-WA',])),]

coh_alpha2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho',])),]
coh_alpha_goa2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$ocean.basin=='GOA',])),]
coh_alpha_bs2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$ocean.basin=='BS',])),]
coh_alpha_wc2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$ocean.basin=='WC',])),]
coh_alpha_ak2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$state2=='AK',])),]
coh_alpha_bc2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$state2=='BC',])),]
coh_alpha_orwa2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$state2=='OR-WA',])),]

pi_alpha2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species2=='Pink',])),]
pi_alpha_goa2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$ocean.basin=='GOA',])),]
pi_alpha_bs2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$ocean.basin=='BS',])),]
pi_alpha_wc2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$ocean.basin=='WC',])),]
pi_alpha_ak2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$state2=='AK',])),]
pi_alpha_bc2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$state2=='BC',])),]
pi_alpha_orwa2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$state2=='OR-WA',])),]

soc_alpha2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species2=='Sockeye',])),]
soc_alpha_goa2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$ocean.basin=='GOA',])),]
soc_alpha_bs2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$ocean.basin=='BS',])),]
soc_alpha_wc2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$ocean.basin=='WC',])),]
soc_alpha_ak2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$state2=='AK',])),]
soc_alpha_bc2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$state2=='BC',])),]
soc_alpha_orwa2=prod_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$state2=='OR-WA',])),]

#umsy matrices
chi_alpha=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook',])),]
chi_alpha_goa=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$ocean.basin=='GOA',])),]
chi_alpha_bs=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$ocean.basin=='BS',])),]
chi_alpha_wc=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$ocean.basin=='WC',])),]
chi_alpha_ak=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$state2=='AK',])),]
chi_alpha_bc=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$state2=='BC',])),]
chi_alpha_orwa=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$state2=='OR-WA',])),]

chu_alpha=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum',])),]
chu_alpha_goa=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$ocean.basin=='GOA',])),]
chu_alpha_bs=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$ocean.basin=='BS',])),]
chu_alpha_wc=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$ocean.basin=='WC',])),]
chu_alpha_ak=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$state2=='AK',])),]
chu_alpha_bc=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$state2=='BC',])),]
chu_alpha_orwa=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$state2=='OR-WA',])),]

coh_alpha=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho',])),]
coh_alpha_goa=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$ocean.basin=='GOA',])),]
coh_alpha_bs=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$ocean.basin=='BS',])),]
coh_alpha_wc=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$ocean.basin=='WC',])),]
coh_alpha_ak=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$state2=='AK',])),]
coh_alpha_bc=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$state2=='BC',])),]
coh_alpha_orwa=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$state2=='OR-WA',])),]

pi_alpha=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species2=='Pink',])),]
pi_alpha_goa=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$ocean.basin=='GOA',])),]
pi_alpha_bs=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$ocean.basin=='BS',])),]
pi_alpha_wc=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$ocean.basin=='WC',])),]
pi_alpha_ak=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$state2=='AK',])),]
pi_alpha_bc=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$state2=='BC',])),]
pi_alpha_orwa=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$state2=='OR-WA',])),]

soc_alpha=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species2=='Sockeye',])),]
soc_alpha_goa=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$ocean.basin=='GOA',])),]
soc_alpha_bs=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$ocean.basin=='BS',])),]
soc_alpha_wc=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$ocean.basin=='WC',])),]
soc_alpha_ak=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$state2=='AK',])),]
soc_alpha_bc=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$state2=='BC',])),]
soc_alpha_orwa=umsy_matr[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$state2=='OR-WA',])),]

#umsy scaled matrices
chi_alpha2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook',])),]
chi_alpha_goa2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$ocean.basin=='GOA',])),]
chi_alpha_bs2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$ocean.basin=='BS',])),]
chi_alpha_wc2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$ocean.basin=='WC',])),]
chi_alpha_ak2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$state2=='AK',])),]
chi_alpha_bc2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$state2=='BC',])),]
chi_alpha_orwa2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chinook'&stock_info_filtered$state2=='OR-WA',])),]

chu_alpha2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum',])),]
chu_alpha_goa2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$ocean.basin=='GOA',])),]
chu_alpha_bs2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$ocean.basin=='BS',])),]
chu_alpha_wc2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$ocean.basin=='WC',])),]
chu_alpha_ak2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$state2=='AK',])),]
chu_alpha_bc2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$state2=='BC',])),]
chu_alpha_orwa2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Chum'&stock_info_filtered$state2=='OR-WA',])),]

coh_alpha2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho',])),]
coh_alpha_goa2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$ocean.basin=='GOA',])),]
coh_alpha_bs2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$ocean.basin=='BS',])),]
coh_alpha_wc2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$ocean.basin=='WC',])),]
coh_alpha_ak2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$state2=='AK',])),]
coh_alpha_bc2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$state2=='BC',])),]
coh_alpha_orwa2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Coho'&stock_info_filtered$state2=='OR-WA',])),]

pi_alpha2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species2=='Pink',])),]
pi_alpha_goa2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$ocean.basin=='GOA',])),]
pi_alpha_bs2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$ocean.basin=='BS',])),]
pi_alpha_wc2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$ocean.basin=='WC',])),]
pi_alpha_ak2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$state2=='AK',])),]
pi_alpha_bc2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$state2=='BC',])),]
pi_alpha_orwa2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Pink'&stock_info_filtered$state2=='OR-WA',])),]

soc_alpha2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species2=='Sockeye',])),]
soc_alpha_goa2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$ocean.basin=='GOA',])),]
soc_alpha_bs2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$ocean.basin=='BS',])),]
soc_alpha_wc2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$ocean.basin=='WC',])),]
soc_alpha_ak2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$state2=='AK',])),]
soc_alpha_bc2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$state2=='BC',])),]
soc_alpha_orwa2=umsy_matr2[as.numeric(rownames(stock_info_filtered[stock_info_filtered$species=='Sockeye'&stock_info_filtered$state2=='OR-WA',])),]

