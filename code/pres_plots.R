#start####
library(here);library(ggplot2);library(dplyr);library(cowplot);library(formattable);library(ggthemes);library(ggspatial);library(sf);library(rnaturalearth);library(rnaturalearthdata)
source(here('code','samEst code','stan_functions.R'))
source(here('code','samEst code','lfo_stan_functions.R'))
source(here('code','samEst code','util_functions.R'))
options(dplyr.summarise.inform = FALSE)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_aug2022.csv'))
stock_info<- read.csv(here('data','filtered datasets','all_stocks_info_aug2022.csv'))
#Remove stocks with less than 15 years of recruitment data
st_info=subset(stock_info,n.years>=18) #242 stocks
st_info$stock.name=gsub('/','_',st_info$stock.name)
st_info$stock.name=gsub('&','and',st_info$stock.name)

st_info$state=factor(st_info$state,levels=c('OR','WA','BC','AK'))
st_info$ocean.basin=factor(st_info$ocean.basin,levels=c('WC','GOA','BS'))

stock_dat2=subset(stock_dat,stock.id %in% st_info$stock.id)
length(unique(stock_dat2$stock.id)) #242

stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)

st_info$state2=as.character(st_info$state)
st_info=st_info %>% mutate(state2=ifelse(state=='OR','OR-WA',state2),
                                                   state2=ifelse(state=='WA','OR-WA',state2))

#plot aesthetics####
sp_cols=cbind(c("darkgray", "darkgreen", "darkblue","sienna","darkred"),c('Chinook','Chum','Coho','Pink','Sockeye'))

mod_col=cbind(c("azure4", "azure3", "cyan4","cadetblue3","deepskyblue4",'seagreen','olivedrab','darkgreen'),seq=1:8)

mc_col=cbind(c("azure3","deepskyblue4",'darkgreen'),c('static','dynamic','regime'))


#Model class ms criterion####
aic=read.csv(here('outputs','ms_rmd','aic_table.csv'))
bic=read.csv(here('outputs','ms_rmd','aic_table.csv'))
tmb_lfo=read.csv(here('outputs','ms_rmd','tmb_lfo_table.csv'))
lfo=read.csv(here('outputs','ms_rmd','lfo_table.csv'))
loo=read.csv(here('outputs','ms_rmd','stan_looic_table_stackweight.csv'))

aic$mod_class=NA
aic= aic %>%
  mutate(
    mod_class=ifelse(best_mod<=2,'static', mod_class),
    mod_class=ifelse(best_mod %in% c(3,4,5),'dynamic', mod_class),
    mod_class=ifelse(best_mod>=6,'regime', mod_class),
  )
aic$mod_class=factor(aic$mod_class,levels=c('static','dynamic','regime'))

tot_aic = aic %>% group_by(mod_class) %>% summarize(n=n())
tot_aic$ms=rep('mle_aic',nrow(tot_aic))
tot_aic$prop=tot_aic$n/sum(tot_aic$n)

bic$mod_class=NA
bic= bic %>%
  mutate(
    mod_class=ifelse(best_mod<=2,'static', mod_class),
    mod_class=ifelse(best_mod %in% c(3,4,5),'dynamic', mod_class),
    mod_class=ifelse(best_mod>=6,'regime', mod_class),
  )
bic$mod_class=factor(bic$mod_class,levels=c('static','dynamic','regime'))

tot_bic = bic %>% group_by(mod_class) %>% summarize(n=n())
tot_bic$ms=rep('mle_bic',nrow(tot_bic))
tot_bic$prop=tot_bic$n/sum(tot_bic$n)

tmb_lfo$mod_class=NA
tmb_lfo= tmb_lfo %>%
  mutate(
    mod_class=ifelse(best_mod<=2,'static', mod_class),
    mod_class=ifelse(best_mod %in% c(3,4,5),'dynamic', mod_class),
    mod_class=ifelse(best_mod>=6,'regime', mod_class),
  )
tmb_lfo$mod_class=factor(tmb_lfo$mod_class,levels=c('static','dynamic','regime'))

tot_tmb_lfo = tmb_lfo %>% group_by(mod_class) %>% summarize(n=n())
tot_tmb_lfo$ms=rep('mle_lfo',nrow(tot_tmb_lfo))
tot_tmb_lfo$prop=tot_tmb_lfo$n/sum(tot_tmb_lfo$n)

loo$mod_class=NA
loo= loo %>%
  mutate(
    mod_class=ifelse(best_mod<=2,'static', mod_class),
    mod_class=ifelse(best_mod %in% c(3,4,5),'dynamic', mod_class),
    mod_class=ifelse(best_mod>=6,'regime', mod_class),
  )
loo$mod_class=factor(loo$mod_class,levels=c('static','dynamic','regime'))

tot_loo = loo %>% group_by(mod_class) %>% summarize(n=n())
tot_loo$ms=rep('mcmc_loo',nrow(tot_loo))
tot_loo$prop=tot_loo$n/sum(tot_loo$n)

lfo$mod_class=NA
lfo= lfo %>%
  mutate(
    mod_class=ifelse(best_mod<=2,'static', mod_class),
    mod_class=ifelse(best_mod %in% c(3,4,5),'dynamic', mod_class),
    mod_class=ifelse(best_mod>=6,'regime', mod_class),
  )
lfo$mod_class=factor(lfo$mod_class,levels=c('static','dynamic','regime'))

tot_lfo = lfo %>% group_by(mod_class) %>% summarize(n=n())
tot_lfo$ms=rep('mcmc_lfo',nrow(tot_lfo))
tot_lfo$prop=tot_lfo$n/sum(tot_lfo$n)

tot_ms=rbind(tot_aic,tot_bic,tot_tmb_lfo,tot_loo,tot_lfo)

tot_ms$ms=factor(tot_ms$ms,levels=c('mle_aic','mle_bic','mle_lfo','mcmc_loo','mcmc_lfo'))

ggplot(tot_ms, aes(fill=factor(mod_class), y=prop, x=ms))+
  scale_fill_manual(values=mc_col[match(levels(factor(tot_ms$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("Top-ranked models")+
  geom_bar(position="stack", stat="identity",color='white') +
  geom_text(label=round(tot_ms$prop,2), position = position_stack(vjust = 0.5),colour='white',size=3)+
  ylim(0,1)+
  theme_minimal() +
  xlab('')+
  ylab('Prop. populations')+ 
  theme(text=element_text(size=14),axis.text=element_text(size=14),axis.line = element_line(colour = "black", 
                                                                                            size = 1, linetype = "solid"),plot.title = element_text(size=14))

tot_aic = aic %>% group_by(best_mod) %>% summarize(n=n())
tot_aic$ms=rep('mle_aic',nrow(tot_aic))
tot_aic$prop=tot_aic$n/sum(tot_aic$n)
tot_bic = bic %>% group_by(best_mod) %>% summarize(n=n())
tot_bic$ms=rep('mle_bic',nrow(tot_bic))
tot_bic$prop=tot_bic$n/sum(tot_bic$n)
tot_tmb_lfo = tmb_lfo %>% group_by(best_mod) %>% summarize(n=n())
tot_tmb_lfo$ms=rep('mle_lfo',nrow(tot_tmb_lfo))
tot_tmb_lfo$prop=tot_tmb_lfo$n/sum(tot_tmb_lfo$n)
tot_loo = loo %>% group_by(best_mod) %>% summarize(n=n())
tot_loo$ms=rep('mcmc_loo',nrow(tot_loo))
tot_loo$prop=tot_loo$n/sum(tot_loo$n)
tot_lfo = lfo %>% group_by(best_mod) %>% summarize(n=n())
tot_lfo$ms=rep('mcmc_lfo',nrow(tot_lfo))
tot_lfo$prop=tot_lfo$n/sum(tot_lfo$n)

tot_ms=rbind(tot_aic,tot_bic,tot_tmb_lfo,tot_loo,tot_lfo)

tot_ms$ms=factor(tot_ms$ms,levels=c('mle_aic','mle_bic','mle_lfo','mcmc_loo','mcmc_lfo'))


ggplot(tot_ms, aes(fill=factor(best_mod), y=prop, x=ms))+
  scale_fill_manual(values=mod_col[match(levels(factor(tot_ms$best_mod)),mod_col[,2])],name='Model')+ 
  ggtitle("Top-ranked models")+
  geom_bar(position="stack", stat="identity",color='white') +
  geom_text(label=round(tot_ms$prop,2), position = position_stack(vjust = 0.5),colour='white',size=3)+
  ylim(0,1)+
  theme_minimal() +
  xlab('')+
  ylab('Prop. populations')+ 
  theme(text=element_text(size=14),axis.text=element_text(size=14),axis.line = element_line(colour = "black", 
  size = 1, linetype = "solid"),plot.title = element_text(size=14))



#by species####
sp = aic %>% group_by(species,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = aic %>% group_by(species) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$species,sp1$species)]

aic_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=species))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MLE AIC")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

sp = bic %>% group_by(species,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = bic %>% group_by(species) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$species,sp1$species)]

bic_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=species))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MLE BIC")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

sp = loo %>% group_by(species,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = loo %>% group_by(species) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$species,sp1$species)]

loo_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=species))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MCMC LOO-CV")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

sp = lfo %>% group_by(species,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = lfo %>% group_by(species) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$species,sp1$species)]

lfo_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=species))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MCMC LFO-CV, L=2/3n")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

sp = tmb_lfo %>% group_by(species,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = tmb_lfo %>% group_by(species) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$species,sp1$species)]

tmb_lfo_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=species))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MLE LFO-CV, L=10")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

plot_grid(aic_pl+ theme(legend.position="none"),
              bic_pl+ theme(legend.position="none"),
              tmb_lfo_pl+ theme(legend.position="none"),
              loo_pl+ theme(legend.position="none"),
               lfo_pl + theme(legend.position="none"),
              legend= get_legend(lfo_pl),
              ncol=3,nrow=2,labels=c("A","B","C","D","E",''))

#by region####
aic$state2=as.character(aic$state)
aic=aic %>% mutate(state2=ifelse(state=='OR','OR-WA',state2),
                                                   state2=ifelse(state=='WA','OR-WA',state2))
aic$state2=factor(aic$state2,levels=c('OR-WA','BC','AK'))

bic$state2=as.character(bic$state)
bic=bic %>% mutate(state2=ifelse(state=='OR','OR-WA',state2),
                   state2=ifelse(state=='WA','OR-WA',state2))
bic$state2=factor(bic$state2,levels=c('OR-WA','BC','AK'))

tmb_lfo$state2=as.character(tmb_lfo$state)
tmb_lfo=tmb_lfo %>% mutate(state2=ifelse(state=='OR','OR-WA',state2),
                   state2=ifelse(state=='WA','OR-WA',state2))
tmb_lfo$state2=factor(tmb_lfo$state2,levels=c('OR-WA','BC','AK'))

loo$state2=as.character(loo$state)
loo=loo %>% mutate(state2=ifelse(state=='OR','OR-WA',state2),
                   state2=ifelse(state=='WA','OR-WA',state2))
loo$state2=factor(loo$state2,levels=c('OR-WA','BC','AK'))

lfo$state2=as.character(lfo$state)
lfo=lfo %>% mutate(state2=ifelse(state=='OR','OR-WA',state2),
                   state2=ifelse(state=='WA','OR-WA',state2))
lfo$state2=factor(lfo$state2,levels=c('OR-WA','BC','AK'))

sp = aic %>% group_by(state2,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = aic %>% group_by(state2) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$state2,sp1$state2)]

aic_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=state2))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MLE AIC")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

sp = bic %>% group_by(state2,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = bic %>% group_by(state2) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$state2,sp1$state2)]

bic_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=state2))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MLE BIC")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

sp = loo %>% group_by(state2,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = loo %>% group_by(state2) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$state2,sp1$state2)]

loo_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=state2))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MCMC LOO-CV")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

sp = lfo %>% group_by(state2,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = lfo %>% group_by(state2) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$state2,sp1$state2)]

lfo_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=state2))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MCMC LFO-CV, L=2/3n")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

sp = tmb_lfo %>% group_by(state2,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = tmb_lfo %>% group_by(state2) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$state2,sp1$state2)]

tmb_lfo_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=state2))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MLE LFO-CV, L = 10")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

plot_grid(aic_pl+ theme(legend.position="none"),
          bic_pl+ theme(legend.position="none"),
          tmb_lfo_pl+ theme(legend.position="none"),
          loo_pl+ theme(legend.position="none"),
          lfo_pl + theme(legend.position="none"),
          legend= get_legend(lfo_pl),
          ncol=3,nrow=2,labels=c("A","B","C","D","E",''))

#by length####
st_info$ts_length_bin=cut(
  st_info$n.years,
  breaks = quantile(st_info$n.years, seq(0,1,by=0.2)),
  labels = c('18-26','27-35','35-45','46-50','50+'))
st_info$ts_length_bin[match(197,st_info$stock.id)]='18-26'
aic$ts_length_bin=st_info$ts_length_bin[match(aic$stock.name,st_info$stock.name)]
bic$ts_length_bin=st_info$ts_length_bin[match(bic$stock.name,st_info$stock.name)]
tmb_lfo$ts_length_bin=st_info$ts_length_bin[match(tmb_lfo$stock.name,st_info$stock.name)]
loo$ts_length_bin=st_info$ts_length_bin[match(loo$stock.name,st_info$stock.name)]
lfo$ts_length_bin=st_info$ts_length_bin[match(lfo$stock.name,st_info$stock.name)]


sp = aic %>% group_by(ts_length_bin,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = aic %>% group_by(ts_length_bin) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$ts_length_bin,sp1$ts_length_bin)]

aic_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=ts_length_bin))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MLE AIC")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

sp = bic %>% group_by(ts_length_bin,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = bic %>% group_by(ts_length_bin) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$ts_length_bin,sp1$ts_length_bin)]

bic_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=ts_length_bin))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MLE BIC")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)

sp = loo %>% group_by(ts_length_bin,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = loo %>% group_by(ts_length_bin) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$ts_length_bin,sp1$ts_length_bin)]

loo_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=ts_length_bin))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MCMC LOO-CV")+
  geom_bar(position="stack", stat="identity",color='white') +
  geom_text(label=sp$n, position = position_stack(vjust = 0.5),colour='white',size=2.2)+
  ylim(0,1)+
  theme_minimal() +
  xlab('Time-Series length')+
  ylab('Prop. populations')+ 
  theme(text=element_text(size=14),axis.text=element_text(size=14),axis.line = element_line(colour = "black", 
                                                                                            size = 1, linetype = "solid"),plot.title = element_text(size=14))+
  annotation_custom(grid::textGrob(sp1$n[1]), xmin = 1, xmax = 1, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[2]), xmin = 2, xmax = 2, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[3]), xmin = 3, xmax = 3, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[4]), xmin = 4, xmax = 4, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)


sp = lfo %>% group_by(ts_length_bin,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = lfo %>% group_by(ts_length_bin) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$ts_length_bin,sp1$ts_length_bin)]

lfo_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=ts_length_bin))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MCMC LFO-CV, L=2/3n")+
  geom_bar(position="stack", stat="identity",color='white') +
  geom_text(label=sp$n, position = position_stack(vjust = 0.5),colour='white',size=2.2)+
  ylim(0,1)+
  theme_minimal() +
  xlab('Time-Series length')+
  ylab('Prop. populations')+ 
  theme(text=element_text(size=14),axis.text=element_text(size=14),axis.line = element_line(colour = "black", 
                                                                                            size = 1, linetype = "solid"),plot.title = element_text(size=14))+
  annotation_custom(grid::textGrob(sp1$n[1]), xmin = 1, xmax = 1, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[2]), xmin = 2, xmax = 2, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[3]), xmin = 3, xmax = 3, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[4]), xmin = 4, xmax = 4, ymin = 1.03, ymax = 1.03)+
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)


sp = tmb_lfo %>% group_by(ts_length_bin,mod_class) %>% summarize(n=n(),.groups='keep')
sp1 = tmb_lfo %>% group_by(ts_length_bin) %>% summarize(n=n())
sp$prop = sp$n/sp1$n[match(sp$ts_length_bin,sp1$ts_length_bin)]

tmb_lfo_pl=ggplot(sp, aes(fill=factor(mod_class), y=prop, x=ts_length_bin))+
  scale_fill_manual(values=mc_col[match(levels(factor(sp$mod_class)),mc_col[,2])],name='Model')+ 
  ggtitle("MLE LFO-CV, L = 10")+
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
  annotation_custom(grid::textGrob(sp1$n[5]), xmin = 5, xmax = 5, ymin = 1.03, ymax = 1.03)


plot_grid(aic_pl+ theme(legend.position="none"),
          bic_pl+ theme(legend.position="none"),
          tmb_lfo_pl+ theme(legend.position="none"),
          loo_pl+ theme(legend.position="none"),
          lfo_pl + theme(legend.position="none"),
          legend= get_legend(lfo_pl),
          ncol=3,nrow=2,labels=c("A","B","C","D","E",''))

##
i=24

m1f=samEst::sr_mod(type='static',ac = FALSE,par='n',lfo =F)
m3f=samEst::sr_mod(type='rw',par='a',lfo=F)

df=data.frame(by=s$broodyear,R=s$recruits,S=s$spawners)
x_new=seq(min(df$S),max(df$S),length.out=200)

f1 = rstan::sampling(m1f, 
                     data = list(N=nrow(s),
                                 L=max(s$broodyear)-min(s$broodyear)+1,
                                 ii=s$broodyear-min(s$broodyear)+1,
                                 R_S=s$logR_S,
                                 S=s$spawners),
                     control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)

f3 = rstan::sampling(m3f, 
                     data = list(N=nrow(s),
                                 L=max(s$broodyear)-min(s$broodyear)+1,
                                 ii=s$broodyear-min(s$broodyear)+1,
                                 R_S=s$logR_S,
                                 S=s$spawners),
                     control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)

post=rstan::extract(mod)
pred_df=data.frame(pred=exp(median(post$log_a)-median(post$b)*x_new)*x_new,x_new=x_new)

sr_plot(df=df,mod=f1,type='static',form='stan',title=paste(stock_info_filtered$stock.name[i]),sr_only=FALSE)

ggplot2::ggplot(df, aes(S, R)) +
  geom_line(data=pred_df,aes(x=x_new,y=pred),linewidth=1.3)+
  geom_point(aes(colour = by),size=2.5) +
  scale_colour_viridis_c(name='Year')+
  ggtitle(title)+
  xlab("Spawners") + 
  ylab("Recruits")+
  xlim(0, max(df$S))+
  ylim(0, max(df$R))+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))

x_new=seq(min(df$S),max(df$S),length.out=200)
by_q=round(quantile(df$by,seq(0,1,by=0.1)))

post=rstan::extract(mod)
pred_df=data.frame(x_new)
for(n in 1:length(by_q)){
  pred_df[,1+n]=exp(median(post$log_a[,match(by_q[n],df$by)])-median(post$b)*x_new)*x_new
}
ggplot2::ggplot(df, aes(S, R)) +
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
  geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
  geom_point(aes(colour = by),size=2.5) +
  scale_colour_viridis_c(name='Year')+
  ggtitle(title)+
  xlab("Spawners") + 
  ylab("Recruits")+
  xlim(0, max(df$S))+
  ylim(0, max(df$R))+
  theme_classic(14)+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))

#weight variance####
aic$max_w=apply(aic[,8:15],1,max)
bic$max_w=apply(bic[,8:15],1,max)
tmb_lfo$max_w=apply(tmb_lfo[,8:15],1,max)
loo$max_w=apply(loo[,8:15],1,max)
lfo$max_w=apply(lfo[,8:15],1,max)

hist(aic$max_w)
hist(bic$max_w)
hist(tmb_lfo$max_w)
hist(loo$max_w)
hist(lfo$max_w)

his_aic=ggplot(aic, aes(x=max_w)) + 
  geom_histogram(binwidth=0.05,color='white')+theme_classic(14)+
  xlab("Max. model weight")+ggtitle("AIC")
his_bic=ggplot(bic, aes(x=max_w)) + 
  geom_histogram(binwidth=0.05,color='white')+theme_classic(14)+
  xlab("Max. model weight")+ggtitle("BIC")
his_tmb_lfo=ggplot(tmb_lfo, aes(x=max_w)) + 
  geom_histogram(binwidth=0.05,color='white')+theme_classic(14)+
  xlab("Max. model weight")+ggtitle("MLE LFO-CV")
his_loo=ggplot(loo, aes(x=max_w)) + 
  geom_histogram(binwidth=0.05,color='white')+theme_classic(14)+
  xlab("Max. model weight")+ggtitle("MCMC LOO-CV")
his_lfo=ggplot(lfo, aes(x=max_w)) + 
  geom_histogram(binwidth=0.05,color='white')+theme_classic(14)+
  xlab("Max. model weight")+ggtitle("MCMC LFO-CV")


plot_grid(his_aic,
          his_bic,
          his_tmb_lfo,
          his_loo,
          his_lfo,
          ncol=3,nrow=2,labels=c("A","B","C","D","E",''))

#weight variance####
sum(aic[,8])/nrow(aic)
sum(aic[,9])/nrow(aic)
sum(aic[,10])/nrow(aic)
sum(aic[,11])/nrow(aic)
sum(aic[,12])/nrow(aic)
sum(aic[,13])/nrow(aic)
sum(aic[,14])/nrow(aic)
sum(aic[,15])/nrow(aic)

sum(bic[,8])/nrow(bic)
sum(bic[,9])/nrow(bic)
sum(bic[,10])/nrow(bic)
sum(bic[,11])/nrow(bic)
sum(bic[,12])/nrow(bic)
sum(bic[,13])/nrow(bic)
sum(bic[,14])/nrow(bic)
sum(bic[,15])/nrow(bic)

sum(loo[,8])/nrow(loo)
sum(loo[,9])/nrow(loo)
sum(loo[,10])/nrow(loo)
sum(loo[,11])/nrow(loo)
sum(loo[,12])/nrow(loo)
sum(loo[,13])/nrow(loo)
sum(loo[,14])/nrow(loo)
sum(loo[,15])/nrow(loo)

sum(lfo[,8])/nrow(lfo)
sum(lfo[,9])/nrow(lfo)
sum(lfo[,10])/nrow(lfo)
sum(lfo[,11])/nrow(lfo)
sum(lfo[,12])/nrow(lfo)
sum(lfo[,13])/nrow(lfo)
sum(lfo[,14])/nrow(lfo)
sum(lfo[,15])/nrow(lfo)
