library(here);library(dplyr);library(ggplot2);library(cowplot)
source(here('code','util_func.R'))

#stock data ####
stk.dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation2024-10-17.csv'))
stk.info=read.csv(here('data','filtered datasets','stock_info2024-10-17.csv'))


#egregik####
i=29

s=subset(stk.dat,stk.dat$stock==stk.info$stock.name[i])
df=data.frame(S=s$spawners,logRS=s$logRS,by=s$broodyear,R=s$recruits)

cmdstanr::set_cmdstan_path(path='C:/Users/greenbergda/Documents/.cmdstan/cmdstan-2.29.2')
library(cmdstanr)
options(mc.cores=6)

file=file.path(cmdstanr::cmdstan_path(),'sr models', "m3f_ip.stan")
m.tv=cmdstanr::cmdstan_model(file)

file=file.path(cmdstanr::cmdstan_path(),'sr models', "m2f_ip.stan")
m.stv=cmdstanr::cmdstan_model(file)


dl=list(N=nrow(df),
        L=max(s$broodyear)-min(s$broodyear)+1,
        ii=s$broodyear-min(s$broodyear)+1,
        R_S=s$logRS,
        S=s$spawners,
        pSmax_mean=0.5*max(s$spawners),
        pSmax_sig=2*max(s$spawners))


rw.f.s=     m.tv$sample(data=dl,
                        chains = 6, 
                        iter_warmup = 500,
                        iter_sampling =1000,
                        refresh = 100,
                        adapt_delta = 0.9999,
                        max_treedepth = 20)

x_new=seq(0,max(df$S),length.out=200)
by_q=round(quantile(df$by,seq(0,1,by=0.1)))

mod=rstan::read_stan_csv(rw.f.s$output_files())

post=rstan::extract(mod)
pred_df=data.frame(x_new)
for(n in 1:length(by_q)){
  pred_df[,1+n]=exp(median(post$log_a[,match(by_q[n],df$by)])-median(post$b)*x_new)*x_new
}



alpha_df=data.frame(by=seq(min(df$by),max(df$by)),med=apply(post$log_a,2,median),l90=apply(post$log_a,2,quantile,0.1),u90=apply(post$log_a,2,quantile,0.9))
plot2=ggplot2::ggplot(alpha_df, aes(by,exp(med))) +
  geom_line(aes(x=by,y=exp(med)),linewidth=1.3)+
  geom_point(aes(colour = by),size=4) +
  scale_colour_viridis_c(name='Year')+
  geom_ribbon(aes(ymin =exp(l90), ymax =exp(u90)), alpha = 0.2)+
  xlab("Year") + 
  ylab(paste0("Productivity, max. R/S"))+
  theme_classic(14)+
  theme(panel.background = element_blank(),text = element_text(size = 20),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=15),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold", size=18),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
plot2

pred_df=pred_df/1e3

plot1=
  ggplot2::ggplot(df, aes(S/1e3, R/1e3)) +
  geom_path(alpha=0.075,linewidth=0.75)+
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

  geom_point(aes(colour = by),size=4) +
  scale_colour_viridis_c(name='Year')+
  xlab("Spawners (thousands)") + 
  ylab("Recruits (thousands)")+
  xlim(0, max(df$S/1e3))+
  ylim(0, max(df$R/1e3))+
  theme_classic(14)+
  theme(panel.background = element_blank(),text = element_text(size = 20),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
plot1

#title <- ggdraw() + 
#  draw_label(
#    'Sockeye - Egregik, Bristol Bay, AK',
#    fontface = 'bold',
#    x = 0.5, y = 0.5, hjust = 0.5)

legend = cowplot::get_legend(plot1)

plot_rw_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                             plot2 + theme(legend.position="none"),
                             ncol=2,nrow=1)
plot=cowplot::plot_grid(plot_rw_a,legend,rel_widths = c(3,.25))
#plot=cowplot::plot_grid(title,plot,ncol=1,rel_heights = c(0.1,1))
plot

ggsave2('BabineLate_SR_RWA.png',plot=plot,dpi=600,width=14,height=6.5)

st.f.s=     m.st$sample(data=dl,
                        chains = 6, 
                        iter_warmup = 500,
                        iter_sampling =1000,
                        refresh = 100,
                        adapt_delta = 0.9999,
                        max_treedepth = 20)

x_new=seq(0,max(df$S),length.out=200)
by_q=round(quantile(df$by,seq(0,1,by=0.1)))

mod=rstan::read_stan_csv(st.f.s$output_files())

post=rstan::extract(mod)
pred_df=data.frame(x_new,pred=exp(median(post$log_a)-median(post$b)*x_new)*x_new)

plotx=
  ggplot2::ggplot(df, aes(S/1e6, R/1e6)) +
  geom_abline(slope = 1, intercept = 0, color='darkred')+
  geom_path(alpha=0.2,linewidth=0.75)+
  geom_point(aes(colour = by),size=4) +
  scale_colour_viridis_c(name='Year')+
  xlab("Spawners (millions)") + 
  ylab("Recruits (millions)")+
  xlim(0, max(df$S/1e6))+
  ylim(0, max(df$R/1e6))+
  geom_line(data=pred_df,aes(x=x_new/1e6,y=pred/1e6),linewidth=2)+
  theme_classic(14)+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
plotx

ggsave2('egregik_SR_static.png',plot=plotx,dpi=600,width=8,height=6)


smsy_df=data.frame(by=seq(min(df$by),max(df$by)),med=apply(post$Smsy,2,median),l90=apply(post$Smsy,2,quantile,0.1),u90=apply(post$Smsy,2,quantile,0.9))
plot3=ggplot2::ggplot(smsy_df, aes(by,med)) +
  geom_line(aes(x=by,y=med),linewidth=1.3)+
  geom_point(aes(colour = by),size=4) +
  scale_colour_viridis_c(name='Year')+
  geom_ribbon(aes(ymin =l90, ymax =u90), alpha = 0.2)+
  xlab("Year") + 
  ylab(paste0("Smsy"))+
  theme_classic(14)+
  theme(panel.background = element_blank(),text = element_text(size = 20),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
plot3
ggsave2('BabineLate_smsy.png',plot=plot3,dpi=600,width=8,height=6.5)

umsy_df=data.frame(by=seq(min(df$by),max(df$by)),med=apply(post$Umsy,2,median),l90=apply(post$Umsy,2,quantile,0.1),u90=apply(post$Umsy,2,quantile,0.9))
plot4=ggplot2::ggplot(umsy_df, aes(by,med)) +
  geom_line(aes(x=by,y=med),linewidth=1.3)+
  geom_point(aes(colour = by),size=4) +
  scale_colour_viridis_c(name='Year')+
  geom_ribbon(aes(ymin =l90, ymax =u90), alpha = 0.2)+
  xlab("Year") + 
  ylab(paste0("Umsy"))+
  theme_classic(14)+
  theme(panel.background = element_blank(),text = element_text(size = 20),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
plot4

ggsave2('BabineLate_umsy.png',plot=plot4,dpi=600,width=8,height=6.5)

plot_rw_a=cowplot::plot_grid(plot3 + theme(legend.position="none"),
                             plot4 + theme(legend.position="none"),
                             ncol=2,nrow=1)
#plot=cowplot::plot_grid(plot_rw_a,legend,rel_widths = c(3,.25))

ggsave2('egregik_refpts.png',plot=plot_rw_a,dpi=600,width=14,height=6.5)


#cowichan####
i=257

s=subset(stk.dat,stk.dat$stock==stk.info$stock.name[i])
df=data.frame(S=s$spawners,logRS=s$logRS,by=s$broodyear,R=s$recruits)

dl=list(N=nrow(df),
        L=max(s$broodyear)-min(s$broodyear)+1,
        ii=s$broodyear-min(s$broodyear)+1,
        R_S=s$logRS,
        S=s$spawners,
        pSmax_mean=0.5*max(s$spawners),
        pSmax_sig=2*max(s$spawners))


rw.f.s=     m.tv$sample(data=dl,
                        chains = 6, 
                        iter_warmup = 500,
                        iter_sampling =1000,
                        refresh = 100,
                        adapt_delta = 0.9999,
                        max_treedepth = 20)

x_new=seq(0,max(df$S),length.out=200)
by_q=round(quantile(df$by,seq(0,1,by=0.1)))

mod=rstan::read_stan_csv(rw.f.s$output_files())

post=rstan::extract(mod)
pred_df=data.frame(x_new)
for(n in 1:length(by_q)){
  pred_df[,1+n]=exp(median(post$log_a[,match(by_q[n],df$by)])-median(post$b)*x_new)*x_new
}



alpha_df=data.frame(by=seq(min(df$by),max(df$by)),med=apply(post$log_a,2,median),l90=apply(post$log_a,2,quantile,0.1),u90=apply(post$log_a,2,quantile,0.9))
plot2=ggplot2::ggplot(alpha_df, aes(by,exp(med))) +
  geom_line(aes(x=by,y=exp(med)),linewidth=1.3)+
  geom_point(aes(colour = by),size=4) +
  scale_colour_viridis_c(name='Year')+
  geom_ribbon(aes(ymin =exp(l90), ymax =exp(u90)), alpha = 0.2)+
  xlab("Year") + 
  ylab(paste0("Productivity, max. R/S"))+
  theme_classic(14)+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
plot2

pred_df=pred_df/1e3

plot1=
  ggplot2::ggplot(df, aes(S/1e3, R/1e3)) +
  geom_abline(slope = 1, intercept = 0, color='darkred')+
  geom_path(alpha=0.2,linewidth=0.75)+
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
  
  geom_point(aes(colour = by),size=4) +
  scale_colour_viridis_c(name='Year')+
  xlab("Spawners (thousands)") + 
  ylab("Recruits (thousands)")+
  xlim(0, max(df$S/1e3))+
  ylim(0, max(df$R/1e3))+
  theme_classic(14)+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
plot1

title <- ggdraw() + 
  draw_label(
    'Chinook - Cowichan River, BC',
    fontface = 'bold',
    x = 0.5, y = 0.5, hjust = 0.5)

legend = cowplot::get_legend(plot1)

plot_rw_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                             plot2 + theme(legend.position="none"),
                             ncol=2,nrow=1)
plot=cowplot::plot_grid(plot_rw_a,legend,rel_widths = c(3,.25))
plot=cowplot::plot_grid(title,plot,ncol=1,rel_heights = c(0.1,1))
plot

ggsave2('Cowichan_SR_RWA.png',plot=plot,dpi=600,width=11,height=6.5)

plotx=
  ggplot2::ggplot(df, aes(S/1e3, R/1e3)) +
  geom_abline(slope = 1, intercept = 0, color='darkred')+
  geom_path(alpha=0.2,linewidth=0.75)+
  geom_point(aes(colour = by),size=4) +
  scale_colour_viridis_c(name='Year')+
  xlab("Spawners (thousands)") + 
  ylab("Recruits (thousands)")+
  xlim(0, max(df$S/1e3))+
  ylim(0, max(df$R/1e3))+
  ggtitle("Chinook - Cowichan River, BC")+
  theme_classic(14)+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
plotx

ggsave2('Cowichan_SRobs.png',plot=plotx,dpi=600,width=8,height=6)

#grays harbour####
i=100

s=subset(stk.dat,stk.dat$stock==stk.info$stock.name[i])
df=data.frame(S=s$spawners,logRS=s$logRS,by=s$broodyear,R=s$recruits)

dl=list(N=nrow(df),
        L=max(s$broodyear)-min(s$broodyear)+1,
        ii=s$broodyear-min(s$broodyear)+1,
        R_S=s$logRS,
        S=s$spawners,
        pSmax_mean=0.5*max(s$spawners),
        pSmax_sig=2*max(s$spawners))


rw.f.s=     m.tv$sample(data=dl,
                        chains = 6, 
                        iter_warmup = 500,
                        iter_sampling =1000,
                        refresh = 100,
                        adapt_delta = 0.9999,
                        max_treedepth = 20)

x_new=seq(0,max(df$S),length.out=200)
by_q=round(quantile(df$by,seq(0,1,by=0.1)))

mod=rstan::read_stan_csv(rw.f.s$output_files())

post=rstan::extract(mod)
pred_df=data.frame(x_new)
for(n in 1:length(by_q)){
  pred_df[,1+n]=exp(median(post$log_a[,match(by_q[n],df$by)])-median(post$b)*x_new)*x_new
}



alpha_df=data.frame(by=seq(min(df$by),max(df$by)),med=apply(post$log_a,2,median),l90=apply(post$log_a,2,quantile,0.1),u90=apply(post$log_a,2,quantile,0.9))
plot2=ggplot2::ggplot(alpha_df, aes(by,exp(med))) +
  geom_line(aes(x=by,y=exp(med)),linewidth=1.3)+
  geom_point(aes(colour = by),size=4) +
  scale_colour_viridis_c(name='Year')+
  geom_ribbon(aes(ymin =exp(l90), ymax =exp(u90)), alpha = 0.2)+
  xlab("Year") + 
  ylab(paste0("Productivity, max. R/S"))+
  theme_classic(14)+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
plot2

pred_df=pred_df/1e3

plot1=
  ggplot2::ggplot(df, aes(S/1e3, R/1e3)) +
  geom_abline(slope = 1, intercept = 0, color='darkred')+
  geom_path(alpha=0.2,linewidth=0.75)+
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
  
  geom_point(aes(colour = by),size=4) +
  scale_colour_viridis_c(name='Year')+
  xlab("Spawners (thousands)") + 
  ylab("Recruits (thousands)")+
  xlim(0, max(df$S/1e3))+
  ylim(0, max(df$R/1e3))+
  theme_classic(14)+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
plot1

title <- ggdraw() + 
  draw_label(
    'Chum - Hood Canal, WA',
    fontface = 'bold',
    x = 0.5, y = 0.5, hjust = 0.5)

legend = cowplot::get_legend(plot1)

plot_rw_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                             plot2 + theme(legend.position="none"),
                             ncol=2,nrow=1)
plot=cowplot::plot_grid(plot_rw_a,legend,rel_widths = c(3,.25))
plot=cowplot::plot_grid(title,plot,ncol=1,rel_heights = c(0.1,1))
plot

ggsave2('hoodcanal_SR_RWA.png',plot=plot,dpi=600,width=11,height=6.5)

plotx=
  ggplot2::ggplot(df, aes(S/1e3, R/1e3)) +
  geom_abline(slope = 1, intercept = 0, color='darkred')+
  geom_path(alpha=0.2,linewidth=0.75)+
  geom_point(aes(colour = by),size=4) +
  scale_colour_viridis_c(name='Year')+
  xlab("Spawners (thousands)") + 
  ylab("Recruits (thousands)")+
  xlim(0, max(df$S/1e3))+
  ylim(0, max(df$R/1e3))+
  ggtitle("Chum - Hood Canal, WA")+
  theme_classic(14)+
  theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(face="bold", size=12),
        axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
plotx

ggsave2('hoodcanal_SRobs.png',plot=plotx,dpi=600,width=8,height=6)



##sockeye all prod. trends####
tv.scy=readRDS('./outputs/model fits/tv_fit_sockeye.RDS')


scy.info=subset(stk.info,species=='Sockeye')
scy.dat=stk.dat[stk.dat$stock.id%in%scy.info$stock.id,]

png('socprod.png',width=8,height=6,res=600,units='in')
prod_ts_plot(tv.scy,info=scy.info,title='',mu.col='#811b0a',log=F,ylim=c(0,10))
dev.off()

scy.info2=scy.info[as.numeric(sample(nrow(scy.info),20)),]

cols=c('#01038f',
  '#492f9c',
  '#6f54a9',
  '#907bb5',
  '#aea3c1',
  '#cccccc',
  '#b7c8d0',
  '#a1c3d3',
  '#88bfd7',
  '#6abada',
  '#40b6de','#6c8f84',
  '#849e95',
  '#9bada7',
  '#b4bcb9',
  '#cccccc',
  '#acc6ae',
  '#8bbf91',
  '#65b874',
  '#33b056')

plot(c(0,25)~c(min(scy.info$begin),max(scy.info$end)),bty='l',type='n',ylab='Productivity, max. R/S',xlab='Brood year',cex.lab=1.3,cex.axis=1.3)
logat=tv.scy$draws(variable='log_a_t',format='draws_matrix')
for(i in 1:nrow(scy.info2)){
  logas=logat[,grepl(paste(',',rownames(scy.info2)[i],']',sep=''),colnames(logat))]
  lines(apply(exp(logas),2,median)~seq(min(scy.info$begin),max(scy.info$end)),col=adjustcolor(cols[i],alpha.f=0.3),lwd=2)
  lines(apply(exp(logas[,seq(scy.info2$begin[i]-min(scy.info$begin)+1,scy.info2$end[i]-min(scy.info$begin)+1)]),2,median)~seq(scy.info2$begin[i],scy.info2$end[i]),col=adjustcolor(cols[i],alpha.f=0.5),lwd=2)
} 
dev.off()


all_lat_lon = scy.info2 %>% group_by(lat,lon,species2) %>% summarize(n=n(),.groups = 'keep')
all_lat_lon$lon=ifelse(all_lat_lon$lon>0,-all_lat_lon$lon,all_lat_lon$lon)


library(dplyr)
world <- ne_countries(scale = "medium", returnclass = "sf")


socmap=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("Longitude") + ylab("Latitude") +
  geom_point(data =  scy.info2, mapping = aes(x = lon, y = lat), color = cols, size = 2.5, alpha = 0.7) +
  annotate(geom = 'point', x = -145, y = 54, size = 2.5*log(1+2)) +
  annotate(geom = 'point', x = -145, y = 53, size = 2.5*log(3+2)) + 
  annotate(geom = 'point', x = -145, y = 52, size = 2.5*log(10+2)) + 
  annotate(geom = 'point', x = -145, y = 51, size = 2.5*log(20+2)) + 
  annotate(geom = 'text', x = -142, y = 55,label='n', size = 6) +
  annotate(geom = 'text', x = -142, y = 54,label='1', size = 5) +
  annotate(geom = 'text', x = -142, y = 53,label='3', size = 5) + 
  annotate(geom = 'text', x = -142, y = 52,label='10', size = 5) + 
  annotate(geom = 'text', x = -142, y = 51,label='20', size = 5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(47, 65), expand = T) +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = 'dashed', linewidth = 0.1), panel.background = element_rect(fill = '#c0e6ea'),legend.title = element_blank())
socmap
cowplot::ggsave2('soc_map.png',plot=socmap,dpi=600,width=8,height=6)






scy.info.wc=subset(scy.info,ocean.basin=='WC')
scy.info.seak=subset(scy.info,ocean.basin=='SEAK')
scy.info.goa=subset(scy.info,ocean.basin=='GOA')
scy.info.bs=subset(scy.info,ocean.basin=='BS')

library(ggspatial);library(sf);library(rnaturalearth);library(ggplot2)

sp_cols=cbind(c('#6f9ba7',
                '#316024',
                '#295f8b',
                '#f8a64c',
                '#811b0a'),c('Chinook','Chum','Coho','Pink','Sockeye'))

basin_cols=cbind(c('#016c59',
                   '#8eae8f',
                   '#ffa600',
                   '#f8732a')
                 ,c('WC','SEAK','GOA','BS'))


mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

stk.info$species2=gsub("-.*", "", stk.info$species)

library(dplyr)
world <- ne_countries(scale = "medium", returnclass = "sf")

all_lat_lon = stk.info %>% group_by(lat,lon,species2) %>% summarize(n=n(),.groups = 'keep')
all_lat_lon$lon=ifelse(all_lat_lon$lon>0,-all_lat_lon$lon,all_lat_lon$lon)
all_lat_lon$ocean.basin=stk.info$ocean.basin[match(all_lat_lon$lat,stk.info$lat)]
ll_p1=subset(all_lat_lon,species2=='Pink')
ll_ch1=subset(all_lat_lon,species2=='Chum')
ll_coh1=subset(all_lat_lon,species2=='Coho')
ll_soc1=subset(all_lat_lon,species2=='Sockeye')
ll_chi1=subset(all_lat_lon,species2=='Chinook')

socmap=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("Longitude") + ylab("Latitude") +
  geom_point(data =  ll_soc1, mapping = aes(x = lon, y = lat), color = sp_cols[5], size = 2.5*log(ll_soc1$n+2), alpha = 0.7) +
  annotate(geom = 'point', x = -145, y = 54, size = 2.5*log(1+2)) +
  annotate(geom = 'point', x = -145, y = 53, size = 2.5*log(3+2)) + 
  annotate(geom = 'point', x = -145, y = 52, size = 2.5*log(10+2)) + 
  annotate(geom = 'point', x = -145, y = 51, size = 2.5*log(20+2)) + 
  annotate(geom = 'text', x = -142, y = 55,label='n', size = 6) +
  annotate(geom = 'text', x = -142, y = 54,label='1', size = 5) +
  annotate(geom = 'text', x = -142, y = 53,label='3', size = 5) + 
  annotate(geom = 'text', x = -142, y = 52,label='10', size = 5) + 
  annotate(geom = 'text', x = -142, y = 51,label='20', size = 5) + 
  annotation_scale(location = 'bl', width_hint = 0.5) + 
  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(47, 65), expand = T) +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = 'dashed', linewidth = 0.1), panel.background = element_rect(fill = 'white'),legend.title = element_blank())

cowplot::ggsave2('soc_map.png',plot=socmap,dpi=600,width=8,height=6)

ll_soc1=subset(all_lat_lon,species2=='Sockeye'&ocean.basin=='WC')
ll_soc2=subset(all_lat_lon,species2=='Sockeye'&ocean.basin=='SEAK')
ll_soc3=subset(all_lat_lon,species2=='Sockeye'&ocean.basin=='GOA')
ll_soc4=subset(all_lat_lon,species2=='Sockeye'&ocean.basin=='BS')

socmap=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  ll_soc4, mapping = aes(x = lon, y = lat), color = basin_cols[4], size = 2.5*log(ll_soc4$n+2), alpha = 0.7) +
  geom_point(data =  ll_soc3, mapping = aes(x = lon, y = lat), color =basin_cols[3], size = 2.5*log(ll_soc3$n+2), alpha = 0.7) +
  geom_point(data =  ll_soc2, mapping = aes(x = lon, y = lat), color = basin_cols[2], size = 2.5*log(ll_soc2$n+2), alpha = 0.7) +
   geom_point(data =  ll_soc1, mapping = aes(x = lon, y = lat), color =basin_cols[1], size = 2.5*log(ll_soc1$n+2), alpha = 0.7) +
  annotate(geom = 'point', x = -145, y = 54, size = 2.5*log(1+2)) +
  annotate(geom = 'point', x = -145, y = 53, size = 2.5*log(3+2)) + 
  annotate(geom = 'point', x = -145, y = 52, size = 2.5*log(10+2)) + 
  annotate(geom = 'point', x = -145, y = 51, size = 2.5*log(20+2)) + 
  annotate(geom = 'text', x = -142, y = 55,label='n', size = 6) +
  annotate(geom = 'text', x = -142, y = 54,label='1', size = 5) +
  annotate(geom = 'text', x = -142, y = 53,label='3', size = 5) + 
  annotate(geom = 'text', x = -142, y = 52,label='10', size = 5) + 
  annotate(geom = 'text', x = -142, y = 51,label='20', size = 5) + 
  annotate(geom = 'text', x = -134, y = 50.5, label = 'West Coast', color = basin_cols[1], size = 7) + 
  annotate(geom = 'text', x = -136.5, y = 54.5, label = 'SE AK', color = basin_cols[2], size = 7) + 
  annotate(geom = 'text', x = -146, y = 57.5, label = 'Gulf of Alaska', color = basin_cols[3], size = 7) + 
  annotate(geom = 'text', x = -164, y = 58, label = 'Bering Sea', color = basin_cols[4], size = 7) + 
  #annotation_scale(location = 'bl', width_hint = 0.5) + 
  #annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = 'dashed', linewidth = 0.1), panel.background = element_rect(fill = 'white'),legend.title = element_blank())
socmap
cowplot::ggsave2('soc_map2.png',plot=socmap,dpi=600,width=8,height=6)


loga=tv.scy$draws(variables='log_a_t',format='draws_matrix')

L=c(max(scy.info$end)-min(scy.info$begin)+1)

wc.info=subset(scy.info,ocean.basin=='WC')
mu_loga.wc=matrix(nrow=nrow(loga),ncol=L)
subset.rows=paste(',',rownames(wc.info),']',sep='')
subset.rows2=paste(subset.rows, collapse = "|")
wc.loga=loga[,grepl(subset.rows2,colnames(loga))]
logat=wc.loga
for(t in 1:L){
    colnames(logat)=gsub("]", ".", colnames(logat))
    colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
    
    logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
    mu_loga.wc[,t]=apply(logas,1,mean)
  }    

seak.info=subset(scy.info,ocean.basin=='SEAK')
mu_loga.seak=matrix(nrow=nrow(loga),ncol=L)
subset.rows=paste(',',rownames(seak.info),']',sep='')
subset.rows2=paste(subset.rows, collapse = "|")
seak.loga=loga[,grepl(subset.rows2,colnames(loga))]
logat=seak.loga
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.seak[,t]=apply(logas,1,mean)
}    

head(mu_loga.seak)

goa.info=subset(scy.info,ocean.basin=='GOA')
mu_loga.goa=matrix(nrow=nrow(loga),ncol=L)
subset.rows=paste(',',rownames(goa.info),']',sep='')
subset.rows2=paste(subset.rows, collapse = "|")
goa.loga=loga[,grepl(subset.rows2,colnames(loga))]
logat=goa.loga
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.goa[,t]=apply(logas,1,mean)
}    

head(mu_loga.goa)

bs.info=subset(scy.info,ocean.basin=='BS')
mu_loga.bs=matrix(nrow=nrow(loga),ncol=L)
subset.rows=paste(',',rownames(bs.info),']',sep='')
subset.rows2=paste(subset.rows, collapse = "|")
bs.loga=loga[,grepl(subset.rows2,colnames(loga))]
logat=bs.loga
for(t in 1:L){
  colnames(logat)=gsub("]", ".", colnames(logat))
  colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
  
  logas=logat[,grepl(paste('log_a_t.',t,',',sep=''),colnames(logat))]
  mu_loga.bs[,t]=apply(logas,1,mean)
}    

head(mu_loga.bs)

png('reg_socprod.png',width=8,height=6,res=600,units='in')
par(mar=c(5,5,1,1))
ylim=c(2,8)
plot(c(ylim[1],ylim[2])~c(min(scy.info$begin),max(scy.info$end)),bty='l',type='n',ylab=expression(paste("Mean Productivity, max. R/S")),xlab='Brood year',cex.axis=1.3,cex.lab=1.3)
lines(apply(exp(mu_loga.bs),2,median)~seq(min(scy.info$begin),max(scy.info$end)),lwd=3,col=adjustcolor(basin_cols[4],alpha.f=0.9))
lines(apply(exp(mu_loga.goa),2,median)~seq(min(scy.info$begin),max(scy.info$end)),lwd=3,col=adjustcolor(basin_cols[3],alpha.f=0.9))
lines(apply(exp(mu_loga.seak),2,median)~seq(min(scy.info$begin),max(scy.info$end)),lwd=3,col=adjustcolor(basin_cols[2],alpha.f=0.9))
lines(apply(exp(mu_loga.wc),2,median)~seq(min(scy.info$begin),max(scy.info$end)),lwd=3,col=adjustcolor(basin_cols[1],alpha.f=0.9))
dev.off()


#species trends
tv.chk=readRDS('./outputs/model fits/tv_fit_chinook.RDS')
tv.cho=readRDS('./outputs/model fits/tv_fit_coho.RDS')
tv.chm=readRDS('./outputs/model fits/tv_fit_chum.RDS')
tv.pke=readRDS('./outputs/model fits/tv_fit_pinkeven.RDS')
tv.pko=readRDS('./outputs/model fits/tv_fit_pinkodd.RDS')
tv.scy=readRDS('./outputs/model fits/tv_fit_sockeye.RDS')


mu_loga.chk=tv.chk$draws(variables='mu_log_a',format='draws_matrix')
mu_loga.cho=tv.cho$draws(variables='mu_log_a',format='draws_matrix')
mu_loga.chm=tv.chm$draws(variables='mu_log_a',format='draws_matrix')
mu_loga.pke=tv.pke$draws(variables='mu_log_a',format='draws_matrix')
mu_loga.pko=tv.pko$draws(variables='mu_log_a',format='draws_matrix')
mu_loga.scy=tv.scy$draws(variables='mu_log_a',format='draws_matrix')

scy.info=subset(stk.info,species=='Sockeye')
chk.info=subset(stk.info,species=='Chinook')
cho.info=subset(stk.info,species=='Coho')
chm.info=subset(stk.info,species=='Chum')
pke.info=subset(stk.info,species=='Pink-Even')
pko.info=subset(stk.info,species=='Pink-Odd')

png('sp_prod.png',width=8,height=6,res=600,units='in')
par(mar=c(5,5,1,1))
ylim=c(2,6.5)
plot(c(ylim[1],ylim[2])~c(min(stk.info$begin),max(stk.info$end)),bty='l',type='n',ylab=expression(paste("Mean Productivity, max. R/S")),xlab='Brood year',cex.axis=1.3,cex.lab=1.3)
lines(apply(exp(mu_loga.chk),2,median)~seq(min(chk.info$begin),max(chk.info$end)),lwd=3,col=adjustcolor(sp_cols[1],alpha.f=0.9))
lines(apply(exp(mu_loga.chm),2,median)~seq(min(chm.info$begin),max(chm.info$end)),lwd=3,col=adjustcolor(sp_cols[2],alpha.f=0.9))
lines(apply(exp(mu_loga.cho),2,median)~seq(min(cho.info$begin),max(cho.info$end)),lwd=3,col=adjustcolor(sp_cols[3],alpha.f=0.9))
lines(apply(exp(mu_loga.pke),2,median)~seq(min(pke.info$begin),max(pke.info$end)),lwd=3,col=adjustcolor(sp_cols[4],alpha.f=0.9))
lines(apply(exp(mu_loga.pko),2,median)~seq(min(pko.info$begin),max(pko.info$end)),lwd=3,col=adjustcolor(sp_cols[4],alpha.f=0.9))
lines(apply(exp(mu_loga.scy),2,median)~seq(min(scy.info$begin),max(scy.info$end)),lwd=3,col=adjustcolor(sp_cols[5],alpha.f=0.9))
dev.off()

ll_soc1=subset(all_lat_lon,species2=='Sockeye')
ll_chk1=subset(all_lat_lon,species2=='Chinook')
ll_chm1=subset(all_lat_lon,species2=='Chum')
ll_cho1=subset(all_lat_lon,species2=='Coho')
ll_pke1=subset(all_lat_lon,species2=='Pink')

map=ggplot(data = world) + 
  geom_sf(fill= 'antiquewhite') +  xlab("") + ylab("") +
  geom_point(data =  ll_soc1, mapping = aes(x = lon, y = lat), color =sp_cols[5], size = 2.5*log(ll_soc1$n+2), alpha = 0.7) +
  geom_point(data =  ll_chk1, mapping = aes(x = lon, y = lat), color = sp_cols[1], size = 2.5*log(ll_chk1$n+2), alpha = 0.7) +
  geom_point(data =  ll_chm1, mapping = aes(x = lon, y = lat), color =sp_cols[2], size = 2.5*log(ll_chm1$n+2), alpha = 0.7) +
  geom_point(data =  ll_cho1, mapping = aes(x = lon, y = lat), color = sp_cols[3], size = 2.5*log(ll_cho1$n+2), alpha = 0.7) +
  geom_point(data =  ll_pke1, mapping = aes(x = lon, y = lat), color = sp_cols[4], size = 2.5*log(ll_pke1$n+2), alpha = 0.7) +
  annotate(geom = 'point', x = -145, y = 54, size = 2.5*log(1+2)) +
  annotate(geom = 'point', x = -145, y = 53, size = 2.5*log(3+2)) + 
  annotate(geom = 'point', x = -145, y = 52, size = 2.5*log(10+2)) + 
  annotate(geom = 'point', x = -145, y = 51, size = 2.5*log(20+2)) + 
  annotate(geom = 'text', x = -142, y = 55,label='n', size = 6) +
  annotate(geom = 'text', x = -142, y = 54,label='1', size = 5) +
  annotate(geom = 'text', x = -142, y = 53,label='3', size = 5) + 
  annotate(geom = 'text', x = -142, y = 52,label='10', size = 5) + 
  annotate(geom = 'text', x = -142, y = 51,label='20', size = 5) + 
  annotate(geom = 'text', x = -155, y = 54, label = 'Chinook', color = sp_cols[1], size = 7) + 
  annotate(geom = 'text', x = -155, y = 52, label = 'Chum', color = sp_cols[2], size = 7) + 
  annotate(geom = 'text', x = -155, y = 50, label = 'Coho', color = sp_cols[3], size = 7) + 
  annotate(geom = 'text', x = -155, y = 48, label = 'Pink', color = sp_cols[4], size = 7) + 
  annotate(geom = 'text', x = -155, y = 46, label = 'Sockeye', color = sp_cols[5], size = 7) + 
#  annotation_scale(location = 'bl', width_hint = 0.5) + 
#  annotation_north_arrow(location = 'bl', which_north = 'true', pad_x = unit(0.35, 'in'), pad_y = unit(0.25, 'in'), style = north_arrow_fancy_orienteering) + 
  coord_sf(xlim = c(-170, -124), ylim = c(46, 65), expand = T) +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = 'dashed', linewidth = 0.1), panel.background = element_rect(fill = 'white'),legend.title = element_blank())
map
cowplot::ggsave2('map2.png',plot=map,dpi=600,width=8,height=6)
