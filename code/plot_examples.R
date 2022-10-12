rm(list=ls())
library(here);library(dplyr);library(rstan)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_aug2022.csv'))
stock_info<- read.csv(here('data','filtered datasets','all_stocks_info_aug2022.csv'))

source(here('code','functions.R'))
source(here('code','lfo_functions.R'))

#Remove stocks with less than 15 years of recruitment data
stock_info_filtered=subset(stock_info,n.years>=18) #242 stocks

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock.id)) #242

stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)

i=23 #nass meziadin
s<- subset(stock_dat2,stock.id==stock_info_filtered$stock.id[i])
s<- s[complete.cases(s$spawners),]


##baseline
plot(recruits~spawners,data=s,bty='l',pch=21,cex=1.5,bg=adjustcolor('black',alpha.f=0.7),xlab='Spawners',ylab='Recruits')
m<-lm(logR_S~spawners,data=s)
a<- exp(m$coefficients[1]);b=-m$coefficients[2]
rmax<- a/b*exp(-1)
x=seq(0,max(s$spawners),by=1000)
pred_r=x*exp(m$coefficients[1]+m$coefficients[2]*x)

plot(recruits~spawners,data=s,bty='l',ylab='Recruits (R)',xlab='Spawners (S)',type='n',cex.axis=1.5,cex.lab=1.5)
lines(pred_r~x,lwd=2)
lines(c(0,rmax)~c(1/b,1/b),lty=5)
points(recruits~spawners,data=s,bg=adjustcolor('black',alpha.f=0.7),cex=1.7,pch=21,col='white')


###tv
m3f=sr_mod(type='tv',par='a',loglik=F)

r = rstan::sampling(m3f, 
                    data = list(N=nrow(s),
                                L=max(s$broodyear)-min(s$broodyear)+1,
                                ii=s$broodyear-min(s$broodyear)+1,
                                R_S=s$logR_S,
                                S=s$spawners),
                    control = list(adapt_delta = 0.99), warmup = 200, chains = 6, iter = 700)

params=extract(r)

a_mat=params$log_a;
med_a=apply(a_mat,2,median)

plot(recruits~spawners,data=s,bty='l',ylab='Recruits (R)',xlab='Spawners (S)',type='n',cex.axis=2,cex.lab=2)
cols=wes_palette("Zissou1",length(seq(0,1,by=0.02)), type="continuous")
br=quantile(s$broodyear,probs=seq(0,1,by=0.02))
col_points1=cols[findInterval(s$broodyear,br)]
col_points2=adjustcolor(cols[findInterval(s$broodyear,br)],alpha.f=0.8)
points(recruits~spawners,data=s,bg=col_points2,cex=3,pch=21,col='white')

r_pred_l=list()
for(i in 1:length(med_a)){
   x=seq(0,max(s$spawners),by=1000)
   r_pred_l[[i]]=x*exp(a_mat[i]-median(params$b)*x)
   
 lines(r_pred_l[[i]]~x,col=col_points2[i],lwd=2)
}

#Legend##
legend_image <- as.raster(matrix(rev(wes_palette("Zissou1",length(seq(0,1,by=0.02)), type="continuous")), ncol=1))
rasterImage(legend_image, 0, 0, 1,1)
op <- par(  ## set and store par
  fig=c(0.9,0.99, 0.7, 0.95),    ## set figure region , 
  mar=c(1, 1, 1, 3),                                  ## set margins
  new=TRUE)                                ## set new for overplot w/ next plot

plot(c(0, 2), c(0, 1), type='n', axes=F, xlab='', ylab='')  ## ini plot2
rasterImage(legend_image, 0, 0, 1, 1)                       ## the gradient
lbsq <- seq.int(0, 1, l=5)                                  ## seq. for labels ## axis ticks
mtext(round(seq(min(s$broodyear),max(s$broodyear),l=5)), 4, -.5, at=lbsq, las=2, cex=.9,font=2)                      ## tick labels
mtext('Year', 3, 0.125, cex=1.2, adj=-0.1, font=2)              ## title


##a trend
a_mat=params$log_a;
med_a=apply(a_mat,2,median)

plot(exp(med_a)~s$broodyear,type='n',xlab='Cohort Year',ylab='Max productivity (R/S)',bty='l',cex.lab=2,cex.axis=2)
lines(exp(med_a)~s$broodyear)
points(exp(med_a)~s$broodyear,col='white',pch=21,bg=col_points2,cex=3)
