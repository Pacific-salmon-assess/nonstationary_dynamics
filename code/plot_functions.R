#functions
library(wesanderson);library(RColorBrewer);library(cowplot)
#Colour assignment function - x - vector of values (broodyear here), cols = palette of choice, breaks - set the division ie. 1/10 (0.1); 1/5 (0.2)


whichColour<-function(x,cols,breaks){
  b=quantile(x,probs=seq(0,1,by=breaks)) 
  i<-1
  while(x>=b[i]&&x>b[i+1]) i<-i+1
  cols[i]
}

assignCol<- function(x,cols,breaks){
  sapply(x$broodyear,whichColour,cols=cols,breaks=breaks)
}

#Plot function - recruit per spawner & SR plot by time
plot_SR=function(x,m,path){ #x = dataset, m = S-R model fit, path = desired output folder
  x_new<- seq(min(x$spawners),max(x$spawners))
  pred_r<- exp(m$coefficients[1]+m$coefficients[2]*x_new)*x_new
  a<- exp(m$coefficients[1]);b=-m$coefficients[2]
  rmax<- a/b*exp(-1)

  cols=wes_palette("Zissou1",length(seq(0,1,by=0.02)), type="continuous")
  
  br=quantile(x$broodyear,probs=seq(0,1,by=0.02))
  col_points1=cols[findInterval(x$broodyear,br)]
  col_points2=adjustcolor(cols[findInterval(x$broodyear,br)],alpha.f=0.8)

 pdf(file.path(path,paste(paste(x$species,gsub('/','_',x$stock),x$stock.id,sep='_'),'.pdf',sep='')),width=14,height=7)
 par(mfrow=c(1,2))
 plot(exp(logR_S)~broodyear,data=x,bty='l',ylab='Recruits per Spawner',xlab='Year',type='n',cex.axis=1.2,cex.lab=1.2)
 lines(exp(logR_S)~broodyear,lwd=2,data=x)
 points(exp(logR_S)~broodyear,data=x,bg=col_points1,cex=1.7,pch=21,col='white')
 mtext(paste(unique(x$stock),unique(x$species),sep=' '),side=3,font=2,adj=1.5,line=2,cex=1.5)
 
  plot(recruits~spawners,data=x,bty='l',ylab='Recruits (R)',xlab='Spawners (S)',type='n',cex.axis=1.2,cex.lab=1.2)
  lines(pred_r~x_new,lwd=2)
  lines(c(0,rmax)~c(1/b,1/b),lty=5)
  points(recruits~spawners,data=x,bg=col_points2,cex=1.7,pch=21,col='white')
  
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
  mtext(round(seq(min(x$broodyear),max(x$broodyear),l=5)), 4, -.5, at=lbsq, las=2, cex=.6,font=2)                      ## tick labels
  mtext('Year', 3, 0.125, cex=.8, adj=-0.1, font=2)              ## title
  ##
  dev.off()  ## reset par

}


#plot acf
plot_acf<- function(x){
  plot(c(1:nrow(x))~rep(1,nrow(x)),xlim=c(-1,1),xlab='Autocorrelation [t-1]',ylab='',bty='l',type='n',yaxt='n',main=unique(x$species))
  abline(v=0,lty=5)
  for(i in 1:nrow(x)){
    lines(rep(i,2)~c(x$ar1.l95[i],x$ar1.u95[i]),lwd=2)
    points(i~x$ar1[i],pch=21,cex=1.2,bg='black')
  }
}

#plot Stan output of SR
plot_SR=function(x,m,params){ #x = dataset, m = S-R model fit, path = desired output folder
  x_new<- seq(min(x$spawners),max(x$spawners))
  pred_r<- matrix(nrow = length(x_new),ncol=length(params$a))
  for(i in 1:length(params$a)){
    pred_r[,i]<-exp(params$a[i]*x_new*exp(-params$b[i]*x_new))
  }
  a<- exp(m$coefficients[1]);b=-m$coefficients[2]
  rmax<- a/b*exp(-1)
  
  cols=wes_palette("Zissou1",length(seq(0,1,by=0.02)), type="continuous")
  
  br=quantile(x$broodyear,probs=seq(0,1,by=0.02))
  col_points1=cols[findInterval(x$broodyear,br)]
  col_points2=adjustcolor(cols[findInterval(x$broodyear,br)],alpha.f=0.8)
  
  pdf(file.path(path,paste(paste(x$species,gsub('/','_',x$stock),x$stock.id,sep='_'),'.pdf',sep='')),width=14,height=7)
  par(mfrow=c(1,2))
  plot(exp(logR_S)~broodyear,data=x,bty='l',ylab='Recruits per Spawner',xlab='Year',type='n',cex.axis=1.2,cex.lab=1.2)
  lines(exp(logR_S)~broodyear,lwd=2,data=x)
  points(exp(logR_S)~broodyear,data=x,bg=col_points1,cex=1.7,pch=21,col='white')
  mtext(paste(unique(x$stock),unique(x$species),sep=' '),side=3,font=2,adj=1.5,line=2,cex=1.5)
  
  plot(recruits~spawners,data=x,bty='l',ylab='Recruits (R)',xlab='Spawners (S)',type='n',cex.axis=1.2,cex.lab=1.2)
  lines(pred_r~x_new,lwd=2)
  lines(c(0,rmax)~c(1/b,1/b),lty=5)
  points(recruits~spawners,data=x,bg=col_points2,cex=1.7,pch=21,col='white')
  
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
  mtext(round(seq(min(x$broodyear),max(x$broodyear),l=5)), 4, -.5, at=lbsq, las=2, cex=.6,font=2)                      ## tick labels
  mtext('Year', 3, 0.125, cex=.8, adj=-0.1, font=2)              ## title
  ##
  dev.off()  ## reset par
  
}




###Maybe nix this
mass_plot_prod=function(x,m){ # x = list of datasets for each each stock
  x_d<- do.call(rbind.data.frame, x)
  x_d<- x_d[complete.cases(x_d$recruits),]
  plot(logR_S~broodyear,data=x_d,type='n',lwd=2,bty='l',xlim=c(min(na.omit(x_d$broodyear)),max(na.omit(x_d$broodyear))),ylim=c(min(na.omit(x_d$logR_S)),max(na.omit(x_d$logR_S))))
  for(i in 1:length(x)){
    lines(logR_S~broodyear,data=x[[i]],lwd=2)
    points(logR_S~broodyear,data=x[[i]],cex=1.2,bg='black',pch=21)
  }
    
}

