make_design_matrix=function(x,grp){
    x2=matrix(nrow=length(x),ncol=length(unique(grp)))
    for(i in 1:length(unique(grp))){
      x2[,i]=ifelse(grp==levels(factor(grp))[i],1,0)*x
    }
  return(x2)
}

rag_n=function(x){
  rle=rle(x) #get run lengths of strings
  start_n=c(1,1+cumsum(rle$lengths)) #starting points
  start_n=start_n[-length(start_n)] #drop last value
  end_n=c(start_n[-1]-1) #end points
  end_n[length(start_n)]=length(x) #final obs.
  
  return(cbind(start_n,end_n))
}

umsyCalc <- function(a) {
  #gives the min Ricker log-likelihood
  umsy <- (1 - gsl::lambert_W0(exp(1 - a)))
  
  return(umsy)
}

corr_heatmap=function(x,info,init=FALSE,k=NA,title,cols='pal1'){
  L=max(info$end)-min(info$begin)+1
  
  logas=x$draws(variables='log_a_t',format='draws_matrix')
  logas2=apply(logas,2,median)
  logas.mat=matrix(logas2,nrow=L,ncol=nrow(info))
 
#  for(i in 1:ncol(logas.mat)){
#    logas.mat[!(seq_along(logas.mat[,i]) %in% seq(info$begin[i]-min(info$begin)+1,info$end[i]-min(info$begin)+1)),i]=NA
#  }
  
  corr=cor(logas.mat,method='pearson',use='pairwise.complete.obs')
#  corr[is.na(corr)]=0
#  diag(corr)=0
  rownames(corr)=gsub(paste('-',unique(info$species),sep=''),'',info$stock.name)
  rownames(corr)=gsub("_",' ',rownames(corr))

  #initial fit - no specified group number
  e1=pheatmap::pheatmap(corr,clustering_distance_rows='correlation',clustering_distance_cols='correlation',color=RColorBrewer::brewer.pal(n = 9, name =
                                                                                                                                                    "PuBu"), main = title,display_numbers=T,angle_col =0,fontsize=5)
  abline(c(0,1),col='darkred')
  if(init==TRUE){
    return(e1)
  }
  if(init==FALSE){
    if(is.na(k)==TRUE){print('must provide integer k, number of groups, if unknown use init==TRUE to see hiarchical clustering')}else{
      cls=cutree(e1$tree_row, k = k)
      plot(e1$tree_row)
      clsts=rect.hclust(e1$tree_row,k=k)
      clsts2=unlist(clsts)
      cls=rep(1:length(clsts), lengths(clsts))
      names(cls)=names(clsts2)
      colnames(corr)=cls[match(rownames(corr),names(cls))]
    }
    
    e=pheatmap::pheatmap(corr,clustering_distance_rows='correlation',clustering_distance_cols='correlation',color=RColorBrewer::brewer.pal(n = 9, name = "PuBu"), main = title,display_numbers=T,angle_col =0,fontsize=5)
    
    if(cols=='pal1'){
      grp_cols=list()
      grp_cols[[1]]=c('#000647','#ffa600')
      grp_cols[[2]]=c('#000647','#b91153','#ffa600')
      grp_cols[[3]]=c('#000647','#87005b','#e04144','#ffa600')
      grp_cols[[4]]=c('#000647','#6b005a','#b91153','#ee5939','#ffa600')
      grp_cols[[5]]=c('#000647','#590059','#9c0059','#d22e4b','#f46832','#ffa600')
      grp_cols[[6]]=c('#000647','#4d0257','#87005b','#b91153','#e04144','#f8722c','#ffa600')
      grp_cols[[7]]=c('#000647','#440355','#77005b','#a50058','#cb264e','#e94f3e','#fa7928','#ffa600')
    }
    if(k==2){
      cols=grp_cols[[1]][as.numeric(e$gtable$grobs[[5]]$label)]
    }
    if(k==3){
      cols=grp_cols[[2]][as.numeric(e$gtable$grobs[[5]]$label)]
    }
    if(k==4){
      cols=grp_cols[[3]][as.numeric(e$gtable$grobs[[5]]$label)]
    }
    if(k==5){
      cols=grp_cols[[4]][as.numeric(e$gtable$grobs[[5]]$label)]
    }
    if(k==6){
      cols=grp_cols[[5]][as.numeric(e$gtable$grobs[[5]]$label)]
    }
    if(k==7){
      cols=grp_cols[[6]][as.numeric(e$gtable$grobs[[5]]$label)]
    }
    if(k==8){
      cols=grp_cols[[7]][as.numeric(e$gtable$grobs[[5]]$label)]
    }
    e$gtable$grobs[[5]]$gp=grid::gpar(col=cols)
    e$gtable$grobs[[6]]$gp=grid::gpar(col=cols)
    return(e)
  }
}

scaler=function(x){return(c((x - mean(x)) / sd(x)))}

prod_ts_plot=function(x,info,title,subset.rows=seq(1:nrow(info)),log=T,umsy=F,scaled=F,ylim=c(-1,4),pdf=F,filename=NA,dim=c(6,8),mu.col='black'){
  loga=x$draws(variables='log_a_t',format='draws_matrix')
  mu_loga=x$draws(variables='mu_log_a',format='draws_matrix')
  L=c(max(info$end)-min(info$begin)+1)

  if(length(subset.rows)<nrow(info)){
    info=info[subset.rows,]
    mu_loga=matrix(nrow=nrow(loga),ncol=L)
    logat=loga
    for(t in 1:L){
      colnames(logat)=gsub("]", ".", colnames(logat))
      colnames(logat)=gsub("log_a_t[", "log_a_t.", colnames(logat),fixed=T)
      
      logas=loga[,grepl(paste('log_a_t.',t,',',sep=''),colnames(loga))]
      mu_loga[,t]=apply(logas,1,mean)
    }    
  }
  
  par(mfrow=c(1,1))
  if(log==F){
    loga=exp(loga)
    mu_loga=exp(mu_loga)
  }
  
  if(pdf==T){
    pdf(paste(filename,'_mvprodseries.pdf',sep=''),width=dim[2],height=dim[1])
  }
  if(umsy==T){
    dumsy=x$draws(variables='Umsy',format='draws_matrix')
    mu_umsy=apply(mu_loga,1:2,umsyCalc)
   
    plot(c(0,1)~c(min(info$begin),max(info$end)),bty='l',type='n',ylab='Umsy',xlab='Brood year',main=title)
    for(i in 1:nrow(info)){
      dumsy_j=dumsy[,grepl(paste(',',subset.rows[i],']',sep=''),colnames(dumsy))]
      lines(apply(dumsy_j,2,median)~seq(min(info$begin),max(info$end)),col=adjustcolor('darkgray',alpha.f=0.5),lty=5,lwd=0.8)
      lines(apply(dumsy_j[,seq(info$begin[i]-min(info$begin)+1,info$end[i]-min(info$begin)+1)],2,median)~seq(info$begin[i],info$end[i]),col=adjustcolor('darkgray',alpha.f=0.5))
    }
    lines(apply(mu_umsy,2,median)~seq(min(info$begin),max(info$end)),lwd=3,col=adjustcolor(mu.col,alpha.f=0.9))
    lines(apply(mu_umsy,2,quantile,0.025)~seq(min(info$begin),max(info$end)),lty=5,lwd=2,col=adjustcolor(mu.col,alpha.f=0.5))
    lines(apply(mu_umsy,2,quantile,0.975)~seq(min(info$begin),max(info$end)),lty=5,lwd=2,col=adjustcolor(mu.col,alpha.f=0.5))
  }else{
  if(scaled==T){
    plot(c(0,ylim[2])~c(min(info$begin),max(info$end)),bty='l',type='n',ylab='Productivity - scaled SD from median',xlab='Brood year',main=title)
    abline(h=0,lwd=0.5)
  }else{
    plot(c(ylim[1],ylim[2])~c(min(info$begin),max(info$end)),bty='l',type='n',ylab=if(log==T){expression(paste("Productivity - log(", alpha, ")"))}else{expression(paste("Productivity, max. R/S "))},xlab='Brood year',main=title)
    }
  for(i in 1:nrow(info)){
    loga_j=loga[,grepl(paste(',',subset.rows[i],']',sep=''),colnames(loga))]
    if(scaled==T){
      
      
      loga_js=apply(loga_j,1,scaler) 
      lines(apply(loga_js,1,median)~seq(min(info$begin),max(info$end)),col=adjustcolor('darkgray',alpha.f=0.5),lty=5,lwd=0.8)
      lines(apply(loga_js[seq(info$begin[i]-min(info$begin)+1,info$end[i]-min(info$begin)+1),],1,median)~seq(info$begin[i],info$end[i]),col=adjustcolor('darkgray',alpha.f=0.5))
      
    }else{
      lines(apply(loga_j,2,median)~seq(min(info$begin),max(info$end)),col=adjustcolor('darkgray',alpha.f=0.5),lty=5,lwd=0.8)
      lines(apply(loga_j[,seq(info$begin[i]-min(info$begin)+1,info$end[i]-min(info$begin)+1)],2,median)~seq(info$begin[i],info$end[i]),col=adjustcolor('darkgray',alpha.f=0.5))
    }
  }
  if(scaled==T){
    mu_logas=apply(mu_loga,1,scaler)
    lines(apply(mu_logas,1,median)~seq(min(info$begin),max(info$end)),lwd=3,col=adjustcolor(mu.col,alpha.f=0.9))
    lines(apply(mu_logas,1,quantile,0.025)~seq(min(info$begin),max(info$end)),lty=5,lwd=2,col=adjustcolor(mu.col,alpha.f=0.5))
    lines(apply(mu_logas,1,quantile,0.975)~seq(min(info$begin),max(info$end)),lty=5,lwd=2,col=adjustcolor(mu.col,alpha.f=0.5))
  }else{
    lines(apply(mu_loga,2,median)~seq(min(info$begin),max(info$end)),lwd=3,col=adjustcolor(mu.col,alpha.f=0.9))
    lines(apply(mu_loga,2,quantile,0.025)~seq(min(info$begin),max(info$end)),lty=5,lwd=2,col=adjustcolor(mu.col,alpha.f=0.5))
    lines(apply(mu_loga,2,quantile,0.975)~seq(min(info$begin),max(info$end)),lty=5,lwd=2,col=adjustcolor(mu.col,alpha.f=0.5))
  }
  }
  if(pdf==T){dev.off()}
}


mv_prod_timeseries_clust=function(x,info,heatmap,k,ylim=c(-1,4)){
  clusters=data.frame(clst=heatmap$gtable$grobs[[5]]$label,stk=heatmap$gtable$grobs[[6]]$label)
  clusters$mat.id=match(clusters$stk,gsub(paste('-',unique(info$species),sep=''),'',info$stock.name))
  plot(heatmap$tree_row)
  r=rect.hclust(heatmap$tree_row, k = k) #selects K groups based on hierarchical clustering
  length_r=unlist(lapply(r,length)) #number of stocks in each group
  
  loga=x$draws(variables='log_a_t',format='draws_matrix')
  
  L=max(info$end)-min(info$begin)+1
  
  plot(c(ylim[1],ylim[2])~c(min(info$begin),max(info$end)),bty='l',type='n',ylab='log(alpha)',xlab='Brood year',main='West Coast (Pink - Even)')
  mu_log_a=list()
  log_a_agg=list()
  prod_mat=list()
  for(i in k){
    id=r[[k]][1]
    log_a_agg[[k]]=loga[,grepl(paste(',',r[[k]][1],']',sep=''),colnames(loga))]  
    prod_mat[[k]]=matrix(nrow=length_r[k],ncol=L)
    rownames(prod_mat[[k]])=clusters$stk[cluswters$clst==k]
    
    prod_mat[[k]][1,]=apply(log_a_agg[[k]],2,median)
    for(n in 2:length_r[k]){
      id=r[[k]][n]
      log_a_j=loga[,grepl(paste(',',r[[k]][n],']',sep=''),colnames(loga))]
      log_a_agg[[k]]=cbind(log_a_agg,log_a_j)
      prod_mat[[k]][n,]=apply(log_a_j,2,median)
      lines(prod_mat[[k]][n,]~seq(min(info$begin),max(info$end)),col=adjustcolor('darkgray',alpha.f=0.5))
    }
    
    colnames(loga)=paste('log_a.',expand.grid(seq(1,L),seq(1,nrow(info)))[,1],',',expand.grid(seq(1,L),seq(1,nrow(info)))[,2],'.',sep='')
    mu_log_a[[k]]=matrix(nrow=nrow(loga),ncol=L)
    for(l in 1:L){
      log_a_t=log_a_agg[[k]][,grepl(paste('log_a.',l,',',sep=''),colnames(log_a_agg[[k]]))]  
      mu_log_a[[k]][,l]=apply(log_a_t,1,mean)
    }
    
    lines(apply(mu_log_a[[k]],2,median)~seq(min(info$begin),max(info$end)),lwd=3)
    lines(apply(mu_log_a[[k]],2,quantile,0.025)~seq(min(info$begin),max(info$end)),lwd=2,lty=3)
    lines(apply(mu_log_a[[k]],2,quantile,0.975)~seq(min(info$begin),max(info$end)),lwd=2,lty=3)
  }
}


prod_change_plot=function(x,years=5,info,title,trend=F,umsy=F,pdf=F,filename=NA,dim=c(5,8)){
  if(umsy==F){
    loga=x$draws(variables='log_a_t',format='draws_matrix')
    loga=exp(loga)
    d=matrix(nrow=nrow(info),ncol=7)
    for(i in 1:nrow(info)){
      loga_j=loga[,grepl(paste(',',i,']',sep=''),colnames(loga))]
      #remove interpolated productivity - just to observations in time-series:
      loga_j=loga_j[,seq(info$begin[i]-min(info$begin)+1,info$end[i]-min(info$begin)+1)]
      loga_j_b=apply(loga_j[,1:years],1,mean)
      d[i,1]=median(loga_j_b);d[i,2]=quantile(loga_j_b,0.025);d[i,3]=quantile(loga_j_b,0.975)
      loga_j_f=apply(loga_j[,c(ncol(loga_j)-years):ncol(loga_j)],1,mean)
      d[i,4]=median(loga_j_f);d[i,5]=quantile(loga_j_f,0.025);d[i,6]=quantile(loga_j_f,0.975)
      d[i,7]=(d[i,4]/d[i,1]-1)*100        }  
  }else{
    dumsy=x$draws(variables='Umsy',format='draws_matrix')
    d=matrix(nrow=nrow(info),ncol=7)
    for(i in 1:nrow(info)){
      dumsy_j=dumsy[,grepl(paste(',',i,']',sep=''),colnames(dumsy))]
      #remove interpolated productivity - just to observations in time-series:
      dumsy_j=dumsy_j[,seq(info$begin[i]-min(info$begin)+1,info$end[i]-min(info$begin)+1)]
      dumsy_j_b=apply(dumsy_j[,1:years],1,mean)
      d[i,1]=median(dumsy_j_b);d[i,2]=quantile(dumsy_j_b,0.025);d[i,3]=quantile(dumsy_j_b,0.975)
      dumsy_f=apply(dumsy_j[,c(ncol(dumsy_j)-years):ncol(dumsy_j)],1,mean)
      d[i,4]=median(dumsy_f);d[i,5]=quantile(dumsy_f,0.025);d[i,6]=quantile(dumsy_f,0.975)
      #Umsy below zero set to 0
      d[i,4]=ifelse(d[i,4]<0,0.001,d[i,4])
      d[i,1]=ifelse(d[i,1]<0,0.001,d[i,1])
      d[i,7]=(d[i,4]/d[i,1]-1)*100        }  
  }
  if(pdf==T){
    pdf(paste(filename,'_changeplot.pdf',sep=''),width=dim[2],height=dim[1])
  }
par(mfrow=c(1,2))
plot(d[,4]~d[,1],type='n',bty='l',ylim=if(umsy==T){c(0,1)}else{c(0,max(d[,1:2])*1.2)},xlim=if(umsy==T){c(0,1)}else{c(0,max(d[,1:2])*1.2)},ylab=if(umsy==T){paste('Recent (last ',years,' years) Umsy',sep='')}else{paste('Recent (last ',years,' years) productivity (R/S)',sep='')},xlab=if(umsy==T){paste('Initial (first ',years,' years) Umsy',sep='')}else{paste('Initial (first ',years,' years) productivity (R/S)',sep='')},cex.axis=0.8,cex.lab=0.8)  
abline(c(0,1))
p_cols1=RColorBrewer::brewer.pal(n=10,name='RdYlBu')
cuts=cut(d[,7], breaks = c(seq(-100,100,length.out=9),Inf))
p_cols2=p_cols1[as.numeric(cuts)]

for(i in 1:nrow(d)){
  lines(c(d[i,5],d[i,6])~rep(d[i,1],2),col=adjustcolor('darkgray',alpha.f = 0.15))
  lines(rep(d[i,4],2)~c(d[i,2],d[i,3]),col=adjustcolor('darkgray',alpha.f = 0.15))
  lines(c(d[i,4],d[i,1])~rep(d[i,1],2),col=adjustcolor(p_cols2[i],alpha.f = 0.5))
}
points(d[,4]~d[,1],bg=adjustcolor(p_cols2,alpha.f=0.9),col=adjustcolor('black',alpha.f=0.8),pch=21,cex=1.15)
break_cols=p_cols1[cut(c(seq(-99,100,by=10),Inf),breaks = c(seq(-100,100,length.out=9),Inf))]
break_cols[22:23]=rep(break_cols[21],2)
dcut=cut(d[,7],c(seq(-100,100,by=10),Inf))
num.dcut=as.numeric(dcut)
num.dcut=ifelse(num.dcut==21,23,num.dcut)
hist(num.dcut,breaks=seq(1:max(num.dcut)),xlab='Median % change',main='',col=adjustcolor(break_cols,alpha.f=0.9),xaxt='n',cex.axis=0.8,cex.lab=0.8,xlim=c(0,23))
axis(side=1,at=c(1,5,10,15,20,23),labels=c(seq(-100,100,by=50),'>100'),cex.axis=0.7,cex.lab=0.7)
abline(v=10)
mtext(title,side=3,outer=T,line=-3)
if(pdf==T){
 dev.off()
}
}
