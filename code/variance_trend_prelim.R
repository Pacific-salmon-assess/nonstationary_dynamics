#Preliminary batched model support for varying dynamics in sockeye
library(here);library(dplyr)
sock_dat<- read.csv(here('data','filtered datasets','sockeye_final.csv'))
sock_info<- read.csv(here('data','filtered datasets','filtered_sockeye_info.csv'))

source(here('code','plot_functions.R'))


resid_trends=data.frame(stock=sock_info$stock,m.resid=NA,b.resid=NA,b.resid.l95=NA,b.resid.u95=NA,b.resid2=NA,b.resid2.l95=NA,b.resid2.u95=NA)
for(i in 1:nrow(sock_info)){
  s<- subset(sock_dat,stock.id==sock_info$stock.id[i])
  
  resids<- resid(lm(logR_S~spawners,data=s))
  s$sq_resids<- resids^2
  resid_t<- lm(s$sq_resids~seq(1:length(unique(s$broodyear))))
  logresid_t<- lm(log(s$sq_resids)~seq(1:length(unique(s$broodyear))))
  
  resid_trends[i,2]=mean(s$sq_resids)
  resid_trends[i,3]=resid_t$coefficients[2]
  resid_trends[i,4]=confint(resid_t)[2,1]
  resid_trends[i,5]=confint(resid_t)[2,2]
  resid_trends[i,6]=logresid_t$coefficients[2]
  resid_trends[i,7]=confint(logresid_t)[2,1]
  resid_trends[i,8]=confint(logresid_t)[2,2]
  
}

plot_resid_t(resid_trends,m.col=3, l95.col = 4,u95.col=5,sp='Sockeye')
