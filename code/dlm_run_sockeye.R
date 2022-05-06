#Preliminary batched model support for varying dynamics in sockeye

library(here)

sock_dat<- read.csv(here('data','filtered datasets','sockeye_final.csv'))
sock_info<- read.csv(here('data','sockeye','sockeye_info.csv'))

sock_info<- subset(sock_info, Stock.ID %in% sock_dat$stock.id)

source(here("code","dlm-wrapper.R"))


AICtab<-matrix(NA, ncol=4, nrow=nrow(sock_info))

for(i in seq_len(nrow(sock_info))){
  
  s <- subset(sock_dat,stock.id==sock_info$Stock.ID[i])
  
  SRdata <- data.frame(byr=s$broodyear,
    spwn=s$spawners,
    rec=s$recruits)
  
  #Model 1 - static a & b
  simple<-fitDLM(SRdata)
  
  
  #Model 2 - tv a and static b
  avary<-fitDLM(SRdata, alpha_vary = TRUE, beta_vary = FALSE)
  
  #Model 3 - tv b and static a  
  bvary<-fitDLM(SRdata, alpha_vary = FALSE, beta_vary = TRUE)


  #Model 4 - tv b and static a  
  abvary<-fitDLM(SRdata, alpha_vary = TRUE, beta_vary = TRUE)
  
  AICtab[i,]<-c(simple$AIC,avary$AIC,bvary$AIC,abvary$AIC)

  plotDLM(abvary)
}
warnings()
table( c("static","avary","bvary","abvary")[apply(AICtab,1,which.min)])/nrow(sock_info)









  
  