#multipop S-R curves
#stock characteristics
library(here);library(dplyr);library(rstan);library(ggplot2)
stock_dat<- read.csv(here('data','filtered datasets','salmon_productivity_compilation_aug2022.csv'))
stock_info<- read.csv(here('data','filtered datasets','all_stocks_info_aug2022.csv'))

#source(here('code','samEst code','stan_functions.R'))
#source(here('code','samEst code','lfo_stan_functions.R'))
#source(here('code','samEst code','lfo_TMB_functions.R'))
##source(here('code','samEst code','TMB_functions.R'))
#source(here('code','samEst code','util_functions.R'))
library(samEst)
options(mc.cores = parallel::detectCores())
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')
source(here('code','util_func.R'))

###Load in data####
#Remove stocks with less than 15 years of recruitment data
stock_info_filtered=subset(stock_info,n.years>=18) #242 stocks
stock_info_filtered$stock.name=gsub('/','_',stock_info_filtered$stock.name)
stock_info_filtered$stock.name=gsub('&','and',stock_info_filtered$stock.name)

stock_dat2=subset(stock_dat,stock.id %in% stock_info_filtered$stock.id)
length(unique(stock_dat2$stock.id)) #242
stock_info_filtered$stock.id2=seq(1:nrow(stock_info_filtered))
if(any(stock_dat2$spawners==0)){stock_dat2$spawners=stock_dat2$spawners+1}
if(any(stock_dat2$recruits==0)){stock_dat2$recruits=stock_dat2$recruits+1}

stock_dat2$logR_S=log(stock_dat2$recruits/stock_dat2$spawners)
stock_dat2=stock_dat2[complete.cases(stock_dat2$logR_S),]

#Skeena/Nass stocks
skeena_stocks=subset(stock_info_filtered,lat==54.2237)
skeena_dat=subset(stock_dat2,stock.id %in% skeena_stocks$stock.id)
length(unique(skeena_dat$stock.id))

stock_year=expand.grid(unique(skeena_dat$stock),unique(skeena_dat$broodyear))
stock_year=stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")

#Find the position of each stock-year combination in all possible stock-year combinations (ie. where in the total vector should we expect this estimate)
skeena_dat$stock_yr=match(paste(skeena_dat$stock,skeena_dat$broodyear,sep='_'),stock_year[,3])

X_s=make_design_matrix(skeena_dat$spawners,grp=skeena_dat$stock)

mc="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int J_i[N];//index of stocks
  int J_ii[N];//index of stock-year combination
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  }
parameters{
  real log_a0; //inital productivity (global)
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<upper = 0>[J] log_b; // rate capacity - fixed in this

 //variance components  
 vector<lower = 0>[J] sigma;

//MVN parameters  
  cholesky_factor_corr[J] Lcorr;
  vector[L-1] z_dev; //average deviation in stock productivity among stocks
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  real<lower = 0> sigma_a; //variance in average productivity among years
  vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years
  
}
transformed parameters{
  vector[J] b; //capacity rate
  //global state process
  vector[L] log_a; //global (ie. average) productivity among stocks over time
  //stock state process
  matrix[L,J] log_a_s; //stock deviation from global log_a over time
  matrix[L-1,J] a_dev; //stock-level deviations in year-to-year productivity
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  
  b=exp(log_b);
  
  //initial productivities
  log_a[1] = log_a0; //average productivity among stocks at t = 1
  log_a_s[1,] = to_row_vector(log_a0_s); //stock deviation in productivity at t =1
  log_a_t[1,] = log_a[1]+to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  
  for(t in 1:L-1){
   a_dev[t,] = (diag_pre_multiply(sigma_a_s, Lcorr) * to_vector(z_dev_s[t,]))';
  }
  
  for(t in 2:L){
    log_a[t] = log_a[t-1] + z_dev[t-1]*sigma_a; //global prod. random walk
    log_a_s[t,] = log_a_s[t-1,] + a_dev[t-1,]; //stock-specific deviation in random walk
    log_a_t[t,] = log_a[t] + log_a_s[t,]; //final estimate of stock prod. (global + stock)
  }
  
}  
model{
  //priors
  log_a0 ~ gamma(3,1.5); //initial globalproductivity - wide prior
  log_a0_s ~ normal(0,1); //stock-level deviations around the global productivity
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  z_dev ~ std_normal(); //standardized stock-level deviances in prod
  to_vector(z_dev_s) ~ std_normal(); //global deviation in prod.
  
  Lcorr ~ lkj_corr_cholesky(2.0); // prior for correlation of process deviances
 
  //variance terms
  target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior  
  target += normal_lpdf(sigma_a_s | 0, 0.5) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior
  target += normal_lpdf(sigma_a| 0, 0.5) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior
 

  R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities {
  vector[N] log_lik;

  for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
  corr_matrix[J] Cor_t = multiply_lower_tri_self_transpose(Lcorr);
}


"

mc2="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int J_i[N];//index of stocks
  int J_ii[N];//index of stock-year combination
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  }
parameters{
  real log_a0; //inital productivity (global)
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<upper = 0>[J] log_b; // rate capacity - fixed in this

 //variance components  
 vector<lower = 0>[J] sigma;

//MVN parameters  
  vector[L-1] z_dev; //average deviation in stock productivity among stocks
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  real<lower = 0> sigma_a; //variance in average productivity among years
  vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years
  
}
transformed parameters{
  vector[J] b; //capacity rate
  //global state process
  vector[L] log_a; //global (ie. average) productivity among stocks over time
  //stock state process
  matrix[L,J] log_a_s; //stock deviation from global log_a over time
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  
  b=exp(log_b);
  
  //initial productivities
  log_a[1] = log_a0; //average productivity among stocks at t = 1
  
  log_a_s[1,] = to_row_vector(log_a0_s); //stock deviation in productivity at t =1
  log_a_t[1,] = log_a[1]+to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  for(t in 2:L){
    log_a[t] = log_a[t-1] + z_dev[t-1]*sigma_a; //global prod. random walk
    log_a_s[t,] = log_a_s[t-1,] + z_dev_s[t-1,].*to_row_vector(sigma_a_s); //stock-specific deviation in random walk
    log_a_t[t,] = log_a[t] + log_a_s[t,]; //final estimate of stock prod. (global + stock)
    }
  
}  
model{
  //priors
  log_a0 ~ gamma(3,1); //initial globalproductivity - wide prior
  log_a0_s ~ normal(0,1); //stock-level deviations around the initial global productivity
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  z_dev ~ std_normal(); //standardized stock-level deviances in prod
  to_vector(z_dev_s) ~ std_normal(); //global deviation in prod.
  
  //variance terms
  target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior  
  target += normal_lpdf(sigma_a_s | 0, 0.5) - normal_lcdf(0 | 0, 0.5); //remove density below zero to accomodate half-normal prior
  target += normal_lpdf(sigma_a| 0, 0.5) - normal_lcdf(0 | 0, 0.5); //remove density below zero to accomodate half-normal prior
 
  R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities{
vector[N] log_lik;

for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
}
"

mc3="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int J_i[N];//index of stocks
  int J_ii[N];//index of stock-year combination
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  }
parameters{
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<upper = 0>[J] log_b; // rate capacity - fixed in this

 //variance components  
 vector<lower = 0>[J] sigma;

//MVN parameters  
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years
  
}
transformed parameters{
  vector[J] b; //capacity rate
  //global state process
  //stock state process
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  
  b=exp(log_b);
  
  //initial productivities
  
  log_a_t[1,] = to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  for(t in 2:L){
    log_a_t[t,] = log_a_t[t-1,] + z_dev_s[t-1,].*to_row_vector(sigma_a_s); //stock-specific deviation in random walk
    }
}  
model{
  //priors
  log_a0_s ~ normal(3,2); //stock-level deviations around the initial global productivity
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  to_vector(z_dev_s) ~ std_normal(); //global deviation in prod.
  
  //variance terms
  target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior  
  target += normal_lpdf(sigma_a_s | 0, 0.5) - normal_lcdf(0 | 0, 0.5); //remove density below zero to accomodate half-normal prior

  R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities{
vector[N] log_lik;

for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
}
"

mc_avg="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int J_i[N];//index of stocks
  int J_ii[N];//index of stock-year combination
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  }
parameters{
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<upper = 0>[J] log_b; // rate capacity - fixed in this

 //variance components  
 vector<lower = 0>[J] sigma;

//MVN parameters  
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years
  
}
transformed parameters{
  vector[J] b; //capacity rate
  //global state process
  //stock state process
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  vector[L] log_a_g; //global (mean) productivity
  
  b=exp(log_b);
  
  //initial productivities
  
  log_a_t[1,] = to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  log_a_g[1] = mean(to_array_1d(log_a_t[1,]));
  for(t in 2:L){
    log_a_t[t,] = log_a_t[t-1,] + z_dev_s[t-1,].*to_row_vector(sigma_a_s); //stock-specific deviation in random walk
    log_a_g[t]= mean(to_array_1d(log_a_t[t,]));
  }
    
    
}  
model{
  //priors
  log_a0_s ~ normal(3,2); //stock-level deviations around the initial global productivity
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  to_vector(z_dev_s) ~ std_normal(); //global deviation in prod.
  
  //variance terms
  target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior  
  target += normal_lpdf(sigma_a_s | 0, 0.5) - normal_lcdf(0 | 0, 0.5); //remove density below zero to accomodate half-normal prior

  R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities{
vector[N] log_lik;

for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
}
"

multi_avg="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int J_i[N];//index of stocks
  int J_ii[N];//index of stock-year combination
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  }
parameters{
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<upper = 0>[J] log_b; // rate capacity - fixed in this

 //variance components  
 vector<lower = 0>[J] sigma;

 
//MVN parameters  
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  cholesky_factor_corr[J] Lcorr;
  vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years
  
}
transformed parameters{
  vector[J] b; //capacity rate
  //global state process
  //stock state process
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  matrix[L-1,J] a_dev; //stock-level deviations in year-to-year productivity
  vector[L] log_a_g; //global (mean) productivity
  
  b=exp(log_b);
  
  //initial productivities
  
  log_a_t[1,] = to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  log_a_g[1] = mean(to_array_1d(log_a_t[1,]));
  
    
  for(t in 1:L-1){
   a_dev[t,] = (diag_pre_multiply(sigma_a_s, Lcorr) * to_vector(z_dev_s[t,]))';
  }
  
  
  for(t in 2:L){
    log_a_t[t,] = log_a_t[t-1,] + a_dev[t-1,]; //stock-specific deviation in random walk
    log_a_g[t]= mean(to_array_1d(log_a_t[t,]));
  }
    
    
}  
model{
  //priors
  log_a0_s ~ normal(3,2); //stock-level deviations around the initial global productivity
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  to_vector(z_dev_s) ~ std_normal(); //global deviation in prod.
  
  Lcorr ~ lkj_corr_cholesky(2.0); // prior for correlation of process deviances
 
  //variance terms
  target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior  
  target += normal_lpdf(sigma_a_s | 0, 0.5) - normal_lcdf(0 | 0, 0.5); //remove density below zero to accomodate half-normal prior

  R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities{
vector[N] log_lik;

for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
corr_matrix[J] Cor_t = multiply_lower_tri_self_transpose(Lcorr);
}
"

#Multi-var with global avg
mv_avg_m=rstan::stan_model(model_code=multi_avg)

test_mv = rstan::sampling(mv_avg_m, 
                     data = list(N=nrow(skeena_dat),
                                 L=max(skeena_dat$broodyear)-min(skeena_dat$broodyear)+1,
                                 J_i=as.numeric(factor(skeena_dat$stock)),
                                 J_ii=skeena_dat$stock_yr,
                                 J=length(unique(skeena_dat$stock)),
                                 R_S=skeena_dat$logR_S,
                                 S=X_s),
                     control = list(adapt_delta = 0.99,max_treedepth=20), warmup = 100, chains = 4, iter = 300)
d=extract(test_mv)

plot(apply(d$log_a_g,2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2,xlab='year')
for(i in 1:13){
  lines(apply(d$log_a_t[,,i],2,median),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}

sum(apply(d$log_lik,2,log_mean_exp))

loo1=loo(test_mv)

#Multivar with evolving global trend
mv_ag_m=rstan::stan_model(model_code=mc)

test_mv2 = rstan::sampling(mv_ag_m, 
                           data = list(N=nrow(skeena_dat),
                                       L=max(skeena_dat$broodyear)-min(skeena_dat$broodyear)+1,
                                       J_i=as.numeric(factor(skeena_dat$stock)),
                                       J_ii=skeena_dat$stock_yr,
                                       J=length(unique(skeena_dat$stock)),
                                       R_S=skeena_dat$logR_S,
                                       S=X_s),
                           control = list(adapt_delta = 0.99,max_treedepth=20), warmup = 100, chains = 4, iter = 300)
d=extract(test_mv2)

plot(apply(d$log_a,2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2,xlab='year')
for(i in 1:13){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}

loo2=loo(test_mv2)
sum(apply(d$log_lik,2,log_mean_exp))

#Independent fit w/ average
test3=rstan::stan_model(model_code=mc_avg)

test_run3 = rstan::sampling(test3, 
                            data = list(N=nrow(skeena_dat),
                                        L=max(skeena_dat$broodyear)-min(skeena_dat$broodyear)+1,
                                        J_i=as.numeric(factor(skeena_dat$stock)),
                                        J_ii=skeena_dat$stock_yr,
                                        J=length(unique(skeena_dat$stock)),
                                        R_S=skeena_dat$logR_S,
                                        S=X_s),
                            control = list(adapt_delta = 0.99,max_treedepth=20), warmup = 100, chains = 4, iter = 300)
d=extract(test_run3)

plot(apply(d$log_a_g,2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:13){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}


loo3=loo(test_run3)

sum(apply(d$log_lik,2,log_mean_exp))

test4=rstan::stan_model(model_code=mc4)

test_run4 = rstan::sampling(test4, 
                            data = list(N=nrow(sk_nass_dat),
                                        L=max(sk_nass_dat$broodyear)-min(sk_nass_dat$broodyear)+1,
                                        J_i=as.numeric(factor(sk_nass_dat$stock)),
                                        J_ii=sk_nass_dat$stock_yr,
                                        J=length(unique(sk_nass_dat$stock)),
                                        R_S=sk_nass_dat$logR_S,
                                        S=X_s),
                            control = list(adapt_delta = 0.99,max_treedepth=20), warmup = 100, chains = 4, iter = 300)
d=extract(test_run4)

plot(apply(d$log_a_t,2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:16){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}

##Compare
loo_mlt=loo(test)


#Fit individually

m3f=samEst::sr_mod(type='rw',par='a',lfo=F)

prod_mat=matrix(ncol=max(sk_nass_stocks$end)-min(sk_nass_stocks$begin)+1,nrow=nrow(sk_nass_stocks))
colnames(prod_mat)=seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end))
rownames(prod_mat)=sk_nass_stocks$stock.name
for(i in 1:16){
  s=subset(sk_nass_dat,stock.id==sk_nass_stocks$stock.id[i])

  f3 = rstan::sampling(m3f, 
                       data = list(N=nrow(s),
                                   L=max(s$broodyear)-min(s$broodyear)+1,
                                   ii=s$broodyear-min(s$broodyear)+1,
                                   R_S=s$logR_S,
                                   S=s$spawners),
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 600)
  
  d_3=rstan::extract(f3)
  med_log_a=apply(d_3$log_a,2,median)
  names(med_log_a)=seq(sk_nass_stocks$begin[i],sk_nass_stocks$end[i])
  
  prod_mat[i,c(match(names(med_log_a),colnames(prod_mat)))]=med_log_a
  }

prod_avg=apply(prod_mat,2,mean,na.rm=TRUE)

plot(prod_avg~colnames(prod_mat),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:16){
  st_pr=prod_mat[i,];st_pr=st_pr[complete.cases(st_pr)]
  lines(st_pr~names(st_pr),col=adjustcolor('black',alpha.f=0.2))
  points(st_pr~names(st_pr),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}



test2=rstan::stan_model(model_code=mc2)

test_run2 = rstan::sampling(test2, 
                           data = list(N=nrow(sk_nass_dat),
                                       L=max(sk_nass_dat$broodyear)-min(sk_nass_dat$broodyear)+1,
                                       J_i=as.numeric(factor(sk_nass_dat$stock)),
                                       J_ii=sk_nass_dat$stock_yr,
                                       J=length(unique(sk_nass_dat$stock)),
                                       R_S=sk_nass_dat$logR_S,
                                       S=X_s),
                           control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 100, chains = 4, iter = 300)
d=extract(test_run2)

plot(apply(d$log_a,2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:16){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
}

test3=rstan::stan_model(model_code=mc_avg)

test_run3 = rstan::sampling(test3, 
                           data = list(N=nrow(skeena_dat),
                                       L=max(skeena_dat$broodyear)-min(skeena_dat$broodyear)+1,
                                       J_i=as.numeric(factor(skeena_dat$stock)),
                                       J_ii=skeena_dat$stock_yr,
                                       J=length(unique(skeena_dat$stock)),
                                       R_S=skeena_dat$logR_S,
                                       S=X_s),
                           control = list(adapt_delta = 0.99,max_treedepth=20), warmup = 100, chains = 4, iter = 300)
d=extract(test_run3)

plot(apply(d$log_a_g,2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:13){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(skeena_stocks$begin),max(skeena_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}

sum(apply(d$log_lik,2,log_mean_exp))

test4=rstan::stan_model(model_code=mc4)

test_run4 = rstan::sampling(test4, 
                            data = list(N=nrow(sk_nass_dat),
                                        L=max(sk_nass_dat$broodyear)-min(sk_nass_dat$broodyear)+1,
                                        J_i=as.numeric(factor(sk_nass_dat$stock)),
                                        J_ii=sk_nass_dat$stock_yr,
                                        J=length(unique(sk_nass_dat$stock)),
                                        R_S=sk_nass_dat$logR_S,
                                        S=X_s),
                            control = list(adapt_delta = 0.99,max_treedepth=20), warmup = 100, chains = 4, iter = 300)
d=extract(test_run4)

plot(apply(d$log_a_t,2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:16){
  lines(apply(d$log_a_t[,,i],2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),col=adjustcolor('black',alpha.f=0.2))
  points(apply(d$log_a_t[,,i],2,median)~seq(min(sk_nass_stocks$begin),max(sk_nass_stocks$end)),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}


loo_ind=loo(test_run2)
loo_multi=loo(test_run)
loo_ind_2=loo(test_run3)



#Chinook
chi_stocks=subset(stock_info_filtered,species=='Chinook')
chi_dat=subset(stock_dat2,stock.id %in% chi_stocks$stock.id)
length(unique(chi_dat$stock.id))

stock_year=expand.grid(unique(chi_dat$stock),unique(chi_dat$broodyear))
stock_year=stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")

chi_dat$stock_yr=stock_year[match(paste(chi_dat$stock,chi_dat$broodyear,sep='_'),stock_year[,3]),3]

X_s=make_design_matrix(chi_dat$spawners,grp=chi_dat$stock)

test_run3 = rstan::sampling(test2, 
                            data = list(N=nrow(chi_dat),
                                        L=max(chi_dat$broodyear)-min(chi_dat$broodyear)+1,
                                        J_i=as.numeric(factor(chi_dat$stock)),
                                        J_ii=as.numeric(factor(chi_dat$stock_yr)),
                                        J=length(unique(chi_dat$stock)),
                                        R_S=chi_dat$logR_S,
                                        S=X_s),
                            control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 100, chains = 4, iter = 300,init=1)
d=extract(test_run3)

summary(test_run3,pars=c('sigma'))

plot(apply(d$log_a,2,median),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:16){
  lines(apply(d$log_a_t[,,i],2,median),col=adjustcolor('black',alpha.f=0.2))
}



##Chinook
#Fit individually
chi_stocks=subset(stock_info_filtered,species=='Chinook')
chi_dat=subset(stock_dat2,stock.id %in% chi_stocks$stock.id)
length(unique(chi_dat$stock.id))

m3f=samEst::sr_mod(type='rw',par='a',lfo=F)

prod_mat=matrix(ncol=max(chi_stocks$end)-min(chi_stocks$begin)+1,nrow=nrow(chi_stocks))
colnames(prod_mat)=seq(min(chi_stocks$begin),max(chi_stocks$end))
rownames(prod_mat)=chi_stocks$stock.name
for(i in 1:nrow(chi_stocks)){
  s=subset(chi_dat,stock.id==chi_stocks$stock.id[i])
  
  f3 = rstan::sampling(m3f, 
                       data = list(N=nrow(s),
                                   L=max(s$broodyear)-min(s$broodyear)+1,
                                   ii=s$broodyear-min(s$broodyear)+1,
                                   R_S=s$logR_S,
                                   S=s$spawners),
                       control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 600)
  
  d_3=rstan::extract(f3)
  med_log_a=apply(d_3$log_a,2,median)
  names(med_log_a)=seq(chi_stocks$begin[i],chi_stocks$end[i])
  
  prod_mat[i,c(match(names(med_log_a),colnames(prod_mat)))]=med_log_a
}

prod_avg=apply(prod_mat,2,mean,na.rm=TRUE)

plot(prod_avg~colnames(prod_mat),ylim=c(-1,4),type='l',bty='l',ylab='Productivity - log alpha',lwd=2)
for(i in 1:16){
  st_pr=prod_mat[i,];st_pr=st_pr[complete.cases(st_pr)]
  lines(st_pr~names(st_pr),col=adjustcolor('black',alpha.f=0.2))
  points(st_pr~names(st_pr),col=adjustcolor('black',alpha.f=0.15),cex=0.9)
}



#By species
mc="data{
  int<lower=1> N;//number of annual samples 
  int L; //years covered by time-series
  int J;//number of stocks
  int J_i[N];//index of stocks
  int J_ii[N];//index of stock-year combination
  vector[N] R_S; //log(recruits per spawner)
  matrix[N,J] S; //design matrix of spawners in time T
  }
parameters{
  real log_a0; //inital productivity (global)
  vector[J] log_a0_s;// initial productivity (stock-specific)
  vector<upper = 0>[J] log_b; // rate capacity - fixed in this

 //variance components  
 vector<lower = 0>[J] sigma;

//MVN parameters  
  cholesky_factor_corr[J] Lcorr;
  vector[L-1] z_dev; //average deviation in stock productivity among stocks
  matrix[L-1,J] z_dev_s;  // deviations in stock productivity from year-to-year
  real<lower = 0> sigma_a; //variance in average productivity among years
  vector<lower = 0>[J] sigma_a_s; //variance in stock productivity among years
  
}
transformed parameters{
  vector[J] b; //capacity rate
  //global state process
  vector[L] log_a; //global (ie. average) productivity among stocks over time
  //stock state process
  matrix[L,J] log_a_s; //stock deviation from global log_a over time
  matrix[L-1,J] a_dev; //stock-level deviations in year-to-year productivity
  matrix[L,J] log_a_t; //stock productivity at time t (local + global)
  
  b=exp(log_b);
  
  //initial productivities
  log_a[1] = log_a0; //average productivity among stocks at t = 1
  log_a_s[1,] = to_row_vector(log_a0_s); //stock deviation in productivity at t =1
  log_a_t[1,] = log_a[1]+to_row_vector(log_a0_s); //global prod + stock deviation from global prod.
  
  for(t in 1:L-1){
   a_dev[t,] = (diag_pre_multiply(sigma_a_s, Lcorr) * to_vector(z_dev_s[t,]))';
  }
  
  for(t in 2:L){
    log_a[t] = log_a[t-1] + z_dev[t-1]*sigma_a; //global prod. random walk
    log_a_s[t,] = log_a_s[t-1,] + a_dev[t-1,]; //stock-specific deviation in random walk
    log_a_t[t,] = log_a[t] + log_a_s[t,]; //final estimate of stock prod. (global + stock)
  }
  
}  
model{
  //priors
  log_a0 ~ gamma(3,1.5); //initial globalproductivity - wide prior
  log_a0_s ~ normal(0,1); //stock-level deviations around the global productivity
  log_b ~ normal(-12,3); //per capita capacity parameter - wide prior
  
  z_dev ~ std_normal(); //standardized stock-level deviances in prod
  to_vector(z_dev_s) ~ std_normal(); //global deviation in prod.
  
  Lcorr ~ lkj_corr_cholesky(2.0); // prior for correlation of process deviances
 
  //variance terms
  target += normal_lpdf(sigma | 0, 1) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior  
  target += normal_lpdf(sigma_a_s | 0, 0.5) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior
  target += normal_lpdf(sigma_a| 0, 0.5) - normal_lcdf(0 | 0, 1); //remove density below zero to accomodate half-normal prior
 

  R_S ~ normal(to_vector(log_a_t)[J_ii] - S*b, sigma[J_i]); 
  
}
generated quantities {
  vector[N] log_lik;

  for(i in 1:N) log_lik[i] = normal_lpdf(R_S[i]|to_vector(log_a_t)[J_ii[i]] - S[i]*b, sigma[J_i[i]]);
  corr_matrix[J] Cor_t = multiply_lower_tri_self_transpose(Lcorr);
}


"

#Chinook####

chi_stocks=subset(stock_info_filtered,species=='Chinook')
chi_dat=subset(stock_dat2,stock.id %in% chi_stocks$stock.id)
length(unique(chi_dat$stock.id))

test=rstan::stan_model(model_code=multi_avg)


stock_year=expand.grid(unique(chi_dat$stock),unique(chi_dat$broodyear))
stock_year=stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")

chi_dat$stock_yr=match(paste(chi_dat$stock,chi_dat$broodyear,sep='_'),stock_year[,3])

X_s=make_design_matrix(chi_dat$spawners,grp=chi_dat$stock)

chi_indfit_avg = rstan::sampling(test, 
                            data = list(N=nrow(chi_dat),
                                        init=0,
                                        L=max(chi_dat$broodyear)-min(chi_dat$broodyear)+1,
                                        J_i=as.numeric(factor(chi_dat$stock)),
                                        J_ii=chi_dat$stock_yr,
                                        J=length(unique(chi_dat$stock)),
                                        R_S=chi_dat$logR_S,
                                        S=X_s),
                            control = list(adapt_delta = 0.99,max_treedepth=20),init=2, warmup = 100, chains = 1, iter = 300)

coh_stocks=subset(stock_info_filtered,species=='Coho')
coh_dat=subset(stock_dat2,stock.id %in% coh_stocks$stock.id)
length(unique(coh_dat$stock.id))

test=rstan::stan_model(model_code=multi_avg)


stock_year=expand.grid(unique(coh_dat$stock),unique(coh_dat$broodyear))
stock_year=stock_year[order(stock_year[,1]),]
stock_year[,3]=paste(stock_year[,1],stock_year[,2],sep="_")

coh_dat$stock_yr=match(paste(coh_dat$stock,coh_dat$broodyear,sep='_'),stock_year[,3])

X_s=make_design_matrix(coh_dat$spawners,grp=coh_dat$stock)

coh_indfit_avg = rstan::sampling(test, 
                                 data = list(N=nrow(coh_dat),
                                             init=0,
                                             L=max(coh_dat$broodyear)-min(coh_dat$broodyear)+1,
                                             J_i=as.numeric(factor(coh_dat$stock)),
                                             J_ii=coh_dat$stock_yr,
                                             J=length(unique(coh_dat$stock)),
                                             R_S=coh_dat$logR_S,
                                             S=X_s),
                                 control = list(adapt_delta = 0.99,max_treedepth=20),init=2, warmup = 100, chains = 4, iter = 300)

