#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

 // dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres;
  logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  if(give_log)return logres; else return exp(logres);
}

template <class Type>
Type dstudent(Type x, Type mean, Type sigma, Type df, int give_log = 0) {
  // from metRology::dt.scaled()
  // dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - log(sd)
  Type logres = dt((x - mean) / sigma, df, true) - log(sigma);
  if (give_log)
    return logres;
  else
    return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs_logRS);   // observed log recruitment
  DATA_VECTOR(obs_S);    // observed  Spawner
  
  PARAMETER(alpha);
  PARAMETER(logbeta);
  PARAMETER(logsigobs);
  
  int timeSteps=obs_logRS.size();

  //priors - based on evaluation done with the prior predictive check
  Type ans= Type(0);
  ans -=dnorm(alpha,Type(0.0),Type(2.5),true);
  ans -=dnorm(logbeta,Type(-12.0),Type(3.0),true);
  
  ans -= dnorm(logsigobs,Type(0.0),Type(2.0),true);
  //ans -= dexp(sigobs,Type(2.0),true);
  //ans -= dt(sigobs,Type(3.0),true);
  
  //model
  Type beta = exp(logbeta);
  Type sigobs = exp(logsigobs);
  Type Smax  = Type(1.0)/beta;
  
  vector<Type> pred_logRS(timeSteps), pred_logR(timeSteps), residuals(timeSteps); 
   
  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logRS(i))){
      pred_logRS(i) = alpha - beta * obs_S(i) ; 
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));
      residuals(i) = obs_logRS(i) - pred_logRS(i);
      ans+=-dnorm(obs_logRS(i),pred_logRS(i),sigobs,true);
    }
  
  }
  Type umsy = Type(.5) * alpha - Type(0.07) * (alpha * alpha);
  Type Smsy = alpha/beta * (Type(0.5) -Type(0.07) * alpha);

  SIMULATE {
   vector<Type> R_Proj(timeSteps);
   for(int i=0; i<timeSteps; ++i){
     R_Proj(i) = exp(rnorm(pred_logR(i), sigobs));
   }
   REPORT(R_Proj);
  }

  REPORT(pred_logR)
  REPORT(pred_logRS)
  REPORT(residuals)
  REPORT(alpha)  
  REPORT(beta)
  REPORT(sigobs)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  
  return ans;
}

