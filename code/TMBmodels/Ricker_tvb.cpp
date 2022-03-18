#include <TMB.hpp>

// Code taken from Brooke
// Set up Lambert's W function to use to calculate SMSY
// Code taken from https://kaskr.github.io/adcomp/lambert_8cpp_source.html
// Step 1: Code up a plain C version
// Double version of Lambert W function
double LambertW(double x) {
  double logx = log(x);
  double y = (logx > 0 ? logx : 0);
  int niter = 100, i=0;
  for (; i < niter; i++) {
    if ( fabs( logx - log(y) - y) < 1e-9) break;
    y -= (y - exp(logx - y)) / (1 + y);
  }
  if (i == niter) Rf_warning("W: failed convergence");
  return y;
}

TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  LambertW
  ,
  // OUTPUT_DIM
  1,
  // ATOMIC_DOUBLE
  ty[0] = LambertW(tx[0]); // Call the 'double' version
,
// ATOMIC_REVERSE
Type W  = ty[0];                    // Function value from forward pass
Type DW = 1. / (exp(W) * (1. + W)); // Derivative
px[0] = DW * py[0];                 // Reverse mode chain rule
)
  
  // Scalar version
  template<class Type>
  Type LambertW(Type x){
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return LambertW(tx)[0];
  }
  

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

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(obs_logR);   // observed recruitment
  DATA_VECTOR(obs_S);    // observed  Spawner
  
  DATA_SCALAR(prbeta1); //beta prior parameter
  DATA_SCALAR(prbeta2); //beta prior parameter
  
  //logbeta     -> log of beta from ricker curve
  //alphao      -> initial alpha value
  //rho         -> Proportion of total variance associated with obs error.
  //varphi      -> Total precision
  //alpha       -> Time-varying alpha

  PARAMETER(logbetao);
  PARAMETER(alpha);
  PARAMETER(rho);
  PARAMETER(logvarphi);


  PARAMETER_VECTOR(logbeta);
  
  
  int timeSteps=obs_logR.size();

  
  
  
  //theta       -> total standard deviation
  //sig         -> obs error std
  //tau         -> proc error (beta) sd
  
  Type varphi     = exp(logvarphi);
  Type theta     = sqrt(Type(1.0)/varphi);
  Type sig       = sqrt(rho) * theta;
  Type tau        = sqrt(Type(1.0)-rho) * theta ;


  vector<Type> pred_logR(timeSteps), pred_logRS(timeSteps), Smsy(timeSteps), residuals(timeSteps);
  vector<Type> Srep(timeSteps), beta(timeSteps), Smax(timeSteps);

  

  //prior on observation and process variance ratio
  Type ans= -dbeta(rho,prbeta1,prbeta2,true);  
  
  ans+= -dnorm(logbeta(0),logbetao,tau,true);

  // Use the Hilborn approximations for Smsy and umsy
  beta(0) = exp(logbeta(0));
  Type umsy  = Type(.5) * alpha - Type(0.07) * (alpha * alpha);
  Smsy(0)  =  alpha/beta(0) * (Type(0.5) -Type(0.07) * alpha);  
  Srep(0)  = alpha/beta(0);

  
  for(int i=1;i<timeSteps;i++){
  
    ans+= -dnorm(logbeta(i),logbeta(i-1),tau,true);
  
  }

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logR(i))){
      beta(i) = exp(logbeta(i));
      pred_logRS(i) = alpha - beta(i) * obs_S(i) ;
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));

      Smsy(i) =  alpha/beta(i) * (Type(0.5) -Type(0.07) * alpha);
      Srep(i) = alpha/beta(i);

      residuals(i) = obs_logR(i) - pred_logR(i);
      ans+=-dnorm(obs_logR(i),pred_logR(i),sig,true);
    }
  
  }

  REPORT(pred_logR)
  REPORT(alpha)
  REPORT(sig)
  REPORT(tau)
  REPORT(rho)
  REPORT(theta)
  REPORT(residuals)
  REPORT(beta)
  REPORT(varphi)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  REPORT(Srep)
  return ans;
}

