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

  DATA_VECTOR(obs_logRS);   // observed recruitment
  DATA_VECTOR(obs_S);    // observed  Spawner
  
  //DATA_SCALAR(prbeta1); //beta prior parameter
  //DATA_SCALAR(prbeta2); //beta prior parameter
  //DATA_SCALAR(prbeta3); //beta prior parameter
  //DATA_SCALAR(prbeta4); //beta prior parameter
  
  //logbeta     -> log of beta from ricker curve
  //alphao      -> initial alpha value
  //rho         -> Proportion of total variance associated with obs error.
  //varphi      -> Total precision
  //alpha       -> Time-varying alpha

  PARAMETER(logbetao);
  PARAMETER(alphao);
  
  PARAMETER(logsigobs);
  PARAMETER(logsiga);
  PARAMETER(logsigb);

  //PARAMETER(rho);
  //PARAMETER(logvarphi);
  //PARAMETER(kappa);


  PARAMETER_VECTOR(logbeta);
  PARAMETER_VECTOR(alpha);
  
  
  int timeSteps=obs_logRS.size();

  Type sigobs   = exp(logsigobs);
  Type siga   = exp(logsiga);  
  Type sigb   = exp(logsigb); 
  
  
  //theta       -> total standard deviation
  //sig         -> obs error std
  //tau         -> proc error (alpha and beta) sd
  
  //Type varphi = exp(logvarphi);
  //Type theta  = sqrt(Type(1.0)/varphi);
  //Type sig    = sqrt(rho) * theta;
  //Type tau    = sqrt(Type(1.0)-rho) * theta ;
  //Type siga   = sqrt(kappa) * tau;  
  //Type sigb   = sqrt(Type(1.0)-kappa) * tau; 



  vector<Type> pred_logR(timeSteps), pred_logRS(timeSteps), Smsy(timeSteps), residuals(timeSteps);
  vector<Type> Srep(timeSteps), beta(timeSteps), Smax(timeSteps), umsy(timeSteps);

  
  Type ans= Type(0); 

  //prior on parameters
  ans -=dnorm(alphao,Type(0.0),Type(5.0),true);
  ans -=dstudent(logbetao,Type(-8.0),Type(10.0),Type(4.0),true);
  //prior on observation and process variance ratio
  //Type ans= -dbeta(rho,prbeta1,prbeta2,true); 
  //ans= -dbeta(kappa,prbeta3,prbeta4,true);  
  ans -= dexp(sigobs,Type(2.0),true);
  ans -= dexp(siga,Type(2.0),true);
  ans -= dexp(sigb,Type(2.0),true);
  
  //ans -= dnorm(sigobs,Type(0.0),Type(2.0),true);
  //ans -= dnorm(siga,Type(0.0),Type(2.0),true);
  //ans -= dnorm(sigb,Type(0.0),Type(2.0),true);

  

  
  ans+= -dnorm(alpha(0),alphao,siga,true);
  ans+= -dnorm(logbeta(0),logbetao,sigb,true);

  // Use the Hilborn approximations for Smsy and umsy  

  
  for(int i=1;i<timeSteps;i++){
  
    ans+= -dnorm(alpha(i),alpha(i-1),siga,true);
    ans+= -dnorm(logbeta(i),logbeta(i-1),sigb,true);
  
  }

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logRS(i))){
      beta(i) = exp(logbeta(i));
      pred_logRS(i) = alpha(i) - beta(i) * obs_S(i) ;
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));

      Smsy(i) =  alpha(i)/beta(i) * (Type(0.5) -Type(0.07) * alpha(i));
      Srep(i) = alpha(i)/beta(i);
      umsy(i)  = Type(.5) * alpha(i) - Type(0.07) * (alpha(i) * alpha(i));

      residuals(i) = obs_logRS(i) - pred_logRS(i);
      ans+=-dnorm(obs_logRS(i),pred_logRS(i),sigobs,true);
    }
  
  }

  REPORT(pred_logRS)
  REPORT(alpha)
  REPORT(sigobs)
  REPORT(siga)
  REPORT(sigb)
  //REPORT(kappa)
  //REPORT(theta)
  REPORT(residuals)
  REPORT(beta)
  //REPORT(varphi)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  REPORT(Srep)
  return ans;
}

