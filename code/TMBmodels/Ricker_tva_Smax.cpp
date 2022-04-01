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
  
  //logbeta     -> log of beta from ricker curve
  //alphao      -> initial alpha value

  PARAMETER(alphao);
  PARAMETER(logbeta);
  PARAMETER(logsigobs);
  PARAMETER(logsiga);
 
  PARAMETER_VECTOR(alpha);
    
  int timeSteps=obs_logRS.size();

  Type beta = exp(logbeta);
  Type Smax = Type(1.0)/beta;


  //Type varphi     = exp(logvarphi);
  //Type theta     = sqrt(Type(1.0)/varphi);
  //Type sigobs       = sqrt(rho) * theta;
  //Type siga        = sqrt(Type(1.0)-rho) * theta ;

  Type sigobs = exp(logsigobs);
  Type siga = exp(logsiga);


  vector<Type> pred_logR(timeSteps), pred_logRS(timeSteps),umsy(timeSteps), Smsy(timeSteps), residuals(timeSteps);
  vector<Type> Srep(timeSteps); //, alpha(timeSteps)

 
  //priors on parameters
 
  Type ans= Type(0);
  ans -=dnorm(alphao,Type(0.0),Type(5.0),true);
  ans -=dstudent(logbeta,Type(-8.0),Type(10.0),Type(4.0),true);
  
  ans -= dnorm(sigobs,Type(0.0),Type(2.0),true);
  ans -= dnorm(siga,Type(0.0),Type(2.0),true);
  //ans -= dexp(sigobs,Type(2.0),true);
  //ans -= dexp(siga,Type(2.0),true);
  //ans -= dt(sigobs,Type(3.0),true);
  //ans -= dt(siga,Type(3.0),true);
  //ans -= dbeta(rho,Type(3.0),Type(3.0),true);  


  
  // Use the Hilborn approximations for Smsy and umsy
  //alpha(0) = alphao;
  ans+= -dnorm(alpha(0),alphao,siga,true);
  //umsy(0)  = Type(.5) * alpha(0) - Type(0.07) * (alpha(0) * alpha(0));
  //Smsy(0)  =  alpha(0)/beta * (Type(0.5) -Type(0.07) * alpha(0));  
  //Srep(0)  = alpha(0)/beta;
  
  for(int i=1;i<timeSteps;i++){
  
    ans+= -dnorm(alpha(i),alpha(i-1),siga,true);
    //alpha(i) = alpha(i-1) + alphadev(i-1)*siga;
    //ans+= -dnorm(alphadev(i-1),Type(0.0),Type(1.0),true);
  
  }

  for(int i=0;i<timeSteps;i++){
    if(!isNA(obs_logRS(i))){
      
      pred_logRS(i) = alpha(i) - beta * obs_S(i) ;
      pred_logR(i) = pred_logRS(i) + log(obs_S(i));

      umsy(i) = Type(.5) * alpha(i) - Type(0.07) * (alpha(i) * alpha(i)); 
      Smsy(i) =  alpha(i)/beta * (Type(0.5) -Type(0.07) * alpha(i));
      Srep(i) = alpha(i)/beta;

      residuals(i) = obs_logRS(i) - pred_logRS(i);
      ans+=-dnorm(obs_logRS(i),pred_logRS(i),sigobs,true);
    }
  
  }

  REPORT(pred_logRS)
  REPORT(alpha)
  REPORT(alpha)
  REPORT(sigobs)
  REPORT(siga)
  REPORT(residuals)
  REPORT(beta)
  REPORT(alphao)
  REPORT(Smax)
  REPORT(umsy)
  REPORT(Smsy)
  REPORT(Srep)

  return ans;
}

