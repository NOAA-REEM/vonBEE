#include <TMB.hpp>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Eigen>
  // ------------------------------------------------------------------------- //
  //                    vonBT version 1.0.1                                    //
  //       vonB growth model with environmental covariates                     //
  //                                                                           //
  // AUTHORS:   K Holsman                                                      //
  // CITATIONS:                                                                //
  // 0. TBA et al. 2019                                                        //
  // 1. Holsman, et al. 2016                                                   //
  // ------------------------------------------------------------------------- //
  // 
  //  INDEX 
  //  0. Load data
  //  1. Estimated Parameters
  //  2. Parameters
  //  3. Init. calcs
  //  4. vonB model
  //  5. Calc obj fun
  //  6. Report section
  //  7. End
  
template<class Type>
Type objective_function<Type>::operator() ()
{

  // 0. LOAD DATA
  // ------------------------------------------------------------------------- //

    DATA_IVECTOR(YrIndex);        // sample year index
    DATA_VECTOR(age);             // fish age
    DATA_VECTOR(weight);          // weight
    DATA_IVECTOR(cohortYr);       // cohort yr
    DATA_INTEGER( ncov );         // number of environmental covariates
    DATA_INTEGER( nobs );         // number of observations
    DATA_INTEGER( nyrs );         // number of years
    DATA_INTEGER(d_offset);       // default to .1
    DATA_MATRIX(vonb_cov);        // environmental covariates
    DATA_MATRIX(sdvonb_cov);      // stdev environmental covariates
    // int ncov = vonb_cov.size();   // integer of cov
    // int nobs = age.size();        // number of observations
    
  // 1. DEFINE ESTIMATED PARAMETERS
  // ------------------------------------------------------------------------- //
    PARAMETER(logH);          
    PARAMETER(logK);
    PARAMETER(t0);
    PARAMETER(mu_d);
    PARAMETER(log_sigma_proc);     // log(process SD)
    PARAMETER(log_sigma_obs);      // log(observation SD)
    PARAMETER_VECTOR(beta_c);      // vector of ncovs length
    PARAMETER(beta_1);             // AR1 beta
    PARAMETER_VECTOR(beta_0);      // cohort beta
    PARAMETER_VECTOR(u_y);         // unobserved state vector length years
    
  // 2. DEFINE PARAMETERS
  // ------------------------------------------------------------------------- //
    vector<Type> d(nobs);          // declare local vector
    vector<Type> Winf(nobs);       // declare local vector
    vector<Type> logWhat(nobs);    // declare local vector
    vector<Type> What(nobs);       // declare local vector
    vector<Type> logWobs(nobs);    // declare local vector
    matrix<Type> cov( ncov,nobs ); // declare local matrix
    
  // 3. Init. Calcs
  // ------------------------------------------------------------------------- //
  // procedures: (transformed parameters)
     Type sigma_proc = exp(log_sigma_proc);
     Type sigma_obs  = exp(log_sigma_obs);
     Type H  = exp(logH);
     Type K  = exp(logK);
    
     for (int i = 0; i < nobs; i++)
      logWobs(i) = log(weight(i));
  
  // reports on transformed parameters:
    ADREPORT(sigma_proc);
    ADREPORT(sigma_obs);
    
  // 4. Fit model
  // ------------------------------------------------------------------------- //
    Type nll_proc = 0.0; // initialize negative log likelihood
    Type nll_obs  = 0.0; // initialize negative log likelihood
    Type nll      = 0.0; // initialize negative log likelihood
    
    for(int y = 1; y < nyrs; y++){
      Type m = beta_1*u_y(y-1); // Gompertz
      // process model:
      nll_proc -= dnorm(u_y(y), m, sigma_proc, true);
    }
    
    for(int i = 0; i < nobs; i++)
    {
        Type covars  = Type(0.0);
        for ( int c = 0; c < ncov; c++ ){
          cov(c,i) = vonb_cov( c,i );
          SIMULATE {
            cov( c,i ) = rnorm(vonb_cov( c,i ) , sdvonb_cov( c,i)) ; // Simulate env
          }
          covars+= beta_c( c )* cov( c,i );
        }
      Type x        = mu_d + u_y(YrIndex(i)) + beta_0(cohortYr(i)) + covars;
      d(i)          = ( 1 - d_offset )/( 1 + exp(-x) );
      Winf(i)       = pow((H/K),1.0/(1.0 - d(i)) );
      logWhat(i)    = log(Winf(i)) + (1.0/(1.0 - d(i)))*log(1.0 - exp(-K * (1.0 - d(i)) * (age(i) - t0))) ;
      What(i)       = exp(logWhat(i));
     
      
      // observation model:
      nll_obs -= dnorm(logWobs(i), logWhat(i), sigma_obs, true);
    }
    
    // 5. calc obj fun
    // ------------------------------------------------------------------------- //
    nll = nll_obs + nll_proc;

  // 6. report section
  // ------------------------------------------------------------------------- //
    
    REPORT(mu_d);
    REPORT(H);
    REPORT(K);
    REPORT(t0);
    REPORT(beta_c);
    REPORT(beta_1);
    REPORT(beta_0);
    REPORT(u_y);
    REPORT(log_sigma_obs);
    REPORT(log_sigma_proc);
    REPORT(nll_obs);
    REPORT(nll_proc);
    REPORT(nll);
    REPORT(Winf);
    REPORT(d);
    REPORT(cov);
    REPORT(logWhat);
    REPORT(logWobs);
    REPORT(What);
    REPORT(weight);
    REPORT(age);
    
    //this for any pars to estimate sd 
    ADREPORT(mu_d);
    ADREPORT(H);
    ADREPORT(K);
    ADREPORT(beta_c);
    ADREPORT(t0);
    ADREPORT(log_sigma_obs);
    ADREPORT(log_sigma_proc);
    
    return nll;
    
    // Quoting from ?sdreport: If ignore.parm.uncertainty=TRUE
    // then the Hessian calculation is omitted and a zero-matrix is used in
    //   place of V(theta).
    // 
    // If you want parameters without sds you can get them directly from the
    // optimizer output. Or you can grab them from the object environment
    // after running the optimizer (see ?MakeADFun): obj$env$last.par.best
    // 
    // Get the parameter estimates in the original list structure using: obj$env$parList(par=obj$env$last.par.best)
    //   
    //   The REPORT / ADREPORT macros are for derived quantities. REPORT is
    //   fastest and can output general objects (matrices, arrays,etc). 
    // ADREPORT is different. It joins all output to a single vector
    //   and performs AD on the result.
    // 
}




