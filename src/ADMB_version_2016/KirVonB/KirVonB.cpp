#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  
    DATA_VECTOR(age);
    DATA_VECTOR(weight);  
    DATA_VECTOR(Temp);
    // DATA_VECTOR(R);
    // DATA_VECTOR(Sp);
    // DATA_VECTOR(Sex);
    // DATA_VECTOR(L);

    int nobs = age.size();  // integer of y

    PARAMETER(H);
    PARAMETER(K);
    PARAMETER(t0);
    PARAMETER(log_mean_d);
    PARAMETER(logSigma);
    PARAMETER(Tcoef);

    vector<Type> d(nobs);  // declare local integer
    vector<Type> Winf(nobs);  // declare local integer
    vector<Type> logWhat(nobs);  // declare local integer

    Type neglogL = 0.0;

    for (int i=0;i<nobs;i++)
    {
      d(i)          = exp(log_mean_d+Tcoef*Temp(i)); //mfexp(log_mean_d+Pcoef*Pval_lag1(i)+Tcoef*Temp(i));
      Winf(i)       = pow((H/K),1.0/(1.0 - d(i)) );
      logWhat(i)    =  log(Winf(i)) + (1.0/(1.0 - d(i)))*log(1.0 - exp(-K * (1.0 - d(i)) * (age(i) - t0))) ;
    }

    neglogL = -sum(dnorm(vector<Type>(log(weight)), logWhat, exp(logSigma), true));


    REPORT(log_mean_d);
    REPORT(H);
    REPORT(K);
    REPORT(Tcoef);
    REPORT(t0);
    REPORT(logSigma);
    // REPORT(Winf);
    // REPORT(d);
    // REPORT(logWhat);

    ADREPORT(log_mean_d);
    ADREPORT(H);
    ADREPORT(K);
    ADREPORT(Tcoef);
    ADREPORT(t0);
    ADREPORT(logSigma);
    // ADREPORT(Winf);
    // ADREPORT(d);
    // ADREPORT(logWhat);
    return neglogL;

}




