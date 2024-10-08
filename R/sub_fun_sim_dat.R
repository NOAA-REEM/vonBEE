#'sim_dat()
#'
#'simulate the TMB model in for VonBT in R
#'
#' @param par =dlistIN_0$par
#' @param data  = datIN
#' @param simulateIT =FALSE
#'
#' @export
#'
#' @returns list(out = cbind(data.frame(W_obs = exp(logWobs),W_hat = What, age=age),t(vonb_cov)),What=What,Wobs=exp(logWobs),age=age, vonb_cov=vonb_cov,d=d, Winf=Winf,nll=nll)
#'
sim_dat <- function(par=dlistIN_0$par,data = datIN,simulateIT=FALSE){

  for(nm in names(par))
    eval(parse(text=paste0(nm,"<-par$",nm)))
  for(nm in names(data))
    eval(parse(text=paste0(nm,"<-data$",nm)))
  sigma_proc <- exp(log_sigma_proc)
  sigma_obs  <- exp(log_sigma_obs)
  cov     <- vonb_cov*0
  d       <- rep(0,nobs)
  Winf    <- rep(0,nobs)
  logWhat <- rep(0,nobs)
  What    <- rep(0,nobs)
  logWobs <- rep(0,nobs)
  logWobs = log(weight)
  H   <- exp(logH)
  K   <- exp(logK)


  nll <- 0

  for(y in 2:nyrs){

    m = 0 + beta_1*u_y[y-1] # Gompertz

    # process model:
    #nll<- nll - dnorm(u_y[y], m, sigma_proc, TRUE)
  }

  for( i in 1:nobs){

    covars = 0 # Initialize

    for ( c in 1:ncov ){
      cov[c,i] <- vonb_cov[ c,i ]
      if(simulateIT) {
        cov[c,i] <- rnorm(1,  vonb_cov[ c,i ] , sdvonb_cov[ c,i ]) # Simulate env
      }
      covars<- covars+ beta_c[ c ]* cov[ c,i ]
    }
    x             = mu_d + u_y[YrIndex[i]+1] + beta_0[cohortYr[i]+1] + covars
    d[i]          =  ( 1 - d_offset )/( 1 + exp(-x) )
    Winf[i]       = (H/K)^(1.0/(1.0 - d[i]) )
    logWhat[i]    = log(Winf[i]) + ( 1.0/(1.0 - d[i]) )*log(1.0 - exp(-K * (1.0 - d[i]) * (age[i] - t0)))
    What[i]       = exp(logWhat[i])

    # observation model:
    # nll = nll- dnorm(logWobs[i], logWhat[i], sigma_obs, TRUE);

  }
  return(list(out = cbind(data.frame(W_obs = exp(logWobs),W_hat = What, age=age),t(vonb_cov)),
              What=What,Wobs=exp(logWobs),age=age, vonb_cov=vonb_cov,d=d, Winf=Winf,nll=nll))
  detach()
}
