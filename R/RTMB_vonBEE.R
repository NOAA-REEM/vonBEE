vonBEE <- function(pars){
  require(RTMB)
  RTMB::getAll(pars, data_list)

  # Parameter transform ----
  sigmaProc = exp(log_sigma_proc)
  sigmaObs = exp(log_sigma_obs)
  H = exp(logH)
  K = exp(logK)
  rho = 1/(1+exp(-logit_rho))

  # Data transform ----
  logWeightObs = log(weight)

  # Model ----
  x = mu_d + x_mat %*% beta_x + u_y[YrIndex] + beta_c[cohortYr]
  d          = ( 1 - d_offset )/( 1 + exp(-x) )
  Winf       = (H/K) ^ (1.0/(1.0 - d))
  logWhat    = log(Winf) + (1.0/(1.0 - d))*log(1.0 - exp(-K * (1.0 - d) * (age - t0)))
  What       = exp(logWhat)

  # Likelihood ----
  nll_proc = -dnorm(u_y[1], 0, sigmaProc, TRUE)
  for(i in 2:nyrs){
    nll_proc = nll_proc - dnorm(u_y[i], rho*u_y[i-1], sigmaProc, TRUE)
  }

  nll_obs = -sum(dnorm(logWeightObs, logWhat, sigmaObs, TRUE))
  nll = nll_obs + nll_proc;

  # Report ----
  RTMB::REPORT(sigmaProc)
  RTMB::REPORT(sigmaObs)
  RTMB::REPORT(H)
  RTMB::REPORT(K)
  RTMB::REPORT(rho)
  RTMB::REPORT(x)
  RTMB::REPORT(d)
  RTMB::REPORT(Winf)
  RTMB::REPORT(What)

  return(nll)
}



FitVonBEE <- function(dat){

  # - Rearrange data
  dat <- dat %>%
    mutate(
      YrIndex = YEAR - min(YEAR) + 1,
      cohortYr   = YEAR-age,
      cohortIndex = cohortYr - min(cohortYr) + 1,
      cohortIndexFactor = as.numeric(factor(cohortIndex))) %>%
    group_by(YEAR) %>%
    mutate(mnTemp = mean(Temp,na.rm=T), # Probably want to scale these
           seTemp = sd(Temp,na.rm=T),
           mnTemp2 = mean(Temp^2, na.rm =T),
           seTemp2 = sd(Temp^2, na.rm =T)) %>%
    ungroup() %>%
    group_by(cohortYr) %>%
    mutate(chort_mnTemp = mean(Temp,na.rm=T),
           chort_seTemp = sd(Temp,na.rm=T),
           chort_mnTemp2 = mean(Temp^2, na.rm =T),
           chort_seTemp2 = sd(Temp^2, na.rm =T)) %>%
    ungroup()

  data_list <- list(
    YrIndex    = factor(dat$YrIndex),
    age        = dat$age,
    weight     = dat$weight,
    cohortYr   = factor(dat$cohortIndex),
    ncov       = 3,
    nobs       = nrow(dat),
    nyrs       = length(unique(dat$YEAR)),
    d_offset   = 0 ,#0.01,  # prevents WInf --> 0
    x_mat   = dat %>%
      dplyr::select(mnTemp, mnTemp2, chort_mnTemp) %>%
      as.matrix()
  )


  # Parameters ----
  par_list <- list(
    logH   = log(12.0),
    logK   = log(0.15),
    t0  = -0.25,
    mu_d = -.005,
    log_sigma_proc   = -.6,
    log_sigma_obs   = -.3,
    logit_rho     = 0,
    beta_c     = rep(0,length(unique(data_list$cohortYr))),
    beta_x     = rep(0, data_list$ncov) ,
    u_y        = rep(0,data_list$nyrs)
  )

  # Build and fit ----
  obj <- MakeADFun(vonBEE, par_list, map = map, silent = TRUE, random = "u_y")
  fit <- optim(par = obj$par,
               fn = obj$fn,
               gr = obj$gr,
               control = list(maxit = 1e6))
  report <- obj$report(obj$env$last.par.best)

  # Return ----
  return(list(obj = obj, data = dat, fit = fit, report = report))
}
