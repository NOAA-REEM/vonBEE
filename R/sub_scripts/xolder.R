
mod <- m4 
sim_4 <- sim_dat(par = list(
  logH   = mod$sdreport$value["logH"],
  logK   = mod$sdreport$value["logK"],
  t0  = mod$sdreport$value["t0"],
  mu_d = mod$sdreport$value["mu_d"],
  log_sigma_proc   = mod$sdreport$value["log_sigma_proc"],
  log_sigma_obs   = mod$sdreport$value["log_sigma_obs"],
  beta_c     =  mod$sdreport$value[names(mod$sdreport$value)=="beta_c"],
  beta_1     = 0,
  beta_0     = rep(length(unique(dlistIN_4$data$cohortYr))),
  u_y        = rep(0,nyrs)),data = dlistIN_4$data)
names(simd_T2$out)[-(1:3)]<-paste0("cov",names(simd_T2$out)[-(1:3)])

x<-1:20  
y <-0.063*x+-0.0094*x^2
plot(x,y)
ggplot()+
  geom_point(data=simd$out,aes(x=age,y=W_obs,color="observed"))+
  geom_point(data=simd_T2$out,aes(x=age,y=W_hat,color=factor(round(cov2))))+
  geom_point(data=simd$out,aes(x=age,y=W_hat,color="predicted"))



# run base model without covars:
model <- MakeADFun(data       =  dlistIN_0$data,
                   parameters =  dlistIN_0$parameters,
                   random     =  dlistIN_0$random,
                   DLL        =  'VonBTv1',
                   map        =  dlistIN_0$map, 
                   hessian    =  TRUE)

hessian$hessian <- FALSE
fit <- nlminb(model$env$last.par.best,model$fn,model$gr) 
for(i in 1:3)
  fit <- nlminb(model$env$last.par.best,model$fn,model$gr) 



rep <- sdreport(model)

# extract estimated process:
u <- summary(rep, "random")[, "Estimate"]
u_se <- summary(rep, "random")[, "Std. Error"]
# extract fixed effects:
fixed <- summary(rep, "fixed")


m0 <- suppressWarnings(
  runMod(dlistIN=dlistIN_0,
         src_fldr   = getwd(),
         version    = "VonBTv1",
         nrep       = 3,
         hessian    = FALSE,
         silentIN   = TRUE,
         se.fit     = TRUE,
         simulate   = FALSE,
         sim_nitr   = 1000,
         recompile  = FALSE,
         maxitr     = 10000,
         maxeval    = 10000))



tmpmod  <- list()
tmpmod$model    <-  model
tmpmod$fit      <-  fit
if(tmpmod$fit$objective==Inf) 
  stop("Problem with objective function (Inf)")
if(is.na(tmpmod$fit$objective)) 
  stop("Problem with objective function (NaN)")
tmpmod$objFun   <-  model$fn()
tmpmod$report   <-  model$report()
tmpmod$sdreport <-  sdreport(model,getJointPrecision=TRUE)
tmpmod$rep      <-  sdreport(model,getJointPrecision=TRUE)
tmpmod$input    <-  dlistIN$data
tmpmod$out      <-  summary(tmpmod$rep)
tmpmod$mle      <-  model$env$last.par.best
tmpmod$Hessian  <-  with(model,optimHess(tmpmod$mle,model$fn,model$gr))
tmpmod$LnDet    <-  sum(determinant(tmpmod$Hessian, logarithm=TRUE)$modulus[[1]])
lp              <-  model$env$last.par


rep       <- ADreport(model)
rep       <- sdreport(model,getJointPrecision=TRUE)
est       <- data.frame(t(model$env$last.par.best))
est$Tcoef <- 0
est       <- data.frame(t(data.frame(val=rep$value,sd=rep$sd)))



What_fun <- function( estIn = est, newdat = data.frame(age=seq(0,30,.1),Temp=0),randn = 1000){
  nobs<-length(newdat$age)
  logH<-estIn$logH
  logK<-estIn$logK
  t0<-estIn$t0
  Tcoef<-estIn$Tcoef
  mu_d<-estIn$mu_d
  logSigma<-estIn$logSigma
  mnd<-mnWinf<-mnlogWhat<-rep(0,nobs)
  mcmc_val<-NA
  if(is.na(randn)!=TRUE){
    lgmnd<-rnorm(randn,mu_d[1],mu_d[2])
    Tcoef_tmp<-rnorm(randn,Tcoef[1],Tcoef[2])
    H_tmp<-rnorm(randn,logH[1],logH[2])
    K_tmp<-rnorm(randn,logK[1],logK[2])
    t0_tmp<-rnorm(randn,t0[1],t0[2])
    d<-Winf<-logWhat<-matrix(0,nobs,randn)
  }
  
  for (i in 1:nobs)
  {
    mnd[i]			<-exp(mu_d[1]+Tcoef[1]*newdat$Temp[i]) # #mfexp(mu_d+Pcoef*Pval_lag1(i)+Tcoef*Temp(i));
    mnWinf[i]       <-(logH[1]/logK[1])^(1.0/(1.0 - mnd[i] ))
    mnlogWhat[i]    <-log(mnWinf[i]) + (1.0/(1.0 - mnd[i]))*log(1.0 - exp(-logK[1] * (1.0 - mnd[i]) * (newdat$age[i] - t0[1]) )) 
    
    if(is.na(randn)!=TRUE){
      
      d[i,]			<-exp(lgmnd+Tcoef_tmp*newdat$Temp[i]) # #mfexp(mu_d+Pcoef*Pval_lag1(i)+Tcoef*Temp(i));
      Winf[i,]       <-(H_tmp/K_tmp)^(1.0/(1.0 - d[i,] ))
      logWhat[i,]    <-log(Winf[i,]) + (1.0/(1.0 - d[i,]))*log(1.0 - exp(-K_tmp* (1.0 - d[i,]) * (newdat$age[i] - t0_tmp) )) 
      
    }
  }
  mn_val<-data.frame(d=mnd,Winf=mnWinf,logWhat=mnlogWhat)
  if(is.na(randn)!=TRUE) mcmc_val<-list(d=d,Winf=Winf,logWhat=logWhat)
  
  return(list(mn_val=mn_val,mcmc_val=mcmc_val))
  
}
tmpdat<-data.frame(age=seq(0,30,.1),Temp=0)
W_hat<-What_fun(estIn=est,newdat=tmpdat)
mn<-tapply(log(dat$weight),dat$age,mean)  # means of each age data
plot(dat$age,dat$weight,pch=16)
points(as.numeric(names(mn)),exp(mn),col="red",pch=16) # means
quantile(W_hat[[2]]$logWhat)

tmp<-apply(W_hat[[2]]$logWhat,1,quantile,prob=c(0,.05,.10,.50,.90,.95,1))

points(tmpdat$age,exp(tmp[4,]),col="gray",pch=16) # W at Age function
points(tmpdat$age,exp(tmp[2,]),col="gray",pch=16) # W at Age function
points(tmpdat$age,exp(tmp[6,]),col="gray",pch=16) # W at Age function

out<-summary(rep)
model$he()
cov2cor(solve(model$he()))
print(rep)
#-----------------------------------


