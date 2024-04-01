
################################################################################

## Mixed effect temperature specific VonB  (sensu Holsman et al. 2016)
## Kirstin Holsman
## 2016
################################################################################

rm(list=ls())
setwd("/Users/kkari/GitHub/TMB_models/KirVonB")
load("LWA_ALL.Rdata")
dat.all<-na.omit(data.frame(
	age=LWA_ALL$AGE,
	weight=LWA_ALL$WEIGHT,
	Temp=LWA_ALL$GEAR_TEMPERATURE,
	Sp=factor(LWA_ALL$GOAPOLL_PRED)))
	# R=factor(LWA_ALL$REGION),
	# Sp=factor(LWA_ALL$GOAPOLL_PRED),
	# Sex=factor(LWA_ALL$SEX),
	# L=LWA_ALL$LENGTH,
	# )

dat<-dat.all[which(dat.all$Sp=="WALLEYE POLLOCK"),]
# dat<-dat[dat$age<=25,]
	par<-list(
		H=9.0,
		K=0.5,
		t0=-0.25,
		log_mean_d= -.3,
		logSigma= -.2,
		Tcoef=0)
	mapp<-list(
		H=factor("H"),
		K=factor("K"),
		t0=factor("T0"),
		log_mean_d=factor("log_mean_d"),
		logSigma= factor("logSigma"),
		Tcoef=factor(NA))

## TMB Model for W @ A:
#-----------------------------------
	require(TMB) 
	# precompile()
	compile('KirVonB.cpp',"-o0 -g") 
	dyn.load(dynlib('KirVonB')) 
	model<-MakeADFun(data=dat,parameters=par,DLL='KirVonB',map=mapp, hessian = TRUE)#,map=list(H=factor(1),K=factor(2),c=factor(NA))) 
	fit<-nlminb(model$env$last.par.best,model$fn,model$gr) 
	for(i in 1:3)
		fit<-nlminb(model$env$last.par.best,model$fn,model$gr) 
	
	rep<-sdreport(model,getJointPrecision=TRUE)
	est<-data.frame(t(model$env$last.par.best))
	est$Tcoef<-0
	est<-data.frame(t(data.frame(val=rep$value,sd=rep$sd)))
	What_fun<-function(estIn=est,newdat=data.frame(age=seq(0,30,.1),Temp=0),randn=1000){
		nobs<-length(newdat$age)
		H<-estIn$H
		K<-estIn$K
		t0<-estIn$t0
		Tcoef<-estIn$Tcoef
		log_mean_d<-estIn$log_mean_d
		logSigma<-estIn$logSigma
		mnd<-mnWinf<-mnlogWhat<-rep(0,nobs)
		mcmc_val<-NA
 		if(is.na(randn)!=TRUE){
				lgmnd<-rnorm(randn,log_mean_d[1],log_mean_d[2])
				Tcoef_tmp<-rnorm(randn,Tcoef[1],Tcoef[2])
				H_tmp<-rnorm(randn,H[1],H[2])
				K_tmp<-rnorm(randn,K[1],K[2])
				t0_tmp<-rnorm(randn,t0[1],t0[2])
				d<-Winf<-logWhat<-matrix(0,nobs,randn)
		}

		for (i in 1:nobs)
	    {
	      mnd[i]			<-exp(log_mean_d[1]+Tcoef[1]*newdat$Temp[i]) # //mfexp(log_mean_d+Pcoef*Pval_lag1(i)+Tcoef*Temp(i));
	      mnWinf[i]       <-(H[1]/K[1])^(1.0/(1.0 - mnd[i] ))
	      mnlogWhat[i]    <-log(mnWinf[i]) + (1.0/(1.0 - mnd[i]))*log(1.0 - exp(-K[1] * (1.0 - mnd[i]) * (newdat$age[i] - t0[1]) )) 
	    
		    if(is.na(randn)!=TRUE){
			 
		      d[i,]			<-exp(lgmnd+Tcoef_tmp*newdat$Temp[i]) # //mfexp(log_mean_d+Pcoef*Pval_lag1(i)+Tcoef*Temp(i));
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


