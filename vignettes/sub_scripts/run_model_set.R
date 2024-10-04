#'
#'run_model_set.R
#' this code loops through the model set for each species
#' K Holsman Aug 2024
#' 
#'


# A. Run model set:
#1. Run null model
#2. Run model with SST and SST^2
#3. Run model with BT and BT^2
#4. Run model with RFR
#5. Run model with subset of covariates for EBS

# B. Run AIC selection

# C. Select top model for use in CEATTLE

  #source("R/make.R")       # loads packages, data, setup, etc.

  print(cat("now running models and plots for: ",spIN,"in Region", regIN, "\n"))
 
  # select subset for EBS pollock
  dat <- dat.all%>%filter(Sp==spIN,REGION==regIN)%>%
    mutate(cohortYr   = YEAR-age,
           cohortYrfac = (factor(YEAR-age)))
  
  cov_dat <- dat%>%select(Temp, YEAR)%>%group_by(YEAR)%>%
    summarize(mnTemp = mean(Temp,na.rm=T),
              seTemp = se (Temp,na.rm=T))%>%
    mutate(factYr = factor(YEAR),
           YrIndex = as.numeric(factYr))%>%
    relocate(YrIndex,YEAR, mnTemp, seTemp,factYr)%>%arrange(YrIndex)%>%
    ungroup()
  
  subNA <- function(x,sub =0){
    if(any(is.na(x)))
      x[is.na(x)] <-sub
    return(x)
  }
  cohort_cov_dat <- dat%>%select(Temp, cohortYr)%>%
    group_by(cohortYr)%>%
    summarize(chort_mnTemp = mean(Temp,na.rm=T),
              chort_seTemp = se(Temp,na.rm=T))%>%
    mutate_at(c("chort_mnTemp","chort_seTemp"),~subNA(.x))%>%
    ungroup()
   
  
  dat     <- dat%>%left_join(cov_dat)%>%left_join(cohort_cov_dat)
  
  # add in diet covars too:
  #load()
  
  
  covs    <- as.matrix(data.frame(blank=0,mnTemp=dat$mnTemp,mnTemp2=dat$mnTemp^2,chort_mnTemp=dat$chort_mnTemp))
  covs_se <- as.matrix(data.frame(blank=0,seTemp=dat$seTemp,seTemp2 =dat$seTemp^2,chort_seTemp = dat$chort_seTemp))
  ncov    <- dim(covs)[2]
  nyrs    <- dim(cov_dat)[1]
  nobs    <- dim(dat)[1]
  nyrs_cohort<- length(levels(dat$cohortYrfac))
  
  datIN <- list(YrIndex    = dat$YrIndex-1, 
                age        = dat$age, 
                weight     = dat$weight, 
                cohortYr   = as.numeric(dat$cohortYrfac)-1, 
                ncov       = ncov,
                nobs       = nobs,
                nyrs       = nyrs,
                d_offset   = 0 ,#0.01,  # prevents WInf --> 0
                #covYr      = cov_dat$YrIndex,
                vonb_cov   = matrix(t(covs),ncov,nobs),
                sdvonb_cov = matrix(t(covs_se),ncov,nobs))
  
  save(dat,file = file.path(out_fn,"dat.Rdata"))
  save(covs,file = file.path(out_fn,"covs.Rdata"))
  save(covs_se,file = file.path(out_fn,"covs_se.Rdata"))
  
  # define starting par
  par <-list(
    logH   = log(12.0),
    logK   = log(0.15),
    t0  = -0.25,
    mu_d = -.005,
    log_sigma_proc   = -.6,
    log_sigma_obs   = -.3,
    beta_c     = (c(0,0,0,0)),
    beta_1     = (0),
    beta_0     = (rep(0,nyrs_cohort)),
    u_y        = (rep(0,nyrs)))
  
  # create map for base model
  dlistIN_0 <- list(data = datIN, 
                    parameters  = par, 
                    random = NULL, #random = "u",
                    map = makeMap(param=par,estpar=list(
                      logH   = TRUE,
                      logK   = TRUE,
                      t0  = TRUE,
                      mu_d = TRUE,
                      log_sigma_proc  = FALSE,
                      log_sigma_obs   = TRUE,
                      beta_c     = FALSE,
                      beta_1     = FALSE,
                      beta_0     = FALSE,
                      u_y        = FALSE)))
  cat("  -- running m0 \n")
  m0 <- tryCatch({
    suppressWarnings(
    runMod(dlistIN    = dlistIN_0,
           src_fldr   = getwd(),
           version    = "VonBTv1",
           nrep       = 3,
           hessianIN  = TRUE,
           silentIN   = silent_tmp,
           se.fit     = TRUE,
           simulate   = FALSE,
           sim_nitr   = 1000,
           recompile  = FALSE,
           maxitr     = 10000,
           maxeval    = 10000))},
  error=function(cond) {
    message("model m4 failed")
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(list(fit = NA, objFun=NA, convergence_test=NA,fitted = data.frame(age=NA,W_hat=NA,W_obs=NA)))
  },
  warning=function(cond) {
    message("m4 caused a warning:")
    message("Here's the original warning message:")
    message(cond)
    return(NULL)
  },
  finally={
    message("model complete")
  })
  
  #covariate effects
  tmp_map <- makeMap(param=par,estpar=list(
    logH   = TRUE,
    logK   = TRUE,
    t0  = TRUE,
    mu_d = TRUE,
    log_sigma_proc  = FALSE,
    log_sigma_obs   = TRUE,
    beta_c     = TRUE,
    beta_1     = FALSE,
    beta_0     = FALSE,
    u_y        = FALSE))
  tmp_map$beta_c<-factor(c(NA,6,7,FALSE),levels=c(6,7))
  
  dlistIN_1 <- list(data = datIN, 
                    parameters  = par, 
                    random = NULL, #random = "u",
                    map = tmp_map)
  
  cat("  -- running m1 \n")
  m1 <- tryCatch({
    suppressWarnings(
    runMod(dlistIN    = dlistIN_1,
           src_fldr   = getwd(),
           version    = "VonBTv1",
           nrep       = 5,
           hessianIN  = TRUE,
           silentIN   = silent_tmp,
           se.fit     = TRUE,
           simulate   = FALSE,
           sim_nitr   = 1000,
           recompile  = FALSE,
           maxitr     = 50000,
           maxeval    = 50000))},
  error=function(cond) {
    message("model m4 failed")
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(list(fit = NA, objFun=NA, convergence_test=NA,fitted = data.frame(age=NA,W_hat=NA,W_obs=NA)))
  },
  warning=function(cond) {
    message("m4 caused a warning:")
    message("Here's the original warning message:")
    message(cond)
    return(NULL)
  },
  finally={
    message("model complete")
  })
  
  # covariates in Temp_y but also cohort yearTemp
  dlistIN_2<- list(data = datIN, 
                   parameters  = par, 
                   random = NULL, #random = "u",
                   map = makeMap(param=par,estpar=list(
                     logH   = TRUE,
                     logK   = TRUE,
                     t0  = TRUE,
                     mu_d = TRUE,
                     log_sigma_proc  = FALSE,
                     log_sigma_obs   = TRUE,
                     beta_c     = TRUE,
                     beta_1     = FALSE,
                     beta_0     = FALSE,
                     u_y        = FALSE)))
  cat("  -- running m2 \n")
  m2 <- tryCatch({
    suppressWarnings(
    runMod(dlistIN    = dlistIN_2,
           src_fldr   = getwd(),
           version    = "VonBTv1",
           nrep       = 5,
           hessianIN  = TRUE,
           silentIN   = silent_tmp,
           se.fit     = TRUE,
           simulate   = FALSE,
           sim_nitr   = 1000,
           recompile  = FALSE,
           maxitr     = 50000,
           maxeval    = 50000))},
  error=function(cond) {
    message("model m4 failed")
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(list(fit = NA, objFun=NA, convergence_test=NA,fitted = data.frame(age=NA,W_hat=NA,W_obs=NA)))
  },
  warning=function(cond) {
    message("m4 caused a warning:")
    message("Here's the original warning message:")
    message(cond)
    return(NULL)
  },
  finally={
    message("model complete")
  })
  
  # covariates in Temp_y but also cohort yearTemp and fixed effect of cohort year
  dlistIN_3<- list(data = datIN, 
                   parameters  = par, 
                   random = "beta_0", #random = "u",
                   map = makeMap(param=par,estpar=list(
                     logH   = TRUE,
                     logK   = TRUE,
                     t0  = TRUE,
                     mu_d = TRUE,
                     log_sigma_proc  = FALSE,
                     log_sigma_obs   = TRUE,
                     beta_c     = TRUE,
                     beta_1     = FALSE,
                     beta_0     = TRUE,
                     u_y        = FALSE)))
  dlistIN_3$random<-NULL
  cat("  -- running m3 \n")
  m3 <- tryCatch({
    suppressWarnings(
    runMod(dlistIN    = dlistIN_3,
           src_fldr   = getwd(),
           version    = "VonBTv1",
           nrep       = 5,
           hessianIN  = TRUE,
           silentIN   = silent_tmp,
           se.fit     = TRUE,
           simulate   = FALSE,
           sim_nitr   = 1000,
           recompile  = FALSE,
           maxitr     = 50000,
           maxeval    = 50000))},
  error=function(cond) {
    message("model m4 failed")
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(list(fit = NA, objFun=NA, convergence_test=NA,fitted = data.frame(age=NA,W_hat=NA,W_obs=NA)))
  },
  warning=function(cond) {
    message("m4 caused a warning:")
    message("Here's the original warning message:")
    message(cond)
    return(NULL)
  },
  finally={
    message("model complete")
  })
  
  # covariates in Temp_y but also cohort year Temp and ramdom effect of cohort year
  dlistIN_3$random<-"beta_0"
  cat("  -- running m3_rnd \n")
  m3_rnd <- tryCatch(
    {suppressWarnings(
    runMod(dlistIN    = dlistIN_3,
           src_fldr   = getwd(),
           version    = "VonBTv1",
           nrep       = NULL,
           hessianIN  = FALSE,
           silentIN   = silent_tmp,
           se.fit     = TRUE,
           simulate   = FALSE,
           sim_nitr   = 1000,
           recompile  = FALSE,
           maxitr     = 50000,
           maxeval    = 50000))},
  error=function(cond) {
    message("model m4 failed")
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(list(fit = NA, objFun=NA, convergence_test=NA,fitted = data.frame(age=NA,W_hat=NA,W_obs=NA)))
  },
  warning=function(cond) {
    message("m4 caused a warning:")
    message("Here's the original warning message:")
    message(cond)
    return(NULL)
  },
  finally={
    message("model complete")
  })
  
  # state space with auto regressive process error, no cohort yr effect
  dlistIN_4<- list(data = datIN, 
                   parameters  = par, 
                   random = "u_y", #random = "u",
                   map = makeMap(param=par,estpar=list(
                     logH   = TRUE,
                     logK   = TRUE,
                     t0  = TRUE,
                     mu_d = TRUE,
                     log_sigma_proc  = TRUE,
                     log_sigma_obs   = TRUE,
                     beta_c     = TRUE,
                     beta_1     = TRUE,
                     beta_0     = FALSE,
                     u_y        = TRUE)))
  
  cat("  -- running m4 \n")
  m4 <- tryCatch(
    {
    suppressWarnings(
    runMod(dlistIN    = dlistIN_4,
           src_fldr   = getwd(),
           version    = "VonBTv1",
           nrep       = NULL,
           hessianIN  = FALSE,
           silentIN   = silent_tmp,
           se.fit     = TRUE,
           simulate   = FALSE,
           sim_nitr   = 1000,
           recompile  = FALSE,
           maxitr     = 50000,
           maxeval    = 50000))
  },
  error=function(cond) {
    message("model m4 failed")
    message("Here's the original error message:")
    message(cond)
    # Choose a return value in case of error
    return(list(fit = NA, objFun=NA, convergence_test=NA,fitted = data.frame(age=NA,W_hat=NA,W_obs=NA)))
  },
  warning=function(cond) {
    message("m4 caused a warning:")
    message("Here's the original warning message:")
    message(cond)
    return(NULL)
  },
  finally={
    message("model complete")
  }
  )   
  
  # state space with auto regressive process error, no cohort yr effect
  dlistIN_5<- list(data = datIN, 
                   parameters  = par, 
                   random = c("u_y","beta_0"), #random = "u",
                   map = makeMap(param=par,estpar=list(
                     logH   = TRUE,
                     logK   = TRUE,
                     t0  = TRUE,
                     mu_d = TRUE,
                     log_sigma_proc  = TRUE,
                     log_sigma_obs   = TRUE,
                     beta_c     = TRUE,
                     beta_1     = TRUE,
                     beta_0     = TRUE,
                     u_y        = TRUE)))
  
  cat("  -- running m5 \n")
  m5 <- tryCatch(
    {
      suppressWarnings(
        runMod(dlistIN    = dlistIN_5,
               src_fldr   = getwd(),
               version    = "VonBTv1",
               nrep       = NULL,
               hessianIN  = FALSE,
               silentIN   = silent_tmp,
               se.fit     = TRUE,
               simulate   = FALSE,
               sim_nitr   = 1000,
               recompile  = FALSE,
               maxitr     = 50000,
               maxeval    = 50000))
    },
    error=function(cond) {
      message("model m4 failed")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(list(fit = NA, objFun=NA, convergence_test=NA,fitted = data.frame(age=NA,W_hat=NA,W_obs=NA)))
    },
    warning=function(cond) {
      message("m4 caused a warning:")
      message("Here's the original warning message:")
      message(cond)
      return(NULL)
    },
    finally={
      message("model complete")
    }
  )   
  cat("  -- running m5 (slow) \n")
  # state space with auto regressive process error, no cohort yr effect; no covars
  dlistIN_0_rnd<- list(data = datIN, 
                   parameters  = par, 
                   random = "u_y", #random = "u",
                   map = makeMap(param=par,estpar=list(
                     logH   = TRUE,
                     logK   = TRUE,
                     t0  = TRUE,
                     mu_d = TRUE,
                     log_sigma_proc  = TRUE,
                     log_sigma_obs   = TRUE,
                     beta_c     = FALSE,
                     beta_1     = TRUE,
                     beta_0     = FALSE,
                     u_y        = TRUE)))
  m0_rnd <- tryCatch({
    suppressWarnings(
    runMod(dlistIN    = dlistIN_0_rnd,
           src_fldr   = getwd(),
           version    = "VonBTv1",
           nrep       = NULL,
           hessianIN  = FALSE,
           silentIN   = silent_tmp,
           se.fit     = TRUE,
           simulate   = FALSE,
           sim_nitr   = 1000,
           recompile  = FALSE,
           maxitr     = 50000,
           maxeval    = 50000))},
    error=function(cond) {
      message("model m4 failed")
      message("Here's the original error message:")
      message(cond)
      # Choose a return value in case of error
      return(list(fit = NA, objFun=NA, convergence_test=NA,fitted = data.frame(age=NA,W_hat=NA,W_obs=NA)))
    },
    warning=function(cond) {
      message("m4 caused a warning:")
      message("Here's the original warning message:")
      message(cond)
      return(NULL)
    },
    finally={
      message("model complete")
    })
 
  pdat <- rbind(m0$fitted%>%mutate(type ="m0: base",model="m0" ),
                m1$fitted%>%mutate(type ="m1: Temp+Temp2",model="m1"  ),
                m2$fitted%>%mutate(type ="m2: Temp+Temp2; cohort yr Temp" ,model="m2" ),
                m3$fitted%>%mutate(type ="m3: Temp+Temp2; cohort yr Temp; cohort yr effect" ,model="m3" ),
                m3_rnd$fitted%>%mutate(type ="m3_rnd: Temp+Temp2; cohort yr Temp; rnd cohort yr effect" ,model="m3_rnd" ),
                m4$fitted%>%mutate(type ="m4: Temp+Temp2; cohort yr Temp; no cohort yr; proc error state space",model="m4"  ),
                m5$fitted%>%mutate(type ="m5: Temp+Temp2; cohort yr Temp; cohort yr; proc error state space",model="m5"  ),
                m0_rnd$fitted%>%mutate(type ="m0_rnd: proc error state space only (no covars)",model="m0_rnd"  ))
  
  # now plot model results
  plots<-list()
  cat("  -- plotting results \n")
  p1_models <- ggplot()+
    geom_point(data=m0$fitted,aes(x=age,y=W_obs,color="observed"),alpha=.3)+
    geom_point(data=pdat,aes(x=age,y=W_hat,color=type))+
    facet_wrap(model~.,nrow = 2)+
    scale_colour_manual( values = c(col_ramp(9)[c(1,3,4,5,6,7,8,9)],"gray"))+
    labs(title = spIN)+
    theme(plot.title = element_text(size = rel(.8)))+
    theme_minimal()
  
  p1_models
  plots[["p1_models"]]<-p1_models
  
  sclr<-1
  jpeg(filename = file.path(plot_fn,"model_plots.jpg"), res = 350, width = sclr*10, height = sclr*6,units = "in")
  print(p1_models)
  dev.off()
  
  rm(age)
  age     = (1:25)
  
  newData <- list(
    YrIndex = age*0,
    age     = age,
    weight  = age*0,
    cohortYr = age*0,
    ncov     = 4,
    nobs     = length(age),
    nyrs     = 1,
    d_offset = 0,
    vonb_cov  = matrix(rep(c(0,1,1,1),length(age)),4,length(age)),
    sdvonb_cov = matrix(rep(c(0,1,1,1),length(age)),4,length(age))*0)
  
  fit_to_list <- function(modelpar){
    out <- list()
    for(nm in names(modelpar))
      out[[nm]]<- as.numeric(modelpar[names(modelpar)==nm])
    return(out)
  }
  
  # if(!is.na(m5$fit)){
  #   Temp_m_par <- fit_to_list(m5$fit$par)
  #   Temp_m_par$beta_0<-par$beta_0
  #   Temp_m_par$beta_1<-par$beta_1
  #   Temp_m_par$u_y<-par$u_y
  # }else{
    Temp_m_par <- fit_to_list(m4$fit$par)
    for(p in names(par)){
      if(!any(names(Temp_m_par)%in%p)){
        Temp_m_par[[p]]<-par[[p]]
      }
    }
# 
#   }

  
  tt <- 0
  Tlist <- c(-4,-2,0,2,4,6)
  for(Temp in seq(-5,10,.2) ){
    tt <- tt+1
    
    newData <- list(
      YrIndex = age*0,
      age     = age,
      weight  = age*0,
      cohortYr = age*0,
      ncov     = 4,
      nobs     = length(age),
      nyrs     = 1,
      d_offset = 0,
      vonb_cov  = matrix(rep(c(0,Temp,Temp^2,Temp),length(age)),4,length(age)),
      sdvonb_cov = matrix(rep(c(0,Temp,Temp^2,Temp),length(age)),4,length(age))*0)
    
    
    tmp <- sim_dat (par = Temp_m_par,data = newData)
    tmpout <- tmp$out
    colnames(tmpout)[-(1:3)]<-c("blank","Temp","Temp2","CohortT")
    tmpout <- tmpout%>%mutate(type =paste(Temp," deg" ))
    if(tt ==1)
      pdat_out <- tmpout
    if(tt!=1)
      pdat_out<- rbind(pdat_out,tmpout)
    rm(tmpout)
    rm(newData)
    rm(tmp)
    
  }
  
  tt <- 0
  Tlist <- c(-4,-2,0,2,4,6)
  for(Temp in seq(-5,10,.2) ){
    tt <- tt+1
    
    newData <- list(
      YrIndex = age*0,
      age     = age,
      weight  = age*0,
      cohortYr = age*0,
      ncov     = 4,
      nobs     = length(age),
      nyrs     = 1,
      d_offset = 0,
      vonb_cov  = matrix(rep(c(0,0,0,Temp),length(age)),4,length(age)),
      sdvonb_cov = matrix(rep(c(0,0,0,Temp),length(age)),4,length(age))*0)
    
    
    tmp <- sim_dat (par = Temp_m_par,data = newData)
    tmpout <- tmp$out
    colnames(tmpout)[-(1:3)]<-c("blank","Temp","Temp2","CohortT")
    tmpout <- tmpout%>%mutate(type =paste(Temp," deg" ))
    if(tt ==1)
      pdat_out2 <- tmpout
    if(tt!=1)
      pdat_out2<- rbind(pdat_out2,tmpout)
    rm(tmpout)
    rm(newData)
    rm(tmp)
    
  }
  
  tt <- 0
  Tlist <- c(-4,-2,0,2,4,6)
  for(Temp in seq(-5,10,.2) ){
    tt <- tt+1
    
    newData <- list(
      YrIndex = age*0,
      age     = age,
      weight  = age*0,
      cohortYr = age*0,
      ncov     = 4,
      nobs     = length(age),
      nyrs     = 1,
      d_offset = 0,
      vonb_cov  = matrix(rep(c(0,Temp,Temp^2,0),length(age)),4,length(age)),
      sdvonb_cov = matrix(rep(c(0,0,0,0),length(age)),4,length(age))*0)
    
    
    tmp     <- sim_dat (par = Temp_m_par,data = newData)
    tmpout  <- tmp$out
    colnames(tmpout)[-(1:3)]<-c("blank","Temp","Temp2","CohortT")
    tmpout <- tmpout%>%mutate(type =paste(Temp," deg" ))
    if(tt ==1)
      pdat_out3 <- tmpout
    if(tt!=1)
      pdat_out3<- rbind(pdat_out3,tmpout)
    rm(tmpout)
    rm(newData)
    rm(tmp)
    
  }
  
  T_models <- ggplot()+
    geom_point(data=m0$fitted,aes(x=age,y=W_obs,color="observed"),alpha=.3)+
    geom_point(data=pdat_out%>%filter(Temp%in%Tlist),aes(x=age,y=W_hat,color=type))+
    scale_colour_manual( values = c(col_ramp(8)[c(1,3,4,5,6,7)],"gray"))+
    labs(title = spIN)+
    theme(plot.title = element_text(size = rel(.8)))+
    theme_minimal()
  
  plots[["T_models"]]<-T_models

  pdat_out3$type2 <- "fitted (T_y)"
  pdat_out$type2  <- "fitted (T_y+cohort_T)"
  
  ageIN <-1:20
  T_models2 <- ggplot()+
    geom_point(data=dat%>%filter(age%in%ageIN),aes(x=Temp,y=log(weight),color="observed"),alpha=.3)+
    geom_point(data=pdat_out3%>%filter(age%in%ageIN, Temp<6,Temp>-3),aes(x=Temp,y=log(W_hat),color=type2))+
    #geom_point(data=pdat_out2%>%filter(age%in%ageIN, Temp<6,Temp>-3),aes(x=CohortT,y=W_hat,color="fitted (cohort_T)"))+
    geom_point(data=pdat_out%>%filter(age%in%ageIN, Temp<6,Temp>-3),aes(x=Temp,y=log(W_hat),color=type2))+
    scale_colour_manual( values = c(col_ramp(8)[c(5,7)],"gray"))+
    facet_wrap(age~.,ncol=2, scales="free_y")+
    labs(title = spIN)+
    theme(plot.title = element_text(size = rel(.8)))+
    theme_minimal()
  
  plots[["T_models2"]]<-T_models2
  
  T_models2
  
  ageIN <-1:20
  T_models3 <- ggplot()+
    geom_point(data=dat%>%filter(age%in%ageIN),aes(x=Temp,y=log(weight)),color="gray",alpha=.3)+
    geom_point(data=pdat_out3%>%filter(age%in%ageIN, Temp<6,Temp>-3),aes(x=Temp,y=log(W_hat),color=factor(age,levels=1:20)))+
    #geom_point(data=pdat_out2%>%filter(age%in%ageIN, Temp<6,Temp>-3),aes(x=CohortT,y=W_hat,color="fitted (cohort_T)"))+
    geom_point(data=pdat_out%>%filter(age%in%ageIN, Temp<6,Temp>-3),aes(x=Temp,y=log(W_hat),color=factor(age,levels=1:20)))+
    scale_colour_manual( values = c(col_ramp(20),"gray"))+
    facet_wrap(type2~.,ncol=2, scales="free_y")+
    labs(title = spIN)+
    theme(plot.title = element_text(size = rel(.8)))+
    theme_minimal()
  
  plots[["T_models3"]]<-T_models3
  T_models3
  
  sclr<-1
  jpeg(filename = file.path(plot_fn,"model_Temp.jpg"), res = 350, width = sclr*6, height = sclr*5,units = "in")
  print(T_models)
  dev.off()
  
  sclr<-1
  jpeg(filename = file.path(plot_fn,"model_Temp_byage.jpg"), res = 350, width = sclr*9, height = sclr*10,units = "in")
  print(T_models2)
  dev.off()
  
  sclr<-1.5
  jpeg(filename = file.path(plot_fn,"model_Temp_byage2.jpg"), res = 350, width = sclr*6, height = sclr*4,units = "in")
  print(T_models3)
  dev.off()
  
  models <- list(m0=m0,m1=m1,m2=m2,m3=m3,m3_rnd=m3_rnd,m4=m4,m5=m5,species=spIN)
  save(models,file=file.path(out_fn,"models.Rdata"))
  save(plots,file=file.path(out_fn,"plots.Rdata"))
  rmlist <- c("models","plots","T_models3","T_models2","T_models","datIN","pdat_out2",
              "pdat_out3","pdat_out","pdat","m0","m1","m2","m3","m3_rnd","m4","m5")
  print("run_Model_set complete\n ----------------------------")
