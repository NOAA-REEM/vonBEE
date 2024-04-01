#'
#' ============================================
#' Mixed effect temperature specific VonB  (sensu Holsman et al. 2016)
#' Kirstin Holsman
#' 2023
#' ============================================
#' 


  source("R/make.R")

  # LOAD data
  # ------------------------------------
  load("data/in/NotShared/2023_09_26/LWA_ALL.Rdata")
    dat.all <- LWA_ALL%>%filter(!is.na(AGE), !is.na(WEIGHT),WEIGHT>0, AGE>0,!is.na(GEAR_TEMPERATURE),!is.na(SPECIES_NAME))%>%
      mutate(YEAR  = as.numeric(substr(START_TIME,1,4)),
             MONTH = as.numeric(substr(START_TIME,6,7)),
             DAY   = as.numeric(substr(START_TIME,9,10)), START_DATE = as.Date(substr(START_TIME,1,10)))%>%
      filter(CRUISE_TYPE =="Race_Groundfish")%>%
      select(REGION,YEAR,CRUISE_TYPE,HAULJOIN,STATIONID, STRATUM,age=AGE,weight=WEIGHT,Temp = GEAR_TEMPERATURE, SPECIES_CODE,CN = COMMON_NAME, SN = SPECIES_NAME,CRUISE_TYPE)%>%
      mutate(Sp = factor(CN))
  
    species <- (unique(dat.all$Sp))
    species <- c("walleye pollock","Pacific cod", "arrowtooth flounder","yellowfin sole","sablefish","Pacific ocean perch")
    regions <- unique(dat.all$REGION)

    # TMB Model for W @ A:
    #-----------------------------------
    require(TMB) 
    # precompile()
    tmpdir <- getwd()
    src_dir<- "src/TMB/models/vonBTv1/"
    setwd(src_dir)
    #TMB::compile('VonBTv1.cpp',"-o0 -g") 
    TMB::compile('VonBTv1.cpp') 
    #recompile("VonBTv1")
    dyn.load(dynlib('VonBTv1')) 
    
    # test the TMB model
    source("../../../../R/sub_scripts/test_TMB_model.R")
    
    for(regIN in regions){
      out_fn <- file.path("../../../../data/out",regIN)
      if(!dir.exists(out_fn))
        dir.create(out_fn)
      
      for(spIN in species[1:4]){  # skip sablefish and POP for now
        out_fn <- file.path("../../../../data/out",regIN,spIN)
        if(!dir.exists(out_fn))
          dir.create(out_fn)
        
        plot_fn<- file.path(out_fn,"Figs")
        if(!dir.exists(plot_fn))
          dir.create(plot_fn)
        
        silent_tmp <- TRUE # change to FALSE to see print out
        source("../../../../R/sub_scripts/run_model_set.R")
        
        # run aic on models
        nms    <- (names(models[[1]]))
        fits   <- lapply(models[-length(models)],"[[",which(nms =="fit"))
        if( any(is.na(fits)) ) models <- models[-which(is.na(fits))]
        
        nll    <- rep(0, length(models)-1)
        
        for(i in 1:(length(models)-1))
          nll[i] <- models[[i]]$objFun
        fits   <- lapply(models[-length(models)],"[[",which(nms =="fit"))
        input  <- lapply(models[-length(models)],"[[",which(nms =="input"))
        fitted <- r2 <- list()
        for(i in 1:length(models[-length(models)])){
          fitted[[names(models)[i]]] <- tmpfit <- models[[i]]$fitted
          r2[[names(models)[i]]]     <- cor(tmpfit$W_hat, tmpfit$W_obs,method ="pearson")^2
        }

        estpar <- lapply(fits,"[[",which(names(fits[[1]])=="par"))
        convergence <- lapply(fits,"[[",which(names(fits[[1]])=="convergence"))
        convergence_optim <- rep(NA,length(models)-1)
        for(i in 1:(length(models)-1))
          if(any(names(models[[i]])=="convergence_test"))
            convergence_optim[i]<- models[[i]]$convergence_test[2]
     
        message <- lapply(fits,"[[",which(names(fits[[1]])=="message"))
        
        nobs   <- lapply(input,"[[",which(names(input[[1]])=="nobs"))
        
        aic_table <- data.frame(model   = names(fits),
                                r2      = as.numeric(r2),
                                nll     = as.numeric(nll), 
                                k       = lengths(estpar),
                                n       = as.numeric(nobs),
                                converg_nlminb = as.numeric(convergence),
                                convergence_optim = as.numeric(convergence_optim),
                                message_nlminb = unlist(message))
        aic_table<- aic_table%>%mutate(aic = aic(k, nll),
                           aicc = aicc(k, nll, n),
                           delta_aic = aicc-min(aicc), 
                           aic_L= exp(-delta_aic/2))%>%arrange(delta_aic)
        top_model <- aic_table$aic_L[1]
        aic_table$sum_aic_L <- sum(aic_table$aic_L)
        aic_table<- aic_table%>%
          mutate(aic_wt = aic_L/sum_aic_L, 
                 sum_aic_wt = round( cumsum(aic_wt),4),
                 evidence_wt_top = round(top_model/aic_L,4), top_set = aic_wt>0.9 )
        # create mean weight at age table
        top_model <- models[[aic_table$model[1]]]
        top_model$fit$par
        top_model_par <- fit_to_list(top_model$fit$par)
        missing<- (names(par)%in%names(top_model_par))
        if(any(missing==FALSE))
          for(p in names(par)[which(!missing)])
            top_model_par[[p]]<-par[[p]]
        
          tmp_age <- (1:25)
          Temp    <- mean(dat$Temp,na.rm=T)
          newData <- list(
            YrIndex = tmp_age*0,
            age     = tmp_age,
            weight  = tmp_age*0,
            cohortYr = tmp_age*0,
            ncov     = 4,
            nobs     = length(tmp_age),
            nyrs     = 1,
            d_offset = 0,
            vonb_cov  = matrix(rep(c(0,Temp,Temp^2,Temp),length(tmp_age)),4,length(tmp_age)),
            sdvonb_cov = matrix(rep(c(0,Temp,Temp^2,Temp),length(tmp_age)),4,length(tmp_age))*0)
          
          
          tmp      <- sim_dat (par = top_model_par,data = newData)
          mn_wtAge <- data.frame(age = tmp$age, weight = tmp$What, Temp = Temp)
          years    <- sort(unique(dat$YEAR))
          nyrs     <- length(years)
          for(y in 1:nyrs){
            sub <- dat%>%filter(YEAR==years[y])
            Temp    <- mean(sub$Temp,na.rm=T)
            newData <- list(
              YrIndex = tmp_age*0,
              age     = tmp_age,
              weight  = tmp_age*0,
              cohortYr = tmp_age*0,
              ncov     = 4,
              nobs     = length(tmp_age),
              nyrs     = 1,
              d_offset = 0,
              vonb_cov  = matrix(rep(c(0,Temp,Temp^2,Temp),length(tmp_age)),4,length(tmp_age)),
              sdvonb_cov = matrix(rep(c(0,Temp,Temp^2,Temp),length(tmp_age)),4,length(tmp_age))*0)
            tmp      <- sim_dat (par = top_model_par,data = newData)
            if(y ==1)
            mn_wtAge_yr <- data.frame(Year = years[y],age = tmp$age, weight = tmp$What, Temp = Temp)
            if(y>1)
              mn_wtAge_yr <- rbind(mn_wtAge_yr,data.frame(Year = years[y],age = tmp$age, weight = tmp$What, Temp = Temp))
          }
          save(mn_wtAge,file=file.path(out_fn,"mn_wtAge.Rdata"))
          save(mn_wtAge_yr,file=file.path(out_fn,"mn_wtAge_yr.Rdata"))
          save(top_model_par, file=file.path(out_fn,"top_model_par.Rdata"))
          save(aic_table, file=file.path(out_fn,"aic_table.Rdata"))
          write.csv(reshape2::acast(mn_wtAge_yr%>%select(-Temp), Year ~ age)/1000,file=file.path(out_fn,"wt_age_yr_kg.csv"))
          write.csv(aic_table,file=file.path(out_fn,"aic_table.csv"))
          rm(list = rmlist)
      }
    }
    
    setwd(tmpdir)
    
    

    
        
    