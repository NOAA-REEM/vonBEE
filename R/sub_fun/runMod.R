#' run the vonBT()  model
#'
#' runMod() will run the vonBT() recruitment model
#' @import vonBT
#' @email For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#' @weblink 
#' @param optimizer default is "nlminb"; can also be set to "optim"
#' @param methodIN  default is NULL, can be set to "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", or "Brent"
#' @param convergenceTest Test for failed convergence by also running the model using optim
#' @param dlistIN   list of input data and parameters used for the makeADfun() dependency
#' * dlistIN$parameters  list of parameters for the TMB vonBT model:
#' * dlistIN$randonm   list of which parms should be radnom effects
#' * dlistIN$data      list of data to read into the model
#' * dlistIN$map       list of the mapp for the TMB model
#' * dlistIN$phases    (defunct) phases for estimating parameters
#' * dlistIN$inputs    archive of intputs for this function
#' @param version     recruitment model to run. Default is 'vonBT'
#' @param src_fldr    subfolder where model is located. Default is 'src'
#' @param silentIN    T/F , default = 'TRUE'
#' @param se.fit      T/F , default = 'TRUE'
#' @param simulate    T/F , simulate the model? default = 'FALSE'
#' @param sim_nitr    number of iterations for the simulate function
#' @param recompile   recompile the .cpp file? default is 'FALSE'
#' @param maxitr      10000  Max iterations for fitting the objective function
#' @param maxeval     10000  Max evaluations for fitting the objective function
#' @return returns  summary of the model including the mle
#' 
#' @examples
#' datlist <- readMake_vonBT_data("data/in/vonBT_Inputs.xlsx" )
#' mm      <- runRecMod(dlistIN   = datlist, version   = 'vonBT',recompile = FALSE,simulate  = TRUE,sim_nitr  = 1000)  
#' @export
runMod<-function(
    dlistIN,
    src_fldr   = "src",
    version    = "vonBT",
    nrep       = 3,
    optimizer  = "nlminb",
    convergenceTest  = TRUE,
    methodIN   = NULL, # method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent")
    hessianIN  = TRUE,
    silentIN   = TRUE,
    se.fit     = TRUE,
    simulate   = FALSE,
    sim_nitr   = 1000,
    recompile  = FALSE,
    maxitr     = 10000,
    maxeval    = 10000){
  
  wd0 <- getwd()
  if(!is.null(src_fldr))
    setwd(src_fldr)
  
  if(recompile){  
    if(file.exists(paste0(version,".o")))
      file.remove(paste0(version,".o"))
    if(file.exists(paste0(version,".so")))
      file.remove(paste0(version,".so"))
    compile(paste0(version, ".cpp")) 
  }
  
  dyn.load( dynlib(version) ) 
  
  model  <- MakeADFun(
    data                 =  dlistIN$data,
    parameters           =  dlistIN$parameters,
    random               =  dlistIN$random,
    method               =  methodIN,
    DLL                  =  version,
    checkParameterOrder  =  TRUE,
    hessian              =  hessianIN,
    map                  =  dlistIN$map,
    silent               =  silentIN)
  
  tmpmod  <- list()
  tmpmod$model    <-  model
  
  if(optimizer == "nlminb"){
   tmpmod$fit      <- nlminb(start = model$par, obj = model$fn, gr = model$gr, 
                              control=list(iter.max=maxitr,eval.max=maxeval), method = methodIN)
   if(!is.null(nrep))
     for(i in 1:nrep)
       tmpmod$fit      <-  nlminb(start   = model$env$last.par.best,
                                  obj     = model$fn,
                                  gr      = model$gr, 
                                  control = list(iter.max=maxitr,eval.max=maxeval), 
                                  method  = methodIN)
    if(tmpmod$fit$objective==Inf) 
       stop("Problem with objective function (Inf)")
    if(is.na(tmpmod$fit$objective)) 
       stop("Problem with objective function (NaN)")
     if(convergenceTest){
       if(tmpmod$fit$convergence!=0){
         cat(" problem with nlminb convergence != 0 \n model reports",tmpmod$fit$message, "\n trying optim() to test for convergence")
         model2  <- MakeADFun(
           data                 =  dlistIN$data,
           parameters           =  dlistIN$parameters,
           random               =  dlistIN$random,
           method               =  NULL,
           DLL                  =  version,
           checkParameterOrder  =  TRUE,
           hessian              =  TRUE,
           map                  =  dlistIN$map,
           silent               =  silentIN)
         tmpmod$fit_optim     <- optim(par    = model2$par,
                                 fn     = model2$fn, 
                                 gr     = model2$gr, 
                                 method = methodIN,
                                 control = list(maxit=maxitr))
         if(tmpmod$fit_optim$convergence==0)
           cat("... optim was able to converge; convergence == 0\n")
         if(tmpmod$fit_optim$convergence!=0)
           cat("...optim also failed; convergence = ",tmpmod$fit_optim$convergence,"\n" )
         tmpmod$convergence_test = c(nlminb = tmpmod$fit$convergence, optim = tmpmod$fit_optim$convergence)
      }
     }
   }
  if(optimizer == "optim"){
   tmpmod$fit     <- optim(par    = model$par,
                           fn     = model$fn, 
                           gr     = model$gr, 
                           method = methodIN)
   if(!is.null(nrep))
     for(i in 1:nrep)
       tmpmod$fit      <-  optim(par    = model$env$last.par.best,
                                 fn     = model$fn,
                                 gr     = model$gr, 
                                 method = methodIN)
   }
  
  
 
  tmpmod$objFun   <-  model$fn()
  tmpmod$report   <-  model$report()
  tmpmod$sdreport <-  sdreport(model,getJointPrecision=TRUE)
  tmpmod$rep      <-  sdreport(model,getJointPrecision=TRUE)
  tmpmod$input    <-  dlistIN$data
  tmpmod$out      <-  summary(tmpmod$rep)
  tmpmod$mle      <-  model$env$last.par.best
  tmpmod$Hessian  <- NULL
  if(hessianIN){
    tmpmod$Hessian  <-  with(model,optimHess(tmpmod$mle,model$fn,model$gr))
    tmpmod$LnDet    <-  sum(determinant(tmpmod$Hessian, logarithm=TRUE)$modulus[[1]])
  }
    
  lp              <-  model$env$last.par
  
  if (!se.fit) {
    pred          <- unlist(model$report(lp))
    tmpmod$pred   <- data.frame(def= names(pred),pred=pred,pred.se=NA)
  } else {
    tmpmod$sdr<-NULL
    if(hessianIN)
      tmpmod$sdr    <- sdreport(model,tmpmod$mle,hessian.fixed=tmpmod$Hessian,getReportCovariance=F)
    #tmpmod$sdrsum <- TMB:::summary.sdreport(tmpmod$sdr, "report")
    tmpmod$sdrsum <- summary(sdreport(model,getJointPrecision=TRUE))
    pred          <- tmpmod$sdrsum[,"Estimate"]
    se            <- tmpmod$sdrsum[,"Std. Error"]
    tmpmod$pred   <- data.frame(def= names(pred),pred=pred,pred.se=se)
  }
  
  if(simulate){
    #tmpmod$simdat <- array(sim_nitr)
    sim <- replicate(sim_nitr, {
      simdata <- model$simulate(par = model$par, complete=TRUE)
      
      obj2    <- MakeADFun(data       = simdata,
                           parameters = dlistIN$parameters,
                           DLL        = version,
                           checkParameterOrder=TRUE,
                           hessian    = TRUE,
                           map        = dlistIN$map,
                           profile    = NULL, # TODO: Optionally "beta"
                           silent     = TRUE)
      nlminb(obj2$par, obj2$fn, obj2$gr)$par
    })
    tmpmod$sim_df <- data.frame(estimate=as.vector(sim), parameter=names(model$par)[row(sim)])
    tmpmod$sim    <- sim
    
  }
  setwd(wd0 )
  
  
  tmpmod$fitted <- 
    data.frame(age = (tmpmod$report$age),
             W_hat = exp(tmpmod$report$logWhat),
             W_obs = exp(tmpmod$report$logWobs))
  
  return(tmpmod)
}