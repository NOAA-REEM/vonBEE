# ----------------------------------------
# setup.R
# subset of Yasumiishi et al. in prep
# EBS salmon code
# 
# updated 2022
# ----------------------------------------

  for(d in dir("R/sub_fun")) 
    source(paste0("R/sub_fun/",d))
  
logit <- function(p){
  log(p/(1-p))
}
inv.logit <- function(x){
  1/(1+exp(-x))
}


se <- function(x,na.rm=T){
  if(na.rm) x<- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

  

aic <- function(k, nll){
  (2*k-2*-nll)
}
aicc <- function(k, nll,n){
  aic <- (2*k-2*-nll)
  return(aic+(2*k^2 + 2*k)/(n-k-1))
}
 
  
  
  

  
  
  
  