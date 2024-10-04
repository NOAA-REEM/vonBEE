#' makeMap will create a map list for reading into the runmod() function
#'
#' This function prepares the estimated data for futR
#' This function is a subset of makeDat()
#' For more information contact author Kirstin Holsman (kirstin.holsman@noaa.gov)
#'
#' @param parameters  list of starting values and parameters for the TMB futR() model
#' @param estpar      vector of T/F of which parameters will be estimated
#' @export            maplist a list with map values
#'
#' @examples
#'
#' datlist   <- readMake_futR_data("data/in/futR_Inputs.xlsx" )
#'
#' estparams <- datlist$estparams
#' estparams
#'
#' parameters <- atlist$parameters
#' parameters
#' maplist <- makeMap(param=parameters,estpar=estparams)
#' @export
makeMap<-function(param=parameters,estpar=estparams){
  mapseq  <- 1:length(unlist(param))
  ixx     <- 0
  maplist <- param
  for(inm in names(param)){
    if(inm %in% names(estpar)){
      if(estpar[[inm]]){

        ixx             <-  1:length(param[[inm]])+rev(ixx)[1]
        maplist[[inm]]  <- factor(ixx)

      }else{
        maplist[[inm]]  <- factor(rep(NA,length(maplist[[inm]])))
      }
    }else{
      maplist[[inm]]    <- factor(rep(NA,length(maplist[[inm]])))
    }
  }
  return(maplist)

}
