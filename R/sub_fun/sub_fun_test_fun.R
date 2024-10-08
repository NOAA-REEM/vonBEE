#' test_fun
#'
#' maybe defunct?
#'
#' @param x
#' @param age =5
#' @param t0 =-.25
#' @param K = 12
#' @param H = .25
#' @param offset =0
#'
#' @export
#'
#' @returns list(d=d, Winf = Winf,logWhat = logWhat, What=What, parta = parta , partb=partb)
#'

test_fun <- function(x,age=5,t0=-.25, K = 12,H = .25,offset=0){
  d          =  ((1-offset)/(1+exp(-x)))
  Winf       = (H/K)^(1.0/(1.0 - d) );
  parta      = (1.0/(1.0 - d))
  partb      = log(1.0 - exp(-K * (1.0 - d) * (age - t0)))
  logWhat    = log(Winf) + (1.0/(1.0 - d))*log(1.0 - exp(-K * (1.0 - d) * (age - t0))) ;
  What       = exp(logWhat);
  return(list(d=d, Winf = Winf,logWhat = logWhat, What=What, parta = parta , partb=partb))
}
