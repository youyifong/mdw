#' Bandwidth Selection
#'
#' get.bw applies a specified bandwidth selection method to the dataset subject-wisely and return the median of the n selected bandwidths as the choice of bandwidth for entropy.weight.
#' @param x n by p maxtrix containing observations of p biomarkers of n subjects.
#' @param bw bandwidth selectors of nrd, ucv, bcv, and SJ corresponding to R functions bw.nrd, bw.ucv, bw.bcv, and bw.SJ.
#' @keywords bandwidth selection
#' @param nb number of bins to use, 'na' if bw='nrd' 
#' @export
#' @examples
#' library(MASS)
#' # a ten biomarkers dataset generated from independent normal(0,1)
#' x = MASS::mvrnorm(n = 100, mu=rep(0,10), Sigma=diag(10), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#' get.bw(x,bw='ucv',nb=100)
#' get.bw(x,bw='nrd',nb='na')

get.bw=function(x,bw=c('nrd','ucv','bcv','SJ'),nb){
  h=0
  bw=match.arg(bw)
  for(i in 1:nrow(x)){
    if(nb=='na'){
      if(bw=='nrd'){
        h[i]=bw.nrd(x[i,]) 
      } 
    } else{
      nb=as.numeric(nb)
      if(bw=='ucv'){
        h[i]=bw.ucv(x[i,],nb=nb,tol=1e-10)
      } else if(bw=='bcv'){
        h[i]=bw.bcv(x[i,],nb=nb,tol=1e-10)
      } else if(bw=='SJ.ste'){
        h[i]=bw.SJ(x[i,],nb=nb,method='ste',tol=1e-10)
      } else if(bw=='SJ.dpi'){
        h[i]=bw.SJ(x[i,],nb=nb,method='dpi',tol=1e-10)
      }
    }
  }
  median(h)
}
