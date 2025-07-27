#' Asymptotic variance for maximum entropy weights
#'
#' asym.v.e produces estimated asymptotic covariance matrix of the first p-1 maximum entropy weights (because the p weights sum to 1). 
#' @param X n by p maxtrix containing observations of p biomarkers of n subjects.
#' @param h bandwidth for kernel density estimation.
#' @param w maximum entropy weights for dateset X with bandwidth h used 
#' @keywords weighting
#' @export
#' @importFrom "stats" "dnorm" "integrate"
#' @examples
#' library(MASS)
#' # a three biomarkers dataset generated from independent normal(0,1)
#' X = mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#' h = 1
#' w <- entropy.weight(X,h)
#' asym.v.e(X,w,h)


asym.v.e <- function(X,w,h){
  M_theta.e <- function(X,w,h){
    out <- matrix(0, nrow=nrow(X), ncol=ncol(X)-1)
    for(i in 1:nrow(X)){
      for(k in 1:(ncol(X)-1)){
        int.m <- function(x){
          w.sum <- 0
          for(j in 1:ncol(X)){
            w.sum <- w.sum + w[j]*dnorm(x, mean = X[i,j])
          }
          (dnorm(x, mean = X[i,k]) - dnorm(x, mean = X[i,ncol(X)]))*(1 + log(w.sum))
        }
        out[i,k] <- -integrate(int.m, range(X[i,])[1] - 3*h, range(X[i,])[2] + 3*h)$value
      }
    }
    t(out)%*%out/nrow(X)
  }
  
  V_theta.e <- function(X,w,h){
    temp <- array(0, c(ncol(X)-1, ncol(X)-1, nrow(X)))
    for(i in 1:nrow(X)){
      for(j in 1:(ncol(X)-1)){
        for(l in 1:(ncol(X)-1)){
          int.v <- function(x){
            w.sum <- 0
            for(k in 1:ncol(X)){
              w.sum <- w.sum + w[k]*dnorm(x, mean = X[i,k])
            }
            (dnorm(x, mean = X[i,j]) - dnorm(x, mean = X[i,ncol(X)]))*(dnorm(x, mean = X[i,l]) - dnorm(x, mean = X[i,ncol(X)]))/w.sum  
          }
          temp[j,l,i] <- -integrate(int.v, range(X[i,])[1] - 3*h, range(X[i,])[2] + 3*h)$value
        }
      }
    }
    apply(temp, c(1,2), sum)/nrow(X)
  }
  
    solve(V_theta.e(X,w,h))%*%M_theta.e(X,w,h)%*%solve(V_theta.e(X,w,h))
}

