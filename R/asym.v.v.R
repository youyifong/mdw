#' Asymptotic variance for maximum variance weights
#'
#' asym.v.v produces estimated asymptotic covariance matrix of the first p-1 maximum variance weights (because the p weights sum to 1). 
#' @param X n by p maxtrix containing observations of p biomarkers of n subjects.
#' @param w maximum variance weights for dateset X
#' @keywords weighting
#' @export
#' @examples
#' library(MASS)
#' # a three biomarkers dataset generated from independent normal(0,1)
#' X = mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#' w <- var.weight(X)
#' asym.v.v(X,w)


asym.v.v <- function(X,w){
  M_theta.v <- function(X,w){
    out <- matrix(0, nrow=nrow(X), ncol=ncol(X)-1)
    for(i in 1:nrow(X)){
      for(k in 1:(ncol(X)-1)){
        out[i,k] <- X[i,k]^2-X[i,ncol(X)]^2 - 2*(X[i,k]-X[i,ncol(X)])*(t(w)%*%X[i,])
      }
    }
    t(out)%*%out/nrow(X)
  }
  
  V_theta.v <- function(X){
    temp <- array(0, c(ncol(X)-1, ncol(X)-1, nrow(X)))
    for(i in 1:nrow(X)){
      for(j in 1:(ncol(X)-1)){
        for(l in 1:(ncol(X)-1)){
          temp[j,l,i] <- -2*(X[i,j]-X[i,ncol(X)])*(X[i,l]-X[i,ncol(X)]) 
        }
      }
    }
    apply(temp, c(1,2), sum)/nrow(X)
  }
  
  solve(V_theta.v(X))%*%M_theta.v(X,w)%*%solve(V_theta.v(X))
}

