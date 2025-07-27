#' Maximum entropy weights
#'
#' entropy.weight produces a set of weights that maximizes the total weighted entropy of the distribution of different biomarkers within each subject, values of biomarkers can be either continuous or categorical. 
#' @param X n by p maxtrix containing observations of p biomarkers of n subjects.
#' @param h bandwidth for kernel density estimation，if data is categorical, set to 'na'.
#' @keywords weighting
#' @export
#' @importFrom "stats" "dnorm" "integrate"
#' @examples
#' library(MASS)
#' # a three biomarkers dataset generated from independent normal(0,1)
#' set.seed(1)
#' X = mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#' entropy.weight(X, h=1)
#' ###
#' # a three categorical biomarkers dataset
#' set.seed(1)
#' tmp=mvrnorm(n=10,mu=c(0,0,0),Sigma = diag(3))
#' dat=t(apply(tmp, 1, function(x) cut(x,c(-Inf,-0.5,0.5,Inf),labels=1:3)))
#' entropy.weight(dat,h='na')

entropy.weight <-function(X, h){
  p=ncol(X)
  n=nrow(X)
  if(h=='na'){
    m=length(unique(as.vector(X)))
    initial<-rep(1/p, p-1)
    Ci <- c(rep(c(0,-1),p-1),-1)
    Ui <- matrix(0, ncol =  p-1, nrow = 2*(p-1))
    for(j in 1:(p-1)){
      Ui[(2*j-1):(2*j), j] <- c(1, -1)
    }
    Ui <- rbind(Ui, rep(-1, p-1))
    I=array(0,c(n,p,m))
    for(j in 1:m){
      for(i in 1:n){
        for(k in 1:p){
          I[i,k,j]=ifelse(test = X[i,k]==j,yes = 1,no=0)
        }
      }
    }
    cp <- function(w){
      w <- c(w[1:(p-1)],1-sum(w))
      pr=matrix(0,nrow = n,ncol=m)
      ent=matrix(0,nrow = n,ncol=m)
      for(j in 1:m){
        for(i in 1:n){
          pr[i,j]=sum(w*I[i,,j])
          ent[i,j]=ifelse(test = pr[i,j]==0,yes = 0,no=-pr[i,j]*log(pr[i,j]))
        }
      }
      return(sum(ent))
    }
    
    res <- constrOptim(initial, cp, grad = NULL, ui = Ui, ci = Ci, control=list(fnscale=-1))$par

  } else{
    initial<-rep(1/p, p-1)
    Ci <- c(rep(c(0,-1),p-1),-1)
    Ui <- matrix(0, ncol =  p-1, nrow = 2*(p-1))
    for(j in 1:(p-1)){
      Ui[(2*j-1):(2*j), j] <- c(1, -1)
    }
    Ui <- rbind(Ui, rep(-1, p-1))
    
    X <- t(X)
    cp <- function(w){
      w <- c(w[1:(p-1)],1-sum(w))
      temp <- c(0, ncol(X))
      
      for (i in 1:ncol(X)) {
        
        cd <- function(x) {
          c <- 0
          for(j in 1:nrow(X)){
            c <- c + w[j]*dnorm(x, mean = X[j,i], sd = h)
          }
          return(c*log(c))
        }
        if(p>=100){
          temp[i] <- integrate(cd, range(X[,i])[1], range(X[,i])[2])$value 
        }
        temp[i] <- integrate(cd, range(X[,i])[1] - 3*h, range(X[,i])[2] + 3*h)$value
        
      }
      return(-sum(temp))
    }
    res <- constrOptim(initial, cp, grad = NULL, ui = Ui, ci = Ci, control=list(fnscale=-1))$par
  }
  res=c(res,  1- sum(res))
  eval(eval(substitute(expression( res.len <<- length(res) ))))     # set res.len to be used outside the current environment
  res
}
