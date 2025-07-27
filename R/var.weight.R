#' Maximum variance weights
#'
#' var.weight produces a set of weights that maximizes the total weighted variance of the distribution of different biomarkers within each subject.
#' @param X n by p maxtrix containing observations of p biomarkers of n subjects.
#' @param method optim (default) using R constrOptim function from stats package for optimization, mosek using mosek function from Rmosek package for optimization
#' @keywords weighting
#' @export 
#' @import "stats" "Rmosek" "MASS"
#' @importFrom "Matrix" "Matrix"
#' @examples
#' library(MASS)
#' library(Rmosek)
#' # a three biomarkers dataset generated from independent normal(0,1)
#' X = mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#' # compute maximum variance weights using constrOptim for optimization
#' var.weight(X)
#' # compute maximum variance weights using mosek for optimization
#' var.weight(X,'mosek')

var.weight <-function(X,method=c('optim','mosek')){
  p=ncol(X)
  n=nrow(X)
  method=match.arg(method)
  if(method=='optim'){
    initial<-rep(1/p, p-1)
    Ci <- c(rep(c(0,-1),p-1),-1)
    Ui <- matrix(0, ncol = p-1, nrow =  2*(p-1))
    for(j in 1:(p-1)){
      Ui[(2*j-1):(2*j), j] <- c(1, -1)
    }
    Ui <- rbind(Ui, rep(-1, p-1))
    
    cp <- function(w){
      w<-c(w[1:(p-1)], 1-sum(w))
      t1 <- t(w)%*%t(X)^2
      t2 <- t(w)%*%t(X)
      return(sum(t1 - (t2)^2))
    }
    
    res <- constrOptim(theta = initial, f = cp, grad = NULL, ui = Ui, ci = Ci, control=list(fnscale=-1))$par
    res=c(res,  1- sum(res))
  } else if(method=='mosek'){
    Q0=0
    for(i in 1:n){
      Q0=Q0+X[i,]%*%t(X[i,])
    }
    Q0=-2*Q0
    var.w <- list(sense = "max")
    var.w$c <- colSums(X^2)
    var.w$A <- Matrix::Matrix( rep(1,p), nrow=1, byrow=TRUE)
    var.w$bc <- rbind(blc=1,buc=1)
    var.w$bx <- rbind(blx = rep(0, p),
                      bux = rep(1, p))
    ind.row=NULL
    for(i in 1:p){
      ind.row=c(ind.row,rep(i,i))
    }
    ind.col=NULL
    for(i in 1:p){
      ind.col=c(ind.col,1:i)
    }
    values=0
    for(v in 1:length(ind.row)){
      values[v]=Q0[ind.row[v],ind.col[v]]
    }
    var.w$qobj <- list(i = ind.row,
                       j = ind.col,
                       v = values)
    fit.mosek <- Rmosek::mosek(var.w,opts=list(verbose=0))
    if ( inherits (fit.mosek$sol$itr$xx, "NULL") ){
      res=rep(NA,p)
    } else{
      res=fit.mosek$sol$itr$xx
    }
  }

  
  eval(eval(substitute(expression( res.len <<- length(res) ))))     # set res.len to be used outside the current environment

  res
}
