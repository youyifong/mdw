# copied from lcmix on Rforge
rmvgamma <- function(n, shape=1, rate=1, corr=diag(length(shape)))
{
    ## extract parameters, do sanity checks, deal with univariate case
    
    if(!is.matrix(corr) || !isSymmetric(corr))
        stop("'corr' must be a symmetric matrix")
    D = ncol(corr)
    
    Ds = length(shape)
    if(Ds > D)
        warning("'shape' longer than width of 'corr', truncating to fit")
    if(Ds != D)
        shape = rep(shape, length.out=D)
    
    Dr = length(rate)
    if(Dr > D)
        warning("'rate' longer than width of 'corr', truncating to fit")
    if(Dr != D)
        rate = rep(rate, length.out=D)
    
    if(D == 1) rgamma(n, shape, rate)
    
    ## generate standard multivariate normal matrix, convert to CDF
    
    Z = mvtnorm::rmvnorm(n, cov=corr)
    cdf = pnorm(Z)
    
    ## convert to gamma, return
    
    sapply(1:D, function(d) qgamma(cdf[,d], shape[d], rate[d]))
}

sim.dat <- function(dist,nlv,n,p,rho,scenario='na'){
  if(nlv>0){
    if(dist=='normal'){
      if(rho==0){
        X = MASS::mvrnorm(n = n, mu=rep(0,p), Sigma=diag(p), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
      } else{
        eps=sqrt((1-rho)/rho)
        L = MASS::mvrnorm(n = n, mu=rep(0,nlv), Sigma=diag(nlv), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
        X=matrix(0,ncol=p,nrow=n)
        for(i in 1:(.2*p)){
          X[,i]=L[,1] + rnorm(n, 0, eps)
        }
        for(i in (.2*p+1):p){
          X[,i]=L[,2] + rnorm(n, 0, eps)
        }
      }
    } else if(dist=='lognormal'){
      eps=sqrt(exp(1)*(exp(1)-1))*sqrt((1-rho)/rho)
      L=matrix(0,ncol=nlv,nrow=n)
      for(i in 1:nlv){
        L[,i]=rlnorm(n=n,meanlog = 0,sdlog = 1)
      }
      X=matrix(0,ncol=p,nrow=n)
      for(i in 1:(.2*p)){
        X[,i]=L[,1] + rnorm(n, 0, 1/3*sqrt(exp(1)*(exp(1)-1)))
      }
      for(i in (.2*p+1):p){
        X[,i]=L[,2] + rnorm(n, 0, eps)
      }
    }
  } else if(nlv==0){
    if(dist=='lognormal'){
      if(p==5){
        Mu=exp(rep(exp(1/2),p))
        Sigma=rep((exp(1)-1)*(exp(1)),p)
        if(scenario=='a'){
          R=diag(p)
          R[1,2]=R[2,1]=R[4,5]=R[5,4]=0.93
          R[1,3]=R[2,3]=R[3,1]=R[3,2]=rho
        } else if(scenario=='b')
          R=diag(p)
        R[1,2]=R[2,1]=R[4,5]=R[5,4]=rho
      }
      # X=MethylCapSig::mvlognormal(n=n,Mu=Mu,Sigma=Sigma,R=R)
      stop("stopped since MethylCapSig is not on CRAN anymore")
    } else if(dist=='t'){
      if(p==5){
        if(scenario=='a'){
          Sigma=diag(p)
          Sigma[1,2]=Sigma[2,1]=Sigma[4,5]=Sigma[5,4]=0.9
          Sigma[1,3]=Sigma[2,3]=Sigma[3,1]=Sigma[3,2]=rho
        } else if(scenario=='b'){
          Sigma=diag(p)
          Sigma[1,2]=Sigma[2,1]=Sigma[4,5]=Sigma[5,4]=rho
        }
        X=mvtnorm::rmvt(n=n,sigma=Sigma,df=5)
      }
    } else if(dist=='gamma'){
      if(p==5){
        if(scenario=='a'){
          Sigma=diag(p)
          Sigma[1,2]=Sigma[2,1]=Sigma[4,5]=Sigma[5,4]=0.9
          Sigma[1,3]=Sigma[2,3]=Sigma[3,1]=Sigma[3,2]=rho
        } else if(scenario=='b'){
          Sigma=diag(p)
          Sigma[1,2]=Sigma[2,1]=Sigma[4,5]=Sigma[5,4]=rho
        }
        X=rmvgamma(n=n,shape=rep(2,p),rate=rep(.5,p),corr=Sigma)
      }
    } else{
      sigma=rbind(cbind(matrix(data=rho,ncol=.2*p,nrow=.2*p),matrix(data=0,ncol=.8*p,nrow=.2*p)),cbind(matrix(data=0,ncol=.2*p,nrow=.8*p),matrix(data=rho,ncol=.8*p,nrow=.8*p)))
      diag(sigma)=1
      X=MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=sigma)
    }
  }
  scale(X)
}

sim.cat=function(n, seed, k=100) {
  set.seed(seed)
  m1=c(1/3,1/3,1/3); m2=c(1/3,1/3,1/3)
  # m1=c(.6, .2, .2)
  # m2=c(.2, .2, .6)
  #    m1=c(.8, .2)
  #    m2=c(.2, .8)
  t(sapply (1:n, function(i) {
    if (k!=Inf){
      p1<-p2<-p3<-gtools::rdirichlet(1, m1*k)
      p4<-p5<-gtools::rdirichlet(1, m2*k)
    } else {
      p1<-p2<-p3<-m1
      p4<-p5<-m2
    }
    sapply(1:5, function(p) rbern(1, get("p"%.%p), generalized = T) )
  } ) )
}
