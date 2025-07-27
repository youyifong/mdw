\name{entropy.weight}
\alias{entropy.weight}
\title{Maximum entropy weights}
\usage{
entropy.weight(X, h)
}
\arguments{
\item{X}{n by p maxtrix containing observations of p biomarkers of n subjects.}

\item{h}{bandwidth for kernel density estimation. if data is categorical, set to 'na'.}
}
\description{
entropy.weight produces a set of weights that maximizes the total weighted entropy of the distribution of different biomarkers within each subject, values of biomarkers can be either continuous or categorical.
}
\examples{
library(MASS)
# a three biomarkers dataset generated from independent normal(0,1)
set.seed(1)
X = mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
entropy.weight(X, h=1)
###
# a three categorical biomarkers dataset
set.seed(1)
tmp=mvrnorm(n=10,mu=c(0,0,0),Sigma = diag(3))
dat=t(apply(tmp, 1, function(x) cut(x,c(-Inf,-0.5,0.5,Inf),labels=1:3)))
entropy.weight(dat,h='na')
}
\keyword{weighting}
