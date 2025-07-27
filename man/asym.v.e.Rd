\name{asym.v.e}
\alias{asym.v.e}
\title{Asymptotic variance for maximum entropy weights}
\usage{
asym.v.e(X, w, h)
}
\arguments{
\item{X}{n by p maxtrix containing observations of p biomarkers of n subjects.}

\item{w}{maximum entropy weights for dateset X with bandwidth h used}

\item{h}{bandwidth for kernel density estimation.}
}
\description{
asym.v.e produces estimated asymptotic covariance matrix of the first p-1 maximum entropy weights (because the p weights sum to 1).
}
\examples{
library(MASS)
# a three biomarkers dataset generated from independent normal(0,1)
X = mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
h = 1
w <- entropy.weight(X,h)
asym.v.e(X,w,h)
}
\keyword{weighting}
