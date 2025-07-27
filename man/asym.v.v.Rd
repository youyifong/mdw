\name{asym.v.v}
\alias{asym.v.v}
\title{Asymptotic variance for maximum variance weights}
\usage{
asym.v.v(X, w)
}
\arguments{
\item{X}{n by p maxtrix containing observations of p biomarkers of n subjects.}

\item{w}{maximum variance weights for dateset X}
}
\description{
asym.v.v produces estimated asymptotic covariance matrix of the first p-1 maximum variance weights (because the p weights sum to 1).
}
\examples{
library(MASS)
# a three biomarkers dataset generated from independent normal(0,1)
X = mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
w <- var.weight(X)
asym.v.v(X,w)
}
\keyword{weighting}
