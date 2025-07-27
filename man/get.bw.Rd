\name{get.bw}
\alias{get.bw}
\title{Bandwidth Selection}
\usage{
get.bw(x, bw = c("nrd", "ucv", "bcv", "SJ"), nb)
}
\arguments{
\item{x}{n by p maxtrix containing observations of p biomarkers of n subjects.}

\item{bw}{bandwidth selectors of nrd, ucv, bcv, and SJ corresponding to R functions bw.nrd, bw.ucv, bw.bcv, and bw.SJ.}

\item{nb}{number of bins to use, 'na' if bw='nrd'}
}
\description{
get.bw applies a specified bandwidth selection method to the dataset subject-wisely and return the median of the n selected bandwidths as the choice of bandwidth for entropy.weight.
}
\examples{
library(MASS)
# a ten biomarkers dataset generated from independent normal(0,1)
x = mvrnorm(n = 100, mu=rep(0,10), Sigma=diag(10), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
get.bw(x,bw='ucv',nb=100)
get.bw(x,bw='nrd',nb='na')
}
\keyword{bandwidth}
\keyword{selection}
