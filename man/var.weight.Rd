\name{var.weight}
\alias{var.weight}
\title{Maximum variance weights}
\usage{
var.weight(X, method = c("optim", "mosek"))
}
\arguments{
\item{X}{n by p maxtrix containing observations of p biomarkers of n subjects.}

\item{method}{optim (default) using R constrOptim function from stats package for optimization, mosek using mosek function from Rmosek package for optimization}
}
\description{
var.weight produces a set of weights that maximizes the total weighted variance of the distribution of different biomarkers within each subject.
}
\examples{
library(MASS)
# a three biomarkers dataset generated from independent normal(0,1)
X = mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
# compute maximum variance weights using constrOptim for optimization
var.weight(X)

\dontrun{
# need mosek installed
# compute maximum variance weights using mosek for optimization
library(Rmosek)
var.weight(X,'mosek')
}
}
\keyword{weighting}
