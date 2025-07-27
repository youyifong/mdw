\name{pca.weight}
\alias{pca.weight}
\title{Weights based on PCA}
\usage{
pca.weight(emp.cor)
}
\arguments{
\item{emp.cor}{empirical correlation matrix of the dataset}
}
\description{
pca.weight produce the coefficients of the first principal compoment
}
\examples{
library(MASS)
# a three biomarkers dataset generated from independent normal(0,1)
X = mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
emp.cor <- cor(X)
pca.weight(emp.cor)
}
\keyword{PCA}
