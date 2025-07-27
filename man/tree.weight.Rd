\name{tree.weight}
\alias{tree.weight}
\title{ Weights based on GSC Tree Method }
\description{ 
  \code{tree.weight} Produce a set of weights for different end points based on a correlation matrix using the GSC tree method
}

\usage{
tree.weight (cor.mat, method="GSC", clustering.method="average", plot=TRUE, 
    orientation=c("vertical","horizontal"), ...) 

}

\arguments{

  \item{cor.mat}{ a matrix, correlation matrix }
  \item{method}{ a string. GSC, implementation of Gerstein et al., is the only implemented currently }
  \item{clustering.method}{ a string, how the bottom-up hierarchical clustering tree is built, is passed to hclust as the method parameter }
  \item{plot}{ a Boolean, whether to plot the tree }
  \item{orientation}{vertical or horizontal}
  \item{\dots}{additional args}

}


\value{
A vector of weights that sum to 1.
}

\author{ Youyi Fong \email{yfong@fhcrc.org} }

\references{
Gerstein, M., Sonnhammer, E., and Chothia, C. (1994), Volume changes in protein evolution. J Mol Biol, 236, 1067-78.
}

\examples{

cor.mat=diag(rep(1,3))
cor.mat[1,2]<-cor.mat[2,1]<-0.9
cor.mat[1,3]<-cor.mat[3,1]<-0.1
cor.mat[2,3]<-cor.mat[3,2]<-0.1
tree.weight(cor.mat)    

}
\keyword{ distribution }
