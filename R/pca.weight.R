#' Weights based on PCA
#'
#' pca.weight produce the coefficients of the first principal compoment
#' @param emp.cor empirical correlation matrix of the dataset
#' @keywords PCA
#' @export
#' @examples
#' library(MASS)
#' # a three biomarkers dataset generated from independent normal(0,1)
#' X = mvrnorm(n = 100, mu=rep(0,3), Sigma=diag(3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#' emp.cor <- cor(X)
#' pca.weight(emp.cor)


pca.weight=function(emp.cor){
  p=nrow(emp.cor)
  tmp=eigen(emp.cor)
  tmp$vectors[,1]
}