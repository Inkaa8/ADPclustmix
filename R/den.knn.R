#' @title K-th Nearest Neighbor Density Calculation
#'
#' @description This is the function to calculate KNN density.
#'
#' @param dist a numeric distance matrix or data frame object.
#' @param k k-th nearest neighbors to be checked to estimate the density of each point.
#'
#' @details
#' If \code{k} is \code{NULL}, k=n^(4/5), where \code{n} is the number of rows in the \code{dist}.
#'
#' @return a vector of density value.
#'
knn.density <- function(dist, k=NULL){
  distmat <- as.matrix(dist)
  n <- nrow(distmat)
  if (k >= n) {stop("k is not a reasonable and valid number!")}
  v.d <- pi^(1/2) /gamma(1/2+1)
  r.k <- apply(distmat, 1, sort)[k+1,]
  den <- k / (n * v.d * r.k)
  return(den)
}













