#' @title Fast Clustering Using Adaptive Density Peak Detection For Mixed Type Data
#'
#' @description Clustering of mixed-type data by finding cluster centers from estimated density peaks. ADPclust is a non-iterative procedure that incorporates multivariate Gaussian density estimation. The number of clusters as well as bandwidths can either be selected by the user or selected automatically through an internal clustering criterion.
#'
#' @details Given n data points x's in p dimensions, adpclustmix() calculates f(x) and delta(x) for each data point x, where f(x) is the local density for mixed-type data at x, and delta(x) is the shortest distance for mixed-type data between x and y for all y such that f(x) <= f(y). Data points with large f and large delta values are labeled class centroids. In other words, they appear as isolated points in the upper right corner of the f vs. delta plot (the decision plot). After cluster centroids are determined, other data points are clustered according to their distances to the closes centroids.
#'
#' This package is built as an extension of \code{"ADPclust"} to handle mixed-type data. See \code{"adpclust"} for more details on parameter settings.
#'
#' For mixed type data, there are three methods available to calculate distance: \code{gower}, \code{rf}, \code{kp}, two methods to perform dimension reduction: \code{FAMD}, \code{mixedPCA} and three methods to estimate local density: \code{knn}, \code{kde}, code{NB}.
#'
#' \code{gower}
#'
#' Gower distance is the most common distance for a mixed variable data set. See Gower (1971).
#'
#' \code{rf}
#'
#' A distance measure based on the proximity induced by the given random forest model. The proximity of two rows is the number of trees in which the two rows end up in the same leaf divided by the total number of trees.
#'
#' \code{kp}
#'
#' K-prototype distance. It measures distance between numerical variables using Euclidean distance (like K-means) and also measures the distance between categorical variables using the number of matching categories. See Huang (1998).
#'
#' \code{FAMD}
#'
#'  Factor analysis of mixed data. FAMD algorithm can be seen as a mixed between principal component analysis (PCA) and multiple correspondence analysis (MCA). See Pagès Jérôme (2014).
#'
#' \code{mixedPCA}
#'
#' It performs principal component analysis of a set of individuals (observations) described by a mixture of qualitative and quantitative variables. It includes ordinary principal component analysis (PCA) and multiple correspondence analysis (MCA) as special cases. See Chavent M. (2014).
#'
#' \code{knn}
#'
#' Nonparametric K-nearest neighbors (k-NN) density estimation.
#'
#' \code{kde}
#'
#'  Nonparametric kernel density estimation, which computes kernel unconditional density estimates on evaluation data, given a set of training data and a bandwidth specification (a bandwidth object or a bandwidth vector, bandwidth type, and kernel type). See Li (2003).
#'
#' \code{NB}
#'
#' Naïve-bayes density, which estimates the density of mixed data assuming independence of variables.
#'
#'
#' @param data a data frame where rows are observations and columns are variables. Mixed type variables are allowed. Categorical variables should be factorized.
#' @param dist.method method to calculate distance matrix. Valid options are \code{"gower"}, \code{"rf"}, \code{"kp"}, \code{"FAMD"} and \code{"mixedPCA"}.
#' @param dens.method method to estimate density. Valid options are \code{"knn"}, \code{"kde"} and \code{"NB"}.
#' @param ... arguments passed to \code{adpclust()}, e.g. \code{nclust=2:6}, or arguments specific to selected distance/density measure, see Details for more information.
#'
#'
#' @return an 'adpclust' object that contains the list of the following items.
#' \itemize{
#' \item{clusters:}{ Cluster assignments. A vector of the same length as the number of observations.}
#' \item{centers:}{ Indices of the clustering centers.}
#' \item{silhouette:}{ Silhouette score from the final clustering result.}
#' \item{nclust:}{ Number of clusters.}
#' \item{h:}{ Final bandwidth.}
#' \item{f:}{ Final density vector f(x).}
#' \item{delta:}{ Final delta vector delta(x).}
#' \item{selection.type:}{ 'user' or 'auto'.}
#' }
#'
#' @import ADPclust
#' @import kmed
#' @import randomForest
#' @import np
#' @import FactoMineR
#' @import PCAmixdata
#' @importFrom stats as.formula as.dist
#' @export
#'
#' @examples
#' # Example 1
#' # Load mixed type data with 2 clusters
#' data(dat_demo1)
#'
#' # Select "gower" distance and "NB" density for clustering
#' res <- adpclustmix(dat_demo1, dist.method = "gower", dens.method = "NB",con.den = "locfit",
#' idnum = 1:2, idcat = 3:4)
#' summary(res)
#' plot(res)
#'
#' # Select "rf" distance and "kde" density for clustering
#' res <- adpclustmix(dat_demo1, dist.method = "rf", dens.method = "kde")
#' summary(res)
#' plot(res)
#'
#'
#' # Example 2
#' # Load mixed type data with 4 clusters and 2 noise variables
#' data(dat_demo2)
#' # Select "kp" distance and "knn" density for clustering
#' res <- adpclustmix(dat_demo2, dist.method = "kp", dens.method = "knn",
#' idnum = c(1:5,11), idcat = c(6:10,12))
#' summary(res)
#' plot(res)
#'
#' # Specify "FAMD" as distance measure to do dimension reduction followed by normal adpclust
#' res <- adpclustmix(dat_demo2, dist.method = "FAMD")
#' summary(res)
#' plot(res)
#'
adpclustmix <- function(data, dist.method, dens.method, ...){

  parms <- list(...)

  # if data is all continuous
  if(all(sapply(data, class) == "numeric")){
    adp.res <- adpclust(data, ...)
    return(adp.res)

  } else{ # if data is mixed type

    # -----------------------------distance--------------------------------#
    # distance measure function
    if(dist.method == "gower"){
      dist.gower <- distmix(data, method = "gower",
                            idnum = parms$idnum, idcat = parms$idcat)
      dist <- as.dist(dist.gower)


    } else if(dist.method == "rf"){
      if(!("mtry" %in% names(parms))){
        parms$mtry=3
      }
      if(!("no.tree" %in% names(parms))){
        parms$no.tree=1000
      }

      RF.prox <- randomForest(data, mtry = parms$mtry, ntree = parms$no.tree, proximity = TRUE)$proximity
      dist <- as.dist(sqrt(1 - RF.prox))


    } else if(dist.method == "kp"){
      dist.huang <- distmix(data, method = "huang",
                            idnum = parms$idnum, idcat = parms$idcat)
      dist <- as.dist(dist.huang)


    } else if(dist.method == "FAMD"){

      if(!("ncp" %in% names(parms))) {
        parms$ncp <- 5
      }
      parms_famd <- c(list(base=data), graph=FALSE, ncp=parms$ncp)
      res.famd <- do.call(FAMD, parms_famd)
      pc <- res.famd$ind$coord


    } else if(dist.method == "mixedPCA"){
      if(!("ncp" %in% names(parms))){
        parms$ncp <- 5
      }

      parms_pcamixed <- c(list(X.quanti=data[, parms$idnum],
                               X.quali=data[, parms$idcat]),
                          ndim=parms$ncp,
                          graph=FALSE,
                          rename.level=TRUE)
      res.pcamix <- do.call(PCAmix, parms_pcamixed)
      pc <- res.pcamix$ind$coord


    }else {
      stop("'dist.method' must be one of 'gower','rf','kp','FAMD','mixedPCA'.
             Other methods are currently unavailable.")
    }


    # if chose dimension reduction method (FAMD or mixedPCA), proceed to normal ADPclust procedure
    # on continuous data
    if (dist.method %in% c("FAMD", "mixedPCA")){
      parms_adp <- c(list(x=pc), parms[!(names(parms) %in% c("idnum", "idcat", "ncp"))])
      adp.res <- do.call(adpclust, parms_adp)

      return(adp.res)

    }

    # -----------------------------density--------------------------------#
    # density measure function
    if(dens.method == "knn"){
      n <- nrow(as.matrix(dist))
      if(is.null(parms$k)){
        parms$k <- n^(4/5)
      }

      parms_knn <- c(list(dist=dist), k=parms$k)
      den.knn <- do.call(knn.density, parms_knn)
      den <- den.knn


    } else if(dens.method == "kde"){
      den.form <- as.formula(paste("~", paste(names(data), collapse = "+")))
      den.exact <- npudens(den.form, data = data)$dens
      den <- den.exact


    } else if(dens.method == "NB"){
      den <- NB.cal(data, idnum = parms$idnum, idcat = parms$idcat, con.den = parms$con.den)


    } else{
      stop("'den.method' must be one of 'knn', 'kde' and 'NB'.
                 Other methods are currently unavailable.")

    }


    # ADPclust procedure using distance measure and density measure
    parms_adp <- c(list(distm=dist, f=den),
                   parms[!names(parms) %in% c("idnum", "idcat", "ncp",
                                              "mtry", "no.tree", "no.rep",
                                              "k","con.den","weight")])
    adp.res <- do.call(adpclust, parms_adp)
    return(adp.res)

  }

}




