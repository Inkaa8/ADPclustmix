###### data generating function
########################################################################################
########################################################################################
# NORMAL
library(MixSim)

# Mixsim()
# BarOmega: value of desired average overlap.
# MaxOmega: value of desired maximum overlap.
# K:number of components.
# p: number of dimensions.
# sph:covariance matrix structure (FALSE - non-spherical, TRUE - spherical).
# hom: heterogeneous or homogeneous clusters (FALSE - heterogeneous, TRUE - homogeneous).
# ecc: maximum eccentricity.
# PiLow: value of the smallest mixing proportion (if 'PiLow' is not reachable with respect to K, equal proportions are taken; PiLow = 1.0 implies equal proportions by default).
# int: mean vectors are simulated uniformly on a hypercube with sides specified by int = (lower.bound, upper.bound). default c(0,1)
# resN: maximum number of mixture resimulations.
# eps: error bound for overlap computation.
# lim: maximum number of integration terms (Davies, 1980).

# simdataset()
# n: sample size.
# Pi: vector of mixing proportions (length K).
# Mu: matrix consisting of components' mean vectors (K * p).
# S: set of components' covariance matrices (p * p * K).
# n.noise: number of noise variables.
# n.out: number of outlying observations.
# alpha: level for simulating outliers.
# max.out: maximum number of trials to simulate outliers.
# int: interval for noise and outlier generation.
# lambda: inverse Box-Cox transformation coefficients.

genComplexMixDat <- function(n, K, p.con, p.cat, BarOmega=NULL, MaxOmega=NULL, sph=FALSE, hom=FALSE, Pilow=1.0,
                             p.noise.con, p.noise.cat, n.out){
  Q <- MixSim(BarOmega = BarOmega, MaxOmega = MaxOmega, K = K, p = p.con + p.cat, sph=FALSE, hom=FALSE, PiLow = Pilow)
  A <- simdataset(n = n, Pi = Q$Pi, Mu = Q$Mu, S = Q$S,
                  n.noise = p.noise.con + p.noise.cat, n.out = n.out)

  A$X <- as.data.frame(A$X)
  # scale continuous variable
  A$X[1:p.con] <- as.vector(scale(A$X[1:p.con]))

  # discretize to make categorize variable
  A$X[(p.con+1):(p.con+p.cat)] <- do.call(data.frame, Map(function(x) cut(x, breaks = quantile(x, probs = c(0, 0.2, 0.4, 0.7, 1)), labels = FALSE, include.lowest = TRUE),
                                                          A$X[(p.con+1):(p.con+p.cat)]))
  A$X[(p.con+1):(p.con+p.cat)] <- as.data.frame(lapply(A$X[(p.con+1):(p.con+p.cat)], factor))


  # if noise variables exist
  if (p.noise.con != 0 | p.noise.cat!=0){
    # scale continuous continuous noise variables
    A$X[(p.con+p.cat+1):(p.con+p.cat+p.noise.con)] <- as.vector(scale(A$X[(p.con+p.cat+1):(p.con+p.cat+p.noise.con)]))
    # randomly discretize to make categorical noise variable
    A$X[(p.con+p.cat+p.noise.con+1):(p.con+p.cat+p.noise.con+p.noise.cat)] <- do.call(data.frame, Map(function(x) cut(x, breaks = quantile(x, probs = c(0, 0.2, 0.4, 0.7, 1)), labels = FALSE, include.lowest = TRUE),
                                                                                                      A$X[(p.con+p.cat+p.noise.con+1):(p.con+p.cat+p.noise.con+p.noise.cat)]))
    A$X[(p.con+p.cat+p.noise.con+1):(p.con+p.cat+p.noise.con+p.noise.cat)] <- as.data.frame(lapply(A$X[(p.con+p.cat+p.noise.con+1):(p.con+p.cat+p.noise.con+p.noise.cat)], factor))

  }

  names(A$X) <- paste0("V", 1:(p.con+p.cat+p.noise.con+p.noise.cat))
  output <- list("SimDat" = A$X, "TrueID" = A$id, "Q" = Q)

}





########################################################################################
########################################################################################
# half LOGNORMAL (on continuous variable)
genComplexMixDat_ln <- function(n, K, p.con, p.cat, BarOmega=NULL, MaxOmega=NULL, sph=FALSE, hom=FALSE, Pilow=1.0,
                                p.noise.con, p.noise.cat, n.out){
  Q <- MixSim(BarOmega = BarOmega, MaxOmega = MaxOmega, K = K, p = p.con + p.cat, sph=FALSE, hom=FALSE, PiLow = Pilow)
  A <- simdataset(n = n, Pi = Q$Pi, Mu = Q$Mu, S = Q$S,
                  n.noise = p.noise.con + p.noise.cat, n.out = n.out)

  A$X <- as.data.frame(A$X)

  ### log-normal transformation
  A$X[1:p.con] <- exp(A$X[1:p.con])
  # scale continuous variable
  A$X[1:p.con] <- as.vector(scale(A$X[1:p.con]))

  # discretize to make categorize variable
  A$X[(p.con+1):(p.con+p.cat)] <- do.call(data.frame, Map(function(x) cut(x, breaks = quantile(x, probs = c(0, 0.2, 0.4, 0.7, 1)), labels = FALSE, include.lowest = TRUE),
                                                          A$X[(p.con+1):(p.con+p.cat)]))
  A$X[(p.con+1):(p.con+p.cat)] <- as.data.frame(lapply(A$X[(p.con+1):(p.con+p.cat)], factor))


  # if noise variables exist
  if (p.noise.con != 0 | p.noise.cat!=0){
    ### log-normal transformation
    A$X[(p.con+p.cat+1):(p.con+p.cat+p.noise.con)] <- exp(A$X[(p.con+p.cat+1):(p.con+p.cat+p.noise.con)])
    # scale continuous continuous noise variables
    A$X[(p.con+p.cat+1):(p.con+p.cat+p.noise.con)] <- as.vector(scale(A$X[(p.con+p.cat+1):(p.con+p.cat+p.noise.con)]))

    # discretize to make categorical noise variable
    A$X[(p.con+p.cat+p.noise.con+1):(p.con+p.cat+p.noise.con+p.noise.cat)] <- do.call(data.frame, Map(function(x) cut(x, breaks = quantile(x, probs = c(0, 0.2, 0.4, 0.7, 1)), labels = FALSE, include.lowest = TRUE),
                                                                                                      A$X[(p.con+p.cat+p.noise.con+1):(p.con+p.cat+p.noise.con+p.noise.cat)]))
    A$X[(p.con+p.cat+p.noise.con+1):(p.con+p.cat+p.noise.con+p.noise.cat)] <- as.data.frame(lapply(A$X[(p.con+p.cat+p.noise.con+1):(p.con+p.cat+p.noise.con+p.noise.cat)], factor))

  }

  names(A$X) <- paste0("V", 1:(p.con+p.cat+p.noise.con+p.noise.cat))
  output <- list("SimDat" = A$X, "TrueID" = A$id, "Q" = Q)

}



########################################################################################
########################################################################################
# Inverse BOX-COX transformation (by lambda parm)
genComplexMixDat_bc <- function(n, K, p.con, p.cat, BarOmega=NULL, MaxOmega=NULL, sph=FALSE, hom=FALSE, Pilow=1.0, lambda=NULL,
                                p.noise.con, p.noise.cat, n.out=0){
  Q <- MixSim(BarOmega = BarOmega, MaxOmega = MaxOmega, K = K, p = p.con + p.cat, sph=FALSE, hom=FALSE, PiLow = Pilow)
  # A <- simdataset(n = n, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, lambda = lambda,
  #                 n.noise = p.noise.con + p.noise.cat, n.out = n.out)

  suppressWarnings(A <- simdataset(n = n, Pi = Q$Pi, Mu = Q$Mu, S = Q$S, lambda = lambda,
                                   n.noise = p.noise.con + p.noise.cat, n.out = n.out))

  A$X <- as.data.frame(A$X)

  # deal with NaN produced during inverse box-cox transformation
  A$X <- replace(A$X, is.na(A$X), -1)

  # scale continuous variable
  A$X[1:p.con] <- as.vector(scale(A$X[1:p.con]))
  # discretize to make categorize variable
  A$X[(p.con+1):(p.con+p.cat)] <- do.call(data.frame, Map(function(x) cut(x, breaks = quantile(x, probs = c(0, 0.2, 0.4, 0.7, 1)), labels = FALSE, include.lowest = TRUE),
                                                          A$X[(p.con+1):(p.con+p.cat)]))
  A$X[(p.con+1):(p.con+p.cat)] <- as.data.frame(lapply(A$X[(p.con+1):(p.con+p.cat)], factor))


  # if noise variables exist
  if (p.noise.con != 0 | p.noise.cat!=0){
    # scale continuous continuous noise variables
    A$X[(p.con+p.cat+1):(p.con+p.cat+p.noise.con)] <- as.vector(scale(A$X[(p.con+p.cat+1):(p.con+p.cat+p.noise.con)]))
    # discretize to make categorical noise variable
    A$X[(p.con+p.cat+p.noise.con+1):(p.con+p.cat+p.noise.con+p.noise.cat)] <- do.call(data.frame, Map(function(x) cut(x, breaks = quantile(x, probs = c(0, 0.2, 0.4, 0.7, 1)), labels = FALSE, include.lowest = TRUE),
                                                                                                      A$X[(p.con+p.cat+p.noise.con+1):(p.con+p.cat+p.noise.con+p.noise.cat)]))
    A$X[(p.con+p.cat+p.noise.con+1):(p.con+p.cat+p.noise.con+p.noise.cat)] <- as.data.frame(lapply(A$X[(p.con+p.cat+p.noise.con+1):(p.con+p.cat+p.noise.con+p.noise.cat)], factor))

  }

  names(A$X) <- paste0("V", 1:(p.con+p.cat+p.noise.con+p.noise.cat))
  output <- list("SimDat" = A$X, "TrueID" = A$id, "Q" = Q)

}





