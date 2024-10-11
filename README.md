
# ADPclustmix

<!-- badges: start -->
<!-- badges: end -->

## Introduction
ADPclustmix (Fast Clustering Using Adaptive Density Peak Detection For Mixed-type Data) is a non-iterative procedure that clusters high dimensional data by finding cluster centers from estimated density peaks developed for mixed-type data. It incorporates multivariate local Gaussian density estimation. The number of clusters as well as bandwidths can either be selected by the user or selected automatically through an internal clustering criterion. Multiple distance and density measures are available for clustering.

## Installation

You can install the most recent version version of ADPclustmix from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Inkaa8/ADPclustmix")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ADPclustmix)
res <- adpclustmix(dat_demo1, dist.method = "gower", dens.method = "NB",con.den = "locfit", idnum = 1:2, idcat = 3:4)
summary(res)
plot(res)
```

