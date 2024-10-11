#' @title Categorical Data Density Calculation
#'
#' @description This is the function to calculate Naive Bayes density for categorical data.
#'
#' @param data a data frame where rows are observations and columns are factorized categorical variables.
#'
#' @return a list of density value and original data.
#'
#' @importFrom dplyr filter select count
#'
#'
est.cat.dens <- function(data) {
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }

  total_count <- nrow(data)

  # count the occurrences of each unique combination
  counts <- data %>%
    dplyr::count(dplyr::across(dplyr::everything()))

  # calculate the probabilities by normalizing the counts
  counts$Probability <- counts$n / total_count

  # drop the frequency column
  probabilities <- counts[, -which(names(counts) == "n")]

  # create a key for merging with the original data
  counts$key <- apply(counts[, 1:ncol(data)], 1, paste, collapse = "_")
  data$key <- apply(data, 1, paste, collapse = "_")

  # merge the probabilities back into the original data
  data_with_probabilities <- merge(data, counts[, c("key", "Probability")], by = "key")

  # drop the key column
  data_with_probabilities <- data_with_probabilities[, -which(names(data_with_probabilities) == "key")]

  results <- list(probabilities=probabilities, data_with_probabilities=data_with_probabilities)
  return(results)

}



#' @title Naive Bayes Density Calculation for Mixed-type Data
#'
#' @description This is the function to calculate Naive Bayes density for mixed-type data.
#'
#' @param data a data frame to be clustered, mixed type variables are allowed. Categorical variables should be factorized.
#' @param idnum a vector indicating the position of continuous variables.
#' @param idcat a vector indicating the position of categorical variables.
#' @param con.den method to calculate density of continuous variables. Options: \code{"logspline"} or \code{"locfit"}. Default is ""locfit".
#'
#' @return a vector of density value.
#'
#' @import logspline
#' @import locfit
#' @importFrom magrittr %>%
#'

NB.cal <- function(data, idnum, idcat, con.den="locfit"){
  NB.den.con <- c()
  if (con.den == "logspline"){
    # calculate density of continuous variables using logspline
    for (i in idnum){
      uni.den <- logspline::dlogspline(data[,i], logspline(data[,i]))
      NB.den.con <- cbind(NB.den.con, uni.den)

    }
  }else if(con.den == "locfit"){
    # calculate density of continuous variables using locfit
    for (i in idnum){
      uni.den <- locfit::density.lf(data[,i], ev = data[,i])$y
      NB.den.con <- cbind(NB.den.con, uni.den)

    }
  }else{
    stop("Specify correct method to calcualte density for continuous variables.")
  }


  # calculate density of categorical variables by joint probability distribution
  NB.den.cat <- est.cat.dens(data[,idcat])$data_with_probabilities$Probability

  NB.den <- cbind(NB.den.con, NB.den.cat)

  NB.dens <- apply(NB.den, 1, prod)

  return(NB.dens)

}


