
# loading neccessary libraries
library("tidyverse")
library("readxl")
library("ReIns")
library("psych")

## Define binarization function

#' Binarization of the data
#'
#' Binarizes the data from likert scale to 0s and 1s
#' @param data The dataset, containing only columns needed to be binarised
#' @param cut_off the cutoff to which we would use to binarise (if its eqal or smaller than cutoff than it will be coded as zero)
#' @return A dataframe containing original columns (as v) and binarized columns (as v_bin)
#' @examples
#' data_bin <- binarize_data(df, 2);
#' @export
binarize_data <-
  function(data, cut_off) {
    data_bin <- data

    ## Rename variables to "Q1:QN" for binarization
    colnames(data_bin) <- c(paste0("V", 1:(ncol(data_bin))))

    for (i in 1:(ncol(data_bin))){
      orig <- paste("v_bin", i, sep = "")
      bin <- paste("V", i, sep = "")
      data_bin[orig] <- dplyr::case_when(data_bin[bin] <= cut_off ~ 0,
                                         data_bin[bin] > cut_off ~ 1)
    }
    return(data_bin)
  }



## Function 2 -  original code
#' Count frequency of symptoms' profiles
#'
#' Using only binarized data  - counting the frequency of the each profile
#' @param data The dataset, containing v_bin columns (of binarized set)
#' @return A dataframe containing frequency of each profile and the sum of scores (from binarized data) of each of these profiles
#' @examples
#' data_f <- get_freq(df);
#' @export
get_freq <- function(data){

  # first grab only v_bin columns
  data1 <- select(data, starts_with('v_bin'))

  ## Count frequency of profiles
  data2_counted <- plyr::count(data1[, ])

  # Create sum score of endorsed symptoms
  data2_counted <- data2_counted %>%
    mutate(total_bin = rowSums(data2_counted)-freq)

  return(data2_counted)
}
