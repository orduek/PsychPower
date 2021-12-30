
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
#' data_bin <- binarizze_data(df, 2);
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
