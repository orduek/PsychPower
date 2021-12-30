
# loading neccessary libraries
library("tidyverse")
#library("ReIns")
library(poweRlaw)

## Define binarization function
# Function 1
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



## Function 2
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


## Function 3
#' Plot phenotypes
#'
#' Using the frequency calculated with the get_freq() function, plotts the N most common phenotypes
#' @param data The dataset returned by the get_freq() function (including a "freq" column)
#' @param number_Phenotypes The number of symptom profiles to plott. Deafult = 100
#' @return a plot
#' @examples
#' data_f <- get_freq(df);
#' @export
plot_pheno <- function(data, number_Phenotypes=100) {

  # order the data frame
  data2 <- data %>%
    arrange(desc(freq))
## plotting most common phenotypes
  freq1_top  <- data2 %>%
    top_n(freq, n = number_Phenotypes) %>%
    select(freq)

  # The frequency of the fifty most common symptom combinations
  g <- ggplot(freq1_top, aes(x=as.factor(1:nrow(freq1_top)),y=freq)) +
    geom_hline(yintercept = c((median(data$freq)), (max(freq1_top$freq))), color = "grey", size = 0.3) + #max and median
    geom_bar(stat = "identity",fill = "grey26") +
    xlab(" ") +
    ylab("Number of endorsements") +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(breaks=c(round((median(data$freq)), 0), round((max(freq1_top$freq)), 0)))  #max and median

  return(g)
}



## Function 4
#' Power Law analysis
#'
#' Using the frequency calculated with the get_freq() function test wether this is distributed powerlaw or not
#' @param data The data set returned by the get_freq() function (including a "freq" column)
#' @param nThreads The number cpus for the bootsraping (default=1)
#' @param rSeed the seed for randomization (default = 123)
#' @param plot whether you want the function to return a plot of the calculation as well, default = true
#' @param bootStrap whether to run the bootstrap or not
#' @return a matrix including the different measures of the powerlaw (x-min, alpha, p value and SDs)
#' @details This function basically runs the powerLaw package basic assessment. Using the bootstrap allows to receive a p value of the chances this adheres to powerlaw or not. If p <0.05 we assume it adheres to powerlaw (from specific xmin value)
#' @examples
#' a <- test_pl(df);
#' @export
test_pl <-
  function(data,nSims=1000, nThreads=1, rSeed = 123, plot=T, bootStrap=T) {

  #### Prepare
  Distribution <- data$freq
  ### Power Law
  m_pl = displ$new(Distribution)
  est_pl = estimate_xmin(m_pl)
  m_pl$setXmin(est_pl)

  if (bootStrap==T) {
  ## Bootstrap parameters
  ## Test whether power law is possible
    bs_p = bootstrap_p(m_pl, no_of_sims = nSims, threads = nThreads, seed = rSeed)

    RESULTS <- matrix(ncol = 5, nrow = 1)
    colnames(RESULTS) <- c("Xmin", "alpha", "BootP", "SDxmin", "SDalpha")
    RESULTS[1,1] <- m_pl$xmin
    RESULTS[1,2] <- m_pl$pars
    RESULTS[1,3] <- bs_p$p
    RESULTS[1,4] <- sd(bs_p$bootstraps$xmin)
    RESULTS[1,5] <- sd(bs_p$bootstraps$pars)

    if (plot==T) {
      plot(bs_p)}

  } else {

    RESULTS <- matrix(ncol = 2, nrow = 1)
    colnames(RESULTS) <- c("Xmin", "alpha")
    RESULTS[1,1] <- m_pl$xmin
    RESULTS[1,2] <- m_pl$pars
  }




  return(RESULTS)
}
