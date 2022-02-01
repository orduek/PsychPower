utils::globalVariables(c("freq", "m_pl", "m_ln", "m_ex", "sd", "x", "y", ".x"))

### Function 1 #####################################################################################
#' Binarization
#'
#' Binarizes vectors or columns based on a cut-off into 0 and 1.
#' @param data A dataset containing only the columns which are to be binarized
#' @param cut_off A cutoff used to binarize.
#'  Values ≤ cut-off will be binarized to 0,
#'  Values > cut-off will be binarized to 1.
#' @return A dataframe containing the original columns named VN (with N: 1 -> number of columns)
#' and binarized columns named v_binN (with N: 1 -> number of columns).
#' @importFrom dplyr case_when
#' @examples
#' \dontrun{
#' data_bin <- binarize(df, 2)
#' data_bin
#' }
#' @export
binarize <- function(data, cut_off) {

  # Check whether data contains characters or logical
  if (any(sapply(data, is.character) == TRUE) |
    any(sapply(data, is.logical) == TRUE)) {
    stop("Character or logical variables are not allowed")
  } else {

    ## Rename variables to "V1:VN" for binarization
    colnames(data) <- c(paste0("V", 1:(ncol(data))))

    for (i in 1:(ncol(data))) {
      orig <- paste("v_bin", i, sep = "")
      bin <- paste("V", i, sep = "")
      data[orig] <- case_when(
        data[bin] <= cut_off ~ 0,
        data[bin] > cut_off ~ 1
      )
    }
    return(data)
  }
}



### Function 2 #####################################################################################
#' Determine phenotype frequency
#'
#' Determines the frequency of each unique phenotype in the dataframe
#' @param data The dataset, containing v_bin columns (of binarized set)
#' @param target_columns The columns which include the variables of the phenotype.
#' default = starts_with("v_bin"), using the output of \code{binarize}
#' @return A dataframe containing the frequency of each phenotype
#' and the sum of the values of all the variables constituting the phenotype
#' The output has as many rows as phenotypes.
#'
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom  plyr count
#'
#' @examples
#' \dontrun{
#' data_f <- pheno_frequency(df, target_columns = 1:12)
#' data_f
#' }
#' @export
pheno_frequency <- function(data, target_columns = starts_with("v_bin")) {

  # Check whether data contains characters or logical
  if (any(sapply(data, is.character) == TRUE) |
    any(sapply(data, is.logical) == TRUE)) {
    stop("Character or logical variables are not allowed")
  } else {

    # first grab only v_bin columns
    data1 <- data %>% select(all_of(target_columns))

    ## Count frequency of profiles
    data2_counted <- data1 %>% plyr::count()

    # Create sum score of endorsed symptoms
    data2_counted <- data2_counted %>%
      mutate(total_bin = rowSums(data2_counted) - freq)

    return(data2_counted)
  }
}



### Function 3 #####################################################################################
#' Plot phenotypes
#'
#' \code{plot_Pheno} plots the frequency of each unique phenotype in the sample in descending order
#' @param data a dataframe or matrix
#' @param frequency column indicating phenotype frequency, default = "freq"
#' @param n_phenotypes The number of phenotypes plotted, default = 100
#' @param color The color of the phenotypes
#' @return a plot
#'
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' fig1 <- plot_pheno(df, frequency = "freq", n_phenotypes = 100, color = "grey26")
#' fig1
#' }
#' @export
plot_pheno <- function(data, frequency = "freq", n_phenotypes = 100, color = "grey26") {

  # Check whether frequency is identified
  if (!any(names(data) == frequency)) {
    stop("Frequency colum not identified")
  } else {

    # Identify frequency column
    names(data)[names(data) == frequency] <- "freq"

    # Select & order phenotypes to be plotted
    data <- data %>%
      arrange(desc(freq)) %>%
      top_n(freq, n = n_phenotypes) %>%
      select(freq)

    # Plot
    g <- ggplot(data, aes(x = as.factor(1:nrow(data)), y = freq)) +
      geom_bar(stat = "identity", fill = color) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      xlab("") +
      ylab("") +
      theme_minimal() +
      theme(
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = .2, color = "black"),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)
      )

    return(g)
  }
}



### Function 4 #####################################################################################
#' Describe phenotypes in sample
#'
#' Describes basic information about the distribution of the phenotype frequency.
#' Uses the frequency calculated with the pheno_frequency() function.
#' @param data a dataframe or matrix
#' @param frequency column indicating phenotype frequency, default = "freq"
#' @return a matrix indicating
#' the number of unique phenotypes
#' the frequency of the most common phenotype
#' the median frequency of all phenotypes
#' @importFrom  stats median
#' @examples
#' \dontrun{
#' description <- describe_pheno(df, frequency = "freq")
#' description
#' }
#' @export
describe_pheno <- function(data, frequency = "freq") {

  # Check whether frequency is identified
  if (!any(names(data) == frequency)) {
    stop("Frequency colum not identified")
  } else {

    # Identify frequency column
    names(data)[names(data) == frequency] <- "freq"

    # Prepare output
    res <- matrix(ncol = 1, nrow = 3)
    rownames(res) <- c("Unique Phenotypes", "N most frequent", "Median frequency")
    colnames(res) <- c("Number")

    res[1, 1] <- nrow(data) # Number of unique phenotypes
    res[2, 1] <- max(data$freq) # Number of endorsements of most common phenotype
    res[3, 1] <- median(data$freq) # Number of endorsements of most common phenotype

    return(res)
  }
}



### Function 5 #####################################################################################
#' Show most common phenotypes
#'
#' Prints the value of each variable of the most common phenotypes
#' @param data a dataframe or matrix
#' @param frequency column indicating phenotype frequency, default = "freq"
#' @param n_phenotypes , number of phenotypes shown, default = 5
#' @return print
#' @importFrom dplyr arrange
#'
#' @examples
#' \dontrun{
#' common_pheno(df, frequency = "freq", n_phenotypes = 5)
#' }
#' @export
common_pheno <- function(data, frequency = "freq", n_phenotypes = 5) {

  # Check whether frequency is identified
  if (!any(names(data) == frequency)) {
    stop("Frequency colum not identified")
  } else {

    # Identify frequency column
    names(data)[names(data) == frequency] <- "freq"

    data <- data %>%
      arrange(desc(freq))

    for (i in 1:n_phenotypes) {
      print(data[i, ])
    }
  }
}


