utils::globalVariables(c("freq", "m_pl", "m_ln", "m_ex", "sd", "x", "y", ".x", "."))

### 0. Helper functions ############################################################################
# 0.1 Check for character or logical columns
check_data_type <- function(data) {
  if (any(purrr::map_lgl(data, ~is.character(.x) | is.logical(.x)))) {
    stop("Character or logical variables are not allowed")
  }
}

# 0.2 Set frequency column name to "freq"
set_freq_colname <- function(data, frequency) {
  if (!any(names(data) == frequency)) {
    stop("Frequency column not identified")
  }
  names(data)[names(data) == frequency] <- "freq"
  return(data)
}


### Function 1 #####################################################################################
#' Binarization
#'
#' Binarizes vectors or columns based on a cut-off into 0 and 1.
#' @param data A dataset containing only the columns which are to be binarized
#' @param cut_off A cutoff used to binarize.
#'  Values <= cut-off will be binarized to 0,
#'  Values > cut-off will be binarized to 1.
#' @return A dataframe containing the original columns named VN (with N: 1 -> number of columns)
#' and binarized columns named v_binN (with N: 1 -> number of columns).
#' @importFrom dplyr case_when
#' @examples
#' \dontrun{
#' data_binarized <- binarize(data_test, cut_off = 2)
#' data_binarized
#' }
#' @export
binarize <- function(data, cut_off) {
  # Check whether data contains characters or logical
  check_data_type(data)

  ## Rename variables to "V1:VN" for binarization
  new_colnames <- paste0("V", seq_along(data))
  colnames(data) <- new_colnames

  bin_names <- paste0("v_bin", seq_along(data))
  data[bin_names] <- ifelse(data <= cut_off, 0, 1)

  return(data)
}


### Function 2 #####################################################################################
#' Determine symptom combinations frequency
#'
#' Determines the frequency of each unique symptom combinations in the dataframe
#' @param data The dataset, containing v_bin columns (of binarized set)
#' @param target_columns The columns which include the variables of the symptom combinations
#' default = starts_with("v_bin"), using the output of \code{binarize}
#' @return A dataframe containing the frequency of each symptom combinations
#' and the sum of the values of all the variables constituting the symptom combinations
#' The output has as many rows as symptom combinations.
#'
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom tidyselect starts_with
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' data_frequency <- pheno_frequency(data_binarized, target_columns = tidyselect::starts_with("v_bin"))
#' data_frequency
#' }
#' @export
pheno_frequency <- function(data, target_columns = tidyselect::starts_with("v_bin")) {
  # Check whether data contains characters or logical
  check_data_type(data)

  # First grab only v_bin columns
  data_selected <- data %>% dplyr::select(all_of(target_columns))

  ## Count frequency of profiles
  data_counted <- data_selected %>%
    dplyr::group_by(across(everything())) %>%
    dplyr::summarise(freq = n(), .groups = "drop")

  # Create sum score of endorsed symptoms
  data_counted <- data_counted %>%
    dplyr::mutate(total_bin = rowSums(.) - freq)

  return(data_counted)
}



### Function 3 #####################################################################################
#' Plot symptom combinations
#'
#' \code{plot_pheno} plots the frequency of each unique symptom combinations in the sample in descending order
#' @param data a dataframe or matrix
#' @param frequency column indicating symptom combination frequency, default = "freq"
#' @param n_phenotypes The number of symptom combinations plotted, default = 100
#' @param color The color of the symptom combinations
#' @return a plot
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' fig1 <- plot_pheno(data_frequency, frequency = "freq", n_phenotypes = 100, color = "grey26")
#' fig1
#' }
#' @export
plot_pheno <- function(data, frequency = "freq", n_phenotypes = 100, color = "grey26") {
  # Check whether frequency is identified
  data <- set_freq_colname(data, frequency)

  # Select & order symptom combinations to be plotted
  data <- data %>%
    dplyr::arrange(desc(freq)) %>%
    dplyr::slice_head(n = n_phenotypes) %>%
    dplyr::select(freq)

  max_y <- max(data$freq) * 1.05 # Increase by 5% to simulate expansion

  g <- ggplot2::ggplot(data, aes(x = as.factor(1:nrow(data)), y = freq)) +
    ggplot2::geom_bar(stat = "identity", fill = color) +
    ggplot2::ylim(0, max_y) +  # Set y-axis limits here
    labs(x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = .2, color = "black"),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )
  return(g)
}




### Function 4 #####################################################################################
#' Describe symptom combinations in sample
#'
#' Describes basic information about the distribution of the symptom combinations frequency.
#' Uses the frequency calculated with the pheno_frequency() function.
#' @param data a dataframe or matrix
#' @param frequency column indicating symptom combinations frequency, default = "freq"
#' @return a matrix indicating
#' the number of unique symptom combinations
#' the frequency of the most common symptom combinations
#' the median frequency of all symptom combinations
#'
#' @importFrom  stats median
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' desc_pheno <- describe_pheno(data_frequency, frequency = "freq")
#' desc_pheno
#' }
#' @export
describe_pheno <- function(data, frequency = "freq") {
  # Check whether frequency is identified
  data <- set_freq_colname(data, frequency)

  # Prepare output
  matrix_data <- c(
    nrow(data), # Number of unique symptom combinations
    max(data$freq), # Frequency of most common symptom combinations
    stats::median(data$freq) # Median frequency
  )

  res <- matrix(matrix_data, ncol = 1, nrow = 3)
  rownames(res) <- c("Unique Symptom combinations", "N most frequent", "Median frequency")
  colnames(res) <- c("Number")

  return(res)
}


### Function 5 #####################################################################################
#' Show most common symptom combinations
#'
#' Prints the value of each variable of the most common symptom combinations
#' @param data a dataframe or matrix
#' @param frequency column indicating symptom combinations frequency, default = "freq"
#' @param n_phenotypes , number of symptom combinations shown, default = 5
#' @return print
#'
#' @importFrom dplyr arrange
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' common_pheno(data_frequency, frequency = "freq", n_phenotypes = 5)
#' }
#' @export
common_pheno <- function(data, frequency = "freq", n_phenotypes = 5) {

   # Check whether frequency is identified
  data <- set_freq_colname(data, frequency)

  top_data <- data %>% dplyr::arrange(desc(freq)) %>% dplyr::slice_head(n = n_phenotypes)

  return(top_data)
}
