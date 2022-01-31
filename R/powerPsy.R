utils::globalVariables(c("freq", "m_pl", "m_ln", "m_ex"))

### Function 1 #####################################################################################
#' Binarization
#'
#' Binarizes vectors or columns based on a cut-off into 0 and 1.
#' @param data A dataset containing only the columns which are to be binarized
#' @param cut_off A cutoff used to binarize.
#'  Values ≤ cut-off will be binarized to 0,
#'  Values > cut-off will be binarized to 1.
#' @return A dataframe containing the original columns named VN (with N: 1 -> number of columns)
#' and binarized columns named v_binN (with N: 1 -> number of columns)
#' @examples
#' \dontrun{
#' data_bin <- binarize(df, 2)
#' data_bin
#'  }
#' @export
binarize <- function(data, cut_off) {

  # Check whether data contains characters or logical
  if(any(sapply(data, is.character) == TRUE) |
     any(sapply(data, is.logical) == TRUE)){
    stop("Character or logical variables are not allowed")
  } else {

    ## Rename variables to "V1:VN" for binarization
    colnames(data) <- c(paste0("V", 1:(ncol(data))))

    for (i in 1:(ncol(data))) {
      orig <- paste("v_bin", i, sep = "")
      bin <- paste("V", i, sep = "")
      data[orig] <- dplyr::case_when(
        data[bin] <= cut_off ~ 0,
        data[bin] > cut_off ~ 1)
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
#' @import dplyr
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
  if(any(sapply(data, is.character) == TRUE) |
     any(sapply(data, is.logical) == TRUE)){
    stop("Character or logical variables are not allowed")
  } else {

    # first grab only v_bin columns
    data1 <- data %>% dplyr::select(all_of(target_columns))

    ## Count frequency of profiles
    data2_counted <- data1 %>% plyr::count()

    # Create sum score of endorsed symptoms
    data2_counted <- data2_counted %>%
      dplyr::mutate(total_bin = rowSums(data2_counted) - freq)

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
  if(!any(names(data) == frequency)){
    stop("Frequency colum not identified")
  } else {

    # Identify frequency column
    names(data)[names(data) == frequency] <- "freq"

    # Select & order phenotypes to be plotted
    data  <- data %>%
      dplyr::arrange(desc(freq)) %>%
      dplyr::top_n(freq, n = n_phenotypes) %>%
      dplyr::select(freq)

    # Plot
    g <- ggplot2::ggplot(data, aes(x=as.factor(1:nrow(data)),y=freq)) +
      ggplot2::geom_bar(stat = "identity",fill = color) +
      ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      ggplot2::xlab("") +
      ggplot2::ylab("") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=.2, color="black" ),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

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
  if(!any(names(data) == frequency)){
    stop("Frequency colum not identified")
  } else {

    # Identify frequency column
    names(data)[names(data) == frequency] <- "freq"

    # Prepare output
    res <- matrix(ncol = 1, nrow = 3)
    rownames(res) <- c("Unique Phenotypes", "N most frequent", "Median frequency")
    colnames(res) <- c("Number")

    res[1,1] <- nrow(data) # Number of unique phenotypes
    res[2,1] <- max(data$freq) # Number of endorsements of most common phenotype
    res[3,1] <- median(data$freq) # Number of endorsements of most common phenotype

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
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' common_pheno(df, frequency = "freq", n_phenotypes = 5)
#' }
#' @export
common_pheno <- function(data, frequency = "freq", n_phenotypes = 5) {

  # Check whether frequency is identified
  if(!any(names(data) == frequency)){
    stop("Frequency colum not identified")
  } else {

    # Identify frequency column
    names(data)[names(data) == frequency] <- "freq"

    data <- data %>%
      dplyr::arrange(desc(freq))

    for (i in 1:n_phenotypes) {
      print(data[i,])
    }
  }
}


### Function 6 #####################################################################################
#' Frequency distribution analysis
#'
#' This function performs poweRlaw functions and gives back best fitting distributions.
#' @param data a dataframe or matrix
#' @param frequency column indicating phenotype frequency, default = "freq"
#' @return a list of distribution objects.
#' Includes approximations of the best fitting power-law, log normal, and exponential distribution
#' @details Approximates different distributions to the phenotype frequency distribution
#' using the poweRlaw package applying methods by Clauset et al. (2011).
#  The distributions approximated include: log-normal, power-law, exponential.
#' Uses the phenotype frequency calculated with the pheno_frequency() function.
#' @import poweRlaw
#' @examples
#' \dontrun{
#' distributions <- pheno_distributions(df, frequency = "freq")
#' }
#' @export
pheno_distributions <-
  function(data, frequency = "freq") {

    # Check frequency
    if(!any(names(data) == frequency)){
      stop("Frequency colum not identified")
    } else {

      ### Power Law
      m_pl <- poweRlaw::displ$new(data$freq)
      est_pl <- poweRlaw::estimate_xmin(m_pl)
      m_pl$setXmin(est_pl)

      ## Log normal with Xmin of PL
      m_ln <- poweRlaw::dislnorm$new(data$freq)
      m_ln$setXmin(m_pl$getXmin())
      est_m_ln <- poweRlaw::estimate_pars(m_ln)
      m_ln$setPars(est_m_ln)

      ## Exponential with Xmin of PL
      m_ex <- poweRlaw::disexp$new(data$freq)
      m_ex$setXmin(m_pl$getXmin())
      est_m_ex <- poweRlaw::estimate_pars(m_ex)
      m_ex$setPars(est_m_ex)

      res <- list()
      res[[1]] <- m_pl
      res[[2]] <- m_ln
      res[[3]] <- m_ex
      return(res)
    }
  }










### Function 7 #####################################################################################
#' Parameters of distributions
#'
#' Using the frequency calculated with the get_freq() function test whether this is distributed power law or not
#' @param data The data set returned by the get_freq() function (including a "freq" column)
#' @param nSims The number of simulations used in the bootstrapping
#' @param nThreads The number CPUs for the bootstrapping (default=1)
#' @param bootStrap whether to run the bootstrap or not
#' @param rSeed the seed for randomization (default = 123)
#' @return a matrix with the parameters of the approximations of the best fitting power-law, log normal, and exponential distribution
#' @details This function displays the parameters of the approximations of the best fitting power-law, log normal, and exponential distribution
#' @examples
#' \dontrun{
#' a <- describe_pheno_distr(df)
#' }
#' @export
describe_pheno_distr <-
  function(data, nSims=1000, nThreads=1, bootStrap=T, rSeed = 123) {

    if (bootStrap==T) {
      ## Bootstrap parameters
      bs_p_pl = poweRlaw::bootstrap_p(data[[1]], no_of_sims = nSims, threads = nThreads, seed = rSeed)
      bs_p_ln = poweRlaw::bootstrap_p(data[[2]], no_of_sims = nSims, threads = nThreads, seed = rSeed)
      bs_p_ex = poweRlaw::bootstrap_p(data[[3]], no_of_sims = nSims, threads = nThreads, seed = rSeed)

      res <- matrix(ncol = 9, nrow = 3)
      rownames(res) <- c("PowerLaw", "LogNormal", "Exponential")
      colnames(res) <- c("Xmin", "alpha/exponent", "log(mu)", "log(sigma)", "Boot_P_Value", "SDxmin",
                             "SDalpha/exponent", "log(mu)", "log(sigma)")
      res[1,1] <- data[[1]]$xmin
      res[1,2] <- data[[1]]$pars
      res[1,3] <- ""
      res[1,4] <- ""
      res[1,5] <- bs_p_pl$p
      res[1,6] <- sd(bs_p_pl$bootstraps$xmin)
      res[1,7] <- sd(bs_p_pl$bootstraps$pars)
      res[1,8] <- ""
      res[1,9] <- ""

      res[2,1] <- data[[2]]$xmin
      res[2,2] <- ""
      res[2,3] <- data[[2]]$pars[1]
      res[2,4] <- data[[2]]$pars[2]
      res[2,5] <- bs_p_ln$p
      res[2,6] <- sd(bs_p_ln$bootstraps$xmin)
      res[2,7] <-
        res[2,8] <- sd(bs_p_ln$bootstraps$pars1)
      res[2,9] <- sd(bs_p_ln$bootstraps$pars2)

      res[3,1] <- data[[3]]$xmin
      res[3,2] <- data[[3]]$pars
      res[3,3] <- ""
      res[3,4] <- ""
      res[3,5] <- bs_p_ex$p
      res[3,6] <- sd(bs_p_ex$bootstraps$xmin)
      res[3,7] <- sd(bs_p_ex$bootstraps$pars)
      res[3,8] <- ""
      res[3,9] <- ""

      if (plot==T) {
        plot(bs_p_pl)
        plot(bs_p_ln)
        plot(bs_p_ex)}

    } else {

      res <- matrix(ncol = 4, nrow = 3)
      colnames(res) <- c("Xmin", "alpha/exponent", "log(mu)", "log(sigma)")
      rownames(res) <- c("PowerLaw", "LogNormal", "Exponential")
      res[1,1] <- data[[1]]$xmin
      res[1,2] <- data[[1]]$pars
      res[1,3] <- ""
      res[1,4] <- ""

      res[2,1] <- data[[2]]$xmin
      res[2,2] <- ""
      res[2,3] <- data[[2]]$pars[1]
      res[2,4] <- data[[2]]$pars[2]

      res[3,1] <- data[[3]]$xmin
      res[3,2] <- data[[3]]$pars
      res[3,3] <- ""
      res[3,4] <- ""
    }

    return(res)
}



### Function 8 #####################################################################################
#' Compare distributions
#'
#' @param data Parameters of the approximations of the best fitting power-law, log normal, and exponential distribution
#' @return a matrix with the parameters of the comparisons
#' @details This function compares parameters of the approximations of the best fitting power-law, log normal, and exponential distribution and outputs the p-values as a matrix
#' @examples
#' \dontrun{
#' a <- compare_pheno_distr(df)
#'  }
#' @export
compare_pheno_distr <-
  function(data) {

    res <- matrix(ncol = 3, nrow = 3)
    colnames(res) <- c("PowerLaw", "LogNormal", "Exponential")
    rownames(res) <- c("PowerLaw", "LogNormal", "Exponential")
    res[1,1] <- ""
    res[1,2] <- round(poweRlaw::compare_distributions(m_ln, m_pl)$p_one_sided, 4)
    res[1,3] <- round(poweRlaw::compare_distributions(m_ex, m_pl)$p_one_sided, 4)
    res[2,1] <- round(poweRlaw::compare_distributions(m_pl,    m_ln)$p_one_sided, 4)
    res[2,2] <- ""
    res[2,3] <- round(poweRlaw::compare_distributions(m_ex, m_ln)$p_one_sided, 4)
    res[3,1] <- round(poweRlaw::compare_distributions(m_pl,    m_ex)$p_one_sided, 4)
    res[3,2] <- round(poweRlaw::compare_distributions(m_ln, m_ex)$p_one_sided, 4)
    res[3,3] <- ""

    return(res)
  }



### Function 9 #####################################################################################
#' Plot distributions
#'
#' @param data Parameters of the approximations of the best fitting power-law, log normal, and exponential distribution
#' @param limity limits of the y-axis
#' @param limitx limits of the x-axis
#' @return a plot
#' @details This function plots parameters of the approximations of the best fitting power-law, log normal, and exponential distribution
#'
#' @import ggplot2
#' @import cowplot
#' @import scales
#'
#' @examples
#' \dontrun{
#' a <- plot_pheno_distr(df)
#' }
#' @export
plot_pheno_distr<-
  function(data, limity = 10^-4, limitx = 10^3) {

    res_pl <- plot(m_pl)
    line_pl <- lines(m_pl)
    line_ln <- lines(m_ln)
    line_ex <- lines(m_ex)

    ## plot
    p <- ggplot2::ggplot(res_pl, aes(x=x,y=y)) +
      geom_point(size = 1)+
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", math_format(10^.x)),
                    expand = c(0, 0),
                    limits = c(limity, 1)) +
      scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", math_format(10^.x)),
                    expand = c(0, 0),
                    limits = c(1, limitx)) +
      geom_line(data = line_pl, aes(x=x, y=y), color = "red", size = 1) +
      geom_line(data = line_ln, aes(x=x, y=y), color = "blue", size = 1,linetype = "dashed")+
      geom_line(data = line_ex, aes(x=x, y=y), color = "orange", size = 1,linetype = "twodash")+
      xlab("") +
      ylab("") +
      ggtitle("")+
      theme_minimal_grid() +
      theme(
        plot.title = element_text(size=11),
        axis.title.x = element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size=9, color = "black"),
        axis.text.y = element_text(size=9, color = "black", margin = margin(t = 0, r = 0, b = 0, l = 5)),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_line(size=.2, color="black"),
        panel.grid.major.y = element_line(size=.2, color="black"),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

    return(p)
  }
