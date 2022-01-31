## Define binarization function
# Function 1
#' Binarization of the data
#'
#' Binarizes the data from likert scale to 0s and 1s
#' @param data The dataset, containing only columns needed to be binarised
#' @param cut_off the cutoff to which we would use to binarise (if its eqal or smaller than cutoff than it will be coded as zero)
#' @return A dataframe containing original columns (as v) and binarized columns (as v_bin)
#' @examples
#' # dont run data_bin <- binarize_data(df, 2)
#' @export
binarize_data <-
  function(data, cut_off) {
    data_bin <- data

    ## Rename variables to "Q1:QN" for binarization
    colnames(data_bin) <- c(paste0("V", 1:(ncol(data_bin))))

    for (i in 1:(ncol(data_bin))) {
      orig <- paste("v_bin", i, sep = "")
      bin <- paste("V", i, sep = "")
      data_bin[orig] <- dplyr::case_when(
        data_bin[bin] <= cut_off ~ 0,
        data_bin[bin] > cut_off ~ 1
      )
    }
    return(data_bin)
  }



## Function 2
#' Count phenotype frequency
#'
#' Using only binarized data  - counting the frequency of the each profile
#' @param data The dataset, containing v_bin columns (of binarized set)
#' @return A dataframe containing frequency of each profile and the sum of scores (from binarized data) of each of these profiles
#'
#' @import dplyr
#' @importFrom  plyr count
#' @examples
#' data_f <- get_freq(df)
#' @export
Pheno_frequency <- function(data) {

  # first grab only v_bin columns
  data1 <- dplyr::select(data, starts_with("v_bin"))

  ## Count frequency of profiles
  data2_counted <- plyr::count(data1[, ])

  # Create sum score of endorsed symptoms
  data2_counted <- data2_counted %>%
    dplyr::mutate(total_bin = rowSums(data2_counted) - freq)

  return(data2_counted)
}


## Function 3
#' Plot phenotypes
#'
#' Plotting the frequency of each phenotype in the sample in descending order
#' @param data a dataframe or matrix
#' @param freq column indicating phenotype frequency, default = "freq"
#' @param n_Pheno The number of phenotypes plotted, default = 100
#' @param color The color of the phenotypes
#' @return a plot
#'
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' data_f <- plot_Pheno(df)
#' @export
plot_Pheno <- function(data, n_Pheno = 100, frequency = "freq", color = "grey26") {

  # Check frequency
  if(!any(names(data_bin_f) == frequency)){
    stop("Frequency colum not identified")
  } else {

    # Identify column
    names(data)[names(data) == frequency] <- "freq"

    # order the data frame
    data2 <- data %>%
      dplyr::arrange(desc(freq))

    ## plotting most common phenotypes
    freq1_top  <- data2 %>%
      dplyr::top_n(freq, n = n_Phenotypes) %>%
      dplyr::select(freq)

    g <- ggplot2::ggplot(freq1_top, aes(x=as.factor(1:nrow(freq1_top)),y=freq)) +
      ggplot2::geom_bar(stat = "identity",fill = color) +
      ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      ggplot2::xlab(" ") +
      ggplot2::ylab("") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = element_text(size=12),
        axis.title.x = element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size=9, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=9, color = "black", margin = margin(t = 0, r = 0, b = 0, l = 5)),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=.2, color="black" ),
        panel.grid.minor.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

    return(g)
  }
}


## Function 5
#' Distribution analysis
#'
#' Using the frequency calculated with the get_freq() function test whether this is distributed powerlaw or not
#' @param data The data set returned by the get_freq() function (including a "freq" column)
#' @param freq column indicating phenotype frequency, default = "freq"
#' @return a list including the approximations of the best fitting power-law, log normal, and exponential distribution
#' @details This function performs poweRlaw functions and gives back best fitting distributions
#'
#' @import poweRlaw
#'
#' @examples
#' a <- pheno_distributions(df)
#' @export
Pheno_distributions <-
  function(data, frequency = "freq") {

    # Check frequency
    if(!any(names(data_bin_f) == frequency)){
      stop("Frequency colum not identified")
    } else {

    #### Prepare
    Distribution <- data$freq

    ### Power Law
    m_pl <- poweRlaw::displ$new(Distribution)
    est_pl <- poweRlaw::estimate_xmin(m_pl)
    m_pl$setXmin(est_pl)

    ## Log normal with Xmin of PL
    m_ln_EQ <- poweRlaw::dislnorm$new(Distribution)
    m_ln_EQ$setXmin(m_pl$getXmin())
    est_m_ln_EQ <- poweRlaw::estimate_pars(m_ln_EQ)
    m_ln_EQ$setPars(est_m_ln_EQ)

    ## Exponential with Xmin of PL
    m_ex_EQ <- poweRlaw::disexp$new(Distribution)
    m_ex_EQ$setXmin(m_pl$getXmin())
    est_m_ex_EQ <- poweRlaw::estimate_pars(m_ex_EQ)
    m_ex_EQ$setPars(est_m_ex_EQ)

    RESULTS <- list()
    RESULTS[[1]] <- m_pl
    RESULTS[[2]] <- m_ln_EQ
    RESULTS[[3]] <- m_ex_EQ
    return(RESULTS)
  }
}


## Function 6
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
#' a <- pheno_distributions_parameters(df)
#' @export
describe_Pheno_distr <-
  function(data, nSims=1000, nThreads=1, bootStrap=T, rSeed = 123) {

    if (bootStrap==T) {
      ## Bootstrap parameters
      bs_p_pl = poweRlaw::bootstrap_p(data[[1]], no_of_sims = nSims, threads = nThreads, seed = rSeed)
      bs_p_ln = poweRlaw::bootstrap_p(data[[2]], no_of_sims = nSims, threads = nThreads, seed = rSeed)
      bs_p_ex = poweRlaw::bootstrap_p(data[[3]], no_of_sims = nSims, threads = nThreads, seed = rSeed)

      RESULTS <- matrix(ncol = 9, nrow = 3)
      rownames(RESULTS) <- c("PowerLaw", "LogNormal", "Exponential")
      colnames(RESULTS) <- c("Xmin", "alpha/exponent", "log(mu)", "log(sigma)", "Boot_P_Value", "SDxmin", "SDalpha/exponent", "log(mu)", "log(sigma)")
      RESULTS[1,1] <- data[[1]]$xmin
      RESULTS[1,2] <- data[[1]]$pars
      RESULTS[1,3] <- ""
      RESULTS[1,4] <- ""
      RESULTS[1,5] <- bs_p_pl$p
      RESULTS[1,6] <- sd(bs_p_pl$bootstraps$xmin)
      RESULTS[1,7] <- sd(bs_p_pl$bootstraps$pars)
      RESULTS[1,8] <- ""
      RESULTS[1,9] <- ""

      RESULTS[2,1] <- data[[2]]$xmin
      RESULTS[2,2] <- ""
      RESULTS[2,3] <- data[[2]]$pars[1]
      RESULTS[2,4] <- data[[2]]$pars[2]
      RESULTS[2,5] <- bs_p_ln$p
      RESULTS[2,6] <- sd(bs_p_ln$bootstraps$xmin)
      RESULTS[2,7] <-
        RESULTS[2,8] <- sd(bs_p_ln$bootstraps$pars1)
      RESULTS[2,9] <- sd(bs_p_ln$bootstraps$pars2)

      RESULTS[3,1] <- data[[3]]$xmin
      RESULTS[3,2] <- data[[3]]$pars
      RESULTS[3,3] <- ""
      RESULTS[3,4] <- ""
      RESULTS[3,5] <- bs_p_ex$p
      RESULTS[3,6] <- sd(bs_p_ex$bootstraps$xmin)
      RESULTS[3,7] <- sd(bs_p_ex$bootstraps$pars)
      RESULTS[3,8] <- ""
      RESULTS[3,9] <- ""

      if (plot==T) {
        plot(bs_p_pl)
        plot(bs_p_ln)
        plot(bs_p_ex)}

    } else {

      RESULTS <- matrix(ncol = 4, nrow = 3)
      colnames(RESULTS) <- c("Xmin", "alpha/exponent", "log(mu)", "log(sigma)")
      rownames(RESULTS) <- c("PowerLaw", "LogNormal", "Exponential")
      RESULTS[1,1] <- data[[1]]$xmin
      RESULTS[1,2] <- data[[1]]$pars
      RESULTS[1,3] <- ""
      RESULTS[1,4] <- ""

      RESULTS[2,1] <- data[[2]]$xmin
      RESULTS[2,2] <- ""
      RESULTS[2,3] <- data[[2]]$pars[1]
      RESULTS[2,4] <- data[[2]]$pars[2]

      RESULTS[3,1] <- data[[3]]$xmin
      RESULTS[3,2] <- data[[3]]$pars
      RESULTS[3,3] <- ""
      RESULTS[3,4] <- ""
    }

    return(RESULTS)
}



## Function 7
#' Compare distributions
#'
#' @param data Parameters of the approximations of the best fitting power-law, log normal, and exponential distribution
#' @return a matrix with the parameters of the comparisons
#' @details This function compares parameters of the approximations of the best fitting power-law, log normal, and exponential distribution and outputs the p-values as a matrix
#' @examples
#' a <- compare_pheno_distributions(df)
#' @export
compare_Pheno_distr <-
  function(data) {

    RESULTS <- matrix(ncol = 3, nrow = 3)
    colnames(RESULTS) <- c("PowerLaw", "LogNormal", "Exponential")
    rownames(RESULTS) <- c("PowerLaw", "LogNormal", "Exponential")
    RESULTS[1,1] <- ""
    RESULTS[1,2] <- round(poweRlaw::compare_distributions(m_ln_EQ, m_pl)$p_one_sided, 4)
    RESULTS[1,3] <- round(poweRlaw::compare_distributions(m_ex_EQ, m_pl)$p_one_sided, 4)
    RESULTS[2,1] <- round(poweRlaw::compare_distributions(m_pl,    m_ln_EQ)$p_one_sided, 4)
    RESULTS[2,2] <- ""
    RESULTS[2,3] <- round(poweRlaw::compare_distributions(m_ex_EQ, m_ln_EQ)$p_one_sided, 4)
    RESULTS[3,1] <- round(poweRlaw::compare_distributions(m_pl,    m_ex_EQ)$p_one_sided, 4)
    RESULTS[3,2] <- round(poweRlaw::compare_distributions(m_ln_EQ, m_ex_EQ)$p_one_sided, 4)
    RESULTS[3,3] <- ""

    return(RESULTS)
  }



## Function 8
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
#'
#' @examples
#' a <- plot_pheno_distributions(df)
#' @export
plot_Pheno_distr<-
  function(data, limity = 10^-4, limitx = 10^3) {

    res_pl <- plot(m_pl)
    line_pl <- lines(m_pl)
    line_ln <- lines(m_ln_EQ)
    line_ex <- lines(m_ex_EQ)

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
