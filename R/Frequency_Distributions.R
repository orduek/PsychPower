
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
