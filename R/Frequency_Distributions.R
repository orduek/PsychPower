
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
    if (!any(names(data) == frequency)) {
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
#' Parameters of frequency distributions
#'
#' This function gives back the parameters of the best fitting distributions estimated with
#' \code{pheno_distributions}.
#' @param data A list of three "poweRlaw" objects.
#' @param bootStrap whether to run the bootstrap or not (default = TRUE)
#' @param nBoots Number of bootstraps (default = 1000)
#' @param nCores Number of cores to use in computing results (default = 1).
#'
#' @importFrom poweRlaw bootstrap_p
#'
#' @return a matrix with the parameters of the approximations of the best fitting power-law, log normal, and exponential distribution
#' @examples
#' \dontrun{
#' a <- describe_pheno_distr(df, bootStrap = T, nBoots = 5, nCores = 1)
#' }
#' @export
describe_pheno_distr <-
  function(data, bootStrap = T, nBoots = 1000, nCores = 1) {

    # Check list
    if (length(data) != 3) {
      stop("Inputs must consist of 3 poweRlaw objects")
    } else {

      # Check poweRlaw objects 1
      if (class(data[[1]]) != "displ") {
        stop("First element of input must be poweRlaw object *displ* (power-law)")
      } else {

        # Check poweRlaw objects 2
        if (class(data[[2]]) != "dislnorm") {
          stop("Second element of input must be poweRlaw object *dislnorm* (log-normal)")
        } else {

          # Check poweRlaw objects 3
          if (class(data[[3]]) != "disexp") {
            stop("Thir element of input must be poweRlaw object *disexp* (exponential)")
          } else {
            if (bootStrap == T) {
              ## Bootstrap parameters
              bs_p_pl <- bootstrap_p(data[[1]], no_of_sims = nBoots, threads = nCores)
              bs_p_ln <- bootstrap_p(data[[2]], no_of_sims = nBoots, threads = nCores)
              bs_p_ex <- bootstrap_p(data[[3]], no_of_sims = nBoots, threads = nCores)

              res <- matrix(ncol = 3, nrow = 9)
              colnames(res) <- c("PowerLaw", "LogNormal", "Exponential")
              rownames(res) <- c(
                "Xmin", "SDxmin", "alpha/exponent", "SDalpha/exponent",
                "log(mu)", "SDlog(mu)", "log(sigma)", "SDlog(sigma)",
                "Boot_P_Value"
              )
              res[1,1] <- data[[1]]$xmin
              res[2,1] <- sd(bs_p_pl$bootstraps$xmin)
              res[3,1] <- data[[1]]$pars
              res[4,1] <- sd(bs_p_pl$bootstraps$pars)
              res[5,1] <- NA
              res[6,1] <- NA
              res[7,1] <- NA
              res[8,1] <- NA
              res[9,1] <- bs_p_pl$p

              res[1,2] <- data[[2]]$xmin
              res[2,2] <- sd(bs_p_ln$bootstraps$xmin)
              res[3,2] <- NA
              res[4,2] <- NA
              res[5,2] <- data[[2]]$pars[1]
              res[6,2] <- sd(bs_p_ln$bootstraps$pars1)
              res[7,2] <- data[[2]]$pars[2]
              res[8,2] <- sd(bs_p_ln$bootstraps$pars2)
              res[9,2] <- bs_p_ln$p

              res[1,3] <- data[[3]]$xmin
              res[2,3] <- sd(bs_p_ex$bootstraps$xmin)
              res[3,3] <- data[[3]]$pars
              res[4,3] <- sd(bs_p_ex$bootstraps$pars)
              res[5,3] <- NA
              res[6,3] <- NA
              res[7,3] <- NA
              res[8,3] <- NA
              res[9,3] <- bs_p_ex$p
            } else {
              res <- matrix(ncol = 3, nrow = 9)
              colnames(res) <- c("PowerLaw", "LogNormal", "Exponential")
              rownames(res) <- c(
                "Xmin", "SDxmin", "alpha/exponent", "SDalpha/exponent",
                "log(mu)", "SDlog(mu)", "log(sigma)", "SDlog(sigma)",
                "Boot_P_Value"
              )

              res[1,1] <- data[[1]]$xmin
              res[2,1] <- NA
              res[3,1] <- data[[1]]$pars
              res[4,1] <- NA
              res[5,1] <- NA
              res[6,1] <- NA
              res[7,1] <- NA
              res[8,1] <- NA
              res[9,1] <- NA

              res[1,2] <- data[[2]]$xmin
              res[2,2] <- NA
              res[3,2] <- NA
              res[4,2] <- NA
              res[5,2] <- data[[2]]$pars[1]
              res[6,2] <- NA
              res[7,2] <- data[[2]]$pars[2]
              res[8,2] <- NA
              res[9,2] <- NA

              res[1,3] <- data[[3]]$xmin
              res[2,3] <- NA
              res[3,3] <- data[[3]]$pars
              res[4,3] <- NA
              res[5,3] <- NA
              res[6,3] <- NA
              res[7,3] <- NA
              res[8,3] <- NA
              res[9,3] <- NA
            }
            return(res)
          }
        }
      }
    }
  }



### Function 8 #####################################################################################
#' Compare frequency distributions
#'
#' Testing the fit of different phenotype frequency approximations
#' @param data list of three poweRlaw objects
#' @return a matrix with P values
#' @details This function compares parameters of the approximations of the best fitting
#' power-law, log normal, and exponential distribution and outputs the P values of one-sided tests.
#' P values < 0.05 indicate refusal of HO (no difference in fit) for the distribution listed in the row
#' (i.e., P values < 0.05 indicate that the distribution in the row has worse fit than the #'distribution in the column.)
#'
#' @importFrom  poweRlaw compare_distributions
#'
#' @examples
#' \dontrun{
#' a <- compare_pheno_distr(df)
#' }
#' @export
compare_pheno_distr <-
  function(data) {

    # Check list
    if (length(data) != 3) {
      stop("Inputs must consist of 3 poweRlaw objects")
    } else {

      # Check poweRlaw objects 1
      if (class(data[[1]]) != "displ") {
        stop("First element of input must be poweRlaw object *displ* (power-law)")
      } else {

        # Check poweRlaw objects 2
        if (class(data[[2]]) != "dislnorm") {
          stop("Second element of input must be poweRlaw object *dislnorm* (log-normal)")
        } else {

          # Check poweRlaw objects 3
          if (class(data[[3]]) != "disexp") {
            stop("Thir element of input must be poweRlaw object *disexp* (exponential)")
          } else {
            res <- matrix(ncol = 3, nrow = 3)
            colnames(res) <- c("PowerLaw", "LogNormal", "Exponential")
            rownames(res) <- c("PowerLaw", "LogNormal", "Exponential")
            res[1,1] <- NA
            res[1,2] <- round(compare_distributions(data[[2]], data[[1]])$p_one_sided, 4)
            res[1,3] <- round(compare_distributions(data[[3]], data[[1]])$p_one_sided, 4)
            res[2,1] <- round(compare_distributions(data[[1]], data[[2]])$p_one_sided, 4)
            res[2,2] <- NA
            res[2,3] <- round(compare_distributions(data[[3]], data[[2]])$p_one_sided, 4)
            res[3,1] <- round(compare_distributions(data[[1]], data[[3]])$p_one_sided, 4)
            res[3,2] <- round(compare_distributions(data[[2]], data[[3]])$p_one_sided, 4)
            res[3,3] <- NA

            message("CAUTION: P values of one-sided tests (testing the distribution in the row).
                    Values < 0.05 indicate *refusal* of HO (no difference in fit).")
            return(res)
          }
        }
      }
    }
  }


### Function 9 #####################################################################################
#' Plot distributions
#'
#' This function plots parameters of the approximations of the best fitting power-law, log normal, and exponential distribution
#' @param data Parameters of the approximations of the best fitting power-law, log normal, and exponential distribution
#' @param limity limits of the y-axis (default = 10^-4)
#' @param limitx limits of the x-axis (default = 10^3)
#' @return a ggplot object
#'
#' @import ggplot2
#' @import scales
#'
#' @examples
#' \dontrun{
#' a <- plot_pheno_distr(distributions, limity = 10^-4, limitx = 10^3)
#' }
#' @export
plot_pheno_distr <- function(data, limity = 10^-4, limitx = 10^3) {

  # Check list
  if (length(data) != 3) {
    stop("Inputs must consist of 3 poweRlaw objects")
  } else {

    # Check poweRlaw objects 1
    if (class(data[[1]]) != "displ") {
      stop("First element of input must be poweRlaw object *displ* (power-law)")
    } else {

      # Check poweRlaw objects 2
      if (class(data[[2]]) != "dislnorm") {
        stop("Second element of input must be poweRlaw object *dislnorm* (log-normal)")
      } else {

        # Check poweRlaw objects 3
        if (class(data[[3]]) != "disexp") {
          stop("Thir element of input must be poweRlaw object *disexp* (exponential)")
        } else {

          res_pl <- plot(data[[1]])
          line_pl <- lines(data[[1]])
          line_ln <- lines(data[[2]])
          line_ex <- lines(data[[3]])

          ## plot
          p <- ggplot(res_pl, aes(x = x, y = y)) +
            scale_y_log10(
              breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)),
              limits = c(limity, 1)
            ) +
            scale_x_log10(
              breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)),
              limits = c(1, limitx)
            ) +
            geom_line(data = line_pl, aes(x = x, y = y), color = "red", size = 1) +
            geom_line(data = line_ln, aes(x = x, y = y), color = "blue", size = 1, linetype = "dashed") +
            geom_line(data = line_ex, aes(x = x, y = y), color = "orange", size = 1, linetype = "twodash") +
            geom_point(size = 1) +
            xlab("") +
            ylab("") +
            theme_minimal() +
            theme(
              axis.ticks = element_blank(),
              panel.grid.major.x = element_line(size = .2, color = "black"),
              panel.grid.major.y = element_line(size = .2, color = "black"),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank()
            )

          return(p)
        }
      }
    }
  }
}
