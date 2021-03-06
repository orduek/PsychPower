% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Frequency_Distributions.R
\name{describe_pheno_distr}
\alias{describe_pheno_distr}
\title{Parameters of frequency distributions}
\usage{
describe_pheno_distr(data, bootStrap = T, nBoots = 1000, nCores = 1)
}
\arguments{
\item{data}{A list of three "poweRlaw" objects.}

\item{bootStrap}{whether to run the bootstrap or not (default = TRUE)}

\item{nBoots}{Number of bootstraps (default = 1000)}

\item{nCores}{Number of cores to use in computing results (default = 1).}
}
\value{
a matrix with the parameters of the approximations of the best fitting power-law, log normal, and exponential distribution
}
\description{
This function gives back the parameters of the best fitting distributions estimated with
\code{pheno_distributions}.
}
\examples{
\dontrun{
desc_freq <- describe_pheno_distr(freq_distribution, bootStrap = T, nBoots = 1000, nCores = 1)
}
}
