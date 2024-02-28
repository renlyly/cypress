
\docType{data}
\name{ibd_prop_param}
\alias{ibd_prop_param}

\alias{quickParaIBD}
\alias{ibd_param}

\title{parameter esitmate results from IBD data with cell type proportion}
\format{
Parameter esitmate results from IBD for simFromPara function
  \describe{
  \item{\code{health_lm_mean}}{A numeric vector of the log mean parameter
  estimates for healthy samples.}
  \item{\code{health_lm_mean_d}}{A matrix containing the variance-covariance
  estimates for the healthy sample means.}
  \item{\code{lod_m}}{A numeric vector of log-dispersion  parameter estimates.}
  \item{\code{lod_d}}{A numeric matrix includes the variance-covariance
  estimates for the dispersion parameters.}
  \item{\code{health_alpha}}{A numeric vector for the alpha parameter estimates
  associated with healthy samples.}
  \item{\code{health_alpha}}{A numeric vector for the alpha parameter estimates
  associated with case samples.}
  }
}

\usage{
data('quickParaIBD')
}
\description{
 An S4 object of that stores parameter estimates associated with the IBD
 dataset. This object contains a variety of numerical vectors and matrices
 representing different statistical parameters used in the simulation.
}
\value{
One S4 object.
}
\examples{
data('quickParaIBD')

}
\keyword{datasets}
