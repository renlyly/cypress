
\docType{data}
\name{GSE60424_param}
\alias{GSE60424_param}
\title{parameter esitmate results from GSE60424}
\format{
 Parameter esitmate results from GSE60424 for simFromPara function
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
data('quickParaGSE60424')
}
\description{
 An S4 object of that stores parameter estimates associated with the GSE60424
 dataset. This object contains a variety of numerical vectors and matrices
 representing different statistical parameters used in the simulation.
 The patients were drawn from subjects with a range of immune-related diseases.
}

\value{
One S4 object.
}

\examples{
data('quickParaGSE60424')

}
\keyword{datasets}
