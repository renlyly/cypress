
\docType{data}
\name{ASD_prop}
\alias{ASD_prop}

\title{SimFromData example raw input data}
\format{
SE object.
}
\usage{
data(ASD_prop_se)
}
\description{
  \code{ASD_prop} contains a \code{SummarizedExperiment} object containing the
   following elements:
  \describe{
    \item{\code{counts}}{A gene expression value dataset From ASD, in the form
    of raw read counts, 29674 genes by 48 sample, with 24 cases and 24 controls}
    \item{\code{colData}}{Sample meta-data. The first column is the group status
    (i.e. case/ctrl), the second column is the subject ID. The remaining are the
    cell type proportions of all samples.}
  }

}
\value{
One S4 object.
}

\examples{
data(ASD_prop_se)
}
\keyword{datasets}
