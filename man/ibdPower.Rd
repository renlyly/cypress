

\name{ibd_propPower}
\alias{ibd_propPower}
\docType{data}
\title{Power calculation results from IBD data with cell type proportion}
\format{
A S4 object.
   \describe{
  \item{\code{ct_TDR_bio_smry}}{A data frame related to cell type specific
   target TDR.}
  \item{\code{TDR_bio_smry}}{A data frame related to target TDR.}
  \item{\code{ct_PWR_bio_smry}}{A data frame related to cell type specific
   target power.}
  \item{\code{PWR_bio_smry}}{A data frame related to target power.}
  \item{\code{PWR_strata_bio_smry}}{A data frame related to target power by
   gene expression stratification.}
  \item{\code{PWR_strata_ct_bio_smry}}{A data frame related to cell type
   specific target power by gene expression stratification.}
  \item{\code{ct_FDC_bio_smry}}{A data frame related to cell type specific
   target FDC.}
  \item{\code{FDC_bio_smry}}{A data frame related to target FDC.}
  }
}

\usage{
data('quickPowerIBD')
}
\description{
pre-calculated power evaluation results from IBD datasets. The obtained results
can be used to create plots from plot functions.
}
\value{
One S4 object.
}
\examples{
data('quickPowerIBD')

}
\keyword{datasets}
