\name{quickPower}
\alias{quickPower}
\title{Obtain pre-calculated results from three available datasets}
\usage{
quickPower(data = "IAD")
}
\arguments{
  \item{data}{A character string specifying the dataset to be retrieved.
  Options are 'IAD', 'IBD', and 'ASD'.}
}
\value{
  \item{\code{ct_TDR_bio_smry}}{ A data frame related to cell type specific
   target TDR.}
  \item{\code{TDR_bio_smry}}{ A data frame related to target TDR.}
  \item{\code{ct_PWR_bio_smry}}{ A data frame related to cell type specific
   target power.}
  \item{\code{PWR_bio_smry}}{A data frame related to target power.}
  \item{\code{PWR_strata_bio_smry}}{A data frame related to target power by
   gene expression stratification.}
  \item{\code{PWR_strata_ct_bio_smry}}{ A data frame related to cell type
   specific target power by gene expression stratification.}
  \item{\code{ct_FDC_bio_smry}}{ A data frame related to cell type specific
  target FDC.}
  \item{\code{FDC_bio_smry}}{ A data frame related to target FDC.}
}
\description{
This function quickly outputs pre-calculated power evaluation results from three
datasets: (IAD, IBD, and ASD). The obtained results can be used to create
plots from plot functions.}

\details{
{IAD: }{whole transcriptome signatures of 6 immune cell subsets. The
patients were drawn from subjects with a range of immune-related diseases.}\cr
{IBD: }{Inflammatory Bowel Disease}\cr
{ASD: }{Autism Spectrum Disorder.}\cr
}


\examples{
# library(cypress)
Quick_power <- quickPower(data = "IAD")
}
\author{
Shilin Yu <sy597@georgetown.edu>
Guanqun Meng <gxm324@case.edu>
}
