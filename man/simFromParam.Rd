\name{simFromParam}
\alias{simFromParam}
\title{Power calculation results on pre-calculated parameters}

\usage{
simFromParam(n_sim = 3, n_gene = 30000, DE_pct = 0.05,
             ss_group_set = c(10, 20, 50, 100),
             lfc_set = c(0, 0.5, 1, 1.5, 2),
             sim_param = "IAD",
             lfc_target = 0.5, fdr_thred = 0.1, BPPARAM=bpparam())
}
\arguments{
  \item{n_sim}{Total number of iterations users wish to conduct. Default to 3}
  \item{n_gene}{Total number of genetic features users with to conduct. Default
  to 30000. Must be greater than or equal to 1000.}
  \item{DE_pct}{Percentage of DEG on each cell type. Default to 0.05.}
  \item{ss_group_set}{Sample sizes per group users wish to simulate. The length
  should be less than or equal to 5. Default to 10,20,50.}
  \item{lfc_set}{effect sizes users wish to simulate. The length should be less
  than or equal to 5. Default to 0,0.5,1,1.5.}
  \item{sim_param}{ Users specify which embedded simulation parameters they wish
  to use. By default set to 'GSE60424_param', which is a cell line specific
  data. Other options include 'IBD' data, and 'ASD' data}
  \item{lfc_target}{Target effect size, should be greater than or equal to 0.
    The absolute LFC lower than this value will be treated as None-DEGs.
    Default to 0.5}
  \item{fdr_thred}{Adjusted p value threshold, should not be below 0 or exceed 1.
    default to 0.1}
  \item{BPPARAM}{ An instance of \code{BiocParallelParam} class, e.g. MulticoreParam, SnowParam, SerialParam, to facilitate parallel computing. If using Unix, MulticoreParam is recommended. Customized options within \code{BiocParallelParam} class is allowed. If not specified, the default back-end is retrieve.
}
}

\value{
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
\description{
This function conducts simulations with various user-defined study design
parameters, including but not limited to sample size, and log fold change.)
}
\details{
  \itemize{
    \item \code{GSE60424:} Whole transcriptome signatures of 6 immune cell subsets. The patients were drawn from subjects with a range of immune-related diseases.
    \item \code{IBD:} Inflammatory Bowel Disease
    \item \code{ASD:} Autism Spectrum Disorder.
  }
}
\examples{
data(quickParaGSE60424)
result <- simFromParam(sim_param="IAD",n_sim = 2,DE_pct = 0.05,n_gene = 1000,
                         ss_group_set = c(8, 10),
                         lfc_set = c(1, 1.5),
                         lfc_target = 0.5, fdr_thred = 0.1)
}
\author{
Shilin Yu <sy597@georgetown.edu>
Guanqun Meng <gxm324@case.edu>
}
