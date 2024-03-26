\name{simFromData}
\alias{simFromData}

\title{Power calculation results on input data}

\usage{
simFromData(INPUTdata = NULL, CT_index = NULL, CT_unk = FALSE,
            n_sim = 3, n_gene = 30000, DE_pct = 0.05,
            ss_group_set = c(10, 20, 50, 100),
            lfc_set = c(0, 0.5, 1, 1.5, 2),
            lfc_target = 0.5, fdr_thred = 0.1, BPPARAM=bpparam())
}
\arguments{
\item{INPUTdata}{The input SE (SummarizedExperiment) object should contain a
count matrix, study design, and an optional cell type proportion matrix. The
study design should have a column named 'disease', where the control by 1 and
the case is indicated by 2. If provided, the cell type proportion matrix should
sum to 1 for each sample. The cell type proportion matrix is optional,
the CT_unk should be True if the user did not provide this matrix}
  \item{CT_index}{Column index for cell types proportion matrix in Coldata, the
  input can be a single number (>3) when the \code{CT_unk} is T.}
  \item{CT_unk}{Logical flag indicating whether unknown cell types are present.
  Defaults to \code{FALSE}.}
  \item{n_sim}{Total number of iterations users wish to conduct. Default to 3}
  \item{n_gene}{Total number of genetic features users with to conduct. Default
  to 30000. Must be greater than or equal to 1000.}
  \item{DE_pct}{Percentage of DEG on each cell type. Default to 0.05.}
  \item{ss_group_set}{Sample sizes per group users wish to simulate. The length
  should be less than or equal to 5. Default to 10,20,50.}
  \item{lfc_set}{Effect sizes users wish to simulate. The length should be less
  than or equal to 5. Default to 0,0.5,1,1.5.}
  \item{lfc_target}{Target effect size, should be greater than or equal to 0.
  The absolute LFC lower than this values will be treated as None-DEGs.
  Default to 0.5}
  \item{fdr_thred}{Adjusted p value threshold, should not below 0 or exceed 1.
  default to 0.1}
   \item{BPPARAM}{ An instance of \code{BiocParallelParam} class, e.g. MulticoreParam, SnowParam, SerialParam, to facilitate parallel computing. If using Unix, MulticoreParam is recommended. Customized options within \code{BiocParallelParam} class is allowed. If not specified, the default back-end is retrieve.
}
}

\value{
  \item{\code{ct_TDR_bio_smry}}{A data frame related to cell type specific target TDR.}
  \item{\code{TDR_bio_smry}}{A data frame related to target TDR.}
  \item{\code{ct_PWR_bio_smry}}{A data frame related to cell type specific target power.}
  \item{\code{PWR_bio_smry}}{A data frame related to target power.}
  \item{\code{PWR_strata_bio_smry}}{A data frame related to target power by gene expression stratification.}
  \item{\code{PWR_strata_ct_bio_smry}}{A data frame related to cell type specific target power by gene expression stratification.}
  \item{\code{ct_FDC_bio_smry}}{A data frame related to cell type specific  target FDC.}
  \item{\code{FDC_bio_smry}}{A data frame related to target FDC.}
}
\description{
This function conducts simulations with various user-defined study design
parameters, Users will need to provide SE object bulk data for parameter
estimation purposes.
}
\details{
  One \code{SummarizedExperiment} object containing the following elements:
  \describe{
    \item{\code{counts}}{A gene expression value dataset}
    \item{\code{colData}}{Sample meta-data. The first column is the group status
    (i.e. case/ctrl) named as 'disease', and the second column is the subject ID.
    The remaining are the cell type proportions of all samples. The user could
    also input the Column index for cell types proportion matrix in Coldata.
    Example: CT_index= 3:8}
  }

}
\examples{
\donttest{
data(ASD_prop_se)
result <- simFromData(INPUTdata = ASD_prop, CT_index = (seq_len(6) + 2), CT_unk = FALSE,
                      n_sim = 2,n_gene = 1000, DE_pct = 0.05,
                      ss_group_set = c(8,10),
                      lfc_set = c(1, 1.5))
                      }
}
\author{
Shilin Yu <sy597@georgetown.edu>
Guanqun Meng <gxm324@case.edu>
}
