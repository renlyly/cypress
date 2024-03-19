\name{cypress-class}

\alias{cypress-class}
\alias{AllClasses}
\alias{est_out-class}
\alias{cypress_out-class}


\title{
Classes for cypress Package
}
\description{
The \code{cypress_out} and \code{est_out} classes are custom S4 classes in the \pkg{cypress} package. both classes are designed as a comprehensive container for various types of analysis results.
}


\section{cypress_out}{
  \subsection{Description}{
The \code{cypress_out} class is a S4 class in the \pkg{cypress} package. This class is customized to better present results and use for cypress plot functions.
  }
  \subsection{Slots}{
    \describe{
      \item{\code{ct_TDR_bio_smry}}{Cell type specific target TDR}
      \item{\code{TDR_bio_smry}}{Target TDR}
      \item{\code{ct_PWR_bio_smry}}{Cell type specific target power}
      \item{\code{PWR_bio_smry}}{Target power}
      \item{\code{PWR_strata_bio_smry}}{Target power by gene expression stratification}
      \item{\code{PWR_strata_ct_bio_smry}}{Cell type specific target power by gene expression stratification}
      \item{\code{ct_FDC_bio_smry}}{Cell type specific target FDC}
      \item{\code{FDC_bio_smry}}{Target FDC.}
    }
  }
}

\section{est_out}{
  \subsection{Description}{
    The \code{est_out} class is designed to output the parameter estimated results, providing a structured representation of results.
  }
  \subsection{Slots}{
    \describe{
      \item{\code{health_alpha}}{Control group proportion simulation parameter.}
      \item{\code{case_alpha}}{Case group proportion simulation parameter.}
      \item{\code{health_lmean_m}}{Mean of genetic distribution mean for each cell.}
      \item{\code{health_lmean_d}}{Var/cov of genetic distribution mean among cell types}
      \item{\code{lod_m}}{Mean of genetic distribution dispersion for each cell.}
      \item{\code{lod_d}}{Var/cov of genetic distribution dispersion among cell types.}
      \item{\code{sample_CT_prop}}{Matrix of sample Cell type proportions.}
      \item{\code{genename}}{Gene Name.}
      \item{\code{samplename}}{Sample Name.}
      \item{\code{CTname}}{Cell type names.}
      \item{\code{dimensions_Z_hat_ary}}{dimensions for the Z hat array.}
    }
  }
}

\value{
  Class for cypress.
}

\author{
Shilin Yu <sy597@georgetown.edu>
}


\examples{
data(quickParaGSE60424)
}
