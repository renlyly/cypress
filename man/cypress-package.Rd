\name{cypress-package}
\alias{cypress}
\alias{cypress-package}

\docType{package}

\title{cypress: cell-type-specific power assessment}
\description{
  The \pkg{cypress}  package is specifically designed to perform comprehensive cell-type-specific power assessment for differential expression using RNA-sequencing experiments. It accepts real Bulk RNAseq data as input for parameter estimation and simulation. The tool provides flexibility by allowing users to customize sample sizes, number of cell types, and log2 fold change values. Additionally, it computes statistical power, true discovery rate (TDR), and false discovery cost (FDC) under different scenarios as results.
}


\details{

cypress is the first statistical tool to evaluate the power in cell-type-specific Differentially Expressed (csDE) genes detection experiments from a prospective way by letting researchers be flexible in tuning sample sizes, effect sizes, csDE genes percentage, total number of genetic features, type I error control, etc.



}

\value{

cypress offers 3 functions for simulation and power evaluation: \code{simFromData} ,\code{simFromParam} and ``quickPower()``. If users have their own bulk RNA-seq count data, they can use the ``simFromData()`` function, otherwise, they can use ``simFromParam()`` estimating parameters from 3 existing studies to perform power evaluation on customized the simulation settings. If users prefer to quickly examine the power evaluation results using the built-in datasets, they can use the ``quickPower()`` function. The output of these 3 functions is a S4 object with a list of power measurements under various experimental settings, such as Statistical Power, TDR, and FDC.

Once users have obtained an S4 object with a list of power measurements using either \code{simFromData}, \code{simFromParam} or \code{quickPower()}, they can use the following functions to generate basic line plots: \code{plotPower()}, \code{plotTDR()} and \code{plotFDC()}.

}


\author{
Shilin Yu <sy597@georgetown.edu>
Guanqun Meng <gxm324@case.edu>
Wen Tang <wxt175@case.edu>
}

\keyword{package}
