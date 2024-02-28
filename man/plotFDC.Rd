

\name{plotFDC}
\alias{plotFDC}

\title{plotFDC: Generate FDC results figure}
\arguments{
\item{simulation_results}{A list of results produced by power evaluation functions.}

\item{sample_size}{A numerical value indicating which sample size to be fixed.
For example, 10 means when plotting the relationship between FDC and effect size, we fixed the scenario of sample size at 10. Default is 10.}
}

\description{
Plot false discovery cost results.
This function plots false discovery cost results in a 2x1 panel. The illustration
of each plot from left to right:

1: False discovery cost(FDC) by effect size, each line represents cell type.
Sample size to be fixed at 10 if sample_size=10.

2: False discovery cost(FDC) by top effect size, each line represents sample size.
FDC was the average value across cell types.
}

\value{
This function does not return a value. It generates a two-panel plot visualizing
the false discovery cost (FDC) results. The first panel shows the FDC by effect
size for each cell type at a fixed sample size (default is 10). The second panel
illustrates the FDC by the top effect sizes, with each line representing a
different sample size, averaged across cell types.
}
\examples{

 data(quickPowerGSE60424)
### Plot power results
 plotFDC(GSE60424Power,sample_size=10)

}

\author{
Wen Tang <wxt175@case.edu>
Shilin Yu <sy597@georgetown.edu>
}
