

\name{plotTDR}
\alias{plotTDR}
\title{plotTDR:Generate TDR results figure}
\arguments{
\item{simulation_results}{A list of results produced by power evaluation functions.}

\item{effect.size}{A numerical value indicating which effect size is to be fixed.
For example, 1 means when plotting the relationship between TDR and top rank genes for cell types or sample size, we fixed the scenario of log2 fold change at 1. The default is 1.}

\item{sample_size}{A numerical value indicating which sample size to be fixed.
For example, 10 means when plotting the relationship between TDR and top rank genes for cell types or effect size, we fixed the scenario of sample size at 10. The default is 10.}
}

\description{
This function plots all true discovery rate measurements in a 2x2 panel.The illustration of each plot is as follow:

1: True discovery rate(TDR) by top rank genes, each line represent cell type.
Sample size to be fixed at 10 and effect size to be fixed at 1 if sample_size=10 and effect.size=1.

2: True discovery rate(TDR) by top rank genes, each line represent effect size.
TDR was average value across cell types and sample size to be fixed at 10 if sample_size=10.

3: True discovery rate(TDR) by top rank genes, each line represent sample size.
TDR was average value across cell types and effect size to be fixed at 1 if effect.size=1.

4: True discovery rate(TDR) by effect size, each line represent sample size.
TDR was calculated under the scenario of top rank gene equals to 350.


}
\value{
This function creates a 2x2 panel plot showcasing various true discovery rate
(TDR) measurements but does not return any values for further programmatic use.
Each panel displays TDR analyses based on top rank genes, with lines
representing different cell types, effect sizes, or sample sizes under specific
conditions.
}
\examples{

 data(quickPowerGSE60424)
### plot power results
 plotTDR(GSE60424Power,effect.size=1,sample_size=10)

}

\author{
Wen Tang <wxt175@case.edu>
Shilin Yu <sy597@georgetown.edu>
}
