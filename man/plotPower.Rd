
\name{plotPower}
\alias{plotPower}
\title{plotPower: Generate statistical power results figure}
\arguments{
\item{simulation_results}{A list of results produced by power evaluation functions.}

\item{effect.size}{A numerical value indicating which effect size is to be fixed.
For example, 1 means when plotting the relationship between power and strata, we fixed the scenario of log fold change at 1. The default is 1.}

\item{sample_size}{A numerical value indicating which sample size to be fixed.
For example, 10 means when plotting the relationship between power and strata, we fixed the scenario of sample size at 10. The default is 10.}
}

\description{
This function plots all statistical power measurements in a 2x3 panel. The illustration of each plot from left to right and from up to bottom is as follows:

1: Statistical power by effect size, each line represents sample size. Statistical power was the average value across cell types.

2: Statistical power by effect size, each line represents cell type. Sample size is fixed at 10 if sample_size=10.

3: Statistical power by sample size, each line represents cell type. Effect size is fixed at 1 if effect.size=1.

4: Statistical power by strata, each line represents cell type. Sample size is fixed at 10 and effect size is fixed at 1 if sample_size=10 and effect.size=1.

5: Statistical power by strata, each line represents sample size. Statistical
power was the average value across cell types and effect size is fixed at 1 if effect.size=1.

6: Statistical power by strata, each line represents effect size. Statistical
power was the average value across cell types and sample size is fixed at 10 if sample_size=10.
}
\value{
This function generates a 2x3 panel plot visualizing various statistical power
measurements, but does not return a programmable value. Each panel displays
power metrics under different conditions such as effect size, sample size, and
stratification, with lines representing either sample size, cell type, or
effect size.
}

\examples{

data(quickPowerGSE60424)
### Plot power results
plotPower(GSE60424Power,effect.size=1,sample_size=10)


}
\author{
Wen Tang <wxt175@case.edu>
Shilin Yu <sy597@georgetown.edu>
}
