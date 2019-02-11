\name{summMCMC}
\alias{summMCMC}
\title{Summarizes Markov chain Monte Carlo (MCMC) samples both graphically and numerically}
\description{Given a set of MCMC for possibly several parameters the following summaries are produced: 1. a trace (time series) plot of each MCMC, 2. a plot of each MCMC sample against the 1-lagged sample, 3. estimated autocorrelatio function (acf), 4. Brooks-Gelman-Rubin (BGR) diagnostic plot in cases where multiple chains are inputted, 5. kernel ddensity estimate of the posterior density function based on the MCMC sample, 6. numerical summary consisting of the MCMC-based estimates of posterior means and credible sets for each parameter.
}
\usage{
summMCMC(xList,EPSfileName,PDFfileName,plotInd = 1,parNames,columnHeadings,
         colourVersion = TRUE,credLevel = 0.95,columnHeadCex = NULL,paletteNum = 1,
         numerSummCex = NULL,BGRsttPos = 10,BGRyRange = c(0.95,1.25),BGRtickPos = 1.2,
         BGRlogTransf,BGRlogitTransf,KDExlim,KDEvertLine = TRUE,KDEvertLineCol = "black",
         addTruthToKDE = NULL)
}
\arguments{
\item{xList}{list of matrices, where each matrix corresponds to a different chain, and the columns of each matrix correspond to different parameters. The matrices each have dimension "numMCMC" by "numParms"; where "numMCMC" is the size of the MCMC sample and "numParms" is the number of parameters being summarized.}
\item{EPSfileName}{filename if the summary is to be saved as a (encapsulated) Postscript file. If this argument and 'PDFfileName' are both not specified then the summary is printed to the screen.}
\item{PDFfileName}{filename if the summary is to be saved as a PDF file. If this argument and 'EPSfileName' are both not specified then the summary is printed to the screen.}
\item{plotInd}{if "numChains" exceeds 1 then this indicates which chain is summarized in the non-BGR panels. The BGR panels are Brooks-Gelman-Rubin diagnostic plots are use all chains. The default value is 1.}
\item{parNames}{list containing a vector of character strings for the parameter names. The maximum length of the vector is 3.}
\item{columnHeadings}{vector containing column headings. The default is: c("parameter","trace","lag 1","acf","BGR","density","summary").}
\item{colourVersion}{logical flag indicating if summary should be in colour. The default is TRUE.}
\item{credLevel}{number between 0 and 1 specifying the credible set level. The default is 0.95.}
\item{columnHeadCex}{positive number specifying character expansion factor for the column headings.}
\item{paletteNum}{either 1 or 2 to specify the colour palette to be used.}
\item{numerSummCex}{positive number specifying character expansion factor for the numerical summary (last column).}
\item{BGRsttPos}{starting position for the Brooks-Gelman-Rubin plots. The default value is 10.}
\item{BGRyRange}{vertical axis limits for the Brooks-Gelman-Rubin plots. The default value is  c(0.95,1.25).}
\item{BGRtickPos}{position of tick mark on vertical axis for the Brooks-Gelman-Rubin plots. The default value is 1.2.}
\item{BGRlogTransf}{vector containing indices of those parameters for which the Brooks-Gelman-Rubin plots should be done on a logarithmic scale.}
\item{BGRlogitTransf}{vector containing indices of those parameters for which the Brooks-Gelman-Rubin plots should be done on a logit scale.}
\item{KDExlim}{list of vectors of length 2 specifying the horizontal axis limits for the kernel density estimates.}
\item{KDEvertLine}{logical flag indicating if a vertical line at zero should be added to the kernel density estimates. The default value is TRUE.}
\item{KDEvertLineCol}{Colour of the vertical line at zero for kernel density estimates. The default value is "black".}
\item{addTruthToKDE}{Vector indicating `true' values of parameters. The default value is NULL. If 'addTruthToKDE' is non-NULL then dashed vertical lines corresponding to true values (if known from simulation) are added.}
}
\author{Matt Wand
\email{matt.wand@uts.edu.au} 
}
\examples{
library(HRW)
xListSingleChain <- list(cbind(rnorm(100),rnorm(100),rnorm(100),rnorm(100)))
summMCMC(xListSingleChain,parNames = list("par1","par2","par3","par4"))
xListMultipleChains <- list(chain1 = cbind(rnorm(100),rnorm(100),rnorm(100),rnorm(100)),
                            chain2 = cbind(rnorm(100),rnorm(100),rnorm(100),rnorm(100)))
summMCMC(xListMultipleChains,parNames = list("par1","par2","par3","par4"))
}
