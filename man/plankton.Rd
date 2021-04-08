\name{plankton}
\alias{plankton}
\usage{data(plankton)}
\title{Flow cytometric measurements on plankton organisms}
\description{
The \code{plankton} data frame has data on six flow cytometric measurements for 400 plankton organisms categorized into 5 different species. The data are synthetic and were generated to test various machine learning algorithms for plankton species classification. More details are given in the 2001 `Cytometry' article by Boddy, Wilkins and Morris (full reference below).
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{timeFlight}}{time of flight.}
\item{\code{forwScatt}}{forward-scatter.}
\item{\code{sideScatt}}{side-scatter.}
\item{\code{redFluorBlueLight}}{red fluorescence under blue light.}
\item{\code{greenFluorBlueLight}}{green fluorescence under blue light.}
\item{\code{redFluorRedLight}}{red fluorescence under red light.}
\item{\code{species}}{name of the plankton species, which
	 is either Dunaliella, Hemiselmis, Isochrysis,
	 Pavlova or Pyramimonas.}
}
}	   
\source{
Boddy, L., Wilkins, M.F. and Morris, C.W. (2001). Pattern recognition in flow cytometry. \emph{Cytometry}, \bold{44}, 195-209.
}
\examples{
library(HRW) ; data(plankton)
pointCols <- c("red","blue","green","orange","purple")
pairs(plankton[,1:6],col = pointCols[plankton$species],pch = ".")
}
