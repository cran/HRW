\name{scallop}
\alias{scallop}
\usage{data(scallop)}
\title{Scallop abundance off Long Island, U.S.A.}
\description{
The \code{scallop} data frame has 148 triplets concerning scallop abundance; based on a 1990 survey cruise in the Atlantic continental shelf off Long Island, New York, U.S.A.}
\format{
This data frame contains the following columns:
\describe{
\item{\code{latitude}}{degrees latitude (north of the Equator).}
\item{\code{longitude}}{degrees longitude (west of Greenwich).}
\item{\code{totalCatch}}{total size of scallop catch at location specified by 'latitude' and 'longitude'.}
}
}
\source{Ecker, M.D. and Heltshe, J.F. (1994).
Geostatistical estimates of scallop abundance. In \emph{Case Studies in Biometry.} Lange, N., Ryan, L., Billard, L., Brillinger, D., Conquest, L. and Greenhouse, J. (eds.) New York: John Wiley & Sons, 107-124.
}
\references{
Ruppert, D., Wand, M.P. and Carroll, R.J. (2003). \emph{Semiparametric Regression.} Cambridge University Press.
}
\examples{
library(HRW) ; data(scallop)
pairs(scallop)
}
