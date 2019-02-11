\name{lidar}
\alias{lidar}
\usage{data(lidar)}
\title{Light detection and ranging}
\description{
The \code{lidar} data frame has 221 pairs from a LIght Detection And Ranging (LIDAR) experiment.
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{range}}{distance traveled before the light is reflected back to its source.}
\item{\code{logratio}}{logarithm of the ratio of received light from two laser sources.}
}
}
\source{Sigrist, M. (Ed.) (1994). \emph{Air Monitoring by Spectroscopic Techniques} (Chemical Analysis Series, vol. 197). New York: Wiley.
}
\examples{
library(HRW) ; data(lidar)
plot(lidar$range,lidar$logratio)
}

