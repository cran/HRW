\name{ragweed}
\alias{ragweed}
\usage{data(ragweed)}
\title{Ragweed pollen in Kalamazoo, U.S.A.}
\description{
The \code{ragweed} data frame has data on ragweed levels and meteorological variables for 334 days in Kalamazoo, Michigan, U.S.A.
}
\format{
  This data frame contains the following columns:
  \describe{
       \item{\code{pollenCount}}{ragweed pollen count (grains per cubic metre).}
       \item{\code{year}}{one of 1991, 1992, 1993 or 1994.}
       \item{\code{dayInSeason}}{day number in the current ragweed pollen season.}
       \item{\code{temperature}}{temperature for the corresponding day (degrees Fahrenheit).}
       \item{\code{temperatureResidual}}{residual from fitting a 5 effective degrees of freedom smoothing splines to temperature versus day number for each annual ragweed pollen season (degrees Fahrenheit).}
       \item{\code{rain}}{indicator of significant rain on the corresponding day:\cr
                   1=at least 3 hours of steady or brief but intense rain,\cr
                   0=otherwise.}
       \item{\code{windSpeed}}{wind speed for the corresponding day (knots).}

  }
}
\source{
Stark, P. C., Ryan, L. M., McDonald, J. L. and Burge, H. A. (1997). Using meteorologic data to model and predict daily ragweed pollen levels. \emph{Aerobiologia}, \bold{13}, 177-184.
}
\references{
Ruppert, D., Wand, M.P. and Carroll, R.J. (2003). \emph{Semiparametric Regression} Cambridge University Press. \url{http://stat.tamu.edu/~carroll/semiregbook/}
}
\examples{
library(HRW) ; data(ragweed)
pairs(ragweed,pch = ".")
}

