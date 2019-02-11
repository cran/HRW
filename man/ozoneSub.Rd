\name{ozoneSub}
\alias{ozoneSub}
\usage{data(ozoneSub)}
\title{Ozone levels in midwest U.S.A.}
\description{This dataset is a subset of the ozone2 dataset in the fields package. It contains the 8-hour average ozone concentration at 147 sites in the midwest region of the U.S.
}
\format{
  A data frame with 147 observations on the following 3 variables:
  \describe{
    \item{\code{longitude}}{observation longitude.}
    \item{\code{latitude}}{observation latitude.}
    \item{\code{ozone}}{ozone level.}
  }
}
\source{
Aerometric Information Retrieval System, the U.S. Environmental Protection Agency air quality data base.
}
\references{
Nychka, D., Furrer, R., Paige, J. and Sain, S. (2017). fields: Tools for spatial data. R package version 9.0. \url{https://www.r-project.org}.
}
\examples{
library(HRW) ; data(ozoneSub)
if (require("mgcv"))
{
   fit.ozone.mgcv.tp <- gam(ozone ~ s(longitude,latitude,bs = "tp"),
                            data = ozoneSub,method = "REML")
   plot(fit.ozone.mgcv.tp,scheme = 2,
        main = "ozone concentration",bty = "l")
   points(ozoneSub$longitude,ozoneSub$latitude)
}
if (require("fields"))
   US(add = TRUE,lwd = 2)
}

