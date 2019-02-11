\name{SydneyRealEstateBdry}
\alias{SydneyRealEstateBdry}
\usage{data(SydneyRealEstateBdry)}
\title{Polygonal boundary of Sydney}
\description{A two-column data frame containing the (longitude,latitude) pairs for the vertices of a 202-sided polygon. The polygon was created manually using the HRW package function createBdry(). The polygon tightly encompasses the majority of the (longitude,latitude) data of the HRW package `SydneyRealEstate` data frame and approximately corresponds to the residential parts of Sydney, Australia. 
}
\format{A data frame with 203 observations on the following 2 variables (note that the first vertex is repeated at the end of the data frame):
  \describe{
    \item{\code{longitude}}{longitudinal position of a vertex.}
    \item{\code{latitude}}{latitudinal position of the same vertex.}
  }
}
\seealso{SydneyRealEstate}
\examples{
library(HRW) ; data(SydneyRealEstateBdry) ; data(SydneyRealEstate)
plot(SydneyRealEstate$longitude,SydneyRealEstate$latitude,cex = 0.1)
lines(SydneyRealEstateBdry,lwd = 5,col = "red")
}

