\name{WarsawApts}
\alias{WarsawApts}
\usage{data(WarsawApts)}
\title{Apartment prices in Warsaw, Poland}
\description{
'WarsawApts' is a subset of the data set 'apartments' in the R package PBImisc. This dataset contains the prices of the apartments which were sold in Warsaw, Poland, during the calendar years 2007 to 2009.}
\format{
  A data frame with 409 observations on the following 6 variables:
  \describe{
    \item{\code{surface}}{area of the apartment in square meters.}
    \item{\code{district}}{a factor corresponding to the district of Warsaw with levels \code{Mokotow}, \code{Srodmiescie}, \code{Wola} and \code{Zoliborz}.}
    \item{\code{n.rooms}}{number of rooms in the apartment.}
    \item{\code{floor}}{floor on which the apartment is located.}
    \item{\code{construction.date}}{year that the apartment was constructed.}
    \item{\code{areaPerMzloty}}{area in square meters per million zloty.}
  }
}
\source{
The Polish real estate web-site \url{www.oferty.net}.
}
\references{
Biecek, P. (2016). PBImisc: A set of datasets in My Classes or in the Book 'Modele Liniowe i Mieszane w R, Wraz z Przykladami w Analizie Danych' 1.0. \url{http://www.biecek.pl/R/}.
}
\examples{
library(HRW) ; data(WarsawApts)
x <- WarsawApts$construction.date
y <- WarsawApts$areaPerMzloty
plot(x,y,bty = "l",col = "dodgerblue")
if (require("mgcv"))
{
   fitGAMcr <- gam(y ~ s(x,bs = "cr",k = 30))
   xg <- seq(min(x),max(x),length = 1001)
   fHatgGAMcr <- predict(fitGAMcr,newdata = data.frame(x = xg))
   lines(xg,fHatgGAMcr,col = "darkgreen")
}
}

