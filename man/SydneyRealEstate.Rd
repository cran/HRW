\name{SydneyRealEstate}
\alias{SydneyRealEstate}
\usage{data(SydneyRealEstate)}
\title{Sydney real estate}
\description{
The \code{SydneyRealEstate} data frame has data on several variables concerning houses sold in Sydney, Australia, 
during 2001.
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{logSalePrice}}{natural logarithm of sale price in Australian dollars.} 
\item{\code{lotSize}}{lot size in square meters but with some imputation.} 
\item{\code{longitude}}{degrees longitude of location of house.}  
\item{\code{latitude}}{degrees latitude of location of house.}  
\item{\code{saleDate}}{sale date in dd/mm/yyy format.} 
\item{\code{saleQtr}}{financial quarter in which sale took place.} 
\item{\code{infRate}}{inflation rate measure as a percentage.} 
\item{\code{postCode}}{four-digit post code of the suburb in which the house located.}  
\item{\code{crimeDensity}}{crime density measure for the suburb in which the house is located.} 
\item{\code{crimeRate}}{crime rate measure for the suburb in which the house is located.}  
\item{\code{income}}{average weekly income of the suburb in which the house is located.} 
\item{\code{distToBusStop}}{distance from house to the nearest bus stop (kilometers).}  
\item{\code{distToCoastline}}{distance from house to the nearest coastline location (kilometers).} 
\item{\code{distToNatPark}}{distance from house to the nearest national park (kilometers).} 
\item{\code{distToPark}}{distance from house to the nearest park (kilometers).} 
\item{\code{distToRailLine}}{distance from house to the nearest railway line (kilometers).}  
\item{\code{distToRailStation}}{distance from house to the nearest railway station (kilometers).} 
\item{\code{distToHighway}}{distance from house to the nearest highway (kilometers).} 
\item{\code{distToFreeway}}{distance from house to the nearest freeway (kilometers).} 
\item{\code{distToTunnel}}{distance from house to the Sydney Harbour Tunnel (kilometers).} 
\item{\code{distToMainRoad}}{distance from house to the nearest main road (kilometers).} 
\item{\code{distToSealedRoad}}{distance from house to the nearest sealed road (kilometers).}  
\item{\code{distToUnsealedRoad}}{distance from house to the nearest sealed road (kilometers).} 
\item{\code{airNoise}}{aircraft noise exposure measure.} 
\item{\code{foreignerRatio}}{proportion of foreigners in the suburb in which the house is located.} 
\item{\code{distToGPO}}{distance from the house to the General Post Office in Sydney's central business district (kilometers).} 
\item{\code{NO}}{nitrous oxide level recorded at the air pollution monitoring station nearest to the house.}  
\item{\code{NO2}}{nitrogen dioxide level recorded at the air pollution monitoring station nearest to the house.}  
\item{\code{ozone}}{ozone level recorded at the air pollution monitoring station nearest to the house.}   
\item{\code{neph}}{nephelometer suspended matter measurement recorded at the air pollution monitoring station nearest to the house.}   
\item{\code{PM10}}{particulate matter with a diameter of under 10 micrometers leve recorded at the air pollution monitoring station nearest to the house.} 
\item{\code{SO2}}{sulphur dioxide level recorded at the air pollution monitoring station nearest to the house.}  
\item{\code{distToAmbulance}}{distance from house to the nearest ambulance station (kilometers).}
\item{\code{distToFactory}}{distance from house to the nearest factory (kilometers).}
\item{\code{distToFerry}}{distance from house to the nearest ferry wharf (kilometers).}  
\item{\code{distToHospital}}{distance from house to the nearest hospital (kilometers).} 
\item{\code{distToMedical}}{distance from house to the nearest medical services (kilometers).}  
\item{\code{distToSchool}}{distance from house to the nearest school (kilometers).}  
\item{\code{distToUniversity}}{distance from house to the nearest university (kilometers).} 
}       
}
\references{Chernih, A. and Sherris, M. (2004). Geoadditive hedonic pricing models. Unpublished manuscript. University of New South Wales, Australia.}
\examples{
library(HRW) ; data(SydneyRealEstate)
\dontrun{
for (colNum in setdiff((2:39),c(5,8)))
{
   plot(jitter(SydneyRealEstate[,colNum]),SydneyRealEstate$logSalePrice,pch = ".",
        xlab = names(SydneyRealEstate)[colNum],ylab = "log(sale price)",col = "blue")
   readline("Hit Enter to continue.\n")
}
}
}

