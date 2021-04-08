\name{carAuction}
\usage{data(carAuction)}
\alias{carAuction}
\title{Cars purchased at auction}
\description{
The \code{carAuction} data frame has data on several variables concerning cars purchased at automobile auctions by automobile dealerships in the United States of America. The origin of these data is a classification competition titled ``Don't Get Kicked!'' that ran on the `kaggle' platform (\url{https://www.kaggle.com}) during 2011-2012. Many of the variables in this data frame have been derived from those in the original data set from \url{https://www.kaggle.com}.
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{RefId}}{unique number assigned to each vehicles.}              
\item{\code{IsBadBuy}}{indicator that the vehicle purchased at auction by an automobile dealership has serious problems that hinder or prevent it being sold - a "bad buy":\cr
                   1 = the vehicle is a "bad buy"\cr
                   0 = the vehicle is a "good buy".\cr
All other indicator variables are defined in this way.}           
\item{\code{purchIn2010}}{indicator that vehicle was purchased in 2010.}        
\item{\code{aucEqAdesa}}{indicator that the auction provider at which the vehicle was purchased was Adesa.}         
\item{\code{aucEqManheim}}{indicator that the auction provider at which the vehicle was purchased was Manheim.}
\item{\code{vehYearEq03}}{indicator that the manufacturer's year of the vehicle is 2003.}        
\item{\code{vehYearEq04}}{indicator that the manufacturer's year of the vehicle is 2004.}        
\item{\code{vehYearEq05}}{indicator that the manufacturer's year of the vehicle is 2005.}        
\item{\code{vehYearEq06}}{indicator that the manufacturer's year of the vehicle is 2006.}        
\item{\code{vehYearEq07}}{indicator that the manufacturer's year of the vehicle is 2007.}        
\item{\code{ageAtSale}}{age of the vehicle in years when sold.}          
\item{\code{makeEqChevrolet}}{indicator that the vehicle's manufacturer is Chevrolet.}    
\item{\code{makeEqFord}}{indicator that the vehicle's manufacturer is Ford.}         
\item{\code{makeEqDodge}}{indicator that the vehicle's manufacturer is Dodge.}        
\item{\code{makeEqChrysler}}{indicator that the vehicle's manufacturer is Chrysler.}     
\item{\code{trimEqBas}}{indicator that the trim level of the vehicle is 'Bas'.}          
\item{\code{trimEqLS}}{indicator that the trim level of the vehicle is 'LS'.}           
\item{\code{trimEqSE}}{indicator that the trim level of the vehicle is 'SE'.}          
\item{\code{subModelEq4DSEDAN}}{indicator that the submodel of the vehicle is '4DSedan'.}  
\item{\code{subModelEq4DSEDANLS}}{indicator that the submodel of the vehicle is '4DSedanLS'.}
\item{\code{subModelEq4DSEDANSE}}{indicator that the submodel of the vehicle is '4DSedanSE'.}
\item{\code{colourEqSilver}}{indicator that the vehicle color is silver.}     
\item{\code{colourEqWhite}}{indicator that the vehicle color is white.}      
\item{\code{colourEqBlue}}{indicator that the vehicle color is blue.}       
\item{\code{colourEqGrey}}{indicator that the vehicle color is grey.}       
\item{\code{colourEqBlack}}{indicator that the vehicle color is black.}      
\item{\code{colourEqRed}}{indicator that the vehicle color is red.}        
\item{\code{colourEqGold}}{indicator that the vehicle color is gold.}       
\item{\code{colourEqOrange}}{indicator that the vehicle color is orange.}     
\item{\code{transEqManual}}{indicator that the vehicle has manual transmission.}      
\item{\code{wheelEqAlloy}}{indicator that the vehicle has alloy wheels.}       
\item{\code{wheelEqCovers}}{indicator that the vehicle has covered wheels.}      
\item{\code{odomRead}}{the vehicle's odometer reading in miles.}           
\item{\code{AmericanMade}}{indicator that the vehicle was manufactured in the United States of America.}       
\item{\code{otherAsianMade}}{indicator that the vehicle was manuctured in an Asian nation other than Japan or South Korea.}     
\item{\code{sizeEqTruck}}{indicator that the size category of the vehicle is 'truck'.}        
\item{\code{sizeEqMedium}}{indicator that the size category of the vehicle is 'medium'.}      
\item{\code{sizeEqSUV}}{indicator that the size category of the vehicle is 'SUV'.}          
\item{\code{sizeEqCompact}}{indicator that the size category of the vehicle is 'compact'.}      
\item{\code{sizeEqVan}}{indicator that the size category of the vehicle is 'van'.}          
\item{\code{price}}{acquisition price for this vehicle in average condition at time of purchase in U.S. dollars.}              
\item{\code{purchInTexas}}{indicator that the vehicle was purchased in Texas.}        
\item{\code{purchInFlorida}}{indicator that the vehicle was purchased in Florida.}      
\item{\code{purchInCalifornia}}{indicator that the vehicle was purchased in California.}   
\item{\code{purchInNorthCarolina}}{indicator that the vehicle was purchased in North Carolina.}
\item{\code{purchInArizona}}{indicator that the vehicle was purchased in Arizona.}           
\item{\code{purchInColorado}}{indicator that the vehicle was purchased in Colorado.}           
\item{\code{purchInSouthCarolina}}{indicator that the vehicle was purchased in South Carolina.}           
\item{\code{costAtPurch}}{acquisition cost paid for the vehicle at time of purchase.}        
\item{\code{onlineSale}}{indicator that the vehicle was purchased online.}         
\item{\code{warrantyCost}}{warranty cost in U.S. dollars.}
}       
}
\source{The ``Don't Get Kicked'' competition, \url{https://www.kaggle.com}.}
\examples{
library(HRW) ; data(carAuction)
\dontrun{
for (colNum in 3:10) 
{
   plot(jitter(carAuction[,colNum]),jitter(carAuction$IsBadBuy),pch = ".",
        xlab = names(carAuction)[colNum],ylab = "is car a bad buy?",col = "blue")
   readline("Hit Enter to continue.\n")
}
for (colNum in 11:51) 
{
   plot(jitter(carAuction[,colNum]),jitter(carAuction$IsBadBuy),pch = ".",
        xlab = names(carAuction)[colNum],ylab = "is car a bad buy?",col = "blue")
   readline("Hit Enter to continue.\n")
}
}
}
